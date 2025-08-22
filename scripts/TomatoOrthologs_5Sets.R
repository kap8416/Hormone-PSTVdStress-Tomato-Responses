## TomatoOrthologs_5Sets.R
## Authors: Octavio Zambada Moreno & Katia Aviña Padilla
## Maintainer (this revision): Katia Aviña Padilla
##
## Purpose:
## Build ortholog tables for five hormone-related gene sets in tomato, oriented
## so that Gen_a is always the tomato gene, and Gen_b holds ortholog IDs for each
## species. Co-orthologs for the same species are concatenated in Gen_b using a
## SAFE delimiter `;` (semicolon) — NEVER tabs — to keep TSV column integrity.
##
## Outputs (written under ./HormoneOutputs_5sets):
##  - ABA_OrthoSpecies.txt, Auxin_OrthoSpecies.txt, Ethylene_OrthoSpecies.txt, MYC2_OrthoSpecies.txt
##  - *_numOrthoLikeCounts.txt
##  - duplicated_genesSpecies_*.txt
##  - HormoneOrthologs.txt                             # normalized master table
##  - HormoneOrthologs_5sets.txt                       # (Gen_a, Gen_b, Hormone, Species_with_Ortholog, Duplicated_Ortholog)
##  - HormoneOrthologs_Alluvial.(png|svg|pdf)

# -------------------- Setup --------------------
setwd("/Users/katiaavinapadilla/Downloads")

options(repos = c(CRAN = "https://cloud.r-project.org"))
pkgs <- c("dplyr","tidyr","readr","stringr","ggplot2","ggalluvial","tibble","forcats")
inst <- setdiff(pkgs, rownames(installed.packages()))
if (length(inst)) install.packages(inst, quiet = TRUE, dependencies = TRUE)
invisible(lapply(pkgs, function(p) suppressPackageStartupMessages(library(p, character.only = TRUE))))

out_dir <- file.path(getwd(), "HormoneOutputs_5sets")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -------------------- Helpers --------------------
`%||%` <- function(a,b) if (!is.null(a)) a else b
strip_version <- function(x) sub("\\.[0-9]+$", "", x %||% "")
clean_id <- function(x) { x <- gsub("[[:space:]\u00A0]+", "", x %||% ""); strip_version(x) }

normalize_ortholog_cols <- function(df) {
  nm <- tolower(names(df)); names(df) <- nm
  ren <- list(
    a="^a$", b="^b$",
    species_a="^species[_]*a$", species_b="^species[_]*b$",
    og="^og$", normalized_bit_score="^normalized[_]*bit[_]*score$"
  )
  for (new in names(ren)) {
    hit <- grep(ren[[new]], nm, value = TRUE)
    if (length(hit) >= 1) names(df)[names(df) == hit[1]] <- new
  }
  stopifnot(all(c("a","b","species_a","species_b") %in% names(df)))
  df$a <- clean_id(df$a); df$b <- clean_id(df$b)
  df
}

read_gene_list_1col <- function(path) {
  x <- read.delim(path, header = FALSE, sep = "\t", quote = "", comment.char = "", check.names = FALSE)
  if (ncol(x) > 1) x <- x[,1,drop=FALSE]
  names(x)[1] <- "V1"
  x$V1 <- clean_id(x$V1)
  dplyr::filter(x, !is.na(V1) & V1 != "") |> dplyr::distinct()
}

# Orient to tomato side: Gen_a = tomato, Gen_b = ortholog (other species)
orient_to_tomato <- function(df) {
  r1 <- df |>
    dplyr::filter(species_b == "tomato") |>
    dplyr::transmute(Gen_a = b, Gen_b = a, Species_a = species_b, Species_b = species_a,
                     og = og %||% NA_character_, normalized_bit_score = normalized_bit_score %||% NA_real_)
  r2 <- df |>
    dplyr::filter(species_a == "tomato") |>
    dplyr::transmute(Gen_a = a, Gen_b = b, Species_a = species_a, Species_b = species_b,
                     og = og %||% NA_character_, normalized_bit_score = normalized_bit_score %||% NA_real_)
  dplyr::bind_rows(r1, r2) |> dplyr::distinct()
}

# Add "No Ortholog" rows for tomato genes absent from the oriented table
add_no_ortholog_rows <- function(tomato_list, oriented_tbl, hormone_label) {
  present <- unique(oriented_tbl$Gen_a)
  missing <- setdiff(tomato_list$V1, present)
  if (length(missing) == 0) {
    return(tibble::tibble(
      en_a = character(), Gen_b = character(), Ortholog_Status = character(),
      Species_with_Ortholog = character(), Hormone = character(), Duplicated_Ortholog = logical()
    ))
  }
  tibble::tibble(
    en_a = missing,
    Gen_b = NA_character_,
    Ortholog_Status = "No Ortholog",
    Species_with_Ortholog = "None",
    Hormone = hormone_label,
    Duplicated_Ortholog = FALSE
  )
}

# Collapse per species; mark co-orthologs; IMPORTANT: use ';' (NO tabs) inside Gen_b
collapse_by_species <- function(oriented_tbl, hormone_label) {
  if (nrow(oriented_tbl) == 0) {
    return(tibble::tibble(
      en_a=character(), Gen_b=character(), Ortholog_Status=character(),
      Species_with_Ortholog=character(), Hormone=character(), Duplicated_Ortholog=logical()
    ))
  }
  oriented_tbl |>
    dplyr::group_by(Gen_a, Species_b) |>
    dplyr::reframe(
      genes = sort(unique(Gen_b)),
      Gen_b_concat = paste(genes, collapse = ";"),        # << SAFE delimiter
      n_genes = dplyr::n_distinct(Gen_b)
    ) |>
    dplyr::mutate(
      en_a = Gen_a,
      Gen_b = Gen_b_concat,
      Ortholog_Status = ifelse(is.na(Gen_b) | Gen_b == "", "No Ortholog", "Has Ortholog"),
      Species_with_Ortholog = ifelse(Ortholog_Status == "Has Ortholog", Species_b, "None"),
      Hormone = hormone_label,
      Duplicated_Ortholog = n_genes > 1
    ) |>
    dplyr::select(en_a, Gen_b, Ortholog_Status, Species_with_Ortholog, Hormone, Duplicated_Ortholog) |>
    dplyr::ungroup()
}

# Per-hormone pipeline to produce a normalized, tomato-oriented table
build_for_hormone <- function(orth, list_path, hormone_label) {
  tom <- read_gene_list_1col(list_path)
  sub <- orth |> dplyr::filter(a %in% tom$V1 | b %in% tom$V1)
  oriented <- orient_to_tomato(sub)
  tbl_ok <- collapse_by_species(oriented, hormone_label)
  no_orth <- add_no_ortholog_rows(tom, oriented, hormone_label)
  dplyr::bind_rows(tbl_ok, no_orth)
}

# -------- Utilities for writing consistent TSVs --------
write_tsv <- function(df, fname) {
  dst <- file.path(out_dir, fname)
  readr::write_tsv(df, dst, na = "NA")
  message("Written: ", dst); dst
}

write_tsv_safe_5cols <- function(df, fname) {
  # Force exactly the five canonical columns for the 5sets export
  df2 <- df |>
    dplyr::transmute(
      Gen_a = en_a,
      Gen_b = Gen_b,
      Hormone = Hormone,
      Species_with_Ortholog = Species_with_Ortholog,
      Duplicated_Ortholog = Duplicated_Ortholog
    )
  write_tsv(df2, fname)
}

count_per_hormone <- function(df) {
  if (nrow(df) == 0) return(tibble::tibble(Gen = character(0), Count = integer(0)))
  df |> dplyr::count(Gen_a, name = "Count") |> dplyr::rename(Gen = Gen_a) |> dplyr::arrange(dplyr::desc(Count))
}

dup_detector <- function(df) {
  if (nrow(df) == 0) return(df)
  df |> dplyr::group_by(Gen_a, Species_b) |> dplyr::summarise(count = dplyr::n(), .groups = "drop") |> dplyr::filter(count > 1)
}

# ---------- join_with_orthologs (for raw per-species dumps) ----------
join_with_orthologs <- function(gene_df, orth) {
  genes <- dplyr::rename(gene_df, V1 = V1)
  sub <- orth |> dplyr::filter(a %in% genes$V1 | b %in% genes$V1)
  orient_to_tomato(sub) |>
    dplyr::select(Gen_a, Gen_b, Species_a, Species_b, og, normalized_bit_score) |>
    dplyr::distinct()
}

# -------------------- Inputs --------------------
ortho_file <- "orthologs_Solanaceae.tsv"
stopifnot(file.exists(ortho_file))
files <- list(
  ABA   = "ABA_interactorsEdit_sinPunto.txt",
  Auxin = "Auxin_interactorsEdit_sinPunto.txt",
  Eth1  = "Ethylene_Solyc02g077370.1_interactorsEdit_sinPunto.txt",
  Eth2  = "Ethylene_Solyc02g093130.1_interactorsEdit_sinPunto.txt",
  MYC2  = "MYCinteractorsEdit.txt"
)
stopifnot(all(file.exists(unlist(files))))

# -------------------- Load and build --------------------
orthologs <- read.delim(ortho_file, sep = "\t", header = TRUE, quote = "", comment.char = "", check.names = FALSE) |>
  normalize_ortholog_cols()

# Raw oriented dumps (useful for diagnostics)
lst_ABA   <- read_gene_list_1col(files$ABA)
lst_Aux   <- read_gene_list_1col(files$Auxin)
lst_Eth   <- dplyr::bind_rows(read_gene_list_1col(files$Eth1), read_gene_list_1col(files$Eth2)) |> dplyr::distinct()
lst_MYC2  <- read_gene_list_1col(files$MYC2)

ABA_tbl_raw <- join_with_orthologs(lst_ABA,  orthologs) |> dplyr::mutate(Hormone = "ABA")
Aux_tbl_raw <- join_with_orthologs(lst_Aux,  orthologs) |> dplyr::mutate(Hormone = "Auxin")
Eth_tbl_raw <- join_with_orthologs(lst_Eth,  orthologs) |> dplyr::mutate(Hormone = "Ethylene")
MYC_tbl_raw <- join_with_orthologs(lst_MYC2, orthologs) |> dplyr::mutate(Hormone = "Jasmonic Acid")

write_tsv(ABA_tbl_raw, "ABA_OrthoSpecies.txt")
write_tsv(Aux_tbl_raw, "Auxin_OrthoSpecies.txt")
write_tsv(Eth_tbl_raw, "Ethylene_OrthoSpecies.txt")
write_tsv(MYC_tbl_raw, "MYC2_OrthoSpecies.txt")
write_tsv(count_per_hormone(ABA_tbl_raw), "ABA_numOrthoLikeCounts.txt")
write_tsv(count_per_hormone(Aux_tbl_raw), "Auxin_numOrthoLikeCounts.txt")
write_tsv(count_per_hormone(Eth_tbl_raw), "Ethylene_numOrthoLikeCounts.txt")
write_tsv(count_per_hormone(MYC_tbl_raw), "MYC2_numOrthoLikeCounts.txt")
write_tsv(dup_detector(ABA_tbl_raw), "duplicated_genesSpecies_ABA.txt")
write_tsv(dup_detector(Aux_tbl_raw), "duplicated_genesSpecies_Auxin.txt")
write_tsv(dup_detector(Eth_tbl_raw), "duplicated_genesSpecies_Ethylene.txt")
write_tsv(dup_detector(MYC_tbl_raw), "duplicated_genesSpecies_MYC2.txt")

# -------------------- Normalized table (for plotting) --------------------
out_ABA <- build_for_hormone(orthologs, files$ABA,   "ABA")
out_Aux <- build_for_hormone(orthologs, files$Auxin, "Auxin")
out_Eth <- dplyr::bind_rows(
  build_for_hormone(orthologs, files$Eth1, "Ethylene"),
  build_for_hormone(orthologs, files$Eth2, "Ethylene")
) |> dplyr::distinct(en_a, Gen_b, Ortholog_Status, Species_with_Ortholog, Hormone, Duplicated_Ortholog, .keep_all = TRUE)
out_MYC <- build_for_hormone(orthologs, files$MYC2,  "Jasmonic Acid")

final_tbl <- dplyr::bind_rows(out_ABA, out_Aux, out_Eth, out_MYC) |>
  dplyr::arrange(Hormone, en_a, Species_with_Ortholog)

# --------- Friendly species names ----------
nice_names <- c(
  "Tomatopimp"        = "S. pimp.",
  "TomatoPennelli"    = "S. pennellii",
  "TomatoCerasiforme" = "Tomato cherry",
  "Potato"            = "Potato",
  "CaChiltepin"       = "Wild pepper",
  "CaZL1"             = "Dom. pepper",
  "None"              = "None"
)
final_tbl <- final_tbl |>
  dplyr::mutate(
    Species_with_Ortholog = gsub("\\.pep$", "", Species_with_Ortholog),
    Species_with_Ortholog = dplyr::recode(Species_with_Ortholog, !!!nice_names, .default = Species_with_Ortholog)
  )

# -------------------- Write normalized outputs --------------------
write_tsv(final_tbl, "HormoneOrthologs.txt")
write_tsv_safe_5cols(final_tbl, "HormoneOrthologs_5sets.txt")

# -------------------- Alluvial plot (fixed: FillClass after count) --------------------
plot_df <- final_tbl |>
  dplyr::mutate(
    Ortholog_Status = ifelse(is.na(Gen_b) | Gen_b == "NA" | Gen_b == "", "No Ortholog", "Has Ortholog"),
    Species_with_Ortholog = ifelse(Ortholog_Status == "Has Ortholog", Species_with_Ortholog, "None"),
    Duplicated_Ortholog = factor(Duplicated_Ortholog, levels = c(TRUE, FALSE)),
    Hormone = forcats::fct_relevel(Hormone, "ABA","Auxin","Ethylene","Jasmonic Acid"),
    Ortholog_Status = forcats::fct_relevel(Ortholog_Status, "Has Ortholog","No Ortholog"),
    Species_with_Ortholog = forcats::fct_relevel(
      Species_with_Ortholog,
      "S. pimp.","S. pennellii","Tomato cherry","Potato","Dom. pepper","Wild pepper","None"
    )
  ) |>
  dplyr::count(Hormone, Ortholog_Status, Species_with_Ortholog, Duplicated_Ortholog, name = "Count") |>
  dplyr::mutate(
    FillClass = ifelse(Ortholog_Status == "No Ortholog", "No ortholog", as.character(Hormone))
  )

p <- ggplot(
  plot_df,
  aes(axis1 = Hormone, axis2 = Ortholog_Status, axis3 = Species_with_Ortholog, axis4 = Duplicated_Ortholog, y = Count)
) +
  ggalluvial::geom_alluvium(aes(fill = FillClass), width = 1/6, alpha = 0.9) +
  ggalluvial::geom_stratum(width = 1/6, fill = "white", color = "black") +
  ggalluvial::geom_stratum(
    data = dplyr::filter(plot_df, Ortholog_Status == "No Ortholog"),
    width = 1/6, fill = "grey85", color = "black", alpha = 0.9
  ) +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 3, fill = "white", label.size = 0.2) +
  scale_x_discrete(
    limits = c("Hormone","Ortholog_Status","Species_with_Ortholog","Duplicated_Ortholog"),
    labels = c("Hormone","Orthologs Presence","Species","Co-ortholog")
  ) +
  scale_fill_manual(
    values = c("ABA"="gold","Auxin"="hotpink2","Ethylene"="lightblue","Jasmonic Acid"="palegreen2","No ortholog"="grey20"),
    breaks = c("ABA","Auxin","Ethylene","Jasmonic Acid","No ortholog"),
    name = "Hormone"
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1))

ggsave(file.path(out_dir, "HormoneOrthologs_Alluvial.png"), p, width = 12, height = 6.5, dpi = 300)
ggsave(file.path(out_dir, "HormoneOrthologs_Alluvial.svg"), p, width = 12, height = 6.5)
ggsave(file.path(out_dir, "HormoneOrthologs_Alluvial.pdf"), p, width = 12, height = 6.5)

# Repro info
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
message("🎉 Done! Files written to: ", out_dir)
