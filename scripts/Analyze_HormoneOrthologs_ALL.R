# Analyze_HormoneOrthologs_ALL.R
# Inputs : HormoneOutputs_5sets/HormoneOrthologs.txt
# Outputs: HormoneOutputs_5sets/Analysis/* (tables + figures)
# Author: Katia Aviña Padilla

# -------------------- Setup --------------------
options(repos = c(CRAN = "https://cloud.r-project.org"))
pkgs <- c("dplyr","tidyr","readr","stringr","ggplot2","pheatmap","ComplexUpset","patchwork","tibble")
inst <- setdiff(pkgs, rownames(installed.packages()))
if (length(inst)) install.packages(inst, quiet = TRUE, dependencies = TRUE)
invisible(lapply(pkgs, function(p) suppressPackageStartupMessages(library(p, character.only = TRUE))))

root_dir <- getwd()
in_file  <- file.path(root_dir, "HormoneOutputs_5sets", "HormoneOrthologs.txt")
stopifnot(file.exists(in_file))
out_dir  <- file.path(root_dir, "HormoneOutputs_5sets", "Analysis")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -------------------- Read & clean --------------------
df <- readr::read_tsv(in_file, show_col_types = FALSE)

df <- df |>
  dplyr::mutate(
    Ortholog_Status = ifelse(is.na(Gen_b) | Gen_b == "NA", "No Ortholog", "Has Ortholog"),
    Species_with_Ortholog = ifelse(Ortholog_Status == "Has Ortholog", Species_with_Ortholog, "None")
  )

nice_names <- c(
  "Tomatopimp"        = "S. pimp.",
  "TomatoPennelli"    = "S. pennellii",
  "TomatoCerasiforme" = "Tomato cherry",
  "Potato"            = "Potato",
  "CaChiltepin"       = "Wild pepper",
  "CaZL1"             = "Dom. pepper",
  "None"              = "None"
)

df <- df |>
  dplyr::mutate(
    Species_with_Ortholog = gsub("\\.pep$", "", Species_with_Ortholog),
    Species_with_Ortholog = dplyr::recode(Species_with_Ortholog, !!!nice_names, .default = Species_with_Ortholog),
    Hormone = factor(Hormone, levels = c("ABA","Auxin","Ethylene","Jasmonic Acid")),
    Ortholog_Status = factor(Ortholog_Status, levels = c("Has Ortholog","No Ortholog"))
  )

# -------------------- Palette (global) --------------------
cols_h <- c(
  "ABA"            = "gold",
  "Auxin"          = "hotpink2",
  "Ethylene"       = "lightblue",
  "Jasmonic Acid"  = "palegreen2",
  "No Ortholog"    = "grey20"
)

species_levels <- c("S. pimp.","S. pennellii","Tomato cherry","Potato","Dom. pepper","Wild pepper")
species_levels <- species_levels[species_levels %in% unique(df$Species_with_Ortholog)]

# -------------------- 1) % Conservation by species × hormone --------------------
species_conservation <- df |>
  dplyr::group_by(Hormone) |>
  dplyr::summarise(TotalGenes = dplyr::n_distinct(en_a), .groups="drop") |>
  dplyr::left_join(
    df |>
      dplyr::filter(Ortholog_Status == "Has Ortholog") |>
      dplyr::distinct(Hormone, en_a, Species_with_Ortholog) |>
      dplyr::count(Hormone, Species_with_Ortholog, name = "GenesWithOrtholog"),
    by = "Hormone"
  ) |>
  dplyr::mutate(Percent = 100 * GenesWithOrtholog / pmax(TotalGenes, 1))

readr::write_tsv(species_conservation, file.path(out_dir, "species_conservation_percent.tsv"))

mat_heat <- species_conservation |>
  dplyr::select(Hormone, Species_with_Ortholog, Percent) |>
  tidyr::pivot_wider(names_from = Hormone, values_from = Percent, values_fill = 0) |>
  tibble::column_to_rownames("Species_with_Ortholog") |>
  as.matrix()

pheatmap::pheatmap(
  mat_heat, cluster_rows = TRUE, cluster_cols = FALSE,
  display_numbers = TRUE, number_format = "%.1f",
  filename = file.path(out_dir, "Heatmap_percent_conservation.png"),
  width = 6, height = 4
)

# -------------------- 2) Global conservation per hormone --------------------
hormone_conservation <- df |>
  dplyr::group_by(Hormone) |>
  dplyr::summarise(
    TotalGenes = dplyr::n_distinct(en_a),
    Conserved  = dplyr::n_distinct(en_a[Ortholog_Status == "Has Ortholog"]),
    Percent    = 100 * Conserved / pmax(TotalGenes, 1),
    .groups = "drop"
  )
readr::write_tsv(hormone_conservation, file.path(out_dir, "hormone_conservation_global.tsv"))

p1 <- ggplot2::ggplot(hormone_conservation,
                      ggplot2::aes(x = Hormone, y = Percent, fill = Hormone)) +
  ggplot2::geom_col(width = 0.65, color = "grey20") +
  ggplot2::scale_fill_manual(values = cols_h, breaks = names(cols_h), name = "Hormone", guide = "none") +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", Percent)),
                     vjust = -0.6, size = 3.7) +
  ggplot2::coord_cartesian(ylim = c(0, 105), expand = FALSE) +
  ggplot2::labs(x = NULL, y = "% genes with an ortholog") +
  ggplot2::theme_minimal(base_size = 12)
ggplot2::ggsave(file.path(out_dir, "Bar_percent_conservation_by_hormone.png"),
                p1, width = 7, height = 4, dpi = 300, bg = "white")

# -------------------- 3) Duplications --------------------
dup_by_species <- df |>
  dplyr::filter(Ortholog_Status == "Has Ortholog") |>
  dplyr::group_by(Hormone, Species_with_Ortholog) |>
  dplyr::summarise(
    N_genes = dplyr::n_distinct(en_a),
    N_dup   = sum(Duplicated_Ortholog, na.rm = TRUE),
    Percent_dup = 100 * N_dup / pmax(N_genes, 1),
    .groups = "drop"
  )
readr::write_tsv(dup_by_species, file.path(out_dir, "duplications_by_species.tsv"))

p2 <- ggplot2::ggplot(dup_by_species,
                      ggplot2::aes(x = Species_with_Ortholog, y = Percent_dup, fill = Hormone)) +
  ggplot2::geom_col(width = 0.65, color = "grey20") +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", Percent_dup)),
                     vjust = -0.4, size = 3) +
  ggplot2::facet_wrap(~ Hormone, ncol = 2) +
  ggplot2::scale_fill_manual(values = cols_h, breaks = names(cols_h), name = "Hormone", guide = "none") +
  ggplot2::labs(x = "Species", y = "% genes with co-orthologs") +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
ggplot2::ggsave(file.path(out_dir, "Bar_percent_duplications_by_species.png"),
                p2, width = 10, height = 6, dpi = 300, bg = "white")

# -------------------- 4) Chi-square tests --------------------
tab_h <- df |>
  dplyr::group_by(Hormone, Ortholog_Status) |>
  dplyr::summarise(n = dplyr::n_distinct(en_a), .groups = "drop") |>
  tidyr::pivot_wider(names_from = Ortholog_Status, values_from = n, values_fill = 0)
chisq_h <- chisq.test(as.matrix(tab_h[, c("Has Ortholog","No Ortholog")]))
sink(file.path(out_dir, "chisq_conservation_by_hormone.txt"))
print(tab_h); print(chisq_h); sink()

total_genes <- df |> dplyr::summarise(TotalGenes = dplyr::n_distinct(en_a)) |> dplyr::pull()
tab_s <- df |>
  dplyr::filter(Ortholog_Status == "Has Ortholog") |>
  dplyr::distinct(en_a, Species_with_Ortholog) |>
  dplyr::count(Species_with_Ortholog, name = "ConservedGenes") |>
  dplyr::mutate(NotConserved = pmax(total_genes - ConservedGenes, 0))
chisq_s <- chisq.test(as.matrix(tab_s[, c("ConservedGenes","NotConserved")]))
sink(file.path(out_dir, "chisq_conservation_by_species.txt"))
print(tab_s); print(chisq_s); sink()

# -------------------- 5) Core & species-specific --------------------
species_set <- setdiff(unique(df$Species_with_Ortholog), "None")
core_genes <- df |>
  dplyr::filter(Ortholog_Status == "Has Ortholog") |>
  dplyr::distinct(Hormone, en_a, Species_with_Ortholog) |>
  dplyr::group_by(Hormone, en_a) |>
  dplyr::summarise(n_sp = dplyr::n_distinct(Species_with_Ortholog), .groups = "drop") |>
  dplyr::filter(n_sp == length(species_set))
readr::write_tsv(core_genes, file.path(out_dir, "core_genes_by_hormone.tsv"))

species_specific <- df |>
  dplyr::filter(Ortholog_Status == "Has Ortholog") |>
  dplyr::distinct(Hormone, en_a, Species_with_Ortholog) |>
  dplyr::group_by(Hormone, en_a) |>
  dplyr::summarise(n_sp = dplyr::n_distinct(Species_with_Ortholog),
                   only_sp = paste(sort(unique(Species_with_Ortholog)), collapse = ","),
                   .groups = "drop") |>
  dplyr::filter(n_sp == 1)
readr::write_tsv(species_specific, file.path(out_dir, "species_specific_genes_by_hormone.tsv"))

# -------------------- 6) UpSet (2x2; puntos vacíos grises cuando sea posible) --------------------
presence_long_has <- df |>
  dplyr::filter(Ortholog_Status == "Has Ortholog") |>
  dplyr::distinct(Hormone, en_a, Species_with_Ortholog) |>
  dplyr::mutate(present = 1L)

# Detecta si tu ComplexUpset soporta 'empty_intersections'
has_empty_arg <- "empty_intersections" %in% names(formals(ComplexUpset::intersection_matrix))

make_upset_hormone <- function(horm) {
  sub <- presence_long_has[presence_long_has$Hormone == horm, ] |>
    dplyr::select(en_a, Species_with_Ortholog, present) |>
    tidyr::pivot_wider(id_cols = en_a, names_from = Species_with_Ortholog,
                       values_from = present, values_fill = 0L)
  if (nrow(sub) == 0) return(NULL)
  for (sp in species_levels) if (!(sp %in% names(sub))) sub[[sp]] <- 0L
  sub <- sub[, species_levels, drop = FALSE]
  
  col_h <- cols_h[[as.character(horm)]]
  
  if (has_empty_arg) {
    # Ruta moderna: puntos conectados al color hormonal y vacíos grises
    p <- ComplexUpset::upset(
      sub, species_levels, name = "Genes", min_size = 1,
      base_annotations = list(
        "Counts" = ComplexUpset::intersection_size(
          mapping = ggplot2::aes(fill = I(col_h)), text = list()
        )
      ),
      set_sizes = ComplexUpset::upset_set_size(
        mapping = ggplot2::aes(fill = I(col_h))
      ),
      matrix = ComplexUpset::intersection_matrix(
        geom = ggplot2::geom_point(size = 2.6, shape = 16, color = col_h),
        empty_intersections = ggplot2::geom_point(size = 2.6, shape = 16, color = "grey70")
      )
    ) +
      ggplot2::scale_fill_identity()
  } else {
    # Ruta compatible: puntos conectados al color hormonal; vacíos quedan con el estilo por defecto
    # (en muchas versiones ya es un gris semi-transparente/contorno)
    p <- ComplexUpset::upset(
      sub, species_levels, name = "Genes", min_size = 1,
      base_annotations = list(
        "Counts" = ComplexUpset::intersection_size(
          mapping = ggplot2::aes(fill = I(col_h)), text = list()
        )
      ),
      set_sizes = ComplexUpset::upset_set_size(
        mapping = ggplot2::aes(fill = I(col_h))
      ),
      matrix = ComplexUpset::intersection_matrix(
        geom = ggplot2::geom_point(size = 2.6, shape = 16, color = col_h)
      )
    ) +
      ggplot2::scale_fill_identity()
  }
  
  p +
    ggplot2::labs(title = paste0(horm)) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold"),
      axis.title.x     = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_blank(),
      axis.ticks.x     = ggplot2::element_blank(),
      axis.title.x.top = ggplot2::element_blank(),
      axis.text.x.top  = ggplot2::element_blank(),
      axis.ticks.x.top = ggplot2::element_blank()
    )
}

plots <- lapply(levels(df$Hormone), make_upset_hormone)
plots <- plots[!vapply(plots, is.null, logical(1))]
if (length(plots)) {
  fig <- patchwork::wrap_plots(plots, ncol = 2)
  ggplot2::ggsave(file.path(out_dir, "UpSet_species_faceted_noXnums.png"),
                  fig, width = 12, height = 8, dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(out_dir, "UpSet_species_faceted_noXnums.pdf"),
                  fig, width = 12, height = 8, bg = "white")
  ggplot2::ggsave(file.path(out_dir, "UpSet_species_faceted_noXnums.svg"),
                  fig, width = 12, height = 8, bg = "white")
}

message("✅ Analysis completed. Files written to: ", out_dir)
