# 🌱 Evolutionary Reconstruction of Master Transcriptional Regulation of Hormone PSTVd-Tomato Responses

**Code Contributors:**  Katia Aviña-Padilla · Octavio Zambada ·Luis Hernández · Manuel Barrios  
**Last Update:** 2025-08-21

This repository hosts code and resources for the **Phylogenomic Insights from PSTVd–Tomato Interactions** study, with a primary focus on identifying the **Evolutionary Reconstruction of GRNs**.  

1) The transcriptomic GRN deconvolution and MRA are implemented in R using **corto**, based on our original script (see *scripts/CORTO_TOMATO-VIROID.R*).

The evolutionary reconstruction of MTRs and their corresponding regulon is performed using REvolutionH-tl [More information here.](https://pypi.org/project/revolutionhtl/)


---

## 📂 Repository Structure

```
Hormone-PSTVdStress-Tomato-Responses/
│
├── data/
│   ├── GSE106912/              # Leaf transcriptome (PSTVd mild, PSTVd severe and Healthy mock-inoculated)
│   ├── GSE111736/              # Root transcriptome (PSTVd mild, PSTVd severe and Healthy mock-inoculated)
│   ├── Tomato_TFs.txt          # TF list (PlantTFDB-derived)
│   └── Tomato_pheno.csv        # Sample metadata (C, M, S)
│
├── scripts/
│   ├── CORTO_TOMATO-VIROID.R   # GRN deconvolution + MRA (corto)
│   ├── MTRs_Network.py         # Radial network viz (5 MTR hubs; shared/unique targets)
│   └── TomatoOrthologs_.R      # Integrates 5 hormone target lists with ortholog tables + alluvial plot
│
├── results/
│   ├── networks/               # net_corto*.txt / .sif (for Cytoscape)
│   ├── mra_plots/              # PDF figures from mraplot()
│   └── tables/                 # Top MTRs per comparison; summary TSVs
│
└── README.md
```

---

## 🔬 Analysis Pipeline (3 modules)

### 1) **PSTVd GRN Construction (Transcriptomics → corto → MRA)**  
Input: GEO datasets **GSE106912** (leaves) and **GSE111736** (roots).  
Method: RMA normalization (affy) → **corto** network inference → **Master Regulator Analysis** (MRA) for C vs S, C vs M, S vs M.  
Output: regulons, ranked MTRs, and PDF plots; exportable network files for Cytoscape.

### 2) **Phylogenomics Analysis (REvolutionH-tl)**  
Input: FASTA sets and (for reconciliation) an **NHX species tree**. 
  -Proteomes from Sol Genomics Network [FTP Site](https://solgenomics.net/) 
We selected three wild and four domesticated hosts
S. pennellii, S. pimpinellifolium, C. annuum var. glabriusculum
S. lycopersicum, S. lycopersicum var. cerasiforme, C. annuum, S. tuberosum

Steps (used here): 1–2 (alignments, best hits & orthogroups) → filter → 3–4 (gene trees, duplication resolution) → 6 (gene–species reconciliation).  
Outputs: best hits, orthogroups, gene trees, resolved trees, reconciled summaries and visual reports.

### 3) **Evolutionary Reconstruction of GRNs**  
### 2) **Evolutionary Reconstruction of GRNs (ortholog integration and comparative analysis)**
This step combines **`scripts/TomatoOrthologs_5Sets.R`** with additional custom Bash and R scripts.

**Input**
- Hormone target lists (ABA, Auxin, Ethylene activator/repressor, MYC2/Jasmonate).
- Ortholog tables (`orthologs*.tsv`).
- Regulatory maps from PlantTFDB:
  - *S. lycopersicum* PlantTFDB (https://planttfdb.gao-lab.org)
  - *S. pennellii* PlantTFDB

**Task**
- Map hormone-specific targets to orthologs and orthogroups across species.
- Use **PlantTFDB regulatory maps** to reconstruct GRNs in domesticated (*S. lycopersicum*) and wild (*S. pennellii*) tomato.
- Identify:
  - **Conserved regulatory elements** across species.
  - **Rewiring events** between healthy vs diseased states.
  - **Conserved regulators** specific to the *S. pennellii* network.

**How to run**
```r
# Ortholog integration (R):
source("TomatoOrthologs_5Sets.R")

# Additional comparative GRN scripts (Bash/R):
bash comparative_grn_analysis.sh
Rscript analyze_rewiring.R


---

## ⚙️ Requirements

### R
- **R ≥ 4.0** and packages: `corto`, `affy`, `GEOquery`, `dplyr`, `ggplot2`
```r
install.packages(c("affy","GEOquery","dplyr","ggplot2"))
# corto from GitHub:
devtools::install_github("federicogiorgi/corto")
```

### Python (for network visualization)
- **Python 3.8+**, `pandas`, `networkx`, `matplotlib`
```bash
python -m pip install pandas networkx matplotlib
```

---

## 🚀 Usage

### A) GRN deconvolution & MRA (R, corto)
From `scripts/`:
```r
source("CORTO_TOMATO-VIROID.R")  # downloads/reads GEO, RMA, builds regulons, runs MRA, writes outputs
```
Key outputs:
- `results/networks/net_corto_prom.txt` 
- `results/mra_plots/cortoMRS_*.pdf`
- `results/tables/results M vs S.txt`

### B) Radial MTR network figure (Python)
The script expects a 2-column edge list (`source`, `target`) named **`net_corto_prom.txt`** (separator auto-detected).  
Run from repo root or `scripts/`:
```bash
python scripts/MTRs_Network.py
```
It produces `MTRs_Network.png/.svg/.pdf` with:
- 5 hub MTRs (hub fill = hormone; border = expression state, ERF_C_5 induced=red, others repressed=blue),
- shared targets colored by hormone combo (EJ, EA, AJ, AEJ),
- unique targets colored by their hub,
- labels for AEJ targets: **Solyc06g008330, Solyc08g075450, Solyc10g044780**.

### C) Hormone target ortholog integration (R)
From `scripts/`:
```r
source("TomatoOrthologs_5Sets.R")
```
Inputs expected in `~/Downloads/` (or adjust `setwd()`):  
`ABA_interactorsEdit*.txt`, `Auxin_interactorsEdit*.txt`, `Ethylene_*interactorsEdit*.txt`, `MYCinteractorsEdit.txt`, and `orthologs*.tsv`.  
Outputs (to `HormoneOutputs_5sets/`):  
- Per-hormone ortholog tables, counts, duplicated hits, merged summary `HormoneOrthologs_5sets.txt`, and an **alluvial plot** `HormoneOrthologs_5sets_alluvial.png`.

### D) REvolutionH-tl quick commands (phylogenomics)
```bash
# Steps 1–2 on FASTA folder:
python3 -m revolutionhtl -steps 1 2 -F fasta_files

# Filter best hits/orthogroups as needed, then:
python3 -m revolutionhtl -steps 3 4 --best_h filtered_best_hits.tsv -D tl_project.distances.tsv

# Reconciliation (Step 6) with NHX species tree:
python3 -m revolutionhtl -steps 6 -T tl_project.resolved_trees.tsv -S species_tree.nhx
```
Tip: use `python3 -m revolutionhtl.plot_summary` and `python3 -m revolutionhtl.plot_reconciliation <orthogroup_ID>` for visual summaries.

---

## 📊 Expected Results
- Tomato GRNs under PSTVd infection, **ranked MTRs** per comparison, and Cytoscape-ready networks.  
- Evolutionary context for regulators (orthogroups, reconciled gene/species histories).  
- Publication-quality figures: MRA PDFs, radial MTR network, and alluvial ortholog plot.

---

## 📖 References
- Mercatelli D., Lopez-Garcia G., Giorgi F. M. (2020). *corto: a lightweight R package for gene network inference and master regulator analysis.* **Bioinformatics**, 36(12):3916–3917. DOI: 10.1093/bioinformatics/btaa223  
- Aviña-Padilla K, Zambada-Moreno O, Herrera-Oropeza GE, et al. Insights into the Transcriptional Reprogramming in Tomato Response to PSTVd Variants Using Network Approaches. Int J Mol Sci. 2022;23(11):5983. Published 2022 May 26. doi:10.3390/ijms23115983  
- GEO datasets: **GSE106912**, **GSE111736**

---

## ✨ Citation
If you use this repository or code, please cite:
> Aviña-Padilla K. et al **Hormone Master Regulators Under PSTVd Stress: A Phylogenomic and Network View in Solanaceae.** (2025).
