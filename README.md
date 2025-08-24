# 🌱 Evolutionary Reconstruction of Hormone-Driven Master Regulators in PSTVd–Tomato Responses

## **Code Contributors:** Katia Aviña-Padilla · Octavio Zambada · Luis Hernández · Manuel Barrios

**Supervision** Katia Aviña-Padilla · Maribel Hernández Rosales

**Last Update:** 2025-08-24  

This repository hosts code and resources for the **phylogenomic and systems-level analysis of PSTVd–tomato interactions**, with a primary focus on identifying the **evolutionary reconstruction of gene regulatory networks (GRNs)** and **Master Transcriptional Regulators (MTRs) linked to hormone signaling**.  

1. Transcriptomic GRN deconvolution and MRA are implemented in R using **corto**, based on our original script (*scripts/corto_PSTVdTomato.R*).  
2. Evolutionary reconstruction of MTRs and their regulons is performed using **REvolutionH-tl** ([More information here](https://pypi.org/project/revolutionhtl/)).  
3. Ortholog integration, comparative GRN analysis, rewiring, and network visualization are implemented in R, Bash, and Python.
![PipelinePlosComputational](https://github.com/user-attachments/assets/5587dd45-9fd5-4c3d-8ba6-b3304d9d304f) 

---

## 📂 Repository Structure

```
Hormone-PSTVdStress-Tomato-Responses/
│
├── data/
│   └── Inputs/
│       └── Hormone_5sets/                 # Input hormone-target interactors
│           ├── ABA_interactorsEdit_sinPunto.txt
│           ├── Auxin_interactorsEdit_sinPunto.txt
│           ├── Ethylene_Solyc02g077370.1_interactors.txt
│           ├── Ethylene_Solyc02g093130.1_interactors.txt
│           ├── MYCinteractorsEdit.txt
│           └── orthologs_Solanaceae.tsv   # Phylogenomic orthologs
│
├── scripts/
│   ├── corto_PSTVdTomato.R                # GRN inference + MRA (corto)
│   ├── extract_regulons.py                # Extract regulons from corto networks
│   ├── MTRs_Network.py                    # Radial visualization of MTRs
│   ├── Analyze_HormoneOrthologs_ALL.R     # Ortholog conservation (heatmaps, UpSet, χ²)
│   ├── Network_visualizer.R               # igraph-based network plotting
│   ├── networkx_pennelli.ipynb            # S. pennellii network visualization (Python + networkx)
│   ├── tomato_healthy.Rmd                 # RMarkdown analysis of healthy tomato
│   ├── int-tomato.Rmd                     # Integrative RMarkdown tomato transcriptomics
│   ├── tomatohealthy_network.ipynb        # Healthy tomato networks (Jupyter)
│   └── healthy-disease.ipynb              # Healthy vs diseased network comparisons (Jupyter)
│
├── results/
│   └── HormoneOutputs_5sets/
│       └── Analysis/
│           ├── HormoneOrthologs.txt / HormoneOrthologs_5sets.txt
│           ├── *_OrthoSpecies.txt / *_numOrthoLikeCounts.txt
│           ├── HormoneOrthologs_Alluvial.(png|pdf|svg)
│           ├── duplicated_genesSpecies_*.txt
│           ├── sessionInfo.txt
│           ├── MTRs_Network_v7.(png|pdf|svg)
│           ├── MTRs_Selected.csv
│           ├── MTRs_orthologs_(heatmap|presence_absence.png)
│           ├── MTRs_orthologs_(summary|table).tsv
│           ├── net_corto_prom.txt
│           ├── regulon_Solyc*.csv / .txt   # Per-MTR regulons
│           └── regulon_summary_five_MTRs.csv
│
└── README.md

```

---
## 🔬 Analysis Pipeline (3 modules)

### **1) PSTVd GRN Construction (Transcriptomics → corto → MRA)**  
- **Input:** GEO datasets [GSE106912](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106912) (leaves) and [GSE111736](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111736) (roots).  
- **Method:** RMA normalization (affy) → **corto** network inference → **MRA** for C vs S, C vs M, S vs M.  
- **Main output:** `net_corto_prom.txt` (input for further analyses).  
- **Additional outputs:** regulons (`regulon_Solyc*.txt/csv`), ranked MTRs, MRA plots (PDF), Cytoscape-ready networks.  

---

### **2) Phylogenomics Analysis (REvolutionH-tl)**  
- **Input:** protein FASTA files from the [Sol Genomics Network](https://solgenomics.net/) and an **NHX species tree**.  
- **Species analyzed:**  
  - *Wild:* *S. pennellii*, *S. pimpinellifolium*, *C. annuum* var. *glabriusculum*  
  - *Domesticated:* *S. lycopersicum*, *S. cerasiforme*, *C. annuum*, *S. tuberosum*  

- **Steps used:**  
  1–2 (alignments, best hits & orthogroups) → filtering → 3–4 (gene trees, duplication resolution) → 6 (gene–species reconciliation).  

- **Outputs:** best hits, orthologs, orthogroups, gene trees, resolved trees, reconciled summaries and visual reports.  
- **Main output:** `orthologs_Solanaceae.tsv`.  

---

### **3) Evolutionary Reconstruction of GRNs (ortholog integration, comparative analysis, rewiring & plotting)**  
This module integrates multi-species orthologs with GRN inference results to reconstruct and compare hormone-driven regulatory programs.  

- **Inputs:**  
  - Hormone target lists (ABA, Auxin, Ethylene activator/repressor, MYC2/Jasmonate).  
  - Ortholog tables (`orthologs*.tsv`).  
  - PlantTFDB regulatory maps (*S. lycopersicum*, *S. pennellii*).  
  - Species/condition-specific edge lists (TSV: `src, tgt, src_mode, tgt_mode`).  

- **Tasks:**  
  1. **Ortholog integration (R):** map hormone targets to orthologs/orthogroups; detect duplicates; generate alluvial and UpSet summaries.  
  2. **Comparative GRN construction (Bash + R):** build per-species/per-condition TSV networks and concatenate by phylogeny.  
  3. **Rewiring analysis (R):** quantify edge/node/global changes between states and species:  
     - *Edge-level:* gained/lost/sign-switch events.  
     - *Node-level (MTRs):* Δ-outdegree, Δ-betweenness, Δ-regulon size.  
     - *Global:* Jaccard similarity and edit distance.  
  4. **Visualization (R):** plot comparable networks (PNG/SVG) using **Regulatory-Network-Plotter** with a unified layout.  


**How to run**
```r
# 1) Ortholog integration
source("scripts/TomatoOrthologs_5Sets.R")
```

```bash
# 2) Comparative networks
```

```bash
# 3) Rewiring analysis 

```r
# 4) Visualization (edit folder_path inside script)
Rscript scripts/Regulatory-Network-Plotter.R
```

**Output**
- Per-hormone ortholog tables (`*_OrthoSpecies.txt`).  
- Merged summary (`HormoneOrthologs_5sets.txt`).  
- Alluvial plot (`HormoneOrthologs_5sets_alluvial.png`).  
- Species/condition TSV networks (per hormone) + concatenated phylogeny networks.  
- **Rewiring reports:** Jaccard, gained, lost, sign-switch edges, centrality deltas.  
- Publication-ready network figures (`.png`, `.svg`).  

---

## ⚙️ Requirements

### R
- **R ≥ 4.0**  
- Packages: `corto`, `affy`, `GEOquery`, `dplyr`, `ggplot2`, `igraph`, `readr`  
```r
install.packages(c("affy","GEOquery","dplyr","ggplot2","igraph","readr"))
devtools::install_github("federicogiorgi/corto")
```

### Python (for MTR network visualization)
- **Python 3.8+**, `pandas`, `networkx`, `matplotlib`  
```bash
python -m pip install pandas networkx matplotlib
```

---

## 🚀 Usage

- **GRN deconvolution (R, corto):**  
  ```r
  source("scripts/CORTO_TOMATO-VIROID.R")
  ```  

- **Radial MTR network (Python):**  
  ```bash
  python scripts/MTRs_Network.py
  ```  

- **Ortholog integration (R):**  
  ```r
  source("scripts/TomatoOrthologs_5Sets.R")
  ```  

- **Rewiring metrics (R):**  
  ```bash
  Rscript scripts/analyze_rewiring.R --netA net1.tsv --netB net2.tsv --out results/rewiring/comparison
  ```  

- **Network plotting (R):**  
  ```r
  Rscript scripts/Regulatory-Network-Plotter.R
  ```  

- **Phylogenomics (REvolutionH-tl):**  
  ```bash
  python3 -m revolutionhtl -steps 1 2 -F fasta_files
  python3 -m revolutionhtl -steps 3 4 --best_h filtered_best_hits.tsv -D tl_project.distances.tsv
  python3 -m revolutionhtl -steps 6 -T tl_project.resolved_trees.tsv -S species_tree.nhx
  ```

---

## 📊 Expected Outputs
- Ranked **MTRs** under PSTVd stress.  
- **Ortholog conservation** across Solanaceae: alluvial, UpSet, heatmaps, χ².  
- **Comparative networks**: per hormone, per species, per condition.  
- **Rewiring analysis**: gained/lost edges, sign-switches, Δcentrality, Jaccard.  
- **Publication-quality figures**: radial MTR networks, ortholog integration plots, rewiring visualizations etc.  

---

## 📖 References
- Mercatelli D., Lopez-Garcia G., Giorgi F. M. (2020). *corto: a lightweight R package for gene network inference and master regulator analysis.* **Bioinformatics**, 36(12):3916–3917.  
- Aviña-Padilla K., Zambada-Moreno O., Herrera-Oropeza G. E., et al. (2022). *Insights into the Transcriptional Reprogramming in Tomato Response to PSTVd Variants Using Network Approaches.* **Int J Mol Sci**, 23(11):5983.  
- Ramírez-Rafael J. A., Korchmaros A., Aviña-Padilla K., et al. (2024). *REvolutionH-tl: Reconstruction of Evolutionary Histories tool.* In **RECOMB-CG 2024**, Springer.  
- Aviña-Padilla K., Zambada-Moreno O., Bustamante Castillo M., Barrios-Izás M. A., Hernández-Rosales M. (2025). *Evolutionary Reconstruction of Hormone-bHLH Regulatory Networks in Solanaceae: Phylogenomic Insights from PSTVd-Tomato Interactions.* **bioRxiv** 2025.03.14.643413. doi: [https://doi.org/10.1101/2025.03.14.643413](https://doi.org/10.1101/2025.03.14.643413)  
- GEO datasets: [GSE106912](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106912), [GSE111736](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111736).  


---

## ✨ Citation
If you use this repository or code, please cite:  
> Aviña-Padilla K., et al. **Hormone Master Regulators Under PSTVd Stress: A Phylogenomic and Network View in Solanaceae.** (2025).  
