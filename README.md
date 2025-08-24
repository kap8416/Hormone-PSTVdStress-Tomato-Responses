# 🌱 Evolutionary Reconstruction of Hormone-Driven Master Regulators in PSTVd–Tomato Responses

**Code Contributors:** Katia Aviña-Padilla · Octavio Zambada · Luis Hernández · Manuel Barrios
  **Supervision** Katia Aviña-Padilla · Maribel Hernández Rosales
**Last Update:** 2025-08-24  

This repository hosts code and resources for the **phylogenomic and systems-level analysis of PSTVd–tomato interactions**, with a primary focus on identifying the **evolutionary reconstruction of gene regulatory networks (GRNs)** and **Master Transcriptional Regulators (MTRs) linked to hormone signaling**.  

1. Transcriptomic GRN deconvolution and MRA are implemented in R using **corto**, based on our original script (*scripts/corto_PSTVdTomato.R*).  
2. Evolutionary reconstruction of MTRs and their regulons is performed using **REvolutionH-tl** ([More information here](https://pypi.org/project/revolutionhtl/)).  
3. Ortholog integration, comparative GRN analysis, rewiring, and network visualization are implemented in R, Bash, and Python.  

---

## 📂 Repository Structure

```
Hormone-PSTVdStress-Tomato-Responses/
│
├── data/
│   ├── GSE106912/              # Leaf transcriptome (PSTVd mild, severe, and healthy)
│   ├── GSE111736/              # Root transcriptome (PSTVd mild, severe, and healthy)
│   ├── Tomato_TFs.txt          # TF list (PlantTFDB-derived)
│   └── Tomato_pheno.csv        # Sample metadata (C, M, S)
│
├── scripts/
│   ├── corto_PSTVdTomato.R            # GRN deconvolution + MRA (corto)
│   ├── MTRs_Network.py                # Radial visualization of MTRs (shared vs unique targets)
│   ├── TomatoOrthologs_5Sets.R        # Integrates hormone target lists with ortholog tables + alluvial plot
│   ├── Regulatory-Network-Plotter.R   # Multi-TSV visualization of tomato-PSTVd networks with igraph
│   └── analyze_rewiring.R             # Rewiring metrics between networks (edge/node/global)
│
├── results/
│   ├── networks/               # net_corto_prom.txt / .sif (for Cytoscape)
│   ├── mra_plots/              # PDF figures from mraplot()
│   ├── tables/                 # Top MTRs per comparison; summary TSVs
│   └── rewiring/               # Gained/lost/signswitch/centrality metrics
│
└── README.md
```

---
![PipelinePlosComputational](https://github.com/user-attachments/assets/5587dd45-9fd5-4c3d-8ba6-b3304d9d304f)


## 🔬 Analysis Pipeline (3 modules)

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

## 📊 Expected Results
- Regulatory networks of tomato under PSTVd infection with ranked **MTRs**.  
- Comparative phylogenomics and ortholog integration of hormone targets.  
- **Rewiring analysis**: conserved vs rewired edges, centrality changes, and global similarity metrics.  
- Publication-quality outputs: MRA plots, radial MTR network, alluvial plots, and rewiring visualizations.  

---

## 📖 References
- Mercatelli D., Lopez-Garcia G., Giorgi F. M. (2020). *corto: a lightweight R package for gene network inference and master regulator analysis.* **Bioinformatics**, 36(12):3916–3917. doi:10.1093/bioinformatics/btaa223  
- Aviña-Padilla K., Zambada-Moreno O., Herrera-Oropeza G. E., et al. (2022). *Insights into the Transcriptional Reprogramming in Tomato Response to PSTVd Variants Using Network Approaches.* **Int J Mol Sci**, 23(11):5983. doi:10.3390/ijms23115983
-  Ramírez-Rafael J. A., Korchmaros A., Aviña-Padilla K., López-Sánchez A., España-Tinajero A. A., Hellmuth M., Stadler P. F., and Hernandez-Rosales M. (2024) REvolutionH-tl: Reconstruction of Evolutionary Histories tool. In Comparative Genomics: 21st International Conference, RECOMB-CG 2024, Boston, MA, USA, April 27–28, 2024, Proceedings. Springer-Verlag, Berlin, Heidelberg, 89–109. https://doi.org/10.1007/978-3-031-58072-7_5
-  
- GEO datasets: **GSE106912**, **GSE111736**  

---

## ✨ Citation
If you use this repository or code, please cite:  
> Aviña-Padilla K., et al. **Hormone Master Regulators Under PSTVd Stress: A Phylogenomic and Network View in Solanaceae.** (2025).  
