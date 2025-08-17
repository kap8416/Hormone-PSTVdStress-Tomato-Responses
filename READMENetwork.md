# Five MTRs Radial Network – Shared Axis Targets

This repository contains the script **`MTRs_network.py`**, which generates a radial network visualization of **five master transcriptional regulators (MTRs)** in *Solanum lycopersicum* under PSTVd infection.  
The figure highlights shared-axis targets regulated across different hormonal pathways (ABA, Auxin, Ethylene, and Jasmonate).

## Features
- Radial layout with **five MTRs** arranged as hubs.
- Regulons automatically extracted from `net_corto_prom.txt`.
- Targets colored by **shared regulatory axis**:
  - E–J shared (blue)  
  - E–A shared (green)  
  - A–J shared (orange)  
  - E–J–A triple shared (red, with **black border**)  
- Specific triple-shared genes are **labeled directly**.
- Legends for both **hormone-specific MTRs** and **shared target groups**.
- Output in **PNG**, **SVG**, and **PDF** formats for publication-ready figures.

## Dependencies
The script requires Python ≥ 3.8 and the following packages:
```bash
pip install networkx matplotlib pandas
```

## Usage
1. Place your regulatory edge list file in the repository root (default: `net_corto_prom.txt`).
2. Run the script:
   ```bash
   python Five_MTRs_radial_EJA_combos_v5.py
   ```
3. The following outputs will be generated:
   - `Five_MTRs_radial_EJA_combos_v5.png`
   - `Five_MTRs_radial_EJA_combos_v5.svg`
   - `Five_MTRs_radial_EJA_combos_v5.pdf`

## Input File Format
The script expects a **tab- or comma-separated file** with at least two columns:  
- **Source** (TF/MTR)  
- **Target** (regulated gene)

Example:
```text
source  target
Solyc08g076930.1  Solyc06g008330
Solyc02g093130.1  Solyc08g075450
Solyc04g009440.2  Solyc10g044780
```

## Example Output
The generated figure shows:
- Five MTRs arranged in a circle.  
- Targets around each regulator.  
- Shared-axis targets colored by overlap.  
- Triple targets highlighted with a **black border** and labels.  

---
