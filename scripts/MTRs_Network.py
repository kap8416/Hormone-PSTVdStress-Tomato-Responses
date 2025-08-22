# MTRs_Network.py
# Script to generate radial network visualization of five MTRs and their regulons in tomato under PSTVd infection.
# Author: Katia Aviña Padilla (2025)

# -*- coding: utf-8 -*-
"""
MTRs_Network.py
--------------------------------
Radial visualization of five MTRs with:
1) Hubs (MTRs): FILL = hormone; BORDER = expression (ERF_C_5 induced=red, others repressed=blue).
2) Edges:
   - If the target is SHARED → color according to hormonal combo (E–J brown, E–A pink, A–J purple, E–J–A turquoise).
   - If the target is UNIQUE → color of the originating HUB (same as its hormonal fill).
3) Shared targets: the NODE is colored with the same combo color (same as the edge).
4) Labels for the three E–J–A genes if present: Solyc06g008330, Solyc08g075450, Solyc10g044780.

Input: net_corto_prom.txt (two columns: source, target; auto-detected separator)
Outputs:  MTRs_Network.png / .svg / .pdf
Requirements: Python 3.8+, networkx, matplotlib, pandas
"""

import re, math, os
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from collections import defaultdict

# -------------------- I/O --------------------
IN_FILE   = "net_corto_prom.txt"
OUT_PNG   = "MTRs_Network.png"
OUT_SVG   = "MTRs_Network.svg"
OUT_PDF   = "MTRs_Network.pdf"

# -------------------- Palettes --------------------
# Hub fill colors by HORMONE
HORM_FILL = {
    "ABA":   "#f2c744",
    "A":     "#6cc36c",   # Auxin
    "E_act": "#f28c28",   # Ethylene activator
    "E_rep": "#ef6a6a",   # Ethylene repressor
    "J":     "#4d9de0",   # Jasmonate
}

# Edge and shared target colors by HORMONAL COMBO
COMBO_COLORS = {
    "EJ":  "#795548",   # brown
    "EA":  "#e91e63",   # pink
    "AJ":  "#7e57c2",   # purple
    "AEJ": "#40E0D0",   # turquoise (changed from red)
}

# -------------------- Definition of MTRs --------------------
MTRS = {
    "Solyc08g076930.1": {"abbr":"MYC2 (bHLH)",       "horm_logic":"J", "fill": HORM_FILL["J"]},
    "Solyc02g077370.1": {"abbr":"ERF_C_5",           "horm_logic":"E", "fill": HORM_FILL["E_act"]},  # ET activator
    "Solyc02g093130.1": {"abbr":"AP2/EREBP TF1",     "horm_logic":"E", "fill": HORM_FILL["E_rep"]},  # ET repressor
    "Solyc06g060490.2": {"abbr":"ABI5-like (bZIP)",  "horm_logic":"B", "fill": HORM_FILL["ABA"]},    # ABA
    "Solyc04g009440.2": {"abbr":"NAC",               "horm_logic":"A", "fill": HORM_FILL["A"]},
}

# Borders by expression: ERF_C_5 induced=red, others repressed=blue
INDUCED  = {"Solyc02g077370.1"}  
RINGCOL  = {k: ("#d32f2f" if k in INDUCED else "#1e88e5") for k in MTRS}

# -------------------- Utilities --------------------
def base_id(x: str) -> str:
    """Remove isoform suffix (.1, .2, …)."""
    return re.sub(r"\.\d+$", "", str(x).strip())

MTRS_BASE = {base_id(k): k for k in MTRS}

def load_edges(path: str) -> pd.DataFrame:
    """Robust load; detects source/target columns and normalizes to base-id."""
    df = pd.read_csv(path, sep=None, engine="python")
    cols = [str(c).strip() for c in df.columns]
    df.columns = cols
    out = df[[cols[0], cols[1]]].copy()
    out.columns = ["source", "target"]
    out["source_base"] = out["source"].apply(base_id)
    out["target_base"] = out["target"].apply(base_id)
    return out

# -------------------- Load data --------------------
edges = load_edges(IN_FILE)
edges = edges[edges["source_base"].isin(MTRS_BASE)].copy()
if edges.empty:
    raise SystemExit("No edges with source among the 5 MTRs found in input file.")

# -------------------- Graph --------------------
G = nx.DiGraph()
for b, full in MTRS_BASE.items():
    meta = MTRS[full]
    G.add_node(b, kind="MTR", label="", fill=meta["fill"],  # <- NO label for hubs
               ring=RINGCOL[full], logic=meta["horm_logic"], full_id=full)

for _, row in edges.iterrows():
    mtr_b = row["source_base"]
    tgt_b = row["target_base"]
    if tgt_b not in G:
        G.add_node(tgt_b, kind="target")
    G.add_edge(mtr_b, tgt_b)

# -------------------- Shared & Unique --------------------
targets = [n for n,d in G.nodes(data=True) if d["kind"]=="target"]
deg_by_mtr = {t: sum(1 for m in MTRS_BASE if G.has_edge(m,t)) for t in targets}
shared_targets = [t for t,c in deg_by_mtr.items() if c >= 2]
unique_targets = [t for t,c in deg_by_mtr.items() if c == 1]

axes = defaultdict(set)
for t in targets:
    for m in MTRS_BASE:
        if G.has_edge(m, t):
            ax = G.nodes[m]["logic"]
            if ax in {"A", "E", "J"}:
                axes[t].add(ax)

def combo_code(s):
    code = "".join(sorted(s))
    if code == "AE": code = "EA"
    if code in {"AEJ", "EJA"}: code = "AEJ"
    return code

combo_per_target = {t: combo_code(ax) for t, ax in axes.items()}

# -------------------- Layout --------------------
n_mtrs = len(MTRS_BASE)
order = list(MTRS_BASE.keys())
angles = {m: (i * 2*math.pi / n_mtrs, (i+1)*2*math.pi / n_mtrs) for i, m in enumerate(order)}
pos = {}
R_mtr, R_outer = 1.6, 3.6

for m in order:
    a0, a1 = angles[m]
    theta = 0.5*(a0 + a1)
    pos[m] = (R_mtr*math.cos(theta), R_mtr*math.sin(theta))

for m in order:
    a0, a1 = angles[m]
    m_targets = [t for t in targets if G.has_edge(m, t)]
    if not m_targets: continue
    uniq = sorted(set(m_targets))
    n = len(uniq)
    for j, t in enumerate(uniq):
        frac  = j / max(1, n-1)
        theta = a0 + frac*(a1 - a0)
        r     = R_outer + 0.5*math.sin(frac*math.pi)
        pos[t]= (r*math.cos(theta), r*math.sin(theta))

# -------------------- Edge & Node colors --------------------
edge_colors = []
for u, v in G.edges():
    combo = combo_per_target.get(v, "")
    if (v in shared_targets) and (combo in COMBO_COLORS):
        edge_colors.append(COMBO_COLORS[combo])
    else:
        full_id = G.nodes[u]["full_id"]
        edge_colors.append(MTRS[full_id]["fill"])

TARGET_GREY = "#C8CCD0"
shared_node_colors = [COMBO_COLORS.get(combo_per_target.get(t, ""), "#000000") for t in shared_targets]

# -------------------- Plot --------------------
plt.figure(figsize=(11, 11))

# Directed edges with arrows
nx.draw_networkx_edges(G, pos, arrows=True, arrowsize=12, width=0.9, alpha=0.85, edge_color=edge_colors)

# Unique targets
nx.draw_networkx_nodes(G, pos, nodelist=unique_targets, node_color=TARGET_GREY, node_size=22, edgecolors="none")

# Shared targets (default size)
nx.draw_networkx_nodes(G, pos, nodelist=[t for t in shared_targets if base_id(t) not in {"Solyc06g008330","Solyc08g075450","Solyc10g044780"}],
                       node_color=[COMBO_COLORS.get(combo_per_target.get(t, ""), "#000000") for t in shared_targets if base_id(t) not in {"Solyc06g008330","Solyc08g075450","Solyc10g044780"}],
                       node_size=30, edgecolors="black", linewidths=1.2)

# Triple AEJ targets in turquoise and bigger
special_nodes = {"Solyc06g008330","Solyc08g075450","Solyc10g044780"}
nx.draw_networkx_nodes(G, pos, nodelist=[t for t in shared_targets if base_id(t) in special_nodes],
                       node_color="#40E0D0", node_size=80, edgecolors="black", linewidths=1.5)

# Hubs (no label)
nx.draw_networkx_nodes(G, pos, nodelist=order,
                       node_color=[G.nodes[m]["fill"] for m in order],
                       node_size=1400, edgecolors=[G.nodes[m]["ring"] for m in order], linewidths=3.0)

# Labels only for special AEJ nodes
label_dict = {t: base_id(t) for t in shared_targets if base_id(t) in special_nodes}
nx.draw_networkx_labels(G, pos, labels=label_dict, font_size=8, font_weight="bold")

plt.axis('off')
plt.tight_layout()
plt.savefig(OUT_PNG, dpi=450, bbox_inches="tight")
plt.savefig(OUT_SVG, bbox_inches="tight")
plt.savefig(OUT_PDF, bbox_inches="tight")
print("Saved:", OUT_PNG, OUT_SVG, OUT_PDF)
