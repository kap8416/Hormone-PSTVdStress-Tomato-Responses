# Five_MTRs_radial_EJA_combos_v5.py
# Script to generate radial network visualization of five MTRs and their regulons in tomato under PSTVd infection.
# Author: Katia Aviña Padilla (2025)

import re, math
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as mpatches
from collections import defaultdict

IN_FILE = "net_corto_prom.txt"
OUT_PNG = "Five_MTRs_radial_EJA_combos_v5.png"
OUT_SVG = "Five_MTRs_radial_EJA_combos_v5.svg"
OUT_PDF = "Five_MTRs_radial_EJA_combos_v5.pdf"

MTRS = {
    "Solyc06g060490.2": {"hormone": "ABA",      "role": "Activator", "color": "#f2c744", "abbr": "ABI5-like (bZIP)"},
    "Solyc04g009440.2": {"hormone": "Auxin",    "role": "Activator", "color": "#6cc36c", "abbr": "NAC"},
    "Solyc02g093130.1": {"hormone": "Ethylene", "role": "Repressor", "color": "#ef6a6a", "abbr": "AP2/EREBP TF1"},
    "Solyc02g077370.1": {"hormone": "Ethylene", "role": "Activator", "color": "#f28c28", "abbr": "ERF_C_5"},
    "Solyc08g076930.1": {"hormone": "Jasmonate","role": "Activator", "color": "#4d9de0", "abbr": "MYC2 (bHLH)"},
}

def base_id(x: str) -> str:
    return re.sub(r"\.\d+$", "", str(x).strip())

MTRS_BASE = {base_id(k): k for k in MTRS.keys()}

def load_edges(path):
    df = pd.read_csv(path, sep=None, engine="python")
    cols = [str(c).strip().lower() for c in df.columns]
    df.columns = cols
    src_candidates = [c for c in cols if any(k in c for k in ["src","source","from","tf","regulator","head","u"])]
    tgt_candidates = [c for c in cols if any(k in c for k in ["tgt","target","to","gene","tail","v"]) and c not in src_candidates]
    if not src_candidates or not tgt_candidates:
        src_col, tgt_col = cols[0], cols[1]
    else:
        src_col, tgt_col = src_candidates[0], tgt_candidates[0]
    out = df[[src_col, tgt_col]].copy()
    out.columns = ["source", "target"]
    out["source_base"] = out["source"].apply(base_id)
    out["target_base"] = out["target"].apply(base_id)
    return out

edges = load_edges(IN_FILE)

# Regulons
regulons = {k: [] for k in MTRS.keys()}
for _, row in edges.iterrows():
    src_base = row["source_base"]
    tgt_full = row["target"]
    if src_base in MTRS_BASE:
        mtr_key = MTRS_BASE[src_base]
        regulons[mtr_key].append(tgt_full)

# Axis sets and combos
target_axes = defaultdict(set)
for mtr, targets in regulons.items():
    axis = MTRS[mtr]["hormone"]
    for t in set(targets):
        if axis == "Ethylene":  target_axes[t].add("E")
        elif axis == "Jasmonate": target_axes[t].add("J")
        elif axis == "Auxin":   target_axes[t].add("A")

def combo_of(axset):
    s = "".join(sorted(axset))
    if s == "AE": s = "EA"
    if s in {"EJA","AEJ"}: s = "AEJ"
    return s

target_combo = {t: combo_of(ax) for t, ax in target_axes.items()}

combo_colors = {"EJ":"#1f77b4", "EA":"#2ca02c", "AJ":"#ff7f0e", "AEJ":"#d62728"}
default_target_color = "#C8CCD0"

# Build graph
G = nx.DiGraph()
for mtr, meta in MTRS.items():
    G.add_node(mtr, kind="MTR", color=meta["color"], size=1400, label=meta["abbr"])

for mtr, targets in regulons.items():
    for t in set(targets):
        if t not in G:
            combo = target_combo.get(t, "")
            col = combo_colors.get(combo, default_target_color)
            is_combo = combo in combo_colors
            G.add_node(t, kind="target", color=col, combo=combo,
                       size=180 if is_combo else 120, is_combo=is_combo)
        style = "dashed" if MTRS[mtr]["role"] == "Repressor" else "solid"
        ecol  = "#d33f49" if MTRS[mtr]["role"] == "Repressor" else "#2a9d8f"
        G.add_edge(mtr, t, style=style, color=ecol, width=1.0)

# Layout
n_mtrs = len(MTRS)
angles = {mtr: (i * 2*math.pi / n_mtrs, (i+1)*2*math.pi / n_mtrs) for i, mtr in enumerate(MTRS.keys())}
pos = {}
R_mtr = 1.6
for mtr in MTRS.keys():
    a0, a1 = angles[mtr]
    theta = (a0 + a1) / 2.0
    pos[mtr] = (R_mtr*math.cos(theta), R_mtr*math.sin(theta))

R1 = 3.4
for mtr, targets in regulons.items():
    a0, a1 = angles[mtr]
    uniq_targets = sorted(list(set(targets)))
    n = len(uniq_targets)
    if n == 0: continue
    for j, t in enumerate(uniq_targets):
        frac = j / max(1, n-1)
        theta = a0 + frac*(a1 - a0)
        r = R1 + 0.6*math.sin(frac*math.pi)
        pos[t] = (r*math.cos(theta), r*math.sin(theta))

plt.figure(figsize=(11, 11))

# Draw edges
for style in ["solid", "dashed"]:
    edges_list = [(u, v) for u, v, d in G.edges(data=True) if d["style"] == style]
    colors = [G.edges[e]["color"] for e in edges_list]
    widths = [G.edges[e]["width"] for e in edges_list]
    nx.draw_networkx_edges(G, pos, edgelist=edges_list, edge_color=colors,
                           arrows=True, arrowstyle='-|>', arrowsize=10,
                           width=widths, style=style, alpha=0.65)

# Base targets (grey)
base_targets = [n for n, d in G.nodes(data=True) if d["kind"]=="target" and not d.get("is_combo", False)]
if base_targets:
    nx.draw_networkx_nodes(G, pos, nodelist=base_targets,
                           node_color=[G.nodes[n]["color"] for n in base_targets],
                           node_size=[G.nodes[n]["size"] for n in base_targets],
                           alpha=0.95, linewidths=0)

# MTRs
mtr_nodes = [n for n, d in G.nodes(data=True) if d["kind"]=="MTR"]
nx.draw_networkx_nodes(G, pos, nodelist=mtr_nodes,
                       node_color=[G.nodes[n]["color"] for n in mtr_nodes],
                       node_size=[G.nodes[n]["size"] for n in mtr_nodes],
                       edgecolors="#2F3B45", linewidths=1.5)

# Combo nodes without border first
combo_nodes = [n for n, d in G.nodes(data=True) if d.get("is_combo", False) and G.nodes[n]["combo"] != "AEJ"]
if combo_nodes:
    nx.draw_networkx_nodes(G, pos, nodelist=combo_nodes,
                           node_color=[G.nodes[n]["color"] for n in combo_nodes],
                           node_size=[G.nodes[n]["size"] for n in combo_nodes],
                           alpha=1.0, linewidths=0)

# Triple nodes (AEJ) with black border
triple_nodes = [n for n, d in G.nodes(data=True) if d.get("is_combo", False) and G.nodes[n]["combo"] == "AEJ"]
if triple_nodes:
    nx.draw_networkx_nodes(G, pos, nodelist=triple_nodes,
                           node_color=[G.nodes[n]["color"] for n in triple_nodes],
                           node_size=[max(200, G.nodes[n]["size"]) for n in triple_nodes],
                           alpha=1.0, edgecolors="black", linewidths=1.6)

# Labels for the three specific triples (if present)
label_ids = { "Solyc06g008330", "Solyc08g075450", "Solyc10g044780",
              "Solyc06g008330.1", "Solyc06g008330.2",
              "Solyc08g075450.1", "Solyc08g075450.2",
              "Solyc10g044780.1", "Solyc10g044780.2" }
present = {n: n for n in triple_nodes if n in label_ids or base_id(n) in {base_id(x) for x in label_ids}}
if present:
    nx.draw_networkx_labels(G, pos, labels=present, font_size=8, font_weight='bold')

# Legends
hormone_legend = [
    mpatches.Patch(color=MTRS["Solyc06g060490.2"]["color"], label="ABA"),
    mpatches.Patch(color=MTRS["Solyc04g009440.2"]["color"], label="Auxin"),
    mpatches.Patch(color=MTRS["Solyc02g093130.1"]["color"], label="Ethylene (Repressor)"),
    mpatches.Patch(color=MTRS["Solyc02g077370.1"]["color"], label="Ethylene (Activator)"),
    mpatches.Patch(color=MTRS["Solyc08g076930.1"]["color"], label="Jasmonate"),
    mpatches.Patch(color="#C8CCD0", label="Targets (others)"),
]
combo_legend = [
    mpatches.Patch(color="#1f77b4", label="E–J shared"),
    mpatches.Patch(color="#2ca02c", label="E–A shared"),
    mpatches.Patch(color="#ff7f0e", label="A–J shared"),
    mpatches.Patch(color="#d62728", label="E–J–A shared (triple)"),
]
leg1 = plt.legend(handles=hormone_legend, loc="lower left", bbox_to_anchor=(0.0, -0.02),
                  ncol=3, frameon=False, fontsize=9)
leg2 = plt.legend(handles=combo_legend, loc="lower right", bbox_to_anchor=(1.0, -0.02),
                  ncol=2, frameon=False, fontsize=9)
plt.gca().add_artist(leg1)

plt.title("Five MTRs and Their Regulons in PSTVd Tomato\nShared-axis targets highlighted (E–J, E–A, A–J, and E–J–A)", fontsize=12, pad=12)
plt.axis('off')
plt.tight_layout()

plt.savefig(OUT_PNG, dpi=450, bbox_inches="tight")
plt.savefig(OUT_SVG, bbox_inches="tight")
plt.savefig(OUT_PDF, bbox_inches="tight")
plt.close()

print("Generated files:", OUT_PNG, OUT_SVG, OUT_PDF)
