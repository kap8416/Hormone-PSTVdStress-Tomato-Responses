#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_mtr_orthologs.py
--------------------------------
Filter and summarize orthologs of a given set of MTRs from a Solanaceae TSV file.

Expected TSV input columns:
- a : gene ID in species A (e.g., Sopen*, SPI*, SLYcer*, PGSC*...)
- b : gene ID in tomato (Solyc*)
- species_a : short name for species A (TomatoPennelli, Tomatopimp, Potato, etc.)
- species_b : usually "tomato"
- OG : (optional) ortholog group ID
- Normalized_bit_score : (optional) score

Outputs:
- <outprefix>_table.tsv            -> complete table of ortholog pairs for the MTRs
- <outprefix>_summary.tsv          -> summary per MTR with ortholog list by species
- <outprefix>_presence_absence.tsv -> binary matrix MTR x species
- <outprefix>_heatmap.png          -> (optional, with --heatmap) presence/absence heatmap (matplotlib)
"""

import argparse
import pandas as pd
import re
from pathlib import Path

def base_id(x: str) -> str:
    """Remove transcript version, e.g., Solyc06g060490.2 -> Solyc06g060490."""
    return re.sub(r"\.\d+$", "", str(x).strip())

def load_orthologs(tsv_path: str) -> pd.DataFrame:
    df = pd.read_csv(tsv_path, sep="\t")
    required = {"a", "b", "species_a"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in TSV: {missing}")
    df["a_base"] = df["a"].astype(str).apply(base_id)
    df["b_base"] = df["b"].astype(str).apply(base_id)
    return df

def filter_by_mtrs(df: pd.DataFrame, mtrs: list[str]) -> pd.DataFrame:
    target_set = set(map(base_id, mtrs))
    out = df[df["b_base"].isin(target_set)].copy()
    out = out.sort_values(["b_base", "species_a", "a_base"]).reset_index(drop=True)
    return out

def make_summary(df_mtrs: pd.DataFrame) -> pd.DataFrame:
    def _join_group(g):
        items = [f"{row['a']} ({row['species_a']})" for _, row in g.iterrows()]
        return pd.Series({
            "MTR": g["b_base"].iloc[0],
            "Orthologs": ", ".join(items)
        })
    return df_mtrs.groupby("b_base", as_index=False).apply(_join_group).reset_index(drop=True)

def presence_absence(df_mtrs: pd.DataFrame) -> pd.DataFrame:
    mat = (
        df_mtrs.groupby(["b_base", "species_a"])
               .size().unstack(fill_value=0)
    )
    mat = (mat > 0).astype(int).reset_index().rename(columns={"b_base": "MTR"})
    return mat

def plot_heatmap(mat_df, out_png):
    import matplotlib.pyplot as plt
    import numpy as np

    df = mat_df.set_index("MTR")
    data = df.values
    fig, ax = plt.subplots(figsize=(max(6, 1.2*data.shape[1]), 0.8*data.shape[0]+2))
    im = ax.imshow(data, aspect="auto")

    # Manual colormap (greens), no seaborn
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap(["#f7fcf5", "#00441b"])
    im.set_cmap(cmap)
    im.set_clim(0,1)

    # Ticks/labels
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    ax.set_xticklabels(df.columns, rotation=45, ha="right")
    ax.set_yticklabels(df.index)

    # Grid and annotations
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            ax.text(j, i, f"{int(data[i, j])}", ha="center", va="center",
                    color="white" if data[i,j] > 0.5 else "black", fontsize=9)

    ax.set_title("Presence/Absence of MTR Orthologs in Solanaceae Species", fontsize=12)
    ax.set_xlabel("Species")
    ax.set_ylabel("MTR")
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tsv", required=True, help="Path to orthologs_Solanaceae.tsv")
    p.add_argument("--mtrs", nargs="+", required=True, help="List of Solyc IDs of MTRs (with or without version)")
    p.add_argument("--outprefix", default="MTRs_orthologs", help="Output file prefix")
    p.add_argument("--heatmap", action="store_true", help="Generate heatmap PNG (matplotlib)")
    args = p.parse_args()

    outprefix = Path(args.outprefix)
    outprefix.parent.mkdir(parents=True, exist_ok=True)

    df = load_orthologs(args.tsv)
    df_mtrs = filter_by_mtrs(df, args.mtrs)

    # Save full table
    table_path = f"{outprefix}_table.tsv"
    df_mtrs.drop(columns=["a_base","b_base"], errors="ignore").to_csv(table_path, sep="\t", index=False)

    # Consolidated summary
    summary_df = make_summary(df_mtrs)
    summary_path = f"{outprefix}_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)

    # Presence/absence matrix
    pa_df = presence_absence(df_mtrs)
    pa_path = f"{outprefix}_presence_absence.tsv"
    pa_df.to_csv(pa_path, sep="\t", index=False)

    # Optional heatmap
    if args.heatmap:
        heatmap_path = f"{outprefix}_heatmap.png"
        plot_heatmap(pa_df, heatmap_path)

    print("Saved:")
    print("  -", table_path)
    print("  -", summary_path)
    print("  -", pa_path)
    if args.heatmap:
        print("  -", heatmap_path)

if __name__ == "__main__":
    main()
