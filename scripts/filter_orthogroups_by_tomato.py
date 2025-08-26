#!/usr/bin/env python3
"""
Filter orthogroups by Tomato IDs (ignoring version after the first dot).

Reads:
  - Tomato IDs file (CSV/TSV/any delimiter). All cells are scanned for tokens.
  - Orthogroups TSV (e.g., "tl_project.orthogroups.tsv").

Matches:
  - A row is kept if ANY token in any species column has a base ID (before the
    first '.') matching one of the Tomato base IDs from the Tomato file.

Writes:
  - A filtered TSV preserving the full identifiers exactly as they appear in
    the orthogroups file.

Example:
  python filter_orthogroups_by_tomato.py \
      --tomato TomatoABA.csv \
      --orthogroups tl_project.orthogroups.tsv \
      --out filtered_orthogroups.tsv

Requirements:
  - Python 3.8+
  - pandas
    pip install pandas
"""

import argparse
import re
import sys
from typing import List, Set

import pandas as pd


def load_id_bases_from_csv(path_csv: str) -> Set[str]:
    """Load IDs from a CSV/TSV (delimiter auto-detected) and return base IDs before first dot."""
    try:
        df = pd.read_csv(path_csv, sep=None, engine="python", dtype=str)
    except Exception:
        # Fallback to comma if autodetect fails
        df = pd.read_csv(path_csv, dtype=str)

    # Flatten all values to strings and split on commas/whitespace/semicolons
    raw_vals: List[str] = []
    for col in df.columns:
        col_vals = df[col].dropna().astype(str).tolist()
        raw_vals.extend(col_vals)

    tokens: List[str] = []
    for s in raw_vals:
        tokens.extend([t for t in re.split(r"[,\s;]+", s.strip()) if t])

    # Keep token base before the first dot
    bases: Set[str] = set()
    for tok in tokens:
        base = tok.split(".", 1)[0]
        if base:
            bases.add(base)
    return bases


def row_matches_bases(row: pd.Series, base_set: Set[str]) -> bool:
    """Return True if any gene token (base part) in the row matches base_set."""
    for val in row:
        if pd.isna(val):
            continue
        s = str(val)
        # Tokens can be comma/space/semicolon separated
        tokens = [t for t in re.split(r"[,\s;]+", s.strip()) if t]
        for t in tokens:
            base = t.split(".", 1)[0]
            if base in base_set:
                return True
    return False


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="Filter orthogroups by Tomato IDs (ignoring version after the first dot)."
    )
    p.add_argument("--tomato", required=True, help="Path to Tomato IDs file (e.g., TomatoABA.csv)")
    p.add_argument("--orthogroups", required=True, help="Path to orthogroups TSV (e.g., tl_project.orthogroups.tsv)")
    p.add_argument("--out", default="filtered_orthogroups.tsv", help="Output TSV path (default: filtered_orthogroups.tsv)")
    args = p.parse_args(argv)

    # Load Tomato base IDs
    tomato_bases = load_id_bases_from_csv(args.tomato)
    if not tomato_bases:
        print("[ERROR] No Tomato base IDs parsed from input. Check the file format.", file=sys.stderr)
        return 2

    # Load orthogroups table
    try:
        df_og = pd.read_csv(args.orthogroups, sep="\t", dtype=str)
    except Exception as e:
        print(f"[ERROR] Could not read orthogroups TSV: {e}", file=sys.stderr)
        return 3

    if df_og.empty:
        print("[ERROR] Orthogroups table is empty.", file=sys.stderr)
        return 4

    # Detect gene columns (skip 'Orthogroup' column if present at index 0)
    cols_for_match = df_og.columns.tolist()
    if cols_for_match and cols_for_match[0].lower().startswith("orthogroup"):
        cols_for_match = cols_for_match[1:]
    if not cols_for_match:
        # Fallback: match across all columns
        cols_for_match = df_og.columns.tolist()

    # Build mask & filter
    match_mask = df_og[cols_for_match].apply(lambda row: row_matches_bases(row, tomato_bases), axis=1)
    df_filtered = df_og[match_mask].copy()

    # Save preserving original cell values
    df_filtered.to_csv(args.out, sep="\t", index=False)

    # Summary
    print(f"[OK] Wrote: {args.out}")
    print(f"  - Total orthogroups:      {len(df_og)}")
    print(f"  - Matched orthogroups:    {len(df_filtered)}")
    print(f"  - Unique Tomato base IDs: {len(tomato_bases)}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
