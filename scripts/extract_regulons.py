# Load and filter the network to extract regulons for the requested MTRs
import pandas as pd
from pathlib import Path
from caas_jupyter_tools import display_dataframe_to_user

# Path to the uploaded network
path = Path("/mnt/data/net_corto_prom.txt")

# Read the tab-separated file (source, target, weight)
df = pd.read_csv(path, sep=r"\s+", header=None, names=["MTR","Target","Weight"], engine="python")

# MTRs of interest
mtrs = {
    "NAC – Solyc04g009440.2": "Solyc04g009440.2",
    "AP2/EREBP TF1 – Solyc02g093130.1": "Solyc02g093130.1",
    "ERF_C_5 – Solyc02g077370.1": "Solyc02g077370.1",
    "MYC2 – Solyc08g076930.1": "Solyc08g076930.1",
    "bZIP (ABI5-like) – Solyc06g060490.2": "Solyc06g060490.2"
}

# Filter regulons
regulons = {label: df[df["MTR"] == gene].copy() for label, gene in mtrs.items()}

# Save each regulon to CSV and show a compact preview
downloads = {}
for label, subdf in regulons.items():
    # Sort by absolute weight descending for interpretability
    subdf_sorted = subdf.sort_values(by="Weight", key=lambda s: s.abs(), ascending=False).reset_index(drop=True)
    # Save
    safe_name = label.split("–")[-1].strip().replace(".", "_")
    out_path = Path(f"/mnt/data/regulon_{safe_name}.csv")
    subdf_sorted.to_csv(out_path, index=False)
    downloads[label] = str(out_path)

    # Show to user (preview top 25 for readability)
    preview = subdf_sorted.head(25)
    display_dataframe_to_user(name=f"Regulon preview: {label}", dataframe=preview)

# Also create a summary table with counts and sign breakdown
summary_rows = []
for label, subdf in regulons.items():
    n = len(subdf)
    n_pos = (subdf["Weight"] > 0).sum()
    n_neg = (subdf["Weight"] < 0).sum()
    summary_rows.append({
        "MTR": label,
        "Total_targets": n,
        "Positive_edges": n_pos,
        "Negative_edges": n_neg
    })

summary_df = pd.DataFrame(summary_rows).sort_values("Total_targets", ascending=False).reset_index(drop=True)

# Save and display summary
summary_path = Path("/mnt/data/regulon_summary_five_MTRs.csv")
summary_df.to_csv(summary_path, index=False)
display_dataframe_to_user(name="Resumen de regulones (5 MTRs)", dataframe=summary_df)

print("Archivos generados:")
for label, p in downloads.items():
    print(f"- {label}: {p}")
print(f"- Resumen: {summary_path}")
