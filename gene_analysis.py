# ============================================================
# 🔬 Multi-Cancer CRISPR × Immune Correlation Pipeline
# Author: Guna Kulothungan | Conference Project | 2025
# ============================================================

import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns
import os

# ------------------------------------------------------------
# Step 1 — Load Core Files
# ------------------------------------------------------------
print("Loading data...")
model = pd.read_csv("Model.csv")
crispr = pd.read_csv("CRISPRGeneEffect.csv", index_col=0)
omics_expr = pd.read_csv("OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")
print("Files loaded successfully!")

# Fix orientation (rows = genes, columns = models)
if crispr.index[0].startswith("ACH-"):
    print("Transposing CRISPR matrix...")
    crispr = crispr.transpose()
print(f"CRISPR shape after orientation: {crispr.shape}")

# Clean omics
omics_expr = omics_expr.set_index("ModelID")
omics_expr = omics_expr[~omics_expr.index.duplicated(keep="first")]
print(f"OMICS shape after cleaning: {omics_expr.shape}")

# ------------------------------------------------------------
# Step 2 — Define Immune Proxy Genes
# ------------------------------------------------------------
immune_markers = ["CD274","CXCL9","CXCL10","STAT1","IRF1","HLA-A","HLA-B","B2M","TAP1","TGFB1"]
immune_genes = [col for col in omics_expr.columns if any(marker in col for marker in immune_markers)]
print(f"Matched {len(immune_genes)} immune genes from omics expression.")

# Output directory
os.makedirs("Cancer_Specific_Results", exist_ok=True)

# ------------------------------------------------------------
# Step 3 — Loop Through Each Cancer Type
# ------------------------------------------------------------
unique_cancers = sorted(model["OncotreePrimaryDisease"].dropna().unique())
print(f"\nFound {len(unique_cancers)} unique cancer types.")
print(unique_cancers[:10], "...")  # preview first 10

summary_rows = []

for cancer in unique_cancers:
    print(f"\n🎯 Analyzing {cancer} cancer...")

    # Identify models for this cancer
    cancer_models = model[model["OncotreePrimaryDisease"].str.lower() == cancer.lower()]["ModelID"].unique().tolist()
    common_models = sorted(list(set(crispr.columns).intersection(set(omics_expr.index)).intersection(set(cancer_models))))

    if len(common_models) < 25:
        print(f"⚠️ Skipping {cancer}: only {len(common_models)} models.")
        continue

    crispr_sub = crispr[common_models]
    omics_sub = omics_expr.loc[common_models]

    immune_proxy = omics_sub[immune_genes].mean(axis=1)

    results = []
    for gene in crispr_sub.index:
        x = crispr_sub.loc[gene].astype(float).values
        y = immune_proxy.astype(float).values
        mask = ~np.isnan(x) & ~np.isnan(y)
        if mask.sum() < 5:
            continue
        rho, p = spearmanr(x[mask], y[mask])
        results.append([gene, rho, p])

    if not results:
        print(f"⚠️ No valid results for {cancer}")
        continue

    results_df = pd.DataFrame(results, columns=["Gene", "Spearman_r", "P"])
    results_df["FDR"] = multipletests(results_df["P"], method="fdr_bh")[1]
    results_df["-log10P"] = -np.log10(results_df["P"])
    results_df = results_df.sort_values("FDR")

    # Save results
    out_file = f"Cancer_Specific_Results/CRISPR_TME_{cancer.replace('/', '_')}_Results.csv"
    results_df.to_csv(out_file, index=False)
    print(f"✅ Saved {len(results_df)} results for {cancer}")

    # Volcano plot
    plt.figure(figsize=(7,5))
    sns.scatterplot(data=results_df, x="Spearman_r", y="-log10P", s=10, alpha=0.7)
    plt.axvline(0, color='grey', linestyle='--', linewidth=0.8)
    plt.axhline(-np.log10(0.05), color='red', linestyle='--', linewidth=0.8, label="p=0.05")
    plt.title(f"{cancer} Cancer — CRISPR vs Immune Proxy")
    plt.xlabel("Spearman Correlation (ρ)")
    plt.ylabel("-log10(P-value)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"Cancer_Specific_Results/Volcano_{cancer.replace('/', '_')}.png", dpi=300)
    plt.close()

    # Summary
    top_pos = results_df.sort_values("Spearman_r", ascending=False).head(3)
    top_neg = results_df.sort_values("Spearman_r").head(3)
    summary_rows.append({
        "Cancer": cancer,
        "N_Models": len(common_models),
        "Top_Positive": "; ".join(top_pos["Gene"].tolist()),
        "Top_Negative": "; ".join(top_neg["Gene"].tolist())
    })

# ------------------------------------------------------------
# Step 4 — Summary Table of All Cancers
# ------------------------------------------------------------
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv("Cancer_Specific_Results/Cancerwise_Summary.csv", index=False)
print("\n✅ Multi-cancer analysis complete!")
print(f"Results saved in folder: Cancer_Specific_Results")
print(summary_df.head(10))
