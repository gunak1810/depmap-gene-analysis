# ============================================================
# 🧬 Post-Analysis: Enrichment, Overlap & Heatmap (Fixed Import)
# Author: Guna Kulothungan | 2025
# ============================================================

import os, glob, pandas as pd, numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(context="talk", style="whitegrid")

# ------------------------------------------------------------
# Step 1 — Locate per-cancer result files
# ------------------------------------------------------------
RESULT_DIR = "Cancer_Specific_Results"
OUT_DIR = "Enrichment_and_Overlap_Results"
os.makedirs(OUT_DIR, exist_ok=True)

pattern = os.path.join(RESULT_DIR, "CRISPR_TME_*_Results.csv")
files = sorted(glob.glob(pattern))
if not files:
    raise FileNotFoundError(f"No result CSVs found in {RESULT_DIR}. Run the main analysis first.")

print(f"Found {len(files)} cancer-specific result files.")

TOP_N = 100   # number of top genes to include in overlap/enrichment

# ------------------------------------------------------------
# Step 2 — Load top gene lists for each cancer
# ------------------------------------------------------------
cancer_top_pos, cancer_top_neg, cancer_list = {}, {}, []

for f in files:
    df = pd.read_csv(f)
    df["Gene_clean"] = df["Gene"].astype(str).apply(lambda s: s.split(" (")[0].strip())
    cancer = os.path.basename(f).replace("CRISPR_TME_", "").replace("_Results.csv", "")
    cancer_list.append(cancer)
    pos = df.sort_values("Spearman_r", ascending=False).head(TOP_N)["Gene_clean"].tolist()
    neg = df.sort_values("Spearman_r").head(TOP_N)["Gene_clean"].tolist()
    cancer_top_pos[cancer] = pos
    cancer_top_neg[cancer] = neg

print("Loaded gene lists for:", ", ".join(cancer_list))

# ------------------------------------------------------------
# Step 3 — GO/KEGG enrichment using g:Profiler
# ------------------------------------------------------------
print("\nRunning GO/KEGG enrichment (requires internet)...")

# Fixed import for Python 3.12+/Windows
try:
    from gprofiler_official import GProfiler
except ImportError:
    print("Installing gprofiler-official...")
    os.system("pip install gprofiler-official")
    from gprofiler_official import GProfiler

gp = GProfiler(return_dataframe=True)

for cancer in cancer_list:
    for label, gene_list in {"Positive": cancer_top_pos[cancer],
                             "Negative": cancer_top_neg[cancer]}.items():
        if not gene_list:
            continue
        try:
            enrich = gp.profile(
                organism="hsapiens",
                query=gene_list,
                sources=["GO:BP", "KEGG"]
            )
            out_csv = os.path.join(OUT_DIR, f"Enrichment_{label}_{cancer}.csv")
            enrich.to_csv(out_csv, index=False)
            # Plot top 10 terms
            top = enrich.head(10).copy()
            if not top.empty:
                top["-log10(p)"] = -np.log10(top["p_value"] + 1e-300)
                plt.figure(figsize=(8, 4))
                sns.barplot(x="-log10(p)", y="name", data=top, color="skyblue")
                plt.title(f"{cancer} — Top {label} Enriched Terms")
                plt.tight_layout()
                plt.savefig(os.path.join(OUT_DIR, f"Enrich_{label}_{cancer}.png"), dpi=300)
                plt.close()
        except Exception as e:
            print(f"⚠️ Enrichment failed for {cancer} ({label}):", e)

# ------------------------------------------------------------
# Step 4 — Jaccard similarity (overlap of top positive genes)
# ------------------------------------------------------------
print("\nCalculating Jaccard similarity between cancers...")

def jaccard(a, b):
    A, B = set(a), set(b)
    if not A and not B:
        return 0.0
    return len(A & B) / len(A | B)

n = len(cancer_list)
jmat = np.zeros((n, n))
for i, a in enumerate(cancer_list):
    for j, b in enumerate(cancer_list):
        jmat[i, j] = jaccard(cancer_top_pos[a], cancer_top_pos[b])

jdf = pd.DataFrame(jmat, index=cancer_list, columns=cancer_list)
jdf.to_csv(os.path.join(OUT_DIR, "Jaccard_Positive_Matrix.csv"))

plt.figure(figsize=(10, 8))
sns.heatmap(jdf, cmap="viridis", linewidths=0.3)
plt.title(f"Jaccard Similarity (Top {TOP_N} Positive Genes)")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "Jaccard_Positive_Heatmap.png"), dpi=300)
plt.close()

# ------------------------------------------------------------
# Step 5 — Gene frequency and presence/absence heatmap
# ------------------------------------------------------------
from collections import Counter
pos_counter = Counter()
for c in cancer_list:
    pos_counter.update(cancer_top_pos[c])

pos_df = pd.DataFrame(pos_counter.items(), columns=["Gene", "Count"]).sort_values("Count", ascending=False)
pos_df.to_csv(os.path.join(OUT_DIR, "Global_Top_Positive_Genes.csv"), index=False)

# Plot top 20 most recurrent positive genes
plt.figure(figsize=(8, 6))
sns.barplot(x="Count", y="Gene", data=pos_df.head(20), palette="viridis")
plt.title(f"Top {TOP_N} Positive Genes — Frequency Across Cancers")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "Top20_Positive_Frequency.png"), dpi=300)
plt.close()

# Presence/absence matrix
M = 60
top_genes = pos_df["Gene"].head(M).tolist()
mat = pd.DataFrame(0, index=top_genes, columns=cancer_list)
for c in cancer_list:
    for g in cancer_top_pos[c]:
        if g in mat.index:
            mat.loc[g, c] = 1
mat.to_csv(os.path.join(OUT_DIR, "Presence_Absence_Positive.csv"))

plt.figure(figsize=(12, 8))
sns.heatmap(mat, cmap="Blues", cbar_kws={"label": "Presence (1=Top Gene)"}, linewidths=0.3)
plt.title(f"Presence of Frequent Positive Genes Across Cancers (Top {M})")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "Presence_Positive_Heatmap.png"), dpi=300)
plt.close()

# ------------------------------------------------------------
# Step 6 — Overlap count matrix (shared positive genes)
# ------------------------------------------------------------
overlap = np.zeros((n, n), dtype=int)
for i, a in enumerate(cancer_list):
    for j, b in enumerate(cancer_list):
        overlap[i, j] = len(set(cancer_top_pos[a]) & set(cancer_top_pos[b]))

overlap_df = pd.DataFrame(overlap, index=cancer_list, columns=cancer_list)
overlap_df.to_csv(os.path.join(OUT_DIR, "Overlap_Count_Positive.csv"))

print("\n✅ Post-analysis complete!")
print(f"All files saved to → {OUT_DIR}")
