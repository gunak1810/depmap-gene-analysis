import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ============================
# Load Your Files
# ============================

top_genes = pd.read_csv(
    r"C:\Users\subha\OneDrive\Desktop\Conference\Dep Map Data\Enrichment_and_Overlap_Results\Global_Top_Positive_Genes.csv"
)

compendium = pd.read_csv(
    r"C:\Users\subha\OneDrive\Desktop\Conference\Dep Map Data\Enrichment_and_Overlap_Results\Compendium_Cancer_Genes.tsv",
    sep="\t"
)

unfiltered = pd.read_csv(
    r"C:\Users\subha\OneDrive\Desktop\Conference\Dep Map Data\Enrichment_and_Overlap_Results\Unfiltered_drivers.tsv",
    sep="\t"
)

print("Loaded top genes:", top_genes.shape)
print("Compendium columns:", compendium.columns.tolist())
print("Unfiltered columns:", unfiltered.columns.tolist())

# ============================
# AUTO-DETECT gene column
# ============================

possible_cols = ["SYMBOL", "symbol", "Gene", "gene", "Hugo_Symbol", "Gene_Symbol"]

def detect_gene_column(df, fname):
    for col in possible_cols:
        if col in df.columns:
            print(f"✔ {fname}: Using gene column → {col}")
            return col
    raise ValueError(f"❌ No valid gene column found in {fname}")

col_comp = detect_gene_column(compendium, "Compendium_Cancer_Genes.tsv")
col_unf = detect_gene_column(unfiltered, "Unfiltered_drivers.tsv")

# ============================
# Convert to uppercase for matching
# ============================

top_genes_set = set(top_genes["Gene"].astype(str).str.upper())
compendium_genes = set(compendium[col_comp].astype(str).str.upper())
unfiltered_genes = set(unfiltered[col_unf].astype(str).str.upper())

# ============================
# MATCHING
# ============================

top_hits_in_compendium = top_genes_set.intersection(compendium_genes)
top_hits_in_unfiltered = top_genes_set.intersection(unfiltered_genes)

print("\n==== VALIDATION SUMMARY ====")
print("Top genes checked:", len(top_genes_set))
print("Matches in Compendium:", len(top_hits_in_compendium))
print("Matches in Unfiltered Drivers:", len(top_hits_in_unfiltered))
print("============================")

# Save matched lists
pd.DataFrame(sorted(top_hits_in_compendium), columns=["Gene"]).to_csv(
    "Matched_in_Compendium.csv", index=False
)

pd.DataFrame(sorted(top_hits_in_unfiltered), columns=["Gene"]).to_csv(
    "Matched_in_Unfiltered.csv", index=False
)

# ============================
# BAR CHART FOR VALIDATION
# ============================

labels = ["IntOGen Compendium", "Unfiltered Drivers"]
values = [len(top_hits_in_compendium), len(top_hits_in_unfiltered)]

plt.figure(figsize=(6,4))
sns.barplot(x=labels, y=values, palette="Blues")
plt.title("Validation of CRISPR-Immune Top Genes\nAgainst IntOGen Cancer Driver Databases")
plt.ylabel("Number of Overlapping Genes")
plt.tight_layout()
plt.savefig("Validation_BarChart.png", dpi=300)
plt.show()
