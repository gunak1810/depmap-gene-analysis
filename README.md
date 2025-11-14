# 🧬 In-Silico Pan-Cancer Profiling of Immune-Linked Gene Dependencies

A computational framework integrating **DepMap CRISPR gene-effect** data with **bulk RNA-seq expression** to identify **cancer-type-specific** and **recurrent immune-linked genetic dependencies** across human cancers.

This repository contains all Python scripts required to reproduce the full workflow, including **correlation analysis**, **pan-cancer enrichment**, and **validation against cancer driver databases**.

---

## 📁 Repository Overview

### Main Scripts
| Script | Purpose |
|--------|---------|
| **gene_analysis.py** | Core correlation analysis for each cancer type. |
| **post_analysis_enrichment.py** | Cross-cancer overlap, GO/KEGG enrichment, and global comparison. |
| **post_analysis_validation.py** | Validation of recurrent hits using curated driver gene databases (IntOGen). |

---

## 1️⃣ Requirements

### 🔧 Software
- Python **3.8+**

### 📦 Python Packages
Install dependencies using:

```bash
pip install pandas numpy scipy statsmodels matplotlib seaborn gprofiler-official
```

> **Note:**  
> `post_analysis_enrichment.py` automatically installs `gprofiler-official` if missing.

---

## 2️⃣ Data Acquisition

All required raw datasets must be downloaded from **DepMap (depmap.org)** and placed in the **project root directory**.

### 📄 Required DepMap Files
| Filename | Description |
|----------|-------------|
| **Model.csv** | Cell line metadata including Oncotree cancer type. |
| **CRISPRGeneEffect.csv** | CRISPR knockout viability scores. |
| **OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv** | Log-TPM expression data. |

### 📄 Validation Files (Step 3)
Place these in **Enrichment_and_Overlap_Results/**:

| Filename | Description |
|----------|-------------|
| **Compendium_Cancer_Genes.tsv** | Curated driver gene list. |
| **Unfiltered_drivers.tsv** | Predicted/unfiltered driver genes. |

---

## 3️⃣ Project Structure

```
.
├── gene_analysis.py
├── post_analysis_enrichment.py
├── post_analysis_validation.py
├── Model.csv
├── CRISPRGeneEffect.csv
├── OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv
└── ...
```

Automatic folders:
- `Cancer_Specific_Results/`
- `Enrichment_and_Overlap_Results/`

---

## 4️⃣ Workflow & Execution

### 🔹 Step 1 — Cancer-Specific CRISPR-Immune Correlation

```bash
python gene_analysis.py
```

Outputs:
- `CRISPR_TME_<CancerType>_Results.csv`
- `Volcano_<CancerType>.png`
- `Cancerwise_Summary.csv`

---

### 🔹 Step 2 — Pan-Cancer Overlap & Enrichment

```bash
python post_analysis_enrichment.py
```

Outputs:
- `Global_Top_Positive_Genes.csv`
- `Jaccard_Positive_Heatmap.png`
- `Presence_Positive_Heatmap.png`
- Enrichment result files

---

### 🔹 Step 3 — Driver Gene Validation

```bash
python post_analysis_validation.py
```

Outputs:
- `Matched_in_Compendium.csv`
- `Matched_in_Unfiltered.csv`
- `Validation_BarChart.png`

---

## 5️⃣ Key Findings

- Immune-linked CRISPR dependencies are **highly cancer-specific**.
- Several **noncanonical genes** (PRH2, TFB1M, NSUN4, TAOK1) strongly correlate with immune activity.
- Many recurrent hits are outside curated driver lists → indicating **novel biology** and **potential therapeutic targets**.

## 📫 Contact / Contribution

Feel free to open issues or submit pull requests for improvements or extensions to the analytical pipeline.
