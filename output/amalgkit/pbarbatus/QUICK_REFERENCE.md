# Quick Reference: P. barbatus Expression Data

**Updated**: October 29, 2025  
**Status**: ✅ **Complete** (83/83 samples)

---

## 📊 Expression Matrix (Pre-built)

### Load Complete Dataset
```python
import pandas as pd

# Load TPM matrix (20,672 × 83)
tpm = pd.read_csv("analysis/expression_matrix_tpm.csv", index_col=0)

# Load counts
counts = pd.read_csv("analysis/expression_matrix_counts.csv", index_col=0)

# Load sample statistics
stats = pd.read_csv("analysis/sample_statistics.csv", index_col=0)
```

---

## 🎯 Common Tasks

### 1. Get Top Expressed Genes
```python
# Mean expression across all samples
mean_expr = tpm.mean(axis=1).sort_values(ascending=False)
print(mean_expr.head(20))

# Top in specific sample
print(tpm['SRR14740487'].nlargest(20))
```

### 2. Find Differentially Expressed Genes
```python
# Compare two groups (example)
group1 = tpm[['SRR14740487', 'SRR14740488', 'SRR14740489']]
group2 = tpm[['SRR14740490', 'SRR14740491', 'SRR14740492']]

fold_change = group1.mean(axis=1) / group2.mean(axis=1)
print(fold_change.nlargest(20))
```

### 3. Filter by Expression Level
```python
# Genes expressed in >50% of samples
expressed = (tpm > 1).sum(axis=1) > (tpm.shape[1] / 2)
expressed_genes = tpm[expressed]
print(f"Expressed genes: {len(expressed_genes)}")
```

### 4. Sample Correlation
```python
# Compute pairwise correlations
corr = tpm.corr()
print(f"Mean correlation: {corr.values[corr.values < 1].mean():.3f}")
```

---

## 📈 Visualizations

Pre-generated plots in `analysis/visualizations/`:

1. **sample_correlation_heatmap.png** - Sample similarity
2. **expression_distribution.png** - Expression patterns
3. **sample_qc_metrics.png** - Quality metrics
4. **top_expressed_genes.png** - Top 20 genes

```bash
# View all plots
open analysis/visualizations/*.png
```

---

## 📊 Dataset Summary

| Metric | Value |
|--------|-------|
| **Samples** | 83 |
| **Transcripts** | 20,672 |
| **Expressed** | 17,191 (83.2%) |
| **Mean genes/sample** | 16,372 |
| **Top transcript** | XM_011631231.1 (25,871 TPM) |

---

## 📁 File Locations

```
analysis/
├── expression_matrix_tpm.csv        ← Main data (12.7MB)
├── expression_matrix_counts.csv     ← Count data (11.2MB)
├── sample_statistics.csv            ← QC metrics
├── ANALYSIS_REPORT.md               ← Full report
└── visualizations/                  ← 4 plots

work/
├── quant/SRR*/abundance.tsv         ← Individual samples (83 files)
└── merge/metadata.tsv               ← Sample metadata
```

---

## 🚀 Next Steps

### Immediate
- ✅ Data ready for differential expression
- ✅ Ready for GO/KEGG enrichment
- ✅ Ready for co-expression networks

### Advanced
- Add functional annotations
- Integrate with other ant species
- Perform tissue-specific analysis

---

## 📖 Documentation

**Full analysis report**: `analysis/ANALYSIS_REPORT.md`  
**Complete workflow**: `docs/rna/amalgkit/END_TO_END_WORKFLOW.md`  
**Directory guide**: `README.md`

---

*Complete dataset with analysis ready for downstream use*