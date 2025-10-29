# Quick Reference: Accessing P. barbatus Gene Expression Data

**Location**: `output/amalgkit/pbarbatus/`  
**Updated**: October 28, 2025

---

## 📊 Where to Find Expression Values

### Individual Sample Files

```bash
# Sample 1
output/amalgkit/pbarbatus/work/quant/SRR14740487/abundance.tsv

# Sample 2
output/amalgkit/pbarbatus/work/quant/SRR14740488/abundance.tsv
```

**Format**: TSV with columns: `target_id`, `length`, `eff_length`, `est_counts`, `tpm`  
**Use**: `tpm` column for normalized expression values

---

## 💻 Load in Python

```python
import pandas as pd
from pathlib import Path

# Load one sample
base = Path("output/amalgkit/pbarbatus/work/quant")
df = pd.read_csv(base / "SRR14740487" / "abundance.tsv", sep='\t')

# View top expressed genes
print(df.nlargest(20, 'tpm')[['target_id', 'tpm', 'est_counts']])

# Merge both samples
df1 = pd.read_csv(base / "SRR14740487" / "abundance.tsv", sep='\t')
df2 = pd.read_csv(base / "SRR14740488" / "abundance.tsv", sep='\t')

expr_matrix = df1[['target_id', 'tpm']].merge(
    df2[['target_id', 'tpm']],
    on='target_id',
    suffixes=('_SRR14740487', '_SRR14740488')
)
```

---

## 🎯 Current Status

| Metric | Value |
|--------|-------|
| **Samples quantified** | 2/83 (2.4%) |
| **Transcripts** | 20,672 |
| **Expression rate** | 89.9-90.6% |
| **Mapping rate** | 64.6-65.8% |
| **Index ready** | ✅ Yes |

---

## 🚀 Next Steps

Process remaining 81 samples using the pre-built index at:
```
output/amalgkit/pbarbatus/work/index/Pogonomyrmex_barbatus.idx
```

See `README.md` for full processing instructions.

---

## 📁 Directory Structure

```
pbarbatus/
├── README.md              ← Full documentation
├── QUICK_REFERENCE.md     ← This file
└── work/
    ├── quant/
    │   ├── SRR14740487/abundance.tsv    ← Expression values here
    │   └── SRR14740488/abundance.tsv    ← Expression values here
    ├── index/                            ← Kallisto index (ready to use)
    └── metadata/metadata.tsv             ← 83 samples to process
```

---

**All variant folders removed** - This is now the single consolidated directory for P. barbatus analysis.

