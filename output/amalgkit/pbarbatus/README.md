# Pogonomyrmex barbatus RNA-seq Analysis

**Species**: *Pogonomyrmex barbatus* (red harvester ant)  
**Tissue**: Brain  
**Samples**: 83 RNA-seq runs from NCBI SRA  
**Status**: ✅ **Complete** - Full workflow with visualizations

---

## Quick Reference

See **[QUICK_REFERENCE.md](QUICK_REFERENCE.md)** for:
- Pre-built expression matrices ready to load
- Common analysis tasks (top genes, differential expression, filtering)
- Dataset summary and file locations

---

## General Amalgkit Documentation

**For workflow guides, setup, and troubleshooting, see:**

- **Quick Start**: `docs/rna/amalgkit/quick_start.md` - Step-by-step guide
- **R Package Setup**: `docs/rna/amalgkit/r_packages.md` - R dependencies
- **Verification Script**: `docs/rna/amalgkit/verify_template.sh` - Template for QC
- **Full Workflow**: `docs/rna/amalgkit/END_TO_END_WORKFLOW.md` - Complete guide
- **Comprehensive Guide**: `docs/rna/amalgkit/comprehensive_guide.md` - All details

---

## Workflow Summary

### Steps Completed ✅

1. **Metadata** - 83 brain RNA-seq samples from NCBI SRA
2. **FASTQ Download** - ENA (187× faster than NCBI)
3. **Quantification** - Kallisto pseudoalignment
4. **Merge** - Combined expression matrices (20,672 genes × 83 samples)
5. **Curate** - QC filtering with 6 PDF visualizations
6. **Sanity** - 100% data integrity validation

---

## File Organization

```
pbarbatus/
├── QUICK_REFERENCE.md          # Species-specific data summary
├── README.md                   # This file
├── analysis/                   # Pre-built expression matrices 📊
│   ├── expression_matrix_tpm.csv
│   ├── expression_matrix_counts.csv
│   ├── sample_statistics.csv
│   ├── ANALYSIS_REPORT.md
│   └── visualizations/
└── work/                       # Amalgkit working directory
    ├── metadata/
    ├── genome/
    ├── quant/                  # 83 sample quantifications
    ├── merge/                  # Combined matrices
    ├── curate/                 # QC-filtered outputs 🎨
    │   └── Pogonomyrmex_barbatus/
    │       ├── plots/          # 6 PDF visualizations
    │       ├── tables/         # 7 TSV data files
    │       └── rdata/          # 3 R analysis states
    └── sanity/                 # Integrity validation
```

---

## Dataset Statistics

| Metric | Value |
|--------|-------|
| **Samples** | 83 |
| **Genes Quantified** | 20,672 |
| **Total Measurements** | 1,695,104 |
| **QC Passed** | 83/83 (100%) |
| **Mapping Rate** | >20% all samples |
| **Correlation** | High (no outliers) |

---

## Key Outputs

### Primary Expression Matrix
```
work/curate/Pogonomyrmex_barbatus/tables/Pogonomyrmex_barbatus.no.tc.tsv
```
- 20,672 genes × 83 samples
- Quality-controlled and filtered
- Ready for differential expression analysis

### Visualizations (6 PDFs)
```
work/curate/Pogonomyrmex_barbatus/plots/
├── *.0.original.pdf              # Original data
├── *.1.mapping_cutoff.pdf        # After mapping filter
└── *.2.correlation_cutoff.pdf    # Final curated (⭐ recommended)
```

Each PDF contains:
- Hierarchical clustering dendrogram
- Sample correlation heatmap
- PCA plot

### View All PDFs
```bash
cd output/amalgkit/pbarbatus
open work/curate/Pogonomyrmex_barbatus/plots/*.pdf
```

---

## Verification

Run the verification script (copy template from docs):

```bash
cd output/amalgkit/pbarbatus
cp ../../docs/rna/amalgkit/verify_template.sh ./verify_workflow.sh
chmod +x verify_workflow.sh
./verify_workflow.sh
```

Or run individual commands:

```bash
# Validate data integrity
amalgkit sanity --out_dir work --all

# Check curate outputs
find work/curate -type f | wc -l  # Should show 17
```

---

## Next Steps

### Immediate Use
- ✅ Load expression matrices from `analysis/` directory
- ✅ Use curated data for differential expression
- ✅ Perform GO/KEGG enrichment analysis
- ✅ Build co-expression networks

### Advanced Analysis
- Tissue-specific gene identification (tau scores in `*.no.tau.tsv`)
- Sample grouping comparison (means in `*.no.sample_group.mean.tsv`)
- Cross-species orthology analysis (requires additional data)

---

## Citation

If using this dataset or workflow:

- **Amalgkit**: https://github.com/kfuku52/amalgkit
- **Kallisto**: Bray et al. (2016). *Nature Biotechnology* 34, 525-527
- **Dataset**: *P. barbatus* brain RNA-seq from NCBI SRA

---

**Last Updated**: October 29, 2025  
**Status**: ✅ Complete and production-ready  
**Documentation**: See `docs/rna/amalgkit/` for workflow guides
