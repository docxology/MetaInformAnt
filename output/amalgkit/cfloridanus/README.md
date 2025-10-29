# Camponotus floridanus RNA-seq Analysis

**Species**: *Camponotus floridanus* (Florida Carpenter Ant)  
**NCBI Taxonomy ID**: 104421  
**Assembly**: GCF_003227725.1_Cflo_v7.5  
**Sequencing**: PacBio long-read sequencing, high-contiguity assembly  
**Status**: ⏳ **In Progress** - Workflow queued

---

## Quick Reference

See **[QUICK_REFERENCE.md](QUICK_REFERENCE.md)** for:
- Pre-built expression matrices (when ready)
- Common analysis tasks
- Dataset summary and file locations

---

## General Amalgkit Documentation

**For workflow guides, setup, and troubleshooting, see:**

- **Quick Start**: `docs/rna/amalgkit/quick_start.md` - Step-by-step guide
- **R Package Setup**: `docs/rna/amalgkit/r_packages.md` - R dependencies
- **Verification Script**: `verify_workflow.sh` - QC validation
- **Full Workflow**: `docs/rna/amalgkit/END_TO_END_WORKFLOW.md` - Complete guide

---

## Workflow Summary

### Steps
1. **Metadata** - Retrieve RNA-seq samples from NCBI SRA
2. **FASTQ Download** - Download via ENA (faster than NCBI)
3. **Quantification** - Kallisto pseudoalignment
4. **Merge** - Combine expression matrices
5. **Curate** - QC filtering with visualizations
6. **Sanity** - Data integrity validation

---

## File Organization

```
cfloridanus/
├── README.md                   # This file
├── QUICK_REFERENCE.md          # Data summary (when ready)
├── verify_workflow.sh          # Verification script
└── work/                       # Amalgkit working directory
    ├── metadata/
    ├── genome/
    ├── quant/
    ├── merge/
    ├── curate/
    └── sanity/
```

---

## Verification

Run the verification script after workflow completes:

```bash
cd output/amalgkit/cfloridanus
./verify_workflow.sh
```

---

## Notes

**Biological Significance**: *C. floridanus* is a model organism for epigenetics, caste determination, and aging research in social insects. The species has well-characterized caste polyphenism and has contributed significantly to understanding social insect biology.

**Assembly Quality**: High-contiguity assembly from PacBio long-read sequencing provides excellent reference for transcriptome analysis and gene structure annotation.

---

**Last Updated**: October 29, 2025  
**Status**: ⏳ Workflow queued  
**Documentation**: See `docs/rna/amalgkit/` for guides

