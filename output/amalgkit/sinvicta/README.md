# Solenopsis invicta RNA-seq Analysis

**Species**: *Solenopsis invicta* (Red Fire Ant)  
**NCBI Taxonomy ID**: 13686  
**Assembly**: GCF_016802725.1_UNIL_Sinv_3.0  
**Sequencing**: PacBio HiFi + Hi-C scaffolding, chromosome-level assembly  
**Status**: ⏳ **In Progress** - Workflow running

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
sinvicta/
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
cd output/amalgkit/sinvicta
./verify_workflow.sh
```

---

## Notes

**Biological Significance**: *S. invicta* is an invasive eusocial species extensively studied for behavioral and caste determination. The species shows remarkable plasticity in social organization and has been a model for understanding social insect genomics.

**Assembly Quality**: Chromosome-level assembly with PacBio HiFi + Hi-C provides excellent contiguity for transcriptome analysis.

---

**Last Updated**: October 29, 2025  
**Status**: ⏳ Workflow in progress  
**Documentation**: See `docs/rna/amalgkit/` for guides

