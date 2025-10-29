# Monomorium pharaonis RNA-seq Analysis

**Species**: *Monomorium pharaonis* (Pharaoh Ant)  
**NCBI Taxonomy ID**: 307658  
**Assembly**: GCF_013373865.1_ASM1337386v2  
**Sequencing**: Illumina + PacBio + Hi-C, chromosome-level assembly  
**Status**: ⏳ **Pending** - Workflow ready to run

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
mpharaonis/
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
cd output/amalgkit/mpharaonis
./verify_workflow.sh
```

---

## Notes

**Biological Significance**: *M. pharaonis* is an indoor pest species and model organism for studying caste determination and alternative splicing in social insects. Despite its small size, it shows remarkable behavioral plasticity and has contributed to understanding social insect genomics.

**Assembly Quality**: Chromosome-level assembly combining Illumina, PacBio, and Hi-C technologies provides excellent contiguity and completeness for transcriptome analysis.

**Research Applications**: Alternative splicing, caste polyphenism, indoor adaptation, chemical communication

---

**Last Updated**: October 29, 2025  
**Status**: ⏳ Workflow ready  
**Documentation**: See `docs/rna/amalgkit/` for guides

