# Multi-Species Amalgkit Workflow Status

**Date**: October 29, 2025  
**Species**: *Solenopsis invicta*, *Camponotus floridanus*  
**Status**: ⏳ **Running**

---

## Overview

Full end-to-end amalgkit workflows have been initiated for two ant species:

1. **Solenopsis invicta** (Red Fire Ant)
2. **Camponotus floridanus** (Florida Carpenter Ant)

The workflows are running sequentially in the background and will take many hours to complete (estimated 18-24 hours per species based on P. barbatus experience).

---

## Configuration Validation

### ✅ S. invicta Configuration
- **Config File**: `config/amalgkit_sinvicta.yaml`
- **Taxonomy ID**: 13686
- **Assembly**: GCF_016802725.1_UNIL_Sinv_3.0
- **Sequencing**: PacBio HiFi + Hi-C (chromosome-level)
- **Work Directory**: `output/amalgkit/sinvicta/work`
- **Threads**: 6
- **Auto-install**: Yes
- **Steps**: 11 (metadata → sanity)
- **Filters**: Illumina, RNA-Seq, paired-end, 10M-400M spots

### ✅ C. floridanus Configuration
- **Config File**: `config/amalgkit_cfloridanus.yaml`
- **Taxonomy ID**: 104421
- **Assembly**: GCF_003227725.1_Cflo_v7.5
- **Sequencing**: PacBio long-read (high-contiguity)
- **Work Directory**: `output/amalgkit/cfloridanus/work`
- **Threads**: 6
- **Auto-install**: Yes
- **Steps**: 11 (metadata → sanity)
- **Filters**: Illumina, RNA-Seq, paired-end, 10M-400M spots

---

## Workflow Steps

Each species will go through these steps:

1. **Genome Download** - Download reference genome from NCBI
2. **Metadata Retrieval** - Find RNA-seq samples in SRA
3. **Integrate** - Integrate metadata with configuration
4. **Config** - Generate amalgkit config files
5. **Select** - Select samples based on filters
6. **GetFASTQ** - Download FASTQ files (via ENA for speed)
7. **Quant** - Kallisto quantification
8. **Merge** - Combine expression matrices
9. **CSTMM** - Cross-species TMM normalization
10. **Curate** - QC filtering and visualization
11. **CSCA** - Cross-species correlation analysis
12. **Sanity** - Final integrity validation

---

## Monitoring Progress

### Real-Time Monitoring
```bash
# Watch progress (updates every 60 seconds)
watch -n 60 ./scripts/rna/monitor_amalgkit_progress.sh

# View live log
tail -f output/amalgkit_multi_species_run.log

# Check specific species
ls -lhR output/amalgkit/sinvicta/work/
ls -lhR output/amalgkit/cfloridanus/work/
```

### Log Locations
- **Main log**: `output/amalgkit_multi_species_run.log`
- **S. invicta logs**: `output/amalgkit/sinvicta/logs/`
- **C. floridanus logs**: `output/amalgkit/cfloridanus/logs/`

---

## Expected Outputs

### For Each Species

#### Directory Structure
```
output/amalgkit/<species>/
├── README.md
├── QUICK_REFERENCE.md (created after completion)
├── verify_workflow.sh
└── work/
    ├── metadata/
    │   └── metadata_original.tsv
    ├── genome/
    │   ├── GCF_*_genomic.fna.gz
    │   ├── GCF_*_rna_from_genomic.fna.gz
    │   └── ...
    ├── quant/
    │   └── SRR*/
    │       ├── abundance.tsv
    │       ├── abundance.h5
    │       └── run_info.json
    ├── merge/
    │   ├── <species>/
    │   │   ├── *_tpm.tsv
    │   │   ├── *_est_counts.tsv
    │   │   └── *_eff_length.tsv
    │   └── metadata.tsv
    ├── curate/
    │   └── <species>/
    │       ├── plots/ (6 PDFs)
    │       ├── tables/ (7 TSVs)
    │       └── rdata/ (3 RData)
    └── sanity/
        └── SRA_IDs_without_fastq.txt
```

#### Key Output Files
- **Expression matrix (TPM)**: `work/curate/<species>/tables/<species>.no.tc.tsv`
- **QC visualizations**: `work/curate/<species>/plots/*.pdf`
- **Tissue specificity**: `work/curate/<species>/tables/<species>.no.tau.tsv`
- **Sample metadata**: `work/merge/metadata.tsv`

---

## Verification

After completion, run verification for each species:

```bash
# S. invicta
cd output/amalgkit/sinvicta
./verify_workflow.sh

# C. floridanus
cd output/amalgkit/cfloridanus
./verify_workflow.sh
```

---

## Documentation

### General Workflow Guides
All documentation has been reorganized to `docs/rna/amalgkit/`:

- **Quick Start**: `docs/rna/amalgkit/quick_start.md`
- **R Package Setup**: `docs/rna/amalgkit/r_packages.md`
- **Verification Template**: `docs/rna/amalgkit/verify_template.sh`
- **Full Workflow**: `docs/rna/amalgkit/END_TO_END_WORKFLOW.md`
- **Comprehensive Guide**: `docs/rna/amalgkit/comprehensive_guide.md`

### Species-Specific
- **S. invicta**: `output/amalgkit/sinvicta/README.md`
- **C. floridanus**: `output/amalgkit/cfloridanus/README.md`

---

## Biological Significance

### Solenopsis invicta
- **Common Name**: Red Fire Ant
- **Research Focus**: Invasive species biology, social organization plasticity, caste determination
- **Genome Quality**: Chromosome-level assembly (PacBio HiFi + Hi-C)
- **RNA-Seq Applications**: Behavioral studies, caste-specific expression, venom research

### Camponotus floridanus
- **Common Name**: Florida Carpenter Ant
- **Research Focus**: Epigenetics, aging, caste determination, circadian rhythms
- **Genome Quality**: High-contiguity assembly (PacBio long-read)
- **RNA-Seq Applications**: Caste polyphenism, developmental studies, aging research

---

## Timeline

### Expected Duration
- **S. invicta**: ~18-24 hours (running)
- **C. floridanus**: ~18-24 hours (queued)
- **Total**: ~36-48 hours

### Started
- **S. invicta**: October 29, 2025 12:22 PM

### Estimated Completion
- **S. invicta**: October 30, 2025 6:00 AM - 12:00 PM
- **C. floridanus**: October 31, 2025 12:00 AM - 6:00 AM
- **Both species**: October 31, 2025 morning

---

## Next Steps After Completion

1. **Verify outputs**: Run `verify_workflow.sh` for each species
2. **Create QUICK_REFERENCE.md**: Add species-specific statistics
3. **Analyze expression data**: Load matrices, identify top genes
4. **Cross-species comparison**: Compare expression patterns between species
5. **GO/KEGG enrichment**: Functional annotation analysis
6. **Publication**: Generate figures and results

---

## Troubleshooting

If workflows fail:

1. **Check logs**: `tail -100 output/amalgkit_multi_species_run.log`
2. **Check species logs**: `ls -lh output/amalgkit/<species>/logs/`
3. **Verify genome download**: `ls -lh output/amalgkit/<species>/genome/`
4. **Check disk space**: `df -h output/`
5. **Re-run specific step**: See `docs/rna/amalgkit/quick_start.md`

---

## Success Metrics

Expected for each species:
- ✅ Genome downloaded (8 files)
- ✅ Metadata retrieved (N samples)
- ✅ FASTQs downloaded (N × 2 files)
- ✅ Quantifications complete (N samples)
- ✅ Expression matrices (3 files: TPM, counts, eff_length)
- ✅ QC visualizations (6 PDFs)
- ✅ Curate data tables (7 TSVs)
- ✅ 100% sanity validation

---

**Status**: ⏳ Workflows running in background  
**Monitor**: `scripts/rna/monitor_amalgkit_progress.sh`  
**Log**: `output/amalgkit_multi_species_run.log`

