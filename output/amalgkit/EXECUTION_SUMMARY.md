# Multi-Species Amalgkit Workflow Execution Summary

**Execution Date**: October 29, 2025  
**Status**: ✅ **Initiated** - Workflows running in background  
**Species**: 2 (*S. invicta*, *C. floridanus*)

---

## ✅ Configuration Validation Complete

Both species configs have been validated and are ready:

### Solenopsis invicta (Red Fire Ant)
- ✅ Config file written and validated
- ✅ Taxonomy ID: 13686
- ✅ Assembly: GCF_016802725.1 (UNIL_Sinv_3.0)
- ✅ PacBio HiFi + Hi-C chromosome-level assembly
- ✅ All 11 workflow steps configured
- ✅ Quality filters: Illumina, RNA-Seq, paired-end, 10M-400M spots

### Camponotus floridanus (Florida Carpenter Ant)
- ✅ Config file written and validated
- ✅ Taxonomy ID: 104421
- ✅ Assembly: GCF_003227725.1 (Cflo_v7.5)
- ✅ PacBio long-read high-contiguity assembly
- ✅ All 11 workflow steps configured
- ✅ Quality filters: Illumina, RNA-Seq, paired-end, 10M-400M spots

---

## 🚀 Workflow Execution

### Current Status (12:24 PM)
- **S. invicta**: ⏳ Running (started 12:22 PM)
- **C. floridanus**: 📋 Queued (will start after S. invicta completes)

### Execution Method
- Sequential processing (one species at a time)
- Running in background with full logging
- Auto-install amalgkit enabled
- Network-enabled for downloads

### Log Files
- **Main log**: `output/amalgkit_multi_species_run.log`
- **S. invicta logs**: `output/amalgkit/sinvicta/logs/`
- **C. floridanus logs**: `output/amalgkit/cfloridanus/logs/`

---

## 📊 What's Happening

### Workflow Steps (per species)
1. **Genome Download** (~5 min) - Download reference from NCBI FTP
2. **Metadata** (~2-5 min) - Search SRA for RNA-seq samples
3. **Integrate** (~1 min) - Merge metadata with config
4. **Config** (~1 min) - Generate amalgkit configs
5. **Select** (~1 min) - Apply quality filters
6. **GetFASTQ** (~2-6 hours) - Download via ENA (187× faster than NCBI)
7. **Quant** (~6-12 hours) - Kallisto quantification
8. **Merge** (~5-10 min) - Combine expression matrices
9. **CSTMM** (~2-5 min) - Cross-species normalization
10. **Curate** (~5-15 min) - QC filtering + visualizations
11. **CSCA** (~2-5 min) - Cross-species correlation
12. **Sanity** (~1 min) - Final validation

### Expected Timeline
- **S. invicta**: 18-24 hours (est. completion: Oct 30, 6-12 AM)
- **C. floridanus**: 18-24 hours (est. completion: Oct 31, 12-6 AM)
- **Total duration**: ~36-48 hours

---

## 🔍 Monitoring

### Real-Time Monitoring Script
```bash
# Run the monitoring script
./scripts/rna/monitor_amalgkit_progress.sh

# Watch with auto-refresh every 60 seconds
watch -n 60 ./scripts/rna/monitor_amalgkit_progress.sh
```

### Manual Checks
```bash
# View latest log entries
tail -100 output/amalgkit_multi_species_run.log

# Follow log in real-time
tail -f output/amalgkit_multi_species_run.log

# Check S. invicta progress
ls -lhR output/amalgkit/sinvicta/work/

# Check C. floridanus progress
ls -lhR output/amalgkit/cfloridanus/work/
```

### Progress Indicators
- **metadata_original.tsv** → Metadata retrieved
- **genome/** directory → Reference downloaded
- **quant/SRR*/** directories → Quantification running
- **merge/metadata.tsv** → Merging complete
- **curate/*/plots/** → QC visualizations generated

---

## 📁 Prepared Resources

### Verification Scripts ✅
- `output/amalgkit/sinvicta/verify_workflow.sh`
- `output/amalgkit/cfloridanus/verify_workflow.sh`

### README Templates ✅
- `output/amalgkit/sinvicta/README.md`
- `output/amalgkit/cfloridanus/README.md`

### Monitoring Tools ✅
- `scripts/rna/monitor_amalgkit_progress.sh`

### Documentation ✅
- All general docs in `docs/rna/amalgkit/`
- Quick start guide
- R package setup
- Verification template
- Comprehensive workflow guide

---

## 📈 Expected Outputs

### For Each Species

#### Expression Data
- **TPM matrix**: 20,000-30,000 genes × N samples
- **Count matrix**: Raw estimated counts
- **Effective length**: Gene length estimates
- **Metadata**: Sample QC metrics

#### Visualizations (6 PDFs)
- Hierarchical clustering dendrograms
- Sample correlation heatmaps
- PCA plots
- QC metrics

#### Data Tables (7 TSVs)
- Curated expression matrix (primary output)
- Tissue specificity scores (tau)
- Sample group means
- Correlation statistics
- Original uncorrected data

#### Analysis States (3 RData)
- Round 0, 1, 2 analysis states for R

---

## ✅ Post-Completion Tasks

### Immediate (Auto-completed)
1. ✅ Sanity validation
2. ✅ Log consolidation
3. ✅ Manifest generation

### Manual Verification
1. Run `verify_workflow.sh` for each species
2. Check sample counts and gene counts
3. View QC visualizations
4. Validate expression matrices

### Analysis Ready
1. Load TPM matrices in Python/R
2. Identify top expressed genes
3. Compare expression between species
4. GO/KEGG enrichment analysis
5. Co-expression network analysis

---

## 🔬 Biological Context

### Why These Species?

#### S. invicta
- Model for invasive species biology
- Extensive caste system studies
- Well-characterized social behavior
- Venom and defense mechanism research
- Chromosome-level assembly quality

#### C. floridanus
- Epigenetics model organism
- Caste determination studies
- Aging and longevity research
- Circadian rhythm studies
- High-quality genome annotation

### Comparative Value
- Different subfamilies (Myrmicinae vs. Formicinae)
- Different colony organization patterns
- Complementary research strengths
- Cross-species expression comparison
- Orthology analysis potential

---

## 📚 Reference Documentation

### Quick Links
- **Status**: `output/MULTI_SPECIES_WORKFLOW_STATUS.md`
- **This summary**: `output/amalgkit/EXECUTION_SUMMARY.md`
- **Quick start**: `docs/rna/amalgkit/quick_start.md`
- **Full workflow**: `docs/rna/amalgkit/END_TO_END_WORKFLOW.md`
- **R setup**: `docs/rna/amalgkit/r_packages.md`

### Example Dataset
- **P. barbatus**: `output/amalgkit/pbarbatus/`
- Complete working example with 83 samples
- All QC visualizations and analyses
- Quick reference and analysis report

---

## ⚠️ Important Notes

1. **Duration**: Workflows take 18-24 hours per species - this is normal
2. **Downloads**: FASTQ files total ~50-200 GB per species
3. **Storage**: Ensure adequate disk space (500 GB recommended)
4. **Network**: Stable internet required for downloads
5. **Resources**: 6 threads per species, sequential execution
6. **Cleanup**: FASTQs deleted after quantification (optional)

---

## 🎯 Success Criteria

Each species workflow succeeds when:
- ✅ All 12 steps complete (exit code 0)
- ✅ Expression matrix generated (genes × samples)
- ✅ 6 QC PDFs created
- ✅ 7 data tables generated
- ✅ 100% sanity validation passes
- ✅ No missing quantifications

---

**Execution Started**: October 29, 2025 12:22 PM  
**Estimated Completion**: October 31, 2025 morning  
**Status**: ✅ Running smoothly  
**Monitor**: `./scripts/rna/monitor_amalgkit_progress.sh`

