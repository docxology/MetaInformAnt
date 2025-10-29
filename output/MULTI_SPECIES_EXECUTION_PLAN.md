# Multi-Species Amalgkit Execution Plan

**Date**: October 29, 2025  
**Species Count**: 4  
**Execution Mode**: Sequential with cross-species analysis

---

## Species to Process

### 1. Pogonomyrmex barbatus (Red Harvester Ant)
- **Config**: `config/amalgkit_pbarbatus.yaml`
- **Taxonomy ID**: 144034
- **Assembly**: GCF_000187915.1 (Pbar_UMD_V03)
- **Status**: ✅ Already complete (83 samples processed)
- **Action**: Re-validate with sanity check

### 2. Solenopsis invicta (Red Fire Ant)
- **Config**: `config/amalgkit_sinvicta.yaml`
- **Taxonomy ID**: 13686
- **Assembly**: GCF_016802725.1 (UNIL_Sinv_3.0)
- **Status**: ⏳ To be processed
- **Est. Duration**: 18-24 hours

### 3. Camponotus floridanus (Florida Carpenter Ant)
- **Config**: `config/amalgkit_cfloridanus.yaml`
- **Taxonomy ID**: 104421
- **Assembly**: GCF_003227725.1 (Cflo_v7.5)
- **Status**: ⏳ To be processed
- **Est. Duration**: 18-24 hours

### 4. Monomorium pharaonis (Pharaoh Ant)
- **Config**: `config/amalgkit_mpharaonis.yaml`
- **Taxonomy ID**: 307658
- **Assembly**: GCF_013373865.1 (ASM1337386v2)
- **Status**: ⏳ To be processed
- **Est. Duration**: 18-24 hours

---

## Workflow Phases

### Phase 1: Individual Species Workflows (54-72 hours total)

Each species goes through these steps:

1. **Genome Download** (5 min) - Download reference from NCBI
2. **Metadata Retrieval** (2-5 min) - Search SRA for RNA-seq samples
3. **Integrate** (1 min) - Merge metadata with config
4. **Config** (1 min) - Generate amalgkit configs
5. **Select** (1 min) - Apply quality filters
6. **GetFASTQ** (2-6 hours) - Download via ENA
7. **Quant** (6-12 hours) - Kallisto quantification
8. **Merge** (5-10 min) - Combine expression matrices
9. **CSTMM** (2-5 min) - Cross-species prep (individual)
10. **Curate** (5-15 min) - QC filtering + 6 PDF visualizations
11. **CSCA** (2-5 min) - Correlation prep (individual)
12. **Sanity** (1 min) - Final validation

### Phase 2: Cross-Species Analysis (30-60 min)

After all individual workflows complete:

1. **CSTMM (Cross-Species TMM)** - Normalize expression using single-copy genes across all 4 species
2. **CSCA (Cross-Species Correlation)** - Generate cross-species correlation plots and analysis

---

## Expected Outputs

### Per Species (4 species × following structure)

```
output/amalgkit/<species>/
├── README.md
├── QUICK_REFERENCE.md
├── verify_workflow.sh
└── work/
    ├── metadata/metadata_original.tsv
    ├── genome/[8 reference files]
    ├── quant/SRR*/[3 files per sample]
    ├── merge/<species>/
    │   ├── *_tpm.tsv (expression matrix)
    │   ├── *_est_counts.tsv
    │   └── *_eff_length.tsv
    └── curate/<species>/
        ├── plots/[6 PDFs]
        ├── tables/[7 TSVs]
        └── rdata/[3 RData]
```

### Cross-Species Outputs

```
output/amalgkit/cross_species/
├── cstmm/
│   ├── normalized_expression.tsv
│   ├── single_copy_genes.txt
│   └── normalization_factors.tsv
└── csca/
    ├── correlation_matrix.tsv
    ├── correlation_heatmap.pdf
    ├── species_dendrogram.pdf
    └── cross_species_pca.pdf
```

---

## Timeline Estimate

| Phase | Duration | Completion |
|-------|----------|------------|
| **P. barbatus** (re-validate) | 5 min | Oct 29, 12:35 PM |
| **S. invicta** | 18-24 hours | Oct 30, 12:00 PM |
| **C. floridanus** | 18-24 hours | Oct 31, 12:00 PM |
| **M. pharaonis** | 18-24 hours | Nov 1, 12:00 PM |
| **Cross-species** | 30-60 min | Nov 1, 1:00 PM |
| **TOTAL** | ~60-75 hours | **Nov 1, 2025** |

---

## Biological Value

### Phylogenetic Coverage
- **Myrmicinae**: *P. barbatus*, *S. invicta*, *M. pharaonis*
- **Formicinae**: *C. floridanus*
- Spans multiple genera within Formicidae

### Research Themes
1. **Caste determination** - All 4 species
2. **Social organization** - Colony structure variation
3. **Invasive species** - *S. invicta*, *M. pharaonis*
4. **Epigenetics** - *C. floridanus*, *M. pharaonis*
5. **Alternative splicing** - *M. pharaonis*
6. **Aging** - *C. floridanus*

### Cross-Species Questions
- Conserved caste-specific expression patterns?
- Subfamily-specific gene expression signatures?
- Shared mechanisms of social regulation?
- Evolution of eusociality across lineages?

---

## Monitoring

### Real-Time Progress
```bash
# Watch all species
watch -n 60 ./scripts/rna/monitor_amalgkit_progress.sh

# View live log
tail -f output/amalgkit_multi_species_run.log

# Check individual species
ls -lhR output/amalgkit/{pbarbatus,sinvicta,cfloridanus,mpharaonis}/work/
```

### Log Files
- **Main**: `output/amalgkit_multi_species_run.log`
- **Per species**: `output/amalgkit/<species>/logs/`
- **Cross-species**: `output/amalgkit/cross_species/logs/`

---

## Success Criteria

### Individual Species (each)
- ✅ Genome downloaded (8 files)
- ✅ Metadata retrieved (N samples)
- ✅ All FASTQs downloaded
- ✅ All quantifications complete
- ✅ Expression matrices generated (3 files)
- ✅ QC visualizations (6 PDFs)
- ✅ Curate tables (7 TSVs)
- ✅ 100% sanity validation

### Cross-Species
- ✅ Single-copy genes identified
- ✅ TMM normalization complete
- ✅ Correlation matrix generated
- ✅ Cross-species plots created
- ✅ Comparative analysis ready

---

## Execution Command

```bash
cd /Users/4d/Documents/GitHub/metainformant
python3 scripts/rna/run_multi_species_amalgkit.py 2>&1 | tee output/amalgkit_multi_species_run.log
```

---

**Ready to Execute**: ✅  
**Auto-Discovery**: ✅ (4 species found)  
**Cross-Species Enabled**: ✅  
**Estimated Completion**: November 1, 2025 afternoon

