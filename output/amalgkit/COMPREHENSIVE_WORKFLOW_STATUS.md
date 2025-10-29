# Comprehensive Multi-Species Amalgkit Workflow Status

**Execution Date**: October 29, 2025  
**Status**: ✅ **RUNNING** - All 4 species + cross-species analysis  
**Updated**: Auto-discovery enabled

---

## 🎯 Workflow Overview

### Phase 1: Individual Species (Sequential)
1. **P. barbatus** - Red Harvester Ant (re-validate existing)
2. **S. invicta** - Red Fire Ant (full workflow)
3. **C. floridanus** - Florida Carpenter Ant (full workflow)
4. **M. pharaonis** - Pharaoh Ant (full workflow)

### Phase 2: Cross-Species Analysis
- **CSTMM**: Cross-species TMM normalization using single-copy genes
- **CSCA**: Cross-species correlation analysis with comparative plots

---

## ✅ What's Different Now

### Auto-Discovery
- ✅ Script now automatically finds all `amalgkit_*.yaml` configs
- ✅ Excludes `amalgkit_template.yaml`
- ✅ Processes all discovered species sequentially
- ✅ No hard-coded species list

### Cross-Species Analysis
- ✅ After all individual workflows complete
- ✅ Uses merged expression data from all species
- ✅ Identifies conserved single-copy genes
- ✅ Normalizes across species using TMM
- ✅ Generates cross-species correlation plots

### Species Discovered
```
✅ Cfloridanus  (C. floridanus)
✅ Mpharaonis   (M. pharaonis) ← NEW!
✅ Pbarbatus    (P. barbatus)
✅ Sinvicta     (S. invicta)
```

---

## 📊 Species Details

### 1. Pogonomyrmex barbatus
- **Status**: ✅ Complete (83 samples already processed)
- **Action**: Re-validate with sanity check
- **Assembly**: GCF_000187915.1
- **Data**: 20,672 genes × 83 samples

### 2. Solenopsis invicta
- **Status**: ⏳ Processing
- **Assembly**: GCF_016802725.1 (PacBio HiFi + Hi-C)
- **Significance**: Invasive species, social organization
- **Est. Duration**: 18-24 hours

### 3. Camponotus floridanus
- **Status**: 📋 Queued
- **Assembly**: GCF_003227725.1 (PacBio long-read)
- **Significance**: Epigenetics, aging, caste determination
- **Est. Duration**: 18-24 hours

### 4. Monomorium pharaonis
- **Status**: 📋 Queued
- **Assembly**: GCF_013373865.1 (Illumina + PacBio + Hi-C)
- **Significance**: Alternative splicing, indoor pest, caste studies
- **Est. Duration**: 18-24 hours

---

## 🔬 Cross-Species Analysis Design

### Scientific Questions
1. Are there conserved caste-specific expression patterns across species?
2. What are subfamily-specific gene expression signatures?
3. Which genes show conserved expression in social insects?
4. How does expression variation correlate with phylogenetic distance?

### Methods

#### CSTMM (Cross-Species TMM Normalization)
- Identifies single-copy orthologous genes across all 4 species
- Normalizes expression using TMM method
- Accounts for library size and composition bias
- Enables direct cross-species expression comparison

#### CSCA (Cross-Species Correlation Analysis)
- Computes correlation matrices across species
- Generates hierarchical clustering dendrograms
- Creates cross-species PCA plots
- Identifies co-expressed gene modules

### Expected Outputs
```
output/amalgkit/cross_species/
├── cstmm/
│   ├── normalized_expression_matrix.tsv    # All species, single-copy genes
│   ├── single_copy_gene_list.txt           # Conserved orthologs
│   ├── normalization_factors.tsv           # TMM factors per species
│   └── species_summary_stats.tsv
└── csca/
    ├── correlation_matrix.tsv              # Species × species correlations
    ├── correlation_heatmap.pdf             # Visual correlation matrix
    ├── species_dendrogram.pdf              # Phylogenetic relationships
    ├── cross_species_pca.pdf               # PCA across all samples
    └── coexpression_modules.tsv            # Shared expression modules
```

---

## 📁 Complete Directory Structure

```
output/amalgkit/
├── pbarbatus/                   # P. barbatus (COMPLETE)
│   ├── README.md
│   ├── QUICK_REFERENCE.md
│   ├── verify_workflow.sh
│   ├── analysis/                # Pre-built matrices
│   └── work/                    # Full workflow outputs
│
├── sinvicta/                    # S. invicta (RUNNING)
│   ├── README.md
│   ├── verify_workflow.sh
│   └── work/                    # Workflow in progress
│
├── cfloridanus/                 # C. floridanus (QUEUED)
│   ├── README.md
│   ├── verify_workflow.sh
│   └── work/                    # Will be created
│
├── mpharaonis/                  # M. pharaonis (QUEUED)
│   ├── README.md
│   ├── verify_workflow.sh
│   └── work/                    # Will be created
│
└── cross_species/               # Cross-species analysis (PENDING)
    ├── cstmm/                   # TMM normalization
    └── csca/                    # Correlation analysis
```

---

## ⏱️ Timeline

| Species/Phase | Start | Duration | Est. Complete |
|---------------|-------|----------|---------------|
| **P. barbatus** | 12:30 PM Oct 29 | 5 min | 12:35 PM Oct 29 |
| **S. invicta** | 12:35 PM Oct 29 | 18-24 hrs | 12:00 PM Oct 30 |
| **C. floridanus** | 12:00 PM Oct 30 | 18-24 hrs | 12:00 PM Oct 31 |
| **M. pharaonis** | 12:00 PM Oct 31 | 18-24 hrs | 12:00 PM Nov 1 |
| **Cross-species** | 12:00 PM Nov 1 | 30-60 min | 1:00 PM Nov 1 |
| **TOTAL** | - | **60-75 hours** | **Nov 1, 2025** |

---

## 📈 Monitoring

### Real-Time Monitoring
```bash
# Auto-refreshing monitor
watch -n 60 ./scripts/rna/monitor_amalgkit_progress.sh

# Live log
tail -f output/amalgkit_multi_species_comprehensive.log

# Check specific species
ls -lhR output/amalgkit/{pbarbatus,sinvicta,cfloridanus,mpharaonis}/work/
```

### Progress Indicators
- `work/metadata/metadata_original.tsv` → Samples found
- `work/genome/` → Reference downloaded
- `work/quant/SRR*/` → Quantifications running
- `work/merge/metadata.tsv` → Merging complete
- `work/curate/*/plots/*.pdf` → QC complete

---

## 🎓 Scientific Value

### Phylogenetic Sampling
- **Myrmicinae** (3 species): *P. barbatus*, *S. invicta*, *M. pharaonis*
- **Formicinae** (1 species): *C. floridanus*
- Multiple genera, single family (Formicidae)
- Good phylogenetic spread for comparative analysis

### Research Themes Covered
1. **Caste determination** - All species
2. **Social organization** - Colony structure variation
3. **Invasive biology** - *S. invicta*, *M. pharaonis*
4. **Epigenetics** - *C. floridanus*, *M. pharaonis*
5. **Alternative splicing** - *M. pharaonis*
6. **Aging mechanisms** - *C. floridanus*
7. **Chemical communication** - All species

### Comparative Power
- **Within-subfamily**: Compare 3 Myrmicinae species
- **Between-subfamily**: Compare Myrmicinae vs. Formicinae
- **Invasive vs. native**: Compare ecological strategies
- **Conserved mechanisms**: Identify shared social insect features

---

## ✅ Success Criteria

### Per Species (4 total)
- ✅ Genome downloaded and indexed
- ✅ Metadata retrieved from SRA
- ✅ All FASTQ files downloaded
- ✅ All samples quantified with Kallisto
- ✅ Expression matrices generated (TPM, counts, eff_length)
- ✅ QC visualizations created (6 PDFs)
- ✅ Curate tables generated (7 TSVs)
- ✅ 100% sanity validation passed

### Cross-Species
- ✅ Single-copy genes identified across all 4 species
- ✅ TMM normalization completed
- ✅ Correlation matrix computed
- ✅ Cross-species plots generated
- ✅ Comparative analysis ready for publication

---

## 📝 Key Files

### Execution
- **Script**: `scripts/rna/run_multi_species_amalgkit.py`
- **Log**: `output/amalgkit_multi_species_comprehensive.log`
- **Plan**: `output/MULTI_SPECIES_EXECUTION_PLAN.md`

### Monitoring
- **Monitor script**: `scripts/rna/monitor_amalgkit_progress.sh`
- **Status doc**: This file

### Configs
- `config/amalgkit_pbarbatus.yaml`
- `config/amalgkit_sinvicta.yaml`
- `config/amalgkit_cfloridanus.yaml`
- `config/amalgkit_mpharaonis.yaml`

---

## 🔄 What Happens Next

1. **Individual workflows complete** (60-72 hours)
   - Each species gets full expression data
   - QC plots and curated matrices
   - Sanity validation

2. **Cross-species CSTMM runs** (15-30 min)
   - Identifies ~5,000-10,000 single-copy genes
   - Normalizes expression across species
   - Creates unified expression matrix

3. **Cross-species CSCA runs** (15-30 min)
   - Computes correlations across species
   - Generates comparative plots
   - Identifies co-expression modules

4. **Analysis ready** 
   - Publication-quality figures
   - Normalized cross-species data
   - Comparative gene expression insights

---

**Status**: ✅ Running smoothly  
**Auto-Discovery**: ✅ 4 species found  
**Cross-Species**: ✅ Configured  
**Est. Completion**: November 1, 2025 afternoon

Monitor progress: `tail -f output/amalgkit_multi_species_comprehensive.log`

