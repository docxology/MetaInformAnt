# GWAS Analysis Scripts

Genome-wide association study workflows with comprehensive visualization and real data integration.

## Directory Structure

```
scripts/gwas/
├── preparation/           # Data download and generation
│   ├── download_genome_scale_data.sh
│   ├── download_honeybee_variants.py
│   ├── generate_phenotypes.py
│   ├── generate_synthetic_variants.py
│   └── query_bioproject_metadata.py
├── qc/                    # Quality Control
│   └── run_qc.py
├── structure/             # Population Structure
│   ├── run_pca.py
│   └── run_kinship.py
├── association/           # Association Testing
│   └── run_association.py
├── visualization/         # Visualization
│   ├── visualizations.py
│   └── generate_missing_plots.py
├── pipelines/             # End-to-end pipelines
│   ├── run_genome_scale_gwas.py
│   ├── run_pbarbatus_gwas.py
│   ├── run_pbarbatus_analysis.py
│   └── run_analysis.py
└── README.md
```

## Pipelines (`scripts/gwas/pipelines/`)

### Genome-Scale GWAS Pipeline (`run_genome_scale_gwas.py`)
Complete end-to-end GWAS workflow from SRA data download through association testing and visualization.

**Usage:**
```bash
python3 scripts/gwas/pipelines/run_genome_scale_gwas.py --config config/gwas/gwas_amellifera.yaml
```

### P. barbatus GWAS (`run_pbarbatus_gwas.py`)
Specialized GWAS analysis for harvester ant behavioral phenotypes.

**Usage:**
```bash
python3 scripts/gwas/pipelines/run_pbarbatus_gwas.py --config config/gwas/gwas_pbarbatus.yaml
```

## Modules

### Preparation (`scripts/gwas/preparation/`)
Scripts for acquiring real data from NCBI SRA or generating synthetic datasets for testing.

### Quality Control (`scripts/gwas/qc/`)
Standalone scripts for variant quality control.

### Structure (`scripts/gwas/structure/`)
Scripts for analyzing population structure (PCA, Kinship).

### Association (`scripts/gwas/association/`)
Scripts for running association tests (Linear/Logistic regression).

### Visualization (`scripts/gwas/visualization/`)
Plotting utilities and regeneration scripts.

## Dependencies

**Required bioinformatics tools:**
- `fasterq-dump` (SRA Toolkit)
- `bwa` (alignment)
- `samtools` (BAM processing)
- `bcftools` (variant calling)

**Python packages:**
- `metainformant.gwas`
- `pandas`, `numpy`
- `matplotlib`, `seaborn`
- `scikit-learn`
- `statsmodels`
