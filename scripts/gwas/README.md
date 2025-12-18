# GWAS Analysis Scripts

Genome-wide association study workflows with comprehensive visualization and real data integration.

## Directory Structure

```
scripts/gwas/
├── run_genome_scale_gwas.py       # Complete SRA→variants→GWAS pipeline ⭐
├── run_pbarbatus_gwas.py          # Pogonomyrmex barbatus GWAS analysis
├── run_pbarbatus_analysis.py      # P. barbatus complete workflow
├── run_analysis.py                # General GWAS analysis orchestrator
├── download_genome_scale_data.sh  # Real A. mellifera data acquisition
├── download_real_honeybee_variants.py # Honeybee variants guide
├── query_bioproject_metadata.py   # NCBI BioProject metadata queries
├── generate_phenotypes.py         # Phenotype data generation
├── generate_synthetic_variants.py # Synthetic variant data generation
├── enhanced_visualizations.py     # Advanced plotting utilities
├── generate_missing_plots.py      # Plot regeneration utility
└── README.md                      # This file
```

## Genome-Scale GWAS Pipeline (`run_genome_scale_gwas.py`)

Complete end-to-end GWAS workflow from SRA data download through association testing and visualization.

**Features:**
- SRA data download and processing
- BWA alignment and variant calling (bcftools)
- Comprehensive statistical analysis
- 32 visualization types (Manhattan, Q-Q, PCA, etc.)
- Real data integration (no synthetic data)

**Usage:**
```bash
# Complete GWAS pipeline
python3 scripts/gwas/run_genome_scale_gwas.py --config config/gwas/gwas_amellifera.yaml

# Download only
python3 scripts/gwas/run_genome_scale_gwas.py --download-only

# Skip download, use existing data
python3 scripts/gwas/run_genome_scale_gwas.py --skip-download --config config/gwas/gwas_amellifera.yaml

# Custom threads
python3 scripts/gwas/run_genome_scale_gwas.py --threads 16 --config config/gwas/gwas_amellifera.yaml
```

**Dependencies:** `fasterq-dump`, `bwa`, `samtools`, `bcftools`

## Real Data Acquisition

### Honeybee Genome-Scale Data (`download_genome_scale_data.sh`)

Batch download script for real *Apis mellifera* behavioral phenotype data.

**Features:**
- Downloads from BioProject PRJNA292680
- Scout vs Recruit behavioral phenotypes
- Progress reporting and validation
- Automatic FASTQ conversion

**Usage:**
```bash
# Download real honeybee data
bash scripts/gwas/download_genome_scale_data.sh
```

**Output:** `data/raw/sra/` with FASTQ files

### Honeybee Variants Guide (`download_real_honeybee_variants.py`)

Comprehensive guide for acquiring real honeybee variant datasets.

**Features:**
- Lists key BioProjects (PRJNA292680, PRJNA392242, PRJNA13343)
- SRA Toolkit usage instructions
- Complete workflow documentation
- Variant calling pipeline guidance

**Usage:**
```bash
# Show available honeybee datasets
python3 scripts/gwas/download_real_honeybee_variants.py
```

## Species-Specific Analysis

### Pogonomyrmex barbatus GWAS (`run_pbarbatus_gwas.py`)

Specialized GWAS analysis for harvester ant behavioral phenotypes.

**Features:**
- Queen vs Worker phenotype association
- Multiple behavioral traits analysis
- Comprehensive statistical validation
- Advanced visualization suite

**Usage:**
```bash
# Complete P. barbatus GWAS
python3 scripts/gwas/run_pbarbatus_gwas.py --config config/gwas/gwas_pbarbatus.yaml
```

## Data Generation and Utilities

### Phenotype Generation (`generate_phenotypes.py`)

Generate synthetic phenotype data for GWAS testing and validation.

**Usage:**
```bash
# Generate phenotypes for testing
python3 scripts/gwas/generate_phenotypes.py --n-samples 1000 --n-traits 5
```

### Synthetic Variants (`generate_synthetic_variants.py`)

Generate synthetic genetic variant data for method development.

**Usage:**
```bash
# Generate synthetic variants
python3 scripts/gwas/generate_synthetic_variants.py --n-variants 10000 --n-samples 500
```

### Metadata Queries (`query_bioproject_metadata.py`)

Query NCBI BioProject metadata for sample information and run counts.

**Usage:**
```bash
# Query BioProject metadata
python3 scripts/gwas/query_bioproject_metadata.py --bioproject PRJNA292680
```

## Visualization and Analysis

### Enhanced Visualizations (`enhanced_visualizations.py`)

Advanced plotting utilities for GWAS results.

**Features:**
- Manhattan plots with gene annotations
- Q-Q plots with inflation factor
- PCA plots with population structure
- LD heatmaps and haplotype blocks
- Power analysis plots

### Plot Regeneration (`generate_missing_plots.py`)

Utility to regenerate missing visualization files.

**Usage:**
```bash
# Regenerate missing plots
python3 scripts/gwas/generate_missing_plots.py --results-dir output/gwas/results
```

## General GWAS Analysis (`run_analysis.py`)

General-purpose GWAS analysis orchestrator for existing genotype/phenotype data.

**Usage:**
```bash
# GWAS analysis on existing data
python3 scripts/gwas/run_analysis.py --genotypes variants.vcf --phenotypes traits.tsv --output output/gwas/analysis
```

## Output Structure

```
output/gwas/
├── [species]_[date]_gwas/
│   ├── raw_data/                 # Downloaded SRA/FASTQ files
│   ├── alignments/               # BAM files from BWA alignment
│   ├── variants/                 # VCF files from bcftools
│   ├── gwas_results/             # Association test results
│   │   ├── manhattan_data.json
│   │   ├── qq_data.json
│   │   └── association_stats.json
│   ├── plots/                    # All visualizations
│   │   ├── manhattan_plot.png
│   │   ├── qq_plot.png
│   │   ├── pca_plot.png
│   │   ├── ld_heatmap.png
│   │   └── phenotype_distributions.png
│   └── logs/                     # Processing logs
└── metadata/                     # BioProject metadata queries
```

## Key Features

✅ **Real Data Integration**: Uses actual NCBI SRA datasets (PRJNA292680, etc.)
✅ **32 Visualization Types**: Comprehensive plotting for GWAS analysis
✅ **Multi-Species Support**: Configurable for different species
✅ **Complete Pipeline**: SRA download → alignment → variants → association → plots
✅ **Step-by-Step Execution**: Can skip completed steps, resume from checkpoints
✅ **Dependency Checking**: Validates required bioinformatics tools
✅ **Comprehensive Error Handling**: Clear error messages and recovery guidance

## Integration

Integrates with:
- **metainformant.gwas**: Core GWAS analysis functionality
- **NCBI SRA Toolkit**: Data download and processing
- **BWA/samtools/bcftools**: Alignment and variant calling
- **Core utilities**: I/O, logging, configuration management
- **Visualization**: Plot generation and customization

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

## Related Documentation

- [GWAS Analysis Documentation](../../docs/gwas/README.md)
- [Bioinformatics Setup](../../docs/gwas/setup.md)
- [Real Data Sources](../../docs/gwas/data_sources.md)
- [METAINFORMANT CLI](../../docs/cli.md)

