# GWAS Configuration for Apis mellifera

This document describes the GWAS configuration setup for *Apis mellifera* (Western Honey Bee).

## Configuration File

**Location**: `config/gwas/gwas_amellifera.yaml`

## Species Information

- **Species**: *Apis mellifera* (Western Honey Bee)
- **NCBI Taxonomy ID**: 7460
- **Reference Assembly**: GCF_003254395.2 (Amel_HAv3.1)
- **FTP URL**: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/
- **Genome Size**: ~260 Mbp (16 chromosomes)

## Configuration Overview

The configuration file is optimized for honeybee GWAS with extensive genomic datasets and strong population structure.

### 1. Reference Genome

```yaml
genome:
  accession: GCF_003254395.2
  dest_dir: output/gwas/amellifera/genome
  include:
    - genome         # Genomic sequences (FASTA)
    - gff3          # Gene annotations
  ftp_url: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/
```

### 2. Variant Data Sources

Honeybees have extensive genomic datasets available. The configuration supports three options:

#### Option 1: Pre-called VCF Files
Large-scale honeybee projects often provide pre-called VCF files:
- British and Irish *A. mellifera mellifera* whole genome sequences
- European population datasets
- African lineage datasets

#### Option 2: Download from Databases
- NCBI dbSNP (if variants available)
- Research consortium datasets (i5K, Honey Bee Genomics Initiative)

#### Option 3: Variant Calling
Many honeybee whole-genome sequencing datasets are available on SRA:
- Search terms: `"Apis mellifera"[Organism] AND "whole genome"[Strategy]`
- Population-specific searches by subspecies or trait
- Consider batch calling by population/region for large cohorts

### 3. Quality Control

Optimized parameters for honeybee genomes:
- MAF threshold: 0.01 (appropriate for population-level analysis)
- Max missing: 0.05 (honeybee datasets often have high coverage)
- Min call rate: 0.95
- HWE p-value: 1e-6 (accounts for population structure)
- Exclude indels: true
- Min quality: 30.0

### 4. Population Structure

**Critical for honeybee GWAS**: Strong population structure due to:
- Multiple subspecies (5 major lineages)
- Geographic isolation
- Breeding practices

Configuration:
- PCA computation: Enabled
- **20 principal components** (increased from default 10 due to subspecies diversity)
- Kinship matrix: Enabled (VanRaden method)
- Important covariates: Population/subspecies, sampling location

### 5. Association Testing

- Model: Linear regression (configurable to logistic for binary traits)
- Trait: Configurable (common traits: Varroa resistance, hygienic behavior, honey production)
- Min sample size: 50-100 samples recommended
- Covariates: Include top PCs (PC1-PC10), population, sampling location

### 6. Multiple Testing Correction

- Method: Bonferroni
- Alpha: 0.05
- Genome-wide significance: ~5e-8 (similar to human GWAS, 16 chromosomes)

## Common Honeybee Traits

The configuration supports analysis of various honeybee traits:

### Disease Resistance
- Varroa destructor resistance
- Nosema resistance
- American foulbrood resistance
- European foulbrood resistance

### Behavioral Traits
- Hygienic behavior
- Defensive behavior
- Foraging behavior
- Swarming tendency

### Productivity Traits
- Honey production
- Colony strength
- Brood production
- Pollen collection

### Environmental Adaptation
- Overwintering success
- Temperature tolerance
- Drought resistance
- Pollinator efficiency

## Usage

### Configuration Validation

```bash
python -m metainformant gwas run --config config/gwas/gwas_amellifera.yaml --check
```

### Full Workflow Execution

```bash
python -m metainformant gwas run --config config/gwas/gwas_amellifera.yaml
```

### Python API

```python
from metainformant.gwas import load_gwas_config, execute_gwas_workflow

# Load configuration
cfg = load_gwas_config("config/gwas/gwas_amellifera.yaml")

# Validate configuration
result = execute_gwas_workflow(cfg, check=True)

# Execute workflow
result = execute_gwas_workflow(cfg, check=False)
```

## Required Data Files

Before running the full workflow, ensure you have:

1. **Phenotype file** (`data/phenotypes/amellifera/phenotypes.tsv`):
   - TSV format with `sample_id` and trait columns
   - Example traits: `varroa_resistance`, `honey_yield`, `hygienic_behavior`
   - Example:
     ```
     sample_id	varroa_resistance	honey_yield
     BEE001	0.85	45.2
     BEE002	0.92	52.1
     ```

2. **Variant data** (one of):
   - VCF file(s) specified in `variants.vcf_files`
   - Or BAM/CRAM files for variant calling in `variants.calling.bam_files`
   - Or configure download from public databases

3. **Covariates file** (optional but recommended):
   - TSV format with `sample_id`, population, subspecies, sampling location, PCs
   - Example:
     ```
     sample_id	population	subspecies	sampling_location	PC1	PC2
     BEE001	European	Apis_mellifera_mellifera	UK	0.0234	-0.0123
     ```

## Finding Honeybee Genomic Datasets

### NCBI SRA Search Strategies

1. **Whole genome sequences**:
   ```
   "Apis mellifera"[Organism] AND "whole genome"[Strategy] AND Illumina[Platform]
   ```

2. **Population-specific searches**:
   ```
   "Apis mellifera"[Organism] AND "Apis mellifera mellifera"[Title]
   "Apis mellifera"[Organism] AND "African"[Title]
   ```

3. **Trait-specific searches**:
   ```
   "Apis mellifera"[Organism] AND "Varroa"[Title]
   "Apis mellifera"[Organism] AND "resistance"[Title]
   ```

### Public Datasets

- **Zenodo**: British and Irish *A. mellifera mellifera* whole genome sequences
- **ENA**: European Nucleotide Archive
- **GSA**: Genome Sequence Archive (China)
- **Research consortia**: i5K, Honey Bee Genomics Initiative publications

## Workflow Steps

1. **Genome Download**: Downloads A. mellifera Amel_HAv3.1 reference genome from NCBI
2. **Variant Acquisition**: 
   - Uses pre-existing VCF files, OR
   - Downloads from databases, OR
   - Calls variants from BAM/CRAM files
3. **Quality Control**: Filters variants based on QC parameters optimized for honeybees
4. **Population Structure**: Computes PCA (20 components) and kinship matrices
5. **Association Testing**: Performs GWAS with population structure correction
6. **Multiple Testing Correction**: Applies Bonferroni correction
7. **Visualization**: Generates Manhattan plots, Q-Q plots, regional plots
8. **Results Export**: Writes results to output directories

## Output Directories

- **Work directory**: `output/gwas/amellifera/work`
- **Genome**: `output/gwas/amellifera/genome`
- **Variants**: `output/gwas/amellifera/variants`
- **Results**: `output/gwas/amellifera/results`
- **Plots**: `output/gwas/amellifera/plots`
- **Logs**: `output/gwas/amellifera/logs`

## Special Considerations for Honeybees

### Population Structure
- Honeybees show strong population structure across subspecies
- **Always include PCs as covariates** (recommended: PC1-PC10)
- Consider population/subspecies as covariates

### Sample Size
- Minimum 50-100 samples per variant recommended
- Large cohorts (500+ samples) provide better power for rare variants

### Chromosome Organization
- 16 chromosomes in A. mellifera
- Genome-wide significance threshold: ~5e-8 (0.05 / 1,000,000 tests)

### Colony Structure
- Account for queen-worker relationships in kinship matrices
- Consider colony-level covariates when available
- Inbreeding coefficients may be relevant for some populations

## Testing

Run configuration tests:

```bash
pytest tests/test_gwas_config_amellifera.py -v
```

## Related Documentation

- [GWAS Configuration Template](../../config/gwas/gwas_template.yaml)
- [GWAS Workflow Documentation](./workflow.md)
- [GWAS Configuration Guide](./config.md)

## References

1. Wallberg et al. (2014). A worldwide survey of genome sequence variation provides insight into the evolutionary history of the honeybee Apis mellifera. *Nature Genetics*.
2. The Honey Bee Genome Sequencing Consortium (2006). Insights into social insects from the genome of the honeybee Apis mellifera. *Nature*.
3. Cridland et al. (2017). The population genetics of structural variants in honeybees. *G3: Genes, Genomes, Genetics*.






