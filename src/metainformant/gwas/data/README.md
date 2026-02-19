# GWAS Data

GWAS data retrieval, download utilities, genome mapping, and sample metadata management.

## Contents

| File | Purpose |
|------|---------|
| `config.py` | GWAS data configuration utilities |
| `download.py` | Reference genome, variant database, and SRA data downloads |
| `genome.py` | Apis mellifera chromosome mapping, normalization, GFF3 parsing |
| `metadata.py` | Sample metadata loading, validation, and population labels |
| `sra_download.py` | SRA sequencing data download and extraction |

## Key Functions

| Function | Module | Description |
|----------|--------|-------------|
| `download_reference_genome()` | `download` | Download reference genome from NCBI (Datasets API or FTP) |
| `download_variant_database()` | `download` | Download dbSNP, 1000 Genomes, or HapMap variants |
| `download_sra_run()` | `download` | Download individual SRA run data |
| `download_sra_project()` | `download` | Download all runs from an SRA project |
| `search_sra_for_organism()` | `download` | Search SRA for organism-specific sequencing data |
| `download_annotation()` | `download` | Download genome annotation (GFF/GTF) files |
| `normalize_chromosome_name()` | `genome` | Normalize chromosome names (NCBI acc, chr, Group, roman) |
| `parse_gff3_genes()` | `genome` | Parse gene locations from GFF3 annotation files |
| `load_sample_metadata()` | `metadata` | Load sample metadata from TSV/CSV |
| `get_population_labels()` | `metadata` | Extract population labels from metadata |

## Usage

```python
from metainformant.gwas.data import download, genome, metadata

# Download reference genome
genome_dir = download.download_reference_genome("GCF_003254395.2", "output/genome/")

# Normalize chromosome names
chrom = genome.normalize_chromosome_name("NC_037638.1")  # -> 1

# Load sample metadata
meta = metadata.load_sample_metadata("samples.tsv")
```
