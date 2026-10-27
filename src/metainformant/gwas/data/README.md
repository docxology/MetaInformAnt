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
| `download_variant_database()` | `download` | Download dbSNP or 1000 Genomes variants (`db_name` `"dbsnp"` or `"1000genomes"`; any other name raises `ValueError`) |
| `download_variant_data()` | `download` | Download the association summary table of a GWAS Catalog study. Only `data_type="summary_stats"` is supported: `"vcf"` (the default) and `"both"` raise `NotImplementedError`, and an unknown type raises `ValueError` |
| `download_sra_run()` | `download` | Download individual SRA run data |
| `download_sra_project()` | `download` | Download all runs from an SRA project |
| `search_sra_for_organism()` | `download` | Search SRA for organism-specific sequencing data |
| `download_annotation()` | `download` | Download genome annotation (GFF/GTF) files |
| `normalize_chromosome_name()` | `genome` | Normalize chromosome names (NCBI acc, chr, Group, roman) |
| `parse_gff3_genes()` | `genome` | Parse gene locations from GFF3 annotation files |
| `load_sample_metadata()` | `metadata` | Load sample metadata from TSV/CSV |
| `get_population_labels()` | `metadata` | Extract population labels from metadata |

NCBI-backed SRA helpers require an explicit ``email=`` argument or
``NCBI_EMAIL``. Pass ``allow_anonymous=True`` only for an intentionally
anonymous request; no helper fabricates or persists a contact address. Genome
download functions raise when both acquisition strategies fail rather than
returning an empty directory.

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
