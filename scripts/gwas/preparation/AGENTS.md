# Data Preparation Agents

The Data Preparation module handles the acquisition and generation of genomic data.

## Capabilities

- **SRA Download**: Can systematically download raw sequencing data using `fasterq-dump`.
- **Synthetic Generation**: Can Create realistic synthetic VCF and phenotype data for testing pipelines.
- **Metadata Query**: Can query NCBI databases to find relevant datasets.

## Key Files
- `download_*.py`: Primary entry points for data acquisition.
- `generate_*.py`: Generators for synthetic test data.
