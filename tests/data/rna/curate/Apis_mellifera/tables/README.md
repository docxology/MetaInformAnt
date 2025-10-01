# Apis mellifera RNA Tables

This directory contains tabular data for Apis mellifera (honey bee) RNA analysis testing.

## Files

### Expression Tables
- **`gene_expression.csv`**: Gene expression matrix for testing
- **`transcript_expression.csv`**: Transcript-level expression data

### Metadata Tables
- **`sample_metadata.csv`**: Sample information and conditions
- **`gene_metadata.csv`**: Gene annotations and descriptions

### Quality Tables
- **`quality_metrics.csv`**: RNA-seq quality metrics
- **`batch_effects.csv`**: Batch effect information

## Usage in Tests

These files are used for:
- Expression matrix processing validation
- Metadata integration testing
- Quality control algorithm verification
- Batch correction testing

## Data Source

- Apis mellifera RNA-seq data from public repositories
- Curated for testing purposes
- Maintains biological relevance while being suitable for unit testing

## Maintenance

- Update when testing new RNA analysis features
- Ensure compatibility with current data formats
- Keep file sizes minimal while maintaining biological relevance
- Document any changes to expected results
