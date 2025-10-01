# DNA Test Data

This directory contains test data for DNA sequence analysis functionality, including various sequence formats, genomic data, and analysis fixtures.

## Files

### Sequence Data
- **`toy.fasta`**: Small test sequences for basic functionality
- **`sample_sequences.fasta`**: Representative DNA sequences for testing
- **`test_genome.fasta`**: Mock genome sequence for integration testing

### Genomic Data
- **`reference_annotation.gff`**: Gene annotation file for testing
- **`test_variants.vcf`**: Variant call format file for testing
- **`chromosome_fragments.fasta`**: Chromosome segment data

### Alignment Data
- **`test_alignment.fasta`**: Multiple sequence alignment
- **`pairwise_test.fasta`**: Pairwise alignment test data

## Usage in Tests

These files are used throughout the DNA module tests:
- Sequence I/O validation
- Alignment algorithm testing
- Phylogenetic tree construction
- Population genetics calculations
- Genomic coordinate operations

## Data Sources

- Synthetic sequences generated for testing purposes
- Public domain reference data where applicable
- Minimal representative datasets for comprehensive coverage

## Maintenance

- Update sequences when testing new formats
- Ensure compatibility with current Biopython versions
- Keep file sizes minimal while maintaining test coverage
- Document any changes to test data expectations
