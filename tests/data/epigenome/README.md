# Epigenome Test Data

This directory contains test data for epigenetic modification analysis functionality.

## Files

### Methylation Data
- **`methylation_calls.txt`**: Test methylation call data
- **`differential_methylation.txt`**: Expected differential methylation results

### Chromatin Data
- **`atac_peaks.bed`**: ATAC-seq peak data for testing
- **`chip_seq_peaks.bed`**: ChIP-seq peak data for testing

## Usage in Tests

These files are used throughout the epigenome module tests:
- Methylation calling algorithm validation
- Peak calling and annotation testing
- Genomic coordinate processing
- Statistical analysis verification

## Data Sources

- Synthetic epigenomic data generated for testing
- Public domain epigenomic datasets where applicable
- Mock processing results for validation

## Maintenance

- Update data when testing new epigenomic formats
- Ensure compatibility with current file format standards
- Keep file sizes minimal while maintaining test coverage
- Document any changes to expected processing behavior
