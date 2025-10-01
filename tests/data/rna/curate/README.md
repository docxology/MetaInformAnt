# RNA Curation Test Data

This directory contains test data for RNA data curation and quality control functionality.

## Files

### Curated Datasets
- **`Apis_mellifera/`**: Species-specific curated RNA data
  - **`metadata.json`**: Curated sample metadata
  - **`expression_data.csv`**: Curated expression matrix
  - **`quality_report.json`**: Curation quality metrics

### Curation Test Cases
- **`curation_test_cases.json`**: Structured test cases for curation algorithms
- **`expected_curation_results.json`**: Expected outputs for validation

## Usage in Tests

These files are used throughout the RNA curation module tests:
- Curation algorithm validation
- Quality control metric testing
- Batch effect correction verification
- Data filtering and standardization testing

## Data Sources

- Curated RNA datasets from public repositories
- Synthetic data generated for edge case testing
- Quality-controlled examples for validation

## Maintenance

- Update curated data when testing new curation methods
- Ensure compatibility with current RNA-seq standards
- Keep file sizes minimal while maintaining test coverage
- Document any changes to curation expectations
