# Ontology Test Data

This directory contains test data for ontology processing and functional annotation functionality.

## Files

### Gene Ontology Data
- **`GO_v3/`**: Gene Ontology version 3 test files
  - **`go.obo`**: GO ontology definition file
  - **`gene_association.test`**: Gene to GO term associations
  - **`test_terms.txt`**: Sample GO terms for testing

### Test Cases
- **`ontology_test_cases.json`**: Structured test cases for ontology parsing
- **`enrichment_test_data.json`**: Data for enrichment analysis testing

## Usage in Tests

These files are used throughout the ontology module tests:
- OBO format parsing validation
- GO term traversal and hierarchy testing
- Enrichment analysis algorithm verification
- Semantic similarity calculation testing

## Data Sources

- Official Gene Ontology Consortium data (subset)
- Synthetic test cases for edge condition testing
- Public domain ontology examples

## Maintenance

- Update GO data when testing new ontology versions
- Ensure compatibility with current OBO format specifications
- Keep file sizes minimal while maintaining comprehensive coverage
- Document any changes to expected parsing behavior
