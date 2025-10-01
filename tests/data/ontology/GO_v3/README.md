# GO v3 Test Data

This directory contains Gene Ontology version 3 test data for ontology parsing and analysis validation.

## Files

### Core Ontology Files
- **`go.obo`**: Complete Gene Ontology definition file in OBO format
- **`go-basic.obo`**: Simplified version for basic parsing tests

### Association Files
- **`gene_association.test`**: Test gene-to-GO term associations
- **`gene_association_subset.test`**: Reduced dataset for performance testing

### Validation Data
- **`expected_parsing_results.json`**: Expected results for parsing validation
- **`hierarchy_test_cases.json`**: Test cases for GO hierarchy traversal

## Usage in Tests

These files are used for:
- OBO format parsing validation
- GO term relationship testing
- Gene annotation processing
- Ontology hierarchy navigation verification

## Data Source

- Gene Ontology Consortium (golang.org)
- Subset created for testing purposes
- Maintains structure and relationships of full GO

## Maintenance

- Update when testing new GO versions
- Ensure compatibility with OBO format 1.4 specification
- Validate parsing results against known correct outputs
- Document any changes in expected behavior
