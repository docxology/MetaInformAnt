# RNA Test Data

This directory contains test data for RNA transcriptomic analysis functionality.

## Files

### Expression Data
- **`expression_matrix.csv`**: Test gene expression matrix
- **`sample_metadata.csv`**: Sample metadata for testing

### Annotation Data
- **`transcript_annotations.gtf`**: Test transcript annotation file
- **`gene_annotations.gff`**: Gene annotation test data

### Workflow Test Data
- **`test_workflow_config.yaml`**: Configuration for workflow testing
- **`expected_results.json`**: Expected workflow outputs

## Usage in Tests

These files are used throughout the RNA module tests:
- Expression matrix processing validation
- Annotation file parsing testing
- Workflow configuration and execution
- Integration testing with external tools

## Data Sources

- Synthetic expression data generated for testing
- Public domain RNA datasets where applicable
- Mock workflow results for validation

## Maintenance

- Update expression data when testing new formats
- Ensure compatibility with current annotation standards
- Keep file sizes minimal while maintaining test coverage
- Document any changes to expected processing behavior
