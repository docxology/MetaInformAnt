# Protein Test Data

This directory contains test data for protein sequence and structure analysis functionality.

## Files

### Sequence Data
- **`sample_proteins.fasta`**: Representative protein sequences for testing
- **`test_sequences.fasta`**: Basic protein sequences for unit tests

### Structure Data
- **`test_structure.pdb`**: Sample PDB structure file for testing
- **`reference_structure.pdb`**: Reference structure for validation

### Database Integration
- **`uniprot_test_data.json`**: Mock UniProt API response data
- **`interpro_test_data.json`**: Mock InterPro domain data
- **`alphafold_test_data.json`**: Mock AlphaFold structure data

## Usage in Tests

These files are used throughout the protein module tests:
- Sequence I/O validation
- Structure parsing and analysis
- Database integration testing
- Format compatibility verification

## Data Sources

- Synthetic protein sequences generated for testing
- Public domain structure files where applicable
- Mock API responses for external service testing

## Maintenance

- Update sequences when testing new protein formats
- Ensure compatibility with current Biopython versions
- Keep file sizes minimal while maintaining test coverage
- Document any changes to test data expectations
