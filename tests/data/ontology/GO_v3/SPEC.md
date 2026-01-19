# SPEC: Ontology Test Data Scripts (GO_v3)

Scripts for generating and processing Gene Ontology test datasets.

## Scripts

- `1_uniprot_ID_extract.py`: Extracts UniProt IDs from raw sources.
- `2_genetogo.py`: Maps genes to GO terms.
- `3_genetogotoanno.py`: Generates functional annotation files.
- `4_genetogo_summary.py`: Aggregates mapping statistics.

## Standards

- **Internal Consistency**: Output formats must remain stable for the regression test suite.
- **Reference**: Uses the GO_v3 schema for all mappings.
