# SPEC: DNA External Scripts

Wrappers and integration logic for external bioinformatics tools used in DNA analysis.

## Key Integrations

- **NCBI/Entrez**: Scripts for downloading and validating genome packages.
- **SRA Toolkit**: Handling fastq-dump and fasterq-dump logic (integrated with failover patterns).

## Standards

- **Verification**: Always verify the availability of the external binary before execution.
- **Output**: Standardize external tool output into the MetaInformAnt `output/` structure.
