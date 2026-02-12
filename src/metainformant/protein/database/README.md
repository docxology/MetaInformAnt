# Protein Database

Remote database clients for UniProt and InterPro protein annotations, providing programmatic access to protein records, domain annotations, and cross-references.

## Contents

| File | Purpose |
|------|---------|
| `uniprot.py` | UniProt REST API client for protein records, FASTA retrieval, and ID mapping |
| `interpro.py` | InterPro API client for domain annotations and functional site queries |

## Key Functions

| Function | Description |
|----------|-------------|
| `fetch_uniprot_record()` | Retrieve complete UniProt record (sequence, annotations, PTMs) |
| `fetch_uniprot_fasta()` | Download protein sequence in FASTA format |
| `search_uniprot_proteins()` | Query UniProt by keyword with pagination |
| `batch_fetch_uniprot_records()` | Bulk retrieval for multiple accessions |
| `map_ids_uniprot()` | Cross-reference ID mapping between databases |
| `validate_uniprot_accession()` | Check accession format validity |
| `fetch_interpro_domains()` | Get InterPro domain annotations for a UniProt protein |
| `batch_fetch_interpro_domains()` | Bulk domain annotation retrieval |
| `search_interpro_entries()` | Search InterPro entries by keyword |
| `get_interpro_go_annotations()` | Fetch GO term mappings for an InterPro entry |

## Usage

```python
from metainformant.protein.database.uniprot import fetch_uniprot_record, search_uniprot_proteins
from metainformant.protein.database.interpro import fetch_interpro_domains

record = fetch_uniprot_record("P12345")
domains = fetch_interpro_domains("P12345")
results = search_uniprot_proteins("kinase AND organism_id:7227", max_results=50)
```
