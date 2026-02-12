# UniProt Database Integration

Utilities for querying the UniProt protein knowledge base, retrieving
records, parsing FASTA headers, fetching annotations, searching proteins,
mapping identifiers, and batch operations.

## Key Concepts

**UniProt** is the most comprehensive protein sequence and functional
annotation resource. METAINFORMANT interfaces with the UniProt REST API
(https://rest.uniprot.org) for programmatic access.

**Accession format**: UniProt accessions follow patterns like `P12345`
(6 characters) or `A0A1234567` (10 characters). The `validate_uniprot_accession`
function validates these formats.

**Record structure**: A fetched UniProt record contains accession, entry name,
protein name, organism, taxonomy ID, sequence, gene name, function description,
subcellular locations, domains, and post-translational modifications.

## Function Reference

### fetch_uniprot_record

```python
def fetch_uniprot_record(uniprot_id: str) -> Dict[str, Any]
```

Fetches a complete UniProt record via the JSON API. Returns a dictionary with
keys: `accession`, `entry_name`, `protein_name`, `organism`, `taxon_id`,
`sequence`, `length`, `gene_name`, `function`, `subcellular_location`,
`domains`, `ptms`. Raises `requests.RequestException` on failure.

### fetch_uniprot_fasta

```python
def fetch_uniprot_fasta(uniprot_id: str) -> Optional[str]
```

Returns the protein sequence in FASTA format, or `None` if the request fails.

### parse_uniprot_fasta_header

```python
def parse_uniprot_fasta_header(header: str) -> Dict[str, str]
```

Parses a UniProt FASTA header line (without the leading `>`) into a
dictionary with keys: `database`, `accession`, `entry_name`, `description`,
`organism`, `gene_name`.

```python
>>> parse_uniprot_fasta_header("sp|P12345|PROT_HUMAN Protein OS=Homo sapiens")
{'database': 'sp', 'accession': 'P12345', 'entry_name': 'PROT_HUMAN', ...}
```

### get_uniprot_annotations

```python
def get_uniprot_annotations(uniprot_id: str) -> List[Dict[str, Any]]
```

Retrieves GO term and keyword annotations for a UniProt entry.

### search_uniprot_proteins

```python
def search_uniprot_proteins(
    query: str,
    max_results: int = 100,
) -> List[Dict[str, Any]]
```

Searches UniProt using a query string (supports the UniProt query syntax).
Returns a list of protein summaries (accession, name, organism, gene, length).
The API caps at 500 results per request.

### get_uniprot_taxonomy_info

```python
def get_uniprot_taxonomy_info(taxon_id: int) -> Optional[Dict[str, Any]]
```

Fetches taxonomy information for an NCBI taxonomy ID. Returns scientific name,
common name, rank, lineage, and parent ID.

### batch_fetch_uniprot_records

```python
def batch_fetch_uniprot_records(
    uniprot_ids: List[str],
) -> Dict[str, Dict[str, Any]]
```

Fetches multiple records sequentially. Returns a mapping from ID to record
(or `None` for failures).

### validate_uniprot_accession

```python
def validate_uniprot_accession(accession: str) -> bool
```

Validates accession format against known UniProt patterns.

```python
>>> validate_uniprot_accession("P12345")
True
>>> validate_uniprot_accession("INVALID")
False
```

### map_ids_uniprot

```python
def map_ids_uniprot(
    protein_ids: List[str],
    source_db: str = "auto",
    target_format: str = "accession",
) -> Dict[str, str]
```

Maps protein identifiers from external databases (Ensembl, RefSeq, PDB,
GeneID) to UniProt accessions using the UniProt ID Mapping API. The
`source_db` parameter can be `"auto"` (detects based on ID prefix),
`"ensembl"`, `"refseq"`, `"pdb"`, or `"geneid"`.

## Usage Example

```python
from metainformant.protein.database.uniprot import (
    fetch_uniprot_record,
    search_uniprot_proteins,
    validate_uniprot_accession,
    map_ids_uniprot,
)

# Fetch a single record
record = fetch_uniprot_record("P53_HUMAN")
print(f"{record['protein_name']} ({record['organism']})")
print(f"Sequence length: {record['length']}")

# Search for kinases in human
results = search_uniprot_proteins("kinase AND organism_id:9606", max_results=10)
for r in results:
    print(f"{r['accession']}: {r['protein_name']}")

# Map Ensembl protein IDs to UniProt
mapping = map_ids_uniprot(["ENSP00000389680"], source_db="ensembl")
```

## Configuration

| Environment Variable | Description               | Default |
|---------------------|---------------------------|---------|
| `PROT_TIMEOUT`      | HTTP request timeout (s)  | 30      |

## Related Modules

- `metainformant.protein.structure.alphafold` -- AlphaFold predictions
- `metainformant.protein.sequence.sequences` -- protein sequence analysis
- `metainformant.protein.database.interpro` -- InterPro domain database
- `metainformant.protein.domains` -- domain detection and classification
