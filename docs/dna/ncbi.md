### DNA: NCBI Integration

The `metainformant.dna.external` package provides three modules for interacting with
NCBI services: `entrez` (Entrez E-utilities client), `ncbi` (higher-level NCBI API
client), and `genomes` (genome assembly retrieval via the NCBI Datasets API).

All network-dependent functions require internet access. NCBI enforces rate limits:
3 requests/second without an API key, 10 requests/second with one.

---

## Module Overview

| Module     | Class / Functions                                         | Primary Use Case                   |
|------------|-----------------------------------------------------------|------------------------------------|
| `entrez`   | `EntrezClient`, `search_genbank`, `fetch_fasta_sequence`  | Entrez E-utilities (efetch, esearch, esummary, elink) |
| `ncbi`     | `NCBIClient`, `search_nucleotide`, `fetch_sequence`       | Nucleotide/protein search, taxonomy lookups |
| `genomes`  | `download_genome_package`, `get_genome_metadata`          | Genome assembly downloads via Datasets API  |

---

## Entrez Client (`dna.external.entrez`)

### `EntrezClient`

A session-based client wrapping the NCBI Entrez E-utilities REST API.

```python
from metainformant.dna.external.entrez import EntrezClient

client = EntrezClient(email="user@example.com", api_key="your_ncbi_key")
```

**Constructor Parameters:**

| Parameter | Type          | Description                                 |
|-----------|---------------|---------------------------------------------|
| `email`   | `str | None`  | Required by NCBI for identification         |
| `api_key` | `str | None`  | NCBI API key for higher rate limits (10/s)  |

**Methods:**

| Method      | Description                                          | Returns               |
|-------------|------------------------------------------------------|-----------------------|
| `search`    | Run esearch against any NCBI database                | `Dict[str, Any]`      |
| `fetch`     | Retrieve records via efetch (FASTA, GenBank, etc.)   | `str`                 |
| `summary`   | Get document summaries via esummary                  | `Dict[str, Any]`      |
| `link`      | Find linked records across databases via elink       | `Dict[str, Any]`      |

### Convenience Functions

```python
from metainformant.dna.external.entrez import (
    search_genbank,
    fetch_genbank_record,
    fetch_fasta_sequence,
    get_sequence_features,
    search_protein_records,
    link_nucleotide_to_protein,
)
```

#### `search_genbank`

Search the GenBank nucleotide database and return structured record summaries.

```python
records = search_genbank(
    "Apis mellifera[Organism] AND mitochondrion[Title]",
    max_results=5,
    email="user@example.com",
)
# Returns list of dicts with keys:
# id, accession, title, organism, length, moltype, created, updated
```

#### `fetch_fasta_sequence`

Download a sequence in FASTA format and return the raw nucleotide string
(header stripped).

```python
sequence = fetch_fasta_sequence("NC_001422.1", email="user@example.com")
# Returns the phiX174 genome as a plain string
```

#### `fetch_genbank_record`

Fetch a full GenBank flat-file record.

```python
record_text = fetch_genbank_record("NC_001422.1", email="user@example.com")
# Returns the complete GenBank-format text
```

#### `get_sequence_features`

Parse features from a GenBank record (simplified parser extracting feature type
and location).

```python
features = get_sequence_features("NC_001422.1", email="user@example.com")
# [{"type": "gene", "location": "...", "qualifiers": {}}, ...]
```

#### `link_nucleotide_to_protein`

Find protein records linked to a nucleotide accession via elink.

```python
protein_ids = link_nucleotide_to_protein("NC_001422.1", email="user@example.com")
```

---

## NCBI Client (`dna.external.ncbi`)

### `NCBIClient`

A higher-level client providing typed search across nucleotide, protein, and
taxonomy databases.

```python
from metainformant.dna.external.ncbi import NCBIClient

client = NCBIClient(email="user@example.com", api_key="your_key")
```

**Methods:**

| Method               | Description                                    | Returns                |
|----------------------|------------------------------------------------|------------------------|
| `search_nucleotide`  | Search nucleotide database with summaries      | `List[Dict[str, Any]]` |
| `fetch_sequence`     | Fetch DNA sequence by accession (FASTA)        | `str | None`           |
| `search_protein`     | Search protein database with summaries         | `List[Dict[str, Any]]` |
| `get_taxonomy_info`  | Retrieve taxonomy by NCBI taxonomy ID          | `Dict[str, Any] | None`|

### Convenience Functions

```python
from metainformant.dna.external.ncbi import (
    search_nucleotide,
    fetch_sequence,
    search_protein,
    get_taxonomy_info,
)

# Search nucleotide database
results = search_nucleotide("BRCA1 Homo sapiens", max_results=10)

# Fetch a specific sequence
seq = fetch_sequence("NM_007294.4")

# Taxonomy lookup
tax_info = get_taxonomy_info(7460)  # Apis mellifera
# {
#     "tax_id": 7460,
#     "scientific_name": "Apis mellifera",
#     "common_name": "honey bee",
#     "rank": "species",
#     "division": "Invertebrates",
#     "lineage": "...",
# }
```

---

## Genome Assembly Retrieval (`dna.external.genomes`)

### `download_genome_package`

Download a genome assembly via the NCBI Datasets v2 API. Validates the accession
before making network requests. Downloads a ZIP archive and extracts it.

```python
from metainformant.dna.external.genomes import download_genome_package

# Download human reference genome (large file)
package_dir = download_genome_package("GCF_000001405.39", "output/genomes/")
```

### `download_genome_package_best_effort`

Tries the Datasets API first, then falls back to FTP download if provided.

### `get_genome_metadata`

Retrieve genome metadata (organism, assembly level, genome size, GC content, etc.)
without downloading the full assembly.

```python
from metainformant.dna.external.genomes import get_genome_metadata

metadata = get_genome_metadata("GCF_000001405.39")
# {
#     "accession": "GCF_000001405.39",
#     "organism_name": "Homo sapiens",
#     "assembly_level": "Chromosome",
#     "genome_size": 3088286401,
#     "gc_percent": 40.9,
#     "contig_count": ...,
#     "release_date": "...",
# }
```

### `list_genome_assemblies`

Search for all genome assemblies for a given organism.

```python
from metainformant.dna.external.genomes import list_genome_assemblies

assemblies = list_genome_assemblies("Drosophila melanogaster", max_results=5)
# List of dicts with accession, organism_name, assembly_level, genome_size, release_date
```

### `download_reference_genome`

Convenience function that searches for the primary reference assembly and downloads it.

```python
from metainformant.dna.external.genomes import download_reference_genome

# genome_path = download_reference_genome("Apis mellifera", "output/genomes/")
```

### `get_chromosome_lengths`

Retrieve chromosome/contig lengths for an assembly. Requires BioPython.

```python
from metainformant.dna.external.genomes import get_chromosome_lengths

# lengths = get_chromosome_lengths("GCF_000001405.39")
# {"chr1": 248956422, "chr2": 242193529, ...}
```

### `validate_genome_files`

Check that a downloaded genome directory contains expected files (FASTA, GFF).

```python
from metainformant.dna.external.genomes import validate_genome_files

results = validate_genome_files("output/genomes/GCF_000001405.39/")
# {"valid": True, "fasta_files": 1, "gff_files": 1, "total_files": 5, "issues": []}
```

---

## Configuration and Rate Limits

NCBI requires an email address for all Entrez requests. An API key is optional but
recommended for production use (raises rate limit from 3 to 10 requests/second).

Set credentials via constructor parameters or environment variables:

```python
import os
os.environ["NCBI_EMAIL"] = "user@example.com"
os.environ["NCBI_API_KEY"] = "your_key_here"
```

Without an API key, the client automatically sleeps 0.4 seconds between requests
to respect the 3 requests/second limit.

---

## See Also

- **[Accessions](accessions.md)** -- Accession format validation
- **[DNA Module Overview](../dna/)** -- Full DNA module documentation
- `metainformant.dna.external.entrez` -- Entrez source module
- `metainformant.dna.external.ncbi` -- NCBI client source module
- `metainformant.dna.external.genomes` -- Genome retrieval source module
