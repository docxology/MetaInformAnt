### DNA: NCBI Integration

The `metainformant.dna.external` package provides four modules for interacting with
NCBI services: `entrez` (Entrez E-utilities client), `ncbi` (higher-level NCBI API
client), `ncbi_datasets` (typed Datasets v2 requests, pagination and downloads), and
`genomes` (genome assembly retrieval built on that Datasets client).

All network-dependent functions require internet access. NCBI enforces rate limits:
E-utilities and Datasets are separate services: E-utilities clients pace requests
at 0.4 seconds without a key and 0.1 seconds with one. The Datasets client sends
its API key in a header; it does not automatically pace requests or retry 429s.

---

## Module Overview

| Module     | Class / Functions                                         | Primary Use Case                   |
|------------|-----------------------------------------------------------|------------------------------------|
| `entrez`   | `EntrezClient`, `search_genbank`, `fetch_fasta_sequence`  | Entrez E-utilities (efetch, esearch, esummary, elink) |
| `ncbi`     | `NCBIClient`, `search_nucleotide`, `fetch_sequence`       | Nucleotide/protein search, taxonomy lookups |
| `ncbi_datasets` | `NCBIDatasetsClient`, `NCBIDatasetsError`, `AssemblyNotFoundError` | Datasets v2 transport, report pagination and streamed packages |
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

HTTP/transport failures, malformed JSON and recognized E-utilities error bodies
raise `NCBIRequestError`, a `requests.RequestException` subclass exported by
`entrez` and `ncbi`. `EntrezClient.search` validates integer `max_results` in
1..10,000. Empty result lists represent no matches, not failed requests.

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

HTTP/transport failures, malformed JSON and recognized E-utilities error bodies
raise `NCBIRequestError`, a `requests.RequestException` subclass exported by
`entrez` and `ncbi`. `EntrezClient.search` validates integer `max_results` in
1..10,000. Empty result lists represent no matches, not failed requests.

The high-level search methods validate integer `max_results` in 1..1,000.
`fetch_sequence` returns one FASTA record without the header; a missing body
returns `None`, while multiple records and non-FASTA error bodies raise.
Unreported sequence lengths remain `None`.

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
before making network requests. Assembly retrieval requires a versioned
accession such as `GCF_003254395.2`. The ZIP is streamed to disk and extracted;
unsafe member paths are refused, but there is no uncompressed-size cap.
`include` accepts `genome`, `gff3`, `gtf`, `gbff`, `rna`, `protein`, `cds`,
`seq-report`, or supported API enum names, sent as repeated
`include_annotation_type` parameters. Omitting it requests the default genomic
FASTA. The assembly report accompanies the package.

```python
from metainformant.dna.external.genomes import download_genome_package

# Download human reference genome (large file)
package_dir = download_genome_package("GCF_000001405.39", "output/genomes/")
```

### `download_genome_package_best_effort`

Tries the Datasets API first, then falls back to FTP download if provided. It
raises ``RuntimeError`` when every strategy fails; an empty directory is never
returned as evidence of a genome download.

### `get_genome_metadata`

Retrieve one assembly report without downloading its package. The request uses
`/genome/accession/{accession}/dataset_report` and requires a versioned accession.

```python
from metainformant.dna.external.genomes import get_genome_metadata

metadata = get_genome_metadata("GCF_003254395.2")
print(metadata["organism_name"])  # Apis mellifera
```

Reported fields include organism/taxonomy, assembly name/level/status, source
database and paired accession, assembly statistics, release date and submitter.
`genome_size` is an integer number of bases; unreported fields are `None`.
`AssemblyNotFoundError` means no matching report. Request or response failures
raise `NCBIDatasetsError`; both are available from `dna.external.ncbi_datasets`,
and request errors derive from `requests.RequestException`. Failures are not
returned as metadata dictionaries carrying an `error` field.

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

Search assembly reports by scientific/common name or taxonomy ID using
`/genome/taxon/{taxon}/dataset_report`, following pagination up to positive
integer `max_results`. Optional keyword-only arguments are `reference_only`
(default `False`), `assembly_source` (default `"refseq"`; also `"genbank"` or
`"all"`) and `client`. A GenBank-only species needs `"genbank"` or `"all"`.
Results follow service order. `[]` means no matches (also an unknown organism);
request failures raise `NCBIDatasetsError`.

```python
from metainformant.dna.external.genomes import list_genome_assemblies

assemblies = list_genome_assemblies("Drosophila melanogaster", max_results=5)
# List of dicts with accession, organism_name, assembly_level, genome_size, release_date
```

### `download_reference_genome`

Select a reported reference assembly and download its package, with optional
`include=` and `client=`. `None` means that no reference is listed. Lookup and
download failures raise; the first arbitrary search result is not assumed to be
a reference.

```python
from metainformant.dna.external.genomes import download_reference_genome

# genome_path = download_reference_genome("Apis mellifera", "output/genomes/")
```

### `get_chromosome_lengths`

Retrieve all paginated sequence reports through
`/genome/accession/{accession}/sequence_reports`; Biopython is not required.
The default `names="accession"` keys lengths by sequence accession.
`names="chromosome"` uses reported chromosome names for assembled molecules and
accession for other sequences. A missing assembly or unsuccessful request raises.

```python
from metainformant.dna.external.genomes import get_chromosome_lengths

# lengths = get_chromosome_lengths("GCF_000001405.39")
# Keys are sequence accessions by default; request names="chromosome" for chromosome labels.
```

### `validate_genome_files`

Search recursively for FASTA (`.fna`, `.fa`, `.fasta`) and GFF (`.gff`, `.gff3`),
each optionally gzip-compressed. This includes the extracted
`ncbi_dataset/data/<accession>/` layout. The result contains file counts,
issues and a `valid` flag.

```python
from metainformant.dna.external.genomes import validate_genome_files

results = validate_genome_files("output/genomes/GCF_000001405.39/")
# {"valid": True, "fasta_files": 1, "gff_files": 1, "total_files": 5, "issues": []}
```

---

## Configuration and Rate Limits

E-utilities clients and `get_chromosome_lengths` share one explicit contact policy. An explicit ``email=``
argument takes precedence over ``NCBI_EMAIL``. If neither is available, callers
must pass ``allow_anonymous=True``; anonymous mode warns, leaves the request
identity unset, and is recorded as ``contact_mode="anonymous"`` on the client.
No configuration file or test fixture invents or persists an email.
`NCBIDatasetsClient` accepts optional contact information and does not require
anonymous opt-in. All three clients honor an explicit API key before
`NCBI_API_KEY`. Genome helpers accept a configured `client=`.

Set credentials via constructor parameters or environment variables:

```python
import os
os.environ["NCBI_EMAIL"] = "user@example.com"
os.environ["NCBI_API_KEY"] = "your_key_here"
```

Without an API key, the client automatically sleeps 0.4 seconds between requests
to respect the 3 requests/second limit. All clients accept ``timeout=`` in
seconds (E-utilities default ``30``, Datasets default ``60``); non-positive
values are rejected. Request exceptions redact credential-bearing URLs.

---

## See Also

- **[Accessions](accessions.md)** -- Accession format validation
- **[DNA Module Overview](../dna/)** -- Full DNA module documentation
- `metainformant.dna.external.entrez` -- Entrez source module
- `metainformant.dna.external.ncbi` -- NCBI client source module
- `metainformant.dna.external.genomes` -- Genome retrieval source module
- `metainformant.dna.external.ncbi_datasets` -- Datasets v2 transport and request builders
