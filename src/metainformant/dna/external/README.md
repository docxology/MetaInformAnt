# External

Clients for NCBI Entrez and Datasets APIs to search databases, download genome assemblies, and retrieve sequence records.

**Contact policy.** E-utilities requests (`entrez.py`, `ncbi.py`, and `get_chromosome_lengths()`) require an
explicit ``email=`` argument, ``NCBI_EMAIL``, or the deliberate ``allow_anonymous=True`` opt-in. The Datasets
genome functions need no contact: they send ``From:`` only when ``email=`` or ``NCBI_EMAIL`` is set, and an
API key (``api_key=`` or ``NCBI_API_KEY``) goes in the ``api-key`` header, never in the URL. Requests are
bounded by a per-call timeout and never use a fabricated contact address.

**Failures raise.** A request that fails (no connection, HTTP error status, a body that is not what was asked
for) raises `NCBIRequestError`, a ``requests.RequestException``; the Datasets functions raise its subclass
`NCBIDatasetsError`. Empty results are returned only when NCBI reported none, and a field NCBI leaves out of
a record is ``None``. Error text never contains ``api_key`` or ``email`` values.

## Contents

| File | Purpose |
|------|---------|
| `entrez.py` | Entrez E-utilities client for GenBank search and sequence retrieval |
| `genomes.py` | Genome assembly metadata, listing and download via the NCBI Datasets v2 API |
| `ncbi.py` | General NCBI client for nucleotide, protein, and taxonomy queries |
| `ncbi_datasets.py` | Datasets v2 request builders, response parsers (pure functions) and `NCBIDatasetsClient` |
| `_http.py` | Shared GET helper, `NCBIRequestError` and credential redaction |

## Key Classes and Functions

| Symbol | Description |
|--------|-------------|
| `EntrezClient` | Session-based client with API key and rate limit support |
| `EntrezClient.search()` | Search any Entrez database (nucleotide, protein, etc.) |
| `fetch_genbank_record()` | Retrieve a GenBank record by accession |
| `fetch_fasta_sequence()` | Download a single FASTA sequence by accession |
| `NCBIClient` | General-purpose NCBI API client |
| `NCBIClient.search_nucleotide()` | Search the NCBI Nucleotide database |
| `download_genome_package()` | Download genome assembly files from NCBI Datasets (streamed, `include=` file types) |
| `download_reference_genome()` | Download the reference genome of a species; `None` if NCBI lists none |
| `validate_accession()` | Check whether a genome accession string is valid |
| `get_genome_metadata()` | Retrieve the assembly report of one versioned accession |
| `list_genome_assemblies()` | List assemblies of an organism or taxon (paged) |
| `get_chromosome_lengths()` | Sequence lengths of an assembly from the Datasets sequence report |
| `validate_genome_files()` | Check an extracted package for FASTA and GFF files |
| `get_taxonomy_info()` | Look up taxonomy information by NCBI tax ID |
| `NCBIDatasetsClient` | Datasets v2 client: API key, contact, timeout, paging, streamed download |

## Datasets v2 endpoints

Base URL `https://api.ncbi.nlm.nih.gov/datasets/v2`; specification
<https://www.ncbi.nlm.nih.gov/datasets/docs/v2/openapi3/openapi3.docs.yaml>.

| Function | Endpoint |
|----------|----------|
| `get_genome_metadata()` | `GET /genome/accession/{accession}/dataset_report` |
| `list_genome_assemblies()`, `download_reference_genome()` | `GET /genome/taxon/{taxon}/dataset_report` (`filters.reference_only`, `filters.assembly_source`, `page_size`, `page_token`) |
| `get_chromosome_lengths()` | `GET /genome/accession/{accession}/sequence_reports` |
| `download_genome_package()` | `GET /genome/accession/{accession}/download` (`include_annotation_type`, repeated) |

Accessions must carry a version (`GCF_003254395.2`); NCBI rejects `GCF_003254395`. Fixtures of real responses
and the status codes seen are in `tests/data/ncbi_datasets/`.

## Usage

```python
from metainformant.dna.external.entrez import EntrezClient, fetch_fasta_sequence
from metainformant.dna.external.genomes import (
    download_reference_genome,
    get_genome_metadata,
    list_genome_assemblies,
)
from metainformant.dna.external.ncbi import NCBIClient
from metainformant.dna.external.ncbi_datasets import NCBIDatasetsClient

client = EntrezClient(email="user@example.com", timeout=30)
results = client.search("nucleotide", "Apis mellifera[ORGN]")

metadata = get_genome_metadata("GCF_003254395.2")  # genome_size is an int, absent fields are None
assemblies = list_genome_assemblies("Apis mellifera", max_results=5, reference_only=True)
with NCBIDatasetsClient(api_key="...", timeout=120) as datasets:  # reuse one session and key
    genome_dir = download_reference_genome("Apis mellifera", "output/genomes", include=["genome", "gff3"], client=datasets)
```
