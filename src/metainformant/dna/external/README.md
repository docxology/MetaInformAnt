# External

Clients for NCBI Entrez and Datasets APIs to search databases, download genome assemblies, and retrieve sequence records.

## Contents

| File | Purpose |
|------|---------|
| `entrez.py` | Entrez E-utilities client for GenBank search and sequence retrieval |
| `genomes.py` | Genome assembly download and validation via NCBI Datasets API |
| `ncbi.py` | General NCBI client for nucleotide, protein, and taxonomy queries |

## Key Classes and Functions

| Symbol | Description |
|--------|-------------|
| `EntrezClient` | Session-based client with API key and rate limit support |
| `EntrezClient.search()` | Search any Entrez database (nucleotide, protein, etc.) |
| `fetch_genbank_record()` | Retrieve a GenBank record by accession |
| `fetch_fasta_sequence()` | Download a single FASTA sequence by accession |
| `NCBIClient` | General-purpose NCBI API client |
| `NCBIClient.search_nucleotide()` | Search the NCBI Nucleotide database |
| `download_genome_package()` | Download genome assembly files from NCBI Datasets |
| `validate_accession()` | Check whether a genome accession string is valid |
| `get_genome_metadata()` | Retrieve metadata for a genome assembly |
| `list_genome_assemblies()` | List available assemblies for an organism |
| `get_chromosome_lengths()` | Fetch chromosome lengths for an assembly |
| `get_taxonomy_info()` | Look up taxonomy information by NCBI tax ID |

## Usage

```python
from metainformant.dna.external.entrez import EntrezClient, fetch_fasta_sequence
from metainformant.dna.external.genomes import download_genome_package, list_genome_assemblies
from metainformant.dna.external.ncbi import NCBIClient

client = EntrezClient(email="user@example.com")
results = client.search("nucleotide", "Apis mellifera[ORGN]")
assemblies = list_genome_assemblies("Apis mellifera", max_results=5)
```
