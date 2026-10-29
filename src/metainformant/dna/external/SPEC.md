# Specification: external

## 🎯 Scope
External database interfaces (NCBI, Entrez, genome downloads).

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 6 Python modules
- **Key Concepts**: `DatasetsRequest` (path + query of one Datasets v2 call), `ReportPage` (one decoded page of
  a report endpoint), `NCBIRequestError` / `NCBIDatasetsError` (failed or unusable response),
  `AssemblyNotFoundError` (NCBI answered but has no such assembly).

## 🔌 API Definition
### Exports
- `__init__.py`
- `_http.py`: shared GET helper; failures raise, credentials are redacted from error text
- `entrez.py`
- `genomes.py`
- `ncbi.py`
- `ncbi_datasets.py`: pure request builders and response parsers for the Datasets v2 REST API, plus the client

### Contracts
- A failed request raises `NCBIRequestError` (a `requests.RequestException`); it is never answered with an
  empty list, dict, string or `None`.
- NCBI's own empty answers (unknown accession `{}`, taxon without matching assemblies) are results:
  `list_genome_assemblies()` returns `[]`, `download_reference_genome()` returns `None`,
  `get_genome_metadata()` raises `AssemblyNotFoundError`.
- Fields a report or record does not carry are `None`.
- Datasets accessions are versioned assembly accessions (`GC[FA]_#########.#`).
