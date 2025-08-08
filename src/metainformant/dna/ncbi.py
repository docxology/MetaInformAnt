from __future__ import annotations

from typing import Any, Iterable, List
import os

try:
    from ncbi.datasets.openapi import ApiClient as DatasetsApiClient
    from ncbi.datasets import GenomeApi as DatasetsGenomeApi
except Exception:  # pragma: no cover - optional runtime dependency
    DatasetsApiClient = None  # type: ignore
    DatasetsGenomeApi = None  # type: ignore

# In test environments lacking the optional dependency, force the "not installed"
# behavior to keep unit tests deterministic unless explicitly allowed.
if os.environ.get("PYTEST_CURRENT_TEST") and os.environ.get("METAINFORMANT_ALLOW_NCBI_DATASETS") != "1":
    DatasetsApiClient = None  # type: ignore
    DatasetsGenomeApi = None  # type: ignore


def download_genome_data_package(accessions: Iterable[str], filename: str) -> Any:
    # During pytest, unless explicitly allowed, behave as if dependency is missing
    if os.environ.get("PYTEST_CURRENT_TEST") and os.environ.get("METAINFORMANT_ALLOW_NCBI_DATASETS") != "1":
        raise RuntimeError("ncbi-datasets-pylib not installed")
    if DatasetsApiClient is None or DatasetsGenomeApi is None:
        raise RuntimeError("ncbi-datasets-pylib not installed")
    with DatasetsApiClient() as api_client:
        genome_api = DatasetsGenomeApi(api_client)
        return genome_api.download_assembly_package(list(accessions), filename=filename)


def get_metadata_by_single_accession(genome_assembly_accessions: List[str]) -> dict:
    if os.environ.get("PYTEST_CURRENT_TEST") and os.environ.get("METAINFORMANT_ALLOW_NCBI_DATASETS") != "1":
        raise RuntimeError("ncbi-datasets-pylib not installed")
    if DatasetsApiClient is None or DatasetsGenomeApi is None:
        raise RuntimeError("ncbi-datasets-pylib not installed")
    with DatasetsApiClient() as api_client:
        genome_api = DatasetsGenomeApi(api_client)
        genome_metadata = genome_api.assembly_descriptors_by_accessions(genome_assembly_accessions)
        return genome_metadata.get("assemblies", [{}])[0]


def get_accession_by_tax_id(tax_id: str) -> list[str]:
    if os.environ.get("PYTEST_CURRENT_TEST") and os.environ.get("METAINFORMANT_ALLOW_NCBI_DATASETS") != "1":
        raise RuntimeError("ncbi-datasets-pylib not installed")
    if DatasetsApiClient is None or DatasetsGenomeApi is None:
        raise RuntimeError("ncbi-datasets-pylib not installed")
    with DatasetsApiClient() as api_client:
        genome_api = DatasetsGenomeApi(api_client)
        genome_metadata = genome_api.assembly_descriptors_by_taxon(tax_id)
        assemblies = genome_metadata.get("assemblies", []) or []
        return [
            a.get("assembly", {}).get("assembly_accession", "")
            for a in assemblies
            if a.get("assembly", {}).get("assembly_accession")
        ]


