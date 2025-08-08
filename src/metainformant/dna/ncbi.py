from __future__ import annotations

from typing import Any, Iterable, List

try:
    from ncbi.datasets.openapi import ApiClient as DatasetsApiClient
    from ncbi.datasets import GenomeApi as DatasetsGenomeApi
except Exception:  # pragma: no cover - optional runtime dependency
    DatasetsApiClient = None  # type: ignore
    DatasetsGenomeApi = None  # type: ignore


def download_genome_data_package(accessions: Iterable[str], filename: str) -> Any:
    if DatasetsApiClient is None or DatasetsGenomeApi is None:
        raise RuntimeError("ncbi-datasets-pylib not installed")
    with DatasetsApiClient() as api_client:
        genome_api = DatasetsGenomeApi(api_client)
        return genome_api.download_assembly_package(list(accessions), filename=filename)


def get_metadata_by_single_accession(genome_assembly_accessions: List[str]) -> dict:
    if DatasetsApiClient is None or DatasetsGenomeApi is None:
        raise RuntimeError("ncbi-datasets-pylib not installed")
    with DatasetsApiClient() as api_client:
        genome_api = DatasetsGenomeApi(api_client)
        genome_metadata = genome_api.assembly_descriptors_by_accessions(genome_assembly_accessions)
        return genome_metadata.get("assemblies", [{}])[0]


def get_accession_by_tax_id(tax_id: str) -> list[str]:
    if DatasetsApiClient is None or DatasetsGenomeApi is None:
        raise RuntimeError("ncbi-datasets-pylib not installed")
    with DatasetsApiClient() as api_client:
        genome_api = DatasetsGenomeApi(api_client)
        genome_metadata = genome_api.assembly_descriptors_by_taxon(tax_id)
        total_count = genome_metadata.get("total_count", 0)
        return [
            genome_metadata["assemblies"][i]["assembly"]["assembly_accession"]
            for i in range(0, total_count)
        ]


