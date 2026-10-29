"""Live tests for the NCBI Datasets v2 genome functions.

These talk to ``https://api.ncbi.nlm.nih.gov/datasets/v2`` and are marked ``network``.  A test is skipped when
the service cannot be reached or is failing on its side (no connection, HTTP 5xx, HTTP 429).  An HTTP 4xx
answer fails the test: it means the query itself is wrong, which is exactly what the port had to fix.

The expectations were read from the service on 2026-10-28 (see ``tests/data/ncbi_datasets/README.md``).  Assembly
accessions are versioned and stable; which assembly is the reference for a species can change, so those
tests check structure instead of an accession.
"""

from __future__ import annotations

import json
import os
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator

import pytest
import requests

from metainformant.dna.external import genomes, ncbi_datasets
from metainformant.dna.external.ncbi_datasets import (
    AssemblyNotFoundError,
    DatasetsRequest,
    NCBIDatasetsClient,
    NCBIDatasetsError,
)

pytestmark = pytest.mark.network

AMEL = "GCF_003254395.2"


@contextmanager
def skip_if_service_unavailable() -> Iterator[None]:
    """Skip on transport failures and server-side trouble; let a client error (HTTP 4xx) fail the test."""
    try:
        yield
    except requests.RequestException as exc:
        status = getattr(exc, "status_code", None)
        if status is not None and status < 500 and status != 429:
            raise
        pytest.skip(f"NCBI Datasets unavailable: {exc}")


def test_live_reference_listing_finds_the_honey_bee_reference() -> None:
    with skip_if_service_unavailable():
        assemblies = genomes.list_genome_assemblies("Apis mellifera", max_results=5, reference_only=True)

    assert assemblies, "NCBI lists a reference genome for Apis mellifera"
    assert assemblies[0]["refseq_category"] == "reference genome"
    assert assemblies[0]["accession"].startswith("GCF_") and assemblies[0]["organism_name"] == "Apis mellifera"
    assert isinstance(assemblies[0]["genome_size"], int) and assemblies[0]["genome_size"] > 200_000_000


def test_live_metadata_of_a_versioned_assembly() -> None:
    with skip_if_service_unavailable():
        metadata = genomes.get_genome_metadata(AMEL)

    assert metadata["accession"] == AMEL
    assert metadata["organism_name"] == "Apis mellifera" and metadata["tax_id"] == 7460
    assert metadata["genome_size"] == 225234541
    assert metadata["contig_count"] == 227 and metadata["scaffold_count"] == 176
    assert metadata["assembly_level"] == "Chromosome"


def test_live_unknown_accession_is_a_lookup_error() -> None:
    with skip_if_service_unavailable():
        with pytest.raises(AssemblyNotFoundError):
            genomes.get_genome_metadata("GCF_999999999.9")


def test_live_unknown_organism_lists_nothing() -> None:
    with skip_if_service_unavailable():
        assert genomes.list_genome_assemblies("Notaspecies nonexistens") == []


def test_live_paging_collects_every_report_once() -> None:
    query = ncbi_datasets.build_taxon_report_request("Apis mellifera", returned_content="ASSM_ACC", page_size=2)

    with skip_if_service_unavailable(), NCBIDatasetsClient() as client:
        collected = client.collect_reports(query)

    accessions = [report["accession"] for report in collected.reports]
    assert collected.total_count is not None and collected.total_count > 2
    assert len(accessions) == collected.total_count
    assert len(set(accessions)) == len(accessions), "a page was read twice"


def test_live_query_for_one_result_terminates() -> None:
    """NCBI sends a page token even with nothing left; the walk must end on the empty page that follows."""
    query = ncbi_datasets.build_taxon_report_request("Apis mellifera", reference_only=True, page_size=1)

    with skip_if_service_unavailable(), NCBIDatasetsClient() as client:
        collected = client.collect_reports(query)

    assert len(collected.reports) == collected.total_count >= 1


def test_live_old_listing_endpoint_is_gone_and_reported_loudly() -> None:
    """The endpoint ``list_genome_assemblies`` used to call answers 404, which the client now raises."""
    with NCBIDatasetsClient() as client:
        try:
            client.get_json(DatasetsRequest("/genome", {"taxon": "Apis mellifera", "page_size": "1"}))
        except NCBIDatasetsError as exc:
            if exc.status_code is None or exc.status_code >= 500:
                pytest.skip(f"NCBI Datasets unavailable: {exc}")
            assert exc.status_code == 404
        else:
            pytest.fail("the retired /genome listing endpoint answered")


def test_live_sequence_lengths_of_the_honey_bee_genome() -> None:
    email = os.environ.get("NCBI_EMAIL", "").strip()
    if not email:
        pytest.skip("Set NCBI_EMAIL to run this live test; tests never invent contact data")

    with skip_if_service_unavailable():
        by_accession = genomes.get_chromosome_lengths(AMEL, email=email)
        by_chromosome = genomes.get_chromosome_lengths(AMEL, email=email, names="chromosome")

    assert len(by_accession) == len(by_chromosome) == 177  # 16 linkage groups, MT and unplaced scaffolds
    assert by_accession["NC_037638.1"] == by_chromosome["LG1"] == 27754200


def test_live_reference_genome_download_of_a_small_genome(tmp_path: Path) -> None:
    """Downloads the 580 kb Mycoplasma genitalium reference, a few hundred kilobytes of zip."""
    with skip_if_service_unavailable():
        genome_dir = genomes.download_reference_genome("Mycoplasma genitalium", tmp_path)

    assert genome_dir is not None, "NCBI lists a reference genome for Mycoplasma genitalium"
    assert genome_dir.parent == tmp_path and genome_dir.name.startswith("GCF_")
    assert list(tmp_path.glob("*.zip")) == [], "the archive is removed after extraction"
    validation = genomes.validate_genome_files(genome_dir)
    assert validation["valid"] is True and validation["fasta_files"] >= 1, validation
    report_path = next(genome_dir.rglob("assembly_data_report.jsonl"))
    assert json.loads(report_path.read_text().splitlines()[0])["accession"] == genome_dir.name
