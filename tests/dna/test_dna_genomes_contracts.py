"""Offline contract tests for ``dna.external.genomes`` that use only its long-standing public signatures.

Each case fails on the version of the module that preceded the NCBI Datasets v2 port: a bad argument was
passed on to the service (or, for ``chr1`` and ``NC_000001``, accepted by ``validate_accession``) and the
failure came back as ``[]``, ``None``, ``{}`` or an error dictionary.  Invalid arguments now fail before any
query is built, so these tests make no network call.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
import requests

from metainformant.dna.external import genomes

AMEL = "GCF_003254395.2"
CONTACT = "tests@example.org"


@pytest.fixture
def dead_proxy(monkeypatch: pytest.MonkeyPatch) -> None:
    """Send every outbound connection of this process to a closed local port.

    The failure is real and immediate (``requests`` and ``urllib`` both honour the proxy variables, lower
    case ones first), nothing leaves the machine, and no service is imitated.
    """
    for name in (
        "http_proxy",
        "https_proxy",
        "all_proxy",
        "ftp_proxy",
        "HTTP_PROXY",
        "HTTPS_PROXY",
        "ALL_PROXY",
        "FTP_PROXY",
    ):
        monkeypatch.setenv(name, "http://127.0.0.1:9")
    for name in ("no_proxy", "NO_PROXY"):
        monkeypatch.setenv(name, "")
    monkeypatch.delenv("NCBI_API_KEY", raising=False)


@pytest.mark.parametrize("accession", ["chr1", "NC_000001", "ASM2732", "GCF_003254395", "", "../etc/passwd"])
def test_metadata_rejects_non_assembly_accessions(accession: str) -> None:
    """``validate_accession`` accepts ``chr1`` and ``NC_000001``, which the assembly report cannot answer."""
    with pytest.raises(ValueError, match="Invalid genome accession"):
        genomes.get_genome_metadata(accession)


@pytest.mark.parametrize("max_results", [0, -3, True, 2.5, "5", None])
def test_listing_rejects_bad_result_limits(max_results: Any) -> None:
    with pytest.raises(ValueError, match="max_results"):
        genomes.list_genome_assemblies("Apis mellifera", max_results=max_results)


@pytest.mark.parametrize("organism", ["", "   ", "Apis mellifera, Pogonomyrmex barbatus"])
def test_listing_rejects_blank_or_multiple_taxa(organism: str) -> None:
    with pytest.raises(ValueError):
        genomes.list_genome_assemblies(organism)


@pytest.mark.parametrize("species", ["", "   "])
def test_reference_lookup_rejects_blank_species(species: str, tmp_path: Path) -> None:
    with pytest.raises(ValueError):
        genomes.download_reference_genome(species, tmp_path)
    assert list(tmp_path.iterdir()) == []


def test_package_rejects_non_assembly_accession_before_touching_the_disk(tmp_path: Path) -> None:
    target = tmp_path / "out"

    with pytest.raises(ValueError, match="Invalid genome accession"):
        genomes.download_genome_package("NC_000001", target)

    assert not target.exists()


def test_package_rejects_unknown_file_types_before_touching_the_disk(tmp_path: Path) -> None:
    """The ``include`` parameter used to be sent as ``include=...``, a name the v2 service ignores."""
    target = tmp_path / "out"

    with pytest.raises(ValueError, match="Unsupported genome package file type"):
        genomes.download_genome_package(AMEL, target, include=["bogus"])

    assert not target.exists()


def test_chromosome_lengths_reject_non_assembly_accession() -> None:
    with pytest.raises(ValueError, match="Invalid genome accession"):
        genomes.get_chromosome_lengths("chr1", email="tests@example.org")


def test_validate_genome_files_sees_the_layout_of_a_genome_package(tmp_path: Path) -> None:
    """A package keeps its files under ``ncbi_dataset/data/<accession>/`` and names the FASTA ``*.fna``."""
    data = tmp_path / "ncbi_dataset" / "data"
    (data / AMEL).mkdir(parents=True)
    (data / "assembly_data_report.jsonl").write_bytes(b"{}\n")
    (data / AMEL / f"{AMEL}_Amel_HAv3.1_genomic.fna").write_bytes(b">x\n")
    (data / AMEL / "genomic.gff").write_bytes(b"##gff-version 3\n")

    results = genomes.validate_genome_files(tmp_path)

    assert results == {"valid": True, "fasta_files": 1, "gff_files": 1, "total_files": 3, "issues": []}


def test_validate_genome_files_counts_files_not_top_level_entries(tmp_path: Path) -> None:
    (tmp_path / "a").mkdir()
    (tmp_path / "a" / "b").mkdir()
    (tmp_path / "a" / "b" / "genome.fa").write_bytes(b">x\n")

    results = genomes.validate_genome_files(tmp_path)

    assert results["total_files"] == 1 and results["fasta_files"] == 1


@pytest.mark.parametrize(
    "ftp_url", ["ftp://example.invalid", "ftp://example.invalid/", "ftp://example.invalid/genomes/all/"]
)
def test_ftp_fallback_refuses_a_directory_url_instead_of_inventing_a_file_name(ftp_url: str, tmp_path: Path) -> None:
    """The old code saved whatever came back as ``genome_data.fasta`` (or as the last directory's name)."""
    target = tmp_path / "out"

    with pytest.raises(ValueError, match="must name a file"):
        genomes._download_from_ftp(ftp_url, target)

    assert not target.exists()


def test_listing_failure_raises_rather_than_returning_empty(dead_proxy: None) -> None:
    """The old function turned every failure, including its own 404, into ``[]``."""
    with pytest.raises(requests.RequestException):
        genomes.list_genome_assemblies("Apis mellifera", max_results=2)


def test_metadata_failure_raises_rather_than_returning_an_error_dict(dead_proxy: None) -> None:
    with pytest.raises(requests.RequestException):
        genomes.get_genome_metadata(AMEL)


def test_reference_lookup_failure_raises_rather_than_returning_none(dead_proxy: None, tmp_path: Path) -> None:
    with pytest.raises(requests.RequestException):
        genomes.download_reference_genome("Apis mellifera", tmp_path)
    assert list(tmp_path.iterdir()) == []


def test_sequence_lengths_failure_raises_rather_than_returning_empty(dead_proxy: None) -> None:
    with pytest.raises(requests.RequestException):
        genomes.get_chromosome_lengths(AMEL, email=CONTACT)
