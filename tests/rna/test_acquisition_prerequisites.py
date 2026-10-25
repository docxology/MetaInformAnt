"""Native backend decisions and fail-closed acquisition admission."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from metainformant.rna.engine.acquisition_prerequisites import (
    PrerequisiteError,
    classify_task_prerequisites,
    require_worker_prerequisites,
)


def _metadata(path: Path, *, platform: str, sampled: str = "yes") -> Path:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "run",
                "scientific_name",
                "platform",
                "lib_layout",
                "spot_length",
                "is_sampled",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "run": "SRR32701718",
                "scientific_name": "Apis mellifera",
                "platform": platform,
                "lib_layout": "single",
                "spot_length": "100",
                "is_sampled": sampled,
            }
        )
    return path


@pytest.mark.parametrize(
    "source,selection,title,expected",
    [
        ("METAGENOMIC", "PCR", "16Sseq of bombus terrerist", "unresolved"),
        ("METAGENOMIC", "PCR", "16S-rRNA -Apis mellifera:workers", "unresolved"),
        ("metagenomic", "pcr", "16S  RNA of Apis mellifera", "unresolved"),
        ("GENOMIC", "PCR", "16Sseq of bee gut", "unresolved"),
        ("genomic", "pcr", "16Sseq of bee gut", "unresolved"),
        ("TRANSCRIPTOMIC", "PCR", "Host gut RNA-seq", "ready"),
        ("GENOMIC", "RANDOM", "16S community RNA", "ready"),
        ("GENOMIC", "PCR", "Host gut RNA-seq", "ready"),
        ("METAGENOMIC", "RANDOM", "16S community RNA", "ready"),
        ("METAGENOMIC", "PCR", "Mixed gut metatranscriptome", "ready"),
    ],
)
def test_explicit_amplicon_contract_conflict_stays_visible_without_exclusion(
    tmp_path: Path, source: str, selection: str, title: str, expected: str
) -> None:
    metadata = _metadata(tmp_path / "metadata.tsv", platform="ILLUMINA")
    with metadata.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    rows[0].update(lib_source=source, lib_selection=selection, exp_title=title, lib_strategy="RNA-Seq")
    with metadata.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    before = metadata.read_bytes()
    result = classify_task_prerequisites(
        metadata_path=metadata,
        tasks=[{"task_id": "apis_mellifera/SRR32701718", "accession": "SRR32701718", "batch_index": 1}],
    )[0]
    assert result.status == expected
    assert result.task_id == "apis_mellifera/SRR32701718"
    assert metadata.read_bytes() == before
    if expected == "unresolved":
        assert "16S" in result.reason and source.upper() in result.reason
        with pytest.raises(PrerequisiteError, match="metadata contract"):
            require_worker_prerequisites([result])


@pytest.mark.parametrize(
    "platform,backend,technology,status",
    [
        ("ILLUMINA", "kallisto", "", "ready"),
        ("OXFORD_NANOPORE", "oarfish", "ont-cdna", "unresolved"),
        ("PACBIO_SMRT", "oarfish", "pac-bio", "unresolved"),
    ],
)
def test_native_backend_requirements(tmp_path: Path, platform: str, backend: str, technology: str, status: str) -> None:
    metadata = _metadata(tmp_path / "metadata.tsv", platform=platform)
    results = classify_task_prerequisites(
        metadata_path=metadata,
        tasks=[
            {
                "task_id": "apis_mellifera/SRR32701718",
                "accession": "SRR32701718",
                "batch_index": 1,
            }
        ],
    )
    assert (
        results[0].backend,
        results[0].sequencing_technology,
        results[0].status,
    ) == (backend, technology, status)


def test_worker_rejects_long_reads_before_tool_or_acquisition(tmp_path: Path) -> None:
    metadata = _metadata(tmp_path / "metadata.tsv", platform="OXFORD_NANOPORE")
    results = classify_task_prerequisites(
        metadata_path=metadata,
        tasks=[
            {
                "task_id": "apis_mellifera/SRR32701718",
                "accession": "SRR32701718",
                "batch_index": 1,
            }
        ],
    )
    with pytest.raises(PrerequisiteError, match="amended frozen transcript FASTA/MMI"):
        require_worker_prerequisites(results)


def test_native_batch_selection_mismatch_is_unresolved(tmp_path: Path) -> None:
    metadata = _metadata(tmp_path / "metadata.tsv", platform="ILLUMINA")
    results = classify_task_prerequisites(
        metadata_path=metadata,
        tasks=[{"task_id": "apis_mellifera/SRR999", "accession": "SRR999", "batch_index": 1}],
    )
    assert results[0].status == "unresolved"
    assert "batch index" in results[0].reason


def test_duplicate_metadata_fails_closed(tmp_path: Path) -> None:
    metadata = _metadata(tmp_path / "metadata.tsv", platform="ILLUMINA")
    with metadata.open("a") as handle:
        handle.write(metadata.read_text().splitlines()[1] + "\n")
    with pytest.raises(ValueError, match="duplicate run IDs"):
        classify_task_prerequisites(metadata_path=metadata, tasks=[])


def test_missing_native_kallisto_is_rejected(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    metadata = _metadata(tmp_path / "metadata.tsv", platform="ILLUMINA")
    results = classify_task_prerequisites(
        metadata_path=metadata,
        tasks=[
            {
                "task_id": "apis_mellifera/SRR32701718",
                "accession": "SRR32701718",
                "batch_index": 1,
            }
        ],
    )
    monkeypatch.setenv("PATH", str(tmp_path / "no-executables"))
    with pytest.raises(FileNotFoundError, match="kallisto executable not found"):
        require_worker_prerequisites(results)


@pytest.mark.parametrize("batch", [None, True, False, "", "not-a-batch", "1.5", 0, -1, 2])
def test_invalid_batch_input_is_unresolved(tmp_path: Path, batch: str | int | bool | None) -> None:
    metadata = _metadata(tmp_path / "metadata.tsv", platform="ILLUMINA")
    results = classify_task_prerequisites(
        metadata_path=metadata,
        tasks=[{"task_id": "apis_mellifera/SRR32701718", "accession": "SRR32701718", "batch_index": batch}],
    )
    assert results[0].status == "unresolved"
    assert "batch index" in results[0].reason


def test_string_batch_input_keeps_native_selection(tmp_path: Path) -> None:
    metadata = _metadata(tmp_path / "metadata.tsv", platform="ILLUMINA")
    results = classify_task_prerequisites(
        metadata_path=metadata,
        tasks=[{"task_id": "apis_mellifera/SRR32701718", "accession": "SRR32701718", "batch_index": "1"}],
    )
    assert results[0].status == "ready"
    assert results[0].backend == "kallisto"
