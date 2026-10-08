"""Real-file fallback promotion preserves layout failures and completed mates."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from metainformant.rna.engine.streaming_orchestrator import (
    _promote_validated_fastq_inputs,
    _validated_local_fastq_inputs,
)


def _write_fastq(path: Path, payload: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        handle.write(payload)
    return path


def test_layout_failed_existing_file_does_not_discard_validated_fallback(tmp_path: Path) -> None:
    accession = "SRR15496110"
    destination_dir = tmp_path / "work/getfastq" / accession
    # Header form captured from the archive's HTTP Range response, not an interleaved pair.
    destination = _write_fastq(
        destination_dir / f"{accession}.fastq.gz",
        "@SRR15496110.1 1/2\nACGT\n+\n!!!!\n@SRR15496110.2 2/2\nTGCA\n+\n!!!!\n",
    )
    failed_bytes = destination.read_bytes()
    stage_dir = tmp_path / "extraction" / accession
    staged = _write_fastq(
        stage_dir / destination.name,
        "@SRR15496110.1 1 length=4\nACGT\n+\n!!!!\n@SRR15496110.2 2 length=4\nTGCA\n+\n!!!!\n",
    )
    staged_bytes = staged.read_bytes()
    validated = _validated_local_fastq_inputs(stage_dir, accession, True)
    assert validated == [staged]

    promoted = _promote_validated_fastq_inputs(validated, destination_dir, accession, True)

    assert promoted
    assert destination.read_bytes() == staged_bytes
    assert destination.with_name(destination.name + ".invalid").read_bytes() == failed_bytes
    assert _validated_local_fastq_inputs(destination_dir, accession, True) == [destination]


def test_completed_numbered_mate_is_preserved(tmp_path: Path) -> None:
    accession = "SRR_PAIRED"
    destination_dir = tmp_path / "work/getfastq" / accession
    mate1 = _write_fastq(destination_dir / f"{accession}_1.fastq.gz", "@read/1\nACGT\n+\n!!!!\n")
    before = mate1.read_bytes(), mate1.stat().st_mtime_ns
    stage_dir = tmp_path / "stage" / accession
    _write_fastq(stage_dir / mate1.name, "@read/1\nACGT\n+\n!!!!\n")
    _write_fastq(stage_dir / f"{accession}_2.fastq.gz", "@read/2\nTGCA\n+\n!!!!\n")
    validated = _validated_local_fastq_inputs(stage_dir, accession, True)

    assert _promote_validated_fastq_inputs(validated, destination_dir, accession, True)

    assert (mate1.read_bytes(), mate1.stat().st_mtime_ns) == before
    assert (destination_dir / f"{accession}_2.fastq.gz").is_file()


@pytest.mark.parametrize(
    "payload",
    [
        "@read/1\nACGT\n+\n!!!\n@next/1\nACGT\n+\n!!!!\n",
        "@read/1\nACGT\n+\n!!!!\n@next/1\nACGT\n+\n",
        "@read/1\nACGT\n+\n!!!!\n@next/1\nACGT\n+\n!!!!\n@read/2\nTGCA\n+\n!!!!\n",
    ],
)
def test_malformed_and_nonadjacent_paired_streams_remain_rejected(tmp_path: Path, payload: str) -> None:
    directory = tmp_path / "work/getfastq/SRR_BAD"
    _write_fastq(directory / "SRR_BAD.fastq.gz", payload)
    assert _validated_local_fastq_inputs(directory, "SRR_BAD", True) == []


@pytest.mark.parametrize("corrupt_gzip", ["layout", "bad-magic", "truncated-stream"])
def test_failed_byte_archive_collisions_preserve_every_version(tmp_path: Path, corrupt_gzip: str) -> None:
    accession = "SRR_REPEAT"
    destination_dir = tmp_path / "work/getfastq" / accession
    destination = _write_fastq(destination_dir / f"{accession}.fastq.gz", "@a/2\nACGT\n+\n!!!!\n@b/2\nTGCA\n+\n!!!!\n")
    if corrupt_gzip == "bad-magic":
        destination.write_bytes(b"not a gzip stream")
    if corrupt_gzip == "truncated-stream":
        _write_fastq(destination, "@a/1\nACGT\n+\n!!!!\n@b/1\nTGCA\n+\n!!!!\n")
        destination.write_bytes(destination.read_bytes()[:-8])
    rejected = destination.read_bytes()
    archived = destination.with_name(destination.name + ".invalid")
    archived.write_bytes(b"older rejected bytes")
    stage_dir = tmp_path / "staging" / accession
    staged = _write_fastq(stage_dir / destination.name, "@a/1\nACGT\n+\n!!!!\n@b/1\nTGCA\n+\n!!!!\n")
    valid_bytes = staged.read_bytes()
    inputs = _validated_local_fastq_inputs(stage_dir, accession, True)

    assert _promote_validated_fastq_inputs(inputs, destination_dir, accession, True)

    assert archived.read_bytes() == b"older rejected bytes"
    assert destination.with_name(destination.name + ".invalid.1").read_bytes() == rejected
    assert destination.read_bytes() == valid_bytes


def test_missing_staged_source_never_reports_success_or_loses_failed_bytes(tmp_path: Path) -> None:
    accession = "SRR_MISSING"
    directory = tmp_path / "work/getfastq" / accession
    destination = _write_fastq(directory / f"{accession}.fastq.gz", "@a/2\nACGT\n+\n!!!!\n@b/2\nTGCA\n+\n!!!!\n")
    rejected = destination.read_bytes()
    missing = tmp_path / "staging" / destination.name

    with pytest.raises(FileNotFoundError):
        _promote_validated_fastq_inputs([missing], directory, accession, True)

    assert destination.with_name(destination.name + ".invalid").read_bytes() == rejected
    assert _validated_local_fastq_inputs(directory, accession, True) == []


def test_valid_existing_singleton_retains_its_own_bytes_and_witness(tmp_path: Path) -> None:
    accession = "SRR_VALID"
    directory = tmp_path / "work/getfastq" / accession
    destination = _write_fastq(directory / f"{accession}.fastq.gz", "@a/1\nACGT\n+\n!!!!\n@b/1\nTGCA\n+\n!!!!\n")
    before = destination.read_bytes(), destination.stat().st_mtime_ns
    stage_dir = tmp_path / "stage" / accession
    staged = _write_fastq(stage_dir / destination.name, "@c/1\nAAAA\n+\n!!!!\n@d/1\nCCCC\n+\n!!!!\n")
    inputs = _validated_local_fastq_inputs(stage_dir, accession, True)

    assert _promote_validated_fastq_inputs(inputs, directory, accession, True)

    assert (destination.read_bytes(), destination.stat().st_mtime_ns) == before
    assert not staged.exists()
    assert not destination.with_name(destination.name + ".invalid").exists()
    assert _validated_local_fastq_inputs(directory, accession, True) == [destination]
