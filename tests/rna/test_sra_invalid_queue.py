from __future__ import annotations

import threading
from pathlib import Path

import pytest

from metainformant.rna.engine.streaming_orchestrator import (
    StreamingPipelineOrchestrator,
    _cached_sra_validation_result,
    _sra_source_fingerprint,
    _write_sra_validation_result,
)


@pytest.mark.parametrize("changed", [False, True])
def test_invalid_sra_witness_skips_busy_slot_only_while_current(tmp_path: Path, changed: bool) -> None:
    source = tmp_path / "SRR_INVALID.sra"
    source.write_bytes(b"NCBI.sra" + bytes(64))
    orchestrator = StreamingPipelineOrchestrator(config_dir=tmp_path, log_dir=tmp_path / "logs",
                                                db_path=tmp_path / "progress.db")

    def record_invalid() -> None:
        _write_sra_validation_result(source, fingerprint=_sra_source_fingerprint(source),
                                     valid=False, returncode=3, detail="corrupt SRA archive")

    record_invalid()
    assert orchestrator._validate_local_sra_archive(source, "SRR_INVALID") is False
    original = source.read_bytes()
    if changed:
        source.write_bytes(b"NCBI.sra" + bytes(128))
        assert _cached_sra_validation_result(source) is None
    expected = source.read_bytes()
    semaphore = threading.BoundedSemaphore(1)
    orchestrator._fasterq_semaphore = semaphore
    semaphore.acquire()
    entered = threading.Event()
    completed = threading.Event()
    outcomes: list[bool] = []

    def extract() -> None:
        entered.set()
        try:
            outcomes.append(orchestrator._extract_sra_to_fastq_bounded(
                source, tmp_path / "sample", tmp_path, "SRR_INVALID", False))
        finally:
            completed.set()

    worker = threading.Thread(target=extract)
    worker.start()
    try:
        assert entered.wait(timeout=2)
        assert completed.wait(timeout=1) is (not changed)
    finally:
        if changed:
            record_invalid()
        semaphore.release()
        worker.join(timeout=10)
        orchestrator.db.close()
    assert not worker.is_alive()
    assert outcomes == [False]
    assert source.read_bytes() == expected
    assert (expected != original) is changed
    assert _cached_sra_validation_result(source) is False
