"""Real Bash contracts for opt-in bulk-HD execution without payload staging."""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import Literal, assert_never

import pytest


def inspect_profile(tmp_path: Path, mode: str | None) -> subprocess.CompletedProcess[str]:
    """Run the actual launcher in an isolated operations tree."""
    source = Path(__file__).resolve().parents[2] / "output/hymenoptera_local_20261019"
    operations = tmp_path / "operations with spaces"
    operations.mkdir()
    shutil.copy2(source / "run.sh", operations / "run.sh")
    data = tmp_path / "existing bulk with spaces"
    data.mkdir()
    environment = {"PATH": os.environ["PATH"], "HOME": os.environ["HOME"],
                   "AMALGKIT_DATA_ROOT": str(data), "AMALGKIT_FROZEN_ROOT": str(tmp_path / "frozen"),
                   "AMALGKIT_LOCAL_FASTQ_STAGE": "yes",
                   "AMALGKIT_LOCAL_QUANT_SCRATCH_DIR": str(tmp_path / "stale SSD setting"),
                   "AMALGKIT_LOCAL_FASTERQ_SCRATCH_DIR": str(tmp_path / "stale extraction setting"),
                   "AMALGKIT_LOCAL_ENA_SCRATCH_DIR": str(tmp_path / "stale ENA setting"),
                   "UV_PYTHON_DOWNLOADS": "never", "UV_NO_SYNC": "1"}
    if mode is not None:
        environment["AMALGKIT_LOCAL_STORAGE_MODE"] = mode
    return subprocess.run(["/bin/bash", str(operations / "run.sh"), "print-storage-profile"],
                          env=environment, capture_output=True, text=True, timeout=20, check=False)


@pytest.mark.parametrize("mode", [None, "ssd_compute", "hd_first"])
def test_profile_query_reports_effective_storage_without_writes(
    tmp_path: Path, mode: Literal["ssd_compute", "hd_first"] | None,
) -> None:
    # Given an isolated launcher and conflicting inherited scratch settings.
    # When its read-only inspection command resolves the selected mode.
    result = inspect_profile(tmp_path, mode)
    # Then it reports effective settings without creating cache or scratch trees.
    assert result.returncode == 0, result.stdout + result.stderr
    values = dict(line.split("=", 1) for line in result.stdout.splitlines())
    assert values["storage_mode"] == (mode or "ssd_compute")
    operations = tmp_path / "operations with spaces"
    assert not (operations / "tmp").exists()
    assert not (operations / "cache").exists()
    assert not (operations / "active_runtime.json").exists()
    assert values["index_cache"] == str(operations / "index-cache")
    assert values["scratch_reserve_gib"] == "128"
    assert values["external_floor_gib"] == "256"
    assert values["system_floor_gib"] == "32"
    match mode:
        case "hd_first":
            assert values["quant_scratch"] == values["fasterq_scratch"] == values["ena_scratch"] == ""
            assert values["fastq_stage"] == "no"
        case None | "ssd_compute":
            assert values["quant_scratch"] == str(operations / "quant-scratch")
            assert values["fasterq_scratch"] == str(operations / "fasterq-scratch")
            assert values["fastq_stage"] == "yes"
            assert values["ena_scratch"] == str(tmp_path / "stale ENA setting")
        case unreachable:
            assert_never(unreachable)


@pytest.mark.parametrize("mode", ["", "HD_FIRST", "../hd_first", "hd_first;touch bad"])
def test_invalid_storage_mode_refuses_before_side_effects(tmp_path: Path, mode: str) -> None:
    # Given a malformed explicit mode, including shell metacharacters.
    # When the real launcher is invoked.
    result = inspect_profile(tmp_path, mode)
    # Then it refuses without interpreting the value or creating controls.
    assert result.returncode == 2
    assert "AMALGKIT_LOCAL_STORAGE_MODE" in result.stderr
    assert not (tmp_path / "operations with spaces/tmp").exists()
    assert not (tmp_path / "operations with spaces/active_runtime.json").exists()
    assert not (tmp_path / "bad").exists()


def test_hd_first_profile_selects_existing_bulk_workspace(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    from metainformant.rna.engine.streaming_orchestrator import (
        StreamingPipelineOrchestrator,
    )

    # Given settings emitted by the real launcher and retained bulk input.
    result = inspect_profile(tmp_path, "hd_first")
    assert result.returncode == 0, result.stderr
    values = dict(line.split("=", 1) for line in result.stdout.splitlines())
    for key, value in (
        ("AMALGKIT_LOCAL_QUANT_SCRATCH_DIR", values["quant_scratch"]),
        ("AMALGKIT_LOCAL_FASTERQ_SCRATCH_DIR", values["fasterq_scratch"]),
        ("AMALGKIT_LOCAL_FASTQ_STAGE", values["fastq_stage"]),
    ):
        monkeypatch.setenv(key, value)
    bulk = tmp_path / "existing bulk with spaces"
    sample = bulk / "getfastq/SRR123"
    sample.mkdir(parents=True)
    source = sample / "retained.fastq.gz"
    source.write_bytes(b"retained input must not be moved")
    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "configs", tmp_path / "logs", db_path=tmp_path / "private-progress.db",
    )

    # When the engine selects a quantification workspace.
    workspace = orchestrator._prepare_local_quant_workspace(
        bulk, "test_species", "SRR123", str(bulk / "metadata.tsv"), ["amalgkit", "quant"],
    )

    # Then the canonical bulk path is used without staging or a new SSD copy.
    assert workspace is None
    assert orchestrator.local_quant_scratch_dir is None
    assert orchestrator.local_fasterq_scratch_dir is None
    assert orchestrator.local_fastq_stage is False
    assert source.read_bytes() == b"retained input must not be moved"
    assert not (tmp_path / "operations with spaces/quant-scratch").exists()
