"""Actual Kallisto/Amalgkit execution and portable idempotent manifest replay."""

from __future__ import annotations

import gzip
import json
import os
import shutil
import subprocess
import sys
from dataclasses import asdict
from pathlib import Path

import pytest

from metainformant.rna.engine.acquisition_manifest import sha256_file
from metainformant.rna.engine.acquisition_scheduling import PlanningAssumptions
from metainformant.rna.engine.acquisition_worker import run_manifest
from metainformant.rna.engine.streaming_orchestrator import (
    _promote_validated_fastq_inputs,
    _validated_local_fastq_inputs,
)


@pytest.mark.external_tool
@pytest.mark.parametrize("dedicated_recovery", [False, True])
@pytest.mark.parametrize("promote_fallback", [False, True])
@pytest.mark.parametrize("planning_profile", [False, True])
@pytest.mark.parametrize(
    "reference_target", ["Test species", "Test species subspecies"]
)
def test_real_manifest_quantification_and_second_run_reuses_bytes(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    reference_target: str,
    planning_profile: bool,
    promote_fallback: bool,
    dedicated_recovery: bool,
    native_acquisition_factory,
) -> None:
    monkeypatch.setenv(
        "AMALGKIT_PIPELINE_DEDICATED_RECOVERY", "yes" if dedicated_recovery else "no"
    )
    fixture = native_acquisition_factory(
        reference_target=reference_target, paired_metadata=promote_fallback
    )
    data, work, index = fixture.data, fixture.work, fixture.index
    metadata, config_dir, manifest = (
        fixture.metadata,
        fixture.config_dir,
        fixture.manifest,
    )
    reads = work / "getfastq/SRR123"
    if promote_fallback:
        fallback_dir = tmp_path / "fallback/SRR123"
        fallback_dir.mkdir(parents=True)
        source = fallback_dir / "SRR123.fastq.gz"
        destination = reads / source.name
        shutil.copyfile(destination, source)
        validated_bytes = source.read_bytes()
        with gzip.open(destination, "wt") as handle:
            for i in range(2):
                handle.write(f"@SRR123.{i + 1} {i + 1}/2\nACGT\n+\n!!!!\n")
        rejected_bytes = destination.read_bytes()
        inputs = _validated_local_fastq_inputs(fallback_dir, "SRR123", True)
        assert _promote_validated_fastq_inputs(inputs, reads, "SRR123", True)
        assert destination.read_bytes() == validated_bytes
        assert (
            destination.with_name(destination.name + ".invalid").read_bytes()
            == rejected_bytes
        )
    kwargs = {
        "manifest_path": manifest,
        "data_root": data,
        "config_dir": config_dir,
        "workers": 2,
        "threads": 4,
        "fastq_threads": 1,
        "compression_threads": 1,
        "validation_slots": 1,
        "quant_slots": 2,
        "fasterq_slots": 1,
        "max_in_flight": 1,
    }
    # Refuse both a known task past cutoff and a task with no duration bound.
    original_manifest = manifest.read_bytes()
    original_snapshot = (tmp_path / "snapshot.json").read_bytes()
    for scenario, deadline, planning, reason in (
        (
            "expired-deadline",
            "1",
            {"planning_seconds": 1, "planning_drain_seconds": 1},
            "insufficient",
        ),
        ("missing-estimate", "9999999999", {}, "positive"),
        ("missing-deadline", "", {}, "no hard deadline"),
    ):
        task = json.loads(original_manifest)
        task.update(planning)
        manifest.write_text(json.dumps(task) + "\n")
        envelope = json.loads(original_snapshot)
        envelope["manifest_sha256"] = sha256_file(manifest)
        (tmp_path / "snapshot.json").write_text(json.dumps(envelope))
        environment = dict(
            os.environ,
            AMALGKIT_WORKER_DEADLINE_EPOCH=deadline,
            AMALGKIT_CLOUD_MAX_RAW_BYTES="15000",
        )
        completed = subprocess.run(
            [
                sys.executable,
                "scripts/rna/acquisition_worker.py",
                "--manifest",
                str(manifest),
                "--data-root",
                str(data),
                "--config-dir",
                str(config_dir),
                "--workers",
                "1",
                "--threads",
                "1",
            ],
            env=environment,
            text=True,
            capture_output=True,
            check=False,
        )
        (tmp_path / f"{scenario}-cli.txt").write_text(
            completed.stdout + completed.stderr
        )
        assert completed.returncode == 1
        deferred = json.loads((data / "cloud_worker_result.json").read_text())
        assert deferred["counts"] == {"unresolved": 1}
        assert deferred["results"][0]["status"] == "unresolved_deadline"
        assert not (work / "quant/SRR123").exists()
        assert reason in deferred["results"][0]["error"]
    manifest.write_bytes(original_manifest)
    (tmp_path / "snapshot.json").write_bytes(original_snapshot)
    if planning_profile:
        task = json.loads(original_manifest)
        task["total_bases"] = 6750
        manifest.write_text(json.dumps(task) + "\n")
        envelope = json.loads(original_snapshot)
        envelope["manifest_sha256"] = sha256_file(manifest)
        envelope["planning_assumptions"] = asdict(
            PlanningAssumptions(
                transfer_bytes_per_second=1000000,
                extraction_bases_per_second=0,
                quant_bases_per_second=1000000,
                source="deterministic native fixture scenario",
                quant_slots=1,
                extraction_slots=1,
            )
        )
        del envelope["planning_assumptions"]["quant_slots"]
        del envelope["planning_assumptions"]["extraction_slots"]
        (tmp_path / "snapshot.json").write_text(json.dumps(envelope))
        monkeypatch.setenv("AMALGKIT_WORKER_DEADLINE_EPOCH", "1")
        expired = run_manifest(**kwargs)
        assert expired["counts"] == {"unresolved": 1}
        assert not (work / "quant/SRR123").exists()
        monkeypatch.setenv("AMALGKIT_WORKER_DEADLINE_EPOCH", "9999999999")
    completed_samples: list[tuple[str, str]] = []

    def checkpoint_sample(species: str, accession: str) -> None:
        sample = work / "quant" / accession
        assert (sample / f"{accession}_abundance.tsv").is_file()
        assert (sample / ".metainformant_quant_provenance.json").is_file()
        active_runs = [
            directory
            for directory in (data / "acquisition_runs").iterdir()
            if not (directory / "result.json").exists()
        ]
        assert len(active_runs) == 1
        completed_samples.append((species, accession))

    first = run_manifest(**kwargs, on_quantified=checkpoint_sample)
    assert first["counts"].get("failed", 0) == 0, first["results"]
    assert first["counts"]["newly_quantified"] == 1
    assert first["quant_slots"] == 1
    assert first["quant_threads_per_sample"] == 4
    assert first["resource_profile"]["compression_level"] == 1
    assert completed_samples == [("test_species", "SRR123")]
    abundance = next((work / "quant/SRR123").glob("*abundance.tsv"))
    before = abundance.read_bytes(), abundance.stat().st_mtime_ns
    second = run_manifest(**kwargs, on_quantified=checkpoint_sample)
    assert second["counts"].get("failed", 0) == 0, second["results"]
    assert second["counts"]["reused"] == 1
    assert (abundance.read_bytes(), abundance.stat().st_mtime_ns) == before
    assert first["task_results_journal"] != second["task_results_journal"]
    assert Path(first["task_results_journal"]).is_file()
    assert completed_samples == [("test_species", "SRR123"), ("test_species", "SRR123")]

    summary_before_failure = (data / "cloud_worker_result.json").read_bytes()

    def refuse_checkpoint(species: str, accession: str) -> None:
        raise OSError(f"Checkpoint unavailable for {species}/{accession}")

    with pytest.raises(OSError, match="Checkpoint unavailable"):
        run_manifest(**kwargs, on_quantified=refuse_checkpoint)
    assert (data / "cloud_worker_result.json").read_bytes() == summary_before_failure
    assert (abundance.read_bytes(), abundance.stat().st_mtime_ns) == before
    assert reads.is_dir()

    # Keep immutable source copies intact while tampering only with the files
    # consumed by the worker. Both must fail before quant reuse or acquisition.
    snapshot_path = tmp_path / "snapshot.json"
    frozen = tmp_path / "inputs/test_species/work"
    for source in (index, metadata):
        target = frozen / source.relative_to(work)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
    snapshot = json.loads(snapshot_path.read_text())
    snapshot["input_files"] = [
        {"path": str(p.relative_to(tmp_path)), "sha256": sha256_file(p)}
        for p in (
            frozen / "index/Test_species.idx",
            frozen / "metadata/metadata_selected.tsv",
        )
    ]
    snapshot_path.write_text(json.dumps(snapshot))
    metadata_bytes = metadata.read_bytes()
    metadata.write_bytes(metadata_bytes + b"\n")
    with pytest.raises(ValueError, match="worker metadata differs"):
        run_manifest(**kwargs)
    metadata.write_bytes(metadata_bytes)
    index_bytes = index.read_bytes()
    index.write_bytes(index_bytes + b"changed-consumed-index")
    with pytest.raises(ValueError, match="worker reference differs"):
        run_manifest(**kwargs)
    index.write_bytes(index_bytes)
    assert (abundance.read_bytes(), abundance.stat().st_mtime_ns) == before
