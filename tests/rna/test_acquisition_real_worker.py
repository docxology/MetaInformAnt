"""Actual Kallisto/Amalgkit execution and portable idempotent manifest replay."""

from __future__ import annotations

import csv
import gzip
import json
import os
import random
import shutil
import subprocess
import sys
from dataclasses import asdict
from pathlib import Path

import pytest

from metainformant.rna.engine.acquisition_manifest import sha256_file
from metainformant.rna.engine.acquisition_scheduling import PlanningAssumptions
from metainformant.rna.engine.acquisition_worker import run_manifest


@pytest.mark.external_tool
@pytest.mark.parametrize("planning_profile", [False, True])
@pytest.mark.parametrize("reference_target", ["Test species", "Test species subspecies"])
def test_real_manifest_quantification_and_second_run_reuses_bytes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, reference_target: str, planning_profile: bool
) -> None:
    if not shutil.which("kallisto") or not shutil.which("amalgkit"):
        pytest.skip("real Kallisto and Amalgkit required")
    for name in (
        "AMALGKIT_DURABLE_BUCKET",
        "AMALGKIT_DURABLE_COHORT",
        "AMALGKIT_CLOUD_MAX_RAW_BYTES",
        "AMALGKIT_WORKER_DEADLINE_EPOCH",
    ):
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setenv("AMALGKIT_MIN_EXTERNAL_FREE_GB", "0")
    monkeypatch.setenv("AMALGKIT_MIN_SYSTEM_FREE_GB", "0")
    monkeypatch.setenv("AMALGKIT_RECLAIM_RAW_AFTER_QUANT", "no")
    data = tmp_path / "data"
    # run_manifest sets these process-wide values; register them for teardown
    # before exercising native tools so later resource tests keep their defaults.
    monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(data))
    monkeypatch.setenv("AMALGKIT_PIPELINE_FASTQ_THREADS", "1")
    monkeypatch.setenv("AMALGKIT_PIPELINE_COMPRESSION_THREADS", "1")
    monkeypatch.setenv("AMALGKIT_PIPELINE_VALIDATION_SLOTS", "1")
    monkeypatch.setenv("AMALGKIT_PIPELINE_COMPRESSION_LEVEL", "1")
    work = data / "test_species/work"
    metadata_dir = work / "metadata"
    index_dir = work / "index"
    metadata_dir.mkdir(parents=True)
    index_dir.mkdir()
    rng = random.Random(161)
    sequences = ["".join(rng.choice("ACGT") for _ in range(500)) for _ in range(3)]
    fasta = tmp_path / "transcripts.fa"
    fasta.write_text("".join(f">gene{i}\n{sequence}\n" for i, sequence in enumerate(sequences)))
    index = index_dir / "Test_species.idx"
    subprocess.run(
        ["kallisto", "index", "-i", str(index), str(fasta)],
        check=True,
        capture_output=True,
        text=True,
    )
    reads = work / "getfastq/SRR123"
    reads.mkdir(parents=True)
    with gzip.open(reads / "SRR123.fastq.gz", "wt") as handle:
        for i in range(90):
            sequence = sequences[i % 3][100:175]
            handle.write(f"@SRR123.{i + 1}\n{sequence}\n+\n{'I' * 75}\n")
    metadata = metadata_dir / "metadata_selected.tsv"
    row = {
        "run": "SRR123",
        "scientific_name": reference_target,
        "lib_layout": "single",
        "total_spots": "90",
        "total_bases": "6750",
        "spot_length": "75",
        "size": "15000",
        "sample_group": "test",
        "is_sampled": "yes",
        "is_qualified": "yes",
        "exclusion": "no",
    }
    with metadata.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)
    config_dir = tmp_path / "config"
    config_dir.mkdir()
    config = config_dir / "amalgkit_test_species.yaml"
    config.write_text(
        "reference_aliases:\n  Test species subspecies: Test_species\nspecies_list: [Test_species]\nsteps:\n  quant:\n    index_dir: output/amalgkit/test_species/work/index\n"
    )
    manifest = tmp_path / "manifest.jsonl"
    manifest.write_text(
        json.dumps(
            {
                "schema": "metainformant.rna.acquisition_task.v1",
                "task_id": "test_species/SRR123",
                "species": "test_species",
                "accession": "SRR123",
                "batch_index": 1,
                "config_name": config.name,
                "config_sha256": sha256_file(config),
                "reference_index_sha256": sha256_file(index),
                "expected_paired": False,
                "fastq_bytes": 15000,
            }
        )
        + "\n"
    )
    (tmp_path / "snapshot.json").write_text(
        json.dumps(
            {
                "schema": "metainformant.rna.acquisition_snapshot.v1",
                "cloud_launch_policy": "checkpointed",
                "task_count": 1,
                "manifest_sha256": sha256_file(manifest),
                "input_files": [
                    {"path": str(p.relative_to(tmp_path)), "sha256": sha256_file(p)} for p in (index, metadata)
                ],
            }
        )
    )
    kwargs = dict(
        manifest_path=manifest,
        data_root=data,
        config_dir=config_dir,
        workers=2,
        threads=2,
        fastq_threads=1,
        compression_threads=1,
        validation_slots=1,
        quant_slots=1,
        fasterq_slots=1,
        max_in_flight=1,
    )
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
        environment = dict(os.environ, AMALGKIT_WORKER_DEADLINE_EPOCH=deadline, AMALGKIT_CLOUD_MAX_RAW_BYTES="15000")
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
        (tmp_path / f"{scenario}-cli.txt").write_text(completed.stdout + completed.stderr)
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
        (tmp_path / "snapshot.json").write_text(json.dumps(envelope))
        monkeypatch.setenv("AMALGKIT_WORKER_DEADLINE_EPOCH", "1")
        expired = run_manifest(**kwargs)
        assert expired["counts"] == {"unresolved": 1}
        assert not (work / "quant/SRR123").exists()
        monkeypatch.setenv("AMALGKIT_WORKER_DEADLINE_EPOCH", "9999999999")
    first = run_manifest(**kwargs)
    assert first["counts"].get("failed", 0) == 0, first["results"]
    assert first["counts"]["newly_quantified"] == 1
    assert first["resource_profile"]["compression_level"] == 1
    abundance = next((work / "quant/SRR123").glob("*abundance.tsv"))
    before = abundance.read_bytes(), abundance.stat().st_mtime_ns
    second = run_manifest(**kwargs)
    assert second["counts"].get("failed", 0) == 0, second["results"]
    assert second["counts"]["reused"] == 1
    assert (abundance.read_bytes(), abundance.stat().st_mtime_ns) == before
    assert first["task_results_journal"] != second["task_results_journal"]
    assert Path(first["task_results_journal"]).is_file()

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
