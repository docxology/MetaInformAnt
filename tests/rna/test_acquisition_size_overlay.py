from __future__ import annotations

import copy
import hashlib
import json
import os
import threading
from pathlib import Path

import pytest

from metainformant.rna.engine.acquisition_frontier import run_acquisition_frontier
from metainformant.rna.engine.acquisition_lanes import AdmissionDecision
from metainformant.rna.engine.acquisition_order import plan_acquisition_order
from metainformant.rna.engine.acquisition_size_overlay import (
    effective_fastq_bytes, load_size_overlay,
)


def fixture_files(tmp_path: Path) -> tuple[Path, list[dict], Path, dict]:
    tasks = [
        {"task_id": "species_b/SRR1", "species": "species_b", "accession": "SRR1",
         "fastq_bytes": 0, "expected_paired": True, "reference_index_sha256": "a" * 64},
        {"task_id": "species_a/SRR2", "species": "species_a", "accession": "SRR2",
         "fastq_bytes": 90, "expected_paired": False, "reference_index_sha256": "b" * 64},
    ]
    manifest = tmp_path / "manifest.jsonl"
    manifest.write_text("".join(json.dumps(t) + "\n" for t in tasks))
    document = {"schema": "metainformant.rna.size_overlay.v1",
                "frozen_manifest_sha256": hashlib.sha256(manifest.read_bytes()).hexdigest(),
                "records": [{"task_id": "species_b/SRR1", "accession": "SRR1", "expected_paired": True,
                             "reference_index_sha256": "a" * 64, "fastq_bytes": 30,
                             "files": [{"name": "SRR1_1.fastq.gz", "url": "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR1/SRR1_1.fastq.gz", "bytes": 10, "md5": "1" * 32},
                                       {"name": "SRR1_2.fastq.gz", "url": "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR1/SRR1_2.fastq.gz", "bytes": 20, "md5": "2" * 32}]}]}
    path = tmp_path / "size-overlay.json"
    path.write_text(json.dumps(document))
    return manifest, tasks, path, document


def loaded(manifest: Path, tasks: list[dict], path: Path):
    return load_size_overlay(manifest, tasks, path=str(path),
                             expected_sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def test_overlay_is_immutable_and_preserves_frozen_identity(tmp_path: Path) -> None:
    manifest, tasks, path, _ = fixture_files(tmp_path)
    original = manifest.read_bytes(); before = copy.deepcopy(tasks)
    overlay = loaded(manifest, tasks, path)
    assert overlay.sizes == {"species_b/SRR1": 30}
    assert effective_fastq_bytes(tasks[0], overlay.sizes) == 30
    assert effective_fastq_bytes(tasks[1], overlay.sizes) == 90
    assert manifest.read_bytes() == original and tasks == before
    assert overlay.receipt()["size_overlay_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    with pytest.raises(TypeError):
        overlay.sizes["species_b/SRR1"] = 31
    assert load_size_overlay(manifest, tasks, path=None, expected_sha256=None).receipt() == {}


@pytest.mark.parametrize("case", [
    "schema", "base", "duplicate", "unknown", "positive_overwrite", "run", "reference", "layout",
    "bool_size", "float_size", "negative", "zero", "sum", "missing_mate", "duplicate_mate", "md5",
    "foreign_url", "url_user", "traversal", "extra_record_field", "extra_root_field",
])
def test_bad_overlay_refuses_real_file(tmp_path: Path, case: str) -> None:
    manifest, tasks, path, document = fixture_files(tmp_path)
    row = document["records"][0]; file = row["files"][0]
    if case == "schema": document["schema"] = "unbound"
    elif case == "base": document["frozen_manifest_sha256"] = "0" * 64
    elif case == "duplicate": document["records"].append(copy.deepcopy(row))
    elif case == "unknown": row["task_id"] = "species_b/SRR9"
    elif case == "positive_overwrite": tasks[0]["fastq_bytes"] = 1
    elif case == "run": row["accession"] = "SRR9"
    elif case == "reference": row["reference_index_sha256"] = "0" * 64
    elif case == "layout": row["expected_paired"] = False
    elif case == "bool_size": row["fastq_bytes"] = True
    elif case == "float_size": row["fastq_bytes"] = 30.0
    elif case == "negative": row["fastq_bytes"] = -1
    elif case == "zero": row["fastq_bytes"] = 0
    elif case == "sum": row["fastq_bytes"] = 31
    elif case == "missing_mate": row["files"].pop()
    elif case == "duplicate_mate": row["files"][1] = copy.deepcopy(file)
    elif case == "md5": file["md5"] = "bad"
    elif case == "foreign_url": file["url"] = file["url"].replace("ftp.sra.ebi.ac.uk", "example.invalid")
    elif case == "url_user": file["url"] = file["url"].replace("https://", "https://user@")
    elif case == "traversal": file["url"] = file["url"].replace("/SRR1/", "/%2e%2e/SRR1/")
    elif case == "extra_record_field": row["trust_me"] = True
    elif case == "extra_root_field": document["trust_me"] = True
    else: raise AssertionError(case)
    path.write_text(json.dumps(document))
    with pytest.raises(ValueError):
        loaded(manifest, tasks, path)


@pytest.mark.parametrize("case", ["missing_hash", "missing_path", "wrong_hash", "relative", "parent", "file_symlink", "directory_symlink", "directory", "fifo", "empty", "large", "duplicate_json"])
def test_bad_overlay_path_or_bytes_refuses(tmp_path: Path, case: str) -> None:
    manifest, tasks, path, _ = fixture_files(tmp_path)
    digest = hashlib.sha256(path.read_bytes()).hexdigest(); supplied: str | None = str(path)
    if case == "missing_hash": digest = None
    elif case == "missing_path": supplied = None
    elif case == "wrong_hash": digest = "0" * 64
    elif case == "relative": supplied = path.name
    elif case == "parent": supplied = str(path.parent / ".." / path.parent.name / path.name)
    elif case == "file_symlink":
        link = tmp_path / "link.json"; link.symlink_to(path); supplied = str(link)
    elif case == "directory_symlink":
        link = tmp_path / "linkdir"; link.symlink_to(tmp_path, target_is_directory=True); supplied = str(link / path.name)
    elif case == "directory": supplied = str(tmp_path)
    elif case == "fifo":
        fifo = tmp_path / "fifo"; os.mkfifo(fifo); supplied = str(fifo)
    elif case == "empty": path.write_bytes(b"")
    elif case == "large": path.write_bytes(b"x" * (1024 * 1024 + 1))
    elif case == "duplicate_json":
        path.write_text('{"schema":"a","schema":"b"}')
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
    else: raise AssertionError(case)
    with pytest.raises(ValueError):
        load_size_overlay(manifest, tasks, path=supplied, expected_sha256=digest)


def test_order_and_reservation_share_overlay_without_rewriting_tasks(tmp_path: Path) -> None:
    manifest, tasks, path, _ = fixture_files(tmp_path); overlay = loaded(manifest, tasks, path)
    original = copy.deepcopy(tasks); outputs = []; acquired = []
    order = plan_acquisition_order(tasks, ["checkpoint/old"], global_size_order=True, size_overrides=overlay.sizes)
    assert order == ["checkpoint/old", "species_b/SRR1", "species_a/SRR2"]
    def acquire(task: dict, defer: bool) -> dict:
        acquired.append(task["task_id"])
        target = tmp_path / (task["accession"] + ".raw")
        target.write_bytes(b"x" * effective_fastq_bytes(task, overlay.sizes))
        return {"ready_for_quant": True}
    def quantify(task: dict) -> dict:
        assert (tmp_path / (task["accession"] + ".raw")).is_file()
        return {"quantified": True}
    summary = run_acquisition_frontier(tasks, acquisition_slots=2, quant_slots=2, recovery_slots=1,
        max_in_flight=2, frontier_size=2, raw_bytes_limit=100, acquire=acquire, quantify=quantify,
        recover=lambda t: False, uses_local_sra=lambda t: False,
        admission=lambda t, live: AdmissionDecision(), finished=lambda t, result: outputs.append(result),
        transition=lambda t, stage: None, size_overrides=overlay.sizes)
    assert len(acquired) == len(outputs) == 2 and all(r["quantified"] for r in outputs)
    assert summary["peak_reserved_raw_bytes"] == 90 and tasks == original


def test_zero_without_overlay_remains_refused(tmp_path: Path) -> None:
    _, tasks, _, _ = fixture_files(tmp_path); outputs = []; dispatched = []
    run_acquisition_frontier(tasks[:1], acquisition_slots=1, quant_slots=1, recovery_slots=1,
        max_in_flight=1, frontier_size=1, raw_bytes_limit=100,
        acquire=lambda t, d: dispatched.append(t), quantify=lambda t: {}, recover=lambda t: False,
        uses_local_sra=lambda t: False, admission=lambda t, live: AdmissionDecision(),
        finished=lambda t, result: outputs.append(result), transition=lambda t, stage: None)
    assert not dispatched and outputs[0]["status"] == "unresolved_deadline"


def test_overlay_reuse_probe_releases_its_actual_reservation(tmp_path: Path) -> None:
    manifest, tasks, path, _ = fixture_files(tmp_path); overlay = loaded(manifest, tasks, path)
    ready = threading.Event(); done = []
    ordered = [tasks[1], tasks[0]]
    def acquire(task: dict, defer: bool) -> dict:
        assert ready.wait(10)
        return {"ready_for_quant": True}
    def reuse(task: dict) -> dict:
        ready.set()
        return {"quantified": True, "skipped": True}
    result = run_acquisition_frontier(ordered, acquisition_slots=1, quant_slots=1, recovery_slots=1,
        max_in_flight=2, frontier_size=2, raw_bytes_limit=120, acquire=acquire,
        quantify=lambda t: {"quantified": True}, recover=lambda t: False,
        uses_local_sra=lambda t: False, admission=lambda t, live: AdmissionDecision(),
        finished=lambda t, value: done.append(value), transition=lambda t, stage: None,
        reuse=reuse, size_overrides=overlay.sizes)
    assert result["reused_without_acquisition"] == 1 and result["peak_reserved_raw_bytes"] == 120
    assert len(done) == 2 and sum(x.get("skipped", False) for x in done) == 1
