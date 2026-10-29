from __future__ import annotations

import json
import threading
from pathlib import Path

import pytest

from metainformant.rna.engine.acquisition_frontier import run_acquisition_frontier
from metainformant.rna.engine.acquisition_lanes import AdmissionDecision, Result, Task


def test_cached_checkpoint_progresses_while_network_slot_is_occupied(tmp_path: Path) -> None:
    entered, cached_finished = threading.Event(), threading.Event()
    completions: list[str] = []
    transitions: list[tuple[str, str]] = []

    def acquire(task: Task, defer: bool) -> Result:
        assert task["task_id"] == "network"
        entered.set()
        assert cached_finished.wait(2), "cached checkpoint starved behind occupied acquisition slot"
        (tmp_path / "network.raw").write_bytes(b"ACGT\n")
        return {"ready_for_quant": True}

    def reuse(task: Task) -> Result | None:
        if task["task_id"] == "cached":
            assert entered.wait(2)
            return {"quantified": True, "skipped": True}
        return None

    def finish(task: Task, result: Result) -> None:
        (tmp_path / (task["task_id"] + ".json")).write_text(json.dumps(result))
        completions.append(task["task_id"])
        if task["task_id"] == "cached":
            cached_finished.set()

    summary = run_acquisition_frontier(
        [{"task_id": "network", "fastq_bytes": 10}, {"task_id": "cached", "fastq_bytes": 10}],
        acquisition_slots=1, quant_slots=1, recovery_slots=1, max_in_flight=3,
        frontier_size=3, raw_bytes_limit=30, acquire=acquire,
        quantify=lambda task: {"quantified": True}, recover=lambda task: False,
        uses_local_sra=lambda task: False, admission=lambda task, live: AdmissionDecision(),
        finished=finish, transition=lambda task, stage: transitions.append((task["task_id"], stage)),
        reuse=reuse,
    )
    assert completions == ["cached", "network"]
    assert summary["reused_without_acquisition"] == 1
    assert summary["peak_reserved_raw_bytes"] == 20
    assert ("cached", "acquiring") not in transitions


def test_reuse_miss_returns_to_bounded_acquisition_once() -> None:
    probes: list[str] = []
    acquired: list[str] = []
    results: list[Result] = []

    def reuse(task: Task) -> Result | None:
        probes.append(task["task_id"])
        return None

    def acquire(task: Task, defer: bool) -> Result:
        acquired.append(task["task_id"])
        return {"ready_for_quant": True}

    summary = run_acquisition_frontier(
        [{"task_id": str(i), "fastq_bytes": 1} for i in range(5)],
        acquisition_slots=1, quant_slots=1, recovery_slots=1, max_in_flight=2,
        frontier_size=3, raw_bytes_limit=2, acquire=acquire,
        quantify=lambda task: {"quantified": True}, recover=lambda task: False,
        uses_local_sra=lambda task: False, admission=lambda task, live: AdmissionDecision(),
        finished=lambda task, result: results.append(result), transition=lambda task, stage: None,
        reuse=reuse,
    )
    assert len(acquired) == len(set(acquired)) == len(results) == 5
    assert len(probes) == len(set(probes))
    assert summary["reused_without_acquisition"] == 0
    assert summary["peak_reserved_raw_bytes"] <= 2


def test_checkpoint_failure_in_reuse_lane_preserves_raw_reservation(tmp_path: Path) -> None:
    acquired: list[str] = []
    probed: list[str] = []
    transitions: list[tuple[str, str]] = []
    release = threading.Event()

    def acquire(task: Task, defer: bool) -> Result:
        acquired.append(task["task_id"])
        assert release.wait(2)
        (tmp_path / "network.raw").write_bytes(b"preserved")
        return {"ready_for_quant": True}

    def reuse(task: Task) -> Result:
        probed.append(task["task_id"])
        return {"quantified": True, "skipped": True}

    def finish(task: Task, result: Result) -> None:
        (tmp_path / "retained.json").write_text(json.dumps(result))
        release.set()
        raise OSError("checkpoint refused")

    with pytest.raises(OSError, match="checkpoint refused"):
        run_acquisition_frontier(
            [{"task_id": "network", "fastq_bytes": 1}, {"task_id": "cached", "fastq_bytes": 1}, {"task_id": "pending", "fastq_bytes": 1}],
            acquisition_slots=1, quant_slots=1, recovery_slots=1, max_in_flight=2,
            frontier_size=3, raw_bytes_limit=3,
            acquire=acquire,
            quantify=lambda task: {"quantified": True}, recover=lambda task: False,
            uses_local_sra=lambda task: False, admission=lambda task, live: AdmissionDecision(),
            finished=finish, transition=lambda task, stage: transitions.append((task["task_id"], stage)),
            reuse=reuse,
        )
    assert "pending" not in acquired
    assert probed == ["cached"]
    assert ("cached", "reusing") in transitions
    assert ("cached", "checkpointed") not in transitions
    assert (tmp_path / "retained.json").is_file()
    assert (tmp_path / "network.raw").read_bytes() == b"preserved"


def test_stop_does_not_probe_or_accept_planned_reuse() -> None:
    completed: list[str] = []
    summary = run_acquisition_frontier(
        [{"task_id": "cached", "fastq_bytes": 1}],
        acquisition_slots=1, quant_slots=1, recovery_slots=1, max_in_flight=1,
        frontier_size=1, raw_bytes_limit=1,
        acquire=lambda task, defer: pytest.fail("acquisition after STOP"),
        quantify=lambda task: pytest.fail("quantification after STOP"), recover=lambda task: False,
        uses_local_sra=lambda task: False, admission=lambda task, live: AdmissionDecision(),
        finished=lambda task, result: completed.append(task["task_id"]), transition=lambda task, stage: None,
        reuse=lambda task: pytest.fail("reuse after STOP"), should_stop=lambda: True,
    )
    assert completed == []
    assert summary["unsubmitted"] == 1


@pytest.mark.parametrize("result", [{"quantified": True}, {"quantified": False, "skipped": True}])
def test_reuse_requires_an_explicit_reused_success(result: Result) -> None:
    release = threading.Event()
    completed: list[str] = []

    def acquire(task: Task, defer: bool) -> Result:
        assert release.wait(2)
        return {"ready_for_quant": True}

    def reuse(task: Task) -> Result:
        release.set()
        return result

    with pytest.raises(ValueError, match="verified reused success"):
        run_acquisition_frontier(
            [{"task_id": "network", "fastq_bytes": 1}, {"task_id": "cached", "fastq_bytes": 1}],
            acquisition_slots=1, quant_slots=1, recovery_slots=1, max_in_flight=2,
            frontier_size=2, raw_bytes_limit=2, acquire=acquire,
            quantify=lambda task: {"quantified": True}, recover=lambda task: False,
            uses_local_sra=lambda task: False, admission=lambda task, live: AdmissionDecision(),
            finished=lambda task, result: completed.append(task["task_id"]), transition=lambda task, stage: None,
            reuse=reuse,
        )
    assert "cached" not in completed


def test_stop_arriving_during_reuse_probe_keeps_cached_work_unsubmitted() -> None:
    stop, release = threading.Event(), threading.Event()
    completed: list[str] = []

    def acquire(task: Task, defer: bool) -> Result:
        assert release.wait(2)
        return {"ready_for_quant": True}

    def reuse(task: Task) -> Result:
        stop.set()
        release.set()
        return {"quantified": True, "skipped": True}

    summary = run_acquisition_frontier(
        [{"task_id": "network", "fastq_bytes": 1}, {"task_id": "cached", "fastq_bytes": 1}],
        acquisition_slots=1, quant_slots=1, recovery_slots=1, max_in_flight=2,
        frontier_size=2, raw_bytes_limit=2, acquire=acquire,
        quantify=lambda task: {"quantified": True}, recover=lambda task: False,
        uses_local_sra=lambda task: False, admission=lambda task, live: AdmissionDecision(),
        finished=lambda task, result: completed.append(task["task_id"]), transition=lambda task, stage: None,
        reuse=reuse, should_stop=stop.is_set,
    )
    assert completed == ["network"]
    assert summary["unsubmitted"] == 1
    assert summary["reused_without_acquisition"] == 0


@pytest.mark.external_tool
@pytest.mark.parametrize("mutation", ["abundance", "species", "config", "symlink"])
def test_reuse_probe_rejects_corrupt_real_quant_without_changing_files(native_acquisition_factory, mutation: str, tmp_path: Path) -> None:
    from metainformant.rna.engine.acquisition_sample import reuse_manifest_quantification
    from metainformant.rna.engine.acquisition_worker import run_manifest

    fixture = native_acquisition_factory()
    run_manifest(manifest_path=fixture.manifest, data_root=fixture.data, config_dir=fixture.config_dir,
                 workers=1, threads=1, fastq_threads=1, compression_threads=1, validation_slots=1)
    task = json.loads(fixture.manifest.read_text())
    result = reuse_manifest_quantification(task, fixture.data, fixture.config_dir)
    assert result is not None and result["quantified"] is True and result["skipped"] is True
    sample = fixture.work / "quant/SRR123"
    provenance = json.loads((sample / ".metainformant_quant_provenance.json").read_text())
    abundance = sample / provenance["quantification_file"]
    match mutation:
        case "abundance":
            abundance.write_text(abundance.read_text() + "corruption\n")
        case "species":
            from metainformant.rna.engine.provenance import quantification_contract_id
            provenance["species"] = "other_species"
            provenance["quant_contract_id"] = quantification_contract_id(provenance)
            (sample / ".metainformant_quant_provenance.json").write_text(json.dumps(provenance))
        case "config":
            fixture.config.write_text(fixture.config.read_text() + "\nthreads: 99\n")
        case "symlink":
            retained = tmp_path / "retained-quant"
            sample.rename(retained)
            sample.symlink_to(retained, target_is_directory=True)
        case _:
            raise AssertionError(mutation)
    before = {p.name: p.read_bytes() for p in sample.iterdir() if p.is_file()}
    assert reuse_manifest_quantification(task, fixture.data, fixture.config_dir) is None
    assert {p.name: p.read_bytes() for p in sample.iterdir() if p.is_file()} == before
