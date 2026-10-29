from __future__ import annotations

import json
import threading
from pathlib import Path

import pytest

from metainformant.rna.engine.acquisition_lanes import AdmissionDecision
from metainformant.rna.engine.acquisition_frontier import run_acquisition_frontier


def test_download_continues_during_quantification_and_retains_reservations(tmp_path: Path) -> None:
    quant_entered, second_downloaded = (threading.Event() for _ in range(2))
    events: list[tuple[str, str]] = []
    completed: list[str] = []
    tasks = [{"task_id": str(i), "fastq_bytes": 10} for i in range(4)]

    def acquire(task: dict, defer: bool) -> dict:
        identity = task["task_id"]
        (tmp_path / (identity + ".raw")).write_bytes(b"ACGT\n" * 2)
        if identity == "1":
            assert quant_entered.wait(10)
            second_downloaded.set()
        return {"ready_for_quant": True}

    def quantify(task: dict) -> dict:
        if task["task_id"] == "0":
            quant_entered.set()
            assert second_downloaded.wait(10)
        return {"quantified": True}

    def finish(task: dict, result: dict) -> None:
        assert result["quantified"]
        identity = task["task_id"]
        (tmp_path / (identity + ".proof")).write_text(json.dumps(result))
        completed.append(identity)

    summary = run_acquisition_frontier(
        tasks, acquisition_slots=1, quant_slots=1, recovery_slots=1,
        max_in_flight=4, frontier_size=48, raw_bytes_limit=20,
        acquire=acquire, quantify=quantify, recover=lambda task: True,
        uses_local_sra=lambda task: False,
        admission=lambda task, live: AdmissionDecision(), finished=finish,
        transition=lambda task, stage: events.append((task["task_id"], stage)),
    )
    assert sorted(completed) == ["0", "1", "2", "3"]
    assert summary["peak_reserved_raw_bytes"] == 20
    assert summary["peak_admitted"] == 2
    assert events.index(("1", "acquiring")) < events.index(("0", "checkpointed"))
    assert all((tmp_path / (str(i) + ".proof")).is_file() for i in range(4))


def test_stop_leaves_planned_items_unsubmitted_and_drains_acquired_inputs(tmp_path: Path) -> None:
    stop = threading.Event()
    completed: list[str] = []

    def acquire(task: dict, defer: bool) -> dict:
        stop.set()
        return {"ready_for_quant": True}

    summary = run_acquisition_frontier(
        [{"task_id": str(i), "fastq_bytes": 1} for i in range(48)],
        acquisition_slots=1, quant_slots=1, recovery_slots=1,
        max_in_flight=48, frontier_size=48, raw_bytes_limit=48,
        acquire=acquire, quantify=lambda task: {"quantified": True},
        recover=lambda task: True, uses_local_sra=lambda task: False,
        admission=lambda task, live: AdmissionDecision(),
        finished=lambda task, result: completed.append(task["task_id"]),
        transition=lambda task, stage: (tmp_path / (task["task_id"] + ".stage")).write_text(stage),
        should_stop=stop.is_set,
    )
    assert completed == ["0"]
    assert summary["unsubmitted"] == 47
    assert summary["peak_frontier"] == 48


@pytest.mark.parametrize("weight", [0, 21])
def test_unknown_or_oversized_reservation_is_explicitly_refused(weight: int) -> None:
    results: list[dict] = []
    summary = run_acquisition_frontier(
        [{"task_id": "held", "fastq_bytes": weight}],
        acquisition_slots=1, quant_slots=1, recovery_slots=1,
        max_in_flight=1, frontier_size=48, raw_bytes_limit=20,
        acquire=lambda task, defer: pytest.fail("Unreserved acquisition ran"),
        quantify=lambda task: pytest.fail("Unreserved quantification ran"),
        recover=lambda task: False, uses_local_sra=lambda task: False,
        admission=lambda task, live: AdmissionDecision(),
        finished=lambda task, result: results.append(result), transition=lambda task, stage: None,
    )
    assert summary["peak_admitted"] == 0
    assert results[0]["status"] == "unresolved_deadline"
    assert results[0]["quantified"] is False
    assert "reservation" in results[0]["error"]


def test_failed_checkpoint_prevents_replacement_acquisition(tmp_path: Path) -> None:
    acquired: list[str] = []

    def acquire(task: dict, defer: bool) -> dict:
        acquired.append(task["task_id"])
        (tmp_path / task["task_id"]).write_bytes(b"preserved")
        return {"ready_for_quant": True}

    def finish(task: dict, result: dict) -> None:
        raise OSError("checkpoint readback failed")

    with pytest.raises(OSError, match="checkpoint readback"):
        run_acquisition_frontier(
            [{"task_id": str(i), "fastq_bytes": 10} for i in range(3)],
            acquisition_slots=1, quant_slots=1, recovery_slots=1,
            max_in_flight=1, frontier_size=48, raw_bytes_limit=10,
            acquire=acquire, quantify=lambda task: {"quantified": True},
            recover=lambda task: False, uses_local_sra=lambda task: False,
            admission=lambda task, live: AdmissionDecision(), finished=finish,
            transition=lambda task, stage: None,
        )
    assert acquired == ["0"]
    assert (tmp_path / "0").read_bytes() == b"preserved"


@pytest.mark.parametrize("recovered", [False, True])
def test_recovery_has_one_attempt_and_keeps_other_downloads_independent(recovered: bool) -> None:
    attempts: list[tuple[str, bool]] = []
    results: list[tuple[str, dict]] = []

    def acquire(task: dict, defer: bool) -> dict:
        attempts.append((task["task_id"], defer))
        return {"fallback_deferred": True} if defer and task["task_id"] == "0" else {"ready_for_quant": True}

    summary = run_acquisition_frontier(
        [{"task_id": str(i), "fastq_bytes": 10} for i in range(2)],
        acquisition_slots=1, quant_slots=1, recovery_slots=1,
        max_in_flight=2, frontier_size=48, raw_bytes_limit=20,
        acquire=acquire, quantify=lambda task: {"quantified": True},
        recover=lambda task: recovered, uses_local_sra=lambda task: False,
        admission=lambda task, live: AdmissionDecision(),
        finished=lambda task, result: results.append((task["task_id"], result)),
        transition=lambda task, stage: None,
    )
    assert len(results) == len({identity for identity, _ in results}) == 2
    assert dict(results)["1"]["quantified"] is True
    assert dict(results)["0"]["quantified"] is recovered
    assert attempts.count(("0", False)) == int(recovered)
    assert summary["deferred"] == 1


def test_busy_local_recovery_does_not_block_later_ena_acquisition() -> None:
    ena_entered = threading.Event()
    completed: list[str] = []

    def acquire(task: dict, defer: bool) -> dict:
        if task["task_id"] == "2":
            ena_entered.set()
        return {"ready_for_quant": True}

    def recover(task: dict) -> bool:
        assert ena_entered.wait(10)
        return True

    run_acquisition_frontier(
        [{"task_id": str(i), "fastq_bytes": 10} for i in range(3)],
        acquisition_slots=1, quant_slots=1, recovery_slots=1,
        max_in_flight=3, frontier_size=48, raw_bytes_limit=30,
        acquire=acquire, quantify=lambda task: {"quantified": True}, recover=recover,
        uses_local_sra=lambda task: task["task_id"] != "2",
        admission=lambda task, live: AdmissionDecision(),
        finished=lambda task, result: completed.append(task["task_id"]),
        transition=lambda task, stage: None,
    )
    assert sorted(completed) == ["0", "1", "2"]


def test_refused_frontier_refills_to_reach_later_supported_work() -> None:
    results: list[tuple[str, dict]] = []
    tasks = [{"task_id": str(i), "fastq_bytes": 0} for i in range(50)]
    tasks.append({"task_id": "supported", "fastq_bytes": 1})
    summary = run_acquisition_frontier(
        tasks, acquisition_slots=1, quant_slots=1, recovery_slots=1,
        max_in_flight=1, frontier_size=48, raw_bytes_limit=20,
        acquire=lambda task, defer: {"ready_for_quant": True},
        quantify=lambda task: {"quantified": True}, recover=lambda task: False,
        uses_local_sra=lambda task: False, admission=lambda task, live: AdmissionDecision(),
        finished=lambda task, result: results.append((task["task_id"], result)),
        transition=lambda task, stage: None,
    )
    assert summary["unsubmitted"] == 0
    assert len(results) == 51
    assert dict(results)["supported"]["quantified"] is True


def test_quantification_stage_refuses_changed_mates_without_downloading(
    native_acquisition_factory, monkeypatch: pytest.MonkeyPatch
) -> None:
    import gzip

    from metainformant.rna.engine.streaming_orchestrator import StreamingPipelineOrchestrator

    fixture = native_acquisition_factory(paired_metadata=True)
    raw = fixture.work / "getfastq/SRR123/SRR123.fastq.gz"
    raw.write_bytes(gzip.compress(b"@wrongA/1\nACGT\n+\n!!!!\n@wrongB/2\nACGT\n+\n!!!!\n", mtime=0))
    monkeypatch.setenv("AMALGKIT_PIPELINE_VALIDATE_FASTQ_RECORDS", "yes")
    owner = StreamingPipelineOrchestrator(config_dir=fixture.config_dir, db_path=fixture.data / "progress.db")
    try:
        result = owner.process_single_sample("SRR123", 1, raw.parent.parent, fixture.config,
                                             "test_species", 1, True, quantification_only=True)
        assert result["quantified"] is False
        assert not (fixture.work / "quant/SRR123").exists()
        assert raw.is_file()
    finally:
        owner.db.close()


@pytest.mark.parametrize("mutation", ["replace-inode", "mutate-in-place", "matching-witness"])
def test_quantification_refuses_same_size_mtime_mate_mutation(
    native_acquisition_factory, monkeypatch: pytest.MonkeyPatch, tmp_path: Path, mutation: str,
) -> None:
    import gzip
    import os

    from metainformant.rna.engine.streaming_orchestrator import (
        StreamingPipelineOrchestrator,
        _raw_validation_marker,
        _write_raw_validation_marker,
        _validated_local_fastq_inputs,
    )

    fixture = native_acquisition_factory(paired_metadata=True)
    sample = fixture.work / "getfastq/SRR123"
    (sample / "SRR123.fastq.gz").rename(tmp_path / "original-single.fastq.gz")
    sequences = [line for line in (tmp_path / "transcripts.fa").read_text().splitlines()
                 if not line.startswith(">")]
    for mate in (1, 2):
        records = []
        for n in range(90):
            sequence, start = sequences[n % 3], 30 + n
            read = (sequence[start:start + 75] if mate == 1 else
                    sequence[start + 125:start + 200].translate(str.maketrans("ACGT", "TGCA"))[::-1])
            records.append(f"@SRR123.{n:03d}/{mate}\n{read}\n+\n{'I' * 75}\n")
        (sample / f"SRR123_{mate}.fastq.gz").write_bytes(gzip.compress("".join(records).encode(), compresslevel=0, mtime=0))
    monkeypatch.setenv("AMALGKIT_PIPELINE_VALIDATE_FASTQ_RECORDS", "yes")
    owner = StreamingPipelineOrchestrator(config_dir=fixture.config_dir, db_path=fixture.data / "progress.db")
    try:
        owner.db.init_species("test_species", ["SRR123"])
        acquired = owner.process_single_sample("SRR123", 1, sample.parent, fixture.config,
                                               "test_species", 2, True, acquisition_only=True)
        assert acquired["ready_for_quant"]
        mate = sample / "SRR123_2.fastq.gz"
        before = mate.stat()
        changed = gzip.compress(gzip.decompress(mate.read_bytes()).replace(b"SRR123", b"SRR999"),
                                compresslevel=0, mtime=0)
        assert len(changed) == before.st_size
        if mutation == "replace-inode":
            temporary = tmp_path / "replacement.fastq.gz"
            temporary.write_bytes(changed)
            os.replace(temporary, mate)
        else:
            mate.write_bytes(changed)
        os.utime(mate, ns=(before.st_atime_ns, before.st_mtime_ns))
        assert mate.stat().st_size == before.st_size and mate.stat().st_mtime_ns == before.st_mtime_ns
        assert (mate.stat().st_ino != before.st_ino) is (mutation == "replace-inode")
        if mutation == "matching-witness":
            _write_raw_validation_marker(_raw_validation_marker(sample, "SRR123"),
                                         [sample / "SRR123_1.fastq.gz", mate], True,
                                         layout_validation="full_fastq_record_scan")
            assert len(_validated_local_fastq_inputs(sample, "SRR123", True, require_records=True)) == 2
        result = owner.process_single_sample("SRR123", 1, sample.parent, fixture.config,
                                            "test_species", 2, True, quantification_only=True)
        assert result["quantified"] is False
        assert "invalid" in result["error"]
        assert not list((fixture.work / "quant").rglob("*abundance.tsv"))
        assert _validated_local_fastq_inputs(sample, "SRR123", True, require_records=True,
                                            fresh_record_scan=True) == []
        if mutation != "matching-witness":
            assert _validated_local_fastq_inputs(sample, "SRR123", True, require_records=True) == []
        assert mate.read_bytes() == changed
    finally:
        owner.db.close()


@pytest.mark.parametrize("paired", [False, True])
def test_real_frontier_quantifies_and_replays_without_changing_raw_or_outputs(
    native_acquisition_factory, monkeypatch: pytest.MonkeyPatch, paired: bool
) -> None:
    from metainformant.rna.engine.acquisition_manifest import sha256_file
    from metainformant.rna.engine.acquisition_worker import run_manifest
    from tests.rna.test_ena_staging import transfer_site

    fixture = native_acquisition_factory(accessions=("SRR123", "SRR124", "SRR125"), paired_metadata=paired)
    for key, value in {
        "AMALGKIT_LOCAL_FRONTIER_SIZE": "48", "AMALGKIT_LOCAL_RAW_RESERVATION_BYTES": "45000",
        "AMALGKIT_PIPELINE_ACQUISITION_SLOTS": "1", "AMALGKIT_PIPELINE_DEDICATED_RECOVERY": "yes",
        "AMALGKIT_PIPELINE_VALIDATE_FASTQ_RECORDS": "yes",
    }.items():
        monkeypatch.setenv(key, value)
    kwargs = dict(manifest_path=fixture.manifest, data_root=fixture.data, config_dir=fixture.config_dir,
                  workers=2, threads=2, fastq_threads=1, compression_threads=1,
                  validation_slots=1, quant_slots=1, fasterq_slots=1, max_in_flight=3)
    raw = fixture.work / "getfastq/SRR123/SRR123.fastq.gz"
    payload = raw.read_bytes()
    raw.rename(raw.parent / "preserved-source.bin")
    with transfer_site({"SRR123.fastq.gz": payload}) as site:
        monkeypatch.setenv("AMALGKIT_ENA_API_BASE_URL", f"{site.base}/api")
        first = run_manifest(**kwargs)
        assert any(method == "GET" and path == "/SRR123.fastq.gz" for method, path, _ in site.requests)
    assert first["counts"]["newly_quantified"] == 3
    files = [*fixture.work.glob("getfastq/*/*.gz"), *fixture.work.glob("quant/*/*abundance.tsv")]
    assert len(files) == 6
    before = {str(path): (sha256_file(path), path.stat().st_mtime_ns) for path in files}
    second = run_manifest(**kwargs)
    assert second["counts"]["reused"] == 3
    assert before == {str(path): (sha256_file(path), path.stat().st_mtime_ns) for path in files}
    journals = list((fixture.data / "acquisition_runs").glob("*/frontier_events.jsonl"))
    assert len(journals) == 2
    stages = [[json.loads(line)["stage"] for line in path.read_text().splitlines()] for path in journals]
    assert sorted(entries.count("ready_for_quant") for entries in stages) == [0, 3]
    assert all(entries.count("checkpointed") == 3 for entries in stages)


@pytest.mark.parametrize("quant_slots", [2, 4])
def test_eight_thread_compute_profiles_run_native_quantification_and_reuse(
    native_acquisition_factory, monkeypatch: pytest.MonkeyPatch, quant_slots: int
) -> None:
    from metainformant.rna.engine.acquisition_worker import run_manifest

    fixture = native_acquisition_factory(accessions=("SRR123", "SRR124", "SRR125", "SRR126"))
    for key, value in {
        "AMALGKIT_LOCAL_FRONTIER_SIZE": "48", "AMALGKIT_LOCAL_RAW_RESERVATION_BYTES": "60000",
        "AMALGKIT_PIPELINE_ACQUISITION_SLOTS": "2", "AMALGKIT_PIPELINE_DEDICATED_RECOVERY": "yes",
        "AMALGKIT_PIPELINE_VALIDATE_FASTQ_RECORDS": "yes",
    }.items():
        monkeypatch.setenv(key, value)
    kwargs = dict(manifest_path=fixture.manifest, data_root=fixture.data, config_dir=fixture.config_dir,
                  workers=4, threads=8, fastq_threads=1, compression_threads=1,
                  validation_slots=1, quant_slots=quant_slots, fasterq_slots=1, max_in_flight=4)
    first = run_manifest(**kwargs)
    assert first["counts"]["newly_quantified"] == 4
    assert all(entry["quantified"] for entry in first["results"])
    assert run_manifest(**kwargs)["counts"]["reused"] == 4
