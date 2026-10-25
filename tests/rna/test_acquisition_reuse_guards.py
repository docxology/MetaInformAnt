from __future__ import annotations

import json
import os
import threading
from pathlib import Path

import pytest

from metainformant.rna.engine.acquisition_frontier import run_acquisition_frontier
from metainformant.rna.engine.acquisition_lanes import AdmissionDecision, Result, Task
from metainformant.rna.engine.acquisition_manifest import sha256_file
from metainformant.rna.engine.acquisition_sample import reuse_manifest_quantification
from metainformant.rna.engine.provenance import write_quant_provenance


@pytest.fixture
def portable_quant_files(tmp_path: Path) -> tuple[Task, Path, Path, Path]:
    root, configs = tmp_path / "data", tmp_path / "configs"
    configs.mkdir()
    config = configs / "test.yaml"
    config.write_text("species_list: [Test_species]\n")
    sample = root / "test_species/work/quant/SRR123"
    sample.mkdir(parents=True)
    (sample / "SRR123_abundance.tsv").write_text(
        "target_id\tlength\teff_length\test_counts\ttpm\na\t100\t90\t2\t1000000\n"
    )
    (sample / "SRR123_run_info.json").write_text(json.dumps(dict(n_processed=10, n_pseudoaligned=2, n_targets=1)))
    (sample / "SRR123_abundance.h5").write_bytes(b"preserved optional checkpoint payload")
    write_quant_provenance(
        sample, species="test_species", run_accession="SRR123", config_path=config, command=["amalgkit", "quant"]
    )
    task: Task = dict(
        species="test_species",
        accession="SRR123",
        config_name=config.name,
        config_sha256=sha256_file(config),
        task_id="test_species/SRR123",
        batch_index=1,
    )
    return task, root, configs, sample


def test_regular_quant_payloads_replay_without_mutation(portable_quant_files) -> None:
    task, root, configs, sample = portable_quant_files
    before = {p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in sample.iterdir()}
    for _ in range(2):
        result = reuse_manifest_quantification(task, root, configs)
        assert result is not None and result["quantified"] is True and result["skipped"] is True
    assert before == {p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in sample.iterdir()}


@pytest.mark.parametrize(
    "payload",
    ["SRR123_abundance.tsv", "SRR123_abundance.h5", "SRR123_run_info.json", ".metainformant_quant_provenance.json"],
)
@pytest.mark.parametrize("outside_active_root", [False, True])
def test_checkpoint_payload_symlinks_are_refused_before_reuse(
    portable_quant_files, tmp_path: Path, payload: str, outside_active_root: bool
) -> None:
    task, root, configs, sample = portable_quant_files
    target = (tmp_path if outside_active_root else root) / ("retained-" + payload)
    source = sample / payload
    source.rename(target)
    source.symlink_to(target)
    before = (target.read_bytes(), target.stat().st_mtime_ns)
    assert reuse_manifest_quantification(task, root, configs) is None
    assert source.is_symlink()
    assert before == (target.read_bytes(), target.stat().st_mtime_ns)


def test_fallback_run_info_symlink_is_refused(portable_quant_files, tmp_path: Path) -> None:
    task, root, configs, sample = portable_quant_files
    target = tmp_path / "retained-info.json"
    (sample / "SRR123_run_info.json").rename(target)
    (sample / "run_info.json").symlink_to(target)
    assert reuse_manifest_quantification(task, root, configs) is None


@pytest.mark.parametrize("payload", ["manifest.json", "config.yaml", "reference_manifest.json"])
def test_restored_input_symlink_is_refused_before_content_reads(
    portable_quant_files, tmp_path: Path, payload: str
) -> None:
    task, root, configs, sample = portable_quant_files
    reference = root / "reference.json"
    reference.write_text(
        json.dumps(dict(species="test_species", status="complete", kallisto_index=str(root / "index.idx")))
    )
    write_quant_provenance(
        sample,
        species="test_species",
        run_accession="SRR123",
        config_path=configs / "test.yaml",
        reference_manifest_path=reference,
        command=["amalgkit", "quant"],
    )
    provenance = json.loads((sample / ".metainformant_quant_provenance.json").read_text())
    restored = sample / ".metainformant_restored_inputs"
    restored.mkdir()
    (restored / "config.yaml").write_bytes((configs / "test.yaml").read_bytes())
    (restored / "reference_manifest.json").write_bytes(reference.read_bytes())
    (restored / "manifest.json").write_text(
        json.dumps(
            dict(
                schema="metainformant.rna.restored_inputs.v1",
                contract_id=provenance["quant_contract_id"],
                config_sha256=provenance["config_sha256"],
                reference_manifest_sha256=provenance["reference_manifest_sha256"],
            )
        )
    )
    assert reuse_manifest_quantification(task, root, configs) is not None
    target = tmp_path / ("retained-" + payload)
    (restored / payload).rename(target)
    (restored / payload).symlink_to(target)
    assert reuse_manifest_quantification(task, root, configs) is None


@pytest.mark.parametrize("payload", ["SRR123_abundance.tsv", "SRR123_abundance.h5", "SRR123_run_info.json"])
def test_nonregular_checkpoint_payload_is_refused_without_reading(
    portable_quant_files, tmp_path: Path, payload: str
) -> None:
    task, root, configs, sample = portable_quant_files
    path = sample / payload
    path.rename(tmp_path / ("retained-" + payload))
    os.mkfifo(path)
    assert reuse_manifest_quantification(task, root, configs) is None


def test_ready_quantification_is_drained_between_serial_reuse_misses(tmp_path: Path) -> None:
    acquired, quant_started = threading.Event(), threading.Event()
    probed: list[str] = []

    def acquire(task: Task, defer: bool) -> Result:
        acquired.set()
        (tmp_path / (task["task_id"] + ".raw")).write_bytes(b"ACGT\n")
        return {"ready_for_quant": True}

    def quantify(task: Task) -> Result:
        quant_started.set()
        return {"quantified": True}

    def reuse(task: Task) -> Result | None:
        assert acquired.wait(2)
        if probed:
            assert quant_started.wait(1), "ready native future was delayed behind serial misses"
        else:
            quant_started.wait(0.04)
        probed.append(task["task_id"])
        return None

    completed: list[str] = []
    summary = run_acquisition_frontier(
        [{"task_id": str(i), "fastq_bytes": 1} for i in range(4)],
        acquisition_slots=1,
        quant_slots=1,
        recovery_slots=1,
        max_in_flight=8,
        frontier_size=8,
        raw_bytes_limit=8,
        acquire=acquire,
        quantify=quantify,
        recover=lambda task: False,
        uses_local_sra=lambda task: False,
        admission=lambda task, live: AdmissionDecision(),
        reuse=reuse,
        finished=lambda task, result: completed.append(task["task_id"]),
        transition=lambda task, stage: None,
    )
    assert sorted(completed) == ["0", "1", "2", "3"]
    assert len(probed) == len(set(probed))
    assert summary["reused_without_acquisition"] == 0
    assert summary["peak_reserved_raw_bytes"] <= 8
