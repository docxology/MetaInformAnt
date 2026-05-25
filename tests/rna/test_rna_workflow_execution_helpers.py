"""Tests for private RNA workflow execution helpers.

The helpers are intentionally private but covered here because they hold the
thin-orchestration behavior used by the public execute_workflow entry point.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.rna.engine.workflow import AmalgkitWorkflowConfig
from metainformant.rna.engine.workflow_execution import (
    _is_streaming_step_already_done,
    _process_streaming_sample,
    _split_streaming_steps,
    _streaming_single_sample_params,
)


def _workflow_config(tmp_path: Path) -> AmalgkitWorkflowConfig:
    work_dir = tmp_path / "work"
    work_dir.mkdir(parents=True)
    return AmalgkitWorkflowConfig(work_dir=work_dir, species_list=["Apis_mellifera"])


def test_split_streaming_steps_skips_integrate() -> None:
    """Integrate is intentionally excluded from streaming mode."""
    chunk_steps, post_steps = _split_streaming_steps(
        [
            ("getfastq", {"threads": 4}),
            ("integrate", {}),
            ("quant", {"threads": 4}),
            ("merge", {}),
        ]
    )

    assert [name for name, _ in chunk_steps] == ["getfastq", "quant"]
    assert [name for name, _ in post_steps] == ["merge"]


def test_streaming_single_sample_params_throttle_threads(tmp_path: Path) -> None:
    """Per-sample parameters should reduce thread fan-out inside a streaming pool."""
    meta_file = tmp_path / "sample.tsv"
    params = _streaming_single_sample_params({"threads": 12, "jobs": 12}, meta_file, chunk_size=5)

    assert params["metadata"] == str(meta_file)
    assert params["threads"] == 2
    assert params["jobs"] == 2


def test_streaming_single_sample_params_rejects_non_positive_chunk_size(tmp_path: Path) -> None:
    """Streaming helper should reject chunk sizes that cannot divide worker resources."""
    with pytest.raises(ValueError, match="chunk_size must be >= 1"):
        _streaming_single_sample_params({"threads": 12}, tmp_path / "sample.tsv", chunk_size=0)


def test_streaming_step_already_done_detects_quant_output(tmp_path: Path) -> None:
    """Existing quant outputs bypass expensive per-sample subprocess startup."""
    config = _workflow_config(tmp_path)
    sample_id = "SRR000001"
    quant_dir = config.work_dir / "quant" / sample_id
    quant_dir.mkdir(parents=True)
    (quant_dir / f"{sample_id}_abundance.tsv").write_text("target_id\ttpm\nGENE1\t1.0\n")

    assert _is_streaming_step_already_done(config, sample_id, "quant", {"redo": "no"})
    assert not _is_streaming_step_already_done(config, sample_id, "quant", {"redo": "yes"})


def test_process_streaming_sample_walk_writes_single_metadata(tmp_path: Path) -> None:
    """Dry-run sample processing should create one-row metadata and planned step results."""
    config = _workflow_config(tmp_path)
    sample_row = {"run": "SRR000001", "organism": "Apis mellifera"}

    results = _process_streaming_sample(
        config=config,
        sample_row=sample_row,
        fieldnames=["run", "organism"],
        current_chunk_idx=0,
        sample_idx_in_chunk=0,
        chunk_steps=[("getfastq", {"threads": 4}), ("quant", {"threads": 4})],
        chunk_size=2,
        step_functions={},
        walk=True,
    )

    assert [result.step_name for result in results] == ["getfastq_SRR000001", "quant_SRR000001"]
    assert all(result.success for result in results)
    metadata_file = config.work_dir / "metadata" / "metadata_chunk_0_sample_SRR000001.tsv"
    assert metadata_file.exists()
    assert "SRR000001" in metadata_file.read_text()
