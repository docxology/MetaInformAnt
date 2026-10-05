"""Real budget and partition counterexamples; no AWS resources are created."""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from metainformant.rna.engine.aws_completion import (
    _render_startup,
    budget_allows,
    choose_partition,
    job_timeout,
    verify_locked_campaign,
)


def test_budget_reserves_entire_job_and_storage() -> None:
    assert budget_allows(67.92, 750, 15300, 2)
    assert not budget_allows(735, 750, 15300, 2)
    assert not budget_allows(740, 750, 3600, 2)


def test_large_run_receives_a_transfer_sized_deadline() -> None:
    assert job_timeout(100 * 1024**3, 14400) == 208400
    assert job_timeout(1024**2, 14400) == 14400
    with pytest.raises(ValueError):
        job_timeout(0, 14400)


@pytest.mark.parametrize("value", [math.nan, math.inf, -1])
def test_invalid_budget_fails_closed(value: float) -> None:
    with pytest.raises(ValueError):
        budget_allows(value, 750, 3600, 2)


def test_partition_excludes_completed_and_unknown_sizes() -> None:
    tasks = [
        {"task_id": f"species/SRR{i}", "accession": f"SRR{i}", "fastq_bytes": size}
        for i, size in enumerate([0, 10, 20, 30, 100], start=1)
    ]
    result = choose_partition(tasks, {"species/SRR2"}, max_bytes=50)
    assert [t["accession"] for t in result] == ["SRR3", "SRR4"]
    assert sum(t["fastq_bytes"] for t in result) <= 50


def test_large_single_task_is_isolated() -> None:
    tasks = [{"task_id": "species/SRR1", "accession": "SRR1", "fastq_bytes": 100}]
    assert choose_partition(tasks, set(), max_bytes=50) == tasks


def test_startup_requires_all_bindings_and_shell_quotes(tmp_path: Path) -> None:
    template = tmp_path / "startup.sh"
    template.write_text("KEY=@@KEY@@\nVALUE=@@VALUE@@\n")
    with pytest.raises(ValueError):
        _render_startup(template, {"KEY": "safe"})
    assert (
        _render_startup(template, {"KEY": "a; echo bad", "VALUE": 1})
        == "KEY='a; echo bad'\nVALUE=1\n"
    )


def test_completion_certificate_refuses_empty_inventory(tmp_path: Path) -> None:
    from metainformant.rna.engine.durable_quant import DirectoryStore

    with pytest.raises(ValueError, match="empty or incomplete"):
        verify_locked_campaign({"species": [], "task_count": 0}, DirectoryStore(tmp_path / "store"), "cohort", tmp_path / "result")
