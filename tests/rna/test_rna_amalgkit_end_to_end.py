"""Bounded current Amalgkit workflow smoke tests."""

from __future__ import annotations

from pathlib import Path

from metainformant.rna.amalgkit import amalgkit
from metainformant.rna.engine.workflow import AmalgkitWorkflowConfig, plan_workflow


def test_current_default_plan_is_ordered_and_current(tmp_path: Path) -> None:
    config = AmalgkitWorkflowConfig(
        work_dir=tmp_path / "work", threads=2, species_list=["Apis_mellifera"]
    )
    planned = plan_workflow(config)
    names = [name for name, _ in planned]
    assert names == [
        "metadata",
        "select",
        "getfastq",
        "integrate",
        "quant",
        "merge",
        "wsfilter",
        "finalize",
        "sanity",
    ]
    assert "merge" in names and names.index("merge") < names.index("wsfilter")
    assert names.index("wsfilter") < names.index("finalize")


def test_current_downstream_commands_build_without_unsupported_flags() -> None:
    common = {
        "out_dir": "work",
        "metadata": "work/wsfilter/metadata.tsv",
        "input_dir": "work/wsfilter",
        "norm": "log2p1-fpkm",
        "batch_effect_alg": "no",
        "threads": 2,
    }
    for step in ("wsfilter", "finalize", "csfilter"):
        command = amalgkit.build_amalgkit_command(step, common)
        assert command[:2] == ["amalgkit", step]
        assert "--input_dir" in command
        assert "--work_dir" not in command
