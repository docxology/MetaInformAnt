"""Tests for metainformant.rna.workflow module.

Tests workflow planning, step ordering, and parameter inheritance.
"""

from __future__ import annotations

from pathlib import Path

from metainformant.rna.workflow import AmalgkitWorkflowConfig, plan_workflow


def test_plan_workflow_orders_steps_and_inherits_common_params(tmp_path: Path):
    """Test that plan_workflow returns steps in correct order and inherits common parameters."""
    cfg = AmalgkitWorkflowConfig(work_dir=tmp_path, threads=6, species_list=["Apis_mellifera"])
    steps = plan_workflow(cfg)

    expected_order = [
        "metadata",
        "config",
        "select",
        "getfastq",
        "integrate",  # Moved after getfastq to integrate downloaded FASTQs
        "quant",
        "merge",
        "cstmm",
        "curate",
        "csca",
        "sanity",
    ]

    got_order = [name for name, _ in steps]
    assert got_order == expected_order

    # Each step should include the common params
    for _, params in steps:
        assert params.get("threads") == 6
        assert params.get("species-list") == ["Apis_mellifera"]
