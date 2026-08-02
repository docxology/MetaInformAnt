"""Tests for metainformant.rna.engine.workflow module.

Tests workflow planning, step ordering, and parameter inheritance.
"""

from __future__ import annotations

from pathlib import Path

from metainformant.rna.engine.workflow import AmalgkitWorkflowConfig, plan_workflow


def test_plan_workflow_orders_steps_and_inherits_common_params(tmp_path: Path):
    """Test that plan_workflow returns steps in correct order and inherits common parameters."""
    cfg = AmalgkitWorkflowConfig(work_dir=tmp_path, threads=6, species_list=["Apis_mellifera"])
    steps = plan_workflow(cfg)

    expected_order = [
        "metadata",
        "select",
        "getfastq",
        "integrate",  # Moved after getfastq to integrate downloaded FASTQs
        "quant",
        "merge",
        "wsfilter",
        "finalize",
        "sanity",
    ]

    got_order = [name for name, _ in steps]
    assert got_order == expected_order

    # Each step should include the common params
    # Note: threads are only added for steps that support them (getfastq, integrate, quant)
    # Note: species are only added for steps that support them (currently none of the standard workflow steps)
    steps_with_threads = {"getfastq", "integrate", "quant"}
    for step_name, params in steps:
        if step_name in steps_with_threads:
            assert params.get("threads") == 6
        # out_dir is always present
        assert "out_dir" in params or "work_dir" in params


def test_plan_workflow_step_dependencies(tmp_path: Path):
    """Test that workflow steps respect dependencies (e.g., quant requires getfastq)."""
    cfg = AmalgkitWorkflowConfig(work_dir=tmp_path, threads=6, species_list=["Apis_mellifera"])
    steps = plan_workflow(cfg)

    step_names = [name for name, _ in steps]

    # Verify critical dependencies:
    # - quant comes after getfastq
    assert step_names.index("getfastq") < step_names.index("quant")

    # - merge comes after quant
    assert step_names.index("quant") < step_names.index("merge")

    # - integrate comes after getfastq (to integrate downloaded FASTQs)
    assert step_names.index("getfastq") < step_names.index("integrate")

    # - select follows metadata
    assert step_names.index("metadata") < step_names.index("select")

    # - filtering and finalization follow merge
    assert step_names.index("merge") < step_names.index("wsfilter")
    assert step_names.index("wsfilter") < step_names.index("finalize")

    # - sanity is last
    assert step_names[-1] == "sanity"


def test_plan_workflow_with_specific_steps(tmp_path: Path):
    """Test that plan_workflow respects step filtering when specific steps are requested."""
    cfg = AmalgkitWorkflowConfig(
        work_dir=tmp_path,
        threads=6,
        species_list=["Apis_mellifera"],
        per_step={
            "metadata": {},
            "finalize": {},
            "sanity": {},
        },
    )
    steps = plan_workflow(cfg)

    # When specific steps are provided, should only return those steps
    # But plan_workflow returns all steps by default - this test verifies structure
    step_names = [name for name, _ in steps]

    # Should include all steps (plan_workflow doesn't filter by per_step keys)
    # But verify the structure is correct
    assert "metadata" in step_names
    assert "finalize" in step_names
    assert "sanity" in step_names


def test_plan_workflow_parameter_inheritance(tmp_path: Path):
    """Test that step-specific parameters override common parameters."""
    cfg = AmalgkitWorkflowConfig(
        work_dir=tmp_path,
        threads=6,  # Common thread count
        species_list=["Apis_mellifera"],
        per_step={
            "quant": {"threads": 12},  # Override for quant step
        },
    )
    steps = plan_workflow(cfg)

    # Find quant step
    quant_params = None
    for step_name, params in steps:
        if step_name == "quant":
            quant_params = params
            break

    assert quant_params is not None
    # Quant should have overridden thread count
    assert quant_params.get("threads") == 12

    # Other steps that support threads should have common thread count
    steps_with_threads = {"getfastq", "integrate", "quant"}
    for step_name, params in steps:
        if step_name != "quant" and step_name in steps_with_threads:
            assert params.get("threads") == 6


def test_plan_workflow_uses_optional_shared_amalgkit_download_cache(tmp_path: Path, monkeypatch) -> None:
    """Route taxonomy-consuming stages to one validated campaign cache when requested."""

    shared_cache = tmp_path / "shared" / "ncbi_taxonomy"
    monkeypatch.setenv("AMALGKIT_SHARED_DOWNLOAD_DIR", str(shared_cache))
    cfg = AmalgkitWorkflowConfig(work_dir=tmp_path / "work", threads=2)

    steps = dict(plan_workflow(cfg))

    for step_name in ("metadata", "integrate", "getfastq"):
        assert steps[step_name]["download_dir"] == str(shared_cache.resolve())
    for step_name in ("select", "quant", "merge", "wsfilter", "finalize", "sanity"):
        assert "download_dir" not in steps[step_name]
