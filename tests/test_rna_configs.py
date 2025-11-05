"""Tests for metainformant.rna.configs module.

Tests species profile configuration, run layout generation, and step parameter building.
"""

from __future__ import annotations

from pathlib import Path

from metainformant.rna.configs import AmalgkitRunLayout, SpeciesProfile, build_step_params
from metainformant.rna.workflow import AmalgkitWorkflowConfig
from metainformant.rna.workflow import plan_workflow as plan_basic
from metainformant.rna.workflow import plan_workflow_with_params


def test_build_step_params_includes_species_tissues_and_layout(tmp_path: Path):
    """Test that build_step_params includes species, tissues, and layout information in step parameters."""
    species = SpeciesProfile(
        name="Apis mellifera",
        taxon_id=7460,
        tissues=["brain", "thorax"],
    )
    layout = AmalgkitRunLayout(base_dir=tmp_path)

    params_map = build_step_params(species, layout)

    # metadata/select should include taxon and tissues
    md = params_map["metadata"]
    sel = params_map["select"]
    assert md.get("taxon-id") == 7460
    assert sel.get("taxon-id") == 7460
    assert md.get("tissue") == ["brain", "thorax"]
    assert sel.get("tissue") == ["brain", "thorax"]

    # directories for other steps
    assert params_map["getfastq"]["out-dir"] == layout.fastq_dir
    assert params_map["quant"]["out-dir"] == layout.quant_dir
    assert params_map["merge"]["out"] == layout.merge_table
    assert params_map["cstmm"]["out-dir"] == layout.cstmm_dir
    assert params_map["curate"]["out-dir"] == layout.curate_dir
    assert params_map["csca"]["out-dir"] == layout.csca_dir


def test_plan_workflow_with_params_merges_common_and_specific(tmp_path: Path):
    """Test that plan_workflow_with_params correctly merges common config with step-specific parameters."""
    cfg = AmalgkitWorkflowConfig(work_dir=tmp_path, threads=7)
    species = SpeciesProfile(name="Apis mellifera", taxon_id=7460, tissues=["brain"])
    layout = AmalgkitRunLayout(base_dir=tmp_path)
    params_map = build_step_params(species, layout)

    steps = plan_workflow_with_params(cfg, params_map)

    # common threads + specific params should appear
    for name, params in steps:
        assert params.get("threads") == 7
        if name in ("metadata", "select"):
            assert params.get("taxon-id") == 7460
            assert params.get("tissue") == ["brain"]
