"""Tests for truthful planning/command rendering.

These tests do NOT require the external `amalgkit` tool; they validate our own
planning defaults and CLI-argument sanitization logic.
"""

from __future__ import annotations

from pathlib import Path


def test_plan_applies_defaults_and_sanitizes_internal_only_flags(tmp_path: Path) -> None:
    from metainformant.rna.workflow import AmalgkitWorkflowConfig, apply_step_defaults, plan_workflow, sanitize_params_for_cli

    cfg = AmalgkitWorkflowConfig(
        work_dir=tmp_path / "work",
        threads=4,
        species_list=["Example_species"],
        per_step={
            "getfastq": {"out_dir": str(tmp_path / "fastq"), "num_download_workers": 5, "aws": "yes"},
        },
    )

    apply_step_defaults(cfg)
    planned = plan_workflow(cfg)
    steps = dict(planned)

    # Defaults should exist (these are required for truthful `--plan`)
    assert "config" in steps
    cfg_params = sanitize_params_for_cli("config", steps["config"])
    assert "out_dir" in cfg_params

    # Internal-only keys should never appear in rendered CLI params
    gf_params = sanitize_params_for_cli("getfastq", steps["getfastq"])
    assert "num_download_workers" not in gf_params
    assert "aws" in gf_params



