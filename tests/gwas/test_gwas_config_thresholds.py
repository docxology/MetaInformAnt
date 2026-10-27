"""YAML-loaded GWAS QC thresholds are numeric whatever notation the user wrote."""

from __future__ import annotations

import pytest
import yaml

from metainformant.gwas.workflow.workflow_config import _normalize_config


def test_exponent_thresholds_without_a_decimal_point_become_floats() -> None:
    # PyYAML reads ``1e-6`` and ``5e-2`` as strings; ``0.9`` and ``20`` as numbers.
    config = yaml.safe_load("qc:\n  hwe_pval: 1e-6\n  min_maf: 5e-2\n  min_call_rate: 0.9\n  min_qual: 20\n")
    assert isinstance(config["qc"]["hwe_pval"], str)

    qc = _normalize_config(config)["quality_control"]

    assert qc == {"min_hwe_p": 1e-6, "min_maf": 0.05, "min_call_rate": 0.9, "min_qual": 20.0}
    assert all(isinstance(value, float) for value in qc.values())


def test_directly_supplied_quality_control_is_coerced_and_flags_kept() -> None:
    config = yaml.safe_load("quality_control:\n  min_hwe_p: 1e-8\n  exclude_indels: true\n")

    qc = _normalize_config(config)["quality_control"]

    assert qc == {"min_hwe_p": 1e-8, "exclude_indels": True}


def test_non_numeric_threshold_is_rejected_with_its_key() -> None:
    with pytest.raises(ValueError, match="quality_control.min_maf must be numeric"):
        _normalize_config({"qc": {"min_maf": "five percent"}})
