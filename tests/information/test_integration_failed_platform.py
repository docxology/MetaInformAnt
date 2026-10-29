"""A platform that cannot be analysed has no entropy; ``multiomics_integration`` does not report 0.0 for it.

``multiomics_integration(method="platform_entropy")`` wrapped the analysis of each sequence platform in
``except Exception`` and stored ``<platform>_entropy = 0.0``. A platform whose data are unusable (here a non-string
entry among the sequences) then looked exactly like a genuinely constant platform, whose entropy really is 0.0. The
failure is now ``None`` with the reason under ``errors``. No test doubles.
"""

from __future__ import annotations

import pytest

from metainformant.information.integration.integration import multiomics_integration

NUMERIC = [[1.0, 2.0], [2.0, 3.0], [3.0, 5.0], [4.0, 4.0]]


def test_an_unusable_sequence_platform_is_none_with_its_error_not_zero_entropy() -> None:
    result = multiomics_integration(omics_data={"dna": ["ACGT", 5], "rna": NUMERIC}, method="platform_entropy")

    assert result["dna_entropy"] is None
    assert "dna" in result["errors"] and result["errors"]["dna"]
    assert result["rna_entropy"] > 0.0  # the healthy platform is unaffected
    assert "rna" not in result["errors"]


def test_a_constant_platform_really_has_zero_entropy_and_no_error() -> None:
    result = multiomics_integration(omics_data={"dna": ["AAAA", "AAAA"]}, method="platform_entropy")

    assert result["dna_entropy"] == pytest.approx(0.0, abs=1e-12)
    assert "errors" not in result


def test_a_healthy_run_has_no_errors_key() -> None:
    result = multiomics_integration(omics_data={"dna": ["ACGT", "ACGA"], "rna": NUMERIC}, method="platform_entropy")

    assert result["dna_entropy"] > 0.0 and "errors" not in result
