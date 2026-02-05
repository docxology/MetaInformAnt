"""Tests for GWAS fine-mapping visualization functions."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from metainformant.gwas.visualization.visualization_finemapping import (
    compute_credible_set,
    conditional_analysis_plot,
    credible_set_plot,
    pip_vs_ld_plot,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_assoc_results(n: int = 20, seed: int = 42) -> list[dict]:
    """Generate synthetic association results with one strong signal."""
    rng = np.random.default_rng(seed)
    results = []
    for i in range(n):
        # First variant is the strong signal
        if i == 0:
            p = 1e-10
        else:
            p = float(rng.uniform(0.001, 1.0))
        results.append({"p_value": p, "position": 1000 + i * 100})
    return results


def _make_ld_matrix(n: int = 20, seed: int = 42) -> np.ndarray:
    """Generate a synthetic LD matrix (symmetric, diag=1, values in [0,1])."""
    rng = np.random.default_rng(seed)
    raw = rng.uniform(0, 1, size=(n, n))
    ld = (raw + raw.T) / 2.0
    np.fill_diagonal(ld, 1.0)
    return ld


# ---------------------------------------------------------------------------
# compute_credible_set
# ---------------------------------------------------------------------------


class TestComputeCredibleSet:
    """Tests for compute_credible_set."""

    def test_basic_20_variants(self) -> None:
        """PIPs sum to ~1.0 and credible set contains the lead variant."""
        results = _make_assoc_results(n=20)
        cs = compute_credible_set(results)

        assert cs["status"] == "success"
        pips = cs["pips"]
        assert len(pips) == 20
        assert math.isclose(sum(pips), 1.0, rel_tol=1e-6)

        # Lead variant (index 0, p=1e-10) must be in the credible set
        assert 0 in cs["credible_set_indices"]
        assert cs["credible_set_size"] >= 1
        assert cs["cumulative_probability"] >= 0.95

    def test_credible_level_095(self) -> None:
        """Credible set at 0.95 is smaller or equal to set at 0.99."""
        results = _make_assoc_results(n=20)
        cs_95 = compute_credible_set(results, credible_level=0.95)
        cs_99 = compute_credible_set(results, credible_level=0.99)

        assert cs_95["status"] == "success"
        assert cs_99["status"] == "success"
        assert cs_95["credible_set_size"] <= cs_99["credible_set_size"]
        assert cs_95["cumulative_probability"] >= 0.95
        assert cs_99["cumulative_probability"] >= 0.99

    def test_single_variant(self) -> None:
        """Single variant should have PIP=1.0 and be the entire credible set."""
        results = [{"p_value": 0.001, "position": 5000}]
        cs = compute_credible_set(results)

        assert cs["status"] == "success"
        assert len(cs["pips"]) == 1
        assert math.isclose(cs["pips"][0], 1.0, rel_tol=1e-6)
        assert cs["credible_set_indices"] == [0]
        assert cs["credible_set_size"] == 1

    def test_all_equal_pvalues(self) -> None:
        """Equal p-values should produce equal PIPs."""
        n = 10
        results = [{"p_value": 0.05, "position": i * 100} for i in range(n)]
        cs = compute_credible_set(results)

        assert cs["status"] == "success"
        pips = cs["pips"]
        expected_pip = 1.0 / n
        for pip in pips:
            assert math.isclose(pip, expected_pip, rel_tol=1e-6)

    def test_empty_input(self) -> None:
        """Empty input should return failed status."""
        cs = compute_credible_set([])
        assert cs["status"] == "failed"


# ---------------------------------------------------------------------------
# credible_set_plot
# ---------------------------------------------------------------------------


class TestCredibleSetPlot:
    """Tests for credible_set_plot."""

    def test_basic_plot(self, tmp_path: Path) -> None:
        """Basic credible set plot should be created."""
        results = _make_assoc_results(n=20)
        output = tmp_path / "credible_set.png"
        result = credible_set_plot(results, output_file=output)

        assert result["status"] == "success"
        assert output.exists()
        assert result["credible_set_size"] >= 1
        assert result["output_path"] == str(output)

    def test_plot_with_ld_matrix(self, tmp_path: Path) -> None:
        """Credible set plot with LD matrix coloring."""
        n = 20
        results = _make_assoc_results(n=n)
        ld = _make_ld_matrix(n=n)
        output = tmp_path / "credible_set_ld.png"

        result = credible_set_plot(results, output_file=output, ld_matrix=ld)

        assert result["status"] == "success"
        assert output.exists()
        assert result["credible_set_size"] >= 1

    def test_plot_no_output(self) -> None:
        """Plot without output file should succeed with None output_path."""
        results = _make_assoc_results(n=10)
        result = credible_set_plot(results)

        assert result["status"] == "success"
        assert result["output_path"] is None

    def test_plot_empty_input(self) -> None:
        """Empty input should return failed status."""
        result = credible_set_plot([])
        assert result["status"] == "failed"


# ---------------------------------------------------------------------------
# conditional_analysis_plot
# ---------------------------------------------------------------------------


class TestConditionalAnalysisPlot:
    """Tests for conditional_analysis_plot."""

    def test_two_rounds(self, tmp_path: Path) -> None:
        """Conditional analysis plot with 2 rounds."""
        rng = np.random.default_rng(99)
        round1 = [{"p_value": float(rng.uniform(1e-8, 1.0)), "position": i * 100} for i in range(30)]
        round2 = [{"p_value": float(rng.uniform(0.01, 1.0)), "position": i * 100} for i in range(30)]

        output = tmp_path / "conditional_2rounds.png"
        result = conditional_analysis_plot([round1, round2], output_file=output)

        assert result["status"] == "success"
        assert output.exists()
        assert result["n_rounds"] == 2

    def test_three_rounds_custom_labels(self, tmp_path: Path) -> None:
        """Conditional analysis plot with 3 rounds and custom labels."""
        rng = np.random.default_rng(77)
        rounds = []
        for _ in range(3):
            r = [{"p_value": float(rng.uniform(1e-6, 1.0)), "position": i * 50} for i in range(20)]
            rounds.append(r)

        labels = ["Unconditional", "Cond. on Lead", "Cond. on Top 2"]
        output = tmp_path / "conditional_3rounds.png"
        result = conditional_analysis_plot(rounds, output_file=output, labels=labels, title="Custom Conditional")

        assert result["status"] == "success"
        assert output.exists()
        assert result["n_rounds"] == 3

    def test_empty_rounds(self) -> None:
        """Empty rounds list should return failed status."""
        result = conditional_analysis_plot([])
        assert result["status"] == "failed"

    def test_no_output_file(self) -> None:
        """Plot without output file should succeed."""
        round1 = [{"p_value": 0.05, "position": i} for i in range(5)]
        result = conditional_analysis_plot([round1])

        assert result["status"] == "success"
        assert result["output_path"] is None
        assert result["n_rounds"] == 1


# ---------------------------------------------------------------------------
# pip_vs_ld_plot
# ---------------------------------------------------------------------------


class TestPipVsLdPlot:
    """Tests for pip_vs_ld_plot."""

    def test_basic_plot(self, tmp_path: Path) -> None:
        """Basic PIP vs LD plot."""
        rng = np.random.default_rng(123)
        n = 30
        pips = rng.uniform(0, 1, n).tolist()
        ld_vals = rng.uniform(0, 1, n).tolist()

        output = tmp_path / "pip_vs_ld.png"
        result = pip_vs_ld_plot(pips, ld_vals, output_file=output)

        assert result["status"] == "success"
        assert output.exists()
        assert result["n_variants"] == n

    def test_with_variant_ids(self, tmp_path: Path) -> None:
        """PIP vs LD plot with variant ID annotations."""
        n = 15
        pips = [0.9, 0.7, 0.5] + [0.05] * (n - 3)
        ld_vals = [0.95, 0.8, 0.6] + [float(i) / n for i in range(n - 3)]
        ids = [f"rs{1000 + i}" for i in range(n)]

        output = tmp_path / "pip_vs_ld_annotated.png"
        result = pip_vs_ld_plot(pips, ld_vals, output_file=output, variant_ids=ids)

        assert result["status"] == "success"
        assert output.exists()
        assert result["n_variants"] == n

    def test_empty_input(self) -> None:
        """Empty inputs should return failed status."""
        result = pip_vs_ld_plot([], [])
        assert result["status"] == "failed"

    def test_length_mismatch(self) -> None:
        """Mismatched lengths should return failed status."""
        result = pip_vs_ld_plot([0.5, 0.3], [0.1])
        assert result["status"] == "failed"

    def test_no_output_file(self) -> None:
        """Plot without output file should succeed."""
        pips = [0.8, 0.1, 0.05, 0.05]
        ld_vals = [1.0, 0.9, 0.3, 0.1]
        result = pip_vs_ld_plot(pips, ld_vals)

        assert result["status"] == "success"
        assert result["output_path"] is None
        assert result["n_variants"] == 4
