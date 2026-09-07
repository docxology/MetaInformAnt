"""
test_library.py — Direct unit tests for the src/ library modules.

The thin-orchestrator contract requires all pipeline logic to live in
``src/template_bioinformatics_project/``; these tests import that library
directly and verify the computational behavior (Real-Implementation policy).
"""

import numpy as np
import pandas as pd

from template_bioinformatics_project import analysis, synthetic


# ── Stage 2: analysis functions ────────────────────────────────────────────────

class TestSummaryStatistics:
    def test_includes_cv_column(self) -> None:
        df = pd.DataFrame({"a": [1.0, 2.0, 3.0, 4.0], "b": [10.0, 10.0, 10.0, 10.0]})
        stats = analysis.compute_summary_statistics(df)
        assert "cv" in stats.columns
        assert list(stats.index) == ["a", "b"]
        # cv = std / |mean|; constant column has cv == 0
        assert stats.loc["a", "cv"] > 0
        assert stats.loc["b", "cv"] == 0

    def test_ignores_non_numeric_columns(self) -> None:
        df = pd.DataFrame({"a": [1.0, 2.0], "group": ["x", "y"]})
        stats = analysis.compute_summary_statistics(df)
        assert list(stats.index) == ["a"]


class TestPcaSummary:
    def test_variance_ratios_sum_to_at_most_one(self) -> None:
        rng = np.random.default_rng(7)
        df = pd.DataFrame(rng.normal(size=(50, 3)), columns=["a", "b", "c"])
        result = analysis.compute_pca_summary(df, n_components=2)
        assert result["n_components"] == 2
        assert 0 < result["total_variance_explained"] <= 1.0
        assert len(result["explained_variance_ratio"]) == 2

    def test_insufficient_columns_returns_error(self) -> None:
        df = pd.DataFrame({"a": [1.0, 2.0, 3.0]})
        result = analysis.compute_pca_summary(df, n_components=2)
        assert "error" in result


# ── Stage 99: synthetic data functions ────────────────────────────────────────

class TestSyntheticGeneration:
    def test_generation_is_seed_deterministic(self, tmp_path) -> None:
        """Same seed produces byte-identical CSVs (reproducibility contract)."""
        digests = []
        for out_dir in (tmp_path / "run1", tmp_path / "run2"):
            rng = np.random.default_rng(1234)
            rows = synthetic.generate_sample_dataframe(30, 4, rng)
            csv_path = out_dir / "samples_A.csv"
            synthetic.write_csv(rows, csv_path)
            digests.append(csv_path.read_bytes())
        assert digests[0] == digests[1]

    def test_write_csv_maps_nan_to_empty_field(self, tmp_path) -> None:
        rows = [{"a": 1.0, "b": float("nan")}]
        out = tmp_path / "out.csv"
        synthetic.write_csv(rows, out)
        assert out.read_text() == "a,b\n1.0,\n"
