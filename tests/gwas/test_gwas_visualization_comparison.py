"""Tests for GWAS comparison visualization functions."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from metainformant.gwas.visualization.statistical.comparison import (
    analyze_genetic_architecture,
    compare_analysis_methods,
    compare_gwas_studies,
    concordance_plot,
    cross_cohort_forest,
    miami_plot,
    multi_trait_manhattan,
    plot_qq_overlay,
)


def test_miami_plot(tmp_path: Path) -> None:
    """Test Miami plot (back-to-back Manhattan plot) for two traits."""
    # Create sample data for two traits
    trait1_results = [
        {"CHROM": "1", "POS": 1000, "p_value": 1e-8},
        {"CHROM": "1", "POS": 2000, "p_value": 1e-6},
        {"CHROM": "2", "POS": 1500, "p_value": 1e-5},
    ]

    trait2_results = [
        {"CHROM": "1", "POS": 1200, "p_value": 1e-7},
        {"CHROM": "1", "POS": 2200, "p_value": 1e-4},
        {"CHROM": "2", "POS": 1600, "p_value": 1e-6},
    ]

    output_path = tmp_path / "miami_plot.png"
    result = miami_plot(trait1_results, trait2_results, output_path, trait1_name="Trait A", trait2_name="Trait B")

    assert result["status"] == "success"
    assert output_path.exists()
    assert result["trait1_variants"] == 3
    assert result["trait2_variants"] == 3


def test_multi_trait_manhattan(tmp_path: Path) -> None:
    """Test multi-trait Manhattan plot for multiple phenotypes."""
    # Create sample data for multiple traits
    trait_data = {
        "Trait A": [
            {"CHROM": "1", "POS": 1000, "p_value": 1e-8},
            {"CHROM": "1", "POS": 2000, "p_value": 1e-6},
        ],
        "Trait B": [
            {"CHROM": "1", "POS": 1200, "p_value": 1e-7},
            {"CHROM": "1", "POS": 2200, "p_value": 1e-4},
        ],
        "Trait C": [
            {"CHROM": "1", "POS": 1100, "p_value": 1e-5},
            {"CHROM": "1", "POS": 2100, "p_value": 1e-3},
        ],
    }

    output_path = tmp_path / "multi_trait_manhattan.png"
    result = multi_trait_manhattan(trait_data, output_path, title="Multi-Trait GWAS")

    assert result["status"] == "success"
    assert output_path.exists()
    assert result["num_traits"] == 3


def test_concordance_plot(tmp_path: Path) -> None:
    """Test concordance plot between discovery and replication cohorts."""
    # Create sample data for discovery cohort
    discovery_results = [
        {"CHROM": "1", "POS": 1000, "beta": 0.8, "p_value": 1e-8},
        {"CHROM": "1", "POS": 2000, "beta": -0.6, "p_value": 1e-6},
        {"CHROM": "2", "POS": 1500, "beta": 0.4, "p_value": 1e-5},
    ]

    # Create sample data for replication cohort (matching variants)
    replication_results = [
        {"CHROM": "1", "POS": 1000, "beta": 0.7, "p_value": 1e-7},
        {"CHROM": "1", "POS": 2000, "beta": -0.5, "p_value": 1e-5},
        {"CHROM": "2", "POS": 1500, "beta": 0.3, "p_value": 1e-4},
    ]

    output_path = tmp_path / "concordance_plot.png"
    result = concordance_plot(
        discovery_results, replication_results, output_path, title="Discovery vs Replication Concordance"
    )

    assert result["status"] == "success"
    assert output_path.exists()


def test_cross_cohort_forest(tmp_path: Path) -> None:
    """Test cross-cohort forest plot (meta-analysis placeholder)."""
    # Create sample cohort data
    cohort_results = {
        "Cohort A": [
            {"CHROM": "1", "POS": 1000, "beta": 0.8, "p_value": 1e-8},
        ],
        "Cohort B": [
            {"CHROM": "1", "POS": 1000, "beta": 0.7, "p_value": 1e-7},
        ],
    }

    output_path = tmp_path / "forest_plot.png"
    result = cross_cohort_forest(cohort_results, output_path, variant_id="1:1000")

    # Currently returns skipped status as meta-analysis is not fully implemented
    assert result["status"] == "skipped"
    assert "message" in result


def test_miami_plot_empty(tmp_path: Path) -> None:
    """Test Miami plot with empty data."""
    empty_results = []
    sample_results = [{"CHROM": "1", "POS": 1000, "p_value": 1e-8}]

    output_path = tmp_path / "miami_empty.png"
    result = miami_plot(empty_results, sample_results, output_path)

    # Should handle gracefully
    assert result["status"] == "failed"
    assert "error" in result


def test_multi_trait_manhattan_single_trait(tmp_path: Path) -> None:
    """Test multi-trait Manhattan with single trait."""
    trait_data = {
        "Single Trait": [
            {"CHROM": "1", "POS": 1000, "p_value": 1e-8},
            {"CHROM": "1", "POS": 2000, "p_value": 1e-6},
        ]
    }

    output_path = tmp_path / "single_trait.png"
    result = multi_trait_manhattan(trait_data, output_path)

    # Should fail with < 2 traits
    assert result["status"] == "failed"
    assert "error" in result


def test_concordance_plot_perfect_correlation(tmp_path: Path) -> None:
    """Test concordance plot with perfect correlation."""
    # Create perfectly correlated data
    discovery_results = [
        {"CHROM": "1", "POS": 1000, "beta": 0.5, "p_value": 1e-8},
        {"CHROM": "1", "POS": 2000, "beta": -0.3, "p_value": 1e-6},
    ]

    replication_results = [
        {"CHROM": "1", "POS": 1000, "beta": 0.5, "p_value": 1e-7},  # Perfect match
        {"CHROM": "1", "POS": 2000, "beta": -0.3, "p_value": 1e-5},  # Perfect match
    ]

    output_path = tmp_path / "perfect_concordance.png"
    result = concordance_plot(discovery_results, replication_results, output_path)

    assert result["status"] == "success"
    assert output_path.exists()
    assert result["correlation"] > 0.99  # Should be nearly perfect correlation


def test_analyze_genetic_architecture_uses_p_column() -> None:
    """Architecture metrics derive from the P column; non-DataFrame or P-less inputs are skipped."""
    pd = pytest.importorskip("pandas")

    df = pd.DataFrame({"P": [0.5] * 100 + [1e-9]})
    architecture = analyze_genetic_architecture(
        {
            "PopA": {"results": df},
            "PopB": {"results": {"not": "a dataframe"}},
        }
    )

    assert "PopA" in architecture
    assert "PopB" not in architecture
    assert architecture["PopA"]["lambda_gc"] == pytest.approx(1.0, abs=0.01)
    assert architecture["PopA"]["heritability"] == pytest.approx(0.5, abs=0.01)
    assert architecture["PopA"]["polygenicity"] == pytest.approx(1 / 101)


def test_analyze_genetic_architecture_skips_missing_p_column() -> None:
    """A DataFrame without a P column contributes no entry."""
    pd = pytest.importorskip("pandas")

    architecture = analyze_genetic_architecture({"PopC": {"results": pd.DataFrame({"BETA": [1.0]})}})

    assert architecture == {}


def test_plot_qq_overlay_plots_observed_against_expected_with_zero_clipping() -> None:
    """Observed -log10(p) uses clip(min=1e-300); NaNs are dropped and the rest sorted.

    ``ndarray.clip`` has no ``lower`` keyword (that is the pandas spelling), so the
    previous implementation raised ``TypeError`` for every NumPy input.
    """
    plt = pytest.importorskip("matplotlib.pyplot")

    fig, ax = plt.subplots()
    try:
        plot_qq_overlay(ax, np.array([0.5, 0.0, 0.2, np.nan, 1e-12]), "Study", "tab:blue")

        (collection,) = ax.collections
        offsets = np.asarray(collection.get_offsets())
    finally:
        plt.close(fig)

    sorted_p = [0.0, 1e-12, 0.2, 0.5]
    expected_x = [-math.log10(rank / (len(sorted_p) + 1)) for rank in range(1, len(sorted_p) + 1)]
    expected_y = [-math.log10(max(p, 1e-300)) for p in sorted_p]
    assert offsets[:, 0] == pytest.approx(expected_x)
    assert offsets[:, 1] == pytest.approx(expected_y)
    assert offsets[0, 1] == pytest.approx(300.0)


def _study_frame(seed: int, n: int = 60) -> Any:
    import pandas as pd

    rng = np.random.default_rng(seed)
    return pd.DataFrame(
        {
            "CHR": [1] * (n // 2) + [2] * (n - n // 2),
            "P": rng.uniform(1e-6, 1.0, n),
            "BETA": rng.normal(0.0, 0.2, n),
        }
    )


def test_compare_gwas_studies_builds_figure_for_real_dataframes(tmp_path: Path) -> None:
    """Studies with a P column must produce one Q-Q overlay per study instead of crashing."""
    pytest.importorskip("pandas")
    pytest.importorskip("seaborn")
    plt = pytest.importorskip("matplotlib.pyplot")

    output_path = tmp_path / "study_comparison.png"
    fig = compare_gwas_studies(
        {"StudyA": {"results": _study_frame(1)}, "StudyB": {"results": _study_frame(2)}},
        output_file=output_path,
    )
    try:
        assert fig is not None
        assert output_path.exists() and output_path.stat().st_size > 0
        qq_axis = fig.axes[1]
        assert len(qq_axis.collections) == 2
        n_points = sorted(len(c.get_offsets()) for c in qq_axis.collections)
        assert n_points == [60, 60]
    finally:
        plt.close(fig)


def test_compare_analysis_methods_builds_figure_for_real_dataframes(tmp_path: Path) -> None:
    """Method comparison draws a Q-Q overlay per method (it shares the same overlay helper)."""
    pytest.importorskip("pandas")
    pytest.importorskip("seaborn")
    plt = pytest.importorskip("matplotlib.pyplot")

    output_path = tmp_path / "method_comparison.png"
    fig = compare_analysis_methods(
        {"Linear": {"results": _study_frame(3)}, "Mixed": {"results": _study_frame(4)}},
        output_file=output_path,
    )
    try:
        assert fig is not None
        assert output_path.exists() and output_path.stat().st_size > 0
        qq_axis = fig.axes[1]
        assert len(qq_axis.collections) == 2
    finally:
        plt.close(fig)
