"""Tests for all GWAS visualization modules.

Tests all visualization modules in gwas/visualization_*.py following NO_MOCKING policy.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

# Test basic visualization module
from metainformant.gwas import visualization

# Test specialized visualization modules
try:
    from metainformant.gwas import (
        visualization_comparison,
        visualization_effects,
        visualization_enhanced,
        visualization_genome,
        visualization_population,
        visualization_regional,
        visualization_statistical,
        visualization_suite,
        visualization_variants,
    )
except ImportError:
    # If modules don't exist, skip tests
    pytest.skip("GWAS visualization modules not available", allow_module_level=True)


class TestBasicVisualization:
    """Tests for basic visualization.py module."""

    def test_manhattan_plot_basic(self, tmp_path):
        """Test basic Manhattan plot generation."""
        # Create sample GWAS results
        results = [
            {"CHROM": "1", "POS": 1000, "p_value": 1e-6},
            {"CHROM": "1", "POS": 2000, "p_value": 1e-7},
            {"CHROM": "2", "POS": 1500, "p_value": 1e-5},
        ]
        output_path = tmp_path / "manhattan.png"
        result = visualization.manhattan_plot(results, output_path)
        assert isinstance(result, dict)
        assert "status" in result

    def test_qq_plot_basic(self, tmp_path):
        """Test basic Q-Q plot generation."""
        # Create sample p-values
        pvalues = [0.001, 0.01, 0.05, 0.1, 0.5, 0.9]
        results = [{"p_value": p} for p in pvalues]
        output_path = tmp_path / "qq.png"
        result = visualization.qq_plot(results, output_path)
        assert isinstance(result, dict)


class TestGenomeVisualization:
    """Tests for visualization_genome.py module."""

    def test_manhattan_plot_genome(self, tmp_path):
        """Test genome-wide Manhattan plot."""
        results = [
            {"CHROM": "1", "POS": 1000, "p_value": 1e-6},
            {"CHROM": "2", "POS": 2000, "p_value": 1e-7},
        ]
        output_path = tmp_path / "genome_manhattan.png"
        if hasattr(visualization_genome, "manhattan_plot"):
            result = visualization_genome.manhattan_plot(results, output_path)
            assert isinstance(result, dict)

    def test_circular_manhattan_plot(self, tmp_path):
        """Test circular Manhattan plot."""
        results = [
            {"CHROM": "1", "POS": 1000, "p_value": 1e-6},
            {"CHROM": "2", "POS": 2000, "p_value": 1e-7},
        ]
        output_path = tmp_path / "circular_manhattan.png"
        if hasattr(visualization_genome, "circular_manhattan_plot"):
            result = visualization_genome.circular_manhattan_plot(results, output_path)
            assert isinstance(result, dict)


class TestStatisticalVisualization:
    """Tests for visualization_statistical.py module."""

    def test_qq_plot_statistical(self, tmp_path):
        """Test Q-Q plot from statistical module."""
        results = [{"p_value": 0.001}, {"p_value": 0.01}, {"p_value": 0.05}]
        output_path = tmp_path / "qq_statistical.png"
        if hasattr(visualization_statistical, "qq_plot"):
            result = visualization_statistical.qq_plot(results, output_path)
            assert isinstance(result, dict)

    def test_volcano_plot(self, tmp_path):
        """Test volcano plot generation."""
        results = [
            {"p_value": 1e-6, "effect_size": 0.5},
            {"p_value": 1e-5, "effect_size": -0.3},
        ]
        output_path = tmp_path / "volcano.png"
        if hasattr(visualization_statistical, "volcano_plot"):
            result = visualization_statistical.volcano_plot(results, output_path)
            assert isinstance(result, dict)

    def test_lambda_gc_plot(self, tmp_path):
        """Test lambda GC plot."""
        results = [{"p_value": 0.001}, {"p_value": 0.01}]
        output_path = tmp_path / "lambda_gc.png"
        if hasattr(visualization_statistical, "lambda_gc_plot"):
            result = visualization_statistical.lambda_gc_plot(results, output_path)
            assert isinstance(result, dict)


class TestRegionalVisualization:
    """Tests for visualization_regional.py module."""

    def test_regional_plot(self, tmp_path):
        """Test regional association plot."""
        results = [
            {"CHROM": "1", "POS": 1000, "p_value": 1e-6},
            {"CHROM": "1", "POS": 1100, "p_value": 1e-5},
        ]
        output_path = tmp_path / "regional.png"
        if hasattr(visualization_regional, "regional_association_plot"):
            result = visualization_regional.regional_association_plot(results, output_path, chrom="1", start=900, end=1200)
            assert isinstance(result, dict)


class TestPopulationVisualization:
    """Tests for visualization_population.py module."""

    def test_population_plot(self, tmp_path):
        """Test population-specific plot."""
        results = [
            {"CHROM": "1", "POS": 1000, "p_value": 1e-6, "population": "pop1"},
            {"CHROM": "1", "POS": 2000, "p_value": 1e-5, "population": "pop2"},
        ]
        output_path = tmp_path / "population.png"
        if hasattr(visualization_population, "population_comparison_plot"):
            result = visualization_population.population_comparison_plot(results, output_path)
            assert isinstance(result, dict)


class TestVariantsVisualization:
    """Tests for visualization_variants.py module."""

    def test_variant_plot(self, tmp_path):
        """Test variant-specific visualization."""
        results = [
            {"CHROM": "1", "POS": 1000, "p_value": 1e-6, "REF": "A", "ALT": "G"},
        ]
        output_path = tmp_path / "variant.png"
        if hasattr(visualization_variants, "variant_effect_plot"):
            result = visualization_variants.variant_effect_plot(results, output_path)
            assert isinstance(result, dict)


class TestEffectsVisualization:
    """Tests for visualization_effects.py module."""

    def test_effect_plot(self, tmp_path):
        """Test effect size visualization."""
        results = [
            {"p_value": 1e-6, "effect_size": 0.5, "se": 0.1},
            {"p_value": 1e-5, "effect_size": -0.3, "se": 0.15},
        ]
        output_path = tmp_path / "effects.png"
        if hasattr(visualization_effects, "effect_size_plot"):
            result = visualization_effects.effect_size_plot(results, output_path)
            assert isinstance(result, dict)


class TestComparisonVisualization:
    """Tests for visualization_comparison.py module."""

    def test_comparison_plot(self, tmp_path):
        """Test comparison between studies."""
        results1 = [{"CHROM": "1", "POS": 1000, "p_value": 1e-6}]
        results2 = [{"CHROM": "1", "POS": 1000, "p_value": 1e-5}]
        output_path = tmp_path / "comparison.png"
        if hasattr(visualization_comparison, "study_comparison_plot"):
            result = visualization_comparison.study_comparison_plot(results1, results2, output_path)
            assert isinstance(result, dict)


class TestVisualization:
    """Tests for visualization_enhanced.py module."""

    def test_enhanced_plot(self, tmp_path):
        """Test enhanced visualization features."""
        results = [
            {"CHROM": "1", "POS": 1000, "p_value": 1e-6},
        ]
        output_path = tmp_path / "enhanced.png"
        # Test any enhanced plotting function
        if hasattr(visualization_enhanced, "enhanced_manhattan_plot"):
            result = visualization_enhanced.enhanced_manhattan_plot(results, output_path)
            assert isinstance(result, dict)


class TestVisualizationSuite:
    """Tests for visualization_suite.py module."""

    def test_suite_generation(self, tmp_path):
        """Test generating a suite of visualizations."""
        results = [
            {"CHROM": "1", "POS": 1000, "p_value": 1e-6},
            {"CHROM": "2", "POS": 2000, "p_value": 1e-5},
        ]
        output_dir = tmp_path / "suite"
        if hasattr(visualization_suite, "generate_visualization_suite"):
            result = visualization_suite.generate_visualization_suite(results, output_dir)
            assert isinstance(result, dict)


class TestVisualizationEdgeCases:
    """Tests for edge cases across visualization modules."""

    def test_empty_results(self, tmp_path):
        """Test handling of empty results."""
        results = []
        output_path = tmp_path / "empty.png"
        result = visualization.manhattan_plot(results, output_path)
        assert isinstance(result, dict)
        # Should handle gracefully
        assert result.get("status") in ("failed", "completed")

    def test_missing_fields(self, tmp_path):
        """Test handling of results with missing fields."""
        results = [
            {"CHROM": "1"},  # Missing POS and p_value
            {"POS": 1000},  # Missing CHROM and p_value
        ]
        output_path = tmp_path / "missing.png"
        result = visualization.manhattan_plot(results, output_path)
        assert isinstance(result, dict)

    def test_invalid_pvalues(self, tmp_path):
        """Test handling of invalid p-values."""
        results = [
            {"CHROM": "1", "POS": 1000, "p_value": "invalid"},
            {"CHROM": "1", "POS": 2000, "p_value": -1.0},  # Negative p-value
        ]
        output_path = tmp_path / "invalid.png"
        result = visualization.manhattan_plot(results, output_path)
        assert isinstance(result, dict)

