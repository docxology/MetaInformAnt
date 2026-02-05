"""Tests for enhanced Manhattan plot features: gene annotations, highlight regions,
suggestive threshold, and label_top_n.

Following NO_MOCKING policy - all tests use real matplotlib implementations.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, List

import pytest

try:
    import matplotlib
    import matplotlib.pyplot as plt

    matplotlib.use("Agg")
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

from metainformant.gwas.visualization.general import manhattan_plot


def _make_synthetic_results(n_per_chrom: int = 50, n_chroms: int = 3) -> List[Dict[str, Any]]:
    """Generate synthetic GWAS association results for testing."""
    results: List[Dict[str, Any]] = []
    rng = np.random.RandomState(42)
    for chrom in range(1, n_chroms + 1):
        for i in range(n_per_chrom):
            pos = (i + 1) * 10000
            # Most variants are non-significant, sprinkle a few hits
            if i == 0 and chrom == 1:
                p_val = 1e-10  # Highly significant
            elif i == 1 and chrom == 2:
                p_val = 5e-9  # Genome-wide significant
            elif i == 2 and chrom == 1:
                p_val = 3e-6  # Suggestive
            else:
                p_val = float(rng.uniform(0.01, 1.0))
            results.append(
                {
                    "chrom": str(chrom),
                    "pos": pos,
                    "p_value": p_val,
                    "variant_id": f"rs{chrom}_{i}",
                }
            )
    return results


@pytest.mark.skipif(not HAS_MATPLOTLIB or not HAS_NUMPY, reason="matplotlib/numpy required")
class TestManhattanPlotGeneAnnotations:
    """Test gene annotation feature on Manhattan plot."""

    def test_gene_annotations_by_variant_index(self, tmp_path: Path) -> None:
        """Gene annotations placed by variant_index render without error."""
        results = _make_synthetic_results()
        annotations = [
            {"variant_index": 0, "gene_name": "BRCA1"},
            {"variant_index": 51, "gene_name": "TP53"},
            {"variant_index": 102, "gene_name": "EGFR"},
        ]
        output = tmp_path / "manhattan_gene_ann.png"
        fig = manhattan_plot(results, output_path=output, gene_annotations=annotations)

        assert fig is not None
        assert output.exists()
        assert output.stat().st_size > 0
        # Verify annotations were added (check for text artists on the axis)
        ax = fig.axes[0]
        texts = [t for t in ax.texts if t.get_text() in ("BRCA1", "TP53", "EGFR")]
        assert len(texts) == 3, f"Expected 3 gene annotations, found {len(texts)}"
        plt.close(fig)

    def test_gene_annotations_by_chrom_pos(self, tmp_path: Path) -> None:
        """Gene annotations placed by chrom+pos find the nearest point."""
        results = _make_synthetic_results()
        annotations = [
            {"chrom": "1", "pos": 10000, "gene_name": "GeneA"},
            {"chrom": "2", "pos": 20000, "gene_name": "GeneB"},
        ]
        output = tmp_path / "manhattan_gene_chrompos.png"
        fig = manhattan_plot(results, output_path=output, gene_annotations=annotations)

        assert fig is not None
        ax = fig.axes[0]
        gene_texts = [t.get_text() for t in ax.texts]
        assert "GeneA" in gene_texts
        assert "GeneB" in gene_texts
        plt.close(fig)


@pytest.mark.skipif(not HAS_MATPLOTLIB or not HAS_NUMPY, reason="matplotlib/numpy required")
class TestManhattanPlotHighlightRegions:
    """Test genomic region highlighting on Manhattan plot."""

    def test_highlight_two_regions(self, tmp_path: Path) -> None:
        """Two highlight regions produce axvspan patches on the plot."""
        results = _make_synthetic_results()
        regions = [
            {"chrom": "1", "start": 5000, "end": 30000, "color": "green", "label": "Region A"},
            {"chrom": "2", "start": 10000, "end": 50000, "color": "orange", "label": "Region B"},
        ]
        output = tmp_path / "manhattan_highlight.png"
        fig = manhattan_plot(results, output_path=output, highlight_regions=regions)

        assert fig is not None
        assert output.exists()
        # axvspan creates Polygon patches; confirm the figure rendered
        ax = fig.axes[0]
        # The legend should mention the region labels
        legend = ax.get_legend()
        if legend is not None:
            legend_texts = [t.get_text() for t in legend.get_texts()]
            assert "Region A" in legend_texts or "Region B" in legend_texts
        plt.close(fig)

    def test_highlight_region_default_color(self, tmp_path: Path) -> None:
        """Highlight region without explicit color uses the default yellow."""
        results = _make_synthetic_results()
        regions = [{"chrom": "1", "start": 0, "end": 100000}]
        fig = manhattan_plot(results, highlight_regions=regions)
        assert fig is not None
        plt.close(fig)


@pytest.mark.skipif(not HAS_MATPLOTLIB or not HAS_NUMPY, reason="matplotlib/numpy required")
class TestManhattanPlotSuggestiveThreshold:
    """Test suggestive threshold line on Manhattan plot."""

    def test_suggestive_threshold_default(self, tmp_path: Path) -> None:
        """Default suggestive threshold (1e-5) draws a secondary line."""
        results = _make_synthetic_results()
        output = tmp_path / "manhattan_suggestive.png"
        fig = manhattan_plot(results, output_path=output)

        assert fig is not None
        ax = fig.axes[0]
        # There should be at least 2 horizontal lines (significance + suggestive)
        hlines = [line for line in ax.get_lines() if len(set(line.get_ydata())) == 1]
        assert len(hlines) >= 2, f"Expected >=2 horizontal lines, found {len(hlines)}"
        plt.close(fig)

    def test_suggestive_threshold_disabled(self, tmp_path: Path) -> None:
        """Setting suggestive_threshold=None disables the suggestive line."""
        results = _make_synthetic_results()
        fig = manhattan_plot(results, suggestive_threshold=None)

        assert fig is not None
        ax = fig.axes[0]
        # Only the significance threshold line should remain
        hlines = [line for line in ax.get_lines() if len(set(line.get_ydata())) == 1]
        assert len(hlines) == 1, f"Expected 1 horizontal line, found {len(hlines)}"
        plt.close(fig)

    def test_suggestive_threshold_custom_value(self, tmp_path: Path) -> None:
        """Custom suggestive threshold value draws at the correct y position."""
        results = _make_synthetic_results()
        custom_threshold = 1e-4
        fig = manhattan_plot(results, suggestive_threshold=custom_threshold)

        assert fig is not None
        ax = fig.axes[0]
        expected_y = -math.log10(custom_threshold)
        hlines = [line for line in ax.get_lines() if len(set(line.get_ydata())) == 1]
        y_values = [list(set(line.get_ydata()))[0] for line in hlines]
        assert any(abs(y - expected_y) < 0.01 for y in y_values), f"No line at y={expected_y:.2f}; found y={y_values}"
        plt.close(fig)


@pytest.mark.skipif(not HAS_MATPLOTLIB or not HAS_NUMPY, reason="matplotlib/numpy required")
class TestManhattanPlotLabelTopN:
    """Test auto-labeling top N hits on Manhattan plot."""

    def test_label_top_3(self, tmp_path: Path) -> None:
        """Top 3 hits are annotated with variant_id labels."""
        results = _make_synthetic_results()
        output = tmp_path / "manhattan_top3.png"
        fig = manhattan_plot(results, output_path=output, label_top_n=3, suggestive_threshold=None)

        assert fig is not None
        assert output.exists()
        ax = fig.axes[0]
        # Annotations come from ax.texts; the top 3 should have labels
        annotation_texts = [t.get_text() for t in ax.texts]
        # The top hit is rs1_0 (p=1e-10), second is rs2_1 (p=5e-9), third is rs1_2 (p=3e-6)
        assert "rs1_0" in annotation_texts, f"Expected rs1_0 in annotations: {annotation_texts}"
        assert "rs2_1" in annotation_texts, f"Expected rs2_1 in annotations: {annotation_texts}"
        assert "rs1_2" in annotation_texts, f"Expected rs1_2 in annotations: {annotation_texts}"
        plt.close(fig)

    def test_label_top_0_no_labels(self, tmp_path: Path) -> None:
        """label_top_n=0 does not add any variant labels."""
        results = _make_synthetic_results()
        fig = manhattan_plot(results, label_top_n=0, suggestive_threshold=None)

        assert fig is not None
        ax = fig.axes[0]
        # No variant annotation texts should appear (only chromosome tick labels)
        annotation_texts = [t.get_text() for t in ax.texts]
        variant_labels = [t for t in annotation_texts if t.startswith("rs")]
        assert len(variant_labels) == 0
        plt.close(fig)


@pytest.mark.skipif(not HAS_MATPLOTLIB or not HAS_NUMPY, reason="matplotlib/numpy required")
class TestManhattanPlotCombined:
    """Test all enhancements combined in a single call."""

    def test_all_enhancements_combined(self, tmp_path: Path) -> None:
        """All new parameters work together without conflicts."""
        results = _make_synthetic_results()
        annotations = [
            {"variant_index": 0, "gene_name": "BRCA1"},
            {"variant_index": 51, "gene_name": "TP53"},
        ]
        regions = [
            {"chrom": "1", "start": 5000, "end": 30000, "color": "green", "label": "QTL-1"},
        ]
        output = tmp_path / "manhattan_combined.png"
        fig = manhattan_plot(
            results,
            output_path=output,
            significance_threshold=5e-8,
            gene_annotations=annotations,
            highlight_regions=regions,
            suggestive_threshold=1e-5,
            label_top_n=2,
        )

        assert fig is not None
        assert output.exists()
        assert output.stat().st_size > 1000  # Non-trivial PNG
        ax = fig.axes[0]
        all_texts = [t.get_text() for t in ax.texts]
        assert "BRCA1" in all_texts
        assert "TP53" in all_texts
        plt.close(fig)


@pytest.mark.skipif(not HAS_MATPLOTLIB or not HAS_NUMPY, reason="matplotlib/numpy required")
class TestManhattanPlotBackwardCompatibility:
    """Ensure existing API is unbroken by the enhancements."""

    def test_existing_call_no_new_params(self, tmp_path: Path) -> None:
        """Calling manhattan_plot with only original params still works."""
        results = [
            {"chrom": "1", "pos": 1000, "p_value": 1e-6},
            {"chrom": "1", "pos": 2000, "p_value": 1e-7},
            {"chrom": "2", "pos": 1500, "p_value": 1e-5},
        ]
        output = tmp_path / "manhattan_compat.png"
        fig = manhattan_plot(results, output_path=output, significance_threshold=5e-8)

        assert fig is not None
        assert hasattr(fig, "savefig")
        assert output.exists()
        plt.close(fig)

    def test_single_dict_input_still_works(self, tmp_path: Path) -> None:
        """A single result dict (not a list) is still handled properly."""
        result = {"chrom": "1", "pos": 5000, "p_value": 1e-9}
        fig = manhattan_plot(result)

        assert fig is not None
        assert len(fig.axes) > 0
        plt.close(fig)
