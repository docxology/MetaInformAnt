"""Tests for interactive GWAS visualization functions (Plotly-based)."""

from __future__ import annotations

import random
from pathlib import Path
from typing import Any, Dict, List

import pytest

from metainformant.gwas.visualization.visualization_interactive import (
    interactive_manhattan,
    interactive_pca,
    interactive_volcano,
)

# ---------------------------------------------------------------------------
# Helpers for synthetic data generation
# ---------------------------------------------------------------------------


def _make_assoc_results(n: int, seed: int = 42) -> List[Dict[str, Any]]:
    """Generate *n* synthetic association results across multiple chromosomes."""
    rng = random.Random(seed)
    chroms = [str(c) for c in range(1, 23)]
    results: List[Dict[str, Any]] = []
    for i in range(n):
        chrom = rng.choice(chroms)
        pos = rng.randint(1, 250_000_000)
        pval = 10 ** rng.uniform(-10, 0)
        beta = rng.gauss(0, 0.5)
        results.append(
            {
                "variant_id": f"rs{100000 + i}",
                "chromosome": chrom,
                "position": pos,
                "p_value": pval,
                "beta": beta,
            }
        )
    return results


def _make_pca_data(n_samples: int, n_pcs: int = 5, seed: int = 42) -> Dict[str, Any]:
    """Generate synthetic PCA data with explained variance ratios."""
    rng = random.Random(seed)
    pcs = [[rng.gauss(0, 1) for _ in range(n_pcs)] for _ in range(n_samples)]
    # Create decreasing explained variance ratios that sum to < 1
    raw_var = sorted([rng.random() for _ in range(n_pcs)], reverse=True)
    total = sum(raw_var) * 1.5  # ensure they sum to less than 1
    explained = [v / total for v in raw_var]
    sample_ids = [f"SAMPLE_{i:03d}" for i in range(n_samples)]
    return {
        "pcs": pcs,
        "explained_variance_ratio": explained,
        "sample_ids": sample_ids,
    }


# ---------------------------------------------------------------------------
# interactive_manhattan tests
# ---------------------------------------------------------------------------


class TestInteractiveManhattan:
    """Tests for the interactive_manhattan function."""

    def test_basic_manhattan(self, tmp_path: Path) -> None:
        """Manhattan plot with 200 synthetic results produces valid HTML."""
        results = _make_assoc_results(200)
        output = tmp_path / "manhattan.html"
        result = interactive_manhattan(results, output_file=output)

        assert result["status"] == "success"
        assert result["n_variants"] == 200
        assert output.exists()
        assert output.stat().st_size > 0

        content = output.read_text(encoding="utf-8")
        assert "<html" in content.lower() or "plotly" in content.lower()

    def test_manhattan_with_annotations(self, tmp_path: Path) -> None:
        """Manhattan plot with gene annotations adds annotation markers."""
        results = _make_assoc_results(200)
        # Pick a few known variant IDs from our generated data for annotation
        annotations = [
            {"variant_id": results[0]["variant_id"], "gene": "BRCA1"},
            {"variant_id": results[10]["variant_id"], "gene": "TP53"},
        ]
        output = tmp_path / "manhattan_annotated.html"
        result = interactive_manhattan(results, output_file=output, gene_annotations=annotations)

        assert result["status"] == "success"
        assert result["n_variants"] == 200
        assert output.exists()
        assert output.stat().st_size > 0

    def test_manhattan_empty_results(self, tmp_path: Path) -> None:
        """Empty results return skipped status without error."""
        output = tmp_path / "manhattan_empty.html"
        result = interactive_manhattan([], output_file=output)

        assert result["status"] == "skipped"
        assert result["n_variants"] == 0

    def test_manhattan_no_output_file(self) -> None:
        """Manhattan plot without output file still returns success."""
        results = _make_assoc_results(50)
        result = interactive_manhattan(results, output_file=None)

        assert result["status"] == "success"
        assert result["n_variants"] == 50
        assert result["output_path"] is None


# ---------------------------------------------------------------------------
# interactive_pca tests
# ---------------------------------------------------------------------------


class TestInteractivePCA:
    """Tests for the interactive_pca function."""

    def test_pca_with_metadata(self, tmp_path: Path) -> None:
        """3D PCA with 30 samples and population metadata."""
        pca_data = _make_pca_data(30)
        metadata = {}
        populations = ["PopA", "PopB", "PopC"]
        for i, sid in enumerate(pca_data["sample_ids"]):
            metadata[sid] = {"population": populations[i % len(populations)]}

        output = tmp_path / "pca_3d.html"
        result = interactive_pca(pca_data, metadata=metadata, output_file=output, color_by="population")

        assert result["status"] == "success"
        assert result["n_samples"] == 30
        assert output.exists()
        assert output.stat().st_size > 0

        content = output.read_text(encoding="utf-8")
        assert "<html" in content.lower() or "plotly" in content.lower() or "Sample" in content

    def test_pca_without_metadata(self, tmp_path: Path) -> None:
        """3D PCA without metadata still succeeds."""
        pca_data = _make_pca_data(20)
        output = tmp_path / "pca_no_meta.html"
        result = interactive_pca(pca_data, metadata=None, output_file=output)

        assert result["status"] == "success"
        assert result["n_samples"] == 20
        assert output.exists()
        assert output.stat().st_size > 0

    def test_pca_empty_data(self, tmp_path: Path) -> None:
        """Empty PCA data returns skipped status."""
        output = tmp_path / "pca_empty.html"
        result = interactive_pca({"pcs": [], "explained_variance_ratio": []}, output_file=output)

        assert result["status"] == "skipped"
        assert result["n_samples"] == 0

    def test_pca_insufficient_components(self, tmp_path: Path) -> None:
        """PCA with fewer than 3 components returns failed."""
        pca_data = {
            "pcs": [[1.0, 2.0], [3.0, 4.0]],
            "explained_variance_ratio": [0.6, 0.3],
            "sample_ids": ["S1", "S2"],
        }
        output = tmp_path / "pca_2pc.html"
        result = interactive_pca(pca_data, output_file=output)

        assert result["status"] == "failed"
        assert result["n_samples"] == 2


# ---------------------------------------------------------------------------
# interactive_volcano tests
# ---------------------------------------------------------------------------


class TestInteractiveVolcano:
    """Tests for the interactive_volcano function."""

    def test_basic_volcano(self, tmp_path: Path) -> None:
        """Volcano plot with 100 synthetic results produces valid HTML."""
        results = _make_assoc_results(100)
        output = tmp_path / "volcano.html"
        result = interactive_volcano(results, output_file=output)

        assert result["status"] == "success"
        assert result["n_variants"] == 100
        assert output.exists()
        assert output.stat().st_size > 0

        content = output.read_text(encoding="utf-8")
        assert "<html" in content.lower() or "plotly" in content.lower()

    def test_volcano_empty_results(self, tmp_path: Path) -> None:
        """Empty results return skipped status without error."""
        output = tmp_path / "volcano_empty.html"
        result = interactive_volcano([], output_file=output)

        assert result["status"] == "skipped"
        assert result["n_variants"] == 0

    def test_volcano_custom_thresholds(self, tmp_path: Path) -> None:
        """Volcano plot with custom significance and fold-change thresholds."""
        results = _make_assoc_results(100)
        output = tmp_path / "volcano_custom.html"
        result = interactive_volcano(
            results,
            output_file=output,
            significance_threshold=1e-5,
            fold_change_threshold=1.0,
        )

        assert result["status"] == "success"
        assert result["n_variants"] == 100
        assert output.exists()

    def test_volcano_no_output_file(self) -> None:
        """Volcano plot without output file returns success and None path."""
        results = _make_assoc_results(50)
        result = interactive_volcano(results, output_file=None)

        assert result["status"] == "success"
        assert result["n_variants"] == 50
        assert result["output_path"] is None
