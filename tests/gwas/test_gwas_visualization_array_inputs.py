"""GWAS plots accept NumPy arrays wherever they accept lists.

Zero-mocks: the real plotting functions run on small arrays with the Agg backend.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from metainformant.gwas.visualization.genomic.variants import missingness_plot  # noqa: E402
from metainformant.gwas.visualization.population.population_pca import pca_scree_plot  # noqa: E402


def test_missingness_plot_accepts_a_genotype_array(tmp_path: Path) -> None:
    genotypes = np.array([[0, 1, -1], [2, -1, -1], [1, 1, 0]])  # variants x samples, -1 missing
    result = missingness_plot({"genotypes": genotypes}, output_file=tmp_path / "missing.png")

    assert result["status"] == "success"
    assert (tmp_path / "missing.png").exists()


def test_missingness_plot_reports_an_empty_array(tmp_path: Path) -> None:
    result = missingness_plot({"genotypes": np.empty((0, 0))}, output_file=tmp_path / "missing.png")
    assert result == {"status": "failed", "error": "Empty genotype data"}


def test_scree_plot_accepts_an_explained_variance_array(tmp_path: Path) -> None:
    result = pca_scree_plot(np.array([0.45, 0.25, 0.15, 0.1, 0.05]), output_file=tmp_path / "scree.png")

    assert result["status"] == "success"
    assert (tmp_path / "scree.png").exists()


def test_scree_plot_reports_an_empty_array() -> None:
    result = pca_scree_plot(np.array([]))
    assert result["status"] == "failed"
