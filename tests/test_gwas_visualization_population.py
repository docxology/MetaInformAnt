"""Tests for GWAS population structure visualization functions."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from metainformant.gwas.visualization_population import (
    pca_plot,
    pca_scree_plot,
    kinship_heatmap,
    admixture_plot,
)


def test_pca_plot_basic(tmp_path: Path) -> None:
    """Test basic PCA plot generation."""
    # Create sample PCA data file
    pca_file = tmp_path / "pca_data.tsv"
    with open(pca_file, "w") as f:
        f.write("sample_id\tPC1\tPC2\tPC3\n")
        for i in range(100):
            pc1, pc2, pc3 = np.random.randn(3)
            f.write(f"sample_{i}\t{pc1}\t{pc2}\t{pc3}\n")

    output_path = tmp_path / "pca_plot.png"
    result = pca_plot(pca_file, output_path)

    assert result["status"] == "success"
    assert output_path.exists()
    assert output_path.stat().st_size > 0


def test_pca_scree_plot(tmp_path: Path) -> None:
    """Test PCA scree plot generation."""
    # Create sample PCA eigenvalues file
    eigenval_file = tmp_path / "eigenvals.txt"
    eigenvals = [0.3, 0.2, 0.15, 0.1, 0.05]
    with open(eigenval_file, "w") as f:
        for val in eigenvals:
            f.write(f"{val}\n")

    output_path = tmp_path / "scree_plot.png"
    result = pca_scree_plot(eigenval_file, output_path)

    assert result["status"] == "success"
    assert output_path.exists()


def test_admixture_plot(tmp_path: Path) -> None:
    """Test admixture plot generation."""
    # Create sample admixture file (ADMIXTURE format)
    admixture_file = tmp_path / "admixture.txt"
    with open(admixture_file, "w") as f:
        for i in range(50):
            # Three ancestry components that sum to 1
            q1, q2 = np.random.rand(2)
            q3 = 1 - q1 - q2
            f.write(f"{q1}\t{q2}\t{q3}\n")

    output_path = tmp_path / "admixture.png"
    result = admixture_plot(admixture_file, output_path)

    assert result["status"] == "success"
    assert output_path.exists()


def test_kinship_heatmap_with_title(tmp_path: Path) -> None:
    """Test kinship heatmap with custom title."""
    # Create sample kinship matrix file
    kinship_file = tmp_path / "kinship2.tsv"
    n_samples = 20
    with open(kinship_file, "w") as f:
        sample_ids = [f"s{i}" for i in range(n_samples)]
        f.write("\t".join(sample_ids) + "\n")

        for i in range(n_samples):
            row = [f"{1.0 if i==j else np.random.uniform(0, 0.3)}" for j in range(n_samples)]
            f.write("\t".join(row) + "\n")

    output_path = tmp_path / "kinship_titled.png"
    result = kinship_heatmap(kinship_file, output_path, title="Sample Kinship Matrix")

    assert result["status"] == "success"
    assert output_path.exists()


def test_kinship_heatmap(tmp_path: Path) -> None:
    """Test kinship matrix heatmap generation."""
    # Create sample kinship matrix file
    kinship_file = tmp_path / "kinship.tsv"
    n_samples = 30
    with open(kinship_file, "w") as f:
        # Header with sample IDs
        sample_ids = [f"sample_{i}" for i in range(n_samples)]
        f.write("\t".join(sample_ids) + "\n")

        # Create symmetric kinship matrix
        for i in range(n_samples):
            row = []
            for j in range(n_samples):
                if i == j:
                    row.append("1.0")  # Diagonal
                elif j > i:
                    val = np.random.uniform(0, 0.5)  # Upper triangle
                    row.append(f"{val}")
                else:
                    # For lower triangle, copy from upper triangle to maintain symmetry
                    val = np.random.uniform(0, 0.5)
                    row.append(f"{val}")
            f.write("\t".join(row) + "\n")

    output_path = tmp_path / "kinship.png"
    result = kinship_heatmap(kinship_file, output_path)

    assert result["status"] == "success"
    assert output_path.exists()


def test_kinship_heatmap_with_populations(tmp_path: Path) -> None:
    """Test kinship heatmap with population annotations."""
    n_samples = 40
    kinship_matrix = np.random.rand(n_samples, n_samples)
    kinship_matrix = (kinship_matrix + kinship_matrix.T) / 2
    np.fill_diagonal(kinship_matrix, 1.0)

    sample_labels = [f"S{i}" for i in range(n_samples)]
    population_labels = ["POP_A"] * 20 + ["POP_B"] * 20

    output_path = tmp_path / "kinship_pop.png"
    result = kinship_heatmap(kinship_matrix, sample_labels, output_path,
                           population_labels=population_labels)

    assert result["status"] == "success"
    assert output_path.exists()


def test_admixture_plot_with_k(tmp_path: Path) -> None:
    """Test admixture plot with specified K."""
    # Create sample admixture file
    admixture_file = tmp_path / "admixture_k4.txt"
    with open(admixture_file, "w") as f:
        for i in range(40):
            # Four ancestry components that sum to 1
            q_vals = np.random.rand(4)
            q_vals = q_vals / q_vals.sum()
            f.write("\t".join(f"{q:.4f}" for q in q_vals) + "\n")

    output_path = tmp_path / "admixture_k4.png"
    result = admixture_plot(admixture_file, output_path, k=4)

    assert result["status"] == "success"
    assert output_path.exists()
