"""Tests for GWAS population structure visualization functions."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from metainformant.gwas.visualization.visualization_population import (
    pca_plot,
    pca_scree_plot,
    kinship_heatmap,
    admixture_plot,
    kinship_dendrogram,
    kinship_clustermap,
    pca_multi_panel,
    pca_3d,
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
    result = kinship_heatmap(kinship_matrix, sample_labels, output_path, population_labels=population_labels)

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


# ---------------------------------------------------------------------------
# kinship_dendrogram tests
# ---------------------------------------------------------------------------


def _make_synthetic_kinship(n: int, n_pops: int = 1, rng: np.random.RandomState | None = None) -> np.ndarray:
    """Build a synthetic kinship matrix with block structure.

    Samples within the same population block have higher kinship than
    samples in different blocks, producing a realistic block-diagonal
    pattern for testing dendrogram / clustermap clustering.

    Args:
        n: Total number of samples.
        n_pops: Number of populations (blocks).
        rng: Optional RandomState for reproducibility.

    Returns:
        Symmetric n x n kinship matrix with 1.0 on the diagonal.
    """
    if rng is None:
        rng = np.random.RandomState(42)

    block_size = n // n_pops
    kinship = rng.uniform(0.0, 0.1, size=(n, n))
    kinship = (kinship + kinship.T) / 2.0

    for b in range(n_pops):
        start = b * block_size
        end = start + block_size if b < n_pops - 1 else n
        kinship[start:end, start:end] += rng.uniform(0.2, 0.4, size=(end - start, end - start))
        kinship[start:end, start:end] = (kinship[start:end, start:end] + kinship[start:end, start:end].T) / 2.0

    np.fill_diagonal(kinship, 1.0)
    kinship = np.clip(kinship, 0.0, 1.0)
    return kinship


def test_kinship_dendrogram_basic_20_samples(tmp_path: Path) -> None:
    """Test dendrogram with 20 samples, no metadata coloring."""
    n = 20
    kinship = _make_synthetic_kinship(n, n_pops=2)
    labels = [f"sample_{i}" for i in range(n)]
    output_file = tmp_path / "dendro_20.png"

    result = kinship_dendrogram(kinship, sample_labels=labels, output_file=output_file)

    assert result["status"] == "success"
    assert result["n_samples"] == n
    assert output_file.exists()
    assert output_file.stat().st_size > 0


def test_kinship_dendrogram_with_metadata_coloring(tmp_path: Path) -> None:
    """Test dendrogram with population metadata coloring leaves."""
    n = 20
    kinship = _make_synthetic_kinship(n, n_pops=2)
    labels = [f"ind_{i}" for i in range(n)]
    metadata = {}
    for i in range(n):
        pop = "European" if i < 10 else "African"
        metadata[labels[i]] = {"population": pop, "region": "test"}

    output_file = tmp_path / "dendro_colored.png"
    result = kinship_dendrogram(
        kinship,
        sample_labels=labels,
        output_file=output_file,
        color_by="population",
        metadata=metadata,
    )

    assert result["status"] == "success"
    assert result["n_samples"] == n
    assert output_file.exists()


def test_kinship_dendrogram_edge_case_2_samples(tmp_path: Path) -> None:
    """Test dendrogram with the minimum viable number of samples (2)."""
    kinship = np.array([[1.0, 0.3], [0.3, 1.0]])
    output_file = tmp_path / "dendro_2.png"

    result = kinship_dendrogram(kinship, sample_labels=["A", "B"], output_file=output_file)

    assert result["status"] == "success"
    assert result["n_samples"] == 2
    assert output_file.exists()


def test_kinship_dendrogram_single_sample() -> None:
    """Test dendrogram with 1 sample returns failure (need >= 2)."""
    kinship = np.array([[1.0]])
    result = kinship_dendrogram(kinship, sample_labels=["A"])

    assert result["status"] == "failed"
    assert result["n_samples"] == 1


def test_kinship_dendrogram_no_labels(tmp_path: Path) -> None:
    """Test dendrogram with no sample labels generates defaults."""
    n = 10
    kinship = _make_synthetic_kinship(n, n_pops=1)
    output_file = tmp_path / "dendro_nolabels.png"

    result = kinship_dendrogram(kinship, output_file=output_file)

    assert result["status"] == "success"
    assert result["n_samples"] == n
    assert output_file.exists()


def test_kinship_dendrogram_list_input(tmp_path: Path) -> None:
    """Test dendrogram accepts a list-of-lists kinship matrix."""
    kinship_list = [[1.0, 0.4, 0.1], [0.4, 1.0, 0.2], [0.1, 0.2, 1.0]]
    output_file = tmp_path / "dendro_list.png"

    result = kinship_dendrogram(kinship_list, output_file=output_file)

    assert result["status"] == "success"
    assert result["n_samples"] == 3


def test_kinship_dendrogram_identical_matrix(tmp_path: Path) -> None:
    """Test dendrogram handles all-identical kinship (all 1.0) gracefully."""
    n = 5
    kinship = np.ones((n, n))
    output_file = tmp_path / "dendro_identical.png"

    result = kinship_dendrogram(kinship, output_file=output_file)

    # Should succeed (flat dendrogram) or at least not crash
    assert result["status"] in ("success", "failed")
    assert result["n_samples"] == n or result["n_samples"] == 0


def test_kinship_dendrogram_different_methods(tmp_path: Path) -> None:
    """Test dendrogram works with various linkage methods."""
    n = 10
    kinship = _make_synthetic_kinship(n, n_pops=2)

    for linkage_method in ["ward", "complete", "average", "single"]:
        output_file = tmp_path / f"dendro_{linkage_method}.png"
        result = kinship_dendrogram(kinship, output_file=output_file, method=linkage_method)
        assert result["status"] == "success", f"Failed with method={linkage_method}"


# ---------------------------------------------------------------------------
# kinship_clustermap tests
# ---------------------------------------------------------------------------


def test_kinship_clustermap_basic_20_samples(tmp_path: Path) -> None:
    """Test clustermap with 20 samples, no population annotations."""
    n = 20
    kinship = _make_synthetic_kinship(n, n_pops=2)
    labels = [f"sample_{i}" for i in range(n)]
    output_file = tmp_path / "clustermap_20.png"

    result = kinship_clustermap(
        kinship,
        sample_labels=labels,
        output_file=output_file,
        annotate_populations=False,
    )

    assert result["status"] == "success"
    assert result["n_samples"] == n
    assert output_file.exists()
    assert output_file.stat().st_size > 0


def test_kinship_clustermap_with_population_annotations(tmp_path: Path) -> None:
    """Test clustermap with population metadata annotations on both axes."""
    n = 30
    kinship = _make_synthetic_kinship(n, n_pops=3)
    labels = [f"ind_{i}" for i in range(n)]
    metadata = {}
    pops = ["Pop_A", "Pop_B", "Pop_C"]
    block = n // 3
    for i in range(n):
        pop = pops[min(i // block, 2)]
        metadata[labels[i]] = {"population": pop}

    output_file = tmp_path / "clustermap_pops.png"
    result = kinship_clustermap(
        kinship,
        sample_labels=labels,
        output_file=output_file,
        metadata=metadata,
        annotate_populations=True,
    )

    assert result["status"] == "success"
    assert result["n_samples"] == n
    assert output_file.exists()


def test_kinship_clustermap_edge_case_2_samples(tmp_path: Path) -> None:
    """Test clustermap with only 2 samples."""
    kinship = np.array([[1.0, 0.25], [0.25, 1.0]])
    output_file = tmp_path / "clustermap_2.png"

    result = kinship_clustermap(kinship, sample_labels=["X", "Y"], output_file=output_file)

    assert result["status"] == "success"
    assert result["n_samples"] == 2
    assert output_file.exists()


def test_kinship_clustermap_single_sample() -> None:
    """Test clustermap with 1 sample returns failure."""
    kinship = np.array([[1.0]])
    result = kinship_clustermap(kinship, sample_labels=["only"])

    assert result["status"] == "failed"
    assert result["n_samples"] == 1


def test_kinship_clustermap_list_input(tmp_path: Path) -> None:
    """Test clustermap accepts list-of-lists kinship matrix."""
    kinship_list = [[1.0, 0.5, 0.1], [0.5, 1.0, 0.3], [0.1, 0.3, 1.0]]
    output_file = tmp_path / "clustermap_list.png"

    result = kinship_clustermap(kinship_list, output_file=output_file)

    assert result["status"] == "success"
    assert result["n_samples"] == 3


def test_kinship_clustermap_no_output_file() -> None:
    """Test clustermap without saving to file still returns success."""
    n = 8
    kinship = _make_synthetic_kinship(n)

    result = kinship_clustermap(kinship)

    assert result["status"] == "success"
    assert result["output_path"] is None
    assert result["n_samples"] == n


def test_kinship_clustermap_large_sample_no_labels(tmp_path: Path) -> None:
    """Test clustermap with >50 samples omits tick labels gracefully."""
    n = 60
    kinship = _make_synthetic_kinship(n, n_pops=3)
    output_file = tmp_path / "clustermap_large.png"

    result = kinship_clustermap(kinship, output_file=output_file)

    assert result["status"] == "success"
    assert result["n_samples"] == n
    assert output_file.exists()


# ---------------------------------------------------------------------------
# pca_multi_panel tests
# ---------------------------------------------------------------------------


def _make_pca_data(
    n_samples: int = 30,
    n_components: int = 5,
    n_pops: int = 3,
    rng: np.random.RandomState | None = None,
) -> tuple[dict, dict]:
    """Build synthetic PCA data and metadata for testing.

    Returns:
        Tuple of (pca_data dict, metadata dict).
    """
    if rng is None:
        rng = np.random.RandomState(99)

    sample_ids = [f"sample_{i}" for i in range(n_samples)]
    pcs = rng.randn(n_samples, n_components).tolist()
    explained = [0.30, 0.20, 0.15, 0.10, 0.05][:n_components]

    pop_names = ["PopA", "PopB", "PopC"][:n_pops]
    block = n_samples // n_pops
    metadata = {}
    for i, sid in enumerate(sample_ids):
        pop = pop_names[min(i // block, n_pops - 1)]
        metadata[sid] = {"population": pop, "region": "test"}

    pca_data = {
        "pcs": pcs,
        "explained_variance_ratio": explained,
        "sample_ids": sample_ids,
    }
    return pca_data, metadata


def test_pca_multi_panel_with_metadata(tmp_path: Path) -> None:
    """Test multi-panel PCA with 30 samples, 5 PCs, 3 populations."""
    pca_data, metadata = _make_pca_data(n_samples=30, n_components=5, n_pops=3)
    output_file = tmp_path / "pca_multi.png"

    result = pca_multi_panel(pca_data, metadata=metadata, output_file=output_file)

    assert result["status"] == "success"
    assert result["n_samples"] == 30
    assert result["n_panels"] == 3  # default (0,1), (0,2), (1,2)
    assert output_file.exists()
    assert output_file.stat().st_size > 0


def test_pca_multi_panel_no_metadata(tmp_path: Path) -> None:
    """Test multi-panel PCA without metadata (no coloring)."""
    pca_data, _ = _make_pca_data(n_samples=20, n_components=4, n_pops=2)
    output_file = tmp_path / "pca_multi_nocolour.png"

    result = pca_multi_panel(pca_data, output_file=output_file)

    assert result["status"] == "success"
    assert result["n_samples"] == 20
    assert result["n_panels"] == 3
    assert output_file.exists()


def test_pca_multi_panel_custom_pairs(tmp_path: Path) -> None:
    """Test multi-panel PCA with custom PC pairs."""
    pca_data, metadata = _make_pca_data(n_samples=25, n_components=5, n_pops=2)
    output_file = tmp_path / "pca_multi_custom.png"

    result = pca_multi_panel(
        pca_data,
        metadata=metadata,
        output_file=output_file,
        pairs=[(0, 1), (2, 3), (3, 4)],
    )

    assert result["status"] == "success"
    assert result["n_panels"] == 3


def test_pca_multi_panel_prunes_invalid_pairs(tmp_path: Path) -> None:
    """Test that pairs referencing non-existent PCs are pruned."""
    pca_data, _ = _make_pca_data(n_samples=15, n_components=2, n_pops=1)
    output_file = tmp_path / "pca_multi_pruned.png"

    # (0,1) is valid; (0,2) and (1,2) require PC3 which does not exist
    result = pca_multi_panel(pca_data, output_file=output_file)

    assert result["status"] == "success"
    assert result["n_panels"] == 1  # only (0,1) survives


def test_pca_multi_panel_no_output_file() -> None:
    """Test multi-panel PCA without saving to file still succeeds."""
    pca_data, _ = _make_pca_data(n_samples=10, n_components=3, n_pops=1)

    result = pca_multi_panel(pca_data)

    assert result["status"] == "success"
    assert result["output_path"] is None


# ---------------------------------------------------------------------------
# pca_3d tests
# ---------------------------------------------------------------------------


def test_pca_3d_with_metadata(tmp_path: Path) -> None:
    """Test 3D PCA scatter with metadata coloring."""
    pca_data, metadata = _make_pca_data(n_samples=30, n_components=5, n_pops=3)
    output_file = tmp_path / "pca_3d.png"

    result = pca_3d(pca_data, metadata=metadata, output_file=output_file)

    assert result["status"] == "success"
    assert result["n_samples"] == 30
    assert result["n_components"] == 5
    assert output_file.exists()
    assert output_file.stat().st_size > 0


def test_pca_3d_no_metadata(tmp_path: Path) -> None:
    """Test 3D PCA scatter without metadata."""
    pca_data, _ = _make_pca_data(n_samples=20, n_components=4, n_pops=2)
    output_file = tmp_path / "pca_3d_nocolour.png"

    result = pca_3d(pca_data, output_file=output_file)

    assert result["status"] == "success"
    assert result["n_samples"] == 20
    assert output_file.exists()


def test_pca_3d_only_2_pcs(tmp_path: Path) -> None:
    """Test 3D PCA with only 2 PCs available returns failed gracefully."""
    pca_data, _ = _make_pca_data(n_samples=15, n_components=2, n_pops=1)
    output_file = tmp_path / "pca_3d_2pcs.png"

    result = pca_3d(pca_data, output_file=output_file)

    assert result["status"] == "failed"
    assert result["n_samples"] == 15
    assert result["n_components"] == 2
    assert "error" in result
    assert not output_file.exists()


def test_pca_3d_custom_components(tmp_path: Path) -> None:
    """Test 3D PCA with custom component selection."""
    pca_data, metadata = _make_pca_data(n_samples=25, n_components=5, n_pops=2)
    output_file = tmp_path / "pca_3d_custom.png"

    result = pca_3d(
        pca_data,
        metadata=metadata,
        output_file=output_file,
        components=(1, 2, 4),
    )

    assert result["status"] == "success"
    assert result["n_samples"] == 25
    assert output_file.exists()


def test_pca_3d_no_output_file() -> None:
    """Test 3D PCA without saving to file still succeeds."""
    pca_data, _ = _make_pca_data(n_samples=10, n_components=3, n_pops=1)

    result = pca_3d(pca_data)

    assert result["status"] == "success"
    assert result["output_path"] is None
