"""Tests for single-cell RNA velocity estimation and analysis.

Real implementation testing for velocity computation, embedding projection,
pseudotime estimation, confidence metrics, and dynamical model fitting.
No mocking used - all tests use real computational methods and data.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from metainformant.singlecell.velocity.rna_velocity import (
    compute_velocity,
    fit_dynamical_model,
    velocity_confidence,
    velocity_embedding,
    velocity_pseudotime,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_spliced_unspliced(
    n_cells: int = 50,
    n_genes: int = 30,
    seed: int = 42,
) -> tuple[list[list[float]], list[list[float]], list[str]]:
    """Create synthetic spliced/unspliced count matrices.

    Generates data where unspliced ~= gamma * spliced + noise so that
    a linear fit can recover velocity genes.
    """
    rng = np.random.RandomState(seed)
    gene_names = [f"gene_{i}" for i in range(n_genes)]

    # Base spliced counts
    spliced = rng.exponential(3.0, size=(n_cells, n_genes))

    # Unspliced proportional to spliced with gene-specific gamma
    gammas = rng.uniform(0.1, 0.8, size=n_genes)
    unspliced = np.zeros_like(spliced)
    for g in range(n_genes):
        unspliced[:, g] = gammas[g] * spliced[:, g] + rng.normal(0, 0.2, size=n_cells)
    unspliced = np.clip(unspliced, 0, None)

    # Make a few genes pure noise (low counts) so they get filtered
    spliced[:, -3:] = rng.exponential(0.01, size=(n_cells, 3))
    unspliced[:, -3:] = rng.exponential(0.01, size=(n_cells, 3))

    return spliced.tolist(), unspliced.tolist(), gene_names


def _make_embedding(n_cells: int = 50, seed: int = 42) -> list[list[float]]:
    """Create a 2D embedding (e.g. UMAP coordinates)."""
    rng = np.random.RandomState(seed)
    t = np.linspace(0, 2 * np.pi, n_cells)
    x = np.cos(t) + rng.normal(0, 0.1, n_cells)
    y = np.sin(t) + rng.normal(0, 0.1, n_cells)
    return [[float(x[i]), float(y[i])] for i in range(n_cells)]


# ---------------------------------------------------------------------------
# compute_velocity
# ---------------------------------------------------------------------------


class TestComputeVelocity:
    """Tests for compute_velocity function."""

    def test_returns_expected_keys(self) -> None:
        s, u, gn = _make_spliced_unspliced()
        result = compute_velocity(s, u, gn)
        for key in [
            "velocity_matrix",
            "gamma",
            "r_squared",
            "velocity_genes",
            "n_velocity_genes",
        ]:
            assert key in result

    def test_velocity_matrix_dimensions(self) -> None:
        n_cells = 40
        s, u, gn = _make_spliced_unspliced(n_cells=n_cells)
        result = compute_velocity(s, u, gn)
        assert len(result["velocity_matrix"]) == n_cells
        assert len(result["velocity_matrix"][0]) == result["n_velocity_genes"]

    def test_velocity_genes_subset(self) -> None:
        s, u, gn = _make_spliced_unspliced()
        result = compute_velocity(s, u, gn)
        assert result["n_velocity_genes"] <= len(gn)
        assert result["n_velocity_genes"] > 0
        for vg in result["velocity_genes"]:
            assert vg in gn

    def test_gamma_positive(self) -> None:
        s, u, gn = _make_spliced_unspliced()
        result = compute_velocity(s, u, gn)
        for gene in result["velocity_genes"]:
            assert result["gamma"][gene] > 0

    def test_r_squared_in_range(self) -> None:
        s, u, gn = _make_spliced_unspliced()
        result = compute_velocity(s, u, gn)
        for gene, r2 in result["r_squared"].items():
            assert 0.0 <= r2 <= 1.0

    def test_dimension_mismatch_spliced_unspliced_rows(self) -> None:
        s = [[1.0, 2.0], [3.0, 4.0]]
        u = [[1.0, 2.0]]  # wrong row count
        with pytest.raises(ValueError, match="same number of rows"):
            compute_velocity(s, u, ["a", "b"])

    def test_dimension_mismatch_spliced_unspliced_cols(self) -> None:
        s = [[1.0, 2.0], [3.0, 4.0]]
        u = [[1.0], [3.0]]  # wrong col count
        with pytest.raises(ValueError, match="same number of columns"):
            compute_velocity(s, u, ["a", "b"])

    def test_gene_names_mismatch_raises(self) -> None:
        s = [[1.0, 2.0]]
        u = [[0.5, 1.0]]
        with pytest.raises(ValueError, match="must match"):
            compute_velocity(s, u, ["only_one"])

    def test_invalid_method_raises(self) -> None:
        s, u, gn = _make_spliced_unspliced(n_cells=10, n_genes=5)
        with pytest.raises(ValueError, match="Invalid method"):
            compute_velocity(s, u, gn, method="dynamical")

    def test_min_counts_filter(self) -> None:
        """Genes with very low counts should be filtered out."""
        s, u, gn = _make_spliced_unspliced()
        # Very high min_counts should filter most genes
        result = compute_velocity(s, u, gn, min_counts=100000)
        assert result["n_velocity_genes"] == 0

    def test_r_squared_threshold(self) -> None:
        s, u, gn = _make_spliced_unspliced()
        # Very high threshold should filter most genes
        result_strict = compute_velocity(s, u, gn, r_squared_threshold=0.99)
        result_loose = compute_velocity(s, u, gn, r_squared_threshold=0.001)
        assert result_strict["n_velocity_genes"] <= result_loose["n_velocity_genes"]

    def test_numpy_array_input(self) -> None:
        s, u, gn = _make_spliced_unspliced(n_cells=30)
        result = compute_velocity(np.array(s), np.array(u), gn)
        assert len(result["velocity_matrix"]) == 30


# ---------------------------------------------------------------------------
# velocity_embedding
# ---------------------------------------------------------------------------


class TestVelocityEmbedding:
    """Tests for velocity_embedding projection."""

    def test_returns_expected_keys(self) -> None:
        s, u, gn = _make_spliced_unspliced(n_cells=30)
        vel_result = compute_velocity(s, u, gn)
        emb = _make_embedding(n_cells=30)
        result = velocity_embedding(vel_result["velocity_matrix"], emb, n_neighbors=10)
        for key in [
            "velocity_embedding",
            "transition_matrix",
            "cell_velocities",
        ]:
            assert key in result

    def test_velocity_embedding_dimensions(self) -> None:
        n_cells = 25
        s, u, gn = _make_spliced_unspliced(n_cells=n_cells)
        vel_result = compute_velocity(s, u, gn)
        emb = _make_embedding(n_cells=n_cells)
        result = velocity_embedding(vel_result["velocity_matrix"], emb, n_neighbors=5)
        assert len(result["velocity_embedding"]) == n_cells
        # Each projected velocity should have 2 dims (matching embedding)
        assert len(result["velocity_embedding"][0]) == 2

    def test_transition_matrix_dimensions(self) -> None:
        n_cells = 20
        s, u, gn = _make_spliced_unspliced(n_cells=n_cells)
        vel_result = compute_velocity(s, u, gn)
        emb = _make_embedding(n_cells=n_cells)
        result = velocity_embedding(vel_result["velocity_matrix"], emb, n_neighbors=5)
        tm = result["transition_matrix"]
        assert len(tm) == n_cells
        assert len(tm[0]) == n_cells

    def test_transition_probabilities_non_negative(self) -> None:
        n_cells = 20
        s, u, gn = _make_spliced_unspliced(n_cells=n_cells)
        vel_result = compute_velocity(s, u, gn)
        emb = _make_embedding(n_cells=n_cells)
        result = velocity_embedding(vel_result["velocity_matrix"], emb, n_neighbors=5)
        for row in result["transition_matrix"]:
            for val in row:
                assert val >= 0.0

    def test_cell_velocities_non_negative(self) -> None:
        n_cells = 25
        s, u, gn = _make_spliced_unspliced(n_cells=n_cells)
        vel_result = compute_velocity(s, u, gn)
        emb = _make_embedding(n_cells=n_cells)
        result = velocity_embedding(vel_result["velocity_matrix"], emb, n_neighbors=5)
        for v in result["cell_velocities"]:
            assert v >= 0.0

    def test_dimension_mismatch_raises(self) -> None:
        velocity = [[1.0, 2.0], [3.0, 4.0]]
        embedding = [[0.1, 0.2]]  # wrong n_cells
        with pytest.raises(ValueError, match="same number of rows"):
            velocity_embedding(velocity, embedding)


# ---------------------------------------------------------------------------
# velocity_pseudotime
# ---------------------------------------------------------------------------


class TestVelocityPseudotime:
    """Tests for velocity_pseudotime."""

    def _run_pseudotime(self, n_cells: int = 25, root: int | None = None) -> dict:
        s, u, gn = _make_spliced_unspliced(n_cells=n_cells)
        vel_result = compute_velocity(s, u, gn)
        emb = _make_embedding(n_cells=n_cells)
        return velocity_pseudotime(
            vel_result["velocity_matrix"],
            emb,
            root_cell=root,
            n_neighbors=5,
        )

    def test_returns_expected_keys(self) -> None:
        result = self._run_pseudotime()
        assert "pseudotime" in result
        assert "root_cell" in result
        assert "terminal_cells" in result

    def test_pseudotime_normalized(self) -> None:
        result = self._run_pseudotime()
        pt = result["pseudotime"]
        assert min(pt) >= 0.0
        assert max(pt) <= 1.0

    def test_pseudotime_length(self) -> None:
        n = 30
        result = self._run_pseudotime(n_cells=n)
        assert len(result["pseudotime"]) == n

    def test_root_cell_has_zero_pseudotime(self) -> None:
        result = self._run_pseudotime(root=0)
        assert result["root_cell"] == 0
        assert result["pseudotime"][0] == 0.0

    def test_auto_root_selection(self) -> None:
        result = self._run_pseudotime()
        root = result["root_cell"]
        assert 0 <= root < 25

    def test_terminal_cells_at_high_pseudotime(self) -> None:
        result = self._run_pseudotime()
        for tc in result["terminal_cells"]:
            assert result["pseudotime"][tc] >= 0.95

    def test_invalid_root_raises(self) -> None:
        s, u, gn = _make_spliced_unspliced(n_cells=10)
        vel = compute_velocity(s, u, gn)
        emb = _make_embedding(n_cells=10)
        with pytest.raises(ValueError, match="out of range"):
            velocity_pseudotime(vel["velocity_matrix"], emb, root_cell=999)

    def test_dimension_mismatch_raises(self) -> None:
        with pytest.raises(ValueError, match="same number of rows"):
            velocity_pseudotime([[1.0]], [[0.1, 0.2], [0.3, 0.4]])


# ---------------------------------------------------------------------------
# velocity_confidence
# ---------------------------------------------------------------------------


class TestVelocityConfidence:
    """Tests for velocity_confidence metrics."""

    def test_returns_expected_keys(self) -> None:
        s, u, gn = _make_spliced_unspliced(n_cells=30)
        vel = compute_velocity(s, u, gn)
        result = velocity_confidence(vel["velocity_matrix"], s)
        for key in [
            "gene_confidence",
            "cell_confidence",
            "overall_confidence",
            "n_high_confidence_genes",
            "n_high_confidence_cells",
        ]:
            assert key in result

    def test_gene_confidence_range(self) -> None:
        s, u, gn = _make_spliced_unspliced(n_cells=30)
        vel = compute_velocity(s, u, gn)
        result = velocity_confidence(vel["velocity_matrix"], s)
        for c in result["gene_confidence"]:
            assert 0.0 <= c <= 1.0

    def test_cell_confidence_range(self) -> None:
        s, u, gn = _make_spliced_unspliced(n_cells=30)
        vel = compute_velocity(s, u, gn)
        result = velocity_confidence(vel["velocity_matrix"], s)
        for c in result["cell_confidence"]:
            assert 0.0 <= c <= 1.0

    def test_overall_confidence_is_mean(self) -> None:
        s, u, gn = _make_spliced_unspliced(n_cells=20)
        vel = compute_velocity(s, u, gn)
        result = velocity_confidence(vel["velocity_matrix"], s)
        expected = sum(result["cell_confidence"]) / len(result["cell_confidence"])
        assert abs(result["overall_confidence"] - expected) < 1e-9

    def test_gene_confidence_count(self) -> None:
        s, u, gn = _make_spliced_unspliced(n_cells=30)
        vel = compute_velocity(s, u, gn)
        result = velocity_confidence(vel["velocity_matrix"], s)
        n_vel_genes = vel["n_velocity_genes"]
        assert len(result["gene_confidence"]) == n_vel_genes

    def test_cell_confidence_count(self) -> None:
        n_cells = 25
        s, u, gn = _make_spliced_unspliced(n_cells=n_cells)
        vel = compute_velocity(s, u, gn)
        result = velocity_confidence(vel["velocity_matrix"], s)
        assert len(result["cell_confidence"]) == n_cells


# ---------------------------------------------------------------------------
# fit_dynamical_model
# ---------------------------------------------------------------------------


class TestFitDynamicalModel:
    """Tests for fit_dynamical_model."""

    def test_returns_expected_keys(self) -> None:
        s, u, gn = _make_spliced_unspliced(n_cells=30, n_genes=15)
        result = fit_dynamical_model(s, u, gn, max_iter=20)
        for key in [
            "alpha",
            "beta",
            "gamma",
            "likelihood",
            "velocity_matrix",
            "fitted_genes",
            "n_iterations",
        ]:
            assert key in result

    def test_kinetic_parameters_positive(self) -> None:
        s, u, gn = _make_spliced_unspliced(n_cells=30, n_genes=10)
        result = fit_dynamical_model(s, u, gn, max_iter=30)
        for gene in result["fitted_genes"]:
            assert result["alpha"][gene] > 0
            assert result["beta"][gene] > 0
            assert result["gamma"][gene] > 0

    def test_velocity_matrix_dimensions(self) -> None:
        n_cells = 25
        s, u, gn = _make_spliced_unspliced(n_cells=n_cells, n_genes=10)
        result = fit_dynamical_model(s, u, gn, max_iter=20)
        assert len(result["velocity_matrix"]) == n_cells
        if result["fitted_genes"]:
            assert len(result["velocity_matrix"][0]) == len(result["fitted_genes"])

    def test_iterations_within_max(self) -> None:
        s, u, gn = _make_spliced_unspliced(n_cells=20, n_genes=8)
        max_iter = 15
        result = fit_dynamical_model(s, u, gn, max_iter=max_iter)
        for gene, iters in result["n_iterations"].items():
            assert iters <= max_iter

    def test_dimension_mismatch_raises(self) -> None:
        s = [[1.0, 2.0]]
        u = [[1.0]]  # wrong cols
        with pytest.raises(ValueError, match="same number of columns"):
            fit_dynamical_model(s, u, ["a", "b"])

    def test_gene_names_mismatch_raises(self) -> None:
        s = [[1.0, 2.0]]
        u = [[0.5, 1.0]]
        with pytest.raises(ValueError, match="must match"):
            fit_dynamical_model(s, u, ["only_one"])

    def test_min_counts_filter(self) -> None:
        s, u, gn = _make_spliced_unspliced(n_cells=20, n_genes=10)
        result = fit_dynamical_model(s, u, gn, min_counts=1000000)
        assert len(result["fitted_genes"]) == 0

    def test_convergence_tolerance(self) -> None:
        """Loose tolerance should converge faster."""
        s, u, gn = _make_spliced_unspliced(n_cells=30, n_genes=10)
        result_loose = fit_dynamical_model(s, u, gn, max_iter=100, tol=0.1)
        result_tight = fit_dynamical_model(s, u, gn, max_iter=100, tol=1e-8)
        # Loose should use fewer or equal iterations on average
        if result_loose["n_iterations"] and result_tight["n_iterations"]:
            avg_loose = sum(result_loose["n_iterations"].values()) / len(result_loose["n_iterations"])
            avg_tight = sum(result_tight["n_iterations"].values()) / len(result_tight["n_iterations"])
            assert avg_loose <= avg_tight + 1  # allow small tolerance
