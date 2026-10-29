"""Regression tests: ``sv_population_structure`` returns the same quantities for any shape of data.

With fewer samples than SVs the function worked from the ``n_samples x n_samples`` matrix ``X Xt`` and divided it
by ``n_svs - 1`` instead of ``n_samples - 1``, so its eigenvalues were those of the covariance matrix scaled by
``(n_samples - 1) / (n_svs - 1)``; and ``pcs`` was the matrix of unit-length eigenvectors there (``U``) but the
sample scores ``X V`` (``U S``) when there were at least as many samples as SVs.

Both shapes now return the sample scores and the eigenvalues of the sample covariance matrix. The oracle is
numpy's singular value decomposition of the standardised matrix ``X = U S Vt``: scores ``U S``, eigenvalues
``S**2 / (n_samples - 1)``, loadings ``V``. Eigenvectors have an arbitrary sign (it depends on the LAPACK build),
so components are compared up to sign, with tolerances, and the components beyond the rank of the data (whose
direction is undefined) are checked through the rank, never through whether a solver raised.
"""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.structural_variants.population.sv_population import sv_population_structure


def _genotypes(n_svs: int, n_samples: int, seed: int) -> np.ndarray:
    """Genotypes (0/1/2) of two groups of samples with different allele frequencies."""
    rng = np.random.default_rng(seed)
    group = np.arange(n_samples) % 2
    freq = rng.uniform(0.05, 0.6, size=(n_svs, 2))
    alleles = rng.random((2, n_svs, n_samples)) < freq[:, group][None, :, :]
    return alleles.sum(axis=0).astype(int)


def _standardised(genotypes: np.ndarray) -> np.ndarray:
    """Samples x SVs, each SV centred and scaled to unit population variance (constant SVs left at 0)."""
    data = np.asarray(genotypes, dtype=float).T
    std = data.std(axis=0)
    std[std == 0] = 1.0
    return (data - data.mean(axis=0)) / std


def _svd_oracle(genotypes: np.ndarray) -> dict[str, np.ndarray | int]:
    x = _standardised(genotypes)
    u, s, vt = np.linalg.svd(x, full_matrices=False)
    return {
        "scores": u * s,
        "eigenvalues": s**2 / max(x.shape[0] - 1, 1),
        "loadings": vt.T,
        "rank": int(np.linalg.matrix_rank(x)),
    }


def _match_signs(actual: np.ndarray, expected: np.ndarray) -> np.ndarray:
    """Flip each column of ``actual`` to the sign that agrees with ``expected``."""
    signs = np.sign(np.sum(actual * expected, axis=0))
    signs[signs == 0] = 1.0
    return actual * signs


# (n_svs, n_samples): more SVs than samples takes the sample-by-sample route, the others the SV-by-SV one.
SHAPES = [(60, 12), (30, 8), (12, 60), (25, 25), (7, 7), (9, 4)]


@pytest.mark.parametrize(("n_svs", "n_samples"), SHAPES)
class TestAgainstTheSingularValueDecomposition:
    def _run(self, n_svs: int, n_samples: int) -> tuple[dict[str, object], dict[str, np.ndarray | int], int]:
        genotypes = _genotypes(n_svs, n_samples, seed=100 * n_svs + n_samples)
        oracle = _svd_oracle(genotypes)
        result = sv_population_structure(genotypes, n_components=min(n_svs, n_samples))
        assert isinstance(oracle["rank"], int)
        return result, oracle, oracle["rank"]

    def test_eigenvalues_are_those_of_the_sample_covariance_matrix(self, n_svs: int, n_samples: int) -> None:
        result, oracle, rank = self._run(n_svs, n_samples)
        expected = oracle["eigenvalues"]
        assert isinstance(expected, np.ndarray)
        assert np.asarray(result["eigenvalues"])[:rank] == pytest.approx(expected[:rank], rel=1e-8, abs=1e-10)

    def test_pcs_are_the_sample_scores(self, n_svs: int, n_samples: int) -> None:
        result, oracle, rank = self._run(n_svs, n_samples)
        scores, expected = np.asarray(result["pcs"])[:, :rank], oracle["scores"]
        assert isinstance(expected, np.ndarray)
        assert _match_signs(scores, expected[:, :rank]) == pytest.approx(expected[:, :rank], rel=1e-6, abs=1e-8)

    def test_loadings_are_the_covariance_eigenvectors(self, n_svs: int, n_samples: int) -> None:
        result, oracle, rank = self._run(n_svs, n_samples)
        loadings, expected = np.asarray(result["loadings"])[:, :rank], oracle["loadings"]
        assert isinstance(expected, np.ndarray)
        assert _match_signs(loadings, expected[:, :rank]) == pytest.approx(expected[:, :rank], rel=1e-6, abs=1e-8)
        assert loadings.T @ loadings == pytest.approx(np.eye(rank), abs=1e-8)

    def test_scores_are_the_data_projected_on_the_loadings(self, n_svs: int, n_samples: int) -> None:
        # The identity that ties the three outputs together, whatever the route taken to compute them.
        genotypes = _genotypes(n_svs, n_samples, seed=100 * n_svs + n_samples)
        result = sv_population_structure(genotypes, n_components=min(n_svs, n_samples))
        assert np.asarray(result["pcs"]) == pytest.approx(
            _standardised(genotypes) @ np.asarray(result["loadings"]), rel=1e-6, abs=1e-8
        )

    def test_the_variance_of_a_score_column_is_its_eigenvalue(self, n_svs: int, n_samples: int) -> None:
        result, _, rank = self._run(n_svs, n_samples)
        pcs = np.asarray(result["pcs"])
        variances = (pcs**2).sum(axis=0) / max(n_samples - 1, 1)
        assert variances[:rank] == pytest.approx(np.asarray(result["eigenvalues"])[:rank], rel=1e-8, abs=1e-10)

    def test_scores_of_different_components_are_orthogonal(self, n_svs: int, n_samples: int) -> None:
        result, _, rank = self._run(n_svs, n_samples)
        pcs = np.asarray(result["pcs"])[:, :rank]
        gram = pcs.T @ pcs
        assert gram == pytest.approx(np.diag(np.diag(gram)), abs=1e-7)

    def test_components_beyond_the_rank_have_no_direction(self, n_svs: int, n_samples: int) -> None:
        result, _, rank = self._run(n_svs, n_samples)
        eigenvalues = np.asarray(result["eigenvalues"])
        assert rank == min(n_svs, n_samples - 1)  # centred data lose one dimension
        assert np.all(eigenvalues[:rank] > 0.0)
        assert np.all(eigenvalues[rank:] == 0.0)
        assert np.all(np.asarray(result["pcs"])[:, rank:] == 0.0)
        assert np.all(np.asarray(result["loadings"])[:, rank:] == 0.0)
        assert list(result["variance_explained"][rank:]) == [0.0] * (len(eigenvalues) - rank)  # type: ignore[index]


class TestTheOldDefectsDirectly:
    def test_eigenvalues_with_fewer_samples_than_svs_are_not_rescaled(self) -> None:
        # 40 SVs x 8 samples: the old eigenvalues were (8 - 1) / (40 - 1) = 0.1795 times the true ones.
        genotypes = _genotypes(40, 8, seed=1)
        result = sv_population_structure(genotypes, n_components=3)
        expected = _svd_oracle(genotypes)["eigenvalues"]
        assert isinstance(expected, np.ndarray)
        assert result["eigenvalues"][0] / expected[0] == pytest.approx(1.0, abs=1e-8)  # type: ignore[index]

    def test_the_sum_of_all_eigenvalues_is_the_trace_of_the_covariance_matrix(self) -> None:
        # A shape-free check on the scale: the eigenvalues of Xt X / (n - 1) add up to its trace, ||X||_F^2 / (n - 1).
        for n_svs, n_samples in [(40, 8), (8, 40), (14, 14)]:
            genotypes = _genotypes(n_svs, n_samples, seed=n_svs)
            result = sv_population_structure(genotypes, n_components=min(n_svs, n_samples))
            trace = float(np.sum(_standardised(genotypes) ** 2)) / (n_samples - 1)
            assert sum(result["eigenvalues"]) == pytest.approx(trace, rel=1e-9)

    def test_pcs_are_not_unit_vectors_when_samples_are_fewer_than_svs(self) -> None:
        genotypes = _genotypes(40, 8, seed=2)
        pcs = np.asarray(sv_population_structure(genotypes, n_components=2)["pcs"])
        # U would have unit columns; the scores U S have squared length (n - 1) * eigenvalue.
        assert np.linalg.norm(pcs[:, 0]) > 2.0

    def test_a_prefix_of_the_components_is_the_same_as_the_full_run(self) -> None:
        genotypes = _genotypes(40, 8, seed=3)
        full = sv_population_structure(genotypes, n_components=7)
        few = sv_population_structure(genotypes, n_components=3)
        assert few["eigenvalues"] == pytest.approx(full["eigenvalues"][:3], rel=1e-9)  # type: ignore[index]
        assert _match_signs(np.asarray(few["pcs"]), np.asarray(full["pcs"])[:, :3]) == pytest.approx(
            np.asarray(full["pcs"])[:, :3], rel=1e-6, abs=1e-8
        )

    def test_permuting_the_samples_permutes_the_scores(self) -> None:
        genotypes = _genotypes(30, 10, seed=4)
        order = np.random.default_rng(0).permutation(10)
        original = np.asarray(sv_population_structure(genotypes, n_components=3)["pcs"])
        permuted = np.asarray(sv_population_structure(genotypes[:, order], n_components=3)["pcs"])
        assert _match_signs(permuted, original[order]) == pytest.approx(original[order], rel=1e-6, abs=1e-8)


class TestDegenerateData:
    def test_constant_svs_add_no_dimensions(self) -> None:
        # 20 SVs of which 8 are constant, 30 samples: the standardised matrix has rank at most 12.
        genotypes = _genotypes(20, 30, seed=5)
        genotypes[:8] = 1
        oracle = _svd_oracle(genotypes)
        rank = oracle["rank"]
        assert isinstance(rank, int) and rank <= 12
        result = sv_population_structure(genotypes, n_components=20)
        eigenvalues = np.asarray(result["eigenvalues"])
        assert np.all(eigenvalues[:rank] > 0.0) and np.all(eigenvalues[rank:] == 0.0)
        assert np.all(np.asarray(result["loadings"])[:8, :rank] == 0.0)  # a constant SV carries no loading

    def test_data_without_variance_has_zero_eigenvalues_and_scores(self) -> None:
        result = sv_population_structure(np.ones((10, 8), dtype=int), n_components=3)
        assert result["eigenvalues"] == [0.0, 0.0, 0.0]
        assert result["variance_explained"] == [0.0, 0.0, 0.0]
        assert np.all(np.asarray(result["pcs"]) == 0.0) and np.all(np.asarray(result["loadings"]) == 0.0)

    def test_a_single_sample_has_no_variance(self) -> None:
        result = sv_population_structure([[0], [1], [2]], n_components=2)
        assert result["eigenvalues"] == [0.0] and np.asarray(result["pcs"]).shape == (1, 1)

    def test_components_are_limited_to_the_smaller_dimension(self) -> None:
        result = sv_population_structure(_genotypes(6, 9, seed=6), n_components=50)
        assert np.asarray(result["pcs"]).shape == (9, 6) and np.asarray(result["loadings"]).shape == (6, 6)


class TestInputs:
    @pytest.mark.parametrize("n_components", [0, -1, -5])
    def test_fewer_than_one_component_is_refused(self, n_components: int) -> None:
        # A negative count used to slice off the last components instead of failing.
        with pytest.raises(ValueError, match="n_components"):
            sv_population_structure(_genotypes(10, 6, seed=7), n_components=n_components)

    @pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf])
    def test_non_finite_genotypes_are_refused(self, bad: float) -> None:
        genotypes = _genotypes(10, 6, seed=8).astype(float)
        genotypes[2, 3] = bad
        with pytest.raises(ValueError, match="non-finite"):
            sv_population_structure(genotypes)

    def test_a_matrix_that_is_not_two_dimensional_is_refused(self) -> None:
        with pytest.raises(ValueError, match="2D"):
            sv_population_structure([0, 1, 2, 1, 0])
