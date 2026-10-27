"""Regression tests: ``sv_population_structure`` reports variance explained against the total variance.

``variance_explained`` was each returned eigenvalue divided by the sum of the *returned* eigenvalues, so the
leading ``k`` components always summed to 1.0 whatever share of the variance they really carried (the top
three of twelve components "explained" 100%). It is now each eigenvalue divided by the sum of all eigenvalues
of the covariance matrix, which is what an eigenvalue fraction means.

The oracle is the singular value decomposition of the standardised genotype matrix: the squared singular
values are proportional to the covariance eigenvalues, so their normalised values are the true fractions.
"""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.structural_variants.population.sv_population import sv_population_structure


def _genotypes(n_svs: int, n_samples: int, seed: int) -> np.ndarray:
    """Genotypes (0/1/2) with a population structure: two groups of samples with different allele frequencies."""
    rng = np.random.default_rng(seed)
    group = np.arange(n_samples) % 2
    freq = rng.uniform(0.05, 0.6, size=(n_svs, 2))
    alleles = rng.random((2, n_svs, n_samples)) < freq[:, group][None, :, :]
    return alleles.sum(axis=0).astype(int)


def _true_fractions(genotypes: np.ndarray) -> np.ndarray:
    data = genotypes.T.astype(float)  # samples x SVs
    std = data.std(axis=0)
    std[std == 0] = 1.0
    singular_values = np.linalg.svd((data - data.mean(axis=0)) / std, compute_uv=False)
    variance = singular_values**2
    return variance / variance.sum()


# (n_svs, n_samples): more SVs than samples uses the sample-by-sample covariance, fewer SVs the SV-by-SV one.
SHAPES = [(60, 12), (12, 60), (25, 25)]


class TestVarianceExplained:
    @pytest.mark.parametrize(("n_svs", "n_samples"), SHAPES)
    @pytest.mark.parametrize("n_components", [1, 3, 5])
    def test_fractions_are_shares_of_the_total_variance(self, n_svs: int, n_samples: int, n_components: int) -> None:
        genotypes = _genotypes(n_svs, n_samples, seed=n_svs + n_samples)

        result = sv_population_structure(genotypes, n_components=n_components)

        expected = _true_fractions(genotypes)[:n_components]
        assert result["variance_explained"] == pytest.approx(expected, rel=1e-8, abs=1e-12)
        assert len(result["variance_explained"]) == n_components
        # The leading components do not carry everything.
        assert 0.0 < sum(result["variance_explained"]) < 1.0

    @pytest.mark.parametrize(("n_svs", "n_samples"), SHAPES)
    def test_all_components_together_explain_everything(self, n_svs: int, n_samples: int) -> None:
        genotypes = _genotypes(n_svs, n_samples, seed=7)
        result = sv_population_structure(genotypes, n_components=min(n_svs, n_samples))
        assert sum(result["variance_explained"]) == pytest.approx(1.0)

    def test_the_fractions_do_not_depend_on_how_many_components_are_asked_for(self) -> None:
        genotypes = _genotypes(40, 15, seed=2)
        few = sv_population_structure(genotypes, n_components=2)["variance_explained"]
        more = sv_population_structure(genotypes, n_components=6)["variance_explained"]
        assert more[:2] == pytest.approx(few)
        assert sum(more) > sum(few)

    def test_fractions_are_non_increasing_and_non_negative(self) -> None:
        values = sv_population_structure(_genotypes(50, 20, seed=4), n_components=8)["variance_explained"]
        assert all(v >= 0.0 for v in values)
        assert values == sorted(values, reverse=True)

    def test_eigenvalues_and_fractions_agree(self) -> None:
        # Fractions are eigenvalue ratios, so two components keep their ratio however many are returned.
        result = sv_population_structure(_genotypes(30, 30, seed=9), n_components=4)
        eigenvalues, fractions = result["eigenvalues"], result["variance_explained"]
        assert fractions[0] / fractions[1] == pytest.approx(eigenvalues[0] / eigenvalues[1])

    def test_data_without_variance_explains_nothing(self) -> None:
        result = sv_population_structure(np.ones((10, 8), dtype=int), n_components=3)
        assert result["variance_explained"] == [0.0, 0.0, 0.0]

    def test_more_components_than_the_data_has_are_limited(self) -> None:
        genotypes = _genotypes(6, 9, seed=5)
        result = sv_population_structure(genotypes, n_components=50)
        assert len(result["variance_explained"]) == 6
        assert sum(result["variance_explained"]) == pytest.approx(1.0)
