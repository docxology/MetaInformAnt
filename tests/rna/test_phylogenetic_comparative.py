"""Tests for phylogenetic comparative methods (PGLS) on validated trees."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from metainformant.rna.analysis import phylogenetic_comparative as pc
from metainformant.rna.analysis.statistics_contract import (
    ProvenanceError,
    TreeInvariantError,
)

# ((A,B):1,(C,D):1): every tip depth 2, sister pairs share depth 1,
# cross-pairs share only the root (0).
SYMMETRIC_NEWICK = "((apis:1.0,bombus:1.0):1.0,(ceratina:1.0,megachile:1.0):1.0);"
# Root is bifurcating; all internal branches zero, so the Brownian VCV is
# the identity (star phylogeny) despite the bifurcating root structure.
STAR_NEWICK = "(apis:1,(bombus:1,(ceratina:1,megachile:1):0):0);"

SYMMETRIC_VCV = np.array(
    [
        [2.0, 1.0, 0.0, 0.0],
        [1.0, 2.0, 0.0, 0.0],
        [0.0, 0.0, 2.0, 1.0],
        [0.0, 0.0, 1.0, 2.0],
    ]
)


class TestTreeValidation:
    def test_undeclared_rootedness_fails_closed(self) -> None:
        with pytest.raises(ProvenanceError, match="rootedness"):
            pc.brownian_vcv("(apis:1,bombus:1);")

    def test_declared_unrooted_raises(self) -> None:
        with pytest.raises(TreeInvariantError, match="declared unrooted"):
            pc.brownian_vcv("(apis:1,bombus:1);", rooted=False)

    def test_malformed_newick_raises(self) -> None:
        with pytest.raises(TreeInvariantError, match="unbalanced"):
            pc.brownian_vcv("((apis:1,bombus:1);", rooted=True)

    def test_missing_branch_length_raises(self) -> None:
        with pytest.raises(TreeInvariantError, match="missing a branch length"):
            pc.brownian_vcv("(apis,bombus);", rooted=True)

    def test_negative_branch_length_raises(self) -> None:
        with pytest.raises(TreeInvariantError, match="non-negative"):
            pc.brownian_vcv("(apis:-1,bombus:1);", rooted=True)

    def test_non_finite_branch_length_raises(self) -> None:
        with pytest.raises(TreeInvariantError, match="finite"):
            pc.brownian_vcv("(apis:nan,bombus:1);", rooted=True)

    def test_duplicate_leaf_labels_raise(self) -> None:
        with pytest.raises(TreeInvariantError, match="duplicate"):
            pc.brownian_vcv("(apis:1,apis:1);", rooted=True)

    def test_non_numeric_branch_length_in_dict_raises(self) -> None:
        tree = {
            "name": "root",
            "children": [
                {"name": "apis", "distance": "1.0"},
                {"name": "bombus", "distance": 1.0},
            ],
        }
        with pytest.raises(TreeInvariantError, match="real number"):
            pc.brownian_vcv(tree, rooted=True)

    def test_dict_tree_missing_branch_length_raises(self) -> None:
        tree = {"name": "root", "children": [{"name": "apis"}, {"name": "bombus", "distance": 1.0}]}
        with pytest.raises(TreeInvariantError, match="missing a branch length"):
            pc.brownian_vcv(tree, rooted=True)

    def test_invalid_tree_type_raises(self) -> None:
        with pytest.raises(TreeInvariantError, match="Newick string or a nested-dict"):
            pc.brownian_vcv(["apis", "bombus"], rooted=True)

    def test_nested_dict_tree_accepted(self) -> None:
        tree = {
            "name": "root",
            "distance": 0.0,
            "children": [
                {"name": "apis", "distance": 2.0},
                {"name": "bombus", "distance": 2.0},
            ],
        }
        vcv = pc.brownian_vcv(tree, rooted=True)
        assert list(vcv.index) == ["apis", "bombus"]
        assert vcv.loc["apis", "apis"] == pytest.approx(2.0)
        assert vcv.loc["apis", "bombus"] == pytest.approx(0.0)


class TestBrownianVcv:
    def test_symmetric_tree_covariance(self) -> None:
        vcv = pc.brownian_vcv(SYMMETRIC_NEWICK, rooted=True)
        assert list(vcv.index) == ["apis", "bombus", "ceratina", "megachile"]
        np.testing.assert_allclose(vcv.to_numpy(), SYMMETRIC_VCV)

    def test_sister_pairs_share_mrca_depth(self) -> None:
        vcv = pc.brownian_vcv(SYMMETRIC_NEWICK, rooted=True)
        assert vcv.loc["apis", "bombus"] == pytest.approx(1.0)
        assert vcv.loc["apis", "ceratina"] == pytest.approx(0.0)
        assert vcv.loc["ceratina", "megachile"] == pytest.approx(1.0)

    def test_root_distance_is_ignored(self) -> None:
        vcv = pc.brownian_vcv("(apis:1,bombus:1):99;", rooted=True)
        assert vcv.loc["apis", "apis"] == pytest.approx(1.0)


class TestLambdaAdjustment:
    def test_lambda_zero_is_star(self) -> None:
        vcv = pc.brownian_vcv(SYMMETRIC_NEWICK, rooted=True)
        adjusted = pc.lambda_adjusted_vcv(vcv, 0.0)
        np.testing.assert_allclose(adjusted.to_numpy(), np.diag(np.diag(vcv.to_numpy())))

    def test_lambda_one_is_brownian(self) -> None:
        vcv = pc.brownian_vcv(SYMMETRIC_NEWICK, rooted=True)
        pd.testing.assert_frame_equal(pc.lambda_adjusted_vcv(vcv, 1.0), vcv)

    def test_lambda_out_of_bounds_raises(self) -> None:
        vcv = pc.brownian_vcv(SYMMETRIC_NEWICK, rooted=True)
        with pytest.raises(ValueError, match="lambda"):
            pc.lambda_adjusted_vcv(vcv, 1.5)
        with pytest.raises(ValueError, match="lambda"):
            pc.lambda_adjusted_vcv(vcv, -0.1)


SEED = 20260917


def _simulated_pgls_data() -> tuple[pd.Series, pd.Series, np.ndarray]:
    """Brownian predictor plus y = 2 + 1.5 x + N(0, sigma^2 V) noise.

    Returns (response, predictor, analytic SYMMETRIC_VCV); sigma = 0.5
    for the residual Brownian process.
    """
    x = pc.simulate_brownian_traits(SYMMETRIC_NEWICK, 1.0, seed=SEED, rooted=True).iloc[:, 0].rename("x")
    noise = pc.simulate_brownian_traits(SYMMETRIC_NEWICK, 0.5, seed=SEED + 1, rooted=True).iloc[:, 0]
    y = pd.Series(2.0 + 1.5 * x + noise, index=x.index)
    return y, x, SYMMETRIC_VCV


def _expected_standard_error(sigma2: float, x_values: np.ndarray, v: np.ndarray) -> float:
    design = np.column_stack([np.ones(x_values.size), x_values])
    information = design.T @ np.linalg.inv(v) @ design
    return float(np.sqrt(sigma2 * np.linalg.inv(information)[1, 1]))


class TestPGLSRecovery:
    def test_brownian_simulation_recovers_coefficients_and_sigma(self) -> None:
        y, x, v = _simulated_pgls_data()
        result = pc.fit_pgls(SYMMETRIC_NEWICK, y, x, rooted=True, lambda_=1.0)
        # True parameters: intercept 2.0, slope 1.5, sigma^2 0.25.
        assert result.coefficients["intercept"] == pytest.approx(2.0, abs=0.6)
        assert result.coefficients["x"] == pytest.approx(1.5, abs=0.6)
        assert result.sigma2 == pytest.approx(0.25, rel=0.5)
        # Recovered slope within 2.5 SE of the truth.
        assert abs(result.coefficients["x"] - 1.5) < 2.5 * result.standard_errors["x"]

    def test_standard_errors_match_analytic_scaling(self) -> None:
        y, x, v = _simulated_pgls_data()
        result = pc.fit_pgls(SYMMETRIC_NEWICK, y, x, rooted=True, lambda_=1.0)
        expected = _expected_standard_error(0.25, x.to_numpy(), v)
        assert result.standard_errors["x"] == pytest.approx(expected, rel=0.5)
        expected_at_hat = _expected_standard_error(float(result.sigma2), x.to_numpy(), v)
        assert result.standard_errors["x"] == pytest.approx(expected_at_hat, rel=1e-8)
        # True slope inside the reported 2-SE interval.
        assert abs(result.coefficients["x"] - 1.5) <= 2.0 * expected

    def test_lambda_near_one_for_brownian_residuals(self) -> None:
        y, x, _ = _simulated_pgls_data()
        result = pc.fit_pgls(SYMMETRIC_NEWICK, y, x, rooted=True)
        assert result.lambda_estimated
        assert 0.5 <= result.lambda_ <= 1.0

    def test_lambda_small_for_independent_residuals(self) -> None:
        # Eight taxa: small panels cannot identify lambda; an 8-taxon
        # pectinate tree separates white noise (lambda ~ 0) from Brownian
        # residual structure (see the companion near-one test).
        pectinate = "(a:1,(b:1,(c:1,(d:1,(e:1,(f:1,(g:1,h:1):1):1):1):1):1):1);"
        species = ["a", "b", "c", "d", "e", "f", "g", "h"]
        rng = np.random.default_rng(SEED)
        x = pd.Series(rng.normal(0.0, 1.0, 8), index=species)
        y = pd.Series(2.0 + 1.5 * x + rng.normal(0.0, 0.5, 8), index=species)
        result = pc.fit_pgls(pectinate, y, x, rooted=True)
        assert result.lambda_ < 0.01

    def test_star_phylogeny_recovers_ols(self) -> None:
        rng = np.random.default_rng(SEED)
        species = ["apis", "bombus", "ceratina", "megachile"]
        x = pd.Series(rng.normal(0.0, 1.0, 4), index=species)
        y = pd.Series(2.0 + 1.5 * x + rng.normal(0.0, 0.5, 4), index=species)
        result = pc.fit_pgls(STAR_NEWICK, y, x, rooted=True, lambda_=0.0)
        # Identity Brownian VCV: PGLS equals OLS exactly.
        design = np.column_stack([np.ones(4), x.to_numpy()])
        ols = np.linalg.lstsq(design, y.to_numpy(), rcond=None)[0]
        np.testing.assert_allclose(result.coefficients.to_numpy(), ols, atol=1e-8)

    def test_fixed_lambda_out_of_bounds_raises(self) -> None:
        y, x, _ = _simulated_pgls_data()
        with pytest.raises(ValueError, match="lambda"):
            pc.fit_pgls(SYMMETRIC_NEWICK, y, x, rooted=True, lambda_=1.2)

    def test_result_fields_and_diagnostics(self) -> None:
        y, x, _ = _simulated_pgls_data()
        result = pc.fit_pgls(SYMMETRIC_NEWICK, y, x, rooted=True, lambda_=0.8)
        assert result.residual_df == 2
        assert result.n_obs == 4
        assert result.aic == pytest.approx(2 * (2 + 1) - 2 * result.log_likelihood)
        assert result.diagnostics["covariance_condition_number"] > 0
        assert result.diagnostics["shapiro_wilk_statistic"] is not None
        assert list(result.covariance.index) == ["apis", "bombus", "ceratina", "megachile"]
        expected_t = result.coefficients / result.standard_errors
        np.testing.assert_allclose(result.t_statistics.to_numpy(), expected_t.to_numpy())


class TestPGLSDataFailureModes:
    def test_species_absent_from_tree_raises(self) -> None:
        y, x, _ = _simulated_pgls_data()
        y_extended = pd.concat([y, pd.Series([0.0], index=["vespa"])])
        x_extended = pd.concat([x, pd.Series([0.0], index=["vespa"])])
        with pytest.raises(ValueError, match="absent from the species tree"):
            pc.fit_pgls(SYMMETRIC_NEWICK, y_extended, x_extended, rooted=True)

    def test_non_finite_response_raises(self) -> None:
        y, x, _ = _simulated_pgls_data()
        y.loc["apis"] = np.nan
        with pytest.raises(ValueError, match="non-finite"):
            pc.fit_pgls(SYMMETRIC_NEWICK, y, x, rooted=True)

    def test_non_finite_predictor_raises(self) -> None:
        y, x, _ = _simulated_pgls_data()
        x.loc["bombus"] = np.inf
        with pytest.raises(ValueError, match="non-finite"):
            pc.fit_pgls(SYMMETRIC_NEWICK, y, x, rooted=True)

    def test_missing_species_in_predictors_raises(self) -> None:
        y, x, _ = _simulated_pgls_data()
        with pytest.raises(ValueError, match="missing species"):
            pc.fit_pgls(SYMMETRIC_NEWICK, y, x.drop(index="megachile"), rooted=True)

    def test_duplicate_species_labels_raise(self) -> None:
        y, x, _ = _simulated_pgls_data()
        y_dup = pd.concat([y, pd.Series([0.1], index=["apis"])])
        x_dup = pd.concat([x, pd.Series([0.1], index=["apis"])])
        with pytest.raises(ValueError, match="duplicate"):
            pc.fit_pgls(SYMMETRIC_NEWICK, y_dup, x_dup, rooted=True)

    def test_non_string_species_labels_raise(self) -> None:
        y, x, _ = _simulated_pgls_data()
        y_int = pd.Series(y.to_numpy(), index=[1, 2, 3, 4])
        x_int = pd.Series(x.to_numpy(), index=[1, 2, 3, 4])
        with pytest.raises(ValueError, match="species name strings"):
            pc.fit_pgls(SYMMETRIC_NEWICK, y_int, x_int, rooted=True)

    def test_no_residual_degrees_of_freedom_raises(self) -> None:
        y = pd.Series([1.0, 2.0], index=["apis", "bombus"])
        x = pd.Series([0.0, 1.0], index=["apis", "bombus"])
        with pytest.raises(ValueError, match="residual degree"):
            pc.fit_pgls("(apis:1,bombus:1);", y, x, rooted=True)

    def test_reserved_intercept_name_raises(self) -> None:
        y, _, _ = _simulated_pgls_data()
        x = pd.DataFrame({"intercept": np.linspace(0.0, 1.0, 4)}, index=y.index)
        with pytest.raises(ValueError, match="intercept"):
            pc.fit_pgls(SYMMETRIC_NEWICK, y, x, rooted=True)

    def test_collinear_predictors_raise(self) -> None:
        y, x, _ = _simulated_pgls_data()
        frame = pd.DataFrame({"a": x.to_numpy(), "b": 2.0 * x.to_numpy()}, index=x.index)
        with pytest.raises(ValueError, match="rank deficient"):
            pc.fit_pgls(SYMMETRIC_NEWICK, y, frame, rooted=True)


class TestTreeUncertainty:
    def _tree_set(self) -> list[str]:
        return [
            SYMMETRIC_NEWICK,
            SYMMETRIC_NEWICK.replace(":1.0):1.0,(ceratina", ":1.2):1.0,(ceratina", 1),
            "((apis:1.0,bombus:1.0):1.5,(ceratina:1.0,megachile:1.0):1.5);",
        ]

    def test_seeded_resampling_is_reproducible(self) -> None:
        y, x, _ = _simulated_pgls_data()
        first = pc.fit_pgls_tree_uncertainty(self._tree_set(), y, x, n_resamples=40, seed=7, rooted=True)
        second = pc.fit_pgls_tree_uncertainty(self._tree_set(), y, x, n_resamples=40, seed=7, rooted=True)
        pd.testing.assert_frame_equal(first["draws"], second["draws"])
        assert first["lambda_mean"] == pytest.approx(second["lambda_mean"])

    def test_different_seed_gives_different_draws(self) -> None:
        y, x, _ = _simulated_pgls_data()
        first = pc.fit_pgls_tree_uncertainty(self._tree_set(), y, x, n_resamples=20, seed=7, rooted=True)
        second = pc.fit_pgls_tree_uncertainty(self._tree_set(), y, x, n_resamples=20, seed=8, rooted=True)
        assert not np.allclose(first["draws"].to_numpy(), second["draws"].to_numpy())

    def test_summary_shapes_and_bounds(self) -> None:
        y, x, _ = _simulated_pgls_data()
        result = pc.fit_pgls_tree_uncertainty(self._tree_set(), y, x, n_resamples=50, seed=7, rooted=True)
        assert result["draws"].shape == (2, 50)
        assert list(result["summary"].index) == ["intercept", "x"]
        assert (result["summary"]["ci_low"] <= result["summary"]["ci_high"]).all()

    def test_uncoverable_species_fails_before_resampling(self) -> None:
        y, x, _ = _simulated_pgls_data()
        y_extra = pd.concat([y, pd.Series([0.0], index=["vespa"])])
        x_extra = pd.concat([x, pd.Series([0.0], index=["vespa"])])
        with pytest.raises(ValueError, match="absent from the species tree"):
            pc.fit_pgls_tree_uncertainty(self._tree_set(), y_extra, x_extra, n_resamples=5, seed=7, rooted=True)

    def test_malformed_tree_in_set_fails_closed(self) -> None:
        y, x, _ = _simulated_pgls_data()
        trees = self._tree_set() + ["((apis:1,bombus:1);"]
        with pytest.raises(TreeInvariantError):
            pc.fit_pgls_tree_uncertainty(trees, y, x, n_resamples=5, seed=7, rooted=True)

    def test_invalid_arguments_raise(self) -> None:
        y, x, _ = _simulated_pgls_data()
        with pytest.raises(ValueError, match="n_resamples"):
            pc.fit_pgls_tree_uncertainty(self._tree_set(), y, x, n_resamples=0, rooted=True)
        with pytest.raises(ValueError, match="at least one tree"):
            pc.fit_pgls_tree_uncertainty([], y, x, rooted=True)
        with pytest.raises(ValueError, match="ci"):
            pc.fit_pgls_tree_uncertainty(self._tree_set(), y, x, ci=(0.9, 0.1), rooted=True)


class TestSimulateBrownianTraits:
    def test_shape_index_and_reproducibility(self) -> None:
        first = pc.simulate_brownian_traits(SYMMETRIC_NEWICK, 1.0, n_traits=3, seed=5, rooted=True)
        second = pc.simulate_brownian_traits(SYMMETRIC_NEWICK, 1.0, n_traits=3, seed=5, rooted=True)
        pd.testing.assert_frame_equal(first, second)
        assert first.shape == (4, 3)
        assert list(first.index) == ["apis", "bombus", "ceratina", "megachile"]
        assert list(first.columns) == ["trait_0", "trait_1", "trait_2"]

    def test_zero_sigma_gives_exact_zeros(self) -> None:
        traits = pc.simulate_brownian_traits(SYMMETRIC_NEWICK, 0.0, n_traits=2, rooted=True)
        np.testing.assert_allclose(traits.to_numpy(), 0.0)

    def test_draws_follow_the_brownian_covariance(self) -> None:
        sigma = 1.0
        traits = pc.simulate_brownian_traits(SYMMETRIC_NEWICK, sigma, n_traits=4000, seed=11, rooted=True)
        sample = np.cov(traits.to_numpy())
        np.testing.assert_allclose(sample, SYMMETRIC_VCV, atol=0.35)

    def test_invalid_sigma_raises(self) -> None:
        with pytest.raises(ValueError, match="sigma"):
            pc.simulate_brownian_traits(SYMMETRIC_NEWICK, -1.0, rooted=True)
        with pytest.raises(ValueError, match="sigma"):
            pc.simulate_brownian_traits(SYMMETRIC_NEWICK, np.nan, rooted=True)

    def test_trait_name_mismatch_raises(self) -> None:
        with pytest.raises(ValueError, match="trait_names"):
            pc.simulate_brownian_traits(SYMMETRIC_NEWICK, 1.0, n_traits=2, trait_names=["a"], rooted=True)
