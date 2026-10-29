"""Round-4 regression tests for the phenotype package.

Each test pins a result that was wrong before its fix, against an oracle that is independent of the code under
test (statsmodels OLS, ``scipy.stats``, a closed form, or a restatement of the definition):

- ``analysis.statistical``: a constant predictor was reported as a significant slope (rank-deficient fit
  attributing the intercept to the slope, p = 0.013); integer group labels made every pairwise t-test disappear
  (``str(1) != 1``); a group of one value produced a NaN "result"; column names with spaces crashed the formula.
- ``gwas_integration.phewas``: covariates were removed by sequential one-covariate regressions (not the
  projection onto their span) and their degrees of freedom were ignored; the heritability screen documented a
  one-sided p-value for h2 > 0 but returned the two-sided p of the slope (0.0007 for a slope of the wrong sign
  and h2 = 0); undefined statistics were reported as 0 with p = 1.
- ``integration.cross_omic``: ``covariates`` and ``method`` were ignored; the "GxE interaction" p-value was
  that of the product's correlation with the phenotype, which any main effect makes significant.
- ``analysis.life_course``: ``confidence`` was a constant 0.7 and the confidence intervals were never produced
  (the loop looked for "predicted" in keys that never contain it); window counts double counted events on a
  boundary; ``trait_categories`` was ignored; explicit windows were echoed back unaggregated; unknown outcome
  measures were skipped silently.
- ``behavior.sequence``: the Markov stationarity statistic assumed both halves hold equally many transitions
  from every state.
- smaller ones: pipeline analyzer failures were hidden inside the results of a "successful" run; the regression
  plot jittered with numpy's global random state; ``dominant_frequency`` returned 0 Hz for any recording with a
  DC offset.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
import statsmodels.api as sm
from scipy import stats

from metainformant.phenotype.analysis import life_course, statistical
from metainformant.phenotype.behavior.ethogram import Ethogram
from metainformant.phenotype.behavior.sequence import BehaviorEvent, BehaviorSequence
from metainformant.phenotype.gwas_integration import phewas
from metainformant.phenotype.integration import cross_omic
from metainformant.phenotype.sonic.signal import AcousticSignal
from metainformant.phenotype.workflow.pipeline import PhenotypePipeline, PipelineConfig

matplotlib.use("Agg")

# ---------------------------------------------------------------------------
# analysis.statistical
# ---------------------------------------------------------------------------


class TestLinearRegressionGuards:
    def test_constant_predictor_is_an_error_not_a_significant_slope(self) -> None:
        df = pd.DataFrame({"x": [2.0] * 5, "y": [1.0, 2.0, 3.0, 4.0, 5.0]})
        result = statistical.perform_linear_regression(df, "x", "y")
        assert "error" in result
        assert "slope" not in result and "significant" not in result

    def test_non_numeric_column_is_an_error(self) -> None:
        df = pd.DataFrame({"x": list("abcde"), "y": [1.0, 2.0, 3.0, 4.0, 6.0]})
        assert "error" in statistical.perform_linear_regression(df, "x", "y")

    def test_column_names_with_spaces_are_fitted(self) -> None:
        df = pd.DataFrame(
            {"body length": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0], "wing width": [2.1, 3.9, 6.2, 7.8, 10.3, 11.9]}
        )
        result = statistical.perform_linear_regression(df, "body length", "wing width")
        oracle = stats.linregress(df["body length"], df["wing width"])
        assert result["slope"] == pytest.approx(oracle.slope, rel=1e-9)
        assert result["intercept"] == pytest.approx(oracle.intercept, rel=1e-9)
        assert result["p_value"] == pytest.approx(oracle.pvalue, rel=1e-6)
        assert result["r_squared"] == pytest.approx(oracle.rvalue**2, rel=1e-9)


class TestPairwiseTTests:
    def test_integer_group_labels_are_compared_as_they_are(self) -> None:
        df = pd.DataFrame({"g": [1, 1, 1, 2, 2, 2, 3, 3, 3], "v": [1.0, 2.0, 3.0, 5.0, 6.0, 7.5, 9.0, 10.0, 11.5]})
        results = statistical.get_comprehensive_pairwise_ttests(df, "v", "g")
        assert len(results) == 3  # was 0: every pair was dropped
        by_pair = {(r["group1"], r["group2"]): r for r in results}
        oracle = stats.ttest_ind(df[df.g == 1].v, df[df.g == 3].v, equal_var=False)
        assert by_pair[(1, 3)]["p_value"] == pytest.approx(oracle.pvalue, rel=1e-9)
        assert [r["p_value"] for r in results] == sorted(r["p_value"] for r in results)

    def test_a_group_of_one_value_is_an_error_not_a_nan_result(self) -> None:
        df = pd.DataFrame({"g": ["a", "b", "b"], "v": [1.0, 2.0, 3.0]})
        result = statistical.perform_ttest(df, "v", "g", "a", "b")
        assert "error" in result
        assert statistical.get_comprehensive_pairwise_ttests(df, "v", "g") == []

    def test_every_term_reports_the_sums_of_squares_it_used(self) -> None:
        df = pd.DataFrame(
            {
                "y": [1.0, 2.0, 3.0, 6.0, 7.0, 8.0, 2.0, 3.0, 4.0, 7.0, 8.0, 9.0],
                "f1": ["a"] * 6 + ["b"] * 6,
                "f2": (["x"] * 3 + ["y"] * 3) * 2,
            }
        )
        result = statistical.perform_multifactor_anova(df, "y ~ C(f1) * C(f2)")
        assert {row["anova_type"] for row in result.values()} == {2}


# ---------------------------------------------------------------------------
# gwas_integration.phewas
# ---------------------------------------------------------------------------


def _correlated_covariate_data(seed: int = 5, n: int = 80) -> Dict[str, Any]:
    rng = np.random.default_rng(seed)
    c1 = rng.normal(size=n)
    c2 = 0.9 * c1 + 0.2 * rng.normal(size=n)  # strongly correlated with c1
    g = rng.integers(0, 3, n).astype(float)
    y = 2.0 * c1 - 1.5 * c2 + 0.3 * g + rng.normal(size=n)
    return {"c1": c1, "c2": c2, "g": g, "y": y}


class TestPhewasCovariates:
    def test_genotype_effect_matches_ols_with_all_covariates_in_the_model(self) -> None:
        data = _correlated_covariate_data()
        covariates = np.column_stack([data["c1"], data["c2"]]).tolist()

        result = phewas.run_phewas(data["g"].tolist(), {"y": data["y"].tolist()}, covariates=covariates)[0]

        design = sm.add_constant(np.column_stack([data["c1"], data["c2"], data["g"]]))
        fit = sm.OLS(data["y"], design).fit()
        assert result["beta"] == pytest.approx(fit.params[3], rel=1e-9)
        assert result["se"] == pytest.approx(fit.bse[3], rel=1e-9)
        assert result["p_value"] == pytest.approx(fit.pvalues[3], rel=1e-8)

    def test_degrees_of_freedom_are_reduced_by_the_covariates(self) -> None:
        # same data, same effect: adjusting for a covariate that explains nothing still costs a degree of freedom
        rng = np.random.default_rng(1)
        n = 12
        g = rng.integers(0, 3, n).astype(float)
        y = 0.6 * g + rng.normal(size=n)
        noise_covariates = rng.normal(size=(n, 4))

        adjusted = phewas.run_phewas(g.tolist(), {"y": y.tolist()}, covariates=noise_covariates.tolist())[0]

        fit = sm.OLS(y, sm.add_constant(np.column_stack([noise_covariates, g]))).fit()
        assert fit.df_resid == n - 6
        assert adjusted["p_value"] == pytest.approx(fit.pvalues[-1], rel=1e-8)

    def test_residualise_is_the_projection_onto_the_covariate_span(self) -> None:
        data = _correlated_covariate_data(seed=2, n=50)
        covariates = np.column_stack([data["c1"], data["c2"]])
        residuals = np.array(phewas._residualise(data["y"].tolist(), covariates.tolist()))
        design = np.column_stack([np.ones(50), covariates])
        expected = data["y"] - design @ np.linalg.lstsq(design, data["y"], rcond=None)[0]
        np.testing.assert_allclose(residuals, expected, atol=1e-10)
        # exactly orthogonal to every covariate (the sequential approximation left part of c2's effect)
        assert abs(residuals @ data["c1"]) < 1e-8 and abs(residuals @ data["c2"]) < 1e-8

    def test_more_covariates_than_the_samples_allow_is_an_error(self) -> None:
        with pytest.raises(ValueError, match="Cannot adjust"):
            phewas._residualise([1.0, 2.0, 3.0], [[1.0, 0.0, 1.0], [0.0, 1.0, 1.0], [1.0, 1.0, 0.0]])

    def test_covariate_rows_must_match_the_samples(self) -> None:
        with pytest.raises(ValueError, match="covariates has 2 rows"):
            phewas.run_phewas([0, 1, 2, 1], {"y": [1.0, 2.0, 3.0, 4.0]}, covariates=[[1.0], [2.0]])

    def test_a_sample_with_a_missing_covariate_is_dropped(self) -> None:
        g = [0, 1, 2, 0, 1, 2, 0, 1, 2, 1]
        y = [0.1, 1.2, 2.1, 0.3, 0.9, 2.2, -0.1, 1.1, 1.8, 1.0]
        covariates = [[float(i)] for i in range(10)]
        covariates[4] = [float("nan")]
        result = phewas.run_phewas(g, {"y": y}, covariates=covariates)[0]
        assert result["n_samples"] == 9

    def test_a_genotype_without_variation_is_skipped_not_reported_as_beta_0_p_1(self) -> None:
        result = phewas.run_phewas([1] * 8, {"y": [float(i) for i in range(8)]})
        assert result == []


class TestHeritabilityScreen:
    @staticmethod
    def _kinship() -> List[List[float]]:
        return [[1.0 if i == j else (0.5 if (i < 5 and j < 5) else 0.0) for j in range(10)] for i in range(10)]

    def test_p_value_is_one_sided_for_a_positive_h2(self) -> None:
        # relatives are *less* alike than strangers: the squared-difference slope is positive, h2 is clamped to 0
        # and the evidence for h2 > 0 is nil. The two-sided p of the slope called this significant (7e-4).
        pheno = {"p": [0.0, 10.0, -10.0, 20.0, -20.0, 1.0, 1.1, 0.9, 1.05, 0.95]}
        result = phewas.phenotype_heritability_screen(pheno, self._kinship())[0]
        assert result["h2"] == 0.0
        assert result["p_value"] > 0.9

    def test_p_value_is_one_sided_for_a_genuinely_heritable_trait(self) -> None:
        rng = np.random.default_rng(3)
        n = 30
        kin = np.eye(n)
        for block in range(0, n, 3):  # families of three with kinship 0.5
            for i in range(block, block + 3):
                for j in range(block, block + 3):
                    if i != j:
                        kin[i, j] = 0.5
        family = np.repeat(rng.normal(scale=3.0, size=n // 3), 3)
        trait = family + rng.normal(scale=0.5, size=n)

        result = phewas.phenotype_heritability_screen({"t": trait.tolist()}, kin)[0]

        # re-derive the one-sided p from the pairwise regression with scipy as the oracle
        diffs, kins = [], []
        for i in range(n):
            for j in range(i + 1, n):
                diffs.append((trait[i] - trait[j]) ** 2)
                kins.append(kin[i, j])
        fit = stats.linregress(kins, diffs)
        assert fit.slope < 0
        assert result["p_value"] == pytest.approx(fit.pvalue / 2.0, rel=1e-6)
        assert result["h2"] == pytest.approx(min(1.0, max(0.0, -fit.slope / (2.0 * np.var(trait, ddof=1)))), rel=1e-9)

    def test_undefined_heritability_is_nan_not_zero_with_p_one(self) -> None:
        flat = phewas.phenotype_heritability_screen({"flat": [5.0] * 6}, self._kinship()[:6])
        assert all(math.isnan(flat[0][key]) for key in ("h2", "se", "p_value"))
        # identical kinship for every pair leaves no slope to estimate
        identity = phewas.phenotype_heritability_screen({"y": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]}, np.eye(6))
        assert math.isnan(identity[0]["p_value"])

    def test_undefined_results_sort_after_defined_ones(self) -> None:
        pheno = {"flat": [5.0] * 10, "p": [0.0, 10.0, -10.0, 20.0, -20.0, 1.0, 1.1, 0.9, 1.05, 0.95]}
        results = phewas.phenotype_heritability_screen(pheno, self._kinship())
        assert [r["phenotype"] for r in results] == ["p", "flat"]


class TestPhenotypeCorrelationMatrixUndefinedPairs:
    def test_constant_or_sparse_phenotypes_have_nan_correlations(self) -> None:
        phenotypes = {
            "a": [1.0, 2.0, 3.0, 4.0, 5.0],
            "flat": [7.0, 7.0, 7.0, 7.0, 7.0],
            "sparse": [1.0, None, None, None, 2.0],
        }
        out = phewas.phenotype_correlation_matrix(phenotypes)
        matrix, pvalues = out["correlation_matrix"], out["p_value_matrix"]
        assert math.isnan(matrix[0][1]) and math.isnan(pvalues[0][1])  # constant
        assert math.isnan(matrix[0][2]) and math.isnan(pvalues[0][2])  # two complete pairs
        assert matrix[0][0] == 1.0

    def test_defined_pairs_are_unchanged(self) -> None:
        a = [1.0, 2.0, 3.0, 4.0, 5.0, 6.5]
        b = [2.0, 1.0, 4.5, 3.0, 6.0, 5.5]
        out = phewas.phenotype_correlation_matrix({"a": a, "b": b})
        oracle = stats.pearsonr(a, b)
        assert out["correlation_matrix"][0][1] == pytest.approx(oracle.statistic, rel=1e-9)
        assert out["p_value_matrix"][0][1] == pytest.approx(oracle.pvalue, rel=1e-6)


# ---------------------------------------------------------------------------
# integration.cross_omic
# ---------------------------------------------------------------------------


class TestGenotypeAssociationCovariates:
    def test_covariates_enter_the_model(self) -> None:
        data = _correlated_covariate_data(seed=8, n=60)
        samples = [f"s{i}" for i in range(60)]
        phenotypes = dict(zip(samples, data["y"].tolist()))

        result = cross_omic.phenotype_genotype_association(
            phenotypes,
            {"v": data["g"].astype(int).tolist()},
            covariates={"c1": data["c1"].tolist(), "c2": data["c2"].tolist()},
        )["associations"]["v"]

        design = sm.add_constant(np.column_stack([data["c1"], data["c2"], data["g"]]))
        fit = sm.OLS(data["y"], design).fit()
        assert result["beta"] == pytest.approx(fit.params[3], rel=1e-9)
        assert result["p_value"] == pytest.approx(fit.pvalues[3], rel=1e-8)
        assert result["r_squared"] == pytest.approx(fit.rsquared, rel=1e-9)
        # and the unadjusted answer is a different one
        unadjusted = cross_omic.phenotype_genotype_association(phenotypes, {"v": data["g"].astype(int).tolist()})
        assert unadjusted["associations"]["v"]["beta"] != pytest.approx(result["beta"], rel=1e-3)

    def test_covariate_lengths_are_checked(self) -> None:
        with pytest.raises(ValueError, match="values for 4 samples"):
            cross_omic.phenotype_genotype_association(
                {f"s{i}": float(i) for i in range(4)}, {"v": [0, 1, 2, 1]}, covariates={"c": [1.0, 2.0]}
            )

    def test_covariates_are_refused_for_the_correlation_method(self) -> None:
        with pytest.raises(ValueError, match="only supported with method='linear_regression'"):
            cross_omic.phenotype_genotype_association(
                {f"s{i}": float(i) for i in range(4)},
                {"v": [0, 1, 2, 1]},
                method="correlation",
                covariates={"c": [1.0, 2.0, 3.0, 4.0]},
            )

    def test_unknown_method_raises_instead_of_running_the_regression(self) -> None:
        with pytest.raises(ValueError, match="Unknown method"):
            cross_omic.phenotype_genotype_association({"a": 1.0, "b": 2.0, "c": 3.0}, {"v": [0, 1, 2]}, method="anova")

    def test_untestable_variants_are_listed_not_reported_as_beta_0_p_1(self) -> None:
        phenotypes = {f"s{i}": float(i) + 0.5 * (i % 3) for i in range(8)}
        result = cross_omic.phenotype_genotype_association(
            phenotypes,
            {"mono": [1] * 8, "short": [0, 1, 2], "ok": [0, 1, 2, 0, 1, 2, 1, 0]},
        )
        assert set(result["associations"]) == {"ok"}
        assert set(result["skipped_variants"]) == {"mono", "short"}
        assert "no variation" in result["skipped_variants"]["mono"]


class TestTraitExpressionCorrelationUndefined:
    def test_constant_expression_is_skipped_not_reported_as_r_0_p_1(self) -> None:
        samples = [f"s{i}" for i in range(6)]
        traits = dict(zip(samples, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]))
        expression = {
            "flat": dict(zip(samples, [3.0] * 6)),
            "varies": dict(zip(samples, [1.0, 3.0, 2.0, 5.0, 4.0, 6.0])),
        }
        result = cross_omic.trait_expression_correlation(traits, expression)
        assert set(result["correlations"]) == {"varies"}
        assert "no variation in the expression" in result["skipped_genes"]["flat"]
        assert result["top_genes"] == ["varies"]

    def test_unknown_method_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown method"):
            cross_omic.trait_expression_correlation({"a": 1.0}, {"g": {"a": 1.0}}, method="kendall")


class TestMultiPhenotypeIntegration:
    @staticmethod
    def _phenotypes() -> Dict[str, Dict[str, float]]:
        rng = np.random.default_rng(4)
        n = 25
        base = rng.normal(size=n)
        columns = {
            "a": base + 0.2 * rng.normal(size=n),
            "b": 0.8 * base + 0.5 * rng.normal(size=n),
            "c": rng.normal(size=n),
        }
        return {name: {f"s{i}": float(v) for i, v in enumerate(values)} for name, values in columns.items()}

    def test_pca_is_the_eigenanalysis_of_the_correlation_matrix(self) -> None:
        phenotypes = self._phenotypes()
        result = cross_omic.multi_phenotype_integration(phenotypes, method="pca")

        eigenvalues = np.sort(np.linalg.eigvalsh(np.array(result["correlation_matrix"])))[::-1]
        np.testing.assert_allclose(result["pca"]["eigenvalues"], eigenvalues, atol=1e-9)
        assert sum(result["pca"]["explained_variance_ratio"]) == pytest.approx(1.0)
        assert eigenvalues.sum() == pytest.approx(3.0)  # trace of a 3 x 3 correlation matrix
        # a loading is the correlation between a phenotype and a component score
        scores = np.array([result["pca"]["scores"][f"s{i}"] for i in range(25)])
        values = np.array([[phenotypes[name][f"s{i}"] for name in "abc"] for i in range(25)])
        for k, name in enumerate("abc"):
            assert result["pca"]["loadings"][name][0] == pytest.approx(np.corrcoef(values[:, k], scores[:, 0])[0, 1])

    def test_the_default_method_adds_no_pca(self) -> None:
        assert "pca" not in cross_omic.multi_phenotype_integration(self._phenotypes())

    def test_unknown_method_raises_instead_of_returning_the_correlations(self) -> None:
        with pytest.raises(ValueError, match="Unknown method"):
            cross_omic.multi_phenotype_integration(self._phenotypes(), method="ica")

    def test_constant_phenotype_has_nan_correlations_and_blocks_pca(self) -> None:
        phenotypes = self._phenotypes()
        phenotypes["flat"] = {f"s{i}": 1.0 for i in range(25)}
        result = cross_omic.multi_phenotype_integration(phenotypes)
        flat = result["phenotype_names"].index("flat")
        assert math.isnan(result["correlation_matrix"][0][flat])
        with pytest.raises(ValueError, match="constant"):
            cross_omic.multi_phenotype_integration(phenotypes, method="pca")


class TestGxEInteraction:
    @staticmethod
    def _data(interaction: float) -> Dict[str, Any]:
        rng = np.random.default_rng(11)
        n = 60
        genotype = rng.integers(0, 3, n)
        environment = rng.normal(size=n)
        outcome = 1.0 * genotype + 2.0 * environment + interaction * genotype * environment
        outcome = outcome + rng.normal(scale=0.5, size=n)
        samples = [f"s{i}" for i in range(n)]
        return {
            "phenotypes": dict(zip(samples, outcome.tolist())),
            "genotypes": {"v": genotype.tolist()},
            "environment": dict(zip(samples, environment.tolist())),
            "g": genotype.astype(float),
            "e": environment,
            "y": outcome,
        }

    def test_main_effects_alone_are_not_an_interaction(self) -> None:
        data = self._data(interaction=0.0)
        result = cross_omic.phenotype_environment_interaction(
            data["phenotypes"], data["genotypes"], data["environment"]
        )
        entry = result["interactions"]["v"]
        # The product still correlates strongly with the phenotype (r = 0.75), which is all the old p-value
        # (4e-12) measured; the interaction coefficient itself is not distinguishable from 0.
        assert entry["gxe_correlation"] > 0.7
        assert entry["p_value"] > 0.05
        assert result["significant_interactions"] == []

    def test_a_real_interaction_is_detected_with_the_ols_statistics(self) -> None:
        data = self._data(interaction=1.2)
        result = cross_omic.phenotype_environment_interaction(
            data["phenotypes"], data["genotypes"], data["environment"]
        )
        entry = result["interactions"]["v"]
        design = sm.add_constant(np.column_stack([data["g"], data["e"], data["g"] * data["e"]]))
        fit = sm.OLS(data["y"], design).fit()
        assert entry["gxe_beta"] == pytest.approx(fit.params[3], rel=1e-9)
        assert entry["gxe_se"] == pytest.approx(fit.bse[3], rel=1e-9)
        assert entry["p_value"] == pytest.approx(fit.pvalues[3], rel=1e-6)
        assert entry["p_value"] < 1e-10
        assert [row["variant"] for row in result["significant_interactions"]] == ["v"]

    def test_untestable_variants_are_listed(self) -> None:
        data = self._data(interaction=0.0)
        data["genotypes"]["short"] = [0, 1]
        result = cross_omic.phenotype_environment_interaction(
            data["phenotypes"], data["genotypes"], data["environment"]
        )
        assert "short" in result["skipped_variants"] and "short" not in result["interactions"]


# ---------------------------------------------------------------------------
# analysis.life_course
# ---------------------------------------------------------------------------


def _sequence(person: str, *events: tuple) -> life_course.EventSequence:
    return life_course.EventSequence(
        person_id=person, events=[life_course.Event(timestamp=t, event_type=kind) for t, kind in events]
    )


class TestPredictLifeCourseOutcomes:
    def test_no_constant_confidence_and_real_intervals_for_every_prediction(self) -> None:
        sequence = _sequence("p1", (0, "job_start"), (400, "job_start"), (800, "job_start"), (1500, "health_event"))

        result = life_course.predict_life_course_outcomes([sequence], prediction_horizon=2.0)

        assert all("confidence" not in prediction for prediction in result["outcome_predictions"].values())
        assert set(result["confidence_intervals"]) == {"p1_career", "p1_health", "p1_social"}

        exposure = max(1.0, 1500 / 365.25)
        for key, (current, label) in {"p1_career": (3, "predicted_additional_jobs"), "p1_health": (1, None)}.items():
            interval = result["confidence_intervals"][key]
            lower, upper = (
                int(q) for q in stats.nbinom.ppf([0.025, 0.975], current + 0.5, exposure / (exposure + 2.0))
            )
            assert interval["lower_bound"] == current + lower
            assert interval["upper_bound"] == current + upper
            assert interval["level"] == 0.95
            assert interval["lower_bound"] <= interval["predicted"] <= interval["upper_bound"]
            if label:
                assert result["outcome_predictions"][key][label] == int(current / exposure * 2.0)

    def test_more_history_gives_a_narrower_relative_interval(self) -> None:
        short = _sequence("short", *[(i * 400, "job_start") for i in range(2)])
        long = _sequence("long", *[(i * 400, "job_start") for i in range(20)])
        result = life_course.predict_life_course_outcomes([short, long], prediction_horizon=2.0)["confidence_intervals"]

        def relative_width(key: str) -> float:
            row = result[key]
            return (row["upper_bound"] - row["lower_bound"]) / max(row["predicted"], 1)

        assert relative_width("long_career") < relative_width("short_career")


class TestAggregateTemporalPhenotypes:
    @staticmethod
    def _sequences() -> List[life_course.EventSequence]:
        return [_sequence("a", (0.0, "x"), (5.0, "x"), (10.0, "x")), _sequence("b", (5.0, "y"))]

    def test_an_event_on_a_window_boundary_is_counted_once(self) -> None:
        result = life_course.aggregate_temporal_phenotypes(self._sequences(), time_window_years=5.0)

        counts = [window["n_events"] for window in result["time_windows"]]
        assert counts == [1, 3]  # [0, 5) and [5, 10]; closed windows gave [2, 4], i.e. 6 events for 4
        assert sum(counts) == result["aggregates"]["total_events"] == 4

    def test_explicit_windows_are_aggregated_and_still_echoed(self) -> None:
        result = life_course.aggregate_temporal_phenotypes(self._sequences(), time_windows=[(0.0, 5.0), (5.0, 10.0)])
        assert result["time_windows"] == [(0.0, 5.0), (5.0, 10.0)]
        assert [w["n_events"] for w in result["window_aggregates"]] == [1, 3]
        assert [w["n_people"] for w in result["window_aggregates"]] == [1, 2]

    def test_all_events_at_one_instant_still_form_a_window(self) -> None:
        result = life_course.aggregate_temporal_phenotypes([_sequence("a", (3.0, "x"), (3.0, "y"))])
        assert [w["n_events"] for w in result["time_windows"]] == [2]

    def test_trait_categories_restrict_the_events_by_domain(self) -> None:
        from datetime import datetime

        from metainformant.life_events.core.events import Event, EventSequence

        sequence = EventSequence(
            person_id="p",
            events=[
                Event("diagnosis", datetime(2010, 1, 1), "health"),
                Event("degree", datetime(2012, 1, 1), "education"),
                Event("promotion", datetime(2014, 1, 1), "occupation"),
            ],
        )
        everything = life_course.aggregate_temporal_phenotypes([sequence], time_window_years=1.0)
        only_health = life_course.aggregate_temporal_phenotypes(
            [sequence], time_window_years=1.0, trait_categories=["health", "occupation"]
        )
        assert everything["aggregates"]["total_events"] == 3
        assert only_health["aggregates"]["total_events"] == 2
        assert sum(w["n_events"] for w in only_health["time_windows"]) == 2
        assert all(set(w["domain_counts"]) <= {"health", "occupation"} for w in only_health["time_windows"])
        none_match = life_course.aggregate_temporal_phenotypes([sequence], trait_categories=["housing"])
        assert none_match["aggregates"]["total_events"] == 0


class TestTrajectoryOutcomeMeasures:
    def test_unknown_measure_raises_instead_of_being_skipped(self) -> None:
        with pytest.raises(ValueError, match="Unknown outcome measure"):
            life_course.analyze_life_course_trajectories([_sequence("p", (0, "job_start"))], ["career_succes"])

    def test_the_documented_scores(self) -> None:
        sequence = _sequence(
            "p",
            *[(i, "job_start") for i in range(7)],
            *[(100 + i, "health_event") for i in range(3)],
            (200, "relationship_start"),
            (201, "relationship_start"),
            (202, "relationship_end"),
        )
        result = life_course.analyze_life_course_trajectories([sequence])["outcome_analysis"]
        assert result["career_success"]["distribution"] == [10]  # min(10, 2 * 7)
        assert result["health_outcomes"]["distribution"] == [7]  # 10 - 3
        assert result["social_stability"]["distribution"] == [1]  # 2 - 1


# ---------------------------------------------------------------------------
# behavior.sequence
# ---------------------------------------------------------------------------


class TestMarkovStationarity:
    @staticmethod
    def _sequence(codes: str) -> BehaviorSequence:
        ethogram = Ethogram({"A": "a", "B": "b"})
        events = [BehaviorEvent(timestamp=float(i), code=c, duration=1.0) for i, c in enumerate(codes)]
        return BehaviorSequence(events, ethogram)

    @staticmethod
    def _oracle(codes: str) -> tuple:
        """Sum over states of scipy's chi-square test of homogeneity of the halves' next-behavior counts."""
        half = len(codes) // 2
        first, second = codes[:half], codes[half:]

        def counts(part: str) -> Dict[str, Dict[str, int]]:
            table: Dict[str, Dict[str, int]] = {}
            for a, b in zip(part, part[1:]):
                table.setdefault(a, {}).setdefault(b, 0)
                table[a][b] += 1
            return table

        c1, c2 = counts(first), counts(second)
        statistic, dof = 0.0, 0
        for state in set(c1) & set(c2):
            targets = sorted(set(c1[state]) | set(c2[state]))
            if len(targets) < 2:
                continue
            table = [[c1[state].get(t, 0) for t in targets], [c2[state].get(t, 0) for t in targets]]
            chi2, _, dof_state, _ = stats.chi2_contingency(table, correction=False)
            statistic += chi2
            dof += dof_state
        return statistic, dof

    @pytest.mark.parametrize(
        "codes",
        [
            "ABABABABAB" + "AAAABAAAABAAAAB",
            "AABABBAABBAABABBBAABABAABBABAB",
            "AAAAABAAAABAABBBABBBABBBBABBBA",
        ],
    )
    def test_matches_the_per_state_homogeneity_tests(self, codes: str) -> None:
        result = self._sequence(codes).markov_stationarity_chi2()
        statistic, dof = self._oracle(codes)
        assert result["sufficient_data"] is True
        assert result["chi2"] == pytest.approx(statistic, rel=1e-9)
        assert result["df"] == dof
        assert result["p_value"] == pytest.approx(float(stats.chi2.sf(statistic, dof)), rel=1e-9)

    def test_halves_with_different_numbers_of_transitions_from_a_state_are_compared_by_proportions(self) -> None:
        # State A is followed by B and C in equal numbers in the first half (4 and 4) and in the second half
        # (1 and 1), so its transition probabilities are the same; the old expected counts of "total / 2" per
        # half, however, made the unequal totals (8 against 2) look like a change.
        first = "ABACABACABACABAC"  # A->B x4, A->C x4 (plus B->A, C->A)
        second = "ABACBBBBBBBBBBBB"  # A->B x1, A->C x1, then B->B
        codes = first + second
        ethogram = Ethogram({"A": "a", "B": "b", "C": "c"})
        events = [BehaviorEvent(timestamp=float(i), code=c, duration=1.0) for i, c in enumerate(codes)]
        result = BehaviorSequence(events, ethogram).markov_stationarity_chi2()

        # Per state, the 2 x K table is (first half, second half) by next behavior:
        #   A -> (B, C): [[4, 4], [1, 1]]   identical proportions, contributes exactly 0
        #   B -> (A, B): [[4, 0], [1, 11]]  the second half stays on B
        #   C -> (A, B): [[3, 0], [0, 1]]
        by_state = {
            "A": stats.chi2_contingency([[4, 4], [1, 1]], correction=False),
            "B": stats.chi2_contingency([[4, 0], [1, 11]], correction=False),
            "C": stats.chi2_contingency([[3, 0], [0, 1]], correction=False),
        }
        assert by_state["A"][0] == pytest.approx(0.0, abs=1e-12)
        assert result["chi2"] == pytest.approx(sum(t[0] for t in by_state.values()), rel=1e-9)
        assert result["df"] == sum(t[2] for t in by_state.values())

    def test_nothing_to_compare_is_reported_as_insufficient_with_a_nan_p_value(self) -> None:
        for codes in ("AB", "ABABABABABAB"):  # too short; and a single outcome per state in both halves
            result = self._sequence(codes).markov_stationarity_chi2()
            assert result["sufficient_data"] is False
            assert result["chi2"] == 0.0 and result["df"] == 0
            assert math.isnan(result["p_value"])


# ---------------------------------------------------------------------------
# workflow, plots, signal
# ---------------------------------------------------------------------------


class _FailingBehaviouralPipeline(PhenotypePipeline):
    def _analyze_behavioral(self, data: Any, params: Dict) -> Dict[str, Any]:
        raise RuntimeError("analyzer exploded")


class TestPipelineAnalyzerFailures:
    def test_a_failing_analyzer_fails_the_run_and_the_other_domains_still_run(self) -> None:
        config = PipelineConfig(
            name="p",
            phenotype_types=["behavioral", "morphological"],
            steps=["load", "validate", "analyze", "summarize"],
        )
        result = _FailingBehaviouralPipeline(config).run([{"a": 1}, {"a": 2}])

        assert result.success is False
        assert any("behavioral" in error and "analyzer exploded" in error for error in result.errors)
        results = result.outputs["analyze"]["results"]
        assert results["behavioral"] == {"error": "analyzer exploded"}
        assert results["morphological"]["n_specimens"] == 2
        assert "summarize" in result.outputs  # later steps were not abandoned

    def test_healthy_runs_stay_successful(self) -> None:
        config = PipelineConfig(name="p", phenotype_types=["morphological"], steps=["load", "analyze"])
        assert PhenotypePipeline(config).run([{"a": 1}]).success is True


class TestRegressionScatterJitter:
    @staticmethod
    def _offsets() -> np.ndarray:
        pytest.importorskip("seaborn")
        from metainformant.phenotype.visualization.plots import plot_regression_scatter

        df = pd.DataFrame(
            {"x": [1, 2, 3, 1, 2, 3, 1, 2, 3, 2], "y": [2.0, 3.5, 5.0, 2.2, 3.1, 5.4, 1.9, 3.8, 4.7, 3.3]}
        )
        figure, ax = plt.subplots()
        try:
            plot_regression_scatter(df, "x", "y", ax=ax)
            scatter = [c for c in ax.collections if len(c.get_offsets()) == len(df)][0]
            return np.array(scatter.get_offsets())
        finally:
            plt.close(figure)

    def test_jitter_does_not_depend_on_the_global_random_state(self) -> None:
        np.random.seed(1)
        first = self._offsets()
        np.random.seed(2)
        second = self._offsets()
        np.testing.assert_array_equal(first, second)
        assert not np.array_equal(first[:, 0], [1, 2, 3, 1, 2, 3, 1, 2, 3, 2])  # it did jitter

    def test_global_random_state_is_left_alone(self) -> None:
        np.random.seed(3)
        before = np.random.get_state()[1].copy()
        self._offsets()
        np.testing.assert_array_equal(np.random.get_state()[1], before)


class TestDominantFrequency:
    def test_a_dc_offset_is_not_a_frequency(self) -> None:
        sample_rate = 8000
        t = np.arange(sample_rate) / sample_rate
        signal = AcousticSignal(waveform=5.0 + np.sin(2 * np.pi * 100.0 * t), sample_rate=sample_rate)
        assert signal.dominant_frequency() == pytest.approx(100.0)  # was 0.0 (the offset's bin)

    def test_pure_tone_and_flat_signal(self) -> None:
        assert AcousticSignal.generate_tone(440.0, 0.5).dominant_frequency() == pytest.approx(440.0)
        assert AcousticSignal(waveform=np.ones(16), sample_rate=8000).dominant_frequency() == 0.0
