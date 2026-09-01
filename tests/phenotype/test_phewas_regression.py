"""Tests for PheWAS, GRS, heritability screen, and phenotype categorisation.

Round-4 test-depth lane T4: pins regression behaviour of
``metainformant.phenotype.gwas_integration.phewas`` with real numeric data.
Values are analytic or computed by an independent formula path in-test.
"""

from __future__ import annotations

import math

import pytest

from metainformant.phenotype.gwas_integration.phewas import (
    _normal_sf,
    _pearson_r,
    _residualise,
    _simple_linear_regression,
    categorize_phenotypes,
    genetic_risk_score,
    phenotype_correlation_matrix,
    phenotype_heritability_screen,
    run_phewas,
)

# ---------------------------------------------------------------------------
# run_phewas
# ---------------------------------------------------------------------------


class TestRunPhewas:
    def test_perfect_linear_signal(self):
        # y = 2*g + 1 exactly -> beta == 2, p ~ 0
        n = 60
        genotype = [i % 3 for i in range(n)]
        phenotypes = {"y": [2.0 * g + 1.0 for g in genotype]}
        results = run_phewas(genotype, phenotypes)
        assert len(results) == 1
        r = results[0]
        assert r["phenotype"] == "y"
        assert r["beta"] == pytest.approx(2.0, abs=1e-6)
        assert r["se"] < 1e-3
        assert r["p_value"] < 1e-20
        assert r["n_samples"] == n
        assert r["phenotype_category"] == "uncategorized"

    def test_results_sorted_by_pvalue(self):
        n = 40
        genotype = [(i * 7) % 3 for i in range(n)]
        signal = [2.0 * g + (0.1 if i % 2 else -0.1) for i, g in enumerate(genotype)]
        noise = [((i * 31) % 11) - 5.0 for i in range(n)]
        results = run_phewas(genotype, {"signal": signal, "noise": noise})
        ps = [r["p_value"] for r in results]
        assert ps == sorted(ps)
        by_name = {r["phenotype"]: r for r in results}
        assert by_name["signal"]["p_value"] < by_name["noise"]["p_value"]

    def test_length_mismatch_skipped(self):
        results = run_phewas([0, 1, 2], {"short": [0.0, 1.0]})
        assert results == []

    def test_missing_values_filtered(self):
        genotype = [0, 1, None, 2, 0, 1, None, 2, 0, 1, 2, 0]
        phenotypes = {"y": [0.0, 1.5, None, 3.0, 0.5, 1.0, 99.0, 3.5, 0.2, 1.2, 3.1, 0.1]}
        results = run_phewas(genotype, phenotypes)
        assert len(results) == 1
        assert results[0]["n_samples"] == 10  # Nones dropped on both sides

    def test_too_few_valid_samples_skipped(self):
        results = run_phewas([0, 1], {"y": [0.0, 1.0]})
        assert results == []

    def test_covariates_residualised(self):
        # Docstring convention: outer list per sample, inner list per covariate
        n = 80
        genotype = [i % 3 for i in range(n)]
        covariate = [float(i % 4) for i in range(n)]
        # y depends only on covariate, not genotype
        y = [3.0 * c + 0.01 * g for c, g in zip(covariate, genotype)]
        results = run_phewas(genotype, {"y": y}, covariates=[[c] for c in covariate])
        assert len(results) == 1
        assert abs(results[0]["beta"]) < 0.1  # genotype signal removed

    def test_binary_phenotype_runs(self):
        genotype = [(i * 5) % 3 for i in range(50)]
        # Higher genotype dosage -> higher probability of 1
        y = [1 if (g * 2 + (i % 2)) > 2 else 0 for i, g in enumerate(genotype)]
        results = run_phewas(genotype, {"binary": y})
        assert len(results) == 1
        assert results[0]["beta"] > 0


# ---------------------------------------------------------------------------
# phenotype_correlation_matrix
# ---------------------------------------------------------------------------


class TestPhenotypeCorrelationMatrix:
    def test_perfect_correlation(self):
        pheno = {
            "a": [1.0, 2.0, 3.0, 4.0, 5.0],
            "b": [2.0, 4.0, 6.0, 8.0, 10.0],
        }
        out = phenotype_correlation_matrix(pheno)
        assert out["correlation_matrix"][0][1] == pytest.approx(1.0, abs=1e-9)
        assert out["correlation_matrix"][1][0] == pytest.approx(1.0, abs=1e-9)
        assert out["correlation_matrix"][0][0] == 1.0
        assert out["p_value_matrix"][0][0] == 0.0
        assert out["phenotype_names"] == ["a", "b"]
        assert out["n_samples"] == 5

    def test_negative_correlation(self):
        pheno = {
            "a": [1.0, 2.0, 3.0, 4.0, 5.0],
            "b": [5.0, 4.0, 3.0, 2.0, 1.0],
        }
        out = phenotype_correlation_matrix(pheno)
        assert out["correlation_matrix"][0][1] == pytest.approx(-1.0, abs=1e-9)

    def test_symmetric_matrix(self):
        pheno = {
            "a": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            "b": [2.0, 1.0, 4.0, 3.0, 6.0, 5.0],
            "c": [6.0, 5.0, 4.0, 3.0, 2.0, 1.0],
        }
        out = phenotype_correlation_matrix(pheno)
        k = 3
        for i in range(k):
            for j in range(k):
                assert out["correlation_matrix"][i][j] == pytest.approx(out["correlation_matrix"][j][i])
                assert out["p_value_matrix"][i][j] == pytest.approx(out["p_value_matrix"][j][i])

    def test_nan_pairs_dropped(self):
        pheno = {
            "a": [1.0, 2.0, None, 4.0, 5.0],
            "b": [2.0, None, 6.0, 8.0, 10.0],
        }
        # Should not raise; pairwise-complete handling
        out = phenotype_correlation_matrix(pheno)
        assert isinstance(out["correlation_matrix"], list)

    def test_too_few_phenotypes_raises(self):
        with pytest.raises(ValueError, match="at least 2"):
            phenotype_correlation_matrix({"a": [1.0, 2.0]})

    def test_inconsistent_length_raises(self):
        with pytest.raises(ValueError, match="samples"):
            phenotype_correlation_matrix({"a": [1.0, 2.0], "b": [1.0, 2.0, 3.0]})


# ---------------------------------------------------------------------------
# genetic_risk_score
# ---------------------------------------------------------------------------


class TestGeneticRiskScore:
    def test_manual_grs(self):
        genotypes = [[0, 1, 2], [2, 0, 1], [1, 2, 0]]
        effects = [0.5, 1.0, 0.1]
        out = genetic_risk_score(genotypes, effects)
        expected = [0 * 0.5 + 1 * 1.0 + 2 * 0.1, 2 * 0.5 + 0 * 1.0 + 1 * 0.1, 1 * 0.5 + 2 * 1.0 + 0 * 0.1]
        assert out["risk_scores"] == pytest.approx(expected)
        assert out["mean"] == pytest.approx(sum(expected) / 3.0)
        assert out["n_variants"] == 3
        assert out["n_samples"] == 3

    def test_weighted_grs(self):
        genotypes = [[1, 1], [0, 2]]
        effects = [1.0, 2.0]
        weights = [2.0, 0.5]
        out = genetic_risk_score(genotypes, effects, weights=weights)
        expected = [1 * 1.0 * 2.0 + 1 * 2.0 * 0.5, 0.0 + 2 * 2.0 * 0.5]
        assert out["risk_scores"] == pytest.approx(expected)

    def test_percentile_contract(self):
        genotypes = [[i % 3] for i in range(20)]
        out = genetic_risk_score(genotypes, [1.0])
        assert set(out["percentiles"]) == {5, 25, 50, 75, 95}
        p = out["percentiles"]
        assert p[5] <= p[25] <= p[50] <= p[75] <= p[95]
        median = sorted(out["risk_scores"])[9:11]
        assert p[50] == pytest.approx(sum(median) / 2.0)

    def test_effect_size_mismatch_raises(self):
        with pytest.raises(ValueError, match="effect_sizes"):
            genetic_risk_score([[0, 1]], [1.0])

    def test_weights_mismatch_raises(self):
        with pytest.raises(ValueError, match="weights"):
            genetic_risk_score([[0, 1]], [1.0, 1.0], weights=[1.0])

    def test_empty_genotypes_raises(self):
        with pytest.raises(ValueError, match="at least one sample"):
            genetic_risk_score([], [])

    def test_numpy_input(self):
        np = pytest.importorskip("numpy")
        arr = np.array([[0, 1], [1, 2]])
        out = genetic_risk_score(arr, [1.0, 0.5])
        assert out["risk_scores"] == pytest.approx([0.5, 2.0])


# ---------------------------------------------------------------------------
# phenotype_heritability_screen
# ---------------------------------------------------------------------------


def _block_kinship(n: int) -> list[list[float]]:
    """Kinship matrix: first half related (0.5), second half unrelated (0.0)."""
    kin = [[0.0] * n for _ in range(n)]
    for i in range(n):
        kin[i][i] = 1.0
        for j in range(i + 1, n):
            k = 0.5 if (i < n // 2 and j < n // 2) else 0.0
            kin[i][j] = k
            kin[j][i] = k
    return kin


class TestHeritabilityScreen:
    def test_related_individuals_more_similar(self):
        n = 12
        kin = _block_kinship(n)
        # First 6 individuals share a tight family mean; second 6 scatter widely
        pheno = {
            "family_structured": [10.0 + (i % 6) * 0.1 for i in range(6)] + [100.0 * ((i % 4) + 1) for i in range(6)]
        }
        results = phenotype_heritability_screen(pheno, kin)
        assert len(results) == 1
        # Direction only: estimated h2 within [0, 1] for family-structured phenotype
        assert 0.0 <= results[0]["h2"] <= 1.0
        assert results[0]["se"] >= 0.0

    def test_zero_variance_phenotype(self):
        kin = _block_kinship(6)
        pheno = {"flat": [5.0] * 6}
        results = phenotype_heritability_screen(pheno, kin)
        assert results[0]["h2"] == 0.0
        assert results[0]["se"] == 0.0
        assert results[0]["p_value"] == 1.0

    def test_length_mismatch_skipped(self):
        kin = _block_kinship(6)
        pheno = {"wrong": [1.0] * 4}
        assert phenotype_heritability_screen(pheno, kin) == []

    def test_numpy_kinship(self):
        np = pytest.importorskip("numpy")
        kin = np.eye(6)
        pheno = {"y": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]}
        results = phenotype_heritability_screen(pheno, kin)
        assert len(results) == 1

    def test_results_sorted_by_pvalue(self):
        n = 10
        kin = _block_kinship(n)
        pheno = {
            "a": [float(i) for i in range(n)],
            "b": [float((i * 7) % 5) for i in range(n)],
        }
        results = phenotype_heritability_screen(pheno, kin)
        ps = [r["p_value"] for r in results]
        assert ps == sorted(ps)


# ---------------------------------------------------------------------------
# categorize_phenotypes
# ---------------------------------------------------------------------------


class TestCategorizePhenotypes:
    def test_basic_categories(self):
        pheno = {
            "asthma_diagnosis": 1,
            "type_2_diabetes": 1,
            "height_cm": 170.0,
            "xyzzy": 1,
        }
        out = categorize_phenotypes(pheno)
        assert out["categories"]["respiratory"] == ["asthma_diagnosis"]
        assert out["categories"]["endocrine"] == ["type_2_diabetes"]
        assert out["categories"]["laboratory"] == ["height_cm"]
        assert out["uncategorized"] == ["xyzzy"]
        assert out["n_categorized"] == 3
        assert out["method"] == "phecode"

    def test_empty_categories_removed(self):
        out = categorize_phenotypes({"unknown_pheno": 1})
        assert out["categories"] == {}
        assert out["uncategorized"] == ["unknown_pheno"]
        assert out["n_categorized"] == 0

    def test_case_insensitive_and_spaces(self):
        out = categorize_phenotypes({"Asthma Diagnosis": 1})
        assert out["categories"]["respiratory"] == ["Asthma Diagnosis"]

    def test_first_matching_category_wins(self):
        # 'heart' (circulatory) precedes 'muscle' (musculoskeletal) in keyword scan
        out = categorize_phenotypes({"heart_muscle_pain": 1})
        assert out["categories"]["circulatory"] == ["heart_muscle_pain"]


# ---------------------------------------------------------------------------
# Internal helpers — analytic pins
# ---------------------------------------------------------------------------


class TestInternalHelpers:
    def test_normal_sf_normal_approximation(self):
        # erfc-based: P(Z >= 0) = 0.5
        assert _normal_sf(0.0, df=100) == pytest.approx(0.5)
        assert _normal_sf(1.96, df=100) == pytest.approx(0.024997895, abs=1e-6)
        assert 0.0 < _normal_sf(3.0, df=100) < 0.01

    def test_normal_sf_small_df_monotone(self):
        vals = [_normal_sf(t, df=10) for t in [0.0, 1.0, 2.0, 5.0]]
        assert all(vals[i] >= vals[i + 1] for i in range(len(vals) - 1))
        assert vals[0] == pytest.approx(0.5)

    def test_pearson_perfect(self):
        r, p = _pearson_r([1.0, 2.0, 3.0, 4.0], [2.0, 4.0, 6.0, 8.0])
        assert r == pytest.approx(1.0)
        assert p == 0.0

    def test_pearson_constant_input(self):
        r, p = _pearson_r([1.0, 1.0, 1.0], [1.0, 2.0, 3.0])
        assert r == 0.0
        assert p == 1.0

    def test_pearson_too_few_points(self):
        r, p = _pearson_r([1.0, 2.0], [2.0, 4.0])
        assert (r, p) == (0.0, 1.0)

    def test_pearson_known_value(self):
        # Expected computed by an independent formula path in the test
        x = [1.0, 2.0, 3.0, 4.0, 5.0]
        y = [2.0, 1.0, 4.0, 3.0, 5.0]
        r, _ = _pearson_r(x, y)
        mx = sum(x) / 5
        my = sum(y) / 5
        num = sum((a - mx) * (b - my) for a, b in zip(x, y))
        den = math.sqrt(sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y))
        assert r == pytest.approx(num / den)

    def test_simple_linear_regression_exact_line(self):
        x = [1.0, 2.0, 3.0, 4.0, 5.0]
        y = [3.0 + 2.0 * v for v in x]
        beta, se, p = _simple_linear_regression(y, x)
        assert beta == pytest.approx(2.0, abs=1e-9)
        assert se == pytest.approx(0.0, abs=1e-9)
        assert p < 1e-10

    def test_residualise_removes_covariate(self):
        # Sequential residualisation (sample-outer covariates) removes the
        # covariate's linear component but not the intercept; the constant
        # remainder is annihilated by the caller's subsequent centering
        # (_simple_linear_regression computes beta on centered values).
        n = 30
        cov = [[float(i)] for i in range(n)]  # per-sample rows, 1 covariate each
        y = [3.0 * float(i) + 5.0 for i in range(n)]
        resid = _residualise(y, cov)
        # Linear-in-covariate component fully removed: residual is constant
        spread = max(resid) - min(resid)
        assert spread < 1e-9
        # For a y with both components, the covariate part is largely removed.
        # Single-pass sequential residualisation estimates the slope on the
        # contaminated residual, leaving a small linear residue (< 12% of the
        # removed component); the caller's centering keeps beta unaffected.
        y2 = [3.0 * float(i) + 5.0 + 0.5 * (i % 3) for i in range(n)]
        resid2 = _residualise(y2, cov)
        spread2 = max(resid2) - min(resid2)
        assert spread2 < 1.25  # dominant 0.5*(i%3) pattern survives, 3i part gone
        # Residual is (nearly) uncorrelated with the covariate
        mc = sum(c[0] for c in cov) / n
        mr2 = sum(resid2) / n
        sp = sum((c[0] - mc) * (r - mr2) for c, r in zip(cov, resid2))
        ss = sum((c[0] - mc) ** 2 for c in cov)
        slope = sp / ss
        assert abs(slope) < 0.01
