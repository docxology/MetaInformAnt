"""Statistical audits of the fine-mapping implementations in credible_sets.

These tests exercise the real fine-mapping statistics: approximate Bayes
factors, posterior inclusion probabilities, credible-set construction,
SuSiE regression, conditional analysis, and colocalization posteriors.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from metainformant.gwas.finemapping.credible_sets import (
    annotate_credible_set,
    colocalization,
    compute_bayes_factors,
    compute_credible_set,
    conditional_analysis,
    susie_regression,
)


# ---------------------------------------------------------------------------
# Bayes factors and PIP normalization
# ---------------------------------------------------------------------------


class TestPipNormalization:
    """ABF/PIP audits: posteriors must be a proper probability vector."""

    def test_bayes_factors_positive_and_monotone_in_z(self) -> None:
        """ABF grows with |z| and stays positive."""
        z = [0.0, 1.0, 3.0, 5.0, 8.0]
        bf = compute_bayes_factors(z)

        assert len(bf) == 5
        assert all(b > 0 for b in bf)
        assert bf == sorted(bf)
        # Wakefield ABF: strong z (8) must beat the null variant by a wide margin.
        assert bf[-1] > 2.0 * bf[0]

    def test_pips_sum_to_one(self) -> None:
        """PIP vector must sum to 1 for uniform-prior credible sets."""
        z = [2.0, 5.0, 1.0, 0.5, 3.0]
        result = compute_credible_set(z)

        assert result["status"] == "success"
        assert len(result["pips"]) == 5
        assert math.isclose(sum(result["pips"]), 1.0, rel_tol=1e-9, abs_tol=1e-12)
        assert all(0.0 <= p <= 1.0 for p in result["pips"])

    def test_pips_sum_to_one_with_ld_adjustment(self) -> None:
        """LD down-weighting must not break the PIP normalization."""
        z = [8.0, 7.5]
        ld = [[1.0, 1.0], [1.0, 1.0]]  # Perfect proxies
        result = compute_credible_set(z, ld_matrix=ld)

        assert math.isclose(sum(result["pips"]), 1.0, rel_tol=1e-9, abs_tol=1e-12)

    def test_ld_adjustment_downweights_proxy(self) -> None:
        """A variant in perfect LD with a stronger signal must lose PIP mass."""
        z = [8.0, 7.5]
        ld = [[1.0, 1.0], [1.0, 1.0]]
        result = compute_credible_set(z, ld_matrix=ld)

        assert result["pips"][0] > 0.9
        assert result["pips"][1] < 0.05
        assert result["pips"][1] < 0.1 * result["pips"][0]


# ---------------------------------------------------------------------------
# Credible-set construction
# ---------------------------------------------------------------------------


class TestCredibleSetConstruction:
    """Coverage-honoring credible set audits."""

    def test_single_variant_gets_full_pip_and_own_set(self) -> None:
        """A lone variant is its own credible set with PIP exactly 1."""
        result = compute_credible_set([8.0])

        assert result["status"] == "success"
        assert result["pips"][0] == pytest.approx(1.0, abs=1e-12)
        assert result["snps_in_set"] == [0]
        assert result["n_snps"] == 1
        assert result["coverage_achieved"] == pytest.approx(1.0, abs=1e-12)

    def test_coverage_is_honored(self) -> None:
        """The smallest set that reaches the requested coverage is returned."""
        z = [15.0, 1.0, 1.0, 1.0, 1.0]
        result = compute_credible_set(z, coverage=0.95)

        assert result["status"] == "success"
        assert result["coverage_achieved"] >= 0.95
        assert 0 < result["n_snps"] < 5
        # The dominant variant must anchor the set.
        assert result["snps_in_set"][0] == 0

    def test_two_signals_need_larger_set_than_one(self) -> None:
        """A second comparable signal forces more variants into the set."""
        single = compute_credible_set([15.0], coverage=0.95)
        two = compute_credible_set([12.0, 10.0, 0.5, 0.5], coverage=0.95)

        assert single["n_snps"] == 1
        assert two["n_snps"] >= 2
        assert len(two["snps_in_set"]) == len(set(two["snps_in_set"]))

    def test_invalid_coverage_rejected(self) -> None:
        """Coverage outside (0, 1] must fail loudly."""
        result = compute_credible_set([1.0, 2.0], coverage=1.5)
        assert result["status"] == "error"

        result = compute_credible_set([1.0, 2.0], coverage=0.0)
        assert result["status"] == "error"


# ---------------------------------------------------------------------------
# Conditional analysis: removes its own signal
# ---------------------------------------------------------------------------


class TestConditionalAnalysisSignalRemoval:
    """Stepwise conditioning must not re-report the signal it conditioned."""

    def test_two_independent_signals(self) -> None:
        """Two block-independent signals are both found, once each."""
        z = [8.0, 0.5, 0.3, 0.1, 0.2, 6.0, 0.4, 0.1]
        ld = np.eye(8).tolist()

        signals = conditional_analysis(z, ld, max_signals=10)

        assert len(signals) == 2
        assert [s["index"] for s in signals] == [0, 5]
        assert [s["signal_number"] for s in signals] == [1, 2]
        # The conditioned lead is removed: its conditional z is zeroed.
        assert signals[0]["z_conditional"] == pytest.approx(8.0)
        assert signals[0]["z_score"] == pytest.approx(8.0)

    def test_single_signal_not_rereported(self) -> None:
        """After conditioning, the lead cannot be found a second time."""
        z = [8.0, 0.5, 0.3, 0.1, 0.2, 0.4, 0.4, 0.1]
        ld = np.eye(8).tolist()

        signals = conditional_analysis(z, ld, max_signals=10)

        assert len(signals) == 1
        assert signals[0]["index"] == 0

    def test_proxy_in_perfect_ld_is_absorbed(self) -> None:
        """A near-perfect proxy of the lead must not surface as signal 2."""
        z = [8.0, 7.9, 0.2, 0.1]
        ld = [
            [1.0, 0.999, 0.0, 0.0],
            [0.999, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]

        signals = conditional_analysis(z, ld, max_signals=10)

        assert len(signals) == 1
        assert signals[0]["index"] == 0

    def test_max_signals_cap(self) -> None:
        """max_signals bounds the number of reported signals."""
        z = [8.0, 6.0, 5.5, 0.1]
        ld = np.eye(4).tolist()

        signals = conditional_analysis(z, ld, max_signals=2)

        assert len(signals) == 2
        assert [s["index"] for s in signals] == [0, 1]


# ---------------------------------------------------------------------------
# SuSiE regression audits
# ---------------------------------------------------------------------------


class TestSusieRegressionAudits:
    """Single-effect SuSiE must produce normalized alphas and a sane PIP."""

    def test_single_causal_variant_identified(self) -> None:
        """The simulated causal variant wins the PIP ranking."""
        rng = np.random.default_rng(9)
        n, p = 60, 12
        X = rng.binomial(2, 0.25, size=(n, p)).astype(float)
        X = X - X.mean(axis=0)
        y = 1.2 * X[:, 4] + 0.15 * rng.standard_normal(n)

        result = susie_regression(X, y, L=1, max_iter=100, tol=1e-4)

        assert result["status"] == "success"
        alpha = result["alpha"][0]
        # Posterior over a single effect is a proper distribution.
        assert math.isclose(sum(alpha), 1.0, rel_tol=1e-9, abs_tol=1e-12)
        pip = result["pip"]
        assert len(pip) == p
        assert all(0.0 <= v <= 1.0 for v in pip)
        assert int(np.argmax(pip)) == 4
        assert pip[4] > 0.5
        cs = result["credible_sets"][0]
        assert 4 in cs["indices"]
        assert cs["coverage_achieved"] >= 0.95

    def test_return_contract(self) -> None:
        """The result contract carries converged/elbo/sigma2 of sane types."""
        rng = np.random.default_rng(10)
        n, p = 30, 6
        X = rng.standard_normal((n, p))
        y = X[:, 2] + 0.5 * rng.standard_normal(n)

        result = susie_regression(X, y, L=1, max_iter=20)

        assert result["status"] == "success"
        assert isinstance(result["converged"], bool)
        assert math.isfinite(result["elbo"])
        assert result["sigma2"] > 0
        assert len(result["alpha"]) == 1
        assert len(result["mu"]) == 1

    def test_insufficient_samples(self) -> None:
        """Fewer than 3 samples errors instead of fabricating inference."""
        result = susie_regression([[1.0], [2.0]], [0.0, 1.0])
        assert result["status"] == "error"


# ---------------------------------------------------------------------------
# Colocalization posteriors
# ---------------------------------------------------------------------------


class TestColocalizationPosteriorAudits:
    """Coloc posteriors must be normalized and hypothesis-consistent."""

    def test_posteriors_sum_to_one(self) -> None:
        """The five hypothesis probabilities are a proper distribution."""
        z1 = [2.0, 1.0, 0.5, 0.2, 0.1, 0.3, 0.2, 0.1, 0.4, 0.2]
        z2 = [1.5, 0.8, 0.6, 0.1, 0.2, 0.3, 0.2, 0.1, 0.4, 0.2]

        result = colocalization(z1, z2)

        assert result["status"] == "success"
        total = sum(result[f"PP_H{i}"] for i in range(5))
        assert math.isclose(total, 1.0, rel_tol=1e-9, abs_tol=1e-12)

    def test_shared_strong_signal_favors_h4(self) -> None:
        """Both traits sharing one overwhelming variant must drive PP_H4."""
        strong = [20.0] + [1.0] * 9
        z1 = strong[:]
        z2 = strong[:]

        result = colocalization(z1, z2)

        assert result["status"] == "success"
        assert result["PP_H4"] > 0.9
        assert result["PP_H4"] > result["PP_H3"]
        assert "H4" in result["summary"]

    def test_trait2_only_signal_favors_h2_over_shared(self) -> None:
        """A signal present only in trait 2 must not be called colocalized."""
        z1 = [1.0] * 10
        z2 = [20.0] + [1.0] * 9

        result = colocalization(z1, z2)

        assert result["status"] == "success"
        assert result["PP_H2"] > result["PP_H1"]
        assert result["PP_H2"] > result["PP_H4"]

    def test_input_validation(self) -> None:
        """Empty, mismatched, or prior-invalid inputs return error status."""
        assert colocalization([], [1.0])["status"] == "error"
        assert colocalization([1.0, 2.0], [1.0])["status"] == "error"
        assert (
            colocalization([1.0, 2.0], [1.0, 2.0], prior_p12=0.0)["status"] == "error"
        )


# ---------------------------------------------------------------------------
# Annotation overlay
# ---------------------------------------------------------------------------


class TestAnnotateCredibleSet:
    """Annotation overlay must not fabricate enrichment."""

    def test_no_annotations_no_fabricated_enrichment(self) -> None:
        """Without annotations the enrichment fields stay empty."""
        cs = compute_credible_set([15.0, 1.0, 1.0])
        result = annotate_credible_set(cs, annotations=None)

        assert result["status"] == "success"
        assert result["n_annotated"] == 0
        assert result["enrichment"] == {}

    def test_annotated_enrichment_pinned(self) -> None:
        """With annotations, PIP mass concentrates in the tagged category."""
        cs = compute_credible_set([15.0, 1.0, 1.0, 1.0])
        result = annotate_credible_set(
            cs, annotations={"coding": [0], "intronic": [1, 2, 3]}
        )

        assert result["status"] == "success"
        assert result["n_annotated"] == 1
        coding = result["enrichment"]["coding"]
        # The coding variant holds ~95% of the set's PIP but only 25% of variants.
        assert coding["enrichment_fold"] > 2.0
        assert coding["pip_in_category"] == pytest.approx(cs["pips"][0])
