"""Regression tests for approximate Bayes factors, credible sets and colocalization.

Defect covered: the fine-mapping Bayes factors paired a beta-scale prior
(W = 0.04) with V = 1 on the z-score scale, which made them almost flat. A
z = 7.5 hit among 300 variants then needed a 285-variant 95% credible set and
``eqtl_coloc`` gave PP.H4 = 0.003 for a gene driven by the very same variant.

Independent oracles used here (none of them calls the code under test):

* literal log ABFs computed with 60-digit ``decimal`` arithmetic and cross-checked
  against a base-R transcription of coloc's ``approx.bf.estimates``, ``Var.data``
  and ``Var.data.cc``;
* Wakefield's Bayes factor from its definition, the ratio of the two marginal
  likelihoods N(beta_hat; 0, V + W) / N(beta_hat; 0, V);
* the same ratio by numerical integration over the true effect;
* coloc's five-hypothesis posteriors by brute-force enumeration in ``decimal``,
  and literals from a base-R transcription of ``coloc.bf_bf``.

All inputs are real computations on seeded pseudo-random z-scores; nothing is mocked.
"""

from __future__ import annotations

import itertools
import math
import random
from collections import Counter
from decimal import Decimal, localcontext
from typing import Sequence

import pytest

from metainformant.gwas.finemapping import colocalization as coloc_module, credible_sets

# ---------------------------------------------------------------------------
# Hand-computed oracle values
# ---------------------------------------------------------------------------

#: (z, V, W, log ABF). log ABF = (log(1 - r) + r z^2) / 2 with r = W / (V + W), evaluated
#: with 60-digit decimals and confirmed by a base-R transcription of coloc.
ORACLE_CASES = {
    # beta_hat 0.12, SE 0.03, quantitative trait, W = 0.15^2
    "quantitative_se": (4.0, 0.03**2, 0.15**2, 6.06325942329695128),
    # log OR 0.25, SE 0.05, case-control trait, W = 0.2^2
    "case_control_se": (5.0, 0.05**2, 0.2**2, 10.3480992103248331),
    # N = 5000, MAF 0.2: V = 1 / (2 N f (1 - f)) = 1 / 1600
    "quantitative_n_maf": (6.0, 1.0 / 1600.0, 0.15**2, 15.7080545571914013),
    # N = 1000, standardised genotype: V = 1 / N
    "quantitative_n_standardised": (6.0, 1.0 / 1000.0, 0.15**2, 15.6555423426164327),
    # N = 8000, MAF 0.3, 25% cases: V = 1 / (2 N f (1 - f) s (1 - s)) = 1 / 630
    "case_control_n_maf_cases": (5.5, 1.0 / 630.0, 0.2**2, 12.9148302182805966),
    # z-score scale: V = 1, W = 5.2^2
    "z_scale_5": (5.0, 1.0, 5.2**2, 10.3873922428112418),
    "z_scale_7_5": (7.5, 1.0, 5.2**2, 25.4551525851792875),
}


def _null_z(seed: int, n: int) -> list[float]:
    """Seeded N(0, 1) z-scores."""
    rng = random.Random(seed)
    return [rng.gauss(0.0, 1.0) for _ in range(n)]


def _region_with_signals(seed: int, n: int, signals: dict[int, float]) -> list[float]:
    """Null z-scores with the given ``{index: z}`` signals planted."""
    z = _null_z(seed, n)
    for index, value in signals.items():
        z[index] = value
    return z


def _normal_pdf(x: float, variance: float) -> float:
    return math.exp(-(x * x) / (2.0 * variance)) / math.sqrt(2.0 * math.pi * variance)


def _bayes_factor_from_marginal_likelihoods(beta_hat: float, v: float, w: float) -> float:
    """Wakefield's ABF as the ratio N(beta_hat; 0, V + W) / N(beta_hat; 0, V)."""
    return _normal_pdf(beta_hat, v + w) / _normal_pdf(beta_hat, v)


def _bayes_factor_by_quadrature(beta_hat: float, v: float, w: float, n_steps: int = 40000) -> float:
    """Bayes factor by composite Simpson integration of the likelihood over the true effect."""
    half_width = abs(beta_hat) + 12.0 * math.sqrt(w)
    step = 2.0 * half_width / n_steps

    def integrand(beta: float) -> float:
        return _normal_pdf(beta_hat - beta, v) * _normal_pdf(beta, w)

    total = integrand(-half_width) + integrand(half_width)
    for i in range(1, n_steps):
        total += (4.0 if i % 2 else 2.0) * integrand(-half_width + i * step)
    marginal = total * step / 3.0
    return marginal / _normal_pdf(beta_hat, v)


def _pips_by_enumeration(z: Sequence[float], v: float, w: float) -> list[float]:
    """Uniform-prior PIPs from the closed-form ABF evaluated in 80-digit decimals."""
    with localcontext() as ctx:
        ctx.prec = 80
        dv, dw = Decimal(v), Decimal(w)
        r = dw / (dv + dw)
        bayes_factors = [(((1 - r).ln() + r * Decimal(x) * Decimal(x)) / 2).exp() for x in z]
        total = sum(bayes_factors, Decimal(0))
        return [float(bf / total) for bf in bayes_factors]


def _coloc_posteriors_by_enumeration(
    log_bf_1: Sequence[float],
    log_bf_2: Sequence[float],
    p1: float = 1e-4,
    p2: float = 1e-4,
    p12: float = 1e-5,
    h0_weight: float | None = None,
) -> list[float]:
    """coloc's five hypothesis posteriors by explicit enumeration over variant pairs.

    ``h0_weight`` defaults to ``1 - p1 - p2 - p12``, the convention of this module;
    pass 1.0 for coloc's own convention.
    """
    with localcontext() as ctx:
        ctx.prec = 100
        bf1 = [Decimal(x).exp() for x in log_bf_1]
        bf2 = [Decimal(x).exp() for x in log_bf_2]
        n = len(bf1)
        dp1, dp2, dp12 = Decimal(p1), Decimal(p2), Decimal(p12)
        w0 = Decimal(1) - dp1 - dp2 - dp12 if h0_weight is None else Decimal(h0_weight)
        w1 = dp1 * sum(bf1, Decimal(0))
        w2 = dp2 * sum(bf2, Decimal(0))
        w3 = dp1 * dp2 * sum((bf1[i] * bf2[j] for i in range(n) for j in range(n) if i != j), Decimal(0))
        w4 = dp12 * sum((bf1[i] * bf2[i] for i in range(n)), Decimal(0))
        weights = [w0, w1, w2, w3, w4]
        total = sum(weights, Decimal(0))
        return [float(w / total) for w in weights]


def _posteriors(result: dict) -> list[float]:
    return [result[f"PP_H{i}"] for i in range(5)]


# ---------------------------------------------------------------------------
# Wakefield ABF values
# ---------------------------------------------------------------------------


class TestWakefieldAbfAgainstOracles:
    """The ABF equals Wakefield's closed form for every way of supplying V."""

    def test_standard_error_inputs_match_hand_computed_values(self) -> None:
        z, v, _, expected = ORACLE_CASES["quantitative_se"]
        got = credible_sets.compute_log_bayes_factors([z], standard_errors=math.sqrt(v))
        assert got[0] == pytest.approx(expected, rel=1e-12)

        z, v, _, expected = ORACLE_CASES["case_control_se"]
        got = credible_sets.compute_log_bayes_factors([z], standard_errors=math.sqrt(v), case_fraction=0.5)
        assert got[0] == pytest.approx(expected, rel=1e-12)

    def test_linear_scale_value_is_the_exponential(self) -> None:
        z, v, _, expected = ORACLE_CASES["quantitative_se"]
        abf = credible_sets.compute_bayes_factors([z], standard_errors=math.sqrt(v))
        assert abf[0] == pytest.approx(4.29773971708411494e2, rel=1e-12)
        assert abf[0] == pytest.approx(math.exp(expected), rel=1e-12)

    def test_sample_size_variance_models_match_hand_computed_values(self) -> None:
        z, _, _, expected = ORACLE_CASES["quantitative_n_maf"]
        got = credible_sets.compute_log_bayes_factors([z], sample_size=5000, maf=0.2)
        assert got[0] == pytest.approx(expected, rel=1e-12)

        z, _, _, expected = ORACLE_CASES["quantitative_n_standardised"]
        got = credible_sets.compute_log_bayes_factors([z], sample_size=1000)
        assert got[0] == pytest.approx(expected, rel=1e-12)

        z, _, _, expected = ORACLE_CASES["case_control_n_maf_cases"]
        got = credible_sets.compute_log_bayes_factors([z], sample_size=8000, maf=0.3, case_fraction=0.25)
        assert got[0] == pytest.approx(expected, rel=1e-12)

    def test_z_score_scale_fallback_matches_hand_computed_values(self) -> None:
        for name in ("z_scale_5", "z_scale_7_5"):
            z, _, _, expected = ORACLE_CASES[name]
            assert credible_sets.compute_log_bayes_factors([z])[0] == pytest.approx(expected, rel=1e-12)
        assert credible_sets.compute_bayes_factors([7.5])[0] == pytest.approx(1.13509524333440120e11, rel=1e-12)

    def test_default_priors_are_colocs(self) -> None:
        assert credible_sets.DEFAULT_PRIOR_VARIANCE_QUANTITATIVE == pytest.approx(0.15**2)
        assert credible_sets.DEFAULT_PRIOR_VARIANCE_CASE_CONTROL == pytest.approx(0.2**2)
        # A trait marked case-control picks the 0.2^2 default, a quantitative one 0.15^2.
        z = [5.0]
        cc = credible_sets.compute_log_bayes_factors(z, standard_errors=0.05, case_fraction=0.3)
        explicit_cc = credible_sets.compute_log_bayes_factors(z, 0.2**2, standard_errors=0.05)
        quant = credible_sets.compute_log_bayes_factors(z, standard_errors=0.05)
        explicit_quant = credible_sets.compute_log_bayes_factors(z, 0.15**2, standard_errors=0.05)
        assert cc == pytest.approx(explicit_cc, rel=1e-14)
        assert quant == pytest.approx(explicit_quant, rel=1e-14)
        assert cc != pytest.approx(quant, rel=1e-6)

    @pytest.mark.parametrize(
        ("beta_hat", "se", "prior_sd"),
        [(0.12, 0.03, 0.15), (0.25, 0.05, 0.2), (-0.4, 0.2, 0.15), (0.05, 0.04, 0.3), (1.9, 0.5, 1.0)],
    )
    def test_matches_ratio_of_marginal_likelihoods(self, beta_hat: float, se: float, prior_sd: float) -> None:
        """ABF is the Bayes factor by definition, not just a formula that looks right."""
        z = beta_hat / se
        got = credible_sets.compute_log_bayes_factors([z], prior_sd**2, standard_errors=se)[0]
        expected = math.log(_bayes_factor_from_marginal_likelihoods(beta_hat, se**2, prior_sd**2))
        assert got == pytest.approx(expected, rel=1e-10, abs=1e-12)

    @pytest.mark.parametrize(("beta_hat", "se", "prior_sd"), [(0.12, 0.03, 0.15), (0.25, 0.05, 0.2), (-0.4, 0.2, 0.15)])
    def test_matches_numerical_integration_over_the_effect(self, beta_hat: float, se: float, prior_sd: float) -> None:
        z = beta_hat / se
        got = math.exp(credible_sets.compute_log_bayes_factors([z], prior_sd**2, standard_errors=se)[0])
        assert got == pytest.approx(_bayes_factor_by_quadrature(beta_hat, se**2, prior_sd**2), rel=1e-8)

    def test_unit_changes_leave_the_abf_unchanged(self) -> None:
        """Rescaling beta, SE and the prior SD together (new units) cannot change the evidence."""
        z = [1.0, 3.5, -6.0]
        base = credible_sets.compute_log_bayes_factors(z, 0.15**2, standard_errors=0.02)
        for scale in (0.001, 7.0, 1000.0):
            scaled = credible_sets.compute_log_bayes_factors(z, (0.15 * scale) ** 2, standard_errors=0.02 * scale)
            assert scaled == pytest.approx(base, rel=1e-10)

    def test_sample_size_model_equals_the_implied_standard_error(self) -> None:
        z = [2.0, 4.5, -7.0]
        n, maf = 12000, [0.1, 0.25, 0.4]
        by_n = credible_sets.compute_log_bayes_factors(z, sample_size=n, maf=maf)
        implied_se = [1.0 / math.sqrt(2.0 * n * f * (1.0 - f)) for f in maf]
        by_se = credible_sets.compute_log_bayes_factors(z, standard_errors=implied_se)
        assert by_n == pytest.approx(by_se, rel=1e-12)

    def test_abf_is_symmetric_monotone_and_tends_to_one_for_a_vanishing_prior(self) -> None:
        values = credible_sets.compute_log_bayes_factors([-6.0, 6.0, 0.0, 2.0, 4.0], standard_errors=0.1)
        assert values[0] == pytest.approx(values[1], rel=1e-14)
        assert values[2] < values[3] < values[4] < values[1]
        # A prior with no mass away from the null cannot favour association.
        flat = credible_sets.compute_log_bayes_factors([6.0, 0.0], 1e-14, standard_errors=0.1)
        assert flat == pytest.approx([0.0, 0.0], abs=1e-9)

    def test_old_flat_default_is_gone(self) -> None:
        """Regression: z = 7.5 used to score ABF = 2.9 against about 1 for a null variant."""
        strong, null = credible_sets.compute_bayes_factors([7.5, 0.0])
        assert strong > 1e9
        assert strong / null > 1e9

    def test_huge_z_stays_finite_on_both_scales(self) -> None:
        z = [0.0, 40.0, 5000.0]
        log_values = credible_sets.compute_log_bayes_factors(z, standard_errors=0.01)
        linear = credible_sets.compute_bayes_factors(z, standard_errors=0.01)
        assert all(math.isfinite(v) for v in log_values + linear)
        assert log_values[2] > log_values[1] > log_values[0]
        assert linear[2] == pytest.approx(math.exp(500))  # capped, never inf

    def test_empty_input_gives_empty_output(self) -> None:
        assert credible_sets.compute_log_bayes_factors([]) == []
        assert credible_sets.compute_bayes_factors([], standard_errors=0.1) == []


class TestAbfInputValidation:
    """The scale of V must be stated explicitly and consistently."""

    def test_standard_errors_and_sample_size_conflict(self) -> None:
        with pytest.raises(ValueError, match="not both"):
            credible_sets.compute_log_bayes_factors([1.0], standard_errors=0.1, sample_size=100)

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"standard_errors": 0.0},
            {"standard_errors": -0.1},
            {"standard_errors": [0.1, 0.2, 0.3]},
            {"standard_errors": float("nan")},
            {"sample_size": 0},
            {"sample_size": [10, 20, 30]},
            {"sample_size": 100, "maf": 0.0},
            {"sample_size": 100, "maf": 1.0},
            {"sample_size": 100, "maf": [0.1, 0.2, 0.3]},
            {"sample_size": 100, "case_fraction": 0.0},
            {"standard_errors": 0.1, "case_fraction": 1.0},
        ],
    )
    def test_invalid_scale_inputs_are_rejected(self, kwargs: dict) -> None:
        with pytest.raises(ValueError):
            credible_sets.compute_log_bayes_factors([1.0, 2.0], **kwargs)

    @pytest.mark.parametrize("bad", [0.0, -1.0, float("inf"), float("nan")])
    def test_prior_variance_must_be_positive_and_finite(self, bad: float) -> None:
        with pytest.raises(ValueError, match="prior_variance"):
            credible_sets.compute_log_bayes_factors([1.0], bad, standard_errors=0.1)

    def test_non_finite_z_is_rejected(self) -> None:
        with pytest.raises(ValueError, match="z_scores"):
            credible_sets.compute_log_bayes_factors([1.0, float("nan")])

    def test_scalar_inputs_are_broadcast(self) -> None:
        z = [1.0, 2.0, 3.0]
        scalar = credible_sets.compute_log_bayes_factors(z, standard_errors=0.1)
        vector = credible_sets.compute_log_bayes_factors(z, standard_errors=[0.1, 0.1, 0.1])
        assert scalar == vector


# ---------------------------------------------------------------------------
# Credible sets
# ---------------------------------------------------------------------------


class TestCredibleSetsAreInformative:
    """A strong signal among nulls must yield a tight credible set."""

    def test_single_strong_signal_gives_a_tiny_credible_set(self) -> None:
        """Regression: this region used to need 285 of 300 variants."""
        for seed in range(25):
            z = _region_with_signals(seed, 300, {137: 7.5})
            result = credible_sets.compute_credible_set(z, coverage=0.95)

            assert result["status"] == "success"
            assert result["n_snps"] <= 3, f"seed {seed}: set of {result['n_snps']}"
            assert 137 in result["snps_in_set"]
            assert result["pips"][137] > 0.99
            assert math.isclose(sum(result["pips"]), 1.0, rel_tol=1e-12)

    @pytest.mark.parametrize(
        "scale_kwargs",
        [
            {"standard_errors": 0.02},
            {"sample_size": 20000},
            {"sample_size": 20000, "maf": 0.3},
            {"sample_size": 20000, "maf": 0.3, "case_fraction": 0.4},
        ],
    )
    def test_explicit_scale_inputs_also_give_tiny_sets(self, scale_kwargs: dict) -> None:
        z = _region_with_signals(3, 300, {42: 7.5})
        result = credible_sets.compute_credible_set(z, **scale_kwargs)

        assert result["status"] == "success"
        assert result["snps_in_set"] == [42]
        expected_model = "standard_errors" if "standard_errors" in scale_kwargs else "sample_size"
        assert result["variance_model"] == expected_model

    def test_z_only_calls_report_the_z_scale_model(self) -> None:
        result = credible_sets.compute_credible_set([7.5, 0.1, 0.2])
        assert result["variance_model"] == "z_scale"
        assert result["prior_variance"] == pytest.approx(5.2**2)

    def test_pips_match_independent_enumeration(self) -> None:
        z = [0.4, 3.1, -5.2, 1.7, 4.4, -0.2]
        se = 0.05
        result = credible_sets.compute_credible_set(z, standard_errors=se)
        expected = _pips_by_enumeration(z, se**2, 0.15**2)
        assert result["pips"] == pytest.approx(expected, rel=1e-10, abs=1e-14)

    def test_pips_stay_exact_beyond_double_precision_overflow(self) -> None:
        """|z| above ~37 overflows exp(z^2 / 2); PIPs must still rank and normalise exactly."""
        z = [60.0, 55.0, 0.5]
        result = credible_sets.compute_credible_set(z, standard_errors=0.01)

        assert all(math.isfinite(p) for p in result["pips"])
        assert result["pips"][0] == pytest.approx(1.0)
        assert result["pips"][1] < 1e-100
        assert result["snps_in_set"] == [0]
        assert all(math.isfinite(b) for b in result["bayes_factors"] + result["log_bayes_factors"])

    def test_comparable_signals_share_the_credible_set(self) -> None:
        """Two near-identical hits must both stay in the set (a 2-SD gap in z does not count)."""
        result = credible_sets.compute_credible_set([12.0, 11.9, 0.5, 0.5], standard_errors=0.1, coverage=0.95)
        assert result["snps_in_set"] == [0, 1]

    def test_ld_adjustment_is_applied_in_log_space(self) -> None:
        z = [60.0, 59.5, 0.3]
        ld = [[1.0, 0.99, 0.0], [0.99, 1.0, 0.0], [0.0, 0.0, 1.0]]
        result = credible_sets.compute_credible_set(z, ld_matrix=ld, standard_errors=0.01)

        assert all(math.isfinite(p) for p in result["pips"])
        assert math.isclose(sum(result["pips"]), 1.0, rel_tol=1e-12)
        # The tag variant is penalised by (1 - 0.99^2), a factor of about 50, on top of its lower z.
        unadjusted = credible_sets.compute_credible_set(z, standard_errors=0.01)
        assert result["pips"][1] < unadjusted["pips"][1]

    def test_invalid_scale_inputs_return_error_status(self) -> None:
        both = credible_sets.compute_credible_set([1.0, 2.0], standard_errors=0.1, sample_size=50)
        assert both["status"] == "error"
        assert "not both" in both["message"]
        assert credible_sets.compute_credible_set([1.0, float("nan")])["status"] == "error"
        assert credible_sets.compute_credible_set([1.0, 2.0], standard_errors=[0.1])["status"] == "error"


# ---------------------------------------------------------------------------
# Pairwise colocalization
# ---------------------------------------------------------------------------


class TestColocalizationPosteriors:
    """coloc's hypothesis posteriors from log Bayes factors."""

    CASES = [
        ([0.5, 2.0, 11.0, 3.0, -1.0, 0.0], [0.2, 1.5, 10.0, 8.0, -0.5, 0.3]),
        ([14.0, 0.1, 0.2, 0.3, 0.1, 0.2], [0.2, 0.1, 13.0, 0.3, 0.1, 0.2]),
        ([450.0, 3.0, -2.0, 1.0], [440.0, 2.0, 400.0, 0.5]),
        ([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
    ]

    @pytest.mark.parametrize(("lbf1", "lbf2"), CASES)
    def test_matches_brute_force_enumeration(self, lbf1: list[float], lbf2: list[float]) -> None:
        got = credible_sets.colocalization_from_log_bayes_factors(lbf1, lbf2)
        assert got["status"] == "success"
        assert _posteriors(got) == pytest.approx(_coloc_posteriors_by_enumeration(lbf1, lbf2), rel=1e-10, abs=1e-300)

    def test_matches_non_default_priors(self) -> None:
        lbf1, lbf2 = self.CASES[0]
        got = credible_sets.colocalization_from_log_bayes_factors(lbf1, lbf2, 2e-4, 3e-4, 5e-5)
        expected = _coloc_posteriors_by_enumeration(lbf1, lbf2, 2e-4, 3e-4, 5e-5)
        assert _posteriors(got) == pytest.approx(expected, rel=1e-10)

    def test_agrees_with_coloc_up_to_its_h0_convention(self) -> None:
        """Base-R transcription of coloc.bf_bf (H0 weight 1); this module weighs H0 by 1 - p1 - p2 - p12."""
        r_shared = [
            7.575731828464553e-05,
            4.538214342006801e-04,
            1.895075624682930e-04,
            1.360932759038109e-04,
            9.991448204091417e-01,
        ]
        r_distinct = [
            1.816148926213303e-04,
            2.184119381698570e-02,
            8.034995091925783e-03,
            9.662931264220164e-01,
            3.649069776450899e-03,
        ]
        shared = credible_sets.colocalization_from_log_bayes_factors(*self.CASES[0])
        distinct = credible_sets.colocalization_from_log_bayes_factors(*self.CASES[1])
        assert _posteriors(shared) == pytest.approx(r_shared, rel=5e-4)
        assert _posteriors(distinct) == pytest.approx(r_distinct, rel=5e-4)
        # And exactly, once the H0 convention is aligned.
        coloc_convention = _coloc_posteriors_by_enumeration(*self.CASES[0], h0_weight=1.0)
        assert coloc_convention == pytest.approx(r_shared, rel=1e-12)

    def test_dominant_variant_does_not_lose_precision(self) -> None:
        """sum(BF1) * sum(BF2) - sum(BF1 * BF2) cancels away when one variant dominates both traits."""
        lbf1, lbf2 = self.CASES[2]
        got = credible_sets.colocalization_from_log_bayes_factors(lbf1, lbf2)
        expected = _coloc_posteriors_by_enumeration(lbf1, lbf2)
        assert got["PP_H3"] == pytest.approx(expected[3], rel=1e-10)
        assert got["PP_H4"] == pytest.approx(expected[4], rel=1e-10)
        assert 0.0 < got["PP_H3"] < 1.0

    def test_one_variant_cannot_carry_distinct_causal_variants(self) -> None:
        got = credible_sets.colocalization_from_log_bayes_factors([5.0], [4.0])
        assert got["PP_H3"] == 0.0
        assert math.isclose(sum(_posteriors(got)), 1.0, rel_tol=1e-12)

    def test_snp_level_h4_probabilities_locate_the_shared_variant(self) -> None:
        lbf1, lbf2 = self.CASES[0]
        got = credible_sets.colocalization_from_log_bayes_factors(lbf1, lbf2)
        assert math.isclose(sum(got["snp_pp_h4"]), 1.0, rel_tol=1e-12)
        assert max(range(6), key=lambda i: got["snp_pp_h4"][i]) == 2

    def test_invalid_inputs_return_error_status(self) -> None:
        good = [0.0, 1.0]
        assert credible_sets.colocalization_from_log_bayes_factors([], good)["status"] == "error"
        assert credible_sets.colocalization_from_log_bayes_factors(good, [0.0])["status"] == "error"
        assert credible_sets.colocalization_from_log_bayes_factors(good, [0.0, float("inf")])["status"] == "error"
        assert credible_sets.colocalization_from_log_bayes_factors(good, [0.0, float("nan")])["status"] == "error"
        assert credible_sets.colocalization_from_log_bayes_factors(good, good, prior_p12=0.0)["status"] == "error"
        assert credible_sets.colocalization_from_log_bayes_factors(good, good, 0.5, 0.4, 0.2)["status"] == "error"

    def test_z_score_entry_point_uses_the_same_bayes_factors(self) -> None:
        z1 = _region_with_signals(1, 80, {10: 6.0})
        z2 = _region_with_signals(2, 80, {10: 5.0})
        direct = credible_sets.colocalization(z1, z2, standard_errors_1=0.02, sample_size_2=15000, maf=0.2)
        lbf1 = credible_sets.compute_log_bayes_factors(z1, standard_errors=0.02)
        lbf2 = credible_sets.compute_log_bayes_factors(z2, sample_size=15000, maf=0.2)
        via_bf = credible_sets.colocalization_from_log_bayes_factors(lbf1, lbf2)

        assert _posteriors(direct) == pytest.approx(_posteriors(via_bf), rel=1e-12)
        assert direct["bayes_factor_model"]["trait_1"]["variance_model"] == "standard_errors"
        assert direct["bayes_factor_model"]["trait_2"]["variance_model"] == "sample_size"

    def test_z_score_entry_point_rejects_bad_scale_inputs(self) -> None:
        assert credible_sets.colocalization([1.0, 2.0], [1.0, 2.0], standard_errors_1=[0.1])["status"] == "error"
        assert credible_sets.colocalization([1.0, 2.0], [1.0, 2.0], sample_size_2=-5)["status"] == "error"
        assert credible_sets.colocalization([1.0, float("nan")], [1.0, 2.0])["status"] == "error"


class TestLogSumExp:
    def test_matches_direct_sum_and_survives_extremes(self) -> None:
        assert credible_sets.log_sum_exp([math.log(2.0), math.log(3.0)]) == pytest.approx(math.log(5.0), rel=1e-14)
        assert credible_sets.log_sum_exp([1000.0, 1000.0]) == pytest.approx(1000.0 + math.log(2.0), rel=1e-14)
        assert credible_sets.log_sum_exp([-1000.0, -1001.0]) == pytest.approx(-1000.0 + math.log1p(math.exp(-1.0)))
        assert credible_sets.log_sum_exp([]) == -math.inf
        assert credible_sets.log_sum_exp([-math.inf, -math.inf]) == -math.inf
        assert credible_sets.log_sum_exp([-math.inf, 3.0]) == pytest.approx(3.0)


# ---------------------------------------------------------------------------
# eqtl_coloc: the symptom that exposed the defect
# ---------------------------------------------------------------------------


class TestEqtlColocDiscriminates:
    """A gene driven by the GWAS variant must colocalize; one driven elsewhere must not."""

    def test_shared_causal_variant_has_high_pp_h4(self) -> None:
        """Regression: PP.H4 was 0.003 here (and 0.003 for distinct variants too)."""
        for seed in range(15):
            gwas_z = _region_with_signals(seed, 300, {137: 7.5})
            eqtl_z = _region_with_signals(1000 + seed, 300, {137: -7.5})
            result = coloc_module.eqtl_coloc(gwas_z, eqtl_z, "GENE_SHARED")

            assert result["status"] == "success"
            assert result["pp_h4"] > 0.99, f"seed {seed}: PP.H4 {result['pp_h4']}"
            assert result["is_colocalized"] is True
            assert result["shared_variant"] == 137
            assert result["lead_variant_gwas"] == result["lead_variant_eqtl"] == 137

    def test_distinct_causal_variants_have_low_pp_h4(self) -> None:
        for seed in range(15):
            gwas_z = _region_with_signals(seed, 300, {137: 7.5})
            eqtl_z = _region_with_signals(1000 + seed, 300, {40: 7.5})
            result = coloc_module.eqtl_coloc(gwas_z, eqtl_z, "GENE_DISTINCT")

            assert result["status"] == "success"
            assert result["pp_h4"] < 0.01, f"seed {seed}: PP.H4 {result['pp_h4']}"
            assert result["is_colocalized"] is False
            assert result["coloc_result"]["PP_H3"] > 0.99
            assert result["lead_variant_gwas"] == 137
            assert result["lead_variant_eqtl"] == 40

    def test_explicit_scale_inputs_reach_the_bayes_factors(self) -> None:
        gwas_z = _region_with_signals(7, 200, {50: 7.0})
        eqtl_z = _region_with_signals(8, 200, {50: 9.0})
        maf = [0.05 + 0.4 * (i % 10) / 10.0 for i in range(200)]
        eqtl_se = [0.03 + 0.001 * (i % 7) for i in range(200)]
        result = coloc_module.eqtl_coloc(
            gwas_z,
            eqtl_z,
            "G",
            gwas_sample_size=50000,
            eqtl_standard_errors=eqtl_se,
            maf=maf,
            gwas_case_fraction=0.2,
        )
        lbf_gwas = credible_sets.compute_log_bayes_factors(gwas_z, sample_size=50000, maf=maf, case_fraction=0.2)
        lbf_eqtl = credible_sets.compute_log_bayes_factors(eqtl_z, standard_errors=eqtl_se)
        expected = credible_sets.colocalization_from_log_bayes_factors(lbf_gwas, lbf_eqtl)

        assert result["status"] == "success"
        assert _posteriors(result["coloc_result"]) == pytest.approx(_posteriors(expected), rel=1e-12)
        models = result["coloc_result"]["bayes_factor_model"]
        assert models["trait_1"]["variance_model"] == "sample_size"
        assert models["trait_1"]["prior_variance"] == pytest.approx(0.2**2)
        assert models["trait_2"]["variance_model"] == "standard_errors"
        assert models["trait_2"]["prior_variance"] == pytest.approx(0.15**2)

    def test_explicit_prior_variances_are_used(self) -> None:
        z = _region_with_signals(4, 60, {5: 6.0})
        result = coloc_module.eqtl_coloc(z, z, "G", gwas_prior_variance=0.09, eqtl_prior_variance=0.01)
        models = result["coloc_result"]["bayes_factor_model"]
        assert models["trait_1"]["prior_variance"] == 0.09
        assert models["trait_2"]["prior_variance"] == 0.01

    def test_bad_scale_inputs_return_error_status(self) -> None:
        z = [1.0, 2.0, 3.0]
        result = coloc_module.eqtl_coloc(z, z, "G", gwas_standard_errors=[0.1, 0.1])
        assert result["status"] == "error"
        assert result["gene_id"] == "G"

    @pytest.mark.parametrize(
        "strength", [0.0, 1.0, 3.0, 5.0, 8.0, 12.0, 20.0, 30.0, 37.0, 38.0, 50.0, 100.0, 300.0, 1000.0]
    )
    def test_results_are_finite_across_the_signal_range(self, strength: float) -> None:
        """Correct ABFs reach exp(z^2 / 2); a linear-scale sum overflows to NaN above |z| ~ 26."""
        n, causal = 120, 60
        gwas_z = _region_with_signals(11, n, {causal: strength})
        eqtl_z = _region_with_signals(12, n, {causal: -strength})
        result = coloc_module.eqtl_coloc(gwas_z, eqtl_z, "GENE")

        assert result["status"] == "success"
        posteriors = _posteriors(result["coloc_result"])
        assert all(math.isfinite(p) and 0.0 <= p <= 1.0 for p in posteriors)
        assert math.isclose(sum(posteriors), 1.0, rel_tol=1e-12)
        assert math.isfinite(result["pp_h4"])
        assert all(math.isfinite(p) for p in result["coloc_result"]["snp_pp_h4"])
        if strength >= 8.0:
            assert result["pp_h4"] > 0.99
            assert result["shared_variant"] == causal

    def test_strong_shared_signals_of_30_sd_pick_the_shared_variant(self) -> None:
        gwas_z = _region_with_signals(21, 200, {75: 30.0, 76: 29.0})
        eqtl_z = _region_with_signals(22, 200, {75: -31.0, 76: -30.0})
        result = coloc_module.eqtl_coloc(gwas_z, eqtl_z, "STRONG")

        assert math.isfinite(result["pp_h4"])
        assert result["pp_h4"] == pytest.approx(1.0, abs=1e-6)
        assert result["shared_variant"] == 75  # |z| 30/31 beats 29/30 in both traits
        assert result["is_colocalized"] is True

    def test_strong_distinct_signals_of_30_sd_stay_low(self) -> None:
        gwas_z = _region_with_signals(23, 200, {75: 30.0})
        eqtl_z = _region_with_signals(24, 200, {20: -30.0})
        result = coloc_module.eqtl_coloc(gwas_z, eqtl_z, "DISTINCT")

        assert math.isfinite(result["pp_h4"])
        assert result["pp_h4"] < 1e-6
        assert result["coloc_result"]["PP_H3"] > 0.999


# ---------------------------------------------------------------------------
# multi_trait_coloc and regional_coloc
# ---------------------------------------------------------------------------


def _multi_trait_oracle(log_bfs: dict[str, list[float]], p1: float = 1e-4, p12: float = 1e-5) -> tuple[float, float]:
    """``(P(all traits share one variant), P(every trait has a causal variant))`` by enumerating every configuration.

    A configuration gives each trait no causal variant or one variant; traits at the same variant form a group of
    prior ``p1 * (p12 / p1) ** (k - 1)`` for ``k`` traits (coloc's ``p1`` and ``p12`` for one and two), and the Bayes
    factors of the traits multiply. "No trait associated" weighs one minus the prior mass of every non-empty trait
    set, the convention of ``colocalization``. This is the definition of the model, with no dynamic programme.
    """
    with localcontext() as ctx:
        ctx.prec = 100
        names = list(log_bfs)
        n_traits = len(names)
        n_variants = len(log_bfs[names[0]])
        bayes = [[Decimal(x).exp() for x in log_bfs[name]] for name in names]
        d1, d12 = Decimal(p1), Decimal(p12)

        def group_prior(size: int) -> Decimal:
            return d1 * (d12 / d1) ** (size - 1)

        total = Decimal(1) - sum((math.comb(n_traits, k) * group_prior(k) for k in range(1, n_traits + 1)), Decimal(0))
        shared = all_associated = Decimal(0)
        for config in itertools.product(range(n_variants + 1), repeat=n_traits):  # 0: no causal variant
            if not any(config):
                continue
            weight = Decimal(1)
            for size in Counter(c for c in config if c).values():
                weight *= group_prior(size)
            for trait, variant in enumerate(config):
                if variant:
                    weight *= bayes[trait][variant - 1]
            total += weight
            if all(config):
                all_associated += weight
                if len(set(config)) == 1:
                    shared += weight
        return float(shared / total), float(all_associated / total)


class TestMultiTraitColoc:
    def test_pp_shared_matches_independent_enumeration(self) -> None:
        z = {
            "a": _region_with_signals(31, 6, {2: 4.5}),
            "b": _region_with_signals(32, 6, {2: 4.0}),
            "c": _region_with_signals(33, 6, {2: 3.5, 4: 2.0}),
        }
        result = coloc_module.multi_trait_coloc(z, standard_errors={"a": 0.05, "b": 0.05}, sample_sizes={"c": 3000})
        log_bfs = {
            "a": credible_sets.compute_log_bayes_factors(z["a"], standard_errors=0.05),
            "b": credible_sets.compute_log_bayes_factors(z["b"], standard_errors=0.05),
            "c": credible_sets.compute_log_bayes_factors(z["c"], sample_size=3000),
        }
        shared, all_associated = _multi_trait_oracle(log_bfs)
        assert result["status"] == "success"
        assert result["multi_trait_pp_shared"] == pytest.approx(shared, rel=1e-10)
        assert result["multi_trait_pp_all_associated"] == pytest.approx(all_associated, rel=1e-10)

    def test_z_only_traits_match_the_z_scale_oracle(self) -> None:
        z = {"a": _region_with_signals(34, 8, {3: 6.0}), "b": _region_with_signals(35, 8, {3: 6.5})}
        result = coloc_module.multi_trait_coloc(z)
        log_bfs = {name: credible_sets.compute_log_bayes_factors(values) for name, values in z.items()}
        assert result["multi_trait_pp_shared"] == pytest.approx(_multi_trait_oracle(log_bfs)[0], rel=1e-10)
        assert result["multi_trait_pp_shared"] > 0.99
        assert set(result["pairwise"]) == {"a_vs_b"}

    def test_shared_signals_of_30_sd_across_three_traits_are_finite_and_near_one(self) -> None:
        z = {
            "a": _region_with_signals(41, 150, {80: 30.0}),
            "b": _region_with_signals(42, 150, {80: -30.0}),
            "c": _region_with_signals(43, 150, {80: 31.0}),
        }
        result = coloc_module.multi_trait_coloc(z)

        assert result["status"] == "success"
        assert math.isfinite(result["multi_trait_pp_shared"])
        assert result["multi_trait_pp_shared"] == pytest.approx(1.0, abs=1e-9)
        for pair in result["pairwise"].values():
            assert pair["PP_H4"] == pytest.approx(1.0, abs=1e-6)
            assert max(range(150), key=lambda i: pair["snp_pp_h4"][i]) == 80

    def test_strong_distinct_signals_do_not_colocalize(self) -> None:
        z = {
            "a": _region_with_signals(44, 150, {80: 30.0}),
            "b": _region_with_signals(45, 150, {20: 30.0}),
        }
        result = coloc_module.multi_trait_coloc(z)

        pair = result["pairwise"]["a_vs_b"]
        assert pair["PP_H4"] < 1e-6
        assert pair["PP_H3"] > 0.999
        assert math.isfinite(result["multi_trait_pp_shared"])
        assert result["multi_trait_pp_shared"] < 0.5
        assert "Limited evidence" in result["summary"]

    @pytest.mark.parametrize("strength", [0.0, 2.0, 6.0, 15.0, 30.0, 40.0, 80.0, 400.0])
    def test_results_are_finite_across_the_signal_range(self, strength: float) -> None:
        z = {
            "a": _region_with_signals(51, 60, {10: strength}),
            "b": _region_with_signals(52, 60, {10: strength}),
            "c": _region_with_signals(53, 60, {10: -strength}),
        }
        result = coloc_module.multi_trait_coloc(z)

        assert result["status"] == "success"
        assert math.isfinite(result["multi_trait_pp_shared"])
        assert 0.0 <= result["multi_trait_pp_shared"] <= 1.0
        for pair in result["pairwise"].values():
            assert all(math.isfinite(v) for v in _posteriors(pair))
            assert math.isclose(sum(_posteriors(pair)), 1.0, rel_tol=1e-12)

    def test_unknown_trait_names_are_rejected(self) -> None:
        z = {"a": [1.0, 2.0], "b": [1.0, 2.0]}
        result = coloc_module.multi_trait_coloc(z, standard_errors={"a": 0.1, "typo": 0.1})
        assert result["status"] == "error"
        assert "typo" in result["message"]

    def test_bad_scale_inputs_return_error_status(self) -> None:
        z = {"a": [1.0, 2.0], "b": [1.0, 2.0]}
        assert coloc_module.multi_trait_coloc(z, sample_sizes={"a": -1})["status"] == "error"
        assert coloc_module.multi_trait_coloc(z, standard_errors={"a": [0.1, 0.1, 0.1]})["status"] == "error"

    def test_case_fraction_selects_the_case_control_prior_per_trait(self) -> None:
        z = {"a": _region_with_signals(61, 30, {2: 5.0}), "b": _region_with_signals(62, 30, {2: 5.0})}
        result = coloc_module.multi_trait_coloc(
            z, standard_errors={"a": 0.05, "b": 0.05}, case_fractions={"a": 0.3}, prior_variances={"b": 0.01}
        )
        models = result["pairwise"]["a_vs_b"]["bayes_factor_model"]
        assert models["trait_1"]["prior_variance"] == pytest.approx(0.2**2)
        assert models["trait_2"]["prior_variance"] == 0.01


class TestRegionalColoc:
    def test_per_variant_scale_inputs_are_sliced_per_window(self) -> None:
        n = 250
        z = {"a": _region_with_signals(71, n, {120: 6.0}), "b": _region_with_signals(72, n, {120: 6.0})}
        se_a = [0.02 + 0.0001 * (i % 13) for i in range(n)]
        maf = [0.1 + 0.3 * (i % 5) / 5.0 for i in range(n)]
        region = {"chrom": "1", "start": 0, "end": n * 1000}
        result = coloc_module.regional_coloc(
            z, region, standard_errors={"a": se_a}, sample_sizes={"b": 20000}, maf=maf, case_fractions={"b": 0.5}
        )

        assert result["status"] == "success"
        assert result["n_windows"] > 1
        for window in result["all_windows"]:
            lo, hi = window["window_start_idx"], window["window_end_idx"]
            expected = coloc_module.multi_trait_coloc(
                {name: values[lo:hi] for name, values in z.items()},
                standard_errors={"a": se_a[lo:hi]},
                sample_sizes={"b": 20000},
                maf=maf[lo:hi],
                case_fractions={"b": 0.5},
            )
            assert window["multi_trait_pp_shared"] == pytest.approx(expected["multi_trait_pp_shared"], rel=1e-12)
        best = result["best_window"]
        assert best["window_start_idx"] <= 120 < best["window_end_idx"]
        assert best["multi_trait_pp_shared"] > 0.99

    def test_per_variant_inputs_must_cover_the_region(self) -> None:
        z = {"a": [1.0] * 30, "b": [1.0] * 30}
        region = {"chrom": "1", "start": 0, "end": 30000}
        too_short = coloc_module.regional_coloc(z, region, standard_errors={"a": [0.1] * 29})
        assert too_short["status"] == "error"
        assert "standard_errors" in too_short["message"]
        assert coloc_module.regional_coloc(z, region, maf=[0.2] * 31)["status"] == "error"

    def test_strong_signals_do_not_overflow_in_any_window(self) -> None:
        n = 220
        z = {"a": _region_with_signals(73, n, {110: 40.0}), "b": _region_with_signals(74, n, {110: -45.0})}
        result = coloc_module.regional_coloc(z, {"chrom": "2", "start": 0, "end": n * 500})

        assert result["status"] == "success"
        assert all(math.isfinite(w["multi_trait_pp_shared"]) for w in result["all_windows"])
        assert result["best_window"]["multi_trait_pp_shared"] == pytest.approx(1.0, abs=1e-9)
