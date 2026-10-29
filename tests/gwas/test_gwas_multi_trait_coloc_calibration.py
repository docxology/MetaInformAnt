"""Regression tests: ``multi_trait_pp_shared`` is a calibrated posterior, not a heuristic score.

Defect: ``multi_trait_coloc`` reported ``sum_v prod_t BF / (n + sum_t sum_v BF + sum_v prod_t BF)`` as "the posterior
probability that all traits share a single causal variant". That gives the null, "some trait associated" and "all
share" equal prior weight, lets a trait add its Bayes factors to the "individual" term, and has no hypothesis of
distinct causal variants. A trait with no signal multiplies the product term by a Bayes factor near 1, so when two
traits share a signal a *third trait with nothing in the region* leaves the statistic at about 1.

Measured on 400 simulated regions per scenario (100 variants in AR(1) LD with r = 0.9 between neighbours, signal
strength z = 5 to 9, default z-score scale; threshold 0.8), the heuristic called these "shared" and the posterior
does not:

* three traits, two share and the third has no signal: 99.0% (now 0.0%);
* three traits, two share and the third has another signal far away: 87.3% (now 0.0%);
* two traits with distinct signals far apart: 2.0% (now 0.0%); one trait with a signal: 1.5% (now 0.0%);
* the traits that do share one variant are still called: 100% (now 98%), the cost of calling it a probability.

The posterior is the exact answer of coloc's model generalised to ``N`` traits, and for two traits it *is* coloc's
``PP.H4``. Oracles: ``colocalization`` (two traits), enumeration of every configuration in 100-digit decimals (more
traits, in ``test_gwas_finemapping_abf.py`` and here), and the simulated truth.
"""

from __future__ import annotations

import itertools
import math
from collections import Counter
from decimal import Decimal, localcontext
from typing import Dict, List, Optional

import numpy as np
import pytest

from metainformant.gwas.finemapping import colocalization as coloc_module
from metainformant.gwas.finemapping.credible_sets import colocalization as pairwise_coloc, compute_log_bayes_factors

#: The thresholds the verdict of ``multi_trait_coloc`` uses: coloc's conventional 0.8, and 0.5 for "moderate".
STRONG = 0.8
MODERATE = 0.5


# ---------------------------------------------------------------------------
# A region simulator: AR(1) linkage disequilibrium, one causal variant per trait
# ---------------------------------------------------------------------------

N_VARIANTS = 100
RHO = 0.9  # correlation between neighbouring variants; r = RHO ** distance


def trait_z(rng: np.random.Generator, causal: Optional[int], strength: float) -> List[float]:
    """z-scores of one trait: LD-correlated noise (stationary AR(1)) plus ``strength * r(., causal)``."""
    noise = np.empty(N_VARIANTS)
    noise[0] = rng.standard_normal()
    step = math.sqrt(1.0 - RHO**2)
    for i in range(1, N_VARIANTS):
        noise[i] = RHO * noise[i - 1] + step * rng.standard_normal()
    if causal is not None:
        noise = noise + strength * RHO ** np.abs(np.arange(N_VARIANTS) - causal)
    return noise.tolist()


def far_from(rng: np.random.Generator, taken: List[int], distance: int = 20) -> int:
    while True:
        candidate = int(rng.integers(0, N_VARIANTS))
        if all(abs(candidate - t) >= distance for t in taken):
            return candidate


def causal_variants(rng: np.random.Generator, scenario: str, n_traits: int) -> List[Optional[int]]:
    """Causal variant (or None) of each trait in a named scenario."""
    c = int(rng.integers(10, N_VARIANTS - 10))
    if scenario == "shared":
        return [c] * n_traits
    if scenario == "distinct_far":
        taken = [c]
        while len(taken) < n_traits:
            taken.append(far_from(rng, taken))
        return list(taken)
    if scenario == "two_share_one_distinct":
        return [c, c] + [far_from(rng, [c])] * (n_traits - 2)
    if scenario == "all_but_one_share_one_null":
        return [c] * (n_traits - 1) + [None]
    if scenario == "one_trait_only":
        return [c] + [None] * (n_traits - 1)
    if scenario == "null":
        return [None] * n_traits
    raise ValueError(scenario)


def rate_above(threshold: float, scenario: str, n_traits: int, replicates: int = 80, seed: int = 11) -> float:
    """Share of simulated regions of a scenario whose ``multi_trait_pp_shared`` exceeds ``threshold``."""
    rng = np.random.default_rng(seed)
    flagged = 0
    for _ in range(replicates):
        causals = causal_variants(rng, scenario, n_traits)
        strengths = rng.uniform(5.0, 9.0, size=n_traits)
        z = {f"trait{k}": trait_z(rng, causals[k], float(strengths[k])) for k in range(n_traits)}
        result = coloc_module.multi_trait_coloc(z)
        assert result["status"] == "success", result
        flagged += result["multi_trait_pp_shared"] > threshold
    return flagged / replicates


# ---------------------------------------------------------------------------
# Calibration on simulated truth
# ---------------------------------------------------------------------------


class TestCalibration:
    @pytest.mark.parametrize("n_traits", [2, 3])
    def test_traits_that_share_a_variant_are_called(self, n_traits: int) -> None:
        assert rate_above(STRONG, "shared", n_traits) >= 0.92
        assert rate_above(MODERATE, "shared", n_traits) >= 0.95

    def test_a_trait_without_a_signal_is_not_called_shared(self) -> None:
        """The old statistic stayed near 1 here (99% of regions above 0.8): a null trait multiplies by BF ~ 1."""
        assert rate_above(STRONG, "all_but_one_share_one_null", 3) <= 0.05
        assert rate_above(MODERATE, "all_but_one_share_one_null", 3) <= 0.05

    def test_a_trait_with_a_distinct_signal_is_not_called_shared(self) -> None:
        """The old statistic had no distinct-variants hypothesis (87% of regions above 0.8)."""
        assert rate_above(STRONG, "two_share_one_distinct", 3) <= 0.05
        assert rate_above(MODERATE, "two_share_one_distinct", 3) <= 0.05

    @pytest.mark.parametrize("n_traits", [2, 3])
    def test_distinct_signals_far_apart_are_not_called_shared(self, n_traits: int) -> None:
        assert rate_above(STRONG, "distinct_far", n_traits) <= 0.03
        assert rate_above(MODERATE, "distinct_far", n_traits) <= 0.03

    @pytest.mark.parametrize("n_traits", [2, 3])
    def test_one_trait_with_a_signal_is_not_called_shared(self, n_traits: int) -> None:
        assert rate_above(STRONG, "one_trait_only", n_traits) <= 0.03
        assert rate_above(MODERATE, "one_trait_only", n_traits) <= 0.03

    @pytest.mark.parametrize("n_traits", [2, 3])
    def test_a_region_with_no_signal_is_never_called_shared(self, n_traits: int) -> None:
        assert rate_above(MODERATE, "null", n_traits) == 0.0

    def test_the_thresholds_are_the_documented_ones(self) -> None:
        assert coloc_module.COLOC_STRONG_THRESHOLD == STRONG
        assert coloc_module.COLOC_MODERATE_THRESHOLD == MODERATE

    def test_the_verdict_follows_the_posterior(self) -> None:
        rng = np.random.default_rng(3)
        shared = {f"t{k}": trait_z(rng, 50, 8.0) for k in range(3)}
        with_null = {**dict(list(shared.items())[:2]), "t2": trait_z(rng, None, 0.0)}

        assert "Strong evidence" in coloc_module.multi_trait_coloc(shared)["summary"]
        assert "Strong evidence" not in coloc_module.multi_trait_coloc(with_null)["summary"]


# ---------------------------------------------------------------------------
# Exactness
# ---------------------------------------------------------------------------


def enumerate_posteriors(log_bfs: Dict[str, List[float]], p1: float = 1e-4, p12: float = 1e-5) -> tuple[float, float]:
    """(P(all share), P(all have a causal variant)) from the definition: every configuration, 100 digits."""
    with localcontext() as ctx:
        ctx.prec = 100
        names = list(log_bfs)
        n_traits, n_variants = len(names), len(log_bfs[names[0]])
        bayes = [[Decimal(x).exp() for x in log_bfs[name]] for name in names]
        d1, d12 = Decimal(p1), Decimal(p12)

        def prior(size: int) -> Decimal:
            return d1 * (d12 / d1) ** (size - 1)

        total = Decimal(1) - sum((math.comb(n_traits, k) * prior(k) for k in range(1, n_traits + 1)), Decimal(0))
        shared = associated = Decimal(0)
        for config in itertools.product(range(n_variants + 1), repeat=n_traits):
            if not any(config):
                continue
            weight = Decimal(1)
            for size in Counter(c for c in config if c).values():
                weight *= prior(size)
            for trait, variant in enumerate(config):
                if variant:
                    weight *= bayes[trait][variant - 1]
            total += weight
            if all(config):
                associated += weight
                if len(set(config)) == 1:
                    shared += weight
        return float(shared / total), float(associated / total)


def random_z(seed: int, n_traits: int, n_variants: int, planted: Optional[int] = None) -> Dict[str, List[float]]:
    rng = np.random.default_rng(seed)
    z = {f"t{k}": (rng.normal(size=n_variants) * 1.2).tolist() for k in range(n_traits)}
    if planted is not None:
        for k, values in enumerate(z.values()):
            values[planted] = 5.0 + k * 0.7
    return z


class TestExactness:
    @pytest.mark.parametrize("seed", range(6))
    def test_two_traits_are_coloc_pp_h4(self, seed: int) -> None:
        z = random_z(seed, 2, 40, planted=seed % 40)

        result = coloc_module.multi_trait_coloc(z)
        pair = pairwise_coloc(z["t0"], z["t1"])

        assert result["multi_trait_pp_shared"] == pytest.approx(pair["PP_H4"], abs=1e-12)

    def test_two_traits_are_coloc_pp_h4_with_scale_inputs_and_custom_priors(self) -> None:
        z = random_z(7, 2, 25, planted=4)
        scale = {"standard_errors": {"t0": 0.05, "t1": 0.04}}

        result = coloc_module.multi_trait_coloc(z, prior_p1=2e-4, prior_p12=3e-5, **scale)
        pair = pairwise_coloc(
            z["t0"],
            z["t1"],
            prior_p1=2e-4,
            prior_p2=2e-4,
            prior_p12=3e-5,
            standard_errors_1=0.05,
            standard_errors_2=0.04,
        )

        assert result["multi_trait_pp_shared"] == pytest.approx(pair["PP_H4"], abs=1e-12)
        assert result["pairwise"]["t0_vs_t1"]["PP_H4"] == pytest.approx(pair["PP_H4"], abs=1e-12)

    @pytest.mark.parametrize("n_traits, n_variants", [(3, 5), (4, 4)])
    def test_more_traits_match_enumeration_of_every_configuration(self, n_traits: int, n_variants: int) -> None:
        z = random_z(n_traits, n_traits, n_variants, planted=1)
        log_bfs = {name: compute_log_bayes_factors(values) for name, values in z.items()}

        result = coloc_module.multi_trait_coloc(z)
        shared, associated = enumerate_posteriors(log_bfs)

        assert result["multi_trait_pp_shared"] == pytest.approx(shared, rel=1e-10)
        assert result["multi_trait_pp_all_associated"] == pytest.approx(associated, rel=1e-10)

    def test_custom_priors_match_enumeration(self) -> None:
        z = random_z(9, 3, 4, planted=2)
        log_bfs = {name: compute_log_bayes_factors(values) for name, values in z.items()}

        result = coloc_module.multi_trait_coloc(z, prior_p1=5e-3, prior_p12=2e-3)

        assert result["multi_trait_pp_shared"] == pytest.approx(enumerate_posteriors(log_bfs, 5e-3, 2e-3)[0], rel=1e-10)
        assert result["priors"]["p1"] == 5e-3 and result["priors"]["p12"] == 2e-3

    def test_the_null_weight_is_what_is_left_of_the_prior(self) -> None:
        result = coloc_module.multi_trait_coloc(random_z(1, 3, 4))

        # three singles, three pairs, one triple: p1 = 1e-4, p12 = 1e-5, p123 = p12 ** 2 / p1 = 1e-6
        assert result["priors"]["p0"] == pytest.approx(1 - 3e-4 - 3e-5 - 1e-6, rel=1e-12)

    def test_a_larger_shared_prior_raises_the_probability(self) -> None:
        z = random_z(2, 3, 30, planted=6)

        low = coloc_module.multi_trait_coloc(z, prior_p12=1e-6)["multi_trait_pp_shared"]
        high = coloc_module.multi_trait_coloc(z, prior_p12=5e-5)["multi_trait_pp_shared"]

        assert 0.0 <= low < high <= 1.0


# ---------------------------------------------------------------------------
# What the numbers must satisfy
# ---------------------------------------------------------------------------


class TestInvariants:
    def test_all_associated_bounds_shared_and_per_variant_shares_sum_to_one(self) -> None:
        z = random_z(5, 3, 30, planted=11)

        result = coloc_module.multi_trait_coloc(z)

        assert 0.0 <= result["multi_trait_pp_shared"] <= result["multi_trait_pp_all_associated"] <= 1.0
        assert sum(result["snp_pp_shared"]) == pytest.approx(1.0, abs=1e-12)
        assert max(range(30), key=lambda i: result["snp_pp_shared"][i]) == 11

    def test_adding_a_trait_without_a_signal_lowers_the_probability(self) -> None:
        pair = random_z(4, 2, 40, planted=13)
        triple = {**pair, "t2": random_z(99, 1, 40)["t0"]}

        two = coloc_module.multi_trait_coloc(pair)["multi_trait_pp_shared"]
        three = coloc_module.multi_trait_coloc(triple)["multi_trait_pp_shared"]

        assert two > 0.99
        assert three < 0.05  # the old statistic stayed above 0.99

    def test_a_trait_with_another_signal_far_away_gives_distinct_variants(self) -> None:
        pair = random_z(4, 2, 40, planted=13)
        other = random_z(98, 1, 40, planted=33)["t0"]

        result = coloc_module.multi_trait_coloc({**pair, "t2": other})

        assert result["multi_trait_pp_shared"] < 0.01
        # all three are probably associated, just not at one variant
        assert result["multi_trait_pp_all_associated"] > 0.5
        assert result["multi_trait_pp_all_associated"] > 50 * result["multi_trait_pp_shared"]

    def test_trait_order_does_not_change_the_probability(self) -> None:
        z = random_z(6, 3, 20, planted=3)
        reordered = {name: z[name] for name in ("t2", "t0", "t1")}

        a = coloc_module.multi_trait_coloc(z)["multi_trait_pp_shared"]
        b = coloc_module.multi_trait_coloc(reordered)["multi_trait_pp_shared"]

        assert a == pytest.approx(b, rel=1e-12)

    @pytest.mark.parametrize("strength", [0.0, 6.0, 80.0, 400.0])
    def test_four_traits_stay_finite_whatever_the_signal(self, strength: float) -> None:
        z = random_z(8, 4, 15)
        for values in z.values():
            values[4] = strength

        result = coloc_module.multi_trait_coloc(z)

        assert result["status"] == "success"
        assert math.isfinite(result["multi_trait_pp_shared"]) and 0.0 <= result["multi_trait_pp_shared"] <= 1.0
        if strength >= 80.0:
            assert result["multi_trait_pp_shared"] == pytest.approx(1.0, abs=1e-9)

    def test_ld_matrix_is_accepted_and_changes_nothing(self) -> None:
        z = random_z(3, 2, 6, planted=2)
        identity = np.eye(6).tolist()

        assert coloc_module.multi_trait_coloc(z, ld_matrix=identity) == coloc_module.multi_trait_coloc(z)


class TestInputs:
    @pytest.mark.parametrize("bad", [0.0, -1e-4, float("nan"), float("inf"), True, "1e-4", None])
    def test_a_bad_prior_is_an_error(self, bad) -> None:
        for keyword in ("prior_p1", "prior_p12"):
            result = coloc_module.multi_trait_coloc(random_z(0, 2, 5), **{keyword: bad})

            assert result["status"] == "error"
            assert keyword in result["message"]

    def test_priors_that_sum_to_more_than_one_are_an_error(self) -> None:
        result = coloc_module.multi_trait_coloc(random_z(0, 2, 5), prior_p1=0.4, prior_p12=0.4)

        assert result["status"] == "error"
        assert "more than 1" in result["message"]

    def test_more_traits_than_the_exact_limit_are_an_error_not_an_approximation(self) -> None:
        n = coloc_module.MAX_EXACT_TRAITS + 1

        result = coloc_module.multi_trait_coloc(random_z(0, n, 3))

        assert result["status"] == "error"
        assert str(coloc_module.MAX_EXACT_TRAITS) in result["message"] and "subsets" in result["message"]

    def test_the_limit_itself_is_accepted(self) -> None:
        result = coloc_module.multi_trait_coloc(random_z(0, 4, 3))

        assert result["status"] == "success"

    def test_regional_coloc_reports_the_new_statistic(self) -> None:
        z = random_z(12, 3, 60, planted=30)
        region = {"chrom": "1", "start": 0, "end": 60_000}

        result = coloc_module.regional_coloc(z, region)

        assert result["status"] == "success"
        window = result["best_window"]
        assert window["multi_trait_pp_shared"] == pytest.approx(
            coloc_module.multi_trait_coloc(z)["multi_trait_pp_shared"], rel=1e-12
        )
