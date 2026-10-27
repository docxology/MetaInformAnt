"""Regression tests for the ``bayesian`` model of ``call_methylation_from_signal``.

The option used to be ``(1 + p) / 3`` above 0.5 and ``p`` below it: a capped linear map
(never above 2/3) that ignored every other read at the site. It is now a hierarchical
Bayesian posterior. Reads at one ``(chrom, position, strand)`` share a methylation rate with a
Beta prior, and each read's probability is its posterior given the LLRs of all reads at the site.

The tests check that posterior against an exact enumeration of every methylated/unmethylated
state of the reads (marginalising the rate analytically through the beta-binomial), against an
exact recursion for larger sites, and for the behaviours that follow from the model. Everything
runs on real computation; there are no mocks.
"""

from __future__ import annotations

import itertools
import math
import random
from typing import Any

import pytest

from metainformant.longread.methylation.calling import aggregate_methylation, call_methylation_from_signal


def _read(read_id: str, llr: float, position: int = 100, chrom: str = "chr1", strand: str = "+") -> dict[str, Any]:
    return {
        "read_id": read_id,
        "chrom": chrom,
        "position": position,
        "strand": strand,
        "log_likelihood_ratio": llr,
    }


def _site(llrs: list[float], position: int = 100, chrom: str = "chr1", strand: str = "+") -> list[dict[str, Any]]:
    return [_read(f"{chrom}:{position}:{strand}:{i}", llr, position, chrom, strand) for i, llr in enumerate(llrs)]


def _probability(call: dict[str, Any]) -> float:
    """Recover the probability behind a call made at threshold 0.5: confidence is twice its distance from 0.5."""
    return 0.5 + call["confidence"] / 2 if call["methylated"] else 0.5 - call["confidence"] / 2


def _sigmoid(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x))


def _log_beta(a: float, b: float) -> float:
    return math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)


def _enumerated_posteriors(llrs: list[float], alpha: float, beta: float) -> list[float]:
    """Exact P(read i methylated | all LLRs): sum over every state z of the reads.

    Marginalising theta ~ Beta(alpha, beta) analytically, P(z | LLRs) is proportional to
    ``exp(sum of LLR_i over methylated reads) * B(alpha + k, beta + n - k)`` with k methylated reads.
    """
    n = len(llrs)
    states = list(itertools.product((0, 1), repeat=n))
    log_weights = [
        sum(llr for llr, z in zip(llrs, state) if z) + _log_beta(alpha + sum(state), beta + n - sum(state))
        for state in states
    ]
    peak = max(log_weights)
    weights = [math.exp(w - peak) for w in log_weights]
    total = sum(weights)
    return [sum(w for w, state in zip(weights, states) if state[i]) / total for i in range(n)]


def _recursive_posteriors(llrs: list[float], alpha: float, beta: float) -> list[float]:
    """Exact posteriors through elementary symmetric polynomials of the likelihood ratios (small |LLR| only)."""
    n = len(llrs)
    ratios = [math.exp(llr) for llr in llrs]

    def symmetric(values: list[float]) -> list[float]:
        e = [1.0] + [0.0] * len(values)
        for value in values:
            for k in range(len(values), 0, -1):
                e[k] += value * e[k - 1]
        return e

    prior_k = [math.exp(_log_beta(alpha + k, beta + n - k) - _log_beta(alpha, beta)) for k in range(n + 1)]
    normaliser = sum(e * w for e, w in zip(symmetric(ratios), prior_k))
    posteriors = []
    for i in range(n):
        others = symmetric(ratios[:i] + ratios[i + 1 :])
        posteriors.append(ratios[i] * sum(others[k - 1] * prior_k[k] for k in range(1, n + 1)) / normaliser)
    return posteriors


class TestSingleRead:
    @pytest.mark.parametrize("llr", [-30.0, -3.0, -0.1, 0.0, 0.1, 2.0, 10.0, 45.0])
    def test_flat_prior_gives_the_threshold_probability(self, llr: float) -> None:
        (threshold_call,) = call_methylation_from_signal([_read("r", llr)], model="threshold")
        (bayes_call,) = call_methylation_from_signal([_read("r", llr)], model="bayesian")
        assert bayes_call["methylated"] == threshold_call["methylated"]
        assert bayes_call["confidence"] == pytest.approx(threshold_call["confidence"], abs=1e-12)

    @pytest.mark.parametrize(("alpha", "beta"), [(7.0, 3.0), (1.0, 4.0), (0.5, 0.5), (30.0, 10.0)])
    @pytest.mark.parametrize("llr", [-4.0, -0.5, 0.0, 0.5, 4.0])
    def test_informative_prior_shifts_the_odds_by_the_prior_odds(self, llr: float, alpha: float, beta: float) -> None:
        (call,) = call_methylation_from_signal([_read("r", llr)], model="bayesian", prior_alpha=alpha, prior_beta=beta)
        assert _probability(call) == pytest.approx(_sigmoid(llr + math.log(alpha / beta)), abs=1e-9)

    def test_a_prior_that_favours_methylation_turns_a_weakly_negative_read_methylated(self) -> None:
        read = _read("r", -0.5)
        assert call_methylation_from_signal([read], model="threshold")[0]["methylated"] is False
        assert call_methylation_from_signal([read], model="bayesian")[0]["methylated"] is False
        assert call_methylation_from_signal([read], model="bayesian", prior_alpha=7, prior_beta=3)[0]["methylated"]

    def test_the_old_two_thirds_cap_is_gone(self) -> None:
        (call,) = call_methylation_from_signal([_read("r", 12.0)], model="bayesian", threshold=0.9)
        assert call["methylated"] is True
        assert call["confidence"] > 0.99


class TestPosteriorIsExact:
    @pytest.mark.parametrize(("alpha", "beta"), [(1.0, 1.0), (0.5, 0.5), (2.0, 5.0), (8.0, 1.5), (0.25, 3.0)])
    def test_matches_exhaustive_enumeration_of_the_read_states(self, alpha: float, beta: float) -> None:
        rng = random.Random(f"enumeration-{alpha}-{beta}")
        for _ in range(12):
            n = rng.randint(2, 9)
            scale = rng.choice([0.5, 2.0, 6.0, 15.0])
            llrs = [rng.gauss(rng.choice([-1, 0, 1]) * scale / 2, scale) for _ in range(n)]
            calls = call_methylation_from_signal(_site(llrs), model="bayesian", prior_alpha=alpha, prior_beta=beta)
            expected = _enumerated_posteriors(llrs, alpha, beta)
            assert [_probability(c) for c in calls] == pytest.approx(expected, abs=1e-6)

    def test_matches_the_exact_recursion_on_a_larger_site(self) -> None:
        rng = random.Random(30)
        llrs = [rng.gauss(rng.choice([-2.0, 2.0]), 2.0) for _ in range(30)]
        for alpha, beta in ((1.0, 1.0), (2.0, 5.0)):
            calls = call_methylation_from_signal(_site(llrs), model="bayesian", prior_alpha=alpha, prior_beta=beta)
            assert [_probability(c) for c in calls] == pytest.approx(_recursive_posteriors(llrs, alpha, beta), abs=1e-9)

    def test_two_reads_have_the_closed_form_of_the_two_state_posterior(self) -> None:
        # n = 2, flat prior: P(z) is proportional to exp(sum LLR over methylated reads) * B(1 + k, 3 - k)
        llrs = [1.5, -0.7]
        calls = call_methylation_from_signal(_site(llrs), model="bayesian")
        # B(1, 3) = 1/3, B(2, 2) = 1/6, B(3, 1) = 1/3
        weights = {(0, 0): 1 / 3, (1, 0): math.exp(1.5) / 6, (0, 1): math.exp(-0.7) / 6, (1, 1): math.exp(0.8) / 3}
        total = sum(weights.values())
        expected = [(weights[(1, 0)] + weights[(1, 1)]) / total, (weights[(0, 1)] + weights[(1, 1)]) / total]
        assert [_probability(c) for c in calls] == pytest.approx(expected, abs=1e-9)

    def test_site_with_thousands_of_reads_stays_exact_and_fast(self) -> None:
        rng = random.Random(2)
        llrs = [rng.gauss(rng.choice([-1.5, 2.0]), 2.5) for _ in range(5000)]
        probabilities = [_probability(c) for c in call_methylation_from_signal(_site(llrs), model="bayesian")]
        assert all(0.0 <= p <= 1.0 for p in probabilities)
        # the posterior of a read is increasing in its own LLR
        ordered = [p for _, p in sorted(zip(llrs, probabilities))]
        assert ordered == sorted(ordered)


class TestPoolingAcrossReads:
    def test_a_weak_dissenting_read_is_pulled_toward_the_site(self) -> None:
        site = _site([4.0] * 20 + [-0.4])
        threshold_calls = call_methylation_from_signal(site, model="threshold")
        bayes_calls = call_methylation_from_signal(site, model="bayesian")
        assert threshold_calls[-1]["methylated"] is False
        assert bayes_calls[-1]["methylated"] is True
        assert all(c["methylated"] for c in bayes_calls)

    def test_the_mirror_case_is_pulled_toward_unmethylated(self) -> None:
        site = _site([-4.0] * 20 + [0.4])
        assert call_methylation_from_signal(site, model="threshold")[-1]["methylated"] is True
        assert call_methylation_from_signal(site, model="bayesian")[-1]["methylated"] is False

    def test_strong_evidence_is_not_overridden_by_the_site(self) -> None:
        # one read with overwhelming evidence stays methylated at a site that is otherwise unmethylated
        site = _site([-4.0] * 20 + [14.0])
        assert call_methylation_from_signal(site, model="bayesian")[-1]["methylated"] is True

    def test_a_mixed_site_keeps_reads_near_their_own_evidence(self) -> None:
        site = _site([3.0, 3.0, 3.0, 3.0, -3.0, -3.0, -3.0, -3.0, 0.2, -0.2])
        bayes = call_methylation_from_signal(site, model="bayesian")
        assert [c["methylated"] for c in bayes] == [True] * 4 + [False] * 4 + [True, False]

    def test_pooling_does_not_depend_on_the_order_of_reads(self) -> None:
        rng = random.Random(8)
        site = _site([rng.gauss(rng.choice([-2.0, 2.0]), 2.0) for _ in range(40)])
        shuffled = site[:]
        rng.shuffle(shuffled)
        by_id = {c["read_id"]: c for c in call_methylation_from_signal(site, model="bayesian")}
        for call in call_methylation_from_signal(shuffled, model="bayesian"):
            assert call["methylated"] == by_id[call["read_id"]]["methylated"]
            assert call["confidence"] == pytest.approx(by_id[call["read_id"]]["confidence"], abs=1e-9)

    def test_sites_are_pooled_by_chromosome_position_and_strand(self) -> None:
        methylated = _site([5.0] * 12, position=100, strand="+")
        unmethylated_other_strand = _site([-5.0] * 12, position=100, strand="-")
        unmethylated_other_position = _site([-5.0] * 12, position=200, strand="+")
        unmethylated_other_chrom = _site([-5.0] * 12, position=100, chrom="chr2", strand="+")
        uninformative = [
            _read("u_plus", 0.0, 100, "chr1", "+"),
            _read("u_minus", 0.0, 100, "chr1", "-"),
            _read("u_pos", 0.0, 200, "chr1", "+"),
            _read("u_chr", 0.0, 100, "chr2", "+"),
        ]
        calls = call_methylation_from_signal(
            methylated
            + unmethylated_other_strand
            + unmethylated_other_position
            + unmethylated_other_chrom
            + uninformative,
            model="bayesian",
        )
        by_id = {c["read_id"]: c["methylated"] for c in calls}
        assert by_id["u_plus"] is True  # joins the methylated site
        assert by_id["u_minus"] is False and by_id["u_pos"] is False and by_id["u_chr"] is False

    def test_a_lone_uninformative_read_stays_at_the_prior(self) -> None:
        (call,) = call_methylation_from_signal([_read("r", 0.0)], model="bayesian")
        assert call["confidence"] == pytest.approx(0.0, abs=1e-12)

    def test_entries_without_a_site_are_not_pooled(self) -> None:
        entries = [{"read_id": f"r{i}", "log_likelihood_ratio": 4.0} for i in range(20)]
        entries.append({"read_id": "odd", "log_likelihood_ratio": -0.4})
        calls = call_methylation_from_signal(entries, model="bayesian")
        assert calls[-1]["methylated"] is False
        assert calls[-1]["confidence"] == pytest.approx(call_methylation_from_signal(entries[-1:])[0]["confidence"])

    def test_pooled_calls_aggregate_like_any_other_calls(self) -> None:
        calls = call_methylation_from_signal(_site([4.0] * 20 + [-0.4]), model="bayesian")
        (site,) = aggregate_methylation(calls, min_coverage=5)["sites"]
        assert site["coverage"] == 21 and site["methylation_freq"] == 1.0


class TestCalibration:
    def test_pooling_beats_the_flat_prior_when_most_reads_are_unmethylated(self) -> None:
        # reads at one site are methylated with probability 0.2; LLR = 2x for a signal x ~ N(+-1, 1)
        rng = random.Random(11)
        truth = [rng.random() < 0.2 for _ in range(3000)]
        llrs = [2.0 * rng.gauss(1.0 if z else -1.0, 1.0) for z in truth]
        site = _site(llrs)
        bayes = [_probability(c) for c in call_methylation_from_signal(site, model="bayesian")]
        flat = [_probability(c) for c in call_methylation_from_signal(site, model="threshold")]

        def brier(probabilities: list[float]) -> float:
            return sum((p - float(z)) ** 2 for p, z in zip(probabilities, truth)) / len(truth)

        assert brier(bayes) < brier(flat) - 0.01
        # the flat 50 % prior over-calls methylation, the learned rate does not
        assert sum(bayes) / len(bayes) == pytest.approx(0.2, abs=0.03)
        assert sum(flat) / len(flat) > 0.28


class TestRobustness:
    def test_nan_llr_is_excluded_from_its_site_and_not_called(self) -> None:
        with_nan = _site([4.0, 4.0, -1.0, float("nan")])
        without = _site([4.0, 4.0, -1.0])
        calls = call_methylation_from_signal(with_nan, model="bayesian")
        reference = call_methylation_from_signal(without, model="bayesian")
        for call, ref in zip(calls[:3], reference):
            assert call["methylated"] == ref["methylated"]
            assert call["confidence"] == pytest.approx(ref["confidence"], abs=1e-12)
        assert calls[3]["methylated"] is False and math.isnan(calls[3]["confidence"])

    def test_extreme_llrs_give_finite_probabilities_without_warnings(self) -> None:
        # warnings are errors under this suite's pytest configuration
        llrs = [float("inf"), -float("inf"), 1e300, -1e300, 800.0, -800.0, 0.0]
        calls = call_methylation_from_signal(_site(llrs), model="bayesian")
        confidences = [c["confidence"] for c in calls]
        assert all(math.isfinite(x) and 0.0 <= x <= 1.0 for x in confidences)
        assert [c["methylated"] for c in calls[:6]] == [True, False, True, False, True, False]

    def test_input_fields_are_kept_and_types_are_plain_python(self) -> None:
        site = _site([2.0, 3.0, -1.0])
        for entry in site:
            entry["extra"] = "kept"
        calls = call_methylation_from_signal(site, model="bayesian")
        assert [c["read_id"] for c in calls] == [e["read_id"] for e in site]
        assert all(c["extra"] == "kept" for c in calls)
        assert all(type(c["methylated"]) is bool and type(c["confidence"]) is float for c in calls)

    def test_empty_input(self) -> None:
        assert call_methylation_from_signal([], model="bayesian") == []

    @pytest.mark.parametrize("bad", [0.0, -1.0, 0.24, float("nan")])
    def test_prior_parameters_below_the_supported_range_raise(self, bad: float) -> None:
        with pytest.raises(ValueError, match="prior_alpha and prior_beta"):
            call_methylation_from_signal(_site([1.0, 2.0]), model="bayesian", prior_alpha=bad)
        with pytest.raises(ValueError, match="prior_alpha and prior_beta"):
            call_methylation_from_signal(_site([1.0, 2.0]), model="bayesian", prior_beta=bad)

    def test_prior_parameters_are_ignored_by_the_threshold_model(self) -> None:
        site = _site([1.0, -2.0])
        assert call_methylation_from_signal(site, model="threshold", prior_alpha=0.0) == call_methylation_from_signal(
            site, model="threshold"
        )
