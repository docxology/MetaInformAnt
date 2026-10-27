"""Regression tests: ``phenotype.analysis.life_course.analyze_life_course`` performs a real outcome analysis.

Defect: the outcome branch was marked "Outcome analysis (placeholder)" and "This would integrate with predictive
models". It counted the sequences whose ``metadata["outcomes"]`` contained the outcome name and divided by *all*
sequences, so a sequence with no recorded outcomes (status unknown) was silently counted as a negative, and when
no sequence recorded outcomes at all the result was a fabricated prevalence of 0.0. Nothing related the outcome to
the trajectories.

The analysis now (a) uses only sequences with a recorded outcome status, reports how many were left out and gives
no prevalence when there are none, and (b) compares the trajectories of people with and without the outcome: three
trajectory features with a Mann-Whitney U test and every event type with a Fisher exact test.

Oracles, all independent of the module and of scipy: the dataset below is small enough to count by hand (the tables
are worked out in comments); the Fisher p-values are re-derived with exact rational arithmetic from the
hypergeometric distribution; the Mann-Whitney p-values with tied data are re-derived from the textbook
tie-corrected normal approximation with continuity correction, and the exact p-value for a tie-free tiny sample by
enumerating every assignment of the ranks.
"""

from __future__ import annotations

import itertools
import json
import math
from collections import Counter
from fractions import Fraction
from typing import Dict, List, Optional, Tuple

import pytest

from metainformant.core.utils import errors
from metainformant.phenotype.analysis.life_course import Event, EventSequence, analyze_life_course

OUTCOME = "relapse"


def _sequence(person: str, events: List[Tuple[float, str]], outcomes: Optional[object] = None) -> EventSequence:
    metadata = {} if outcomes is None else {"outcomes": outcomes}
    return EventSequence(person, [Event(timestamp=t, event_type=e) for t, e in events], metadata=metadata)


# Four people with the outcome, six without, two whose status is unknown (no "outcomes" entry at all).
WITH_OUTCOME = [
    _sequence("P1", [(0, "a"), (10, "b"), (20, "c")], [OUTCOME]),
    _sequence("P2", [(0, "a"), (5, "a"), (15, "c"), (40, "d")], [OUTCOME, "other"]),
    _sequence("P3", [(0, "a"), (30, "b")], [OUTCOME]),
    _sequence("P4", [(0, "a"), (10, "c"), (20, "c"), (30, "c"), (50, "d")], [OUTCOME]),
]
WITHOUT_OUTCOME = [
    _sequence("P5", [(0, "b"), (5, "c")], []),
    _sequence("P6", [(0, "b")], ["other"]),
    _sequence("P7", [(0, "b"), (2, "b"), (4, "c")], []),
    _sequence("P8", [(0, "c")], []),
    _sequence("P9", [(0, "b"), (6, "c"), (9, "d")], []),
    _sequence("P10", [(0, "b"), (8, "b")], []),
]
UNKNOWN_STATUS = [
    _sequence("U1", [(0, "a"), (7, "d")]),
    _sequence("U2", [(0, "c")]),
]
ALL_SEQUENCES = WITH_OUTCOME + WITHOUT_OUTCOME + UNKNOWN_STATUS


def _analysis(sequences: Optional[List[EventSequence]] = None, outcome: str = OUTCOME) -> Dict:
    return analyze_life_course(sequences or ALL_SEQUENCES, outcomes=[outcome])["outcome_analysis"][outcome]


# ---------------------------------------------------------------------------
# Independent oracles
# ---------------------------------------------------------------------------


def fisher_exact_two_sided(a: int, b: int, c: int, d: int) -> Fraction:
    """Sum of the hypergeometric probabilities that do not exceed that of the observed table (exact)."""
    row1, row2, col1 = a + b, c + d, a + c
    total = math.comb(row1 + row2, col1)

    def probability(x: int) -> Fraction:
        return Fraction(math.comb(row1, x) * math.comb(row2, col1 - x), total)

    observed = probability(a)
    support = range(max(0, col1 - row2), min(row1, col1) + 1)
    return sum((probability(x) for x in support if probability(x) <= observed), Fraction(0))


def average_ranks(values: List[float]) -> List[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        for k in range(i, j + 1):
            ranks[order[k]] = (i + j) / 2 + 1  # average of the 1-based ranks i + 1 .. j + 1
        i = j + 1
    return ranks


def mann_whitney_asymptotic_p(x: List[float], y: List[float]) -> float:
    """Two-sided p of the tie-corrected normal approximation with continuity correction."""
    n1, n2 = len(x), len(y)
    n = n1 + n2
    u1 = sum(average_ranks(x + y)[:n1]) - n1 * (n1 + 1) / 2
    mean = n1 * n2 / 2
    tie_term = sum(t**3 - t for t in Counter(x + y).values())
    sigma = math.sqrt(n1 * n2 / 12 * ((n + 1) - tie_term / (n * (n - 1))))
    z = (max(u1, n1 * n2 - u1) - mean - 0.5) / sigma
    return min(1.0, math.erfc(z / math.sqrt(2)))


def mann_whitney_exact_p(x: List[float], y: List[float]) -> Fraction:
    """Two-sided exact p for tie-free data: the share of rank assignments at least as far from the mean U."""
    n1, n = len(x), len(x) + len(y)
    ranks = average_ranks(x + y)
    observed = sum(ranks[:n1]) - n1 * (n1 + 1) / 2
    mean = Fraction(n1 * (n - n1), 2)
    extreme = total = 0
    for chosen in itertools.combinations(range(1, n + 1), n1):
        u = sum(chosen) - n1 * (n1 + 1) / 2
        total += 1
        extreme += abs(u - mean) >= abs(observed - mean)
    return Fraction(extreme, total)


# ---------------------------------------------------------------------------
# Prevalence, and what counts as "no outcome"
# ---------------------------------------------------------------------------


def test_counts_and_prevalence_use_only_sequences_with_a_recorded_status() -> None:
    analysis = _analysis()

    assert analysis["status"] == "ok"
    assert analysis["n_sequences"] == 12
    assert analysis["n_with_outcome_data"] == 10
    assert analysis["n_missing_outcome_data"] == 2
    assert analysis["n_with_outcome"] == 4
    assert analysis["n_without_outcome"] == 6
    # 4 of the 10 people whose outcomes are known; counting the two unknowns as negatives would give 4 / 12.
    assert analysis["prevalence"] == pytest.approx(0.4)
    assert analysis["prevalence"] != pytest.approx(4 / 12)


def test_no_recorded_outcomes_at_all_gives_no_prevalence_instead_of_zero() -> None:
    analysis = _analysis(UNKNOWN_STATUS)

    assert analysis["status"] == "no_outcome_data"
    assert analysis["prevalence"] is None
    assert (analysis["n_with_outcome_data"], analysis["n_missing_outcome_data"]) == (0, 2)
    assert (analysis["n_with_outcome"], analysis["n_without_outcome"]) == (0, 0)
    assert analysis["trajectory_comparison"] == {}
    assert analysis["event_type_association"] == {}


def test_an_outcome_nobody_with_recorded_data_has_is_a_real_zero() -> None:
    analysis = _analysis(outcome="never_recorded")

    assert analysis["status"] == "ok"
    assert (analysis["n_with_outcome"], analysis["n_without_outcome"]) == (0, 10)
    assert analysis["prevalence"] == 0.0
    assert analysis["trajectory_comparison"] == {}  # nobody has it: there is no group to compare with
    assert analysis["event_type_association"] == {}


def test_an_outcome_everybody_has_has_prevalence_one_and_no_comparison() -> None:
    sequences = [_sequence(f"q{i}", [(0, "a"), (1, "b")], [OUTCOME]) for i in range(4)]

    analysis = _analysis(sequences)

    assert analysis["prevalence"] == 1.0
    assert analysis["trajectory_comparison"] == {}
    assert analysis["event_type_association"] == {}


def test_a_single_outcome_name_stored_as_text_is_one_outcome_not_its_letters() -> None:
    sequences = [
        _sequence("t1", [(0, "a"), (4, "b")], OUTCOME),  # "relapse", not ["relapse"]
        _sequence("t2", [(0, "a")], "recovery"),
    ]

    assert _analysis(sequences, OUTCOME)["n_with_outcome"] == 1
    assert _analysis(sequences, "re")["n_with_outcome"] == 0  # a substring of "relapse" is not an outcome


def test_malformed_outcome_metadata_is_rejected() -> None:
    sequences = [_sequence("bad", [(0, "a"), (1, "b")], 5)]

    with pytest.raises(errors.ValidationError, match="collection of outcome names"):
        analyze_life_course(sequences, outcomes=[OUTCOME])


def test_several_outcomes_are_analysed_independently() -> None:
    result = analyze_life_course(ALL_SEQUENCES, outcomes=[OUTCOME, "other", "never_recorded"])["outcome_analysis"]

    assert set(result) == {OUTCOME, "other", "never_recorded"}
    assert result[OUTCOME]["n_with_outcome"] == 4
    assert result["other"]["n_with_outcome"] == 2  # P2 and P6
    assert result["other"]["prevalence"] == pytest.approx(0.2)
    assert result["never_recorded"]["n_with_outcome"] == 0


def test_without_outcomes_there_is_no_outcome_analysis() -> None:
    assert "outcome_analysis" not in analyze_life_course(ALL_SEQUENCES)


# ---------------------------------------------------------------------------
# Trajectory comparison (hand-counted features; independent test oracles)
# ---------------------------------------------------------------------------

# n_events   with outcome [3, 4, 2, 5] (mean 3.5)           without [2, 1, 3, 1, 3, 2] (mean 2.0)
# duration   with         [20, 40, 30, 50] (mean 35)        without [5, 0, 4, 0, 9, 8] (mean 26/6)
# complexity with  [1, 3/4, 1, 3/5] (mean 0.8375)           without [1, 1, 2/3, 1, 1, 1/2] (mean 31/36)
FEATURES = {
    "n_events": ([3.0, 4.0, 2.0, 5.0], [2.0, 1.0, 3.0, 1.0, 3.0, 2.0]),
    "duration": ([20.0, 40.0, 30.0, 50.0], [5.0, 0.0, 4.0, 0.0, 9.0, 8.0]),
    "complexity": ([1.0, 0.75, 1.0, 0.6], [1.0, 1.0, 2 / 3, 1.0, 1.0, 0.5]),
}


@pytest.mark.parametrize("feature", list(FEATURES))
def test_trajectory_features_are_compared_between_the_groups(feature: str) -> None:
    with_values, without_values = FEATURES[feature]

    comparison = _analysis()["trajectory_comparison"][feature]

    assert comparison["n_with_outcome"] == 4
    assert comparison["n_without_outcome"] == 6
    assert comparison["mean_with_outcome"] == pytest.approx(sum(with_values) / 4)
    assert comparison["mean_without_outcome"] == pytest.approx(sum(without_values) / 6)
    assert comparison["difference"] == pytest.approx(sum(with_values) / 4 - sum(without_values) / 6)
    assert comparison["p_value"] == pytest.approx(mann_whitney_asymptotic_p(with_values, without_values), rel=1e-6)
    assert 0.0 <= comparison["p_value"] <= 1.0


def test_the_hand_counted_means_match_the_comments_above() -> None:
    comparison = _analysis()["trajectory_comparison"]

    assert comparison["n_events"]["mean_with_outcome"] == 3.5
    assert comparison["n_events"]["mean_without_outcome"] == 2.0
    assert comparison["duration"]["mean_with_outcome"] == 35.0
    assert comparison["duration"]["mean_without_outcome"] == pytest.approx(26 / 6)
    assert comparison["complexity"]["mean_with_outcome"] == pytest.approx(0.8375)
    assert comparison["complexity"]["mean_without_outcome"] == pytest.approx(31 / 36)


def test_a_tie_free_small_sample_gets_the_exact_rank_sum_p_value() -> None:
    """With 2 + 3 people and no tied values scipy tests exactly; the oracle enumerates all ten rank assignments."""
    sequences = [
        _sequence("e1", [(i, "a") for i in range(2)], [OUTCOME]),
        _sequence("e2", [(i, "a") for i in range(5)], [OUTCOME]),
        _sequence("e3", [(i, "a") for i in range(3)], []),
        _sequence("e4", [(i, "a") for i in range(7)], []),
        _sequence("e5", [(i, "a") for i in range(4)], []),
    ]

    p_value = _analysis(sequences)["trajectory_comparison"]["n_events"]["p_value"]

    # Event counts 2, 5 | 3, 7, 4: U = 2 of a mean 3, and 8 of the 10 assignments are at least that far out.
    assert mann_whitney_exact_p([2.0, 5.0], [3.0, 7.0, 4.0]) == Fraction(4, 5)
    assert p_value == pytest.approx(0.8)


def test_identical_feature_values_give_p_one_instead_of_nan() -> None:
    sequences = [
        _sequence(f"s{i}", [(0, "a"), (2, "b")], [OUTCOME] if i < 3 else []) for i in range(6)
    ]  # every sequence has two events, a duration of 2 and a complexity of 1

    comparison = _analysis(sequences)["trajectory_comparison"]

    for feature in ("n_events", "duration", "complexity"):
        assert comparison[feature]["difference"] == 0.0
        assert comparison[feature]["p_value"] == 1.0


def test_duration_and_complexity_ignore_sequences_without_events() -> None:
    """A sequence with no events has no duration or complexity; it still counts as zero events."""
    sequences = [
        _sequence("w1", [(0, "a"), (10, "b")], [OUTCOME]),
        _sequence("w2", [], [OUTCOME]),
        _sequence("n1", [(0, "a"), (4, "a")], []),
        _sequence("n2", [(0, "b")], []),
    ]

    comparison = _analysis(sequences)["trajectory_comparison"]

    assert comparison["n_events"]["n_with_outcome"] == 2
    assert comparison["n_events"]["mean_with_outcome"] == 1.0  # (2 + 0) / 2
    assert comparison["duration"]["n_with_outcome"] == 1  # only w1 has a duration (10)
    assert comparison["duration"]["mean_with_outcome"] == 10.0
    assert comparison["complexity"]["n_with_outcome"] == 1


def test_a_feature_no_sequence_of_a_group_defines_is_none() -> None:
    sequences = [
        _sequence("w1", [], [OUTCOME]),
        _sequence("w2", [], [OUTCOME]),
        _sequence("n1", [(0, "a"), (4, "b")], []),
        _sequence("n2", [(0, "b")], []),
    ]

    comparison = _analysis(sequences)["trajectory_comparison"]

    assert comparison["n_events"] is not None
    assert comparison["duration"] is None
    assert comparison["complexity"] is None


# ---------------------------------------------------------------------------
# Event-type association (hand-counted 2x2 tables; exact Fisher oracle)
# ---------------------------------------------------------------------------

# Table per event type: [[exposed & outcome, exposed & no outcome], [unexposed & outcome, unexposed & no outcome]]
#   a [[4, 0], [0, 6]]  p = 1/210       b [[2, 5], [2, 1]]  OR 0.2, p = 1/2
#   c [[3, 4], [1, 2]]  OR 1.5, p = 1   d [[2, 1], [2, 5]]  OR 5,   p = 1/2
TABLES = {"a": (4, 0, 0, 6), "b": (2, 5, 2, 1), "c": (3, 4, 1, 2), "d": (2, 1, 2, 5)}


def test_the_hand_counted_tables_have_the_stated_exact_p_values() -> None:
    assert fisher_exact_two_sided(*TABLES["a"]) == Fraction(1, 210)
    assert fisher_exact_two_sided(*TABLES["b"]) == Fraction(1, 2)
    assert fisher_exact_two_sided(*TABLES["c"]) == Fraction(1)
    assert fisher_exact_two_sided(*TABLES["d"]) == Fraction(1, 2)


def test_every_event_type_is_tested_with_the_exact_fisher_p_value() -> None:
    association = _analysis()["event_type_association"]

    assert set(association) == {"a", "b", "c", "d"}
    for event_type, (exposed_with, exposed_without, unexposed_with, unexposed_without) in TABLES.items():
        entry = association[event_type]
        assert entry["n_exposed"] == exposed_with + exposed_without
        assert entry["n_exposed_with_outcome"] == exposed_with
        assert entry["n_unexposed_with_outcome"] == unexposed_with
        assert entry["p_value"] == pytest.approx(float(fisher_exact_two_sided(*TABLES[event_type])), rel=1e-9)


def test_odds_ratios_are_the_sample_odds_ratios_and_none_when_a_cell_is_zero() -> None:
    association = _analysis()["event_type_association"]

    assert association["a"]["odds_ratio"] is None  # no exposed person without the outcome: infinite
    assert association["b"]["odds_ratio"] == pytest.approx(0.2)  # (2 * 1) / (5 * 2)
    assert association["c"]["odds_ratio"] == pytest.approx(1.5)  # (3 * 2) / (4 * 1)
    assert association["d"]["odds_ratio"] == pytest.approx(5.0)  # (2 * 5) / (1 * 2)


def test_event_types_are_ordered_by_p_value_with_ties_broken_by_name() -> None:
    association = _analysis()["event_type_association"]

    assert list(association) == ["a", "b", "d", "c"]  # 1/210, then b and d tied at 1/2, then c at 1


def test_an_event_type_everyone_or_no_one_has_is_a_null_table() -> None:
    sequences = [
        _sequence("m1", [(0, "common"), (1, "rare")], [OUTCOME]),
        _sequence("m2", [(0, "common")], [OUTCOME]),
        _sequence("m3", [(0, "common")], []),
        _sequence("m4", [(0, "common")], []),
    ]

    association = _analysis(sequences)["event_type_association"]

    assert association["common"]["p_value"] == 1.0
    assert association["common"]["odds_ratio"] is None  # nobody is unexposed
    assert association["rare"]["n_exposed"] == 1


def test_the_whole_analysis_is_plain_json() -> None:
    """No NaN, no infinity, no numpy scalars: the analysis can be written out as it is."""
    result = analyze_life_course(ALL_SEQUENCES, outcomes=[OUTCOME, "other", "never_recorded"])

    json.dumps(result["outcome_analysis"], allow_nan=False)


def test_the_rest_of_the_life_course_analysis_is_unchanged() -> None:
    result = analyze_life_course(ALL_SEQUENCES, outcomes=[OUTCOME])

    assert result["n_sequences"] == 12
    assert result["event_frequency"]["a"] == 1 + 2 + 1 + 1 + 1  # P1, P2 (twice), P3, P4, U1
    assert "trajectory_patterns" in result and "duration_stats" in result
