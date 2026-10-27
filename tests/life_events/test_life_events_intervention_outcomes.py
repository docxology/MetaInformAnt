"""Regression tests: ``intervention_analysis(outcomes=...)`` states what ``outcomes`` must look like.

Defect: ``outcomes`` was split as ``outcomes[:n_pre]`` / ``outcomes[n_pre:]`` without checking that the list holds one
value per pre- and post-intervention sequence. A list of the wrong length failed later, inside ``analyze_life_course``,
with "Number of outcomes must match number of sequences" (which says nothing about the pre/post layout), and a list
with exactly ``n_pre`` values was accepted silently and analysed as pre-intervention outcomes only. The length is now
checked up front and the message names the expected layout and both counts.

Oracle: the pre- and post-intervention means are worked out by hand for the outcomes below.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import pytest

from metainformant.life_events.core.events import Event, EventSequence
from metainformant.life_events.workflow.workflow import intervention_analysis

INTERVENTION = datetime(2017, 1, 1)


def _sequence(person: str, *items: tuple[str, datetime]) -> EventSequence:
    return EventSequence(person, [Event(event_type, moment, "life") for event_type, moment in items])


def _sequences() -> list[EventSequence]:
    return [
        # p1 has events on both sides of the intervention, so it contributes a pre and a post sequence.
        _sequence("p1", ("a", datetime(2010, 1, 1)), ("b", datetime(2020, 1, 1))),
        _sequence("p2", ("a", datetime(2012, 1, 1)), ("c", datetime(2013, 1, 1))),  # pre only
        _sequence("p3", ("b", datetime(2021, 1, 1)), ("a", datetime(2022, 1, 1))),  # post only
    ]


T = INTERVENTION.timestamp()


def test_the_split_has_two_pre_and_two_post_sequences() -> None:
    """Guard the fixture: 2 + 2 = 4 outcomes are expected below."""
    result = intervention_analysis(_sequences(), T)

    assert (result["n_pre_sequences"], result["n_post_sequences"]) == (2, 2)


def test_outcomes_for_every_pre_then_every_post_sequence_are_split_correctly(tmp_path: Path) -> None:
    result = intervention_analysis(_sequences(), T, output_dir=tmp_path, outcomes=[0.2, 0.4, 0.9, 0.7])

    # pre: (0.2 + 0.4) / 2 = 0.3; post: (0.9 + 0.7) / 2 = 0.8; change 0.5
    assert result["outcome_change"]["pre_mean"] == pytest.approx(0.3)
    assert result["outcome_change"]["post_mean"] == pytest.approx(0.8)
    assert result["outcome_change"]["mean"] == pytest.approx(0.5)


def test_outcomes_are_equivalent_to_the_explicit_pre_and_post_arguments(tmp_path: Path) -> None:
    combined = intervention_analysis(_sequences(), T, output_dir=tmp_path / "a", outcomes=[0.2, 0.4, 0.9, 0.7])
    explicit = intervention_analysis(
        _sequences(),
        T,
        output_dir=tmp_path / "b",
        pre_intervention_outcomes=[0.2, 0.4],
        post_intervention_outcomes=[0.9, 0.7],
    )

    assert combined["outcome_change"] == explicit["outcome_change"]


@pytest.mark.parametrize("outcomes", [[0.2, 0.4, 0.9], [0.2, 0.4, 0.9, 0.7, 0.1], []])
def test_a_wrong_number_of_outcomes_is_reported_with_the_expected_layout(tmp_path: Path, outcomes: list) -> None:
    with pytest.raises(ValueError, match=r"2 \+ 2 = 4\), got %d" % len(outcomes)):
        intervention_analysis(_sequences(), T, output_dir=tmp_path, outcomes=outcomes)


def test_the_message_explains_that_pre_comes_before_post(tmp_path: Path) -> None:
    with pytest.raises(ValueError) as excinfo:
        intervention_analysis(_sequences(), T, output_dir=tmp_path, outcomes=[0.2, 0.4, 0.9])

    message = str(excinfo.value)
    assert "pre-intervention sequence followed by one per post-intervention sequence" in message


def test_outcomes_for_the_pre_sequences_only_are_no_longer_accepted_silently(tmp_path: Path) -> None:
    """Two outcomes for two pre-intervention sequences used to be taken as pre outcomes with no post outcomes."""
    with pytest.raises(ValueError, match="got 2"):
        intervention_analysis(_sequences(), T, output_dir=tmp_path, outcomes=[0.2, 0.4])


def test_one_outcome_per_person_is_not_the_expected_layout(tmp_path: Path) -> None:
    """Three people, but four sequences once p1 is split: per-person outcomes cannot be assigned to pre and post."""
    with pytest.raises(ValueError, match=r"2 \+ 2 = 4\), got 3"):
        intervention_analysis(_sequences(), T, output_dir=tmp_path, outcomes=[0.2, 0.4, 0.9])


def test_validation_happens_before_anything_is_analysed_or_written(tmp_path: Path) -> None:
    out = tmp_path / "never_written"

    with pytest.raises(ValueError):
        intervention_analysis(_sequences(), T, output_dir=out, outcomes=[1.0])

    assert not (out / "pre").exists()
    assert not (out / "post").exists()


def test_with_no_post_intervention_sequences_one_outcome_per_pre_sequence_is_enough(tmp_path: Path) -> None:
    pre_only = [
        _sequence("q1", ("a", datetime(2010, 1, 1)), ("b", datetime(2011, 1, 1))),
        _sequence("q2", ("a", datetime(2012, 1, 1)), ("c", datetime(2013, 1, 1))),
    ]

    result = intervention_analysis(pre_only, T, output_dir=tmp_path, outcomes=[0, 1])

    assert (result["n_pre_sequences"], result["n_post_sequences"]) == (2, 0)
    assert result["post_intervention"] == {}
    assert result["outcome_change"] == {}  # nothing to compare against


def test_explicit_pre_and_post_outcomes_take_precedence_over_outcomes(tmp_path: Path) -> None:
    """``outcomes`` is documented as ignored when an explicit pre/post argument is given, so it is not validated."""
    result = intervention_analysis(
        _sequences(),
        T,
        output_dir=tmp_path,
        outcomes=[999],
        pre_intervention_outcomes=[0.2, 0.4],
        post_intervention_outcomes=[0.9, 0.7],
    )

    assert result["outcome_change"]["mean"] == pytest.approx(0.5)


def test_no_outcomes_still_works(tmp_path: Path) -> None:
    result = intervention_analysis(_sequences(), T, output_dir=tmp_path)

    assert result["outcome_change"] == {}
    assert result["n_pre_sequences"] == 2
