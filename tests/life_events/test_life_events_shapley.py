"""Regression tests for ``shapley_values`` and for the tie-breaking hack in ``event_importance``.

``shapley_values`` used to return ``np.random.normal(0, 0.1)`` for every event. The tests below check
the defining properties of Shapley values (efficiency, null player, symmetry, additivity) against
games whose exact answer is known, the seeded Monte Carlo path used for long sequences, and the real
``EventSequencePredictor``. ``event_importance`` used to add ``norm(embedding) * 1e-6`` to every
score, which made tied events rank by embedding norm.
"""

from __future__ import annotations

import itertools
from datetime import datetime
from typing import Any, Dict, List, Sequence

import numpy as np
import pytest

from metainformant.life_events.analysis.interpretability import event_importance, shapley_values
from metainformant.life_events.core.events import Event, EventSequence
from metainformant.life_events.models.predictor import EventSequencePredictor


class TablePredictor:
    """A deterministic additive model: the score of a sequence is the sum of per-event weights."""

    is_fitted = True

    def __init__(self, weights: Dict[str, float], offset: float = 0.0):
        self.weights = dict(weights)
        self.offset = offset
        self.evaluated = 0

    def predict(self, sequences: Sequence[Sequence[str]]) -> np.ndarray:
        self.evaluated += len(sequences)
        return np.array([self.offset + sum(self.weights.get(token, 0.0) for token in seq) for seq in sequences])


class ConjunctionPredictor:
    """Scores 1 when every required event is present, in the given order, else 0 (a pure interaction)."""

    is_fitted = True

    def __init__(self, required: List[str]):
        self.required = required

    def _score(self, sequence: Sequence[str]) -> float:
        positions = [list(sequence).index(token) for token in self.required if token in sequence]
        ordered = positions == sorted(positions)
        return float(len(positions) == len(self.required) and ordered)

    def predict(self, sequences: Sequence[Sequence[str]]) -> np.ndarray:
        return np.array([self._score(seq) for seq in sequences])


class SubsetGame:
    """A game given by an arbitrary table: v(S) is a fixed pseudo-random number for every subset of events."""

    is_fitted = True

    def __init__(self, tokens: List[str], seed: int, scale: float = 1.0):
        rng = np.random.default_rng(seed)
        self.table = {
            frozenset(subset): float(scale * rng.normal())
            for size in range(len(tokens) + 1)
            for subset in itertools.combinations(tokens, size)
        }

    def predict(self, sequences: Sequence[Sequence[str]]) -> np.ndarray:
        return np.array([self.table[frozenset(seq)] for seq in sequences])


def shapley_by_definition(tokens: List[str], game: SubsetGame) -> Dict[str, float]:
    """Average marginal contribution of every event over all n! orders in which the events can join."""
    contributions = dict.fromkeys(tokens, 0.0)
    orders = list(itertools.permutations(tokens))
    for order in orders:
        present: set = set()
        for token in order:
            before = game.table[frozenset(present)]
            present.add(token)
            contributions[token] += game.table[frozenset(present)] - before
    return {token: value / len(orders) for token, value in contributions.items()}


def total(values: Dict[str, float]) -> float:
    return float(sum(values.values()))


def test_additive_model_gives_each_event_its_own_weight():
    weights = {"a": 1.0, "b": -2.0, "c": 0.5}
    predictor = TablePredictor(weights)

    values = shapley_values(predictor, ["a", "b", "c"], [], n_samples=10)

    assert values == {key: pytest.approx(value, abs=1e-12) for key, value in weights.items()}


def test_efficiency_with_exact_enumeration():
    """The values sum to v(all events) - v(empty sequence), including a non-zero empty baseline."""
    predictor = TablePredictor({"a": 1.5, "b": -0.25, "c": 4.0, "d": 0.75}, offset=3.0)
    sequence = ["a", "b", "c", "d", "a"]

    values = shapley_values(predictor, sequence, [])

    full = predictor.predict([sequence])[0]
    empty = predictor.predict([[]])[0]
    assert empty == 3.0
    assert total(values) == pytest.approx(full - empty, abs=1e-12)


def test_efficiency_holds_for_interactions():
    predictor = ConjunctionPredictor(["a", "b", "c"])
    sequence = ["a", "x", "b", "y", "c"]

    values = shapley_values(predictor, sequence, None)

    assert total(values) == pytest.approx(1.0, abs=1e-12)
    for token in ("a", "b", "c"):
        assert values[token] == pytest.approx(1 / 3, abs=1e-12)


def test_interaction_partners_share_the_credit_equally():
    values = shapley_values(ConjunctionPredictor(["A", "B"]), ["A", "B", "C", "D"], [])

    assert values["A"] == pytest.approx(0.5, abs=1e-12)
    assert values["B"] == pytest.approx(0.5, abs=1e-12)


def test_null_player_gets_exactly_zero():
    predictor = TablePredictor({"a": 2.0, "b": 3.0})

    values = shapley_values(predictor, ["a", "b", "irrelevant", "also_irrelevant"], [])

    assert values["irrelevant"] == 0.0
    assert values["also_irrelevant"] == 0.0


def test_symmetric_events_get_equal_values():
    predictor = TablePredictor({"a": 2.0, "b": 2.0, "c": 7.0})

    values = shapley_values(predictor, ["a", "c", "b"], [])

    assert values["a"] == pytest.approx(values["b"], abs=1e-12)
    assert values["c"] > values["a"]


def test_repeated_events_are_reported_once_with_their_summed_value():
    predictor = TablePredictor({"a": 1.0, "b": 2.0})

    values = shapley_values(predictor, ["a", "a", "b"], [])

    assert values == {"a": pytest.approx(2.0), "b": pytest.approx(2.0)}


def test_the_order_of_the_remaining_events_is_preserved_in_every_coalition():
    predictor = ConjunctionPredictor(["A", "B"])

    in_order = shapley_values(predictor, ["A", "B", "C"], [])
    reversed_order = shapley_values(predictor, ["B", "A", "C"], [])

    assert in_order["A"] == pytest.approx(0.5)
    assert in_order["B"] == pytest.approx(0.5)
    assert all(value == 0.0 for value in reversed_order.values())


def test_exact_enumeration_is_deterministic_and_does_not_depend_on_the_seed():
    predictor = ConjunctionPredictor(["A", "B"])
    sequence = ["A", "B", "C"]

    first = shapley_values(predictor, sequence, [], random_state=1)
    second = shapley_values(predictor, sequence, [], random_state=2)

    assert first == second


def test_exact_enumeration_evaluates_every_coalition_once():
    predictor = TablePredictor({"a": 1.0})

    shapley_values(predictor, ["a", "b", "c", "d", "e"], [])

    assert predictor.evaluated == 2**5


def test_long_sequences_are_estimated_from_sampled_permutations():
    tokens = ["A", "B"] + [f"e{i}" for i in range(10)]  # 12 events: above the exact-enumeration limit
    predictor = ConjunctionPredictor(["A", "B"])

    values = shapley_values(predictor, tokens, [], n_samples=1500, random_state=0)

    assert total(values) == pytest.approx(1.0, abs=1e-12)  # efficiency is exact even when sampling
    assert values["A"] == pytest.approx(0.5, abs=0.06)
    assert values["B"] == pytest.approx(0.5, abs=0.06)
    assert all(values[f"e{i}"] == 0.0 for i in range(10))


def test_sampling_is_seeded():
    tokens = ["A", "B"] + [f"e{i}" for i in range(10)]
    predictor = ConjunctionPredictor(["A", "B"])

    first = shapley_values(predictor, tokens, [], n_samples=40, random_state=3)
    again = shapley_values(predictor, tokens, [], n_samples=40, random_state=3)
    other = shapley_values(predictor, tokens, [], n_samples=40, random_state=4)

    assert first == again
    assert first != other


def test_sampled_additive_model_is_exact_and_uses_far_fewer_evaluations_than_enumeration():
    weights = {f"e{i}": float(i) for i in range(12)}
    predictor = TablePredictor(weights)

    values = shapley_values(predictor, list(weights), [], n_samples=5, random_state=0)

    for token, weight in weights.items():
        assert values[token] == pytest.approx(weight, abs=1e-9)
    assert predictor.evaluated < 2 + 5 * 12  # versus 2**12 = 4096 coalitions


def test_max_exact_events_switches_between_enumeration_and_sampling():
    tokens = ["a", "b", "c", "d", "e", "f"]
    enumerated = TablePredictor({"a": 1.0})
    sampled = TablePredictor({"a": 1.0})

    shapley_values(enumerated, tokens, [], n_samples=3, max_exact_events=6)
    values = shapley_values(sampled, tokens, [], n_samples=3, max_exact_events=3, random_state=0)

    assert enumerated.evaluated == 2**6
    assert sampled.evaluated < 2**6
    assert total(values) == pytest.approx(1.0, abs=1e-12)


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_exact_values_match_the_definition_on_arbitrary_games(seed):
    tokens = [f"e{i}" for i in range(5)]
    game = SubsetGame(tokens, seed)

    values = shapley_values(game, tokens, [])

    expected = shapley_by_definition(tokens, game)
    for token in tokens:
        assert values[token] == pytest.approx(expected[token], abs=1e-12)
    assert total(values) == pytest.approx(game.table[frozenset(tokens)] - game.table[frozenset()], abs=1e-12)


def test_sampling_agrees_with_exact_enumeration_on_an_arbitrary_game():
    tokens = [f"e{i}" for i in range(7)]
    game = SubsetGame(tokens, seed=3, scale=0.25)

    exact = shapley_values(game, tokens, [])
    sampled = shapley_values(game, tokens, [], n_samples=4000, max_exact_events=3, random_state=0)

    for token in tokens:
        assert sampled[token] == pytest.approx(exact[token], abs=0.06)
    assert total(sampled) == pytest.approx(game.table[frozenset(tokens)] - game.table[frozenset()], abs=1e-12)


def test_background_sequences_do_not_change_the_values():
    predictor = TablePredictor({"a": 1.0, "b": 2.0})

    without = shapley_values(predictor, ["a", "b"], None)
    with_background = shapley_values(predictor, ["a", "b"], [["b"], ["a", "a"], []])

    assert without == with_background


def test_event_sequence_objects_are_explained_by_their_tokens():
    sequence = EventSequence(
        "p1",
        [Event("degree", datetime(2010, 1, 1), "education"), Event("diagnosis", datetime(2015, 1, 1), "health")],
    )
    predictor = TablePredictor({"education:degree": 1.0, "health:diagnosis": -3.0})

    values = shapley_values(predictor, sequence, [])

    assert values == {"education:degree": pytest.approx(1.0), "health:diagnosis": pytest.approx(-3.0)}


def test_empty_sequence_has_no_values():
    assert shapley_values(TablePredictor({"a": 1.0}), [], []) == {}


def test_values_for_a_fitted_event_sequence_predictor_satisfy_efficiency():
    sequences = [["a", "b"], ["a", "c"], ["b", "c"], ["a", "b", "c"], ["c"], ["a"], ["b"], ["c", "c"]]
    outcomes = np.array([1.0, 2.0, 3.5, 4.0, 0.5, 1.5, 2.5, 1.0])
    predictor = EventSequencePredictor(model_type="embedding", task_type="regression", embedding_dim=8, random_state=0)
    predictor.fit(sequences, outcomes)

    target = ["a", "b", "c", "a"]
    values = shapley_values(predictor, target, [])

    assert set(values) == {"a", "b", "c"}
    expected = float(predictor.predict([target])[0] - predictor.predict([[]])[0])
    assert total(values) == pytest.approx(expected, abs=1e-9)
    assert shapley_values(predictor, target, []) == values  # no randomness left


def test_values_for_a_simple_length_model_satisfy_efficiency():
    sequences = [["a"] * n for n in (1, 2, 3, 4)] * 2
    outcomes = np.array([1.0, 2.0, 4.0, 8.0] * 2)
    predictor = EventSequencePredictor(model_type="simple", task_type="regression")
    predictor.fit(sequences, outcomes)

    values = shapley_values(predictor, ["a", "b", "c"], [])

    expected = float(predictor.predict([["a", "b", "c"]])[0] - predictor.predict([[]])[0])
    assert total(values) == pytest.approx(expected, abs=1e-9)


@pytest.mark.parametrize(
    "arguments, error, message",
    [
        ({"n_samples": 0}, ValueError, "n_samples"),
        ({"max_exact_events": -1}, ValueError, "max_exact_events"),
        ({"max_exact_events": 21}, ValueError, "max_exact_events"),
    ],
)
def test_invalid_options_are_rejected(arguments, error, message):
    with pytest.raises(error, match=message):
        shapley_values(TablePredictor({"a": 1.0}), ["a", "b"], [], **arguments)


def test_unfitted_predictors_are_rejected():
    with pytest.raises(ValueError, match="must be fitted"):
        shapley_values(EventSequencePredictor(), ["a", "b"], [])


def test_objects_without_predict_are_rejected():
    with pytest.raises(TypeError, match="predict"):
        shapley_values(object(), ["a", "b"], [])


def test_non_numeric_predictions_are_rejected():
    class LabelPredictor:
        is_fitted = True

        def predict(self, sequences: Sequence[Any]) -> List[str]:
            return ["low" for _ in sequences]

    with pytest.raises(ValueError, match="numeric"):
        shapley_values(LabelPredictor(), ["a", "b"], [])


# --- event_importance: ties must stay ties -------------------------------------------------------


def test_event_importance_ties_are_not_broken_by_embedding_norm():
    sequences = [["a", "b", "c"], ["b", "c", "d"], ["a", "d", "c"], ["d", "b", "a"]]
    embeddings = {token: np.full(8, scale) for token, scale in zip("abcd", (1.0, 2.0, 3.0, 4.0))}
    # A single outcome class: the predictor never changes its answer when an event is removed.
    predictor = EventSequencePredictor(
        model_type="embedding", task_type="classification", embedding_dim=8, random_state=0
    )
    predictor.fit(sequences, np.array([1, 1, 1, 1]), event_embeddings=embeddings)

    importance = event_importance(predictor, sequences, embeddings, method="permutation")

    assert set(importance) == {"a", "b", "c", "d"}
    assert set(importance.values()) == {0.0}


def test_event_importance_follows_prediction_changes_not_embedding_norms():
    sequences = [["a", "b", "c"], ["b", "c", "a"], ["c", "a", "b"]]
    predictor = TablePredictor({"a": 3.0, "b": 1.0})
    # "c" never changes the prediction but has by far the largest embedding.
    embeddings = {"a": np.ones(4), "b": np.ones(4) * 2, "c": np.ones(4) * 1e7}

    importance = event_importance(predictor, sequences, embeddings, method="permutation")

    assert importance["a"] == pytest.approx(1.0)
    assert importance["b"] == pytest.approx(1 / 3)
    assert importance["c"] == 0.0
