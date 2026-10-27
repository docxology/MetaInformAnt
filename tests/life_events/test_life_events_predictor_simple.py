"""Regression tests for ``EventSequencePredictor(model_type="simple")`` and model persistence.

Defects covered:

* Classification pushed the length model's mean outcomes (already probabilities) through a sigmoid,
  so every prediction was the positive class and ``predict_proba`` was a constant 0.5.
* ``save_model``/``load_model`` dropped ``sequence_length_model``: a reloaded simple model predicted
  0.0 for everything.
* ``model_type="lstm"`` was documented but always raised ``ValueError``.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pytest

from metainformant.life_events.models.predictor import EventSequencePredictor


def sequences_of_lengths(lengths: List[int]) -> List[List[str]]:
    return [[f"event_{i}" for i in range(n)] for n in lengths]


def length_task() -> Tuple[List[List[str]], np.ndarray]:
    """Outcome 1 exactly when the sequence has at least five events: learnable from the length alone."""
    lengths = [1, 2, 3, 4, 6, 7, 8, 9] * 3
    return sequences_of_lengths(lengths), np.array([int(n >= 5) for n in lengths])


def test_simple_classifier_learns_a_task_that_depends_on_sequence_length():
    sequences, outcomes = length_task()
    predictor = EventSequencePredictor(model_type="simple", task_type="classification").fit(sequences, outcomes)

    predictions = predictor.predict(sequences)

    assert set(predictions.tolist()) == {0, 1}  # not a constant
    np.testing.assert_array_equal(predictions, outcomes)


def test_simple_classifier_probabilities_are_the_positive_rates_per_length():
    # Two-event sequences are positive 1 time in 4; four-event sequences 3 times in 4.
    lengths = [2, 2, 2, 2, 4, 4, 4, 4]
    outcomes = np.array([0, 0, 0, 1, 1, 1, 1, 0])
    predictor = EventSequencePredictor(model_type="simple", task_type="classification")
    predictor.fit(sequences_of_lengths(lengths), outcomes)

    probabilities = predictor.predict_proba(sequences_of_lengths([2, 4]))

    np.testing.assert_allclose(probabilities, [[0.75, 0.25], [0.25, 0.75]])
    np.testing.assert_array_equal(predictor.predict(sequences_of_lengths([2, 4])), [0, 1])


def test_simple_classifier_threshold_is_one_half():
    # Three-event sequences are positive exactly half of the time: ties go to the positive class.
    sequences = sequences_of_lengths([3, 3, 1, 1])
    predictor = EventSequencePredictor(model_type="simple", task_type="classification")
    predictor.fit(sequences, np.array([0, 1, 0, 0]))

    assert predictor.predict_proba(sequences_of_lengths([3]))[0, 1] == pytest.approx(0.5)
    assert predictor.predict(sequences_of_lengths([3, 1])).tolist() == [1, 0]


def test_simple_classifier_supports_arbitrary_binary_labels():
    lengths = [1, 2, 8, 9] * 2
    sequences = sequences_of_lengths(lengths)
    predictor = EventSequencePredictor(model_type="simple", task_type="classification")
    predictor.fit(sequences, np.array([2 if n < 5 else 7 for n in lengths]))

    assert predictor.classes_.tolist() == [2, 7]
    assert predictor.predict(sequences_of_lengths([1, 9])).tolist() == [2, 7]
    np.testing.assert_allclose(predictor.predict_proba(sequences_of_lengths([1, 9])), [[1.0, 0.0], [0.0, 1.0]])


def test_simple_classifier_supports_string_labels():
    lengths = [1, 2, 8, 9] * 2
    sequences = sequences_of_lengths(lengths)
    predictor = EventSequencePredictor(model_type="simple", task_type="classification")
    predictor.fit(sequences, np.array(["short" if n < 5 else "long" for n in lengths]))

    assert predictor.predict(sequences_of_lengths([1, 9])).tolist() == ["short", "long"]


def test_simple_classifier_rejects_more_than_two_classes():
    predictor = EventSequencePredictor(model_type="simple", task_type="classification")

    with pytest.raises(ValueError, match="binary classification"):
        predictor.fit(sequences_of_lengths([1, 2, 3]), np.array([0, 1, 2]))


def test_simple_regressor_predicts_the_mean_outcome_per_length():
    lengths = [1, 1, 3, 3, 3]
    predictor = EventSequencePredictor(model_type="simple", task_type="regression")
    predictor.fit(sequences_of_lengths(lengths), np.array([10.0, 20.0, 30.0, 33.0, 36.0]))

    predictions = predictor.predict(sequences_of_lengths([1, 3]))

    np.testing.assert_allclose(predictions, [15.0, 33.0])


def test_simple_model_uses_the_nearest_observed_length_for_unseen_lengths():
    lengths = [1, 1, 3, 3, 3, 8]
    predictor = EventSequencePredictor(model_type="simple", task_type="regression")
    predictor.fit(sequences_of_lengths(lengths), np.array([10.0, 20.0, 30.0, 33.0, 36.0, 100.0]))

    # 2 is as close to 1 as to 3 (the shorter wins); 0 -> 1; 5 -> 3; 6 -> 8; 40 -> 8.
    predictions = predictor.predict(sequences_of_lengths([2, 0, 5, 6, 40]))

    np.testing.assert_allclose(predictions, [15.0, 15.0, 33.0, 100.0, 100.0])


def test_simple_classifier_unseen_lengths_do_not_default_to_class_zero():
    sequences, outcomes = length_task()
    predictor = EventSequencePredictor(model_type="simple", task_type="classification").fit(sequences, outcomes)

    assert predictor.predict(sequences_of_lengths([20])).tolist() == [1]  # longer than any training sequence
    assert predictor.predict(sequences_of_lengths([0])).tolist() == [0]  # an empty sequence looks like the shortest
    assert predictor.predict(sequences_of_lengths([5])).tolist() == [0]  # equally far from 4 and 6: the shorter wins


def test_simple_classification_round_trip_keeps_the_length_model(tmp_path: Path):
    sequences, outcomes = length_task()
    original = EventSequencePredictor(model_type="simple", task_type="classification").fit(sequences, outcomes)
    path = tmp_path / "simple_classifier.json"
    original.save_model(path)

    loaded = EventSequencePredictor.load_model(path)

    assert loaded.sequence_length_model == original.sequence_length_model
    assert list(loaded.sequence_length_model) == sorted(set(len(seq) for seq in sequences))  # integer lengths again
    np.testing.assert_array_equal(loaded.predict(sequences), original.predict(sequences))
    np.testing.assert_array_equal(loaded.predict(sequences), outcomes)  # and they are not a constant
    np.testing.assert_allclose(loaded.predict_proba(sequences), original.predict_proba(sequences))


def test_simple_regression_round_trip_keeps_the_length_model(tmp_path: Path):
    lengths = [1, 1, 2, 4, 4]
    sequences = sequences_of_lengths(lengths)
    original = EventSequencePredictor(model_type="simple", task_type="regression")
    original.fit(sequences, np.array([1.5, 2.5, 7.0, 11.0, 13.0]))
    path = tmp_path / "simple_regressor.json"
    original.save_model(path)

    loaded = EventSequencePredictor.load_model(path)

    np.testing.assert_allclose(loaded.predict(sequences), original.predict(sequences))
    np.testing.assert_allclose(loaded.predict(sequences_of_lengths([2, 3])), [7.0, 7.0])
    assert json.loads(path.read_text())["sequence_length_model"] == {"1": 2.0, "2": 7.0, "4": 12.0}


def test_loading_a_simple_model_file_without_a_length_table_fails_loudly(tmp_path: Path):
    """Files written before the table was saved (regression models) cannot predict; do not pretend."""
    path = tmp_path / "legacy_regression.json"
    path.write_text(
        json.dumps(
            {
                "model_type": "simple",
                "task_type": "regression",
                "embedding_dim": 100,
                "random_seed": 42,
                "is_fitted": True,
                "event_embeddings": {},
                "vocab": {"event_0": 0},
                "event_vocab": {"event_0": 0},
                "reverse_vocab": {},
            }
        )
    )

    with pytest.raises(ValueError, match="sequence-length table"):
        EventSequencePredictor.load_model(path)


def test_loading_a_legacy_simple_classification_file_recovers_the_table_from_the_classifier(tmp_path: Path):
    path = tmp_path / "legacy_classification.json"
    path.write_text(
        json.dumps(
            {
                "model_type": "simple",
                "task_type": "classification",
                "embedding_dim": 100,
                "random_seed": 42,
                "is_fitted": True,
                "event_embeddings": {},
                "vocab": {"event_0": 0},
                "event_vocab": {"event_0": 0},
                "reverse_vocab": {},
                "classes_": [0, 1],
                "classifier": {"type": "simple", "model": {"1": 0.0, "6": 1.0}},
            }
        )
    )

    loaded = EventSequencePredictor.load_model(path)

    assert loaded.sequence_length_model == {1: 0.0, 6: 1.0}
    assert loaded.predict(sequences_of_lengths([1, 6])).tolist() == [0, 1]


def test_embedding_model_classification_is_unchanged():
    """The logistic embedding model still thresholds a sigmoid of its logit."""
    sequences = [["a", "b"], ["a", "c"], ["b", "c"], ["a", "b", "c"]] * 3
    outcomes = np.array([0, 1, 1, 0] * 3)
    predictor = EventSequencePredictor(
        model_type="embedding", task_type="classification", embedding_dim=8, random_state=3
    )
    predictor.fit(sequences, outcomes)

    probabilities = predictor.predict_proba(sequences)
    predictions = predictor.predict(sequences)

    np.testing.assert_allclose(probabilities.sum(axis=1), 1.0)
    np.testing.assert_array_equal(predictions, np.where(probabilities[:, 1] >= 0.5, 1, 0))


# --- LSTM ------------------------------------------------------------------------------------------


def test_lstm_is_not_an_outcome_model_and_the_error_says_where_to_look():
    sequences, outcomes = length_task()
    predictor = EventSequencePredictor(model_type="lstm", task_type="classification")

    with pytest.raises(ValueError, match="LSTMSequenceModel"):
        predictor.fit(sequences, outcomes)
    with pytest.raises(ValueError, match="'embedding' and 'simple'"):
        predictor.fit(sequences, outcomes)
    assert not predictor.is_fitted


def test_documentation_does_not_promise_an_lstm_model_type():
    init_doc = EventSequencePredictor.__init__.__doc__ or ""
    class_doc = EventSequencePredictor.__doc__ or ""

    assert "'lstm'" not in init_doc
    assert "'embedding' or 'simple'" in init_doc
    assert "LSTM models" not in class_doc  # the old "embedding-based, simple statistical, and LSTM models"
    assert "sequence_models" in class_doc


# --- input validation ----------------------------------------------------------------------------------


@pytest.mark.parametrize("model_type", ["simple", "embedding"])
def test_fit_rejects_sequences_and_outcomes_of_different_lengths(model_type):
    predictor = EventSequencePredictor(model_type=model_type, task_type="regression", embedding_dim=4)

    with pytest.raises(ValueError, match="same length"):
        predictor.fit(sequences_of_lengths([1, 2, 3]), np.array([1.0, 2.0]))


@pytest.mark.parametrize("model_type", ["simple", "embedding"])
def test_fit_rejects_empty_training_data(model_type):
    predictor = EventSequencePredictor(model_type=model_type, task_type="regression", embedding_dim=4)

    with pytest.raises(ValueError, match="cannot be empty"):
        predictor.fit([], np.array([]))


def test_predicting_with_an_unfitted_simple_model_is_an_error():
    with pytest.raises(ValueError, match="must be fitted"):
        EventSequencePredictor(model_type="simple").predict([["a"]])
