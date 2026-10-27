"""Regression tests: ``shapley_values(background_sequences=...)`` is deprecated, not silently ignored.

Defect: ``background_sequences`` was accepted and never used ("Kept for interface compatibility"); a caller who
passed a background dataset expected the Shapley values to be computed against it and got values that ignored it,
with only a debug log line as a hint. A background has no canonical meaning for variable-length event sequences
(an absent event is removed from the sequence, it has no value to be replaced by), so the parameter cannot
meaningfully be implemented. It is deprecated: a non-empty value now emits a ``DeprecationWarning`` that says it is
ignored, and the documentation says so.

Oracle: for an additive model the Shapley value of an event is its own weight (computed by hand below).
"""

from __future__ import annotations

import warnings
from typing import Dict, Sequence

import numpy as np
import pytest

from metainformant.life_events.analysis.interpretability import shapley_values


class AdditivePredictor:
    """A deterministic additive model: the score of a sequence is the sum of per-event weights."""

    is_fitted = True

    def __init__(self, weights: Dict[str, float]):
        self.weights = dict(weights)

    def predict(self, sequences: Sequence[Sequence[str]]) -> np.ndarray:
        return np.array([sum(self.weights.get(token, 0.0) for token in seq) for seq in sequences])


WEIGHTS = {"a": 1.0, "b": 2.0}
BACKGROUND = [["b"], ["a", "a"], []]


def test_a_non_empty_background_emits_a_deprecation_warning_that_says_it_is_ignored() -> None:
    with pytest.warns(DeprecationWarning, match="ignores background_sequences") as caught:
        shapley_values(AdditivePredictor(WEIGHTS), ["a", "b"], BACKGROUND)

    message = str(caught[0].message)
    assert "deprecated" in message
    assert "removed" in message  # tells the caller what happens to absent events instead


def test_the_warning_is_attributed_to_the_caller() -> None:
    with pytest.warns(DeprecationWarning) as caught:
        shapley_values(AdditivePredictor(WEIGHTS), ["a", "b"], BACKGROUND)

    assert caught[0].filename == __file__


def test_the_keyword_form_warns_too() -> None:
    with pytest.warns(DeprecationWarning):
        shapley_values(AdditivePredictor(WEIGHTS), ["a", "b"], background_sequences=BACKGROUND)


def test_the_values_are_unchanged_by_the_background() -> None:
    """Deprecated means ignored: the result is the model's exact Shapley values, here the weights themselves."""
    with pytest.warns(DeprecationWarning):
        with_background = shapley_values(AdditivePredictor(WEIGHTS), ["a", "b"], BACKGROUND)

    assert with_background == pytest.approx({"a": 1.0, "b": 2.0})
    assert with_background == shapley_values(AdditivePredictor(WEIGHTS), ["a", "b"])


@pytest.mark.parametrize("background", [None, [], ()])
def test_no_background_means_no_warning(background: object) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # any warning at all would raise
        values = shapley_values(AdditivePredictor(WEIGHTS), ["a", "b"], background)  # type: ignore[arg-type]

    assert values == pytest.approx({"a": 1.0, "b": 2.0})


def test_the_default_call_does_not_warn() -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        shapley_values(AdditivePredictor(WEIGHTS), ["a", "b"])


def test_the_deprecation_does_not_depend_on_the_sequence_having_events() -> None:
    """The warning concerns the call, so it fires for an empty sequence (which returns {}) as well."""
    with pytest.warns(DeprecationWarning):
        assert shapley_values(AdditivePredictor(WEIGHTS), [], BACKGROUND) == {}


def test_the_docstring_marks_the_parameter_as_deprecated_and_ignored() -> None:
    doc = " ".join((shapley_values.__doc__ or "").split())

    assert "background_sequences: Deprecated and ignored." in doc
    assert "Kept for interface compatibility" not in doc
