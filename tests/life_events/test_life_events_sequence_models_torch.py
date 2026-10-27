"""Regression tests for the optional PyTorch dependency of the LSTM/GRU sequence models.

Without PyTorch, ``fit`` used to log a warning, set ``is_fitted = True`` and return, so the model
looked trained although it never saw the data, and ``predict`` returned a made-up uniform
distribution. The tests make ``import torch`` fail with ``sys.modules["torch"] = None`` (Python's
documented way to block an import), so they behave the same whether or not PyTorch is installed.
"""

from __future__ import annotations

import sys

import numpy as np
import pytest

from metainformant.life_events.models.sequence_models import GRUSequenceModel, LSTMSequenceModel

SEQUENCES = [
    ["health:diagnosis", "occupation:job_change", "income:raise"],
    ["education:degree", "occupation:job_change", "income:raise"],
]

MODELS = [LSTMSequenceModel, GRUSequenceModel]


@pytest.fixture
def torch_is_missing(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "torch", None)


@pytest.mark.parametrize("model_class", MODELS)
def test_fit_without_torch_raises_import_error_with_guidance(torch_is_missing, model_class):
    model = model_class(embedding_dim=8, hidden_dim=4, epochs=1, random_state=0)

    with pytest.raises(ImportError, match="PyTorch") as error:
        model.fit(SEQUENCES, np.array([0, 1]))

    assert model_class.__name__ in str(error.value)
    assert "uv pip install torch" in str(error.value)


@pytest.mark.parametrize("model_class", MODELS)
def test_a_failed_fit_does_not_leave_a_model_that_claims_to_be_fitted(torch_is_missing, model_class):
    model = model_class(embedding_dim=8, hidden_dim=4, epochs=1)

    with pytest.raises(ImportError):
        model.fit(SEQUENCES)

    assert model.is_fitted is False
    assert model.model is None
    assert model.event_vocab is None


@pytest.mark.parametrize("model_class", MODELS)
def test_predict_before_fit_is_an_error_not_a_made_up_distribution(model_class):
    model = model_class(embedding_dim=8, hidden_dim=4)

    with pytest.raises(ValueError, match="must be fitted"):
        model.predict(SEQUENCES)


@pytest.mark.parametrize("model_class", MODELS)
def test_predict_after_a_failed_fit_is_an_error(torch_is_missing, model_class):
    model = model_class(embedding_dim=8, hidden_dim=4, epochs=1)
    with pytest.raises(ImportError):
        model.fit(SEQUENCES)

    with pytest.raises(ValueError, match="must be fitted"):
        model.predict(SEQUENCES)


def test_documentation_states_that_targets_are_not_used():
    doc = LSTMSequenceModel.fit.__doc__ or ""

    assert "Not used" in doc
    assert "ImportError" in doc
