"""Regression tests for ``analyze_life_course``: NumPy outcomes and embedding settings from the config.

``if outcomes:`` raised ``ValueError: The truth value of an array ... is ambiguous`` for a NumPy
array, and the surrounding ``except Exception`` logged a warning and dropped the outcome plot without
failing the run. The embedding section of a workflow config (``method``, ``learning_rate``) was
never passed to the embedding learner.
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import List

import numpy as np
import pytest

from metainformant.core import io
from metainformant.life_events.core.events import Event, EventSequence
from metainformant.life_events.workflow.workflow import _generate_analysis_visualizations, analyze_life_course


def make_sequences() -> List[EventSequence]:
    plans = [
        ["degree", "job_change", "promotion"],
        ["degree", "marriage"],
        ["diagnosis", "recovered", "job_change"],
        ["job_change", "promotion", "retirement"],
        ["marriage", "had_child", "diagnosis"],
        ["degree", "job_change", "marriage", "had_child"],
    ]
    return [
        EventSequence(
            f"person_{i}",
            [Event(name, datetime(2000 + 3 * j, 1 + j, 1), "life") for j, name in enumerate(plan)],
        )
        for i, plan in enumerate(plans)
    ]


def plot_types(results) -> dict:
    return {entry["type"]: entry["file"] for entry in results["visualizations"]}


def test_outcome_distribution_is_plotted_for_numpy_outcomes(tmp_path: Path, caplog: pytest.LogCaptureFixture):
    outcomes = np.array([0, 1, 0, 1, 1, 0])

    with caplog.at_level(logging.WARNING):
        results = analyze_life_course(
            make_sequences(), outcomes=outcomes, output_dir=tmp_path / "out", embedding_dim=8, epochs=1
        )

    plots = plot_types(results)
    assert "outcome_distribution" in plots
    assert Path(plots["outcome_distribution"]).stat().st_size > 0
    assert "Error generating visualizations" not in caplog.text
    # The other plots are unaffected.
    assert {"sequence_length_distribution", "domain_distribution"} <= set(plots)


def test_outcome_distribution_is_plotted_for_list_outcomes(tmp_path: Path):
    results = analyze_life_course(
        make_sequences(), outcomes=[0, 1, 0, 1, 1, 0], output_dir=tmp_path / "out", embedding_dim=8, epochs=1
    )

    assert "outcome_distribution" in plot_types(results)


@pytest.mark.parametrize("outcomes", [np.array([0, 1, 1]), [0, 1, 1], np.array([0.5, 2.5, 1.5])])
def test_visualizations_accept_array_and_list_outcomes(tmp_path: Path, outcomes, caplog: pytest.LogCaptureFixture):
    with caplog.at_level(logging.WARNING):
        visualizations = _generate_analysis_visualizations(make_sequences(), outcomes, tmp_path, {})

    assert "outcome_distribution" in {entry["type"] for entry in visualizations}
    assert (tmp_path / "outcome_distribution.png").exists()
    assert "Error generating visualizations" not in caplog.text


@pytest.mark.parametrize("outcomes", [None, [], np.array([])])
def test_no_outcome_plot_without_outcomes(tmp_path: Path, outcomes, caplog: pytest.LogCaptureFixture):
    with caplog.at_level(logging.WARNING):
        visualizations = _generate_analysis_visualizations(make_sequences(), outcomes, tmp_path, {})

    assert "outcome_distribution" not in {entry["type"] for entry in visualizations}
    assert not (tmp_path / "outcome_distribution.png").exists()
    assert {"sequence_length_distribution", "domain_distribution"} <= {entry["type"] for entry in visualizations}
    assert "Error generating visualizations" not in caplog.text


def learned_embeddings(tmp_path: Path, name: str, **embedding_section) -> dict:
    section = {"embedding_dim": 8, "epochs": 3, **embedding_section}
    results = analyze_life_course(make_sequences(), output_dir=tmp_path / name, config_obj={"embedding": section})
    return io.load_json(results["embeddings"])


def test_embedding_method_in_the_config_selects_the_training_algorithm(tmp_path: Path):
    skipgram = learned_embeddings(tmp_path, "skipgram", method="skipgram")
    cbow = learned_embeddings(tmp_path, "cbow", method="cbow")
    default = learned_embeddings(tmp_path, "default")

    assert skipgram != cbow
    assert default == skipgram  # skip-gram is the default


def test_embedding_learning_rate_in_the_config_is_used(tmp_path: Path):
    slow = learned_embeddings(tmp_path, "slow", learning_rate=0.0005)
    fast = learned_embeddings(tmp_path, "fast", learning_rate=0.1)

    assert slow != fast


def test_an_unknown_embedding_method_in_the_config_is_an_error(tmp_path: Path):
    with pytest.raises(ValueError, match="method must be"):
        learned_embeddings(tmp_path, "bad", method="word2vec")
