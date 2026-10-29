"""Round-4 regression tests for ``metainformant.life_events``.

- ``models.predictor``: ``EventSequencePredictor`` reseeded numpy's *global* generator in its constructor and drew
  its embeddings from it, in the iteration order of a ``set`` of strings. Two predictors built before fitting
  therefore overwrote each other's seed, and the same seed gave different embeddings in different interpreter
  processes (hash randomization). ``load_model`` also ran ``pickle.loads`` on any non-JSON ``classifier`` /
  ``regressor`` string in a model file, i.e. executed arbitrary code from the file.
- ``core.utils``: three synthetic-data generators seeded and consumed the interpreter's global ``random`` module,
  changing the caller's random state (their seeded output is unchanged by the fix).
- ``analysis.interpretability.attention_weights`` returned the string ``"attention_layer_present"`` per encoder
  layer as "attention weights" of a transformer, flagged supported.
- ``survival.recurrent_events`` crashed for a subject with no recorded times; ``visualization.statistical``
  labelled any projection as the method that was asked for, even when PCA was drawn.
"""

from __future__ import annotations

import json
import os
import pickle
import random
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, List

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest

import metainformant
from metainformant.life_events.analysis.interpretability import attention_weights
from metainformant.life_events.core.events import Event, EventSequence
from metainformant.life_events.core.utils import (
    generate_cohort_sequences,
    generate_event_chain,
    generate_synthetic_life_events,
)
from metainformant.life_events.models.predictor import EventSequencePredictor
from metainformant.life_events.survival.time_to_event import recurrent_events

matplotlib.use("Agg")

SOURCE_ROOT = Path(metainformant.__file__).resolve().parents[1]  # the tree under test: the child must import it too


def _sequences() -> List[EventSequence]:
    streams = [
        ["degree", "job", "marriage", "child"],
        ["degree", "job", "illness"],
        ["job", "illness", "retirement"],
        ["marriage", "child", "job", "retirement"],
        ["degree", "marriage", "child"],
        ["job", "retirement"],
    ]
    return [
        EventSequence(f"p{i}", [Event(kind, datetime(2000 + j, 1, 1), "life") for j, kind in enumerate(stream)])
        for i, stream in enumerate(streams)
    ]


OUTCOMES = [1, 0, 0, 1, 1, 0]


class TestPredictorRandomness:
    def test_construction_and_fitting_leave_the_global_numpy_state_alone(self) -> None:
        np.random.seed(123)
        before = np.random.get_state()[1].copy()

        predictor = EventSequencePredictor(model_type="embedding", embedding_dim=8, random_seed=7)
        predictor.fit(_sequences(), OUTCOMES)

        np.testing.assert_array_equal(np.random.get_state()[1], before)

    def test_embeddings_depend_on_the_seed_only(self) -> None:
        def embeddings(seed: int) -> dict:
            predictor = EventSequencePredictor(model_type="embedding", embedding_dim=8, random_seed=seed)
            predictor.fit(_sequences(), OUTCOMES)
            return {event: vector.tolist() for event, vector in predictor.embeddings.items()}

        alone = embeddings(7)
        # building (and fitting) predictors with other seeds in between must not change what seed 7 gives
        first = EventSequencePredictor(model_type="embedding", embedding_dim=8, random_seed=7)
        EventSequencePredictor(model_type="embedding", embedding_dim=8, random_seed=8).fit(_sequences(), OUTCOMES)
        first.fit(_sequences(), OUTCOMES)

        assert {event: vector.tolist() for event, vector in first.embeddings.items()} == alone
        assert embeddings(8) != alone

    def test_vocabulary_indices_are_in_sorted_event_order(self) -> None:
        predictor = EventSequencePredictor(model_type="embedding", embedding_dim=4, random_seed=1)
        predictor.fit(_sequences(), OUTCOMES)
        assert list(predictor.event_vocab) == sorted(predictor.event_vocab)
        assert list(predictor.event_vocab.values()) == list(range(len(predictor.event_vocab)))

    def test_embeddings_do_not_depend_on_the_python_hash_seed(self) -> None:
        script = (
            "import json\n"
            "from datetime import datetime\n"
            "from metainformant.life_events.core.events import Event, EventSequence\n"
            "from metainformant.life_events.models.predictor import EventSequencePredictor\n"
            "streams = [['degree','job','marriage','child'],['degree','job','illness'],['job','illness','retirement'],"
            "['marriage','child','job','retirement'],['degree','marriage','child'],['job','retirement']]\n"
            "seqs = [EventSequence(f'p{i}', [Event(k, datetime(2000 + j, 1, 1), 'life') for j, k in enumerate(s)])"
            " for i, s in enumerate(streams)]\n"
            "model = EventSequencePredictor(model_type='embedding', embedding_dim=6, random_seed=3)\n"
            "model.fit(seqs, [1, 0, 0, 1, 1, 0])\n"
            "print(json.dumps({k: v.tolist() for k, v in sorted(model.embeddings.items())}))\n"
        )
        outputs = []
        for hash_seed in ("1", "2", "31337"):
            env = dict(os.environ, PYTHONHASHSEED=hash_seed, PYTHONPATH=str(SOURCE_ROOT))
            completed = subprocess.run(
                [sys.executable, "-c", script], capture_output=True, text=True, timeout=120, env=env, check=True
            )
            outputs.append(json.loads(completed.stdout))
        assert outputs[0] == outputs[1] == outputs[2]


class _WritesAFile:
    """Unpickling this object calls ``Path(marker).write_text('executed')``."""

    def __init__(self, marker: Path) -> None:
        self.marker = marker

    def __reduce__(self) -> Any:
        return (self.marker.write_text, ("executed",))


class TestModelFilesAreNotUnpickled:
    @staticmethod
    def _saved_model(tmp_path: Path) -> dict:
        predictor = EventSequencePredictor(model_type="embedding", embedding_dim=4, random_seed=1)
        predictor.fit(_sequences(), OUTCOMES)
        path = tmp_path / "model.json"
        predictor.save_model(path)
        return json.loads(path.read_text())

    @pytest.mark.parametrize("component", ["classifier", "regressor"])
    def test_a_pickled_payload_is_rejected_and_never_run(self, tmp_path: Path, component: str) -> None:
        marker = tmp_path / "marker.txt"
        data = self._saved_model(tmp_path)
        data[component] = pickle.dumps(_WritesAFile(marker)).hex()
        path = tmp_path / "evil.json"
        path.write_text(json.dumps(data))

        with pytest.raises(ValueError, match="pickled"):
            EventSequencePredictor.load_model(path)

        assert not marker.exists()

    def test_a_saved_model_still_round_trips(self, tmp_path: Path) -> None:
        predictor = EventSequencePredictor(model_type="embedding", embedding_dim=4, random_seed=1)
        predictor.fit(_sequences(), OUTCOMES)
        path = tmp_path / "model.json"
        predictor.save_model(path)

        loaded = EventSequencePredictor.load_model(path)

        np.testing.assert_allclose(loaded.predict(_sequences()), predictor.predict(_sequences()))


class TestGeneratorsLeaveTheGlobalRandomState:
    @staticmethod
    def _global_state_after(call: Any) -> bool:
        random.seed(99)
        before = random.getstate()
        call()
        return random.getstate() == before

    def test_cohort_sequences(self) -> None:
        assert self._global_state_after(lambda: generate_cohort_sequences(n_sequences=3, random_seed=5))

    def test_event_chain(self) -> None:
        assert self._global_state_after(lambda: generate_event_chain("started_school", n_events=6, random_seed=3))

    def test_synthetic_life_events(self) -> None:
        assert self._global_state_after(lambda: generate_synthetic_life_events(n_sequences=2, random_state=4))

    def test_seeded_output_is_unchanged_by_the_private_generator(self) -> None:
        cohort = generate_cohort_sequences(n_sequences=3, random_seed=5)
        assert [(e.event_type, e.domain, e.timestamp.isoformat()) for e in cohort[0].events][:3] == [
            ("diagnosis", "health", "1980-12-27T00:00:00"),
            ("moved_house", "housing", "1986-11-17T00:00:00"),
            ("started_school", "education", "1993-08-25T00:00:00"),
        ]
        chain = generate_event_chain("started_school", n_events=6, random_seed=3)
        assert [getattr(e, "event_type", e) for e in chain] == [
            "started_school",
            "graduated_high_school",
            "graduated_high_school",
            "started_college",
            "got_first_job",
            "retirement",
        ]
        synthetic, _ = generate_synthetic_life_events(n_sequences=2, random_state=4)
        assert [(e.event_type, e.domain) for e in synthetic[0].events][:3] == [
            ("started_school", "education"),
            ("started_college", "education"),
            ("graduated_high_school", "education"),
        ]

    def test_unseeded_calls_differ_and_seeded_calls_repeat(self) -> None:
        def types(seed: Any) -> list:
            sequences = generate_cohort_sequences(n_sequences=4, random_seed=seed)
            return [e.event_type for s in sequences for e in s.events]

        assert types(11) == types(11)
        assert types(11) != types(12)


class TestAttentionWeights:
    def test_encoder_layers_are_not_reported_as_attention_weights(self) -> None:
        class Layer:
            self_attn = object()

        class Encoder:
            layers = [Layer(), Layer()]

        class TransformerLike:
            transformer = True
            encoder = Encoder()

            def eval(self) -> None:  # pragma: no cover - only reached by the old placeholder branch
                pass

        result = attention_weights(TransformerLike(), [["a", "b"]])
        assert result["supported"] is False
        assert result["attention_weights"] == []
        assert "attention_layer_present" not in str(result)

    def test_a_model_that_provides_weights_is_used(self) -> None:
        class WithWeights:
            attention = True

            def get_attention_weights(self, sequences: List[List[str]]) -> List[List[float]]:
                return [[1.0 / len(seq)] * len(seq) for seq in sequences]

        result = attention_weights(WithWeights(), [["a", "b"], ["c"]])
        assert result["supported"] is True
        assert result["attention_weights"] == [[0.5, 0.5], [1.0]]

    def test_linear_models_report_their_coefficients(self) -> None:
        pytest.importorskip("sklearn")
        from sklearn.linear_model import LogisticRegression

        model = LogisticRegression().fit(np.array([[0.0, 1.0], [1.0, 0.0], [1.0, 1.0], [0.0, 0.0]]), [0, 1, 1, 0])
        result = attention_weights(model, [["a"]])
        assert result["supported"] is True and result["method"] == "coefficients"
        np.testing.assert_allclose(result["attention_weights"], model.coef_.tolist())


class TestRecurrentEventsEdgeCases:
    def test_a_subject_without_times_is_never_at_risk(self) -> None:
        result = recurrent_events(
            times=[[1.0, 2.0, 4.0], []], events=[[1, 1, 0], []], subject_ids=["with_events", "no_records"]
        )
        assert result["total_events"] == 2
        assert result["rate_per_subject"]["no_records"] == 0.0
        # one subject at risk at times 1 and 2, each with one event: the MCF rises by 1 at each
        assert result["mean_cumulative_function"] == [(1.0, 1.0), (2.0, 2.0)]


class TestEmbeddingPlotLabels:
    @staticmethod
    def _embeddings() -> dict:
        rng = np.random.default_rng(0)
        return {"embeddings": rng.normal(size=(40, 6)).tolist()}  # t-SNE needs more samples than its perplexity (30)

    def test_unknown_method_is_an_error_not_a_pca_labelled_with_that_name(self) -> None:
        from metainformant.life_events.visualization.statistical import plot_event_embeddings

        with pytest.raises(ValueError, match="Unknown reduction method"):
            plot_event_embeddings(self._embeddings(), method="mds")

    @pytest.mark.parametrize("method", ["pca", "tsne"])
    def test_axes_carry_the_name_of_the_projection(self, method: str) -> None:
        pytest.importorskip("sklearn")
        from metainformant.life_events.visualization.statistical import plot_event_embeddings

        figure = plot_event_embeddings(self._embeddings(), method=method)
        try:
            assert figure.axes[0].get_xlabel() == f"{method.upper()} Component 1"
        finally:
            plt.close(figure)
