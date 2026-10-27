"""Regression tests for real word2vec training in ``learn_event_embeddings``.

``learn_event_embeddings`` used to return seeded random vectors plus keyword offsets: it never read
``window_size``, ``epochs``, ``sg``/``method`` or the sequences' co-occurrence structure. These tests
pin the behaviour of the NumPy skip-gram / CBOW negative-sampling implementation: deterministic
output for a seed, every hyper-parameter changes the result, and events that co-occur end up closer
than events that never do.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List

import numpy as np
import pytest

from metainformant.life_events.core.events import Event, EventSequence
from metainformant.life_events.models.embeddings import (
    biological_embedding,
    domain_specific_embeddings,
    learn_event_embeddings,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
METHODS = ["skipgram", "cbow"]


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))


def same_embeddings(a: Dict[str, np.ndarray], b: Dict[str, np.ndarray]) -> bool:
    return a.keys() == b.keys() and all(np.array_equal(a[key], b[key]) for key in a)


# Events A and B always occur together; C and D form a separate group that never appears near them.
PAIR_CORPUS: List[List[str]] = [["A", "B"]] * 10 + [["C", "D"]] * 10

# Two disjoint six-event chains: A..F and P..U.
CHAIN_CORPUS: List[List[str]] = [list("ABCDEF")] * 10 + [list("PQRSTU")] * 10

# Events with very different frequencies (A 25 times, B 20, C 10, D 5), so the noise distribution matters.
UNEVEN_CORPUS: List[List[str]] = [["A", "B", "C"]] * 10 + [["A", "B"]] * 10 + [["A", "D"]] * 5


@pytest.mark.parametrize("method", METHODS)
def test_same_seed_gives_identical_embeddings(method):
    first = learn_event_embeddings(PAIR_CORPUS, embedding_dim=16, epochs=5, method=method, random_state=7)
    second = learn_event_embeddings(PAIR_CORPUS, embedding_dim=16, epochs=5, method=method, random_state=7)

    assert same_embeddings(first, second)


def test_seed_keyword_is_an_alias_for_random_state():
    via_random_state = learn_event_embeddings(PAIR_CORPUS, embedding_dim=16, epochs=3, random_state=11)
    via_seed = learn_event_embeddings(PAIR_CORPUS, embedding_dim=16, epochs=3, seed=11)

    assert same_embeddings(via_random_state, via_seed)


def test_different_seeds_give_different_embeddings():
    first = learn_event_embeddings(PAIR_CORPUS, embedding_dim=16, epochs=3, random_state=1)
    second = learn_event_embeddings(PAIR_CORPUS, embedding_dim=16, epochs=3, random_state=2)

    assert not same_embeddings(first, second)


@pytest.mark.parametrize("method", METHODS)
def test_epochs_change_the_embeddings(method):
    # The seed is shared, so only the amount of training can explain a difference.
    short = learn_event_embeddings(PAIR_CORPUS, embedding_dim=16, epochs=1, method=method, random_state=3)
    long = learn_event_embeddings(PAIR_CORPUS, embedding_dim=16, epochs=8, method=method, random_state=3)

    assert not same_embeddings(short, long)


@pytest.mark.parametrize("method", METHODS)
def test_window_size_changes_the_embeddings(method):
    narrow = learn_event_embeddings(
        CHAIN_CORPUS, embedding_dim=16, window_size=1, epochs=3, method=method, random_state=3
    )
    wide = learn_event_embeddings(
        CHAIN_CORPUS, embedding_dim=16, window_size=4, epochs=3, method=method, random_state=3
    )

    assert not same_embeddings(narrow, wide)


def test_window_size_limits_which_events_become_similar():
    """With window 1 the ends of a six-event chain never meet; with window 5 they co-occur."""
    options = dict(embedding_dim=16, epochs=15, learning_rate=0.1, method="skipgram", random_state=0)

    narrow = learn_event_embeddings(CHAIN_CORPUS, window_size=1, **options)
    wide = learn_event_embeddings(CHAIN_CORPUS, window_size=5, **options)

    # Neighbours are much closer than the chain ends when only adjacent events are contexts ...
    assert cosine(narrow["A"], narrow["B"]) - cosine(narrow["A"], narrow["F"]) > 0.5
    # ... and the ends of the chain move together once the window spans the whole sequence.
    assert cosine(wide["A"], wide["F"]) - cosine(narrow["A"], narrow["F"]) > 0.5


@pytest.mark.parametrize("method", METHODS)
def test_cooccurring_events_are_closer_than_events_that_never_meet(method):
    embeddings = learn_event_embeddings(
        PAIR_CORPUS, embedding_dim=16, window_size=2, epochs=10, method=method, random_state=0
    )

    together = cosine(embeddings["A"], embeddings["B"])
    apart = cosine(embeddings["A"], embeddings["C"])

    assert together > apart + 1.0
    assert together > 0.5 > apart


def test_cooccurrence_structure_survives_a_shared_context():
    """A and B always co-occur and share neighbours X and Y; C only ever appears with other events."""
    corpus = [["X", "A", "B", "Y"]] * 12 + [["P", "C", "Q", "D"]] * 12

    embeddings = learn_event_embeddings(
        corpus, embedding_dim=16, window_size=2, epochs=20, method="skipgram", random_state=0
    )

    assert cosine(embeddings["A"], embeddings["B"]) > cosine(embeddings["A"], embeddings["C"]) + 0.5


def test_sg_flag_selects_the_training_algorithm():
    options = dict(embedding_dim=16, epochs=4, random_state=5)
    skipgram = learn_event_embeddings(PAIR_CORPUS, method="skipgram", **options)
    cbow = learn_event_embeddings(PAIR_CORPUS, method="cbow", **options)

    assert not same_embeddings(skipgram, cbow)
    assert same_embeddings(learn_event_embeddings(PAIR_CORPUS, sg=1, **options), skipgram)
    assert same_embeddings(learn_event_embeddings(PAIR_CORPUS, sg=0, **options), cbow)
    # An explicit method wins over sg.
    assert same_embeddings(learn_event_embeddings(PAIR_CORPUS, sg=1, method="cbow", **options), cbow)


def test_learning_rate_is_honoured():
    slow = learn_event_embeddings(PAIR_CORPUS, embedding_dim=16, epochs=3, random_state=2, learning_rate=0.0005)
    fast = learn_event_embeddings(PAIR_CORPUS, embedding_dim=16, epochs=3, random_state=2, learning_rate=0.1)
    aliased = learn_event_embeddings(PAIR_CORPUS, embedding_dim=16, epochs=3, random_state=2, alpha=0.1)

    assert not same_embeddings(slow, fast)
    assert same_embeddings(fast, aliased)


def test_min_count_restricts_the_vocabulary():
    sequences = [["common", "also_common", "rare"], ["common", "also_common"], ["common"]]

    everything = learn_event_embeddings(sequences, embedding_dim=8, epochs=2, random_state=0)
    frequent = learn_event_embeddings(sequences, embedding_dim=8, epochs=2, min_count=2, random_state=0)

    assert set(everything) == {"common", "also_common", "rare"}
    assert set(frequent) == {"common", "also_common"}


def test_min_count_removes_rare_events_from_the_contexts():
    # "R" occurs once. Dropped by min_count, it no longer separates A and B, so training sees (A, B)
    # as neighbours exactly as in a corpus that never contained R.
    with_rare = [["A", "R", "B"]] + [["A", "B"]] * 9 + [["C", "D"]] * 10
    without_rare = [["A", "B"]] * 10 + [["C", "D"]] * 10
    options = dict(embedding_dim=8, window_size=1, epochs=4, random_state=0)

    filtered = learn_event_embeddings(with_rare, min_count=2, **options)
    reference = learn_event_embeddings(without_rare, min_count=2, **options)
    unfiltered = learn_event_embeddings(with_rare, min_count=1, **options)

    assert set(filtered) == set(reference) == {"A", "B", "C", "D"}
    assert same_embeddings(filtered, reference)
    # Kept (min_count=1), the rare event separates A from B in the window and so changes what A learns.
    assert set(unfiltered) == {"A", "B", "C", "D", "R"}
    assert not np.array_equal(unfiltered["A"], reference["A"])


def test_min_count_that_removes_every_event_is_an_error():
    with pytest.raises(ValueError, match="No events found"):
        learn_event_embeddings([["a", "b"], ["c"]], embedding_dim=8, min_count=5)


@pytest.mark.parametrize("dim", [1, 7, 33])
def test_embedding_dim_is_honoured_and_vectors_are_finite(dim):
    embeddings = learn_event_embeddings(PAIR_CORPUS, embedding_dim=dim, epochs=3, random_state=0)

    assert set(embeddings) == {"A", "B", "C", "D"}
    for vector in embeddings.values():
        assert vector.shape == (dim,)
        assert np.all(np.isfinite(vector))


def test_event_sequence_objects_are_tokenised_with_their_domain():
    sequences = [
        EventSequence(
            f"p{i}",
            [
                Event("degree", datetime(2010, 1, 1), "education"),
                Event("job_change", datetime(2012, 1, 1), "occupation"),
            ],
        )
        for i in range(4)
    ]

    embeddings = learn_event_embeddings(sequences, embedding_dim=8, epochs=2, random_state=0)

    assert set(embeddings) == {"education:degree", "occupation:job_change"}


def test_sequences_without_context_still_get_vectors():
    """Single-event sequences have nothing to learn from but must not break the call."""
    embeddings = learn_event_embeddings([["a"], ["b"], ["a"]], embedding_dim=8, epochs=2, random_state=0)

    assert set(embeddings) == {"a", "b"}
    assert all(np.all(np.isfinite(vector)) for vector in embeddings.values())


@pytest.mark.parametrize(
    "options, message",
    [
        ({"epochs": 0}, "epochs must be positive"),
        ({"negative": 0}, "negative must be at least 1"),
        ({"learning_rate": 0.0}, "learning_rate must be positive"),
        ({"combine": "mean"}, "combine must be one of"),
    ],
)
def test_invalid_training_options_are_rejected(options, message):
    with pytest.raises(ValueError, match=message):
        learn_event_embeddings(PAIR_CORPUS, embedding_dim=8, **options)


@pytest.mark.parametrize(
    "method, low, high",
    [
        ("skipgram", {"negative": 1}, {"negative": 8}),
        ("cbow", {"negative": 1}, {"negative": 8}),
        ("skipgram", {"ns_exponent": 0.0}, {"ns_exponent": 1.0}),
        ("cbow", {"cbow_mean": False}, {"cbow_mean": True}),
        ("skipgram", {"min_learning_rate": 0.0}, {"min_learning_rate": 0.02}),
    ],
)
def test_every_tuning_option_changes_the_result(method, low, high):
    options = dict(embedding_dim=16, window_size=2, epochs=3, method=method, random_state=6)

    assert not same_embeddings(
        learn_event_embeddings(UNEVEN_CORPUS, **options, **low),
        learn_event_embeddings(UNEVEN_CORPUS, **options, **high),
    )


@pytest.mark.parametrize("cbow_mean", [True, False])
def test_cbow_learns_with_averaged_and_summed_contexts(cbow_mean):
    embeddings = learn_event_embeddings(
        PAIR_CORPUS,
        embedding_dim=16,
        window_size=2,
        epochs=15,
        method="cbow",
        cbow_mean=cbow_mean,
        learning_rate=0.1,
        random_state=0,
    )

    assert all(np.all(np.isfinite(vector)) for vector in embeddings.values())
    assert cosine(embeddings["A"], embeddings["B"]) > cosine(embeddings["A"], embeddings["C"]) + 0.5


def test_a_string_min_count_from_a_config_is_accepted():
    embeddings = learn_event_embeddings([["a", "b"], ["a"]], embedding_dim=4, epochs=1, min_count="2", random_state=0)

    assert set(embeddings) == {"a"}


def test_input_vectors_alone_are_available_with_combine_input():
    summed = learn_event_embeddings(PAIR_CORPUS, embedding_dim=16, epochs=3, random_state=4)
    input_only = learn_event_embeddings(PAIR_CORPUS, embedding_dim=16, epochs=3, random_state=4, combine="input")

    assert not same_embeddings(summed, input_only)
    assert set(summed) == set(input_only)


@pytest.mark.parametrize("skipgram", [True, False])
def test_training_reduces_the_negative_sampling_loss(skipgram):
    from metainformant.life_events.models.embeddings import _train_word2vec

    vocabulary = sorted({token for sentence in CHAIN_CORPUS for token in sentence})
    index = {token: i for i, token in enumerate(vocabulary)}
    sentences = [np.array([index[token] for token in sentence]) for sentence in CHAIN_CORPUS]
    counts = np.array([sum(sentence.count(token) for sentence in CHAIN_CORPUS) for token in vocabulary], dtype=float)

    _, output_vectors, losses = _train_word2vec(
        sentences,
        counts,
        embedding_dim=16,
        window_size=2,
        skipgram=skipgram,
        epochs=12,
        negative=5,
        learning_rate=0.1,
        min_learning_rate=1e-4,
        ns_exponent=0.75,
        cbow_mean=True,
        rng=np.random.default_rng(0),
    )

    assert len(losses) == 12
    assert losses[-1] < 0.6 * losses[0]
    # The output vectors start at zero in word2vec, so they only move if gradients really flow.
    assert np.linalg.norm(output_vectors) > 0.1


def _naive_word2vec(
    sentences,
    counts,
    *,
    embedding_dim,
    window_size,
    skipgram,
    epochs,
    negative,
    learning_rate,
    min_learning_rate,
    ns_exponent,
    cbow_mean,
    rng,
):
    """Straightforward loop implementation of the same training procedure, for cross-checking.

    It replays the random draws of ``_train_word2vec`` in the same order (initial vectors, sequence
    order, window sizes, noise events) and applies the negative-sampling gradients of one position at
    a time with plain Python loops, so it shares none of the vectorised indexing of the real code.
    """

    def sigmoid(x: float) -> float:
        return 1.0 / (1.0 + np.exp(-max(-35.0, min(35.0, x))))

    vocab_size = len(counts)
    w_in = (rng.random((vocab_size, embedding_dim)) - 0.5) / embedding_dim
    w_out = np.zeros((vocab_size, embedding_dim))
    noise = np.asarray(counts, dtype=float) ** ns_exponent
    noise_cdf = np.cumsum(noise / noise.sum())
    total_steps = max(1, epochs * sum(len(sentence) for sentence in sentences))
    step = 0
    losses = []
    for _ in range(epochs):
        loss_sum, loss_terms = 0.0, 0
        for sentence_index in rng.permutation(len(sentences)):
            sentence = [int(event) for event in sentences[int(sentence_index)]]
            length = len(sentence)
            if length < 2:
                step += length
                continue
            reach = rng.integers(1, window_size + 1, size=length)
            max_context = min(2 * window_size, length - 1)
            draws = rng.random((length, max_context, negative))
            noise_ids = np.minimum(np.searchsorted(noise_cdf, draws, side="right"), vocab_size - 1)
            for position in range(length):
                alpha = max(min_learning_rate, learning_rate - (learning_rate - min_learning_rate) * step / total_steps)
                step += 1
                half = int(reach[position])
                context = sentence[max(0, position - half) : position] + sentence[position + 1 : position + 1 + half]
                centre = sentence[position]
                input_gradient = np.zeros(embedding_dim)
                output_updates = []
                if skipgram:
                    hidden = w_in[centre].copy()
                    for k, target in enumerate(context):
                        examples = [(target, 1.0)] + [(int(n), 0.0) for n in noise_ids[position, k] if int(n) != target]
                        for event, label in examples:
                            score = sigmoid(float(w_out[event] @ hidden))
                            input_gradient += (label - score) * alpha * w_out[event]
                            output_updates.append((event, (label - score) * alpha * hidden))
                            loss_sum -= np.log(max(score if label else 1.0 - score, 1e-12))
                    loss_terms += len(context)
                else:
                    stacked = np.array([w_in[event] for event in context])
                    hidden = stacked.mean(axis=0) if cbow_mean else stacked.sum(axis=0)
                    examples = [(centre, 1.0)] + [(int(n), 0.0) for n in noise_ids[position, 0] if int(n) != centre]
                    for event, label in examples:
                        score = sigmoid(float(w_out[event] @ hidden))
                        input_gradient += (label - score) * alpha * w_out[event]
                        output_updates.append((event, (label - score) * alpha * hidden))
                        loss_sum -= np.log(max(score if label else 1.0 - score, 1e-12))
                    loss_terms += 1
                for event, update in output_updates:
                    w_out[event] += update
                if skipgram:
                    w_in[centre] += input_gradient
                else:
                    for event in context:
                        w_in[event] += input_gradient if cbow_mean else input_gradient / len(context)
        losses.append(loss_sum / loss_terms if loss_terms else 0.0)
    return w_in, w_out, losses


@pytest.mark.parametrize(
    "skipgram, cbow_mean, window_size",
    [(True, True, 1), (True, True, 3), (False, True, 1), (False, True, 3), (False, False, 2)],
)
def test_training_matches_a_naive_reference_implementation(skipgram, cbow_mean, window_size):
    """The vectorised updates must equal plain-loop negative-sampling updates to rounding error."""
    from metainformant.life_events.models.embeddings import _train_word2vec

    corpus = [list("ABCDE"), list("BCDEA"), list("EABCD"), list("AABBC"), list("CCDAB"), ["A"]]
    vocabulary = sorted({token for sentence in corpus for token in sentence})
    index = {token: i for i, token in enumerate(vocabulary)}
    sentences = [np.array([index[token] for token in sentence]) for sentence in corpus]
    counts = np.array([sum(sentence.count(token) for sentence in corpus) for token in vocabulary], dtype=float)
    options = dict(
        embedding_dim=6,
        window_size=window_size,
        skipgram=skipgram,
        epochs=3,
        negative=3,
        learning_rate=0.05,
        min_learning_rate=1e-4,
        ns_exponent=0.75,
        cbow_mean=cbow_mean,
    )

    input_vectors, output_vectors, losses = _train_word2vec(sentences, counts, rng=np.random.default_rng(4), **options)
    expected_in, expected_out, expected_losses = _naive_word2vec(
        sentences, counts, rng=np.random.default_rng(4), **options
    )

    np.testing.assert_allclose(input_vectors, expected_in, rtol=0, atol=1e-12)
    np.testing.assert_allclose(output_vectors, expected_out, rtol=0, atol=1e-12)
    np.testing.assert_allclose(losses, expected_losses, rtol=0, atol=1e-12)


def test_verbose_reports_every_epoch(capsys):
    learn_event_embeddings(PAIR_CORPUS, embedding_dim=8, epochs=3, random_state=0, verbose=True)

    printed = capsys.readouterr().out
    assert "Building embeddings for 4 events" in printed
    assert printed.count("Epoch ") == 3


def test_domain_specific_embeddings_learn_from_events_that_share_a_domain():
    """Events of one domain that co-occur in a sequence must train each other's embeddings.

    Every event used to be fed to the learner as its own one-event sequence, so nothing was learned.
    """
    sequences = [["a1", "a2", "z"]] * 8 + [["a3", "a4", "z"]] * 8
    domains = [["d", "d", "other"]] * 16

    embeddings = domain_specific_embeddings(
        sequences, domains, embedding_dim=16, window_size=2, epochs=15, learning_rate=0.1, random_state=0
    )

    within = embeddings["d"]
    assert set(within) == {"a1", "a2", "a3", "a4"}
    assert cosine(within["a1"], within["a2"]) > cosine(within["a1"], within["a3"]) + 0.5


def test_event_type_embeddings_do_not_depend_on_the_python_hash_seed():
    """biological_embedding('event_type') used the per-process salted built-in hash()."""
    sequences = [
        EventSequence(
            "p",
            [Event("degree", datetime(2010, 1, 1), "education"), Event("diagnosis", datetime(2012, 1, 1), "health")],
        )
    ]
    in_process = biological_embedding(sequences, embedding_type="event_type", embedding_dim=32)["embeddings"]
    assert set(in_process) == {"degree", "diagnosis"}

    script = (
        "import json\n"
        "from datetime import datetime\n"
        "from metainformant.life_events.core.events import Event, EventSequence\n"
        "from metainformant.life_events.models.embeddings import biological_embedding\n"
        "seq = EventSequence('p', [Event('degree', datetime(2010, 1, 1), 'education'),"
        " Event('diagnosis', datetime(2012, 1, 1), 'health')])\n"
        "result = biological_embedding([seq], embedding_type='event_type', embedding_dim=32)['embeddings']\n"
        "print(json.dumps({k: v.tolist() for k, v in sorted(result.items())}))\n"
    )
    outputs = []
    for hash_seed in ("1", "2"):
        env = dict(os.environ, PYTHONHASHSEED=hash_seed, PYTHONPATH=str(REPO_ROOT / "src"))
        completed = subprocess.run(
            [sys.executable, "-c", script], capture_output=True, text=True, timeout=120, env=env, check=True
        )
        outputs.append(json.loads(completed.stdout))

    assert outputs[0] == outputs[1]
    assert outputs[0] == {key: vector.tolist() for key, vector in sorted(in_process.items())}


def test_analysis_module_embeddings_are_trained_by_the_same_word2vec():
    """analysis.interpretability.learn_event_embeddings was a second, simplified implementation."""
    from metainformant.life_events.analysis.interpretability import learn_event_embeddings as analysis_embeddings

    embeddings = analysis_embeddings(PAIR_CORPUS, embedding_dim=16, window_size=2, epochs=10, random_state=0)
    again = analysis_embeddings(PAIR_CORPUS, embedding_dim=16, window_size=2, epochs=10, random_state=0)

    assert set(embeddings) == {"A", "B", "C", "D"}
    assert same_embeddings(embeddings, again)
    assert cosine(embeddings["A"], embeddings["B"]) > cosine(embeddings["A"], embeddings["C"]) + 1.0
    assert same_embeddings(
        embeddings,
        learn_event_embeddings(
            PAIR_CORPUS, embedding_dim=16, window_size=2, epochs=10, min_count=1, method="skipgram", random_state=0
        ),
    )


def test_analysis_module_embeddings_return_nothing_when_no_event_reaches_min_count():
    from metainformant.life_events.analysis.interpretability import learn_event_embeddings as analysis_embeddings

    assert analysis_embeddings([["a", "b"], ["c"]], min_count=3) == {}
