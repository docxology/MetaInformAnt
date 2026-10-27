"""Regression tests for ``biological_embedding('domain')``.

The domain embeddings added ``np.random.normal`` noise from NumPy's global, unseeded generator, so two
calls on the same sequences returned different vectors, nothing could be reproduced with a seed, and
the result changed from one process to the next. The variation is now drawn per domain from a seeded
generator (``random_state``, ``seed`` as an alias, a fixed default seed otherwise).
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

import metainformant
from metainformant.life_events.core.events import Event, EventSequence
from metainformant.life_events.models.embeddings import biological_embedding

# The source tree the imported package lives in, so a subprocess exercises the same code as this process.
SOURCE_ROOT = Path(metainformant.__file__).resolve().parents[1]

# The one-line program the subprocess tests run; it prints the embeddings as JSON.
PROGRAM = (
    "import json\n"
    "from datetime import datetime\n"
    "from metainformant.life_events.core.events import Event, EventSequence\n"
    "from metainformant.life_events.models.embeddings import biological_embedding\n"
    "seq = EventSequence('p', [Event('e1', datetime(2000, 1, 1), 'education'),"
    " Event('e2', datetime(2001, 1, 1), 'career'), Event('e3', datetime(2002, 1, 1), 'health'),"
    " Event('e4', datetime(2003, 1, 1), 'sports')])\n"
    "result = biological_embedding([seq], embedding_type='domain', embedding_dim=8)['embeddings']\n"
    "print(json.dumps({k: v.tolist() for k, v in sorted(result.items())}))\n"
)


def sequences_in(*domains: str) -> List[EventSequence]:
    events = [Event(f"event_{i}", datetime(2000 + i, 1, 1), domain) for i, domain in enumerate(domains)]
    return [EventSequence("person", events)]


def domain_embeddings(*domains: str, **options) -> Dict[str, np.ndarray]:
    result = biological_embedding(sequences_in(*domains), embedding_type="domain", embedding_dim=10, **options)
    return result["embeddings"]


def same(a: Dict[str, np.ndarray], b: Dict[str, np.ndarray]) -> bool:
    return a.keys() == b.keys() and all(np.array_equal(a[key], b[key]) for key in a)


def test_repeated_calls_give_identical_embeddings_by_default():
    first = domain_embeddings("education", "career", "health")
    second = domain_embeddings("education", "career", "health")

    assert same(first, second)


def test_embeddings_do_not_depend_on_numpys_global_random_state():
    np.random.seed(1)
    first = domain_embeddings("education", "career")
    np.random.seed(2)
    second = domain_embeddings("education", "career")

    assert same(first, second)


def test_the_result_does_not_consume_the_global_random_stream():
    np.random.seed(123)
    expected_next = np.random.random()
    np.random.seed(123)
    domain_embeddings("education", "career", "health")

    assert np.random.random() == expected_next


def test_random_state_selects_the_variation():
    default = domain_embeddings("education", "career")
    seeded = domain_embeddings("education", "career", random_state=5)
    again = domain_embeddings("education", "career", random_state=5)
    other = domain_embeddings("education", "career", random_state=6)

    assert same(seeded, again)
    assert not same(seeded, other)
    assert not same(seeded, default)


def test_seed_is_an_alias_for_random_state():
    assert same(domain_embeddings("career", random_state=9), domain_embeddings("career", seed=9))
    assert same(domain_embeddings("career", random_state=9, seed=1), domain_embeddings("career", random_state=9))


def test_the_default_is_a_fixed_seed_that_is_reported_in_the_result():
    sequences = sequences_in("education", "career")
    default = biological_embedding(sequences, embedding_type="domain", embedding_dim=10)
    explicit = biological_embedding(
        sequences, embedding_type="domain", embedding_dim=10, random_state=default["random_state"]
    )

    assert isinstance(default["random_state"], int)
    assert same(default["embeddings"], explicit["embeddings"])
    assert biological_embedding(sequences, embedding_type="domain", random_state=3)["random_state"] == 3


def test_a_domains_embedding_does_not_depend_on_which_other_domains_are_present():
    alone = domain_embeddings("education")
    together = domain_embeddings("career", "education", "health")
    reordered = domain_embeddings("health", "education")

    assert np.array_equal(alone["education"], together["education"])
    assert np.array_equal(alone["education"], reordered["education"])


def test_the_domain_list_is_sorted():
    result = biological_embedding(sequences_in("health", "education", "career"), embedding_type="domain")

    assert result["domains"] == ["career", "education", "health"]


def test_variation_is_small_around_the_domain_centre_and_separates_unknown_domains():
    embeddings = domain_embeddings("education", "alpha", "beta")

    # "education" has its own centre (1, 0, 0); the variation has a standard deviation of 0.1.
    assert np.abs(embeddings["education"][:3] - np.array([1.0, 0.0, 0.0])).max() < 0.5
    assert np.abs(embeddings["education"][3:]).max() < 0.5
    # "alpha" and "beta" share the generic centre, yet get different vectors.
    assert not np.array_equal(embeddings["alpha"], embeddings["beta"])
    assert np.abs(embeddings["alpha"][:3] - 0.5).max() < 0.5


def test_event_type_embeddings_ignore_random_state():
    sequences = sequences_in("education", "career")
    plain = biological_embedding(sequences, embedding_type="event_type", embedding_dim=12)["embeddings"]
    seeded = biological_embedding(sequences, embedding_type="event_type", embedding_dim=12, random_state=4)[
        "embeddings"
    ]

    assert same(plain, seeded)


def run_in_new_process(hash_seed: str) -> Dict[str, List[float]]:
    env = dict(os.environ, PYTHONHASHSEED=hash_seed, PYTHONPATH=str(SOURCE_ROOT))
    completed = subprocess.run(
        [sys.executable, "-c", PROGRAM], capture_output=True, text=True, timeout=120, env=env, check=True
    )
    return json.loads(completed.stdout.strip().splitlines()[-1])


@pytest.mark.parametrize("other_hash_seed", ["1", "2"])
def test_embeddings_are_the_same_in_another_process(other_hash_seed: str):
    """Iterating a set of strings (the domains) changes with PYTHONHASHSEED; the result must not."""
    sequence = EventSequence(
        "p",
        [
            Event("e1", datetime(2000, 1, 1), "education"),
            Event("e2", datetime(2001, 1, 1), "career"),
            Event("e3", datetime(2002, 1, 1), "health"),
            Event("e4", datetime(2003, 1, 1), "sports"),
        ],
    )
    in_process = biological_embedding([sequence], embedding_type="domain", embedding_dim=8)["embeddings"]

    other = run_in_new_process(other_hash_seed)

    assert other == {key: vector.tolist() for key, vector in sorted(in_process.items())}
