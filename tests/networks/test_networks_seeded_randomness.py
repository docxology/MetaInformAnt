"""Regression tests: the randomized network analyses take a seed and leave the global random state alone.

* ``compute_network_motifs`` builds its z-scores from 100 randomized networks but shuffled with the *global*
  ``random`` module, so the z-scores (and the ``significant_motifs`` list derived from them) changed from run to run and
  could not be pinned, and every call consumed the caller's global random stream.
* ``analyze_regulatory_dynamics`` added noise drawn from the global NumPy generator, with the same two problems.

Both now take ``random_state`` (an int or a generator). The same seed reproduces the output, different seeds give
different noise, and the global ``random`` / ``numpy.random`` state is neither seeded nor consumed. The dynamics oracle
is the update rule worked out by hand for the first step, using an independently created NumPy generator.
No test doubles.
"""

from __future__ import annotations

import random

import networkx as nx
import numpy as np
import pytest

from metainformant.networks.interaction.regulatory_core import analyze_regulatory_dynamics
from metainformant.networks.regulatory.grn_inference import compute_network_motifs


def _motif_edges() -> list[dict]:
    rng = random.Random(7)
    nodes = [f"G{i}" for i in range(12)]
    pairs = {(a, b) for a in nodes for b in nodes if a != b and rng.random() < 0.22}
    return [{"source": a, "target": b} for a, b in sorted(pairs)]


def test_motif_z_scores_are_reproducible_for_a_seed() -> None:
    edges = _motif_edges()

    first = compute_network_motifs(edges, motif_size=3, random_state=11)
    second = compute_network_motifs(edges, motif_size=3, random_state=11)

    assert first == second
    assert first["motif_counts"]["cascade"] > 0


def test_motif_z_scores_depend_on_the_seed() -> None:
    edges = _motif_edges()

    z_scores = {tuple(compute_network_motifs(edges, 3, random_state=seed)["z_scores"].values()) for seed in range(8)}

    assert len(z_scores) > 1


def test_a_random_generator_instance_gives_the_same_ensemble_as_its_seed() -> None:
    edges = _motif_edges()

    assert compute_network_motifs(edges, 3, random_state=random.Random(5)) == compute_network_motifs(
        edges, 3, random_state=5
    )


@pytest.mark.parametrize("random_state", [5, None])
def test_motif_analysis_does_not_touch_the_global_random_state(random_state: int | None) -> None:
    edges = _motif_edges()
    random.seed(20261027)
    before = random.getstate()

    compute_network_motifs(edges, 3, random_state=random_state)

    assert random.getstate() == before


def test_observed_motif_counts_do_not_depend_on_the_seed() -> None:
    edges = _motif_edges()

    counts = {tuple(compute_network_motifs(edges, 3, random_state=s)["motif_counts"].items()) for s in range(4)}

    assert len(counts) == 1


def _regulatory_graph() -> nx.DiGraph:
    graph = nx.DiGraph()
    graph.add_edges_from([("TF", "A"), ("A", "B"), ("TF", "B"), ("B", "C")])
    return graph


INITIAL = {"TF": 1.0, "A": 0.2, "B": -0.3, "C": 0.5}


def test_dynamics_are_reproducible_for_a_seed_and_differ_between_seeds() -> None:
    graph = _regulatory_graph()

    first = analyze_regulatory_dynamics(graph, INITIAL, time_steps=12, random_state=4)
    second = analyze_regulatory_dynamics(graph, INITIAL, time_steps=12, random_state=4)
    other = analyze_regulatory_dynamics(graph, INITIAL, time_steps=12, random_state=5)

    assert first == second
    assert first["dynamics"]["C"]["trajectory"] != other["dynamics"]["C"]["trajectory"]


def test_dynamics_first_step_follows_the_documented_update_rule_with_the_seeded_noise() -> None:
    graph = _regulatory_graph()
    nodes = list(graph.nodes())  # ["TF", "A", "B", "C"]
    noise = np.random.default_rng(9)

    result = analyze_regulatory_dynamics(graph, INITIAL, time_steps=1, random_state=9)

    trajectory = {node: result["dynamics"][node]["trajectory"] for node in nodes}
    # TF has no regulator: pure decay, no random draw
    assert trajectory["TF"][1] == pytest.approx(0.9 * INITIAL["TF"])
    # the remaining nodes draw their noise in node order: x1 = 0.9 * mean(regulators at t0) + 0.1 * x0 + noise
    for node in ("A", "B", "C"):
        regulators = list(graph.predecessors(node))
        expected = 0.9 * np.mean([INITIAL[r] for r in regulators]) + 0.1 * INITIAL[node] + noise.normal(0, 0.1)
        assert trajectory[node][1] == pytest.approx(expected)


def test_dynamics_accept_a_numpy_generator() -> None:
    graph = _regulatory_graph()

    from_int = analyze_regulatory_dynamics(graph, INITIAL, time_steps=6, random_state=3)
    from_generator = analyze_regulatory_dynamics(graph, INITIAL, time_steps=6, random_state=np.random.default_rng(3))

    assert from_int == from_generator


@pytest.mark.parametrize("random_state", [3, None])
def test_dynamics_do_not_touch_the_global_numpy_state(random_state: int | None) -> None:
    graph = _regulatory_graph()
    np.random.seed(20261027)
    before = np.random.get_state()

    analyze_regulatory_dynamics(graph, INITIAL, time_steps=6, random_state=random_state)

    after = np.random.get_state()
    assert before[0] == after[0] and before[2:] == after[2:]
    assert np.array_equal(before[1], after[1])
