"""Regression tests: GSEA permutations and the pathway-network layout take a seed and leave the global state alone.

* ``gsea`` estimated its permutation p-values with the *global* ``random`` module and had no way to be seeded, so the
  p-values, NES and FDR changed on every call and every call consumed the caller's global random stream.
* ``pathway_network_plot`` called ``random.seed(42)`` before laying out the nodes: the layout was fixed, but every plot
  silently reset the caller's global ``random`` state to a known point.

``gsea`` now takes ``random_state`` (an int or a ``random.Random``) and ``pathway_network_plot`` a ``layout_seed``; the
same seed reproduces the output and the global state is neither seeded nor consumed. No test doubles.
"""

from __future__ import annotations

import random

import matplotlib
import pytest

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

from metainformant.ontology.pathway_enrichment.enrichment import gsea  # noqa: E402
from metainformant.ontology.visualization.plots import pathway_network_plot  # noqa: E402

RANKED = [(f"G{i:03d}", 5.0 - 0.1 * i) for i in range(50)]  # strongest first
_SET_RNG = random.Random(0)
GENE_SETS = {f"S{j}": _SET_RNG.sample([gene for gene, _ in RANKED], 10) for j in range(5)}


def _run(**kwargs):
    return gsea(RANKED, GENE_SETS, n_permutations=60, min_size=5, **kwargs)


def test_gsea_is_reproducible_for_a_seed() -> None:
    assert _run(random_state=3) == _run(random_state=3)


def test_gsea_p_values_depend_on_the_seed() -> None:
    p_values = {
        tuple(r["p_value"] for r in sorted(_run(random_state=seed), key=lambda r: r["term_id"])) for seed in range(6)
    }

    assert len(p_values) > 1


def test_gsea_accepts_a_random_generator_instance() -> None:
    assert _run(random_state=random.Random(5)) == _run(random_state=5)


@pytest.mark.parametrize("random_state", [3, None])
def test_gsea_does_not_touch_the_global_random_state(random_state: int | None) -> None:
    random.seed(20261027)
    before = random.getstate()

    _run(random_state=random_state)

    assert random.getstate() == before


def test_gsea_observed_enrichment_scores_do_not_depend_on_the_seed() -> None:
    first = {r["term_id"]: r["es"] for r in _run(random_state=1)}
    second = {r["term_id"]: r["es"] for r in _run(random_state=2)}

    assert first == second
    assert set(first) == set(GENE_SETS)


NETWORK = {
    "nodes": [
        {"term_id": f"GO:{i:07d}", "p_value": 0.001 * (i + 1), "adjusted_p": 0.002 * (i + 1), "n_genes": 5 + i}
        for i in range(7)
    ],
    "edges": [
        {"source": "GO:0000000", "target": "GO:0000001", "jaccard": 0.6},
        {"source": "GO:0000001", "target": "GO:0000002", "jaccard": 0.4},
        {"source": "GO:0000003", "target": "GO:0000004", "jaccard": 0.5},
    ],
    "clusters": [],
}


def _node_positions(**kwargs) -> list[tuple[float, float]]:
    fig = pathway_network_plot(NETWORK, **kwargs)
    try:
        return [tuple(collection.get_offsets()[0]) for collection in fig.axes[0].collections]
    finally:
        plt.close(fig)


def test_layout_is_reproducible_for_a_seed_and_differs_between_seeds() -> None:
    first = _node_positions(layout_seed=7)

    assert len(first) == len(NETWORK["nodes"])
    assert first == _node_positions(layout_seed=7)
    assert first != _node_positions(layout_seed=8)


def test_the_default_layout_seed_keeps_the_layout_stable() -> None:
    assert _node_positions() == _node_positions(layout_seed=42)


def test_plotting_does_not_reset_the_global_random_state() -> None:
    random.seed(20261027)
    before = random.getstate()

    _node_positions()

    assert random.getstate() == before
