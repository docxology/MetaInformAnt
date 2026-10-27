"""The multi-omic pathway analyses return the same result in every process.

Randomness audit of ``metainformant.multiomics.pathways.enrichment``: the module has no unseeded ``np.random`` and no
global ``random`` (its two permutation tests use ``random.Random(42)``), but three analyses iterated over a ``set`` of
strings. String hashes are randomised per process (``PYTHONHASHSEED``), so set iteration order, and with it

* ``pathway_topology_analysis``: which node each shuffled score is paired with in the permutation null, so the
  permutation p-value itself (0.937 under one hash seed and 0.270 under another on the graph below), and the order of
  the float sums;
* ``active_module_detection``: which of several equally scored neighbours is added first, and so the last digits of
  ``module_score`` (the float sum is order dependent);
* ``cross_omic_pathway_concordance``: the order of tied pathways in ``heatmap_data`` and the last digits of
  ``concordance_score``,

changed from one run to the next for identical input.

Each test runs the same analyses in child processes under different ``PYTHONHASHSEED`` values (the only way to vary
string hashing, which is fixed at interpreter start) and requires identical output. Real computation in real
subprocesses; no test doubles.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

import metainformant
from metainformant.multiomics.pathways.enrichment import cross_omic_pathway_concordance

SEEDS = ("0", "2", "3")
# The directory holding the ``metainformant`` package this process imported. Children must import the very same
# package: the suite-wide fixture in conftest.py prepends the worktree's ``src`` to PYTHONPATH, which would
# otherwise decide which copy a child sees.
PACKAGE_ROOT = str(Path(metainformant.__file__).resolve().parents[1])

CHILD_PROGRAM = r"""
import json
import random

from metainformant.multiomics.pathways.enrichment import (
    active_module_detection,
    cross_omic_pathway_concordance,
    pathway_topology_analysis,
)


def tie_heavy_network(seed):
    # scores drawn from a few levels, so many nodes tie
    rng = random.Random(seed)
    genes = [f"G{i:02d}" for i in range(30)]
    network = {g: [] for g in genes}
    for _ in range(45):
        a, b = rng.sample(genes, 2)
        network[a].append(b)
        network[b].append(a)
    levels = [0.0005, 0.001, 0.002, 0.004, 0.01, 0.03, 0.5, 0.9]
    return network, {g: rng.choice(levels) for g in genes}


modules = []
for case in range(25):
    network, scores = tie_heavy_network(case)
    found = active_module_detection(network, scores, alpha=0.05, n_permutations=30)
    modules.append([(m["module_genes"], repr(m["module_score"]), m["p_value"]) for m in found])

genes = [f"G{i:02d}" for i in range(24)]
topologies = []
for step in (1, 5, 7, 11):
    graph = {g: [genes[(i + 1) % 24], genes[(i + step) % 24]] for i, g in enumerate(genes)}
    node_p = {g: 0.001 * (1 + (i * step) % 13) for i, g in enumerate(genes)}
    result = pathway_topology_analysis(graph, node_p)
    topologies.append(
        [
            repr(result["impact_factor"]),
            result["p_value"],
            result["perturbed_genes"],
            list(result["pathway_perturbation"]),
        ]
    )

tied = {"rna": {f"PW{i}": 0.5 for i in range(12)}, "prot": {f"PW{i}": 0.5 for i in range(12)}}
untied = {
    "rna": {f"PW{i}": 0.5 / (1 + i) for i in range(12)},
    "prot": {f"PW{i}": 0.2 / (1 + (i * 5) % 12) for i in range(12)},
}
concordance = []
for pathway_results in (tied, untied):
    result = cross_omic_pathway_concordance(pathway_results)
    concordance.append(
        [
            list(result["heatmap_data"]),
            repr(result["concordance_score"]),
            result["concordant_pathways"],
            result["discordant_pathways"],
        ]
    )

print(json.dumps({"modules": modules, "topologies": topologies, "concordance": concordance}))
"""


@pytest.fixture(scope="module")
def runs() -> dict[str, dict[str, Any]]:
    """The analyses' output in one child process per hash seed."""
    results: dict[str, dict[str, Any]] = {}
    for seed in SEEDS:
        python_path = os.pathsep.join(filter(None, [PACKAGE_ROOT, os.environ.get("PYTHONPATH")]))
        completed = subprocess.run(
            [sys.executable, "-c", CHILD_PROGRAM],
            env={**os.environ, "PYTHONHASHSEED": seed, "PYTHONPATH": python_path},
            capture_output=True,
            text=True,
            timeout=600,
        )
        assert completed.returncode == 0, f"PYTHONHASHSEED={seed}: {completed.stderr[-2000:]}"
        results[seed] = json.loads(completed.stdout)
    return results


def _assert_same_for_every_seed(runs: dict[str, dict[str, Any]], key: str) -> None:
    reference = runs[SEEDS[0]][key]
    for seed in SEEDS[1:]:
        assert runs[seed][key] == reference, f"{key} differs between PYTHONHASHSEED={SEEDS[0]} and ={seed}"


def test_pathway_topology_analysis_does_not_depend_on_string_hashing(runs: dict[str, dict[str, Any]]) -> None:
    _assert_same_for_every_seed(runs, "topologies")


def test_active_module_detection_does_not_depend_on_string_hashing(runs: dict[str, dict[str, Any]]) -> None:
    _assert_same_for_every_seed(runs, "modules")


def test_cross_omic_pathway_concordance_does_not_depend_on_string_hashing(runs: dict[str, dict[str, Any]]) -> None:
    _assert_same_for_every_seed(runs, "concordance")


def test_tied_pathways_are_listed_in_identifier_order() -> None:
    # every pathway has the same CV; a set would hand them back in arbitrary order, the contract is by identifier
    ids = [f"PW{i:02d}" for i in range(16)]
    pathway_results = {"rna": {pw: 0.5 for pw in ids}, "prot": {pw: 0.5 for pw in ids}}

    result = cross_omic_pathway_concordance(pathway_results)

    assert list(result["heatmap_data"]) == ids
    assert result["concordant_pathways"] == ids
