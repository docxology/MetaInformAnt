"""Randomness audit of the spatial analysis paths.

Library code must neither read nor advance numpy's global random state, and every draw must be
controllable through a seed. Real data and real library calls only.
"""

from __future__ import annotations

from typing import Any, Callable

import numpy as np

from metainformant.spatial.communication.cell_communication import compute_ligand_receptor_interactions
from metainformant.spatial.deconvolution.spatial_deconvolution import niche_identification

LR_DATABASE = {"pairs": [{"ligand": "LIG1", "receptor": "REC1"}, {"ligand": "LIG2", "receptor": "REC2"}]}
GENE_NAMES = ["LIG1", "REC1", "LIG2", "REC2"]


def _global_state_unchanged_by(call: Callable[[], Any]) -> bool:
    np.random.seed(2024)
    before = np.random.get_state()
    call()
    after = np.random.get_state()
    return before[0] == after[0] and bool(np.array_equal(before[1], after[1])) and before[2:] == after[2:]


def _communication_inputs() -> tuple[np.ndarray, list[str]]:
    expression = np.random.default_rng(0).gamma(2.0, 1.0, size=(30, 4))
    return expression, ["A", "B", "C"] * 10


def _niche_inputs() -> tuple[np.ndarray, list[tuple[float, float]]]:
    proportions = np.random.default_rng(1).dirichlet([1.0, 1.0, 1.0], size=24)
    coordinates = [(float(i % 6), float(i // 6)) for i in range(24)]
    return proportions, coordinates


def _p_values(result: dict[str, Any]) -> list[float]:
    return [interaction["p_value"] for interaction in result["interactions"]]


class TestLigandReceptorPermutations:
    def test_the_same_seed_gives_the_same_p_values(self) -> None:
        expression, cell_types = _communication_inputs()

        first = compute_ligand_receptor_interactions(expression, cell_types, LR_DATABASE, gene_names=GENE_NAMES, seed=5)
        again = compute_ligand_receptor_interactions(expression, cell_types, LR_DATABASE, gene_names=GENE_NAMES, seed=5)

        assert _p_values(first) == _p_values(again)

    def test_different_seeds_draw_different_permutations(self) -> None:
        expression, cell_types = _communication_inputs()

        runs = {
            tuple(
                _p_values(
                    compute_ligand_receptor_interactions(
                        expression, cell_types, LR_DATABASE, gene_names=GENE_NAMES, seed=seed
                    )
                )
            )
            for seed in range(6)
        }

        assert len(runs) > 1

    def test_the_global_random_state_is_neither_read_nor_advanced(self) -> None:
        expression, cell_types = _communication_inputs()

        assert _global_state_unchanged_by(
            lambda: compute_ligand_receptor_interactions(expression, cell_types, LR_DATABASE, gene_names=GENE_NAMES)
        )
        assert _global_state_unchanged_by(
            lambda: compute_ligand_receptor_interactions(
                expression, cell_types, LR_DATABASE, gene_names=GENE_NAMES, seed=5
            )
        )


class TestNicheKMeansInitialisation:
    def test_the_same_seed_gives_the_same_niches(self) -> None:
        proportions, coordinates = _niche_inputs()

        first = niche_identification(proportions, coordinates, n_niches=3, seed=11)
        again = niche_identification(proportions, coordinates, n_niches=3, seed=11)

        np.testing.assert_array_equal(first["niche_labels"], again["niche_labels"])

    def test_the_global_random_state_is_neither_read_nor_advanced(self) -> None:
        proportions, coordinates = _niche_inputs()

        assert _global_state_unchanged_by(lambda: niche_identification(proportions, coordinates, n_niches=3))
        assert _global_state_unchanged_by(lambda: niche_identification(proportions, coordinates, n_niches=3, seed=11))
