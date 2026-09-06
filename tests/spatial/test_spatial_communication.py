"""Tests for metainformant.spatial.communication -- cell-cell communication analysis.

All tests use real implementations (real-implementation policy).
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from metainformant.spatial.communication import cell_communication as cell_communication_module
from metainformant.spatial.communication.cell_communication import (
    _gene_name_to_index,
    build_communication_network,
    communication_pattern_analysis,
    compute_ligand_receptor_interactions,
    default_lr_database,
    spatial_interaction_score,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def expression_data() -> np.ndarray:
    """30 cells x 10 genes with some structure."""
    rng = np.random.RandomState(42)
    return np.abs(rng.standard_normal((30, 10)))


@pytest.fixture()
def cell_types() -> list[str]:
    return ["A"] * 10 + ["B"] * 10 + ["C"] * 10


@pytest.fixture()
def coordinates() -> list[tuple[float, float]]:
    rng = np.random.RandomState(42)
    pts = rng.uniform(0, 100, (30, 2))
    return [(float(pts[i, 0]), float(pts[i, 1])) for i in range(30)]


@pytest.fixture()
def simple_lr_database() -> dict[str, list[dict[str, str]]]:
    """Small LR database using gene indices as names so they map into 10 genes."""
    return {
        "pairs": [
            {"ligand": "0", "receptor": "1"},
            {"ligand": "2", "receptor": "3"},
        ]
    }


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestDefaultLRDatabase:
    def test_returns_pairs(self) -> None:
        db = default_lr_database()
        assert "pairs" in db
        pairs = db["pairs"]
        assert len(pairs) > 50  # the built-in database has ~90 pairs
        assert all("ligand" in p and "receptor" in p for p in pairs)


class TestComputeLigandReceptorInteractions:
    def test_basic_operation(
        self,
        expression_data: np.ndarray,
        cell_types: list[str],
        simple_lr_database: dict[str, list[dict[str, str]]],
    ) -> None:
        result = compute_ligand_receptor_interactions(
            expression_data,
            cell_types,
            lr_database=simple_lr_database,
        )
        assert "interactions" in result
        assert "n_significant" in result
        assert "summary" in result
        # 2 LR pairs x 3 source types x 3 target types = 18 interactions
        assert len(result["interactions"]) == 18
        # Each interaction should have required keys
        for ix in result["interactions"]:
            assert "ligand" in ix
            assert "receptor" in ix
            assert "source_type" in ix
            assert "target_type" in ix
            assert "score" in ix
            assert "p_value" in ix

    def test_with_default_database(
        self,
        expression_data: np.ndarray,
        cell_types: list[str],
    ) -> None:
        result = compute_ligand_receptor_interactions(expression_data, cell_types)
        assert "interactions" in result
        assert isinstance(result["n_significant"], int)


class TestSpatialInteractionScore:
    def test_returns_scores(
        self,
        expression_data: np.ndarray,
        coordinates: list[tuple[float, float]],
    ) -> None:
        lr_pairs = [
            {"ligand_idx": 0, "receptor_idx": 1},
            {"ligand_idx": 2, "receptor_idx": 3},
        ]
        result = spatial_interaction_score(
            expression_data,
            coordinates,
            lr_pairs,
            max_distance=50.0,
        )
        assert "spatial_scores" in result
        assert "distance_decay" in result
        assert "significant_pairs" in result
        assert len(result["spatial_scores"]) == 2
        for s in result["spatial_scores"]:
            assert "score" in s
            assert "n_interacting_pairs" in s


class TestBuildCommunicationNetwork:
    def test_creates_network(
        self,
        expression_data: np.ndarray,
        cell_types: list[str],
        simple_lr_database: dict[str, list[dict[str, str]]],
    ) -> None:
        lr_result = compute_ligand_receptor_interactions(expression_data, cell_types, lr_database=simple_lr_database)
        network = build_communication_network(lr_result["interactions"])
        assert "adjacency_matrix" in network
        assert "cell_types" in network
        assert "edge_list" in network
        assert "hub_types" in network
        assert "pathway_summary" in network
        assert len(network["cell_types"]) == 3
        # adjacency_matrix should be 3x3
        assert len(network["adjacency_matrix"]) == 3
        assert len(network["adjacency_matrix"][0]) == 3

    def test_empty_interactions_raises(self) -> None:
        with pytest.raises(ValueError, match="No interactions"):
            build_communication_network([])


class TestCommunicationPatternAnalysis:
    def test_identifies_patterns(
        self,
        expression_data: np.ndarray,
        cell_types: list[str],
        simple_lr_database: dict[str, list[dict[str, str]]],
    ) -> None:
        lr_result = compute_ligand_receptor_interactions(expression_data, cell_types, lr_database=simple_lr_database)
        patterns = communication_pattern_analysis(lr_result, n_patterns=2)
        assert "patterns" in patterns
        assert "pattern_loadings" in patterns
        assert "dominant_pathways_per_pattern" in patterns
        # Should have at most 2 patterns
        assert len(patterns["patterns"]) <= 2

    def test_empty_interactions_returns_empty(self) -> None:
        empty_result: dict[str, list[dict[str, str]]] = {"interactions": []}
        patterns = communication_pattern_analysis(empty_result, n_patterns=3)
        assert patterns["patterns"] == []


class TestGeneNameToIndex:
    def test_numeric_string_maps_directly(self) -> None:
        assert _gene_name_to_index("3", 10) == 3

    def test_out_of_range_numeric_returns_none(self) -> None:
        assert _gene_name_to_index("10", 10) is None
        assert _gene_name_to_index("-1", 10) is None

    def test_named_gene_maps_to_valid_index(self) -> None:
        idx = _gene_name_to_index("CXCL12", 97)
        assert 0 <= idx < 97

    def test_named_gene_stable_across_hash_seeds(self) -> None:
        # Regression: gene-name indices must not depend on the interpreter
        # hash seed (str hashing is randomized per process). The module file
        # is loaded directly to keep the child interpreter lightweight.
        module_path = Path(cell_communication_module.__file__).resolve()
        code = (
            "import importlib.util\n"
            f"spec = importlib.util.spec_from_file_location('cell_communication_under_test', {str(module_path)!r})\n"
            "mod = importlib.util.module_from_spec(spec)\n"
            "spec.loader.exec_module(mod)\n"
            "print(mod._gene_name_to_index('CXCL12', 97))\n"
        )
        outputs = set()
        for seed in ("0", "12345"):
            proc = subprocess.run(
                [sys.executable, "-c", code],
                capture_output=True,
                text=True,
                check=True,
                env={**os.environ, "PYTHONHASHSEED": seed},
            )
            outputs.add(proc.stdout.strip())
        assert len(outputs) == 1


class TestSeededPermutations:
    def test_seeded_permutations_are_reproducible(
        self,
        expression_data: np.ndarray,
        cell_types: list[str],
        simple_lr_database: dict[str, list[dict[str, str]]],
    ) -> None:
        r1 = compute_ligand_receptor_interactions(expression_data, cell_types, lr_database=simple_lr_database, seed=42)
        r2 = compute_ligand_receptor_interactions(expression_data, cell_types, lr_database=simple_lr_database, seed=42)
        assert r1["interactions"] == r2["interactions"]


class TestBuildCommunicationNetworkMinScore:
    def test_min_score_filters_edges(self) -> None:
        interactions = [
            {"source_type": "A", "target_type": "B", "score": 5.0, "ligand": "L1", "receptor": "R1"},
            {"source_type": "A", "target_type": "C", "score": 0.5, "ligand": "L2", "receptor": "R2"},
        ]
        network = build_communication_network(interactions, min_score=1.0)
        edges = {(e["source"], e["target"]) for e in network["edge_list"]}
        assert ("A", "B") in edges
        assert ("A", "C") not in edges
        assert network["pathway_summary"] == {"L1-R1": {"count": 1, "total_score": 5.0}}
        assert network["hub_types"][0] == "A"
