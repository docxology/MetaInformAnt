"""Tests for metainformant.spatial.communication -- cell-cell communication analysis.

All tests use real implementations (NO mocking policy).
"""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.spatial.communication.cell_communication import (
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
