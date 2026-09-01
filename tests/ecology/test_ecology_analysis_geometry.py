"""Tests for untested ecology analysis/functional geometry helpers.

Covers:
    - convex_hull_2d (Andrew's monotone chain)
    - minimum_spanning_tree (Prim's algorithm)

Uses real implementations only (real-implementation policy).
"""

from __future__ import annotations

from metainformant.ecology.analysis.functional import (
    convex_hull_2d,
    minimum_spanning_tree,
)


class TestConvexHull2D:
    """convex_hull_2d: monotone-chain convex hull."""

    def test_square_hull_drops_interior_and_edge_points(self) -> None:
        # Interior point (0.5, 0.5) and edge point (0.5, 0.0) must not appear.
        points = [(0, 0), (1, 0), (1, 1), (0, 1), (0.5, 0.5), (0.5, 0.0)]
        hull = convex_hull_2d(points)
        assert len(hull) == 4
        assert set(hull) == {(0, 0), (1, 0), (1, 1), (0, 1)}

    def test_single_point_and_empty(self) -> None:
        assert convex_hull_2d([]) == []
        assert convex_hull_2d([(2.0, 3.0)]) == [(2.0, 3.0)]

    def test_collinear_points_collapse_to_endpoints(self) -> None:
        # All collinear: only extreme endpoints survive (strict hull).
        hull = convex_hull_2d([(0, 0), (1, 1), (2, 2), (3, 3)])
        assert set(hull) == {(0, 0), (3, 3)}

    def test_counter_clockwise_order(self) -> None:
        hull = convex_hull_2d([(0, 0), (2, 0), (2, 2), (0, 2)])
        # CCW order starting from lexicographically smallest point.
        assert hull == [(0, 0), (2, 0), (2, 2), (0, 2)]


class TestMinimumSpanningTree:
    """minimum_spanning_tree: Prim's algorithm over a symmetric matrix."""

    def test_triangle_picks_two_cheapest_edges(self) -> None:
        dist = [[0, 1, 3], [1, 0, 2], [3, 2, 0]]
        edges = minimum_spanning_tree(dist)
        total = sum(w for _, _, w in edges)
        assert total == 3
        assert len(edges) == 2

    def test_single_node_and_empty(self) -> None:
        assert minimum_spanning_tree([[0]]) == []
        assert minimum_spanning_tree([]) == []

    def test_four_nodes_known_optimum(self) -> None:
        # Optimum uses the three unit edges (weight 3).
        inf = float("inf")
        dist = [
            [0, 1, 5, 2],
            [1, 0, 1, inf],
            [5, 1, 0, 1],
            [2, inf, 1, 0],
        ]
        edges = minimum_spanning_tree(dist)
        assert sum(w for _, _, w in edges) == 3
        # MST on 4 nodes has exactly 3 edges and connects all vertices.
        assert len(edges) == 3
        nodes = set()
        for i, j, _ in edges:
            nodes.add(i)
            nodes.add(j)
        assert nodes == {0, 1, 2, 3}

    def test_star_topology_all_edges_from_center(self) -> None:
        dist = [
            [0, 2, 2, 2],
            [2, 0, 3, 3],
            [2, 3, 0, 3],
            [2, 3, 3, 0],
        ]
        edges = minimum_spanning_tree(dist)
        # Every greedy edge in this construction has weight 2 (star from node 0).
        assert all(w == 2 for _, _, w in edges)
        assert len(edges) == 3
