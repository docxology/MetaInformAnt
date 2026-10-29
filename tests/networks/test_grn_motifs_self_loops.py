"""Regulatory motifs span distinct nodes: autoregulation and repeated edges are not motif edges.

``compute_network_motifs`` counted a self-loop as an edge of a 3-node motif. ``A->A, A->B`` (two nodes) contained one
feed-forward loop (a = b = A), ``A->B, B->B`` one loop and one co-regulation pair (B regulates itself, so B had "two
regulators"), and a genuine feed-forward loop whose three TFs also regulate themselves was counted 7 times with 4
co-regulated pairs instead of 1 and 1. The 2-node ``single_regulation`` count included each self-loop. Repeated edges
inflated the co-regulation pairs of the observed network. When the randomized counts did not vary at all, the z-score
was reported as 0.0 whatever the observed count was.

The oracles are the hand-counted motifs of the small edge lists below. No test doubles.
"""

from __future__ import annotations

import math

import pytest

from metainformant.networks.regulatory import grn_inference
from metainformant.networks.regulatory.grn_inference import compute_network_motifs


def _counts(pairs: list[tuple[str, str]], size: int = 3) -> dict[str, int]:
    edges = [{"source": s, "target": t} for s, t in pairs]
    return compute_network_motifs(edges, motif_size=size, random_state=1)["motif_counts"]


ZERO = {"feed_forward_loop": 0, "cascade": 0, "co_regulation": 0}


@pytest.mark.parametrize(
    ("pairs", "expected"),
    [
        ([("A", "A"), ("A", "B")], ZERO),  # two nodes cannot form a 3-node motif
        ([("A", "A")], ZERO),
        ([("A", "B"), ("B", "B")], ZERO),
        ([("A", "B"), ("B", "C"), ("C", "C")], {"feed_forward_loop": 0, "cascade": 1, "co_regulation": 0}),
        ([("A", "B"), ("B", "C"), ("A", "C")], {"feed_forward_loop": 1, "cascade": 0, "co_regulation": 1}),
        (
            [("A", "B"), ("B", "C"), ("A", "C"), ("A", "A"), ("B", "B"), ("C", "C")],
            {"feed_forward_loop": 1, "cascade": 0, "co_regulation": 1},
        ),
    ],
)
def test_three_node_motifs_ignore_autoregulation(pairs: list[tuple[str, str]], expected: dict[str, int]) -> None:
    assert _counts(pairs) == expected


def test_repeated_edges_do_not_inflate_co_regulation() -> None:
    # C has two distinct regulators (A, B), however often A->C is listed: one pair, not C(3, 2) = 3
    assert _counts([("A", "C"), ("A", "C"), ("A", "C"), ("B", "C")])["co_regulation"] == 1


def test_a_hub_regulated_by_four_tfs_has_six_co_regulated_pairs() -> None:
    pairs = [(f"T{i}", "HUB") for i in range(4)]

    assert _counts(pairs)["co_regulation"] == math.comb(4, 2)


def test_two_node_motifs_leave_autoregulation_out_of_the_single_regulations() -> None:
    assert _counts([("A", "A"), ("A", "B"), ("B", "A")], size=2) == {"mutual_regulation": 1, "single_regulation": 0}
    assert _counts([("A", "A"), ("A", "B"), ("C", "B")], size=2) == {"mutual_regulation": 0, "single_regulation": 2}


def test_the_same_seed_gives_the_same_z_scores() -> None:
    pairs = [(f"T{i}", f"G{(i * 3 + j) % 9}") for i in range(5) for j in range(3)] + [("T0", "T1"), ("T1", "T2")]
    edges = [{"source": s, "target": t} for s, t in pairs]

    first = compute_network_motifs(edges, motif_size=3, random_state=11)
    second = compute_network_motifs(edges, motif_size=3, random_state=11)

    assert first == second


class TestZScores:
    def test_a_constant_null_that_matches_the_observed_count_is_zero(self) -> None:
        assert grn_inference._motif_z_scores({"ffl": 4}, {"ffl": [4] * 100}) == {"ffl": 0.0}

    def test_a_constant_null_below_or_above_the_observed_count_is_infinite(self) -> None:
        z = grn_inference._motif_z_scores({"up": 5, "down": 0}, {"up": [2] * 100, "down": [3] * 100})

        assert z["up"] == math.inf and z["down"] == -math.inf  # the old code reported 0.0 for both

    def test_an_ordinary_null_gives_the_usual_z_score(self) -> None:
        null = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        mean = sum(null) / len(null)
        std = math.sqrt(sum((v - mean) ** 2 for v in null) / (len(null) - 1))

        assert grn_inference._motif_z_scores({"m": 12}, {"m": null})["m"] == pytest.approx(round((12 - mean) / std, 4))

    def test_a_network_the_shuffle_cannot_change_has_zero_z_scores(self) -> None:
        # one edge: shuffling its target returns the same network every time, and the observed count equals it
        result = compute_network_motifs([{"source": "A", "target": "B"}], motif_size=2, random_state=0)

        assert result["motif_counts"] == {"mutual_regulation": 0, "single_regulation": 1}
        assert result["z_scores"] == {"mutual_regulation": 0.0, "single_regulation": 0.0}
        assert result["significant_motifs"] == []
