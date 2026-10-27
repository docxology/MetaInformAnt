"""Reverse-complement overlap detection tests for long-read assembly.

All sequences are deterministic (seeded). For an antiparallel overlap every
shared canonical minimizer satisfies ``qpos + tpos == const``; the tests pin
that geometry, which the forward-only implementation could never satisfy.
"""

from __future__ import annotations

import random

from metainformant.core.sequence import reverse_complement
from metainformant.longread.assembly.overlap import _chain_minimizer_matches, find_overlaps


def _random_seq(rng: random.Random, length: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


class TestReverseComplementOverlaps:
    def test_rc_overlap_between_flanked_reads(self):
        # read_a = flank_a + core (500 bp); read_b = rc(core + down) (500 bp),
        # so rc(core) occupies B[300:500] and A's suffix core (A[300:500])
        # aligns to it antiparallel. A core k-mer starting at p appears in A
        # at 300+p and in B's forward orientation at 500-10-p, so every
        # shared canonical 10-mer obeys qpos + tpos == 790 and the chained
        # overlap must satisfy
        #   query_start + target_end == query_end + target_start == 790.
        rng = random.Random(33)
        core = _random_seq(rng, 200)
        flank_a = _random_seq(rng, 300)
        down = _random_seq(rng, 300)
        read_a = flank_a + core
        read_b = reverse_complement(core + down)

        overlaps = find_overlaps(
            [read_a, read_b],
            min_overlap=150,
            k=10,
            w=5,
            min_minimizer_matches=3,
            max_overhang=1000,
        )

        assert len(overlaps) == 1
        ov = overlaps[0]
        assert ov.strand == "-"
        assert ov.query_name == "read_0" and ov.target_name == "read_1"
        assert ov.query_start + ov.target_end == 790
        assert ov.query_end + ov.target_start == 790
        assert ov.target_end - ov.target_start == ov.query_end - ov.query_start
        assert ov.overlap_length >= 150
        assert 300 <= ov.query_start <= ov.query_end < 500  # inside A's core suffix
        assert 300 <= ov.target_start <= ov.target_end <= 500  # inside B's rc(core) suffix
        assert not ov.is_contained

    def test_fully_reverse_complemented_read_pair(self):
        # read_b = rc(read_a): the whole reads align antiparallel with
        # qpos + tpos == 500 - 10 == 490 for every shared minimizer.
        rng = random.Random(7)
        read_a = _random_seq(rng, 500)
        read_b = reverse_complement(read_a)

        overlaps = find_overlaps(
            [read_a, read_b],
            min_overlap=150,
            k=10,
            w=5,
            min_minimizer_matches=3,
            max_overhang=1000,
        )

        assert len(overlaps) == 1
        ov = overlaps[0]
        assert ov.strand == "-"
        assert ov.query_start + ov.target_end == 490
        assert ov.query_end + ov.target_start == 490
        assert ov.overlap_length >= 400  # whole-read alignment
        assert ov.is_contained

    def test_identical_reads_still_chain_on_forward_strand(self):
        # Regression guard: the rc frame must not steal same-strand overlaps.
        rng = random.Random(11)
        seq = _random_seq(rng, 500)
        overlaps = find_overlaps(
            [seq, seq],
            min_overlap=150,
            k=10,
            w=5,
            min_minimizer_matches=3,
            max_overhang=1000,
        )

        assert overlaps
        assert all(ov.strand == "+" for ov in overlaps)


class TestChainMinimizerMatchesFrames:
    """The forward and reverse-complement chaining frames and their tie/None outcomes."""

    # Reads of 400 bp / 1000 bp, k-mer size 15, minimizer matches as (query_pos, target_pos).
    _ARGS = {"query_length": 400, "target_length": 1000, "min_overlap": 20, "max_overhang": 150, "k": 15}

    def _chain(self, matches):
        return _chain_minimizer_matches("query", target_name="target", matches=matches, **self._ARGS)

    def test_no_matches_gives_no_overlap(self):
        assert self._chain([]) is None

    def test_single_match_cannot_chain_in_either_frame(self):
        assert self._chain([(100, 500)]) is None

    def test_constant_diagonal_chains_on_forward_strand(self):
        # qpos - tpos is constant (100); in the rc frame every match falls into its own band.
        ov = self._chain([(100, 0), (200, 100), (300, 200)])

        assert ov is not None
        assert ov.strand == "+"
        assert (ov.query_start, ov.query_end, ov.target_start, ov.target_end) == (100, 300, 0, 200)
        assert ov.num_matches == 3

    def test_constant_antidiagonal_chains_on_reverse_strand(self):
        # qpos + tpos is constant (900): the rc frame t' = 1000 - 15 - tpos gives a constant diagonal.
        ov = self._chain([(100, 800), (200, 700), (300, 600)])

        assert ov is not None
        assert ov.strand == "-"
        # Target coordinates are reported in the target's original (forward) orientation.
        assert (ov.query_start, ov.query_end, ov.target_start, ov.target_end) == (100, 300, 600, 800)
        assert ov.num_matches == 3

    def test_equal_chains_prefer_the_forward_strand(self):
        # These matches chain equally well (3 matches) in both frames.
        ov = self._chain([(100, 500), (110, 505), (120, 510)])

        assert ov is not None
        assert ov.strand == "+"
        assert ov.num_matches == 3
