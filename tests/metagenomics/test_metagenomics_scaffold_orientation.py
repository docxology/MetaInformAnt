"""Scaffolds are chains of oriented contigs inferred from the strands of the mates.

The scaffolder indexed both strands of every contig but threw the strand away: contigs were joined as
given, ordered by length, with gaps taken from a formula that never looked at a read, and
``insert_std`` was ignored. Here a genome with known contig boundaries is cut into contigs that are
shuffled and reverse-complemented, paired reads are simulated from it, and the scaffold must put the
contigs back in order and orientation. Real implementations; the oracle is the genome itself.
"""

from __future__ import annotations

import random

import pytest

from metainformant.metagenomics.shotgun.assembly import Contig, Scaffold, scaffold_contigs


def _revcomp(sequence: str) -> str:
    """Independent oracle for the reverse complement."""
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _random_dna(length: int, rng: random.Random) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def _fr_pairs(
    genome: str,
    count: int,
    rng: random.Random,
    insert: int = 500,
    insert_std: float = 0.0,
    read_length: int = 60,
    stop: int | None = None,
) -> dict[str, tuple[str, str]]:
    """Forward/reverse pairs from random fragments of both strands of ``genome`` (standard FR library)."""
    pairs: dict[str, tuple[str, str]] = {}
    limit = len(genome) if stop is None else stop
    for i in range(count):
        length = max(2 * read_length, round(rng.gauss(insert, insert_std))) if insert_std else insert
        start = rng.randint(0, limit - length)
        fragment = genome[start : start + length]
        if rng.random() < 0.5:
            fragment = _revcomp(fragment)
        pairs[f"p{i}"] = (fragment[:read_length], _revcomp(fragment[-read_length:]))
    return pairs


def _presented(
    genome: str, bounds: list[tuple[int, int]], flips: list[bool], order: list[int]
) -> tuple[list[Contig], list[str]]:
    """Cut ``genome`` at ``bounds``, reverse-complement the flipped pieces and list them in ``order``."""
    pieces = [genome[a:b] for a, b in bounds]
    given = [_revcomp(p) if flip else p for p, flip in zip(pieces, flips)]
    return [Contig(f"c{i}", given[i], coverage=10.0) for i in order], pieces


def _expected_sequences(pieces: list[str], gaps: list[int]) -> tuple[str, str]:
    """The scaffold sequence in genome order, and its reverse complement (the other reading direction)."""
    parts = []
    for i, piece in enumerate(pieces):
        parts.append(piece)
        if i < len(gaps):
            parts.append("N" * gaps[i])
    forward = "".join(parts)
    return forward, _revcomp(forward)


class TestOrientationIsInferred:
    BOUNDS = [(0, 800), (900, 1700), (1800, 2600), (2700, 3500), (3800, 4600)]
    TRUE_GAPS = [100, 100, 100, 300]
    FLIPS = [False, True, False, True, True]

    def _scaffold(self, min_links: int, seed: int = 5) -> tuple[list[Scaffold], list[str]]:
        rng = random.Random(seed)
        genome = _random_dna(4800, rng)
        contigs, pieces = _presented(genome, self.BOUNDS, self.FLIPS, order=[3, 0, 4, 2, 1])
        pairs = _fr_pairs(genome, 400, rng, stop=4600)
        scaffolds = scaffold_contigs(contigs, paired_reads=pairs, insert_size=500, insert_std=50, min_links=min_links)
        return scaffolds, pieces

    def test_order_orientation_and_sequence_are_recovered(self) -> None:
        scaffolds, pieces = self._scaffold(3)
        assert len(scaffolds) == 1
        scaffold = scaffolds[0]
        ids = [c.contig_id for c in scaffold.contigs]
        # contigs given flipped must be flipped back ("-"); reading the other way inverts everything
        forward_orientations = ["-" if flip else "+" for flip in self.FLIPS]
        flip_sign = {"+": "-", "-": "+"}
        assert (ids, scaffold.orientations) in [
            (["c0", "c1", "c2", "c3", "c4"], forward_orientations),
            (["c4", "c3", "c2", "c1", "c0"], [flip_sign[o] for o in forward_orientations[::-1]]),
        ]
        reading_forward = ids[0] == "c0"
        assert scaffold.gaps == (self.TRUE_GAPS if reading_forward else self.TRUE_GAPS[::-1])
        forward, reverse = _expected_sequences(pieces, self.TRUE_GAPS)
        assert scaffold.sequence == (forward if reading_forward else reverse)
        assert scaffold.total_length == len(scaffold.sequence) == 4600

    def test_every_junction_reports_its_evidence(self) -> None:
        scaffolds, _ = self._scaffold(3)
        scaffold = scaffolds[0]
        assert len(scaffold.link_counts) == len(scaffold.orientation_support) == len(scaffold.gaps) == 4
        assert all(count >= 3 for count in scaffold.link_counts)
        assert scaffold.orientation_support == scaffold.link_counts  # every pair had both mates on a definite strand

    def test_gaps_follow_the_insert_size_distribution(self) -> None:
        """With fragment lengths ~ N(500, 40) the median of the per-pair estimates stays within a read length."""
        bounds = [(0, 800), (900, 1700), (1780, 2580), (2700, 3500), (3650, 4450)]
        true_gaps = [100, 80, 120, 150]
        rng = random.Random(11)
        genome = _random_dna(4600, rng)
        contigs, _ = _presented(genome, bounds, self.FLIPS, order=[1, 3, 0, 2, 4])
        pairs = _fr_pairs(genome, 1200, rng, insert=500, insert_std=40, stop=4450)
        (scaffold,) = scaffold_contigs(contigs, paired_reads=pairs, insert_size=500, insert_std=40, min_links=5)
        expected = true_gaps if scaffold.contigs[0].contig_id == "c0" else true_gaps[::-1]
        assert scaffold.gaps == pytest.approx(expected, abs=25)

    @pytest.mark.parametrize("flip_first", [False, True])
    @pytest.mark.parametrize("flip_second", [False, True])
    @pytest.mark.parametrize("second_listed_first", [False, True])
    def test_every_pair_of_orientations_is_resolved(
        self, flip_first: bool, flip_second: bool, second_listed_first: bool
    ) -> None:
        rng = random.Random(2)
        genome = _random_dna(1800, rng)
        contigs, pieces = _presented(
            genome, [(0, 800), (900, 1700)], [flip_first, flip_second], order=[1, 0] if second_listed_first else [0, 1]
        )
        pairs = _fr_pairs(genome, 80, rng, stop=1700)
        (scaffold,) = scaffold_contigs(contigs, paired_reads=pairs, insert_size=500, insert_std=50, min_links=3)
        # applying the reported orientations must rebuild the genome (or its reverse complement)
        rebuilt = "N".join(
            (c.sequence if o == "+" else _revcomp(c.sequence)) for c, o in zip(scaffold.contigs, scaffold.orientations)
        )
        forward = pieces[0] + "N" + pieces[1]
        assert rebuilt in (forward, _revcomp(forward))
        assert scaffold.orientation_support[0] >= 3

    def test_pairs_below_min_links_do_not_join(self) -> None:
        rng = random.Random(2)
        genome = _random_dna(1800, rng)
        contigs, _ = _presented(genome, [(0, 800), (900, 1700)], [False, True], order=[0, 1])
        pairs = _fr_pairs(genome, 80, rng, stop=1700)
        joined = scaffold_contigs(contigs, paired_reads=pairs, min_links=3, insert_std=50)
        assert len(joined) == 1
        links = joined[0].link_counts[0]
        assert len(scaffold_contigs(contigs, paired_reads=pairs, min_links=links, insert_std=50)) == 1
        assert len(scaffold_contigs(contigs, paired_reads=pairs, min_links=links + 1, insert_std=50)) == 2


class TestJoinsAreChosenCarefully:
    def test_each_contig_end_is_joined_at_most_once(self) -> None:
        """A weak, spurious A-C link must not displace the well-supported A-B and B-C joins or add a third join."""
        rng = random.Random(8)
        genome = _random_dna(3000, rng)
        contigs, pieces = _presented(genome, [(0, 800), (900, 1700), (1800, 2600)], [False] * 3, order=[0, 1, 2])
        pairs = _fr_pairs(genome, 200, rng, stop=2600)
        a, c = contigs[0].sequence, contigs[2].sequence
        for i in range(6):  # chimeric pairs: forward read at the end of A, mate (reverse) at the start of C
            pairs[f"chimera{i}"] = (a[-120 + i : -60 + i], _revcomp(c[10 + i : 70 + i]))
        (scaffold,) = scaffold_contigs(contigs, paired_reads=pairs, insert_size=500, insert_std=50, min_links=3)
        assert [x.contig_id for x in scaffold.contigs] in (["c0", "c1", "c2"], ["c2", "c1", "c0"])
        assert len(scaffold.gaps) == 2

    def test_a_repeat_shared_by_two_contigs_anchors_nothing(self) -> None:
        """Reads from a segment present in both contigs used to link the contigs to each other."""
        rng = random.Random(4)
        shared = _random_dna(150, rng)
        a = Contig("a", _random_dna(600, rng) + shared, coverage=10.0)
        b = Contig("b", shared + _random_dna(600, rng), coverage=10.0)
        pairs = {f"p{i}": (shared[10 + i : 70 + i], _revcomp(a.sequence[100 + i : 160 + i])) for i in range(5)}
        scaffolds = scaffold_contigs([a, b], paired_reads=pairs, min_links=3)
        assert len(scaffolds) == 2

    def test_mates_inconsistent_with_the_insert_size_are_ignored(self) -> None:
        """Reads ~950 bp from the contig ends cannot come from a 500 bp fragment spanning a short gap."""
        rng = random.Random(6)
        a = Contig("a", _random_dna(1000, rng), coverage=10.0)
        b = Contig("b", _random_dna(1000, rng), coverage=10.0)
        pairs = {f"p{i}": (a.sequence[50 + i : 110 + i], _revcomp(b.sequence[900 + i : 960 + i])) for i in range(5)}
        assert len(scaffold_contigs([a, b], paired_reads=pairs, insert_size=500, insert_std=50, min_links=3)) == 2
        # a wide enough insert distribution makes the same pairs plausible (implied overlap of 1410 bp)
        assert len(scaffold_contigs([a, b], paired_reads=pairs, insert_size=500, insert_std=1000, min_links=3)) == 1

    def test_mates_on_one_contig_and_unplaceable_reads_are_ignored(self) -> None:
        rng = random.Random(3)
        a = Contig("a", _random_dna(600, rng), coverage=10.0)
        b = Contig("b", _random_dna(600, rng), coverage=10.0)
        junk = _random_dna(60, rng)
        pairs = {
            "same": (a.sequence[0:60], _revcomp(a.sequence[300:360])),
            "none": (junk, _revcomp(junk)),
            "short": ("ACGT", "ACGT"),
        }
        assert len(scaffold_contigs([a, b], paired_reads=pairs, min_links=1)) == 2


class TestLimitationIsReported:
    """Contigs that are their own reverse complement cannot show their orientation; the scaffold says so."""

    def test_palindromic_contigs_are_joined_as_given_with_no_orientation_support(self) -> None:
        c1 = Contig("c1", "ATCG" * 250, coverage=10.0)
        c2 = Contig("c2", "GCTA" * 125, coverage=20.0)
        pairs = {"p1": (c1.sequence, c2.sequence), "p2": (c1.sequence, c2.sequence)}
        (scaffold,) = scaffold_contigs([c1, c2], paired_reads=pairs, min_links=1)
        assert [c.contig_id for c in scaffold.contigs] == ["c1", "c2"]
        assert scaffold.orientations == ["+", "+"]  # as given: nothing contradicts it, nothing supports it
        assert scaffold.link_counts == [2]
        assert scaffold.orientation_support == [0]
        assert scaffold.gaps == [10]  # no read position: the minimum gap, not an estimate
        assert scaffold.sequence == c1.sequence + "N" * 10 + c2.sequence

    def test_unjoined_contigs_keep_their_orientation_and_sequence(self) -> None:
        contigs = [Contig("c1", "ATCG" * 250), Contig("c2", "GCTA" * 125)]
        scaffolds = scaffold_contigs(contigs)
        assert [s.orientations for s in scaffolds] == [["+"], ["+"]]
        assert [s.sequence for s in scaffolds] == [c.sequence for c in contigs]
        assert all(s.link_counts == [] and s.orientation_support == [] for s in scaffolds)
