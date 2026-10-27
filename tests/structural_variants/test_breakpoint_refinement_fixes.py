"""Regression tests for ``refine_breakpoints``.

Each test builds real read dictionaries (CIGAR strings, sequences, SA tags) and checks the refined
breakpoints against values worked out by hand, never against the implementation:

* A call turned into a dict with ``dataclasses.asdict`` carries ``chrom2=""`` for an intra-chromosomal
  variant. ``refine_breakpoints`` kept that empty string, so breakpoint 2 was looked up on a chromosome
  called ``""`` and always came back unsupported on ``bp2.chrom == ""``. A translocation's second
  breakpoint is its ``pos2``, not its ``end``.
* Reads were prefiltered on where the alignment *starts* (window + 500 bp), so a read that reaches the
  breakpoint from further away (a long alignment, or a long deletion in its CIGAR) was never seen.
* Microhomology was only ever reported on breakpoint 1, from the first read that carried a sequence,
  wherever that read was. It is now reported on both breakpoints from the reads that support them,
  centred on their own junction.
* Strand was always ``'+'`` and ``inserted_sequence`` was never set.
"""

from __future__ import annotations

import random
from dataclasses import asdict
from typing import Any

from metainformant.structural_variants.detection.breakpoints import refine_breakpoints
from metainformant.structural_variants.detection.sv_calling import StructuralVariant, SVEvidence, SVType


def _bases(seed: int, n: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


# ---------------------------------------------------------------------------
# Variant records: empty chrom2, pos2
# ---------------------------------------------------------------------------


class TestVariantRecords:
    READS = [
        {"chrom": "chr1", "pos": 4950, "cigar": "50M100S", "name": "r1"},
        {"chrom": "chr1", "pos": 4950, "cigar": "50M100S", "name": "r2"},
        {"chrom": "chr1", "pos": 7950, "cigar": "100S50M", "name": "r3"},
    ]

    def test_empty_chrom2_from_asdict_falls_back_to_the_first_chromosome(self) -> None:
        call = StructuralVariant(
            chrom="chr1",
            start=5000,
            end=8000,
            sv_type=SVType.DEL,
            evidence=SVEvidence(split_reads=3, breakpoint1=5000, breakpoint2=8000, chrom1="chr1", chrom2="chr1"),
        )
        record = asdict(call)
        assert record["chrom2"] == ""  # the premise: an intra-chromosomal call has no second chromosome

        pair = refine_breakpoints([record], self.READS)[0]

        assert pair.bp2.chrom == "chr1"
        assert (pair.bp2.position, pair.bp2.support, pair.bp2.read_names) == (7950, 1, ["r3"])
        assert pair.sv_type == "DEL"

    def test_none_chrom2_also_falls_back(self) -> None:
        pair = refine_breakpoints([{"chrom": "chr1", "start": 5000, "end": 8000, "chrom2": None}], self.READS)[0]
        assert (pair.bp2.chrom, pair.bp2.position) == ("chr1", 7950)

    def test_translocation_second_breakpoint_is_pos2_on_chrom2(self) -> None:
        variant = {"chrom": "chr1", "start": 5000, "end": 5000, "sv_type": "TRA", "chrom2": "chr2", "pos2": 3000}
        reads = [
            {"chrom": "chr1", "pos": 4950, "cigar": "50M100S", "name": "a"},
            {"chrom": "chr1", "pos": 4950, "cigar": "50M100S", "name": "b"},
            {"chrom": "chr2", "pos": 3010, "cigar": "100S50M", "name": "c"},
            {"chrom": "chr2", "pos": 3010, "cigar": "100S50M", "name": "d"},
            # Clips near 3010 on the wrong chromosome must not count for breakpoint 2.
            {"chrom": "chr1", "pos": 3010, "cigar": "100S50M", "name": "wrong_chrom"},
        ]

        pair = refine_breakpoints([variant], reads)[0]

        assert (pair.bp1.chrom, pair.bp1.position, pair.bp1.support) == ("chr1", 5000, 2)
        assert (pair.bp2.chrom, pair.bp2.position, pair.bp2.support) == ("chr2", 3010, 2)
        assert pair.bp2.read_names == ["c", "d"]

    def test_translocation_without_a_clip_keeps_pos2(self) -> None:
        variant = {"chrom": "chr1", "start": 5000, "end": 5000, "chrom2": "chr2", "pos2": 3000}
        pair = refine_breakpoints([variant], [])[0]
        assert (pair.bp2.chrom, pair.bp2.position, pair.bp2.support) == ("chr2", 3000, 0)

    def test_end_is_still_the_second_position_of_a_hand_built_translocation(self) -> None:
        # VCF-style record: END holds the position on CHR2 and there is no pos2.
        variant = {"chrom": "chr1", "start": 5000, "end": 3000, "chrom2": "chr2"}
        reads = [{"chrom": "chr2", "pos": 3010, "cigar": "100S50M", "name": "c"}]
        pair = refine_breakpoints([variant], reads)[0]
        assert (pair.bp2.chrom, pair.bp2.position, pair.bp2.support) == ("chr2", 3010, 1)

    def test_explicit_breakpoint_keys_take_precedence(self) -> None:
        variant = {"chrom1": "chr1", "breakpoint1": 5000, "chrom2": "chr1", "breakpoint2": 8000, "pos2": 1, "end": 2}
        pair = refine_breakpoints([variant], self.READS)[0]
        assert (pair.bp1.position, pair.bp2.position) == (5000, 7950)

    def test_cigar_as_operation_tuples(self) -> None:
        reads = [
            {"chrom": "chr1", "pos": 4950, "cigar": [(0, 50), (4, 100)], "name": "t1"},
            {"chrom": "chr1", "pos": 4950, "cigar": [(0, 50), (4, 100)], "name": "t2"},
        ]
        pair = refine_breakpoints([{"chrom": "chr1", "start": 5010, "end": 9000}], reads)[0]
        assert (pair.bp1.position, pair.bp1.support) == (5000, 2)


# ---------------------------------------------------------------------------
# Reads that reach the breakpoint from far away
# ---------------------------------------------------------------------------


class TestReadsFromFarAway:
    def test_long_alignments_clipping_at_the_breakpoint_are_used(self) -> None:
        reads = [
            {"chrom": "chr1", "pos": 3000, "cigar": "2000M100S", "name": "long1"},  # clip at 5000, starts 2 kb away
            {"chrom": "chr1", "pos": 2900, "cigar": "2100M100S", "name": "long2"},  # clip at 5000, starts 2.1 kb away
        ]
        pair = refine_breakpoints([{"chrom": "chr1", "start": 5030, "end": 9000}], reads)[0]
        assert (pair.bp1.position, pair.bp1.support) == (5000, 2)

    def test_a_long_deletion_in_the_cigar_moves_the_clip_away_from_the_read_start(self) -> None:
        # The right clip lies at 1000 + 100 (M) + 5000 (D) + 100 (M) = 6200, 5.2 kb from the read start.
        reads = [{"chrom": "chr1", "pos": 1000, "cigar": "100M5000D100M60S", "name": f"del{i}"} for i in range(3)]
        pair = refine_breakpoints([{"chrom": "chr1", "start": 1100, "end": 6190}], reads)[0]
        assert (pair.bp2.position, pair.bp2.support) == (6200, 3)
        assert pair.bp1.support == 0  # the edge of the D operation is not a soft clip

    def test_reads_are_still_judged_by_the_clip_position(self) -> None:
        # A read that starts next to the breakpoint but clips far outside the window does not count.
        reads = [{"chrom": "chr1", "pos": 4990, "cigar": "1000M100S", "name": "far_clip"}]  # clip at 5990
        pair = refine_breakpoints([{"chrom": "chr1", "start": 5000, "end": 9000}], reads)[0]
        assert (pair.bp1.position, pair.bp1.support) == (5000, 0)

    def test_window_bounds_are_inclusive(self) -> None:
        reads = [{"chrom": "chr1", "pos": 4800, "cigar": "200M50S", "name": "edge"}]  # clip at exactly 5000
        inside = refine_breakpoints([{"chrom": "chr1", "start": 5200, "end": 9000}], reads, window=200)[0]
        outside = refine_breakpoints([{"chrom": "chr1", "start": 5201, "end": 9000}], reads, window=200)[0]
        assert inside.bp1.support == 1
        assert outside.bp1.support == 0


# ---------------------------------------------------------------------------
# Microhomology on both breakpoints, from the supporting reads
# ---------------------------------------------------------------------------

# Junction at offset 20 of a 40-base read. detect_microhomology compares the base k before the junction with
# the base k after it until they differ. For this read (indices 17-19 | 20-23):
#   G|G, C|C, A|A, then T|C differs -> "GCA".
BP1_LEFT, BP1_RIGHT = "C" * 16 + "TACG", "GCA" + "C" + "A" * 16
# This one (17-19 | 20-22): T|T, G|G, then A|C differs -> "TG".
BP2_LEFT, BP2_RIGHT = "A" * 17 + "AGT", "TGC" + "C" * 17
# And this one: G|G, then C|A differs -> "G".
SHORT_LEFT, SHORT_RIGHT = "C" * 19 + "G", "G" + "A" * 19
# A read that is unrelated to both breakpoints but whose own midpoint mirrors: "ACGT|TGCA" -> "TGCA".
DECOY_SEQ = "ACGTTGCA"


def _microhomology_reads(with_sequences: bool = True) -> list[dict[str, Any]]:
    def read(name: str, pos: int, cigar: str, seq: str) -> dict[str, Any]:
        record: dict[str, Any] = {"chrom": "chr1", "pos": pos, "cigar": cigar, "name": name}
        if with_sequences:
            record["seq"] = seq
        return record

    return [
        # The first read that carries a sequence, far from both breakpoints.
        {"chrom": "chr1", "pos": 100, "cigar": "8M", "name": "decoy", "seq": DECOY_SEQ},
        read("short", 4980, "20M20S", SHORT_LEFT + SHORT_RIGHT),  # right clip at 5000, homology "G"
        read("bp1", 4980, "20M20S", BP1_LEFT + BP1_RIGHT),  # right clip at 5000, homology "GCA"
        read("bp2", 7000, "20S20M", BP2_LEFT + BP2_RIGHT),  # left clip at 7000, homology "TG"
    ]


class TestMicrohomology:
    VARIANT = {"chrom": "chr1", "start": 5000, "end": 7000, "sv_type": "DEL"}

    def test_both_breakpoints_get_the_microhomology_of_their_own_supporting_reads(self) -> None:
        pair = refine_breakpoints([self.VARIANT], _microhomology_reads())[0]

        assert (pair.bp1.position, pair.bp1.support) == (5000, 2)
        assert (pair.bp2.position, pair.bp2.support) == (7000, 1)
        assert pair.bp1.microhomology == "GCA"  # the longest of the two supporting reads' ("G", "GCA")
        assert pair.bp2.microhomology == "TG"

    def test_a_read_that_does_not_support_the_breakpoint_is_not_used(self) -> None:
        pair = refine_breakpoints([self.VARIANT], _microhomology_reads())[0]
        assert "TGCA" not in (pair.bp1.microhomology, pair.bp2.microhomology)  # the decoy's own mirror

    def test_no_sequence_on_a_supporting_read_means_no_microhomology(self) -> None:
        reads = _microhomology_reads(with_sequences=False)
        reads.insert(0, {"chrom": "chr1", "pos": 90, "cigar": "8M", "name": "decoy2", "seq": DECOY_SEQ})
        pair = refine_breakpoints([self.VARIANT], reads)[0]
        assert (pair.bp1.microhomology, pair.bp2.microhomology) == ("", "")

    def test_a_supporting_read_whose_sequence_is_none_is_skipped(self) -> None:
        reads = _microhomology_reads()
        reads[2]["seq"] = None  # the read that would give breakpoint 1 its "GCA"
        pair = refine_breakpoints([self.VARIANT], reads)[0]
        assert pair.bp1.microhomology == "G"  # the other supporting read's
        assert pair.bp2.microhomology == "TG"

    def test_only_the_window_around_the_junction_is_examined(self) -> None:
        # A long read whose junction is not at its middle: the junction-centred window still finds "GCA".
        left = "T" * 80 + BP1_LEFT  # 100 bases, junction after the 100th
        right = BP1_RIGHT + "G" * 30  # 50 bases
        reads = [{"chrom": "chr1", "pos": 4900, "cigar": "100M50S", "name": "long", "seq": left + right}]
        pair = refine_breakpoints([self.VARIANT], reads)[0]
        assert pair.bp1.position == 5000
        assert pair.bp1.microhomology == "GCA"


# ---------------------------------------------------------------------------
# Strand
# ---------------------------------------------------------------------------


class TestStrand:
    VARIANT = {"chrom": "chr1", "start": 5000, "end": 8000}

    @staticmethod
    def _read(name: str, pos: int, cigar: str, is_reverse: bool | None) -> dict[str, Any]:
        read: dict[str, Any] = {"chrom": "chr1", "pos": pos, "cigar": cigar, "name": name}
        if is_reverse is not None:
            read["is_reverse"] = is_reverse
        return read

    def test_strand_is_the_majority_orientation_of_the_supporting_reads(self) -> None:
        reads = [
            self._read("f1", 4950, "50M100S", False),
            self._read("r1", 4950, "50M100S", True),
            self._read("r2", 4950, "50M100S", True),
            self._read("f2", 7950, "100S50M", False),
            self._read("f3", 7950, "100S50M", False),
        ]
        pair = refine_breakpoints([self.VARIANT], reads)[0]
        assert pair.bp1.strand == "-"
        assert pair.bp2.strand == "+"

    def test_only_reads_at_the_consensus_position_vote(self) -> None:
        reads = [self._read(f"f{i}", 4950, "50M100S", False) for i in range(3)]
        reads += [self._read(f"r{i}", 4960, "50M100S", True) for i in range(2)]  # clip at 5010, not the consensus
        pair = refine_breakpoints([self.VARIANT], reads)[0]
        assert (pair.bp1.position, pair.bp1.support, pair.bp1.strand) == (5000, 3, "+")

    def test_without_orientation_the_default_is_plus(self) -> None:
        reads = [self._read("a", 4950, "50M100S", None), self._read("b", 4950, "50M100S", None)]
        pair = refine_breakpoints([self.VARIANT], reads)[0]
        assert pair.bp1.strand == "+"
        assert refine_breakpoints([self.VARIANT], [])[0].bp1.strand == "+"

    def test_a_tie_keeps_the_default(self) -> None:
        reads = [self._read("a", 4950, "50M100S", True), self._read("b", 4950, "50M100S", False)]
        assert refine_breakpoints([self.VARIANT], reads)[0].bp1.strand == "+"

    def test_reads_without_orientation_do_not_vote(self) -> None:
        reads = [
            self._read("a", 4950, "50M100S", True),
            self._read("b", 4950, "50M100S", None),
            self._read("c", 4950, "50M100S", None),
        ]
        assert refine_breakpoints([self.VARIANT], reads)[0].bp1.strand == "-"


# ---------------------------------------------------------------------------
# Inserted sequence
# ---------------------------------------------------------------------------


class TestInsertedSequence:
    READ = _bases(1, 100)

    def _split_pair(self, **overrides: Any) -> list[dict[str, Any]]:
        """A read whose first 40 bases align at chr1:1000 and last 30 at chr1:3000, 30 bases inserted between."""
        primary = {
            "name": "split1",
            "chrom": "chr1",
            "pos": 1000,
            "cigar": "40M60S",
            "seq": self.READ,
            "is_reverse": False,
            "sa_tag": "chr1,3001,+,70S30M,60,0;",
        }
        supplementary = {
            "name": "split1",
            "chrom": "chr1",
            "pos": 3000,
            "cigar": "70S30M",
            "seq": self.READ,
            "is_reverse": False,
            "sa_tag": "chr1,1001,+,40M60S,60,0;",
        }
        primary.update(overrides)
        return [primary, supplementary]

    VARIANT = {"chrom": "chr1", "start": 1040, "end": 3000, "sv_type": "DEL"}

    def test_bases_between_the_two_alignments_are_the_inserted_sequence(self) -> None:
        pair = refine_breakpoints([self.VARIANT], self._split_pair())[0]

        assert (pair.bp1.position, pair.bp2.position) == (1040, 3000)
        assert pair.bp1.inserted_sequence == self.READ[40:70]
        assert pair.bp2.inserted_sequence == self.READ[40:70]
        assert len(pair.bp1.inserted_sequence) == 30

    def test_opposite_strand_alignments_are_compared_in_one_orientation(self) -> None:
        # The original read O aligns O[0:40] to the reverse strand at chr1:1000 and O[70:100] to the forward
        # strand at chr1:3000, with O[40:70] inserted. The reverse-strand record stores revcomp(O).
        original = _bases(2, 100)
        primary = {
            "name": "inv1",
            "chrom": "chr1",
            "pos": 1000,
            "cigar": "60S40M",
            "seq": _revcomp(original),
            "is_reverse": True,
            "sa_tag": "chr1,3001,+,70S30M,60,0;",
        }
        supplementary = {
            "name": "inv1",
            "chrom": "chr1",
            "pos": 3000,
            "cigar": "70S30M",
            "seq": original,
            "is_reverse": False,
            "sa_tag": "chr1,1001,-,60S40M,60,0;",
        }
        variant = {"chrom": "chr1", "start": 1000, "end": 3000, "sv_type": "INV"}

        pair = refine_breakpoints([variant], [primary, supplementary])[0]

        # Each breakpoint reports the bases as they read along the forward strand of its own alignment.
        assert pair.bp1.inserted_sequence == _revcomp(original[40:70])
        assert pair.bp2.inserted_sequence == original[40:70]

    def test_the_most_common_sequence_among_the_supporting_reads_wins(self) -> None:
        other = _bases(3, 100)
        reads = [dict(self._split_pair()[0], name=f"a{i}") for i in range(2)]
        reads.append(dict(self._split_pair()[0], name="b", seq=other))
        pair = refine_breakpoints([self.VARIANT], reads)[0]
        assert pair.bp1.support == 3
        assert pair.bp1.inserted_sequence == self.READ[40:70]

    def test_nothing_is_reported_when_it_cannot_be_derived(self) -> None:
        def inserted(**overrides: Any) -> str:
            pair = refine_breakpoints([self.VARIANT], self._split_pair(**overrides))[0]
            assert pair.bp1.support == 1  # the read still supports the breakpoint
            return pair.bp1.inserted_sequence

        assert inserted() != ""  # the baseline does derive it
        assert inserted(sa_tag="") == ""  # no supplementary alignment
        assert inserted(sa_tag="chr1,3001,+,70S30M,60,0;chr1,5001,+,10S30M60S,60,0;") == ""  # more than one
        assert inserted(sa_tag="chr1,9001,+,70S30M,60,0;") == ""  # not at the other breakpoint
        assert inserted(sa_tag="chr2,3001,+,70S30M,60,0;") == ""  # on another chromosome
        assert inserted(sa_tag="chr1,3001,+,30S70M,60,0;") == ""  # the alignments overlap (microhomology)
        assert inserted(sa_tag="chr1,3001,+,40S60M,60,0;") == ""  # the alignments abut: nothing inserted
        assert inserted(sa_tag="chr1,3001,+") == ""  # an SA entry without a CIGAR
        assert inserted(seq=self.READ[:90]) == ""  # the CIGARs do not describe this sequence
        assert inserted(seq="") == ""

    def test_the_orientation_flag_is_required(self) -> None:
        reads = self._split_pair()
        del reads[0]["is_reverse"]
        pair = refine_breakpoints([self.VARIANT], reads)[0]
        assert pair.bp1.inserted_sequence == ""
        assert pair.bp1.strand == "+"  # no orientation either
