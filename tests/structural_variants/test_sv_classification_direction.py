"""Direction-consistency tests for structural variant classification.

Pins the strand/allele direction semantics of ``classify_sv_type``: a call is
only made when the orientation evidence supports it, and weak or
unrecognisable evidence returns ``SVType.UNKNOWN`` instead of a manufactured
deletion.
"""

from __future__ import annotations

from metainformant.structural_variants.detection.sv_calling import (
    SVEvidence,
    SVType,
    classify_sv_type,
)


class TestClassifyDirectionRules:
    def test_translocation(self):
        ev = SVEvidence(
            chrom1="chr1",
            chrom2="chr2",
            breakpoint1=1000,
            breakpoint2=2000,
            strand1="+",
            strand2="-",
        )
        assert classify_sv_type(ev) == SVType.TRA

    def test_insertion_needs_split_read_support(self):
        # Split reads whose sides land ~20 bp apart at the same junction:
        # insertion.
        ev = SVEvidence(
            chrom1="chr1",
            chrom2="chr1",
            breakpoint1=1000,
            breakpoint2=1020,
            strand1="+",
            strand2="-",
            split_reads=3,
        )
        assert classify_sv_type(ev) == SVType.INS

    def test_short_span_fr_pair_without_split_support_is_unknown(self):
        # A normal-looking short-span FR pair carries no insertion evidence:
        # the old heuristic manufactured a deletion here. It must not.
        ev = SVEvidence(
            chrom1="chr1",
            chrom2="chr1",
            breakpoint1=1000,
            breakpoint2=1020,
            strand1="+",
            strand2="-",
        )
        assert classify_sv_type(ev) == SVType.UNKNOWN

    def test_fr_pair_with_span_or_split_support_is_deletion(self):
        # Head-to-head pair with a span beyond a normal fragment: deletion.
        ev = SVEvidence(
            chrom1="chr1",
            chrom2="chr1",
            breakpoint1=1000,
            breakpoint2=5000,
            strand1="+",
            strand2="-",
        )
        assert classify_sv_type(ev) == SVType.DEL
        # ... and with split-read support even at close range.
        ev = SVEvidence(
            chrom1="chr1",
            chrom2="chr1",
            breakpoint1=1000,
            breakpoint2=2100,
            strand1="+",
            strand2="-",
            split_reads=4,
        )
        assert classify_sv_type(ev) == SVType.DEL

    def test_same_strand_pairs_are_inversions(self):
        ev = SVEvidence(
            chrom1="chr1",
            chrom2="chr1",
            breakpoint1=1000,
            breakpoint2=5000,
            strand1="+",
            strand2="+",
        )
        assert classify_sv_type(ev) == SVType.INV
        ev = SVEvidence(
            chrom1="chr1",
            chrom2="chr1",
            breakpoint1=1000,
            breakpoint2=5000,
            strand1="-",
            strand2="-",
        )
        assert classify_sv_type(ev) == SVType.INV

    def test_same_strand_call_independent_of_breakpoint_order(self):
        # Reversed breakpoint order must not change the same-strand call.
        ev = SVEvidence(
            chrom1="chr1",
            chrom2="chr1",
            breakpoint1=5000,
            breakpoint2=1000,
            strand1="+",
            strand2="+",
        )
        assert classify_sv_type(ev) == SVType.INV

    def test_everted_pair_is_duplication_in_either_order(self):
        # Tail-to-tail (upstream '-', downstream '+'): duplication.
        ev = SVEvidence(
            chrom1="chr1",
            chrom2="chr1",
            breakpoint1=1000,
            breakpoint2=5000,
            strand1="-",
            strand2="+",
        )
        assert classify_sv_type(ev) == SVType.DUP
        # Breakpoints given in reversed order: normalisation maps the read at
        # the lower coordinate (strand '-') upstream, so the call is stable.
        ev = SVEvidence(
            chrom1="chr1",
            chrom2="chr1",
            breakpoint1=5000,
            breakpoint2=1000,
            strand1="+",
            strand2="-",
        )
        assert classify_sv_type(ev) == SVType.DUP

    def test_reversed_order_head_to_head_with_span_is_deletion(self):
        # The read at the lower coordinate (1000) carries '+', the read at
        # 5000 carries '-': head-to-head with a large span -> deletion.
        ev = SVEvidence(
            chrom1="chr1",
            chrom2="chr1",
            breakpoint1=5000,
            breakpoint2=1000,
            strand1="-",
            strand2="+",
        )
        assert classify_sv_type(ev) == SVType.DEL

    def test_unknown_strand_characters_never_call(self):
        # Missing orientation evidence must not become an inversion or a
        # deletion.
        ev = SVEvidence(
            chrom1="chr1",
            chrom2="chr1",
            breakpoint1=1000,
            breakpoint2=2100,
            strand1=".",
            strand2=".",
        )
        assert classify_sv_type(ev) == SVType.UNKNOWN
        ev = SVEvidence(
            chrom1="chr1",
            chrom2="chr1",
            breakpoint1=1000,
            breakpoint2=2100,
            strand1=".",
            strand2="+",
        )
        assert classify_sv_type(ev) == SVType.UNKNOWN

    def test_identical_breakpoints_without_evidence_are_unknown(self):
        ev = SVEvidence(
            chrom1="chr1",
            chrom2="chr1",
            breakpoint1=1000,
            breakpoint2=1000,
            strand1="+",
            strand2="-",
        )
        assert classify_sv_type(ev) == SVType.UNKNOWN
