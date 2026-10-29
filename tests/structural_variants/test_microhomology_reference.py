"""Regression tests: ``detect_microhomology`` computes true microhomology from a reference, and labels its heuristic.

Microhomology is a stretch of sequence identical at both breakpoints of a junction, so it cannot be seen in one
sequence. ``detect_microhomology`` only had a mirror-match heuristic on the sequence around one junction (it
compared the base before the midpoint with the base after it, outwards): a palindrome signature, not microhomology,
returned as a bare string. It now takes the reference and the partner breakpoint and returns the longest stretch
identical at both junctions, and every result says which method produced it.

The reference sequences below are built so that the answer is known by construction:

    reference = LEFT + H + MIDDLE + H + RIGHT

with the homology ``H`` repeated at both ends of the segment ``H + MIDDLE`` that a deletion removes. The base
before the first ``H`` (the last base of LEFT) differs from the base before the second (the last base of MIDDLE),
and the base after the first ``H`` (the first base of MIDDLE) differs from the base after the second (the first of
RIGHT), so the microhomology of the deletion is exactly ``H``.
"""

from __future__ import annotations

import copy
import json
import pickle
from typing import Any

import pytest

from metainformant.structural_variants.detection.breakpoints import Breakpoint, detect_microhomology, refine_breakpoints

LEFT, H, MIDDLE, RIGHT = "GGCTTAAG", "TCATG", "CCGAT", "AGTCC"  # LEFT ends in G, MIDDLE in T; MIDDLE starts C, RIGHT A
REFERENCE = LEFT + H + MIDDLE + H + RIGHT
A = len(LEFT)  # junction before the first H: 8
B = len(LEFT) + len(H) + len(MIDDLE)  # junction before the second H: 18


def _bp(position: int, chrom: str = "chr1") -> dict[str, Any]:
    return {"chrom": chrom, "position": position}


class TestAgainstTheReference:
    def test_the_stretch_identical_at_both_junctions_is_found(self) -> None:
        # Deleting REFERENCE[8:18] leaves LEFT + H + RIGHT: the junction can sit anywhere along H.
        result = detect_microhomology(_bp(A), reference=REFERENCE, partner=_bp(B))
        assert result == H
        assert result.method == "reference"

    @pytest.mark.parametrize("shift", [0, 1, 2, 3, 4, 5])
    def test_the_junction_may_be_reported_anywhere_inside_the_homology(self, shift: int) -> None:
        # A caller can place the breakpoints anywhere along H: the junction before its first base (shift 0) through
        # the one after its last base (shift 5). The same five bases are identical at both.
        assert detect_microhomology(_bp(A + shift), reference=REFERENCE, partner=_bp(B + shift)) == H

    def test_a_junction_one_base_past_the_homology_has_none(self) -> None:
        # Shifted by 6, the bases after the junctions are C (MIDDLE[1]) and G (RIGHT[1]), those before them C and A.
        assert detect_microhomology(_bp(A + 6), reference=REFERENCE, partner=_bp(B + 6)) == ""

    def test_a_junction_one_base_before_the_homology_is_a_different_junction(self) -> None:
        # Shifted by -1 the junctions are 7 and 17, which delete a different segment. The bases after them are
        # G (LEFT[-1]) and T (MIDDLE[-1]), different; the bases before them are A and A, then A and G: one base.
        assert detect_microhomology(_bp(A - 1), reference=REFERENCE, partner=_bp(B - 1)) == "A"

    def test_the_two_breakpoints_may_be_given_in_either_order(self) -> None:
        assert detect_microhomology(_bp(B), reference=REFERENCE, partner=_bp(A)) == H

    def test_breakpoint_objects_are_accepted(self) -> None:
        result = detect_microhomology(Breakpoint("chr1", A), reference=REFERENCE, partner=Breakpoint("chr1", B))
        assert result == H and result.method == "reference"

    def test_a_region_of_the_reference_needs_its_start_coordinate(self) -> None:
        region = REFERENCE[4:]  # begins at coordinate 1004
        result = detect_microhomology(
            _bp(1_004 + A - 4), reference=region, reference_start=1_004, partner=_bp(1_004 + B - 4)
        )
        assert result == H

    def test_a_mapping_of_chromosome_sequences_serves_breakpoints_on_two_chromosomes(self) -> None:
        # chrA and chrB are two chromosomes carrying the same homology in different surroundings.
        chr_a = "ACGTTGCA" + H + "CAGGT"  # H starts at 8; before it A, after it C
        chr_b = "TTGCCATG" + H + "GATTC"  # H starts at 8; before it G, after it G
        result = detect_microhomology(
            {"chrom": "chrA", "position": 8},
            reference={"chrA": chr_a, "chrB": chr_b},
            partner={"chrom": "chrB", "position": 8},
        )
        assert result == H

    def test_per_chromosome_region_starts(self) -> None:
        chr_a, chr_b = "ACGTTGCA" + H + "CAGGT", "TTGCCATG" + H + "GATTC"
        result = detect_microhomology(
            {"chrom": "chrA", "position": 108},
            reference={"chrA": chr_a, "chrB": chr_b},
            reference_start={"chrA": 100, "chrB": 500},
            partner={"chrom": "chrB", "position": 508},
        )
        assert result == H

    def test_the_reference_may_be_soft_masked(self) -> None:
        result = detect_microhomology(_bp(A), reference=REFERENCE.lower(), partner=_bp(B))
        assert result == H  # reported in upper case, like the heuristic's

    def test_an_n_never_matches(self) -> None:
        # The same construction with the bases after the first H replaced by N on both sides: N does not extend it.
        ref = "GGCTTAAG" + "TCATG" + "NNNAT" + "TCATG" + "NNNCC"
        assert detect_microhomology(_bp(8), reference=ref, partner=_bp(18)) == H

    def test_the_stretch_is_capped_at_max_homology(self) -> None:
        # A tandem repeat of "AC": shifting the deletion of one unit anywhere along 40 bases gives the same allele.
        repeat = "T" + "AC" * 20 + "G"
        capped = detect_microhomology(_bp(1), reference=repeat, partner=_bp(3), max_homology=6)
        assert capped == "ACACAC" and len(capped) == 6
        # Uncapped, the stretch is the run shared by the two junctions: bases 1..38 against 3..40 (38 bases), the
        # next pair being A (position 39) against G (position 41).
        assert detect_microhomology(_bp(1), reference=repeat, partner=_bp(3)) == "AC" * 19

    def test_bases_beyond_the_ends_of_the_reference_are_not_examined(self) -> None:
        assert detect_microhomology(_bp(0), reference="ACACAC", partner=_bp(2)) == "ACAC"

    def test_a_homology_search_on_a_flanking_sequence_is_not_the_reference_result(self) -> None:
        # The retained allele around the junction, LEFT[-4:] + H + RIGHT[:4]: the mirror heuristic finds nothing
        # in it (C against A across the midpoint), though the junction sits in a 5 bp microhomology.
        junction = LEFT[-4:] + H + RIGHT[:4]
        heuristic = detect_microhomology(_bp(A), junction)
        assert heuristic == "" and heuristic.method == "mirror_heuristic"
        assert detect_microhomology(_bp(A), reference=REFERENCE, partner=_bp(B)) == H


class TestReferenceErrors:
    def test_a_reference_needs_the_partner_breakpoint(self) -> None:
        with pytest.raises(ValueError, match="partner"):
            detect_microhomology(_bp(A), reference=REFERENCE)

    def test_two_breakpoints_at_one_position_have_no_junction(self) -> None:
        with pytest.raises(ValueError, match="no junction"):
            detect_microhomology(_bp(A), reference=REFERENCE, partner=_bp(A))

    @pytest.mark.parametrize("position", [-1, len(REFERENCE) + 1])
    def test_a_breakpoint_outside_the_reference_is_refused(self, position: int) -> None:
        with pytest.raises(ValueError, match="outside the reference"):
            detect_microhomology(_bp(A), reference=REFERENCE, partner=_bp(position))

    def test_a_chromosome_missing_from_the_reference_is_refused(self) -> None:
        with pytest.raises(ValueError, match="chrX"):
            detect_microhomology(_bp(A), reference={"chr1": REFERENCE}, partner=_bp(B, chrom="chrX"))

    def test_a_dict_breakpoint_needs_a_position(self) -> None:
        with pytest.raises(ValueError, match="position"):
            detect_microhomology({"chrom": "chr1"}, reference=REFERENCE, partner=_bp(B))


class TestTheHeuristicSaysWhatItIs:
    def test_the_mirror_match_is_unchanged_and_labelled(self) -> None:
        result = detect_microhomology(Breakpoint(chrom="chr1", position=5), "ATTA")
        assert result == "TA"  # seq[mid-1] == seq[mid] and seq[mid-2] == seq[mid+1]
        assert result.method == "mirror_heuristic"

    def test_an_empty_result_still_says_which_method_found_nothing(self) -> None:
        for flank in ("", "A", "ACGT"):
            result = detect_microhomology({"position": 0}, flank)
            assert result == "" and not result and result.method == "mirror_heuristic"

    def test_a_run_of_n_is_not_homology(self) -> None:
        # "ANNNNA" mirrors N against N and then A against A, but N is not a base: the old result was "NNA".
        assert detect_microhomology({"position": 0}, "ANNNNA") == ""
        assert detect_microhomology({"position": 0}, "CANNNNAC") == ""
        # A run of real bases is reported in full: "CAGG|GGAC" mirrors G/G, G/G, A/A, C/C.
        assert detect_microhomology({"position": 0}, "CAGGGGAC") == "GGAC"

    def test_the_mirror_run_is_read_outwards_from_the_midpoint(self) -> None:
        # "CAG|GAC": G/G, A/A, C/C -> the bases right of the midpoint, "GAC".
        assert detect_microhomology({"position": 0}, "CAGGAC") == "GAC"

    def test_lower_case_is_read_as_upper_case(self) -> None:
        assert detect_microhomology({"position": 0}, "attA") == "TA"

    def test_max_homology_limits_the_run(self) -> None:
        # "CGAT|TAGC" mirrors T/T, A/A, G/G, C/C: "TAGC" in full, "TA" for a limit of 2.
        assert detect_microhomology({"position": 0}, "CGATTAGC") == "TAGC"
        assert detect_microhomology({"position": 0}, "CGATTAGC", max_homology=2) == "TA"

    def test_a_reference_result_is_not_confused_with_the_heuristic(self) -> None:
        heuristic = detect_microhomology({"position": 0}, "ATTA")
        truth = detect_microhomology(_bp(A), reference=REFERENCE, partner=_bp(B))
        assert (heuristic.method, truth.method) == ("mirror_heuristic", "reference")


class TestMicrohomologyValue:
    def test_it_is_a_string(self) -> None:
        from metainformant.structural_variants.detection.breakpoints import Microhomology

        value = Microhomology("TCATG", "reference")
        assert isinstance(value, str) and value == "TCATG" and len(value) == 5 and value + "A" == "TCATGA"
        assert json.dumps(value) == '"TCATG"' and {value: 1}["TCATG"] == 1

    def test_it_says_how_it_was_found_in_its_repr(self) -> None:
        from metainformant.structural_variants.detection.breakpoints import Microhomology

        assert repr(Microhomology("TA", "mirror_heuristic")) == "Microhomology('TA', method='mirror_heuristic')"

    def test_copies_and_pickles_keep_the_method(self) -> None:
        from metainformant.structural_variants.detection.breakpoints import Microhomology

        value = Microhomology("TCATG", "reference")
        for clone in (copy.copy(value), copy.deepcopy(value), pickle.loads(pickle.dumps(value))):
            assert clone == "TCATG" and clone.method == "reference" and type(clone) is Microhomology


class TestRefinedBreakpointsAreLabelled:
    # The read geometry of test_breakpoint_refinement_fixes: a 40-base read clipped after 20 aligned bases.
    LEFT20, RIGHT20 = "C" * 16 + "TACG", "GCA" + "C" + "A" * 16  # junction at 20: G|G, C|C, A|A, then T|C -> "GCA"

    def _read(self, with_sequence: bool = True) -> dict[str, Any]:
        read: dict[str, Any] = {"chrom": "chr1", "pos": 4980, "cigar": "20M20S", "name": "r1"}
        if with_sequence:
            read["seq"] = self.LEFT20 + self.RIGHT20
        return read

    def test_the_microhomology_of_a_read_is_labelled_as_the_heuristic(self) -> None:
        pair = refine_breakpoints([{"chrom": "chr1", "start": 5000, "end": 7000, "sv_type": "DEL"}], [self._read()])[0]
        assert pair.bp1.microhomology == "GCA"
        assert type(pair.bp1.microhomology) is str  # a plain string in the data class
        assert pair.bp1.microhomology_method == "mirror_heuristic"

    def test_no_sequence_means_nothing_was_searched(self) -> None:
        pair = refine_breakpoints(
            [{"chrom": "chr1", "start": 5000, "end": 7000, "sv_type": "DEL"}], [self._read(with_sequence=False)]
        )[0]
        assert pair.bp1.microhomology == "" and pair.bp1.microhomology_method == ""
        assert pair.bp2.microhomology_method == ""  # no read supports it at all

    def test_a_breakpoint_without_a_method_by_default(self) -> None:
        assert Breakpoint("chr1", 10).microhomology_method == ""
