"""Regression tests for dna.sequence fixes (ORFs, repeats, consensus, restriction)."""

import pytest

from metainformant.dna.sequence import restriction
from metainformant.dna.sequence.consensus import find_consensus_breaks, quality_weighted_consensus
from metainformant.dna.sequence.core import find_orfs, find_repeats


class TestFindOrfs:
    def test_forward_frame_offsets_are_absolute(self):
        # ATG at original index 3 (frame +1), AAA, TAA stop ends at 12.
        seq = "GGGATGAAATAATTT"
        orfs = find_orfs(seq, min_length=1)
        assert (3, 12, "+1") in orfs
        # Reported intervals must slice the ORIGINAL sequence.
        for start, end, _label in orfs:
            assert 0 <= start < end <= len(seq)
            assert (end - start) % 3 == 0

    def test_frame_offset_one_start_gets_offset(self):
        # 'CCATGAAATAA': ATG at index 2 -> frame +3 in the old implementation
        # reported search-relative coordinates (0-based within the shifted
        # frame); positions must now be absolute.
        seq = "CCATGAAATAA"
        orfs = find_orfs(seq, min_length=1)
        assert any(start == 2 for start, _end, label in orfs if label == "+3")

    def test_reverse_strand_positions_map_to_original(self):
        from metainformant.dna.sequence.core import reverse_complement

        rc_design = "CCCATGAAATAAGGG"  # forward ORF (3, 12) on this strand
        seq = reverse_complement(rc_design)
        orfs = find_orfs(seq, min_length=1)
        reverse_orfs = [(s, e) for s, e, label in orfs if label.startswith("-")]
        assert reverse_orfs, "reverse-strand ORF must be found"
        for start, end in reverse_orfs:
            interval = seq[start:end]
            # The reverse complement of the reported interval reads ATG..stop.
            assert reverse_complement(interval).startswith("ATG")
            assert (end - start) % 3 == 0

    def test_stops_must_be_in_frame(self):
        # ATG TAA: the in-frame stop ends the ORF at 6.
        seq = "ATGTAATAA"  # codons: ATG TAA TAA
        orfs = find_orfs(seq, min_length=1)
        assert (0, 6, "+1") in orfs
        # The TAA starting at offset 4 is OUT of frame: it must NOT truncate
        # the ORF (no in-frame stop -> ORF extends to the sequence end).
        seq2 = "ATGATAATAA"  # codons: ATG ATA ATA | A
        orfs2 = find_orfs(seq2, min_length=1)
        assert (0, 10, "+1") in orfs2

    def test_min_length_filters(self):
        seq = "ATGTAATAAGGGCCC"
        assert not any(s == 0 for s, _e, _l in find_orfs(seq, min_length=5))


class TestFindRepeats:
    def test_overlap_only_repeat_is_reported(self):
        # 'AAA' occurs 3x overlapping in AAAAA but str.count() sees 1.
        repeats = find_repeats("AAAAA", min_length=3)
        assert repeats.get("AAA") == [0, 1, 2]

    def test_disjoint_repeat_positions(self):
        repeats = find_repeats("ACGTACGTACGT", min_length=3)
        assert repeats.get("ACG") == [0, 4, 8]
        assert repeats.get("CGT") == [1, 5, 9]

    def test_single_occurrence_excluded(self):
        repeats = find_repeats("AAACGT", min_length=3)
        assert "CGT" not in repeats


class TestQualityWeightedConsensus:
    def test_high_quality_base_wins_over_low_quality(self):
        # One read calls A with Q40, one calls T with Q10: the A must win.
        consensus = quality_weighted_consensus(["A", "T"], [[40], [10]])
        assert consensus == "A"

    def test_majority_still_wins_with_equal_quality(self):
        consensus = quality_weighted_consensus(["A", "A", "T"], [[20], [20], [20]])
        assert consensus == "A"

    def test_all_ambiguous_position_yields_n(self):
        consensus = quality_weighted_consensus(["N", "N"], [[30], [30]])
        assert consensus == "N"

    def test_dimension_mismatch_raises(self):
        with pytest.raises(ValueError):
            quality_weighted_consensus(["ATCG"], [[30, 30]])


class TestFindConsensusBreaks:
    def test_window_size_one_does_not_crash(self):
        seqs = ["ATCGATCG", "ATCGATCG", "GCTAGCTA"]
        breaks = find_consensus_breaks(seqs, window_size=1)
        assert len(breaks) == 8
        assert all(isinstance(pos, int) for pos, _score in breaks)


class TestRestrictionDigest:
    def test_digest_fragments_conserve_sequence(self):
        seq = "TTGGAATTCATCGAATTCAA"
        fragments = restriction.virtual_digest(seq, "EcoRI")
        assert "".join(fragments) == seq
        assert len(fragments) == 3  # cuts before each of the two sites

    def test_cut_at_sequence_start(self):
        seq = "GAATTCATCG"
        fragments = restriction.virtual_digest(seq, "EcoRI")
        assert fragments == [seq]  # single boundary cut -> one linear fragment

    def test_calculate_fragment_sizes_match_digest(self):
        seq = "GGGAAGCTTATCGCGAAGCTTCC"
        sizes = restriction.calculate_fragment_sizes(seq, "HindIII")
        assert sum(sizes) == len(seq)

    def test_double_digest_fragments_conserve_sequence(self):
        seq = "GAATTCGGATCCATCG"
        fragments = restriction.double_digest(seq, "EcoRI", "BamHI")
        assert "".join(fragments) == seq
        assert "GAATTC" in fragments[0]

    def test_blunt_cutter_flags_are_accurate(self):
        props = restriction.get_enzyme_properties()
        for enzyme in ("SmaI", "EcoRV", "PvuII", "AluI", "HaeIII", "RsaI"):
            assert props[enzyme]["cuts_blunt"] is True, enzyme
        for enzyme in ("EcoRI", "BamHI", "HindIII", "TaqI"):
            assert props[enzyme]["cuts_blunt"] is False, enzyme

    def test_no_placeholder_enzymes(self):
        assert "Blunt_example" not in restriction.RESTRICTION_ENZYMES
        assert "Blunt_example" not in restriction.get_enzyme_properties()

    def test_unknown_enzyme_raises(self):
        with pytest.raises(ValueError):
            restriction.virtual_digest("ACGT", "NoSuchEnzyme")

    def test_pattern_marks_cut_sites(self):
        seq = "AAGCTTAAGCTT"
        pattern = restriction.find_restriction_pattern(seq, "HindIII")
        # Boundary cut before site 0, cut before site 6, and the terminal
        # marker.
        assert pattern.count("|") == 4
