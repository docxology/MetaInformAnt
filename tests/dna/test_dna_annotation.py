"""Focused real-implementation tests for the dna.annotation surface.

Splice-site detection and ORF prediction are exercised on small
deterministic sequences; every score and coordinate is hand-derived in the
comments (donor/acceptor weights: max total 4.05 for donors, pyrimidine
tract weight 0.6 with a 0.1 branch-point bonus for acceptors, detection
threshold 0.3).
"""

from __future__ import annotations

import pytest

from metainformant.dna.annotation.gene_annotation import find_splice_sites
from metainformant.dna.annotation.gene_finding import predict_orfs


class TestFindSpliceSites:
    def test_strong_donor_site_detected_with_perfect_score(self):
        # "CAG|GTAAGT": GT at 3 has the full consensus
        # C|A|G|G|T|A|A|G|T -> every position contributes its maximum weight
        # (4.05/4.05 = 1.0). The second GT at 7 scores only
        # (0.15 + 0.6 + 0.1) / 4.05 = 0.21 < 0.3 (no downstream exonic
        # bases), and the "AG" at 6 gives an acceptor score of
        # (1/3 pyrimidines * 0.6) = 0.2 < 0.3 -- so exactly one site exists.
        sites = find_splice_sites("CAGGTAAGT")
        assert len(sites) == 1
        site = sites[0]
        assert site["position"] == 3
        assert site["type"] == "donor"
        assert site["dinucleotide"] == "GT"
        assert site["strand"] == "+"
        assert site["score"] == pytest.approx(1.0)

    def test_acceptor_site_detected_with_pyrimidine_tract(self):
        # 22 T's followed by AG: the tract upstream of the AG is 100%
        # pyrimidine (0.6) and the base before AG is a pyrimidine (+0.1)
        # -> score 0.7. No GT donor exists anywhere.
        sites = find_splice_sites("T" * 22 + "AG")
        acceptors = [s for s in sites if s["type"] == "acceptor"]
        donors = [s for s in sites if s["type"] == "donor"]
        assert len(acceptors) == 1 and not donors
        assert acceptors[0]["position"] == 22
        assert acceptors[0]["dinucleotide"] == "AG"
        assert acceptors[0]["score"] == pytest.approx(0.7)

    def test_weak_context_reports_nothing(self):
        # A GT whose flanks all carry minimal weights
        # (0.7/4.05 = 0.17 < 0.3) and no AG anywhere.
        assert find_splice_sites("ATCTGTTCCC") == []

    def test_multiple_sites_are_sorted_by_position(self):
        # "CAGGTAAG" + 12 T's + "CTAG" (24 nt):
        #   donor    at 3  (perfect context, score 1.0)
        #   donor    at 7  (0.444: T,A,A upstream, A,A,G downstream)
        #   acceptor at 22 (the terminal "CTAG": AG at 22-23; tract
        #                   "CAGGTAAGTTTTTTTTT" is 10/17 pyrimidine
        #                   -> 0.588*0.6 = 0.353, plus 0.1 pyrimidine
        #                   before the AG)
        seq = "CAGGTAAG" + "T" * 12 + "CTAG"
        sites = find_splice_sites(seq)
        positions = [s["position"] for s in sites]
        assert positions == sorted(positions)
        types = {(s["position"], s["type"]) for s in sites}
        assert (3, "donor") in types
        assert (7, "donor") in types
        assert (22, "acceptor") in types


class TestPredictOrfs:
    def test_forward_strand_orf(self):
        # GCG|ATG AAA TAA|GCG: ORF spans [3, 12) in frame 1; the protein
        # keeps the stop codon as '*'.
        orfs = predict_orfs("GCGATGAAATAAGCG", min_length=9)
        assert len(orfs) == 1
        orf = orfs[0]
        assert orf["frame"] == 1
        assert (orf["start"], orf["end"], orf["length"]) == (3, 12, 9)
        assert orf["sequence"] == "ATGAAATAA"
        assert orf["protein"] == "MK*"

    def test_reverse_strand_orf_maps_to_original_coordinates(self):
        # rc("ATGAAATAA") = "TTATTTCAT", embedded on the reverse strand of
        # the 15-nt sequence "GCG|TTATTTCAT|GCG". The ORF is found on the
        # reverse-complement strand at strand coordinates [3, 12), which map
        # back to the same [3, 12) interval of the original sequence.
        orfs = predict_orfs("GCG" + "TTATTTCAT" + "GCG", min_length=9)
        rev = [o for o in orfs if o["frame"] < 0]
        assert len(rev) == 1
        orf = rev[0]
        assert orf["frame"] == -1
        assert (orf["start"], orf["end"]) == (3, 12)
        assert orf["sequence"] == "ATGAAATAA"  # reported in the ORF's own orientation
        assert orf["protein"] == "MK*"

    def test_min_length_filters_short_orfs(self):
        assert predict_orfs("GCGATGAAATAAGCG", min_length=12) == []

    def test_multiple_orfs_sorted_longest_first(self):
        # "GG" + ATG AAA AAT TAA + "CC" + ATG AAG TAA + "GG":
        # a 12-nt ORF at [2, 14) in frame 3 and a 9-nt ORF at [16, 25) in
        # frame 2; no ORFs on the reverse strand.
        orfs = predict_orfs(
            "GG" + "ATGAAAAAATAA" + "CC" + "ATGAAGTAA" + "GG", min_length=9
        )
        assert [o["length"] for o in orfs] == [12, 9]
        assert (orfs[0]["frame"], orfs[0]["start"], orfs[0]["end"]) == (3, 2, 14)
        assert (orfs[1]["frame"], orfs[1]["start"], orfs[1]["end"]) == (2, 16, 25)

    def test_invalid_sequence_raises(self):
        with pytest.raises(ValueError):
            predict_orfs("ATGX123")
        with pytest.raises(ValueError):
            predict_orfs("")
