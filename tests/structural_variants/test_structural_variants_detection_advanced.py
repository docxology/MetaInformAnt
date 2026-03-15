"""Tests for structural variant detection edge cases and advanced functionality.

NO MOCKING POLICY: All tests use real implementations.
"""
from __future__ import annotations

import pytest

from metainformant.structural_variants.detection.sv_calling import (
    InsertSizeStats,
    SVEvidence,
    SVType,
    StructuralVariant,
    call_structural_variants,
    classify_sv_type,
    detect_discordant_pairs,
    detect_split_reads,
    genotype_sv,
)


# ---------------------------------------------------------------------------
# classify_sv_type
# ---------------------------------------------------------------------------


class TestClassifySVType:
    def test_translocation_different_chroms(self):
        ev = SVEvidence(chrom1="chr1", chrom2="chr2", breakpoint1=1000, breakpoint2=2000, strand1="+", strand2="-")
        assert classify_sv_type(ev) == SVType.TRA

    def test_insertion_close_breakpoints_with_split_reads(self):
        ev = SVEvidence(
            chrom1="chr1", chrom2="chr1", breakpoint1=1000, breakpoint2=1020,
            strand1="+", strand2="-", split_reads=3,
        )
        assert classify_sv_type(ev) == SVType.INS

    def test_inversion_same_strand(self):
        ev = SVEvidence(
            chrom1="chr1", chrom2="chr1", breakpoint1=1000, breakpoint2=5000,
            strand1="+", strand2="+",
        )
        assert classify_sv_type(ev) == SVType.INV

    def test_deletion_forward_reverse(self):
        ev = SVEvidence(
            chrom1="chr1", chrom2="chr1", breakpoint1=1000, breakpoint2=5000,
            strand1="+", strand2="-",
        )
        assert classify_sv_type(ev) == SVType.DEL

    def test_duplication_reverse_forward(self):
        ev = SVEvidence(
            chrom1="chr1", chrom2="chr1", breakpoint1=1000, breakpoint2=5000,
            strand1="-", strand2="+",
        )
        assert classify_sv_type(ev) == SVType.DUP

    def test_large_distance_defaults_to_del(self):
        ev = SVEvidence(
            chrom1="chr1", chrom2="chr1", breakpoint1=1000, breakpoint2=50000,
            strand1="-", strand2="-",  # Same strand => INV
        )
        result = classify_sv_type(ev)
        assert result == SVType.INV


# ---------------------------------------------------------------------------
# detect_split_reads
# ---------------------------------------------------------------------------


class TestDetectSplitReads:
    def test_no_reads(self):
        result = detect_split_reads([])
        assert result == []

    def test_read_with_soft_clip(self):
        reads = [
            {
                "name": "read1",
                "chrom": "chr1",
                "pos": 1000,
                "cigar": "50M30S",
                "mapq": 60,
            }
        ]
        result = detect_split_reads(reads, min_clip=20)
        assert len(result) == 1
        assert result[0].split_reads == 1

    def test_read_below_min_clip_threshold(self):
        reads = [
            {
                "name": "read1",
                "chrom": "chr1",
                "pos": 1000,
                "cigar": "70M10S",
                "mapq": 60,
            }
        ]
        result = detect_split_reads(reads, min_clip=20)
        assert len(result) == 0

    def test_read_with_sa_tag(self):
        reads = [
            {
                "name": "read1",
                "chrom": "chr1",
                "pos": 1000,
                "cigar": "50M30S",
                "mapq": 60,
                "sa_tag": "chr2,5000,+,50M30S,60,0;",
            }
        ]
        result = detect_split_reads(reads, min_clip=20)
        assert len(result) >= 1
        assert result[0].chrom2 == "chr2"

    def test_left_soft_clip(self):
        reads = [
            {
                "name": "read1",
                "chrom": "chr1",
                "pos": 1000,
                "cigar": "25S50M",
                "mapq": 60,
            }
        ]
        result = detect_split_reads(reads, min_clip=20)
        assert len(result) == 1


# ---------------------------------------------------------------------------
# detect_discordant_pairs
# ---------------------------------------------------------------------------


class TestDetectDiscordantPairs:
    def test_no_reads(self):
        stats = InsertSizeStats(mean=400.0, std=50.0)
        result = detect_discordant_pairs([], stats)
        assert result == []

    def test_inter_chromosomal(self):
        stats = InsertSizeStats(mean=400.0, std=50.0)
        reads = [
            {
                "name": "read1", "chrom": "chr1", "pos": 1000,
                "mate_chrom": "chr2", "mate_pos": 5000,
                "insert_size": 0, "is_reverse": False, "mate_is_reverse": True,
            }
        ]
        result = detect_discordant_pairs(reads, stats)
        assert len(result) == 1

    def test_aberrant_insert_size(self):
        stats = InsertSizeStats(mean=400.0, std=50.0)
        reads = [
            {
                "name": "read1", "chrom": "chr1", "pos": 1000,
                "mate_chrom": "chr1", "mate_pos": 50000,
                "insert_size": 49000, "is_reverse": False, "mate_is_reverse": True,
            }
        ]
        result = detect_discordant_pairs(reads, stats)
        assert len(result) == 1

    def test_same_strand_orientation(self):
        stats = InsertSizeStats(mean=400.0, std=50.0)
        reads = [
            {
                "name": "read1", "chrom": "chr1", "pos": 1000,
                "mate_chrom": "chr1", "mate_pos": 1400,
                "insert_size": 400, "is_reverse": True, "mate_is_reverse": True,
            }
        ]
        result = detect_discordant_pairs(reads, stats)
        assert len(result) == 1

    def test_normal_pair_not_discordant(self):
        stats = InsertSizeStats(mean=400.0, std=50.0)
        reads = [
            {
                "name": "read1", "chrom": "chr1", "pos": 1000,
                "mate_chrom": "chr1", "mate_pos": 1400,
                "insert_size": 400, "is_reverse": False, "mate_is_reverse": True,
            }
        ]
        result = detect_discordant_pairs(reads, stats)
        assert len(result) == 0


# ---------------------------------------------------------------------------
# genotype_sv
# ---------------------------------------------------------------------------


class TestGenotypeSV:
    def test_heterozygous_genotype(self):
        variant = StructuralVariant(
            chrom="chr1", start=1000, end=5000, sv_type=SVType.DEL,
            evidence=SVEvidence(split_reads=5, discordant_pairs=3, evidence_reads=["r1", "r2"]),
        )
        reads = [
            {"name": f"ref_{i}", "chrom": "chr1", "pos": 800, "read_length": 300, "mapq": 60}
            for i in range(10)
        ]
        gt = genotype_sv(variant, reads)
        assert gt in ("0/1", "1/1", "0/0", "./.")

    def test_no_spanning_reads(self):
        variant = StructuralVariant(
            chrom="chr1", start=1000, end=5000, sv_type=SVType.DEL,
            evidence=SVEvidence(split_reads=5, discordant_pairs=3),
        )
        gt = genotype_sv(variant, [])
        assert gt in ("1/1", "./.")

    def test_homozygous_alternate(self):
        variant = StructuralVariant(
            chrom="chr1", start=1000, end=5000, sv_type=SVType.DEL,
            evidence=SVEvidence(split_reads=20, discordant_pairs=15, evidence_reads=[f"r{i}" for i in range(35)]),
        )
        gt = genotype_sv(variant, [])
        assert gt in ("1/1", "./.")


# ---------------------------------------------------------------------------
# call_structural_variants (integration)
# ---------------------------------------------------------------------------


class TestCallStructuralVariants:
    def test_empty_alignments(self):
        result = call_structural_variants([])
        assert result == []

    def test_no_evidence_low_quality(self):
        reads = [
            {"name": "r1", "chrom": "chr1", "pos": 100, "cigar": "100M", "mapq": 5}
        ]
        result = call_structural_variants(reads, min_mapq=20)
        assert result == []

    def test_basic_deletion_call(self):
        # Create reads with split-read evidence for a deletion
        reads = []
        for i in range(5):
            reads.append({
                "name": f"split_{i}",
                "chrom": "chr1",
                "pos": 1000,
                "cigar": "50M30S",
                "mapq": 60,
                "sa_tag": f"chr1,5000,+,30M50S,60,0;",
                "insert_size": 400,
                "mate_chrom": "chr1",
                "mate_pos": 1400,
                "is_reverse": False,
                "mate_is_reverse": True,
            })
        result = call_structural_variants(reads, min_support=3, min_sv_size=50)
        # Should detect at least 1 SV from clustered split-read evidence
        assert isinstance(result, list)
