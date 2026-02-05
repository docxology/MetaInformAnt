"""Comprehensive tests for the structural_variants module.

Tests cover all major functionality: CNV detection, SV calling, breakpoint
refinement, annotation, filtering, merging, and visualization.
Uses real algorithmic implementations (no mocking).
"""

from __future__ import annotations

import math
import os
import tempfile
from pathlib import Path
from typing import Any

import pytest

# Optional dependency checks
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

try:
    import matplotlib

    matplotlib.use("Agg")
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def sample_depth_data() -> dict[str, list[float]]:
    """Simulated read depth data with a deletion and duplication."""
    if not HAS_NUMPY:
        pytest.skip("numpy required")
    rng = np.random.default_rng(42)
    # 500 bins: normal(100,10) with a deletion (bins 100-150) and dup (bins 300-350)
    depths = rng.normal(100, 10, 500).tolist()
    for i in range(100, 150):
        depths[i] = max(0, rng.normal(40, 5))  # ~60% drop => deletion
    for i in range(300, 350):
        depths[i] = rng.normal(160, 10)  # ~60% increase => duplication
    return {"chr1": depths}


@pytest.fixture
def sample_alignments() -> list[dict[str, Any]]:
    """Simulated alignment data with split-read and discordant evidence."""
    alignments: list[dict[str, Any]] = []

    # Normal reads
    for i in range(50):
        alignments.append(
            {
                "name": f"normal_{i}",
                "chrom": "chr1",
                "pos": 1000 + i * 100,
                "mapq": 60,
                "cigar": "150M",
                "mate_chrom": "chr1",
                "mate_pos": 1000 + i * 100 + 400,
                "insert_size": 400,
                "is_reverse": False,
                "mate_is_reverse": True,
                "read_length": 150,
            }
        )

    # Split reads at position ~5000 (deletion evidence)
    for i in range(5):
        alignments.append(
            {
                "name": f"split_{i}",
                "chrom": "chr1",
                "pos": 4900 + i * 10,
                "mapq": 50,
                "cigar": "80M70S",
                "sa_tag": f"chr1,8000,+,70M80S,50,0;",
                "mate_chrom": "chr1",
                "mate_pos": 4900 + i * 10 + 400,
                "insert_size": 400,
                "is_reverse": False,
                "mate_is_reverse": True,
                "read_length": 150,
            }
        )

    # Discordant pairs (large insert size => deletion evidence)
    for i in range(5):
        alignments.append(
            {
                "name": f"discordant_{i}",
                "chrom": "chr1",
                "pos": 4800 + i * 20,
                "mapq": 55,
                "cigar": "150M",
                "mate_chrom": "chr1",
                "mate_pos": 8200 + i * 20,
                "insert_size": 3400,
                "is_reverse": False,
                "mate_is_reverse": True,
                "read_length": 150,
            }
        )

    return alignments


@pytest.fixture
def sample_gene_db() -> list[dict[str, Any]]:
    """Sample gene database."""
    return [
        {"chrom": "chr1", "start": 1000, "end": 5000, "name": "GENE_A", "is_coding": True},
        {"chrom": "chr1", "start": 6000, "end": 9000, "name": "GENE_B", "is_coding": True},
        {"chrom": "chr1", "start": 12000, "end": 15000, "name": "GENE_C", "is_coding": False},
        {"chrom": "chr2", "start": 1000, "end": 3000, "name": "GENE_D", "is_coding": True},
    ]


@pytest.fixture
def sample_variants() -> list[dict[str, Any]]:
    """Sample structural variant calls."""
    return [
        {"chrom": "chr1", "start": 4000, "end": 8000, "sv_type": "DEL", "size": 4000, "quality": 50.0, "support": 8},
        {"chrom": "chr1", "start": 11000, "end": 16000, "sv_type": "DUP", "size": 5000, "quality": 35.0, "support": 5},
        {"chrom": "chr1", "start": 20000, "end": 20100, "sv_type": "INS", "size": 100, "quality": 25.0, "support": 4},
        {"chrom": "chr2", "start": 500, "end": 3500, "sv_type": "DEL", "size": 3000, "quality": 60.0, "support": 12},
    ]


@pytest.fixture
def sample_blacklist() -> list[dict[str, Any]]:
    """Sample blacklist regions."""
    return [
        {"chrom": "chr1", "start": 19000, "end": 21000},
        {"chrom": "chr3", "start": 0, "end": 50000},
    ]


# ---------------------------------------------------------------------------
# Tests: detection/cnv.py
# ---------------------------------------------------------------------------


class TestCNVDetection:
    """Tests for CNV detection from read depth data."""

    @pytest.mark.skipif(not HAS_NUMPY, reason="numpy required")
    def test_detect_cnv_from_depth_basic(self, sample_depth_data: dict[str, list[float]]) -> None:
        from metainformant.structural_variants.detection.cnv import detect_cnv_from_depth

        results = detect_cnv_from_depth(sample_depth_data, window_size=1000)

        assert "chr1" in results
        result = results["chr1"]
        assert len(result.segments) > 0
        assert result.method == "segmentation"
        assert result.bin_size == 1000

        # Check that we detect at least one non-neutral segment
        non_neutral = [s for s in result.segments if s.state != "NEUTRAL"]
        assert len(non_neutral) > 0

    @pytest.mark.skipif(not HAS_NUMPY, reason="numpy required")
    def test_detect_cnv_empty_data(self) -> None:
        from metainformant.structural_variants.detection.cnv import detect_cnv_from_depth

        with pytest.raises(ValueError, match="must not be empty"):
            detect_cnv_from_depth({})

    @pytest.mark.skipif(not HAS_NUMPY, reason="numpy required")
    def test_segment_coverage_constant(self) -> None:
        from metainformant.structural_variants.detection.cnv import segment_coverage

        # Constant array should produce one segment
        data = [0.0] * 100
        segments = segment_coverage(data)
        assert len(segments) == 1
        assert segments[0][0] == 0
        assert segments[0][1] == 100

    @pytest.mark.skipif(not HAS_NUMPY, reason="numpy required")
    def test_segment_coverage_with_changepoint(self) -> None:
        from metainformant.structural_variants.detection.cnv import segment_coverage

        # Clear change-point
        data = [0.0] * 50 + [2.0] * 50
        segments = segment_coverage(data, significance=0.05)
        assert len(segments) >= 2

    @pytest.mark.skipif(not HAS_NUMPY, reason="numpy required")
    def test_call_cnv_states(self) -> None:
        from metainformant.structural_variants.detection.cnv import call_cnv_states

        segments = [(0, 50, -0.8), (50, 100, 0.0), (100, 150, 0.6)]
        result = call_cnv_states(segments)

        assert len(result) == 3
        assert result[0].state == "DEL"
        assert result[1].state == "NEUTRAL"
        assert result[2].state == "DUP"

    @pytest.mark.skipif(not HAS_NUMPY, reason="numpy required")
    def test_call_cnv_states_extreme(self) -> None:
        from metainformant.structural_variants.detection.cnv import call_cnv_states

        segments = [(0, 10, -2.0), (10, 20, 1.5)]
        result = call_cnv_states(segments)

        assert result[0].state == "HOMODEL"
        assert result[1].state == "AMP"

    @pytest.mark.skipif(not HAS_NUMPY, reason="numpy required")
    def test_merge_adjacent_segments(self) -> None:
        from metainformant.structural_variants.detection.cnv import CNVSegment, merge_adjacent_segments

        segs = [
            CNVSegment(chrom="chr1", start=0, end=1000, mean_log2ratio=-0.5, n_bins=10, state="DEL", cn=1),
            CNVSegment(chrom="chr1", start=1000, end=2000, mean_log2ratio=-0.6, n_bins=10, state="DEL", cn=1),
            CNVSegment(chrom="chr1", start=5000, end=6000, mean_log2ratio=0.0, n_bins=10, state="NEUTRAL", cn=2),
        ]

        merged = merge_adjacent_segments(segs, max_gap=1000)
        assert len(merged) == 2
        assert merged[0].start == 0
        assert merged[0].end == 2000
        assert merged[0].state == "DEL"

    @pytest.mark.skipif(not HAS_NUMPY, reason="numpy required")
    def test_calculate_log2_ratio(self) -> None:
        from metainformant.structural_variants.detection.cnv import calculate_log2_ratio

        tumor = [100, 200, 50, 100]
        normal = [100, 100, 100, 100]

        ratios = calculate_log2_ratio(tumor, normal)

        assert len(ratios) == 4
        # Tumor=200, Normal=100 should have positive log2 ratio
        # (after median centering, the relative values matter)
        assert isinstance(ratios, np.ndarray)

    @pytest.mark.skipif(not HAS_NUMPY, reason="numpy required")
    def test_calculate_log2_ratio_gc_correction(self) -> None:
        from metainformant.structural_variants.detection.cnv import calculate_log2_ratio

        tumor = [100, 100, 100, 100, 100]
        normal = [100, 100, 100, 100, 100]
        gc = [0.3, 0.4, 0.5, 0.6, 0.7]

        ratios = calculate_log2_ratio(tumor, normal, gc_content=gc)
        assert len(ratios) == 5

    def test_calculate_log2_ratio_mismatch(self) -> None:
        from metainformant.structural_variants.detection.cnv import calculate_log2_ratio

        with pytest.raises(ValueError, match="mismatch"):
            calculate_log2_ratio([1, 2], [1, 2, 3])


# ---------------------------------------------------------------------------
# Tests: detection/sv_calling.py
# ---------------------------------------------------------------------------


class TestSVCalling:
    """Tests for structural variant calling."""

    def test_call_structural_variants(self, sample_alignments: list[dict[str, Any]]) -> None:
        from metainformant.structural_variants.detection.sv_calling import call_structural_variants

        variants = call_structural_variants(sample_alignments, min_support=2)

        assert isinstance(variants, list)
        # Should find at least one variant from our simulated data
        assert len(variants) >= 1

        for v in variants:
            assert v.chrom == "chr1"
            assert v.quality > 0

    def test_detect_split_reads(self, sample_alignments: list[dict[str, Any]]) -> None:
        from metainformant.structural_variants.detection.sv_calling import detect_split_reads

        evidence = detect_split_reads(sample_alignments, min_clip=20)

        # We have 5 split reads with 70S clips
        assert len(evidence) >= 5

        for ev in evidence:
            assert ev.split_reads == 1
            assert len(ev.evidence_reads) == 1

    def test_detect_discordant_pairs(self, sample_alignments: list[dict[str, Any]]) -> None:
        from metainformant.structural_variants.detection.sv_calling import (
            InsertSizeStats,
            detect_discordant_pairs,
        )

        stats = InsertSizeStats(mean=400.0, std=50.0)
        evidence = detect_discordant_pairs(sample_alignments, stats)

        # We have 5 discordant pairs with large insert sizes
        assert len(evidence) >= 5

        for ev in evidence:
            assert ev.discordant_pairs == 1

    def test_classify_sv_type_deletion(self) -> None:
        from metainformant.structural_variants.detection.sv_calling import SVEvidence, SVType, classify_sv_type

        evidence = SVEvidence(
            split_reads=3,
            breakpoint1=1000,
            breakpoint2=5000,
            chrom1="chr1",
            chrom2="chr1",
            strand1="+",
            strand2="-",
        )
        assert classify_sv_type(evidence) == SVType.DEL

    def test_classify_sv_type_duplication(self) -> None:
        from metainformant.structural_variants.detection.sv_calling import SVEvidence, SVType, classify_sv_type

        evidence = SVEvidence(
            split_reads=3,
            breakpoint1=1000,
            breakpoint2=5000,
            chrom1="chr1",
            chrom2="chr1",
            strand1="-",
            strand2="+",
        )
        assert classify_sv_type(evidence) == SVType.DUP

    def test_classify_sv_type_inversion(self) -> None:
        from metainformant.structural_variants.detection.sv_calling import SVEvidence, SVType, classify_sv_type

        evidence = SVEvidence(
            split_reads=3,
            breakpoint1=1000,
            breakpoint2=5000,
            chrom1="chr1",
            chrom2="chr1",
            strand1="+",
            strand2="+",
        )
        assert classify_sv_type(evidence) == SVType.INV

    def test_classify_sv_type_translocation(self) -> None:
        from metainformant.structural_variants.detection.sv_calling import SVEvidence, SVType, classify_sv_type

        evidence = SVEvidence(
            split_reads=3,
            breakpoint1=1000,
            breakpoint2=5000,
            chrom1="chr1",
            chrom2="chr2",
            strand1="+",
            strand2="-",
        )
        assert classify_sv_type(evidence) == SVType.TRA

    def test_classify_sv_type_insertion(self) -> None:
        from metainformant.structural_variants.detection.sv_calling import SVEvidence, SVType, classify_sv_type

        evidence = SVEvidence(
            split_reads=3,
            breakpoint1=1000,
            breakpoint2=1020,
            chrom1="chr1",
            chrom2="chr1",
            strand1="+",
            strand2="-",
        )
        assert classify_sv_type(evidence) == SVType.INS

    def test_genotype_sv(self) -> None:
        from metainformant.structural_variants.detection.sv_calling import (
            StructuralVariant,
            SVEvidence,
            SVType,
            genotype_sv,
        )

        variant = StructuralVariant(
            chrom="chr1",
            start=5000,
            end=8000,
            sv_type=SVType.DEL,
            evidence=SVEvidence(split_reads=5, discordant_pairs=3, evidence_reads=["r1", "r2"]),
        )

        # Create reads: some spanning the breakpoint (reference), some supporting
        reads = [
            {"name": "r1", "chrom": "chr1", "pos": 4800, "read_length": 500, "mapq": 60},
            {"name": "ref1", "chrom": "chr1", "pos": 4800, "read_length": 500, "mapq": 60},
            {"name": "ref2", "chrom": "chr1", "pos": 4800, "read_length": 500, "mapq": 60},
        ]

        gt = genotype_sv(variant, reads)
        assert gt in ("0/0", "0/1", "1/1")


# ---------------------------------------------------------------------------
# Tests: detection/breakpoints.py
# ---------------------------------------------------------------------------


class TestBreakpoints:
    """Tests for breakpoint detection and refinement."""

    def test_refine_breakpoints(self) -> None:
        from metainformant.structural_variants.detection.breakpoints import refine_breakpoints

        variants = [{"chrom": "chr1", "start": 5000, "end": 8000, "sv_type": "DEL"}]
        reads = [
            {"chrom": "chr1", "pos": 4950, "cigar": "50M100S", "name": "r1"},
            {"chrom": "chr1", "pos": 4950, "cigar": "50M100S", "name": "r2"},
            {"chrom": "chr1", "pos": 7950, "cigar": "100S50M", "name": "r3"},
        ]

        pairs = refine_breakpoints(variants, reads)
        assert len(pairs) == 1
        assert pairs[0].bp1.chrom == "chr1"
        assert pairs[0].sv_type == "DEL"

    def test_detect_microhomology(self) -> None:
        from metainformant.structural_variants.detection.breakpoints import detect_microhomology

        # Sequence with microhomology at the breakpoint (center)
        seq = "AATTCCGGAATTCCGG"  # Palindromic-ish
        bp = {"position": 100}

        result = detect_microhomology(bp, seq)
        assert isinstance(result, str)

    def test_detect_microhomology_empty(self) -> None:
        from metainformant.structural_variants.detection.breakpoints import detect_microhomology

        result = detect_microhomology({"position": 0}, "")
        assert result == ""

    def test_cluster_breakpoints(self) -> None:
        from metainformant.structural_variants.detection.breakpoints import Breakpoint, cluster_breakpoints

        bps = [
            Breakpoint(chrom="chr1", position=1000, support=3),
            Breakpoint(chrom="chr1", position=1050, support=5),
            Breakpoint(chrom="chr1", position=1080, support=2),
            Breakpoint(chrom="chr1", position=5000, support=4),
            Breakpoint(chrom="chr2", position=1000, support=1),
        ]

        clusters = cluster_breakpoints(bps, max_distance=100)

        # Should get 3 clusters: two close chr1 bps, one distant chr1, one chr2
        assert len(clusters) == 3
        # Largest cluster first
        assert len(clusters[0]) == 3

    def test_cluster_breakpoints_dict_input(self) -> None:
        from metainformant.structural_variants.detection.breakpoints import cluster_breakpoints

        bps = [
            {"chrom": "chr1", "position": 1000},
            {"chrom": "chr1", "position": 1050},
        ]

        clusters = cluster_breakpoints(bps, max_distance=100)
        assert len(clusters) == 1

    def test_calculate_breakpoint_confidence(self) -> None:
        from metainformant.structural_variants.detection.breakpoints import calculate_breakpoint_confidence

        # High confidence
        conf_high = calculate_breakpoint_confidence(
            {"support": 15, "total_clips": 20, "position_std": 2.0, "mapq_mean": 55.0}
        )

        # Low confidence
        conf_low = calculate_breakpoint_confidence(
            {"support": 1, "total_clips": 10, "position_std": 50.0, "mapq_mean": 10.0}
        )

        assert 0.0 <= conf_high <= 1.0
        assert 0.0 <= conf_low <= 1.0
        assert conf_high > conf_low

    def test_calculate_breakpoint_confidence_zero(self) -> None:
        from metainformant.structural_variants.detection.breakpoints import calculate_breakpoint_confidence

        conf = calculate_breakpoint_confidence({"support": 0})
        assert conf == 0.0


# ---------------------------------------------------------------------------
# Tests: annotation/overlap.py
# ---------------------------------------------------------------------------


class TestOverlapAnnotation:
    """Tests for gene and regulatory overlap annotation."""

    def test_annotate_gene_overlap(
        self, sample_variants: list[dict[str, Any]], sample_gene_db: list[dict[str, Any]]
    ) -> None:
        from metainformant.structural_variants.annotation.overlap import annotate_gene_overlap

        annotated = annotate_gene_overlap(sample_variants, sample_gene_db)

        assert len(annotated) == 4

        # First variant (chr1:4000-8000) should overlap GENE_A and GENE_B
        first = annotated[0]
        assert "GENE_A" in first["overlapping_genes"]
        assert "GENE_B" in first["overlapping_genes"]
        assert first["n_genes_affected"] == 2

        # chr2 variant should overlap GENE_D
        chr2_var = [v for v in annotated if v["chrom"] == "chr2"][0]
        assert "GENE_D" in chr2_var["overlapping_genes"]

    def test_annotate_regulatory_overlap(self) -> None:
        from metainformant.structural_variants.annotation.overlap import annotate_regulatory_overlap

        variants = [{"chrom": "chr1", "start": 1000, "end": 5000}]
        reg_db = [
            {"chrom": "chr1", "start": 2000, "end": 3000, "name": "enhancer_1", "feature_type": "enhancer"},
            {"chrom": "chr1", "start": 4500, "end": 5500, "name": "promoter_1", "feature_type": "promoter"},
        ]

        annotated = annotate_regulatory_overlap(variants, reg_db)
        assert len(annotated) == 1
        assert "enhancer_1" in annotated[0]["overlapping_regulatory"]
        assert "promoter_1" in annotated[0]["overlapping_regulatory"]
        assert "enhancer" in annotated[0]["regulatory_types"]

    def test_calculate_overlap_fraction(self) -> None:
        from metainformant.structural_variants.annotation.overlap import calculate_overlap_fraction

        # 50% overlap
        frac_a, frac_b = calculate_overlap_fraction((0, 100), (50, 150))
        assert abs(frac_a - 0.5) < 0.01
        assert abs(frac_b - 0.5) < 0.01

        # Complete containment
        frac_a, frac_b = calculate_overlap_fraction((0, 100), (20, 80))
        assert abs(frac_a - 0.6) < 0.01
        assert abs(frac_b - 1.0) < 0.01

        # No overlap
        frac_a, frac_b = calculate_overlap_fraction((0, 100), (200, 300))
        assert frac_a == 0.0
        assert frac_b == 0.0

    def test_find_nearest_gene_overlapping(self, sample_gene_db: list[dict[str, Any]]) -> None:
        from metainformant.structural_variants.annotation.overlap import find_nearest_gene

        variant = {"chrom": "chr1", "start": 2000, "end": 3000}
        result = find_nearest_gene(variant, sample_gene_db)
        assert result["nearest_gene"] == "GENE_A"
        assert result["distance"] == 0
        assert result["direction"] == "overlapping"

    def test_find_nearest_gene_nearby(self, sample_gene_db: list[dict[str, Any]]) -> None:
        from metainformant.structural_variants.annotation.overlap import find_nearest_gene

        variant = {"chrom": "chr1", "start": 5500, "end": 5600}
        result = find_nearest_gene(variant, sample_gene_db)
        assert result["nearest_gene"] in ("GENE_A", "GENE_B")
        assert result["distance"] >= 0

    def test_find_nearest_gene_none_found(self) -> None:
        from metainformant.structural_variants.annotation.overlap import find_nearest_gene

        variant = {"chrom": "chr99", "start": 1000, "end": 2000}
        result = find_nearest_gene(variant, [])
        assert result["nearest_gene"] == ""
        assert result["distance"] == -1

    def test_interval_index_query(self) -> None:
        from metainformant.structural_variants.annotation.overlap import GenomicInterval, IntervalIndex

        intervals = [
            GenomicInterval(chrom="chr1", start=100, end=200, name="A"),
            GenomicInterval(chrom="chr1", start=300, end=400, name="B"),
            GenomicInterval(chrom="chr1", start=350, end=500, name="C"),
        ]
        index = IntervalIndex(intervals)

        # Query overlapping 150-350 should find A and B (and possibly C)
        result = index.query_overlap("chr1", 150, 360)
        names = [iv.name for iv in result]
        assert "A" in names
        assert "B" in names


# ---------------------------------------------------------------------------
# Tests: annotation/functional_impact.py
# ---------------------------------------------------------------------------


class TestFunctionalImpact:
    """Tests for functional impact prediction."""

    def test_predict_functional_impact_deletion(self, sample_gene_db: list[dict[str, Any]]) -> None:
        from metainformant.structural_variants.annotation.functional_impact import predict_functional_impact

        variant = {"chrom": "chr1", "start": 4000, "end": 8000, "sv_type": "DEL"}
        result = predict_functional_impact(variant, sample_gene_db)

        assert result.impact_level in ("HIGH", "MODERATE")
        assert len(result.affected_genes) > 0
        assert result.pathogenicity_score > 0

    def test_predict_functional_impact_intergenic(self, sample_gene_db: list[dict[str, Any]]) -> None:
        from metainformant.structural_variants.annotation.functional_impact import predict_functional_impact

        variant = {"chrom": "chr1", "start": 10000, "end": 11000, "sv_type": "DEL"}
        result = predict_functional_impact(variant, sample_gene_db)

        assert len(result.affected_genes) == 0
        assert result.impact_level in ("LOW", "MODIFIER")

    def test_assess_dosage_sensitivity(self) -> None:
        from metainformant.structural_variants.annotation.functional_impact import assess_dosage_sensitivity

        hi_db = {"BRCA1": 0.9, "NORMAL_GENE": 0.1}
        pli_db = {"BRCA1": 0.99, "NORMAL_GENE": 0.01}

        brca_result = assess_dosage_sensitivity("BRCA1", hi_db, pli_db=pli_db)
        assert brca_result.is_haploinsufficient
        assert brca_result.haploinsufficiency_score == 0.9

        normal_result = assess_dosage_sensitivity("NORMAL_GENE", hi_db, pli_db=pli_db)
        assert not normal_result.is_haploinsufficient

    def test_predict_tad_disruption(self) -> None:
        from metainformant.structural_variants.annotation.functional_impact import predict_tad_disruption

        variant = {"chrom": "chr1", "start": 49000, "end": 60000, "sv_type": "DEL"}
        boundaries = [
            {"chrom": "chr1", "start": 49500, "end": 50500, "genes": ["GENE_X"]},
            {"chrom": "chr1", "start": 100000, "end": 101000, "genes": ["GENE_Y"]},
        ]

        result = predict_tad_disruption(variant, boundaries)
        assert result.n_boundaries_disrupted == 1
        assert "GENE_X" in result.genes_in_affected_tads
        assert result.disruption_score > 0

    def test_predict_tad_no_disruption(self) -> None:
        from metainformant.structural_variants.annotation.functional_impact import predict_tad_disruption

        variant = {"chrom": "chr1", "start": 1000, "end": 2000, "sv_type": "DEL"}
        boundaries = [
            {"chrom": "chr1", "start": 50000, "end": 51000},
        ]

        result = predict_tad_disruption(variant, boundaries)
        assert result.n_boundaries_disrupted == 0
        assert result.disruption_score == 0.0

    def test_score_pathogenicity(self) -> None:
        from metainformant.structural_variants.annotation.functional_impact import score_pathogenicity

        # Pathogenic variant
        score_high = score_pathogenicity(
            {"chrom": "chr1", "start": 0, "end": 500_000, "sv_type": "DEL"},
            {
                "overlapping_genes": ["GENE_A", "GENE_B"],
                "dosage_sensitive_genes": ["GENE_A"],
                "tad_disrupted": True,
                "impact_level": "HIGH",
                "n_coding_genes": 2,
            },
        )

        # Benign variant
        score_low = score_pathogenicity(
            {"chrom": "chr1", "start": 0, "end": 100, "sv_type": "INS"},
            {
                "overlapping_genes": [],
                "dosage_sensitive_genes": [],
                "tad_disrupted": False,
                "impact_level": "MODIFIER",
                "n_coding_genes": 0,
            },
        )

        assert 0.0 <= score_high <= 1.0
        assert 0.0 <= score_low <= 1.0
        assert score_high > score_low


# ---------------------------------------------------------------------------
# Tests: filtering/quality_filter.py
# ---------------------------------------------------------------------------


class TestQualityFiltering:
    """Tests for quality-based filtering."""

    def test_filter_by_quality(self, sample_variants: list[dict[str, Any]]) -> None:
        from metainformant.structural_variants.filtering.quality_filter import filter_by_quality

        passed, stats = filter_by_quality(sample_variants, min_qual=30.0, min_support=3)

        assert len(passed) <= len(sample_variants)
        assert stats.input_count == 4
        assert stats.output_count == len(passed)
        assert stats.filter_name == "quality_filter"

    def test_filter_by_quality_strict(self, sample_variants: list[dict[str, Any]]) -> None:
        from metainformant.structural_variants.filtering.quality_filter import filter_by_quality

        passed, stats = filter_by_quality(sample_variants, min_qual=100.0)
        assert len(passed) == 0

    def test_filter_by_size(self, sample_variants: list[dict[str, Any]]) -> None:
        from metainformant.structural_variants.filtering.quality_filter import filter_by_size

        passed, stats = filter_by_size(sample_variants, min_size=1000)
        assert all(v.get("size", 0) >= 1000 or v.get("sv_type") in ("TRA", "BND") for v in passed)

    def test_filter_by_size_max(self, sample_variants: list[dict[str, Any]]) -> None:
        from metainformant.structural_variants.filtering.quality_filter import filter_by_size

        passed, stats = filter_by_size(sample_variants, min_size=50, max_size=3500)
        sizes = [v.get("size", 0) for v in passed]
        assert all(s <= 3500 for s in sizes)

    def test_filter_by_frequency(self) -> None:
        from metainformant.structural_variants.filtering.quality_filter import filter_by_frequency

        variants = [
            {"chrom": "chr1", "start": 1000, "end": 5000, "sv_type": "DEL", "af": 0.001},
            {"chrom": "chr1", "start": 10000, "end": 15000, "sv_type": "DEL", "af": 0.05},
            {"chrom": "chr1", "start": 20000, "end": 25000, "sv_type": "DUP"},
        ]

        passed, stats = filter_by_frequency(variants, max_af=0.01)
        assert len(passed) == 2  # First (AF=0.001) and third (no AF)

    def test_filter_by_frequency_no_db(self, sample_variants: list[dict[str, Any]]) -> None:
        from metainformant.structural_variants.filtering.quality_filter import filter_by_frequency

        passed, stats = filter_by_frequency(sample_variants, population_db=None)
        assert len(passed) == len(sample_variants)

    def test_apply_blacklist(
        self, sample_variants: list[dict[str, Any]], sample_blacklist: list[dict[str, Any]]
    ) -> None:
        from metainformant.structural_variants.filtering.quality_filter import apply_blacklist

        passed, stats = apply_blacklist(sample_variants, sample_blacklist)

        # The INS at chr1:20000-20100 overlaps blacklist chr1:19000-21000
        assert stats.filtered_count >= 1
        assert all(v.get("sv_type") != "INS" for v in passed if v["chrom"] == "chr1")

    def test_apply_blacklist_empty(self, sample_variants: list[dict[str, Any]]) -> None:
        from metainformant.structural_variants.filtering.quality_filter import apply_blacklist

        passed, stats = apply_blacklist(sample_variants, [])
        assert len(passed) == len(sample_variants)


# ---------------------------------------------------------------------------
# Tests: filtering/merge.py
# ---------------------------------------------------------------------------


class TestMerge:
    """Tests for multi-caller merging and deduplication."""

    def test_merge_callsets(self) -> None:
        from metainformant.structural_variants.filtering.merge import merge_callsets

        callsets = {
            "delly": [
                {"chrom": "chr1", "start": 1000, "end": 5000, "sv_type": "DEL", "quality": 50.0},
                {"chrom": "chr1", "start": 10000, "end": 12000, "sv_type": "DUP", "quality": 30.0},
            ],
            "manta": [
                {"chrom": "chr1", "start": 1050, "end": 5050, "sv_type": "DEL", "quality": 60.0},
                {"chrom": "chr2", "start": 5000, "end": 8000, "sv_type": "INV", "quality": 40.0},
            ],
        }

        merged, stats = merge_callsets(callsets, min_overlap=0.5)

        assert stats.n_input_callsets == 2
        assert stats.n_input_variants == 4
        assert stats.n_output_variants > 0

        # The two DELs should merge
        del_variants = [m for m in merged if m.sv_type == "DEL"]
        assert len(del_variants) >= 1
        if del_variants:
            assert del_variants[0].n_callers == 2

    def test_merge_callsets_no_type_match(self) -> None:
        from metainformant.structural_variants.filtering.merge import merge_callsets

        callsets = {
            "caller1": [{"chrom": "chr1", "start": 1000, "end": 5000, "sv_type": "DEL"}],
            "caller2": [{"chrom": "chr1", "start": 1000, "end": 5000, "sv_type": "DUP"}],
        }

        merged, _ = merge_callsets(callsets, type_match=False)
        assert len(merged) == 1  # Should merge when type_match=False

    def test_calculate_reciprocal_overlap(self) -> None:
        from metainformant.structural_variants.filtering.merge import calculate_reciprocal_overlap

        sv1 = {"chrom": "chr1", "start": 1000, "end": 5000}
        sv2 = {"chrom": "chr1", "start": 1000, "end": 5000}
        assert calculate_reciprocal_overlap(sv1, sv2) == 1.0

        sv3 = {"chrom": "chr1", "start": 3000, "end": 7000}
        ro = calculate_reciprocal_overlap(sv1, sv3)
        assert 0.0 < ro < 1.0

        sv4 = {"chrom": "chr2", "start": 1000, "end": 5000}
        assert calculate_reciprocal_overlap(sv1, sv4) == 0.0

    def test_survivor_merge(self) -> None:
        from metainformant.structural_variants.filtering.merge import survivor_merge

        callsets = [
            [
                {"chrom": "chr1", "start": 1000, "end": 5000, "sv_type": "DEL", "_caller": "a"},
                {"chrom": "chr1", "start": 10000, "end": 12000, "sv_type": "DUP", "_caller": "a"},
            ],
            [
                {"chrom": "chr1", "start": 1100, "end": 4900, "sv_type": "DEL", "_caller": "b"},
            ],
        ]

        merged, stats = survivor_merge(callsets, max_distance=500, min_callers=2)

        # Only the DEL should survive (detected by both callers)
        assert len(merged) >= 1
        del_merged = [m for m in merged if m.sv_type == "DEL"]
        assert len(del_merged) == 1
        assert del_merged[0].n_callers == 2

    def test_deduplicate_variants(self) -> None:
        from metainformant.structural_variants.filtering.merge import deduplicate_variants

        variants = [
            {"chrom": "chr1", "start": 1000, "end": 5000, "sv_type": "DEL", "quality": 50.0},
            {"chrom": "chr1", "start": 1010, "end": 5010, "sv_type": "DEL", "quality": 60.0},
            {"chrom": "chr1", "start": 1020, "end": 4990, "sv_type": "DEL", "quality": 40.0},
            {"chrom": "chr2", "start": 1000, "end": 5000, "sv_type": "DUP", "quality": 30.0},
        ]

        deduped, n_removed = deduplicate_variants(variants)

        # Three chr1 DELs should merge into one, chr2 DUP stays
        assert len(deduped) == 2
        assert n_removed == 2

        # Should keep the highest quality
        chr1_var = [v for v in deduped if v["chrom"] == "chr1"][0]
        assert chr1_var["quality"] == 60.0

    def test_deduplicate_empty(self) -> None:
        from metainformant.structural_variants.filtering.merge import deduplicate_variants

        deduped, n_removed = deduplicate_variants([])
        assert len(deduped) == 0
        assert n_removed == 0


# ---------------------------------------------------------------------------
# Tests: visualization/plots.py
# ---------------------------------------------------------------------------


class TestVisualization:
    """Tests for SV visualization plots."""

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib required")
    def test_plot_sv_type_summary(self, sample_variants: list[dict[str, Any]], tmp_path: Path) -> None:
        from metainformant.structural_variants.visualization.plots import plot_sv_type_summary

        output = tmp_path / "sv_types.png"
        fig = plot_sv_type_summary(sample_variants, output_path=str(output))

        assert output.exists()
        assert output.stat().st_size > 0

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib required")
    @pytest.mark.skipif(not HAS_NUMPY, reason="numpy required")
    def test_plot_sv_size_distribution(self, sample_variants: list[dict[str, Any]], tmp_path: Path) -> None:
        from metainformant.structural_variants.visualization.plots import plot_sv_size_distribution

        output = tmp_path / "sv_sizes.png"
        fig = plot_sv_size_distribution(sample_variants, output_path=str(output))

        assert output.exists()
        assert output.stat().st_size > 0

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib required")
    @pytest.mark.skipif(not HAS_NUMPY, reason="numpy required")
    def test_plot_coverage_track(self, tmp_path: Path) -> None:
        from metainformant.structural_variants.visualization.plots import plot_coverage_track

        coverage = [50 + i * 0.1 for i in range(100)]
        variants = [{"chrom": "chr1", "start": 30000, "end": 60000, "sv_type": "DEL"}]
        region = ("chr1", 0, 100000)

        output = tmp_path / "coverage.png"
        fig = plot_coverage_track(coverage, variants, region, output_path=str(output))

        assert output.exists()

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib required")
    def test_plot_circos(self, sample_variants: list[dict[str, Any]], tmp_path: Path) -> None:
        from metainformant.structural_variants.visualization.plots import plot_circos

        chromosomes = {"chr1": 100_000, "chr2": 80_000}
        output = tmp_path / "circos.png"

        fig = plot_circos(sample_variants, chromosomes, str(output))
        assert output.exists()
        assert output.stat().st_size > 0

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib required")
    def test_plot_cnv_profile(self, tmp_path: Path) -> None:
        from metainformant.structural_variants.visualization.plots import plot_cnv_profile

        segments = [
            {"chrom": "chr1", "start": 0, "end": 30000, "mean_log2ratio": 0.0, "state": "NEUTRAL"},
            {"chrom": "chr1", "start": 30000, "end": 50000, "mean_log2ratio": -0.7, "state": "DEL"},
            {"chrom": "chr1", "start": 50000, "end": 100000, "mean_log2ratio": 0.0, "state": "NEUTRAL"},
            {"chrom": "chr2", "start": 0, "end": 40000, "mean_log2ratio": 0.5, "state": "DUP"},
            {"chrom": "chr2", "start": 40000, "end": 80000, "mean_log2ratio": 0.0, "state": "NEUTRAL"},
        ]
        chromosomes = {"chr1": 100_000, "chr2": 80_000}

        output = tmp_path / "cnv_profile.png"
        fig = plot_cnv_profile(segments, chromosomes, output_path=str(output))

        assert output.exists()

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib required")
    def test_plot_breakpoint_detail(self, tmp_path: Path) -> None:
        from metainformant.structural_variants.visualization.plots import plot_breakpoint_detail

        variant = {
            "chrom": "chr1",
            "start": 5000,
            "end": 8000,
            "sv_type": "DEL",
            "evidence": {"evidence_reads": ["r1", "r2"]},
        }
        reads = [
            {"chrom": "chr1", "pos": 4800 + i * 50, "read_length": 150, "name": f"r{i}", "cigar": "150M", "mapq": 60}
            for i in range(30)
        ]

        output = tmp_path / "breakpoint.png"
        fig = plot_breakpoint_detail(variant, reads, output_path=str(output))

        assert output.exists()


# ---------------------------------------------------------------------------
# Tests: Module-level imports
# ---------------------------------------------------------------------------


class TestModuleImports:
    """Test that the module imports correctly."""

    def test_import_main_module(self) -> None:
        from metainformant import structural_variants

        assert hasattr(structural_variants, "detection")
        assert hasattr(structural_variants, "annotation")
        assert hasattr(structural_variants, "filtering")
        assert hasattr(structural_variants, "visualization")

    def test_import_top_level_functions(self) -> None:
        from metainformant.structural_variants import (
            annotate_gene_overlap,
            call_structural_variants,
            detect_cnv_from_depth,
            filter_by_quality,
            merge_callsets,
            predict_functional_impact,
            refine_breakpoints,
        )

        assert callable(detect_cnv_from_depth)
        assert callable(call_structural_variants)
        assert callable(refine_breakpoints)
        assert callable(annotate_gene_overlap)
        assert callable(predict_functional_impact)
        assert callable(filter_by_quality)
        assert callable(merge_callsets)

    def test_import_sv_types(self) -> None:
        from metainformant.structural_variants.detection.sv_calling import SVType

        assert SVType.DEL.value == "DEL"
        assert SVType.DUP.value == "DUP"
        assert SVType.INV.value == "INV"
        assert SVType.TRA.value == "TRA"
        assert SVType.INS.value == "INS"

    def test_import_dataclasses(self) -> None:
        from metainformant.structural_variants.annotation.functional_impact import FunctionalImpact
        from metainformant.structural_variants.annotation.overlap import GenomicInterval, OverlapResult
        from metainformant.structural_variants.detection.breakpoints import Breakpoint, BreakpointPair
        from metainformant.structural_variants.detection.cnv import CNVResult, CNVSegment
        from metainformant.structural_variants.detection.sv_calling import StructuralVariant, SVEvidence
        from metainformant.structural_variants.filtering.merge import MergedVariant
        from metainformant.structural_variants.filtering.quality_filter import FilterStats

        # Verify they are constructible
        seg = CNVSegment(chrom="chr1", start=0, end=1000, mean_log2ratio=0.0, n_bins=10)
        assert seg.chrom == "chr1"

        gi = GenomicInterval(chrom="chr1", start=0, end=1000, name="test")
        assert gi.name == "test"
