"""Depth tests for structural-variant detection, breakpoint, and plot internals.

Covers CIGAR/SA-tag parsing, insert-size estimation, evidence clustering and
merging, genotyping thresholds, CNV state assignment with environment
overrides, breakpoint consensus internals, nearest-gene direction, and
empty-input plot handling — all with real data structures.
"""

from __future__ import annotations

import pytest

from metainformant.structural_variants.annotation.overlap import (
    GenomicInterval,
    IntervalIndex,
    find_nearest_gene,
)
from metainformant.structural_variants.detection.breakpoints import (
    Breakpoint,
    _consensus_position,
    _position_std,
    detect_microhomology,
    refine_breakpoints,
)
from metainformant.structural_variants.detection.cnv import (
    _assign_cnv_state,
    call_cnv_states,
    merge_adjacent_segments,
)
from metainformant.structural_variants.detection.sv_calling import (
    SVEvidence,
    SVType,
    _cluster_evidence,
    _estimate_insert_size,
    _merge_evidence,
    _parse_cigar_clips,
    _parse_sa_tag,
    genotype_sv,
)

try:
    import matplotlib

    matplotlib.use("Agg")
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


class TestCigarAndTagParsing:
    def test_parse_cigar_clips_string(self) -> None:
        clips = _parse_cigar_clips("25S50M30S", 1000)
        assert clips == [(1000, 25, "left"), (1050, 30, "right")]

    def test_parse_cigar_clips_tuple_form(self) -> None:
        clips = _parse_cigar_clips([(4, 25), (0, 50), (4, 30)], 1000)
        assert clips == [(1000, 25, "left"), (1050, 30, "right")]

    def test_parse_sa_tag_malformed_entries_skipped(self) -> None:
        parsed = _parse_sa_tag("chr2,notanint,+,50M,60,0;chr3,700,-,50M,60,0;")
        assert parsed == [("chr3", 699, "-")]


class TestInsertSizeEstimation:
    def test_empty_reads_use_illumina_defaults(self) -> None:
        stats = _estimate_insert_size([])
        assert stats.mean == 400.0
        assert stats.std == 100.0
        assert stats.median == 400.0

    def test_estimate_from_data(self) -> None:
        reads = [
            {"chrom": "chr1", "mate_chrom": "chr1", "insert_size": s}
            for s in (300, 400, 500)
        ]
        stats = _estimate_insert_size(reads)
        assert stats.mean == pytest.approx(400.0)
        assert stats.std == pytest.approx(100.0)
        assert stats.median == 400.0

    def test_interchromosomal_and_outliers_excluded(self) -> None:
        reads = [
            {"chrom": "chr1", "mate_chrom": "chr2", "insert_size": 400},
            {"chrom": "chr1", "mate_chrom": "chr1", "insert_size": 0},
            {"chrom": "chr1", "mate_chrom": "chr1", "insert_size": 50000},
        ]
        assert _estimate_insert_size(reads).mean == 400.0  # defaults, nothing sampled


class TestEvidenceClustering:
    @staticmethod
    def _ev(bp1: int, bp2: int, name: str, **kw: object) -> SVEvidence:
        return SVEvidence(
            split_reads=1,
            evidence_reads=[name],
            breakpoint1=bp1,
            breakpoint2=bp2,
            chrom1="chr1",
            chrom2="chr1",
            **kw,
        )

    def test_cluster_by_distance(self) -> None:
        evidence = [self._ev(1000, 5000, "a"), self._ev(1100, 5100, "b"), self._ev(20000, 30000, "c")]
        clusters = _cluster_evidence(evidence, max_distance=500)
        assert len(clusters) == 2
        assert sorted(len(c) for c in clusters) == [1, 2]

    def test_merge_sums_support_and_medians_breakpoints(self) -> None:
        e1 = self._ev(1000, 5000, "a", strand1="+", strand2="-")
        e2 = SVEvidence(
            discordant_pairs=1,
            evidence_reads=["b"],
            breakpoint1=1010,
            breakpoint2=4990,
            chrom1="chr1",
            chrom2="chr1",
            strand1="+",
            strand2="-",
        )
        merged = _merge_evidence([e1, e2])
        assert merged.split_reads == 1 and merged.discordant_pairs == 1
        assert merged.breakpoint1 == 1005
        assert merged.breakpoint2 == 4995
        assert set(merged.evidence_reads) == {"a", "b"}
        assert (merged.strand1, merged.strand2) == ("+", "-")

    def test_merge_single_item_is_identity(self) -> None:
        e = self._ev(1000, 5000, "a")
        assert _merge_evidence([e]) is e


class TestGenotypeSVThresholds:
    @staticmethod
    def _variant(split_reads: int, names: list[str]) -> object:
        from metainformant.structural_variants.detection.sv_calling import StructuralVariant

        return StructuralVariant(
            chrom="chr1",
            start=5000,
            end=8000,
            sv_type=SVType.DEL,
            evidence=SVEvidence(split_reads=split_reads, evidence_reads=names),
        )

    def test_low_allele_fraction_is_hom_ref(self) -> None:
        variant = self._variant(1, ["r1"])
        refs = [{"name": f"ref{i}", "chrom": "chr1", "pos": 4700, "read_length": 600, "mapq": 60} for i in range(9)]
        assert genotype_sv(variant, refs) == "0/0"  # AF = 1/10

    def test_high_allele_fraction_is_hom_alt(self) -> None:
        variant = self._variant(20, ["r1"])
        assert genotype_sv(variant, []) == "1/1"

    def test_no_support_is_missing(self) -> None:
        variant = self._variant(0, [])
        assert genotype_sv(variant, []) == "./."

    def test_low_mapq_reads_not_counted(self) -> None:
        variant = self._variant(1, ["r1"])
        refs = [{"name": "ref1", "chrom": "chr1", "pos": 4700, "read_length": 600, "mapq": 5}]
        # Only the single alt read counts: AF = 1/1
        assert genotype_sv(variant, refs, min_mapq=20) == "1/1"


class TestCnvStateAssignment:
    def test_state_boundaries(self) -> None:
        assert _assign_cnv_state(-2.0) == ("HOMODEL", 0)
        assert _assign_cnv_state(-0.5) == ("DEL", 1)
        assert _assign_cnv_state(-0.1) == ("NEUTRAL", 2)
        assert _assign_cnv_state(0.5) == ("DUP", 3)
        assert _assign_cnv_state(1.2) == ("AMP", 5)

    def test_call_cnv_states_preserves_segment_coords(self) -> None:
        from metainformant.structural_variants.detection.cnv import CNVSegment

        seg = CNVSegment(chrom="chr7", start=100, end=400, mean_log2ratio=-0.5, n_bins=3, confidence=0.5)
        result = call_cnv_states([seg])
        assert len(result) == 1
        out = result[0]
        assert (out.chrom, out.start, out.end) == ("chr7", 100, 400)
        assert out.state == "DEL"
        assert out.confidence == 0.0

    def test_env_threshold_overrides(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("SV_CNV_DEL_THRESHOLD", "-0.05")
        monkeypatch.setenv("SV_CNV_DUP_THRESHOLD", "0.4")
        result = call_cnv_states([(0, 10, -0.1), (10, 20, 0.45)])
        assert result[0].state == "DEL"
        assert result[1].state == "DUP"

    def test_merge_adjacent_empty(self) -> None:
        assert merge_adjacent_segments([]) == []


class TestBreakpointInternals:
    def test_consensus_position_mode_and_fallback(self) -> None:
        clips = [(1000, "r1"), (1000, "r2"), (1050, "r3")]
        pos, support, names = _consensus_position(clips, fallback=900)
        assert pos == 1000 and support == 2 and names == ["r1", "r2"]
        assert _consensus_position([], 900) == (900, 0, [])

    def test_position_std(self) -> None:
        assert _position_std([]) == 0.0
        assert _position_std([(1000, "a")]) == 0.0
        assert _position_std([(1000, "a"), (1000, "b")]) == 0.0
        assert _position_std([(1000, "a"), (1010, "b")]) == pytest.approx(7.0711, rel=1e-3)

    def test_refine_breakpoints_falls_back_without_clips(self) -> None:
        pairs = refine_breakpoints([{"chrom": "chr1", "start": 5000, "end": 8000, "sv_type": "DEL"}], [])
        assert len(pairs) == 1
        assert pairs[0].bp1.position == 5000
        assert pairs[0].bp2.position == 8000
        assert pairs[0].bp1.support == 0 and pairs[0].bp1.confidence == 0.0

    def test_detect_microhomology_accepts_breakpoint_object(self) -> None:
        bp = Breakpoint(chrom="chr1", position=5)
        # seq[mid-1] == seq[mid] and seq[mid-2] == seq[mid+1] => "TA"
        assert detect_microhomology(bp, "ATTA") == "TA"


class TestNearestGene:
    GENE_DB = [
        {"chrom": "chr1", "start": 1000, "end": 5000, "name": "GENE_A"},
        {"chrom": "chr1", "start": 6000, "end": 9000, "name": "GENE_B"},
    ]

    def test_downstream_direction(self) -> None:
        result = find_nearest_gene({"chrom": "chr1", "start": 10000, "end": 10100}, self.GENE_DB)
        assert result["nearest_gene"] == "GENE_B"
        assert result["distance"] == 1050
        assert result["direction"] == "downstream"

    def test_upstream_direction(self) -> None:
        result = find_nearest_gene({"chrom": "chr1", "start": 100, "end": 200}, self.GENE_DB)
        assert result["nearest_gene"] == "GENE_A"
        assert result["direction"] == "upstream"

    def test_query_nearest_inside_interval_distance_zero(self) -> None:
        index = IntervalIndex([GenomicInterval(chrom="chr1", start=100, end=200, name="A")])
        iv, dist = index.query_nearest("chr1", 150)
        assert dist == 0 and iv is not None and iv.name == "A"
        assert index.query_nearest("chr1", 1_000_000, max_distance=500) == (None, -1)

    def test_query_overlap_unknown_chromosome(self) -> None:
        index = IntervalIndex([GenomicInterval(chrom="chr1", start=100, end=200, name="A")])
        assert index.query_overlap("chrZ", 0, 1000) == []


@pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib required")
class TestPlotEmptyInput:
    def test_size_distribution_no_variants(self) -> None:
        from metainformant.structural_variants.visualization.plots import plot_sv_size_distribution

        fig = plot_sv_size_distribution([])
        assert fig is not None

    def test_type_summary_no_variants(self) -> None:
        from metainformant.structural_variants.visualization.plots import plot_sv_type_summary

        fig = plot_sv_type_summary([])
        assert fig is not None

    def test_size_distribution_linear_scale(self) -> None:
        from metainformant.structural_variants.visualization.plots import plot_sv_size_distribution

        variants = [{"sv_type": "DEL", "size": 5000}, {"sv_type": "DEL", "size": 8000}]
        fig = plot_sv_size_distribution(variants, log_scale=False)
        assert fig is not None
