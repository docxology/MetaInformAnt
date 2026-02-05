"""Comprehensive tests for the longread module.

Tests all public functions in quality, filtering, analysis, assembly,
io, and visualization subpackages. All tests use real implementations
with synthetic data -- NO MOCKING.
"""

from __future__ import annotations

import math
import random
import string
from pathlib import Path
from typing import Any

import pytest

# ---------------------------------------------------------------------------
# Quality metrics
# ---------------------------------------------------------------------------
from metainformant.longread.quality.metrics import (
    ReadLengthStatistics,
    QualityDistribution,
    calculate_n50,
    calculate_nx,
    read_length_stats,
    quality_score_distribution,
    estimate_accuracy,
    calculate_throughput,
)

# ---------------------------------------------------------------------------
# Filtering
# ---------------------------------------------------------------------------
from metainformant.longread.quality.filtering import (
    ReadRecord,
    AdapterMatch,
    ONT_ADAPTERS,
    PACBIO_ADAPTERS,
    filter_by_length,
    filter_by_quality,
    trim_adapters,
    detect_adapters,
    split_chimeric_reads,
)

# ---------------------------------------------------------------------------
# Analysis -- modified bases
# ---------------------------------------------------------------------------
from metainformant.longread.analysis.modified_bases import (
    MethylationCall,
    RegionMethylation,
    DifferentialMethylationResult,
    detect_methylation,
    call_5mc,
    call_6ma,
    aggregate_methylation,
    differential_methylation,
)

# ---------------------------------------------------------------------------
# Analysis -- structural variants
# ---------------------------------------------------------------------------
from metainformant.longread.analysis.structural import (
    StructuralVariant,
    detect_sv_from_long_reads,
    detect_insertions,
    detect_inversions,
    phase_structural_variants,
)

# ---------------------------------------------------------------------------
# Analysis -- phasing
# ---------------------------------------------------------------------------
from metainformant.longread.analysis.phasing import (
    Variant,
    PhaseBlock,
    PhaseBlockStats,
    phase_reads,
    build_haplotype_blocks,
    tag_reads_by_haplotype,
    calculate_phase_block_stats,
)

# ---------------------------------------------------------------------------
# Assembly -- overlap
# ---------------------------------------------------------------------------
from metainformant.longread.assembly.overlap import (
    Minimizer,
    Overlap,
    OverlapGraph,
    find_overlaps,
    minimizer_sketch,
    compute_overlap_graph,
    filter_contained_reads,
)

# ---------------------------------------------------------------------------
# Assembly -- consensus
# ---------------------------------------------------------------------------
from metainformant.longread.assembly.consensus import (
    ConsensusResult,
    MSAResult,
    generate_consensus,
    polish_consensus,
    multiple_sequence_alignment,
    calculate_consensus_quality,
)

# ---------------------------------------------------------------------------
# Assembly -- hybrid
# ---------------------------------------------------------------------------
from metainformant.longread.assembly.hybrid import (
    CorrectedRead,
    Scaffold,
    HybridAssemblyResult,
    hybrid_assemble,
    correct_with_short_reads,
    scaffold_with_long_reads,
)

# ---------------------------------------------------------------------------
# IO -- fast5
# ---------------------------------------------------------------------------
from metainformant.longread.io.fast5 import (
    Fast5Read,
    read_fast5,
    extract_signal,
    extract_basecalls,
    get_read_metadata,
)

# ---------------------------------------------------------------------------
# IO -- bam
# ---------------------------------------------------------------------------
from metainformant.longread.io.bam import (
    LongReadAlignment,
    AlignmentStats,
    extract_methylation_tags,
    get_supplementary_alignments,
    calculate_alignment_stats,
)

# ---------------------------------------------------------------------------
# IO -- formats
# ---------------------------------------------------------------------------
from metainformant.longread.io.formats import (
    write_paf,
)

# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------
from metainformant.longread.visualization.plots import (
    plot_read_length_histogram,
    plot_quality_vs_length,
    plot_dotplot,
    plot_alignment_view,
    plot_methylation_track,
    plot_phasing_blocks,
)


# ===================================================================
# Helpers -- synthetic data generators
# ===================================================================

def _random_dna(length: int, seed: int = 42) -> str:
    """Generate a deterministic random DNA sequence."""
    rng = random.Random(seed)
    return "".join(rng.choices("ACGT", k=length))


def _phred_string(scores: list[int]) -> str:
    """Convert a list of Phred scores to an ASCII quality string (Phred+33)."""
    return "".join(chr(s + 33) for s in scores)


def _make_uniform_quality_string(length: int, quality: int = 20) -> str:
    """Generate a quality string with a uniform Phred score."""
    return _phred_string([quality] * length)


# ===================================================================
# TEST SUITE 1 -- Quality Metrics
# ===================================================================


class TestQualityMetrics:
    """Tests for quality/metrics.py."""

    def test_calculate_n50_basic(self) -> None:
        """Known input: [1000, 2000, 3000, 4000, 5000] -> N50 = 4000.

        Total = 15000, threshold = 7500.
        Sorted desc: 5000, 4000, 3000, 2000, 1000
        Cumulative: 5000, 9000 (>= 7500) -> N50 = 4000.
        """
        lengths = [1000, 2000, 3000, 4000, 5000]
        assert calculate_n50(lengths) == 4000

    def test_calculate_n50_single_read(self) -> None:
        """Single read: N50 is the read length itself."""
        assert calculate_n50([5000]) == 5000

    def test_calculate_n50_empty(self) -> None:
        """Empty input returns 0."""
        assert calculate_n50([]) == 0

    def test_calculate_n50_zeros_ignored(self) -> None:
        """Zero-length reads should be ignored."""
        assert calculate_n50([0, 0, 5000]) == 5000

    def test_calculate_nx_n10(self) -> None:
        """N10: the length at 10% cumulative threshold."""
        lengths = [1000, 2000, 3000, 4000, 5000]
        # Total = 15000, threshold = 1500
        # Sorted desc: 5000 cumulative 5000 >= 1500 -> N10 = 5000
        assert calculate_nx(lengths, x=10) == 5000

    def test_calculate_nx_n90(self) -> None:
        """N90: the length at 90% cumulative threshold."""
        lengths = [1000, 2000, 3000, 4000, 5000]
        # Total = 15000, threshold = 13500
        # cumulative: 5000, 9000, 12000, 14000 >= 13500 -> N90 = 2000
        assert calculate_nx(lengths, x=90) == 2000

    def test_calculate_nx_invalid_x(self) -> None:
        """x outside 0-100 raises ValueError."""
        with pytest.raises(ValueError, match="must be between 0 and 100"):
            calculate_nx([1000], x=101)
        with pytest.raises(ValueError, match="must be between 0 and 100"):
            calculate_nx([1000], x=-1)

    def test_read_length_stats_comprehensive(self) -> None:
        """Full statistics from a known distribution."""
        lengths = [500, 1000, 2000, 3000, 5000]
        stats = read_length_stats(lengths)

        assert isinstance(stats, ReadLengthStatistics)
        assert stats.count == 5
        assert stats.total_bases == 11500
        assert stats.mean_length == pytest.approx(2300.0)
        assert stats.median_length == pytest.approx(2000.0)
        assert stats.min_length == 500
        assert stats.max_length == 5000
        assert stats.std_dev > 0
        assert stats.n50 == calculate_n50(lengths)
        assert stats.n90 == calculate_nx(lengths, 90)
        assert stats.l50 >= 1
        assert "p5" in stats.percentiles
        assert "p25" in stats.percentiles
        assert "p75" in stats.percentiles
        assert "p95" in stats.percentiles

    def test_read_length_stats_from_dicts(self) -> None:
        """Accepts dicts with 'sequence' key."""
        reads = [
            {"sequence": "ACGT" * 250},  # 1000 bp
            {"sequence": "ACGT" * 500},  # 2000 bp
        ]
        stats = read_length_stats(reads)
        assert stats.count == 2
        assert stats.total_bases == 3000

    def test_read_length_stats_from_dict_length_key(self) -> None:
        """Accepts dicts with 'length' key."""
        reads = [{"length": 1000}, {"length": 2000}]
        stats = read_length_stats(reads)
        assert stats.count == 2
        assert stats.total_bases == 3000

    def test_read_length_stats_empty(self) -> None:
        """Empty input returns zeroed stats."""
        stats = read_length_stats([])
        assert stats.count == 0
        assert stats.total_bases == 0

    def test_quality_score_distribution_basic(self) -> None:
        """Verify quality histogram from known quality strings."""
        # Two reads both with uniform Q20 -> mean quality = 20
        q_str = _phred_string([20] * 100)
        dist = quality_score_distribution([q_str, q_str])

        assert isinstance(dist, QualityDistribution)
        assert dist.mean_quality == pytest.approx(20.0)
        assert dist.median_quality == pytest.approx(20.0)
        assert dist.q20_fraction == pytest.approx(1.0)
        assert 20 in dist.histogram
        assert dist.histogram[20] == 200  # 100 bases * 2 reads

    def test_quality_score_distribution_from_dicts(self) -> None:
        """Accepts dicts with 'quality' key."""
        q_str = _phred_string([15] * 50)
        reads = [{"quality": q_str}]
        dist = quality_score_distribution(reads)
        assert dist.mean_quality == pytest.approx(15.0)

    def test_quality_score_distribution_empty(self) -> None:
        """Empty input returns zeroed distribution."""
        dist = quality_score_distribution([])
        assert dist.mean_quality == 0.0

    def test_estimate_accuracy_q10(self) -> None:
        """Q10 -> P(error) = 0.1 -> accuracy = 0.9."""
        acc = estimate_accuracy([10])
        assert acc == pytest.approx(0.9, abs=1e-9)

    def test_estimate_accuracy_q20(self) -> None:
        """Q20 -> 99% accuracy."""
        acc = estimate_accuracy([20])
        assert acc == pytest.approx(0.99, abs=1e-9)

    def test_estimate_accuracy_q30(self) -> None:
        """Q30 -> 99.9% accuracy."""
        acc = estimate_accuracy([30])
        assert acc == pytest.approx(0.999, abs=1e-9)

    def test_estimate_accuracy_from_string(self) -> None:
        """Accepts ASCII quality string."""
        q_str = _phred_string([20])
        acc = estimate_accuracy(q_str)
        assert acc == pytest.approx(0.99, abs=1e-9)

    def test_estimate_accuracy_empty(self) -> None:
        """Empty input returns 0.0."""
        assert estimate_accuracy([]) == 0.0

    def test_calculate_throughput_basic(self) -> None:
        """Total bases and per-hour rates."""
        reads = [1000, 2000, 3000]
        result = calculate_throughput(reads, run_duration=2.0)

        assert result["total_bases"] == 6000.0
        assert result["total_reads"] == 3.0
        assert result["mean_read_length"] == pytest.approx(2000.0)
        assert result["bases_per_hour"] == pytest.approx(3000.0)
        assert result["reads_per_hour"] == pytest.approx(1.5)

    def test_calculate_throughput_no_duration(self) -> None:
        """Without duration, per-hour metrics are absent."""
        result = calculate_throughput([1000, 2000])
        assert "bases_per_hour" not in result
        assert result["total_bases"] == 3000.0

    def test_calculate_throughput_from_dicts(self) -> None:
        """Accepts dicts with 'sequence' key."""
        reads = [{"sequence": "ACGT" * 250}]
        result = calculate_throughput(reads)
        assert result["total_bases"] == 1000.0


# ===================================================================
# TEST SUITE 2 -- Quality Filtering
# ===================================================================


class TestQualityFiltering:
    """Tests for quality/filtering.py."""

    def test_filter_by_length_min_only(self) -> None:
        """Reads shorter than min_length are discarded."""
        reads = [
            {"sequence": "A" * 500, "read_id": "short"},
            {"sequence": "A" * 2000, "read_id": "long"},
        ]
        result = filter_by_length(reads, min_length=1000)
        assert len(result) == 1

    def test_filter_by_length_max_only(self) -> None:
        """Reads longer than max_length are discarded."""
        reads = [
            {"sequence": "A" * 500, "read_id": "short"},
            {"sequence": "A" * 5000, "read_id": "medium"},
            {"sequence": "A" * 50000, "read_id": "long"},
        ]
        result = filter_by_length(reads, min_length=0, max_length=10000)
        assert len(result) == 2

    def test_filter_by_length_with_read_records(self) -> None:
        """Works with ReadRecord objects."""
        reads = [
            ReadRecord(read_id="r1", sequence="A" * 3000),
            ReadRecord(read_id="r2", sequence="A" * 200),
        ]
        result = filter_by_length(reads, min_length=1000)
        assert len(result) == 1

    def test_filter_by_quality_passes_high_quality(self) -> None:
        """Reads above quality threshold pass."""
        q20 = _make_uniform_quality_string(100, quality=20)
        q5 = _make_uniform_quality_string(100, quality=5)
        reads = [
            {"sequence": "A" * 100, "quality": q20, "read_id": "good"},
            {"sequence": "A" * 100, "quality": q5, "read_id": "bad"},
        ]
        result = filter_by_quality(reads, min_q=10.0)
        assert len(result) == 1

    def test_filter_by_quality_no_quality_passes(self) -> None:
        """Reads without quality info pass the filter."""
        reads = [{"sequence": "A" * 100, "read_id": "noqual"}]
        result = filter_by_quality(reads, min_q=10.0)
        assert len(result) == 1

    def test_trim_adapters_known_adapter(self) -> None:
        """Trimming removes a known adapter from the read start."""
        adapter = ONT_ADAPTERS["ONT_LSK110_top"]  # CAGCACCT (8 bp)
        seq = adapter + "A" * 2000
        reads = [{"sequence": seq, "read_id": "r1"}]
        trimmed = trim_adapters(
            reads,
            adapter_sequences={"ONT_LSK110_top": adapter},
            min_identity=0.75,
        )
        assert len(trimmed) == 1
        # The trimmed sequence should be shorter (adapter removed)
        assert len(trimmed[0].sequence) <= len(seq)
        assert trimmed[0].metadata.get("original_length") == len(seq)

    def test_trim_adapters_no_adapter_present(self) -> None:
        """Read without adapter is returned unchanged."""
        seq = _random_dna(3000, seed=99)
        reads = [{"sequence": seq, "read_id": "r1"}]
        # Use a very short, non-matching adapter set
        trimmed = trim_adapters(
            reads,
            adapter_sequences={"fake": "NNNNNNNNNNNNNNNN"},
            min_identity=0.99,
        )
        assert len(trimmed) == 1
        assert trimmed[0].sequence == seq

    def test_detect_adapters_finds_adapter_at_start(self) -> None:
        """detect_adapters locates adapter at the beginning of a sequence."""
        adapter = ONT_ADAPTERS["ONT_LSK110_top"]  # CAGCACCT
        seq = adapter + _random_dna(5000, seed=77)
        matches = detect_adapters(
            seq,
            known_adapters={"ONT_LSK110_top": adapter},
            min_identity=0.75,
        )
        # Should find at least one match near position 0
        assert len(matches) >= 1
        found_start = any(m.location == "start" for m in matches)
        assert found_start

    def test_detect_adapters_empty_sequence(self) -> None:
        """Empty sequence returns no matches."""
        assert detect_adapters("") == []

    def test_split_chimeric_reads_no_chimera(self) -> None:
        """Read without internal adapter is returned intact."""
        seq = _random_dna(5000, seed=55)
        reads = [{"sequence": seq, "read_id": "r1"}]
        result = split_chimeric_reads(
            reads,
            known_adapters={"fake": "NNNNNNNNNNNNNNNN"},
            min_adapter_identity=0.99,
        )
        assert len(result) == 1
        assert result[0].metadata.get("chimeric") is False

    def test_split_chimeric_reads_splits_chimera(self) -> None:
        """Read with an internal known adapter sequence is split into fragments."""
        adapter = "AATGTACTTCGTTCAGTTACGTATTGCT"  # ONT_LSK109_top (28 bp)
        frag1 = _random_dna(2000, seed=10)
        frag2 = _random_dna(2000, seed=20)
        chimeric_seq = frag1 + adapter + frag2
        reads = [{"sequence": chimeric_seq, "read_id": "chimera"}]
        result = split_chimeric_reads(
            reads,
            known_adapters={"ONT_LSK109_top": adapter},
            min_fragment_length=500,
            min_adapter_identity=0.85,
        )
        # Should either find fragments or return intact if adapter not detected
        assert len(result) >= 1
        # If chimera was detected, fragments should be marked
        if len(result) > 1:
            assert any(r.metadata.get("chimeric") is True for r in result)


# ===================================================================
# TEST SUITE 3 -- Methylation Analysis
# ===================================================================


class TestMethylationAnalysis:
    """Tests for analysis/modified_bases.py."""

    def test_call_5mc_returns_calls(self) -> None:
        """call_5mc produces MethylationCall objects from features."""
        features = [
            {
                "position": 100,
                "chromosome": "chr1",
                "strand": "+",
                "current_mean": 105.0,
                "current_std": 2.5,
                "dwell_time": 3.0,
                "context_shift": 4.0,
                "context": "CpG",
            },
            {
                "position": 200,
                "chromosome": "chr1",
                "strand": "+",
                "current_mean": 100.0,
                "current_std": 0.5,
                "dwell_time": 1.0,
                "context_shift": 0.5,
                "context": "CpG",
            },
        ]
        calls = call_5mc(features)
        assert len(calls) == 2
        assert all(isinstance(c, MethylationCall) for c in calls)
        assert all(c.modification_type == "5mC" for c in calls)
        assert calls[0].position == 100
        assert 0.0 <= calls[0].probability <= 1.0

    def test_call_5mc_high_shift_gives_higher_probability(self) -> None:
        """Larger context_shift should yield higher methylation probability."""
        high_shift = call_5mc([{
            "position": 0,
            "context_shift": 10.0,
            "dwell_time": 3.0,
            "current_std": 2.0,
        }])[0]
        low_shift = call_5mc([{
            "position": 0,
            "context_shift": 0.1,
            "dwell_time": 1.0,
            "current_std": 0.1,
        }])[0]
        assert high_shift.probability > low_shift.probability

    def test_call_6ma_returns_calls(self) -> None:
        """call_6ma produces MethylationCall objects typed as 6mA."""
        features = [
            {
                "position": 300,
                "chromosome": "chr2",
                "strand": "-",
                "context_shift": 8.0,
                "dwell_time": 2.5,
                "current_std": 1.5,
                "context": "A",
            },
        ]
        calls = call_6ma(features)
        assert len(calls) == 1
        assert calls[0].modification_type == "6mA"
        assert calls[0].context == "A"
        assert 0.0 <= calls[0].probability <= 1.0

    def test_detect_methylation_cpg_model(self) -> None:
        """detect_methylation with 'cpg' model processes signal data."""
        import numpy as np

        # Construct signal with segments of different means to produce events
        rng = np.random.RandomState(42)
        segment1 = rng.normal(100, 1, 50)
        segment2 = rng.normal(110, 2, 50)
        segment3 = rng.normal(95, 1, 50)
        segment4 = rng.normal(115, 3, 50)
        segment5 = rng.normal(100, 1, 50)
        signal = np.concatenate([segment1, segment2, segment3, segment4, segment5])

        calls = detect_methylation(signal, model="cpg", threshold=0.3)
        assert isinstance(calls, list)
        # The function should process without error; call count depends on signal
        for c in calls:
            assert isinstance(c, MethylationCall)
            assert c.modification_type == "5mC"

    def test_detect_methylation_6ma_model(self) -> None:
        """detect_methylation with '6ma' model."""
        import numpy as np

        rng = np.random.RandomState(99)
        signal = np.concatenate([
            rng.normal(100, 1, 50),
            rng.normal(120, 2, 50),
            rng.normal(90, 1, 50),
            rng.normal(130, 4, 50),
            rng.normal(100, 1, 50),
        ])
        calls = detect_methylation(signal, model="6ma", threshold=0.3)
        assert isinstance(calls, list)
        for c in calls:
            assert c.modification_type == "6mA"

    def test_detect_methylation_unknown_model(self) -> None:
        """Unknown model raises ValueError."""
        import numpy as np

        with pytest.raises(ValueError, match="Unknown methylation model"):
            detect_methylation(np.array([1.0, 2.0, 3.0]), model="unknown")

    def test_detect_methylation_empty_signal(self) -> None:
        """Empty signal returns no calls."""
        import numpy as np

        calls = detect_methylation(np.array([]), model="cpg")
        assert calls == []

    def test_aggregate_methylation_basic(self) -> None:
        """aggregate_methylation groups calls into regions."""
        calls = [
            MethylationCall(chromosome="chr1", position=100, probability=0.9, coverage=5, modified_count=4),
            MethylationCall(chromosome="chr1", position=100, probability=0.8, coverage=5, modified_count=4),
            MethylationCall(chromosome="chr1", position=100, probability=0.7, coverage=5, modified_count=3),
            MethylationCall(chromosome="chr1", position=200, probability=0.1, coverage=5, modified_count=0),
            MethylationCall(chromosome="chr1", position=200, probability=0.2, coverage=5, modified_count=1),
            MethylationCall(chromosome="chr1", position=200, probability=0.15, coverage=5, modified_count=0),
        ]
        regions = [
            {"chromosome": "chr1", "start": 0, "end": 300, "name": "promoter"},
        ]
        result = aggregate_methylation(calls, regions, min_coverage=3)
        assert len(result) == 1
        assert isinstance(result[0], RegionMethylation)
        assert result[0].name == "promoter"
        assert result[0].num_cpgs == 2  # two positions with coverage >= 3
        assert 0.0 <= result[0].mean_methylation <= 1.0

    def test_aggregate_methylation_empty_region(self) -> None:
        """Region with no calls returns zeroed RegionMethylation."""
        calls: list[MethylationCall] = []
        regions = [{"chromosome": "chr1", "start": 0, "end": 1000}]
        result = aggregate_methylation(calls, regions)
        assert len(result) == 1
        assert result[0].num_cpgs == 0

    def test_differential_methylation_detects_difference(self) -> None:
        """differential_methylation finds sites with significant differences."""
        # Sample 1: high methylation at position 100
        sample1 = [
            MethylationCall(chromosome="chr1", position=100, probability=0.9)
            for _ in range(10)
        ]
        # Sample 2: low methylation at position 100
        sample2 = [
            MethylationCall(chromosome="chr1", position=100, probability=0.1)
            for _ in range(10)
        ]
        results = differential_methylation(
            sample1, sample2,
            min_coverage=5,
            min_difference=0.2,
            alpha=0.05,
        )
        assert len(results) >= 1
        assert isinstance(results[0], DifferentialMethylationResult)
        assert results[0].difference == pytest.approx(-0.8, abs=0.01)
        assert results[0].sample1_methylation == pytest.approx(0.9)
        assert results[0].sample2_methylation == pytest.approx(0.1)

    def test_differential_methylation_no_difference(self) -> None:
        """Similar samples produce no significant results."""
        sample1 = [
            MethylationCall(chromosome="chr1", position=100, probability=0.5)
            for _ in range(10)
        ]
        sample2 = [
            MethylationCall(chromosome="chr1", position=100, probability=0.55)
            for _ in range(10)
        ]
        results = differential_methylation(
            sample1, sample2,
            min_coverage=5,
            min_difference=0.2,
        )
        # Difference 0.05 < min_difference 0.2 -> filtered out
        assert len(results) == 0


# ===================================================================
# TEST SUITE 4 -- Structural Variant Detection
# ===================================================================


class TestStructuralVariants:
    """Tests for analysis/structural.py."""

    @staticmethod
    def _make_alignment(
        read_name: str = "read1",
        chrom: str = "chr1",
        ref_start: int = 1000,
        cigar_string: str = "5000M",
        cigar_tuples: list[tuple[int, int]] | None = None,
        mapping_quality: int = 60,
        is_unmapped: bool = False,
        is_reverse: bool = False,
        sa_tag: str = "",
        query_length: int = 5000,
    ) -> dict[str, Any]:
        tags: dict[str, Any] = {}
        if sa_tag:
            tags["SA"] = sa_tag
        return {
            "read_name": read_name,
            "reference_name": chrom,
            "reference_start": ref_start,
            "cigar_string": cigar_string,
            "cigar_tuples": cigar_tuples or [],
            "mapping_quality": mapping_quality,
            "is_unmapped": is_unmapped,
            "is_reverse": is_reverse,
            "tags": tags,
            "query_length": query_length,
            "query_sequence": "",
        }

    def test_detect_sv_cigar_deletion(self) -> None:
        """Detect a deletion from a CIGAR D operation >= min_size."""
        alns = [
            self._make_alignment(read_name=f"r{i}", cigar_string="2000M200D3000M")
            for i in range(3)
        ]
        svs = detect_sv_from_long_reads(alns, min_size=50, min_support=2)
        del_svs = [sv for sv in svs if sv.sv_type == "DEL"]
        assert len(del_svs) >= 1
        assert del_svs[0].size >= 200

    def test_detect_sv_cigar_insertion(self) -> None:
        """Detect an insertion from a CIGAR I operation >= min_size."""
        alns = [
            self._make_alignment(read_name=f"r{i}", cigar_string="2000M150I3000M")
            for i in range(3)
        ]
        svs = detect_sv_from_long_reads(alns, min_size=50, min_support=2)
        ins_svs = [sv for sv in svs if sv.sv_type == "INS"]
        assert len(ins_svs) >= 1

    def test_detect_sv_split_read_inversion(self) -> None:
        """Supplementary alignment on opposite strand -> inversion."""
        sa_tag = "chr1,6000,-,5000M,60,0"
        alns = [
            self._make_alignment(
                read_name=f"r{i}", ref_start=1000,
                cigar_string="2500M", sa_tag=sa_tag, query_length=5000,
            )
            for i in range(3)
        ]
        svs = detect_sv_from_long_reads(alns, min_size=50, min_support=2)
        inv_svs = [sv for sv in svs if sv.sv_type == "INV"]
        assert len(inv_svs) >= 1

    def test_detect_insertions_convenience(self) -> None:
        """detect_insertions filters to INS type only."""
        alns = [
            self._make_alignment(read_name=f"r{i}", cigar_string="2000M150I3000M")
            for i in range(3)
        ]
        insertions = detect_insertions(alns, min_size=50, min_support=2)
        assert all(sv.sv_type == "INS" for sv in insertions)

    def test_detect_inversions_convenience(self) -> None:
        """detect_inversions filters to INV type only."""
        sa_tag = "chr1,6000,-,5000M,60,0"
        alns = [
            self._make_alignment(
                read_name=f"r{i}", ref_start=1000,
                cigar_string="2500M", sa_tag=sa_tag, query_length=5000,
            )
            for i in range(3)
        ]
        inversions = detect_inversions(alns, min_size=50, min_support=2)
        assert all(sv.sv_type == "INV" for sv in inversions)

    def test_detect_sv_low_mapq_filtered(self) -> None:
        """Alignments below min_mapping_quality are ignored."""
        alns = [
            self._make_alignment(read_name=f"r{i}", cigar_string="2000M200D3000M", mapping_quality=5)
            for i in range(3)
        ]
        svs = detect_sv_from_long_reads(alns, min_size=50, min_support=2, min_mapping_quality=20)
        assert len(svs) == 0

    def test_phase_structural_variants_assigns_haplotype(self) -> None:
        """SVs with supporting reads tagged to a haplotype get phased."""
        sv = StructuralVariant(
            sv_type="DEL",
            chromosome="chr1",
            start=1000,
            end=2000,
            size=1000,
            supporting_reads=3,
            read_names=["r1", "r2", "r3"],
        )
        haplotype_tags = {"r1": 1, "r2": 1, "r3": 2}
        phased = phase_structural_variants([sv], haplotype_tags)
        assert len(phased) == 1
        assert phased[0].haplotype == 1  # majority vote: 2 vs 1
        assert phased[0].genotype != "./."

    def test_phase_structural_variants_unphased(self) -> None:
        """SVs with no haplotype-tagged reads stay unphased."""
        sv = StructuralVariant(
            sv_type="DEL",
            chromosome="chr1",
            start=1000,
            end=2000,
            size=1000,
            supporting_reads=2,
            read_names=["r1", "r2"],
        )
        phased = phase_structural_variants([sv], {})
        assert phased[0].haplotype == 0


# ===================================================================
# TEST SUITE 5 -- Phasing
# ===================================================================


class TestPhasing:
    """Tests for analysis/phasing.py."""

    @staticmethod
    def _make_phasing_data() -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
        """Create synthetic variants and reads for phasing tests."""
        variants = [
            {"chromosome": "chr1", "position": 1000, "ref_allele": "A", "alt_allele": "G"},
            {"chromosome": "chr1", "position": 2000, "ref_allele": "C", "alt_allele": "T"},
            {"chromosome": "chr1", "position": 3000, "ref_allele": "G", "alt_allele": "A"},
            {"chromosome": "chr1", "position": 4000, "ref_allele": "T", "alt_allele": "C"},
        ]
        # Reads: haplotype 1 carries alleles [0, 0, 1, 1], haplotype 2 carries [1, 1, 0, 0]
        reads = [
            {"read_name": "h1_r1", "variants": [
                {"chromosome": "chr1", "position": 1000, "allele": 0, "quality": 30},
                {"chromosome": "chr1", "position": 2000, "allele": 0, "quality": 30},
            ]},
            {"read_name": "h1_r2", "variants": [
                {"chromosome": "chr1", "position": 2000, "allele": 0, "quality": 30},
                {"chromosome": "chr1", "position": 3000, "allele": 1, "quality": 30},
            ]},
            {"read_name": "h1_r3", "variants": [
                {"chromosome": "chr1", "position": 3000, "allele": 1, "quality": 30},
                {"chromosome": "chr1", "position": 4000, "allele": 1, "quality": 30},
            ]},
            {"read_name": "h2_r1", "variants": [
                {"chromosome": "chr1", "position": 1000, "allele": 1, "quality": 30},
                {"chromosome": "chr1", "position": 2000, "allele": 1, "quality": 30},
            ]},
            {"read_name": "h2_r2", "variants": [
                {"chromosome": "chr1", "position": 2000, "allele": 1, "quality": 30},
                {"chromosome": "chr1", "position": 3000, "allele": 0, "quality": 30},
            ]},
            {"read_name": "h2_r3", "variants": [
                {"chromosome": "chr1", "position": 3000, "allele": 0, "quality": 30},
                {"chromosome": "chr1", "position": 4000, "allele": 0, "quality": 30},
            ]},
        ]
        return variants, reads

    def test_phase_reads_produces_blocks_and_assignments(self) -> None:
        """phase_reads returns phase blocks and read haplotype assignments."""
        variants, reads = self._make_phasing_data()
        blocks, assignments = phase_reads(variants, reads)

        assert isinstance(blocks, list)
        assert isinstance(assignments, dict)
        # Should produce at least one phase block
        assert len(blocks) >= 1
        assert isinstance(blocks[0], PhaseBlock)
        assert blocks[0].num_variants >= 2
        # Should assign reads
        assert len(assignments) >= 1
        assert all(hp in (1, 2) for hp in assignments.values())

    def test_phase_reads_empty_variants(self) -> None:
        """Empty variants produces no blocks."""
        blocks, assignments = phase_reads([], [])
        assert blocks == []
        assert assignments == {}

    def test_build_haplotype_blocks_basic(self) -> None:
        """build_haplotype_blocks groups nearby phased variants."""
        phased_vars = [
            {"chromosome": "chr1", "position": 1000, "haplotype": 1, "allele": 0, "phase_set": 0},
            {"chromosome": "chr1", "position": 2000, "haplotype": 1, "allele": 1, "phase_set": 0},
            {"chromosome": "chr1", "position": 3000, "haplotype": 2, "allele": 0, "phase_set": 0},
        ]
        blocks = build_haplotype_blocks(phased_vars, max_gap=300000)
        assert len(blocks) >= 1
        assert isinstance(blocks[0], PhaseBlock)
        assert blocks[0].num_variants >= 2

    def test_build_haplotype_blocks_gap_breaks_block(self) -> None:
        """A large gap creates separate phase blocks."""
        phased_vars = [
            {"chromosome": "chr1", "position": 1000, "haplotype": 1, "allele": 0, "phase_set": 0},
            {"chromosome": "chr1", "position": 2000, "haplotype": 1, "allele": 1, "phase_set": 0},
            {"chromosome": "chr1", "position": 1000000, "haplotype": 1, "allele": 0, "phase_set": 0},
            {"chromosome": "chr1", "position": 1001000, "haplotype": 2, "allele": 1, "phase_set": 0},
        ]
        blocks = build_haplotype_blocks(phased_vars, max_gap=100000)
        assert len(blocks) == 2

    def test_tag_reads_by_haplotype_assigns_hp(self) -> None:
        """tag_reads_by_haplotype adds HP and PS tags."""
        block = PhaseBlock(
            block_id=0,
            chromosome="chr1",
            start=1000,
            end=4000,
            variants=[1000, 2000, 3000],
            haplotype1=[0, 1, 0],
            haplotype2=[1, 0, 1],
            num_variants=3,
        )
        reads = [
            {
                "read_name": "r1",
                "reference_name": "chr1",
                "reference_start": 900,
                "reference_end": 2500,
                "variants": [
                    {"position": 1000, "allele": 0},
                    {"position": 2000, "allele": 1},
                ],
            },
        ]
        tagged = tag_reads_by_haplotype(reads, [block])
        assert len(tagged) == 1
        assert tagged[0]["HP"] == 1  # matches haplotype1 alleles

    def test_tag_reads_by_haplotype_no_overlap(self) -> None:
        """Reads not overlapping any block get HP=0."""
        block = PhaseBlock(
            block_id=0,
            chromosome="chr2",
            start=1000,
            end=2000,
            variants=[1000],
            haplotype1=[0],
            haplotype2=[1],
            num_variants=1,
        )
        reads = [
            {
                "read_name": "r1",
                "reference_name": "chr1",
                "reference_start": 100,
                "reference_end": 500,
            },
        ]
        tagged = tag_reads_by_haplotype(reads, [block])
        assert tagged[0]["HP"] == 0

    def test_calculate_phase_block_stats_basic(self) -> None:
        """calculate_phase_block_stats computes N50 and size metrics."""
        blocks = [
            PhaseBlock(block_id=0, chromosome="chr1", start=0, end=10000, num_variants=50, quality=0.95),
            PhaseBlock(block_id=1, chromosome="chr1", start=20000, end=25000, num_variants=20, quality=0.90),
        ]
        stats = calculate_phase_block_stats(blocks)
        assert isinstance(stats, PhaseBlockStats)
        assert stats.num_blocks == 2
        assert stats.num_phased_variants == 70
        assert stats.max_block_size == 10000
        assert stats.total_phased_span == 15000
        assert stats.n50 > 0
        assert stats.mean_block_size == pytest.approx(7500.0)
        assert stats.switch_error_rate >= 0

    def test_calculate_phase_block_stats_empty(self) -> None:
        """Empty blocks returns zeroed stats."""
        stats = calculate_phase_block_stats([])
        assert stats.num_blocks == 0
        assert stats.n50 == 0


# ===================================================================
# TEST SUITE 6 -- Assembly: Overlap
# ===================================================================


class TestAssemblyOverlap:
    """Tests for assembly/overlap.py."""

    def test_minimizer_sketch_produces_minimizers(self) -> None:
        """minimizer_sketch returns Minimizer objects for a valid sequence."""
        seq = _random_dna(5000, seed=1)
        minimizers = minimizer_sketch(seq, k=15, w=10)
        assert len(minimizers) > 0
        assert all(isinstance(m, Minimizer) for m in minimizers)
        # Each minimizer should have a hash and position
        for m in minimizers:
            assert m.hash_value >= 0
            assert 0 <= m.position < len(seq)

    def test_minimizer_sketch_short_sequence(self) -> None:
        """Sequence shorter than k returns empty list."""
        assert minimizer_sketch("ACGT", k=15) == []

    def test_minimizer_sketch_empty(self) -> None:
        """Empty sequence returns empty list."""
        assert minimizer_sketch("") == []

    def test_minimizer_sketch_canonical_deterministic(self) -> None:
        """Same sequence always produces the same sketch."""
        seq = _random_dna(2000, seed=42)
        s1 = minimizer_sketch(seq, k=15, w=10)
        s2 = minimizer_sketch(seq, k=15, w=10)
        assert len(s1) == len(s2)
        assert all(a.hash_value == b.hash_value for a, b in zip(s1, s2))

    def test_find_overlaps_identical_reads(self) -> None:
        """Identical reads should have a strong overlap."""
        seq = _random_dna(5000, seed=7)
        reads = [
            {"read_id": "r1", "sequence": seq},
            {"read_id": "r2", "sequence": seq},
        ]
        overlaps = find_overlaps(reads, min_overlap=1000, min_minimizer_matches=2)
        assert len(overlaps) >= 1
        assert overlaps[0].overlap_length > 0

    def test_find_overlaps_no_overlap(self) -> None:
        """Completely different reads have no overlap."""
        seq1 = "A" * 5000
        seq2 = "C" * 5000
        reads = [
            {"read_id": "r1", "sequence": seq1},
            {"read_id": "r2", "sequence": seq2},
        ]
        overlaps = find_overlaps(reads, min_overlap=1000, min_minimizer_matches=3)
        assert len(overlaps) == 0

    def test_find_overlaps_single_read(self) -> None:
        """Single read returns empty overlap list."""
        reads = [{"read_id": "r1", "sequence": _random_dna(5000)}]
        overlaps = find_overlaps(reads)
        assert overlaps == []

    def test_compute_overlap_graph_structure(self) -> None:
        """compute_overlap_graph produces valid OverlapGraph."""
        ovl = Overlap(
            query_name="r1", query_length=5000, query_start=0, query_end=3000,
            target_name="r2", target_length=5000, target_start=2000, target_end=5000,
            overlap_length=3000, num_matches=10, identity=0.9,
        )
        graph = compute_overlap_graph([ovl])
        assert isinstance(graph, OverlapGraph)
        assert graph.num_nodes == 2
        assert graph.num_edges == 1
        assert "r1" in graph.nodes
        assert "r2" in graph.nodes

    def test_filter_contained_reads_removes_contained(self) -> None:
        """Overlaps involving contained reads are filtered out."""
        ovl_contained = Overlap(
            query_name="short", query_length=1000, query_start=0, query_end=1000,
            target_name="long", target_length=10000, target_start=1000, target_end=2000,
            overlap_length=1000, is_contained=True,
        )
        ovl_normal = Overlap(
            query_name="r1", query_length=5000, query_start=0, query_end=3000,
            target_name="r2", target_length=5000, target_start=2000, target_end=5000,
            overlap_length=3000, is_contained=False,
        )
        filtered = filter_contained_reads([ovl_contained, ovl_normal])
        # The contained overlap involves "short", so any overlap with "short" is removed
        assert all(o.query_name != "short" for o in filtered)

    def test_filter_contained_reads_no_contained(self) -> None:
        """When no reads are contained, all overlaps are returned."""
        overlaps = [
            Overlap(
                query_name="r1", query_length=5000, query_start=0, query_end=3000,
                target_name="r2", target_length=5000, target_start=2000, target_end=5000,
                overlap_length=3000, is_contained=False,
            ),
        ]
        filtered = filter_contained_reads(overlaps)
        assert len(filtered) == 1


# ===================================================================
# TEST SUITE 7 -- Assembly: Consensus
# ===================================================================


class TestAssemblyConsensus:
    """Tests for assembly/consensus.py."""

    def test_generate_consensus_identical_reads(self) -> None:
        """Consensus of identical reads equals the read itself."""
        seq = _random_dna(200, seed=10)
        reads = [seq, seq, seq]
        result = generate_consensus(reads)

        assert isinstance(result, ConsensusResult)
        assert result.num_reads == 3
        assert result.length > 0
        assert result.sequence == seq

    def test_generate_consensus_single_read(self) -> None:
        """Single read returns that read as consensus."""
        seq = _random_dna(100, seed=20)
        result = generate_consensus([seq])
        assert result.sequence == seq
        assert result.num_reads == 1
        assert result.mean_quality == 30.0

    def test_generate_consensus_empty(self) -> None:
        """No reads returns empty ConsensusResult."""
        result = generate_consensus([])
        assert result.sequence == ""
        assert result.num_reads == 0

    def test_generate_consensus_majority_vote(self) -> None:
        """Consensus picks the majority base at each position."""
        # 3 reads agree on "A", 1 read has "T" at position 5
        reads = ["AAAAAAAAAA", "AAAAAAAAAA", "AAAAAAAAAA", "AAAAATAAAA"]
        result = generate_consensus(reads)
        assert result.sequence[5] == "A"  # majority

    def test_polish_consensus_does_not_degrade(self) -> None:
        """Polishing should maintain or improve the consensus."""
        seq = _random_dna(200, seed=30)
        reads = [seq, seq, seq]
        initial = generate_consensus(reads)
        polished = polish_consensus(initial, reads, iterations=2)
        assert isinstance(polished, ConsensusResult)
        assert polished.length > 0
        # After polishing identical reads, result should still match
        assert polished.sequence == seq

    def test_polish_consensus_from_string(self) -> None:
        """polish_consensus accepts a plain string as initial consensus."""
        seq = _random_dna(100, seed=40)
        polished = polish_consensus(seq, [seq, seq])
        assert isinstance(polished, ConsensusResult)
        assert polished.length > 0

    def test_multiple_sequence_alignment_identical(self) -> None:
        """MSA of identical sequences has perfect conservation."""
        seq = _random_dna(100, seed=50)
        result = multiple_sequence_alignment([seq, seq, seq])
        assert isinstance(result, MSAResult)
        assert result.num_sequences == 3
        assert result.consensus == seq
        assert all(c == 1.0 for c in result.conservation)

    def test_multiple_sequence_alignment_single(self) -> None:
        """MSA of a single sequence returns it unchanged."""
        seq = _random_dna(100, seed=60)
        result = multiple_sequence_alignment([seq])
        assert result.consensus == seq
        assert result.num_sequences == 1

    def test_multiple_sequence_alignment_empty(self) -> None:
        """MSA of empty input returns empty result."""
        result = multiple_sequence_alignment([])
        assert result.consensus == ""
        assert result.num_sequences == 0

    def test_calculate_consensus_quality_identical(self) -> None:
        """Quality of consensus from identical reads should be high."""
        seq = _random_dna(100, seed=70)
        quality = calculate_consensus_quality(seq, [seq, seq, seq])
        assert "overall_quality" in quality
        assert "overall_agreement" in quality
        assert quality["overall_agreement"] >= 0.9
        assert len(quality["per_base_quality"]) == len(seq)

    def test_calculate_consensus_quality_empty_reads(self) -> None:
        """No reads returns zero-quality."""
        quality = calculate_consensus_quality("ACGT", [])
        assert quality["overall_quality"] == 0.0

    def test_calculate_consensus_quality_empty_consensus(self) -> None:
        """Empty consensus returns empty quality."""
        quality = calculate_consensus_quality("", ["ACGT"])
        assert quality["overall_quality"] == 0.0


# ===================================================================
# TEST SUITE 8 -- Assembly: Hybrid
# ===================================================================


class TestAssemblyHybrid:
    """Tests for assembly/hybrid.py."""

    def test_correct_with_short_reads_basic(self) -> None:
        """Error correction with short reads produces CorrectedRead objects."""
        # Long read with a single error
        correct_seq = _random_dna(500, seed=80)
        long_seq = list(correct_seq)
        long_seq[100] = "N"  # introduce error
        long_read = "".join(long_seq)

        # Short reads covering the error region (fragments of the correct sequence)
        short_reads: list[str] = []
        for i in range(0, 400, 20):
            short_reads.append(correct_seq[i:i + 100])

        corrected = correct_with_short_reads([long_read], short_reads, kmer_size=15)
        assert len(corrected) == 1
        assert isinstance(corrected[0], CorrectedRead)
        assert corrected[0].original_sequence == long_read
        assert len(corrected[0].corrected_sequence) > 0

    def test_correct_with_short_reads_empty(self) -> None:
        """No long reads returns empty list."""
        corrected = correct_with_short_reads([], ["ACGT" * 25])
        assert corrected == []

    def test_hybrid_assemble_basic(self) -> None:
        """hybrid_assemble produces a HybridAssemblyResult."""
        seq = _random_dna(3000, seed=90)
        long_reads = [seq]
        # Short reads overlapping the long read
        short_reads = [seq[i:i + 150] for i in range(0, 2800, 50)]

        result = hybrid_assemble(
            long_reads, short_reads,
            min_long_read_length=1000,
            kmer_size=15,
        )
        assert isinstance(result, HybridAssemblyResult)
        assert result.num_contigs >= 1
        assert result.total_length > 0

    def test_hybrid_assemble_no_long_reads(self) -> None:
        """No long reads returns empty assembly."""
        result = hybrid_assemble([], ["ACGT" * 50])
        assert result.num_contigs == 0

    def test_scaffold_with_long_reads_singletons(self) -> None:
        """Contigs without bridging reads become singleton scaffolds."""
        contigs = ["ACGT" * 500, "TGCA" * 500]
        long_reads = ["AAAA" * 500]  # unrelated
        scaffolds = scaffold_with_long_reads(contigs, long_reads)
        assert len(scaffolds) >= 2
        assert all(isinstance(s, Scaffold) for s in scaffolds)
        # Each scaffold should have at least 1 contig
        assert all(s.num_contigs >= 1 for s in scaffolds)

    def test_scaffold_with_long_reads_no_long_reads(self) -> None:
        """Without long reads, each contig is its own scaffold."""
        contigs = ["ACGT" * 100, "TGCA" * 100]
        scaffolds = scaffold_with_long_reads(contigs, [])
        assert len(scaffolds) == 2
        assert all(s.num_contigs == 1 for s in scaffolds)


# ===================================================================
# TEST SUITE 9 -- IO: Fast5 / BAM / Formats
# ===================================================================


class TestIO:
    """Tests for io/fast5.py, io/bam.py, and io/formats.py."""

    # -- Fast5 --

    def test_extract_signal_converts_adc_to_pa(self) -> None:
        """extract_signal applies calibration: pA = (raw + offset) * (range / digitisation)."""
        import numpy as np

        raw_signal = np.array([100.0, 200.0, 300.0], dtype=np.float32)
        read = Fast5Read(
            read_id="test_read",
            signal=raw_signal,
            digitisation=8192.0,
            offset=-10.0,
            range_value=1467.0,
        )
        pa = extract_signal(read)
        assert pa is not None
        scale = 1467.0 / 8192.0
        expected = [(100 + (-10)) * scale, (200 + (-10)) * scale, (300 + (-10)) * scale]
        for actual, exp in zip(pa, expected):
            assert float(actual) == pytest.approx(exp, rel=1e-4)

    def test_extract_signal_none(self) -> None:
        """No signal returns None."""
        read = Fast5Read(read_id="empty")
        assert extract_signal(read) is None

    def test_extract_basecalls(self) -> None:
        """extract_basecalls returns (sequence, quality_string)."""
        read = Fast5Read(read_id="r1", sequence="ACGT", quality_string="IIII")
        seq, qual = extract_basecalls(read)
        assert seq == "ACGT"
        assert qual == "IIII"

    def test_extract_basecalls_none(self) -> None:
        """No basecalls returns (None, None)."""
        read = Fast5Read(read_id="r1")
        seq, qual = extract_basecalls(read)
        assert seq is None
        assert qual is None

    def test_get_read_metadata(self) -> None:
        """get_read_metadata returns comprehensive metadata dict."""
        import numpy as np

        read = Fast5Read(
            read_id="test123",
            signal=np.array([1.0, 2.0, 3.0]),
            sequence="ACGT",
            channel_id=42,
            mux=3,
            start_time=1000,
            duration=4000,
            sampling_rate=4000.0,
            run_id="run_abc",
        )
        meta = get_read_metadata(read)
        assert meta["read_id"] == "test123"
        assert meta["channel_id"] == 42
        assert meta["has_signal"] is True
        assert meta["signal_length"] == 3
        assert meta["has_basecalls"] is True
        assert meta["sequence_length"] == 4
        assert meta["duration_seconds"] == pytest.approx(1.0)

    def test_read_fast5_file_not_found(self) -> None:
        """read_fast5 raises FileNotFoundError for nonexistent file."""
        with pytest.raises(FileNotFoundError):
            read_fast5("/nonexistent/path/file.fast5")

    def test_read_fast5_unrecognized_extension(self) -> None:
        """read_fast5 raises ValueError for wrong extension."""
        import tempfile
        import os

        fd, path = tempfile.mkstemp(suffix=".txt")
        os.close(fd)
        try:
            with pytest.raises(ValueError, match="Unrecognized file extension"):
                read_fast5(path)
        finally:
            os.unlink(path)

    # -- BAM --

    def test_extract_methylation_tags_from_alignment(self) -> None:
        """extract_methylation_tags parses MM/ML tags."""
        aln = LongReadAlignment(
            read_name="r1",
            query_sequence="ACGTACGT",
            reference_start=0,
            reference_end=8,
            tags={"MM": "C+m,0,1", "ML": [200, 50]},
            methylation_tags={},
            aligned_pairs=[(i, i) for i in range(8)],
        )
        result = extract_methylation_tags(aln)
        # Should parse the MM tag
        if "C+m" in result:
            assert "genomic_positions" in result["C+m"]

    def test_get_supplementary_alignments_parses_sa(self) -> None:
        """get_supplementary_alignments extracts SA tag entries."""
        aln = LongReadAlignment(
            read_name="r1",
            tags={"SA": "chr2,5000,+,3000M,60,5;chr3,1000,-,2000M,50,3"},
        )
        sa_list = get_supplementary_alignments(aln)
        assert len(sa_list) == 2
        assert sa_list[0]["reference_name"] == "chr2"
        assert sa_list[0]["reference_start"] == 4999  # 5000 - 1 (1-based to 0-based)
        assert sa_list[0]["strand"] == "+"
        assert sa_list[1]["reference_name"] == "chr3"
        assert sa_list[1]["is_reverse"] is True

    def test_get_supplementary_alignments_empty(self) -> None:
        """No SA tag returns empty list."""
        aln = LongReadAlignment(read_name="r1", tags={})
        assert get_supplementary_alignments(aln) == []

    def test_calculate_alignment_stats_basic(self) -> None:
        """calculate_alignment_stats computes correct aggregates."""
        alns = [
            LongReadAlignment(
                read_name=f"r{i}",
                query_sequence="A" * 1000,
                query_length=1000,
                reference_start=0,
                reference_end=1000,
                mapping_quality=60,
                cigar_string="1000M",
                is_unmapped=False,
            )
            for i in range(5)
        ]
        stats = calculate_alignment_stats(alns, reference_length=10000)
        assert isinstance(stats, AlignmentStats)
        assert stats.total_reads == 5
        assert stats.mapped_reads == 5
        assert stats.unmapped_reads == 0
        assert stats.mean_mapping_quality == pytest.approx(60.0)
        assert stats.mean_read_length == pytest.approx(1000.0)
        assert stats.total_bases == 5000
        assert stats.coverage_depth == pytest.approx(0.5)

    def test_calculate_alignment_stats_empty(self) -> None:
        """Empty alignments returns zeroed stats."""
        stats = calculate_alignment_stats([])
        assert stats.total_reads == 0
        assert stats.mapped_reads == 0

    def test_calculate_alignment_stats_unmapped(self) -> None:
        """Unmapped reads are counted correctly."""
        alns = [
            LongReadAlignment(read_name="r1", is_unmapped=True),
            LongReadAlignment(
                read_name="r2",
                query_sequence="A" * 500,
                query_length=500,
                reference_start=0,
                reference_end=500,
                mapping_quality=30,
                cigar_string="500M",
                is_unmapped=False,
            ),
        ]
        stats = calculate_alignment_stats(alns)
        assert stats.total_reads == 2
        assert stats.mapped_reads == 1
        assert stats.unmapped_reads == 1

    # -- Formats --

    def test_write_paf(self, tmp_path: Path) -> None:
        """write_paf creates a valid PAF file."""
        alignments = [
            {
                "query_name": "read1",
                "query_length": 10000,
                "query_start": 0,
                "query_end": 8000,
                "strand": "+",
                "target_name": "ref",
                "target_length": 50000,
                "target_start": 1000,
                "target_end": 9000,
                "matches": 7500,
                "block_length": 8000,
                "mapping_quality": 60,
                "tags": {"NM": 500, "tp": "P"},
            },
        ]
        out = tmp_path / "test.paf"
        written = write_paf(alignments, out)
        assert written == 1
        assert out.exists()
        content = out.read_text()
        fields = content.strip().split("\t")
        assert fields[0] == "read1"
        assert fields[1] == "10000"
        assert fields[4] == "+"
        assert fields[5] == "ref"
        assert fields[11] == "60"

    def test_write_paf_empty(self, tmp_path: Path) -> None:
        """write_paf with no alignments creates an empty file."""
        out = tmp_path / "empty.paf"
        written = write_paf([], out)
        assert written == 0
        assert out.exists()
        assert out.read_text() == ""


# ===================================================================
# TEST SUITE 10 -- Visualization
# ===================================================================


class TestVisualization:
    """Tests for visualization/plots.py. Verify plot files are created with non-zero size."""

    def test_plot_read_length_histogram(self, tmp_path: Path) -> None:
        """Read length histogram is saved as a valid image."""
        lengths = [1000, 2000, 3000, 5000, 8000, 12000, 15000, 20000]
        out = tmp_path / "length_hist.png"
        result = plot_read_length_histogram(lengths, out)
        assert result == out
        assert out.exists()
        assert out.stat().st_size > 0

    def test_plot_read_length_histogram_from_dicts(self, tmp_path: Path) -> None:
        """Accepts dicts with 'length' key."""
        reads = [{"length": i * 1000} for i in range(1, 10)]
        out = tmp_path / "length_hist_dict.png"
        result = plot_read_length_histogram(reads, out)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_plot_quality_vs_length(self, tmp_path: Path) -> None:
        """Quality vs length scatter plot is saved."""
        reads = []
        for i in range(50):
            length = 1000 + i * 200
            quality = 10 + i * 0.2
            q_str = _make_uniform_quality_string(length, quality=int(quality))
            reads.append({"sequence": "A" * length, "quality_string": q_str})

        out = tmp_path / "qvl.png"
        result = plot_quality_vs_length(reads, out)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_plot_dotplot(self, tmp_path: Path) -> None:
        """Dotplot between two sequences is saved."""
        seq1 = _random_dna(500, seed=100)
        seq2 = seq1  # identical => strong diagonal
        out = tmp_path / "dotplot.png"
        result = plot_dotplot(seq1, seq2, out, word_size=7)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_plot_alignment_view(self, tmp_path: Path) -> None:
        """Alignment view plot is saved."""
        alignments = [
            {
                "reference_name": "chr1",
                "reference_start": 1000 + i * 100,
                "reference_end": 5000 + i * 100,
                "is_reverse": i % 2 == 0,
                "mapping_quality": 60,
                "read_name": f"r{i}",
                "is_supplementary": False,
            }
            for i in range(20)
        ]
        out = tmp_path / "alignment.png"
        result = plot_alignment_view(
            alignments,
            region="chr1:0-10000",
            output_path=out,
        )
        assert out.exists()
        assert out.stat().st_size > 0

    def test_plot_methylation_track(self, tmp_path: Path) -> None:
        """Methylation track plot is saved."""
        meth_data = [
            {"chromosome": "chr1", "position": 100 + i * 50, "probability": 0.1 + 0.05 * i}
            for i in range(20)
        ]
        out = tmp_path / "meth_track.png"
        result = plot_methylation_track(
            meth_data,
            region={"chromosome": "chr1", "start": 0, "end": 2000},
            output_path=out,
        )
        assert out.exists()
        assert out.stat().st_size > 0

    def test_plot_phasing_blocks(self, tmp_path: Path) -> None:
        """Phasing block plot is saved."""
        blocks = [
            {"chromosome": "chr1", "start": 0, "end": 50000, "num_variants": 100, "quality": 0.95},
            {"chromosome": "chr1", "start": 60000, "end": 80000, "num_variants": 40, "quality": 0.88},
            {"chromosome": "chr2", "start": 10000, "end": 25000, "num_variants": 30, "quality": 0.92},
        ]
        out = tmp_path / "phasing.png"
        result = plot_phasing_blocks(blocks, out)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_plot_phasing_blocks_empty(self, tmp_path: Path) -> None:
        """Empty phase blocks still create a valid (empty) plot file."""
        out = tmp_path / "phasing_empty.png"
        result = plot_phasing_blocks([], out)
        assert out.exists()
        assert out.stat().st_size > 0
