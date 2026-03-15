"""Tests for epigenome assay-specific functions.

Tests methylation (compute_beta_values, load_cpg_table,
summarize_beta_by_chromosome, MethylationSite, calculate_methylation_statistics),
ChIP-seq (ChIPPeak, calculate_peak_statistics, filter_peaks_by_score,
find_overlapping_peaks, merge_overlapping_peaks), and
ATAC-seq (ATACPeak, calculate_atac_statistics, calculate_atac_specific_metrics,
compare_atac_conditions) using real implementations. NO MOCKING.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from metainformant.epigenome.assays.atacseq import (
    ATACPeak,
    calculate_atac_specific_metrics,
    calculate_atac_statistics,
    compare_atac_conditions,
    identify_tss_enrichment,
)
from metainformant.epigenome.assays.chipseq import (
    ChIPPeak,
    calculate_peak_enrichment,
    calculate_peak_statistics,
    filter_peaks_by_score,
    find_overlapping_peaks,
    merge_overlapping_peaks,
)
from metainformant.epigenome.assays.methylation import (
    MethylationSite,
    calculate_methylation_statistics,
    compute_beta_values,
    load_cpg_table,
    summarize_beta_by_chromosome,
)

# ---------------------------------------------------------------------------
# MethylationSite
# ---------------------------------------------------------------------------


class TestMethylationSite:
    """Tests for MethylationSite data class."""

    def test_basic_creation(self) -> None:
        site = MethylationSite("chr1", 1000, 8, 10)
        assert site.chromosome == "chr1"
        assert site.position == 1000
        assert site.methylated_reads == 8
        assert site.total_reads == 10

    def test_methylation_level(self) -> None:
        site = MethylationSite("chr1", 100, 7, 10)
        assert site.methylation_level == pytest.approx(0.7)

    def test_zero_coverage_methylation_level(self) -> None:
        # Create with 0 methylated and 0 total - needs total >= methylated >= 0
        site = MethylationSite("chr1", 100, 0, 0)
        assert site.methylation_level == 0.0

    def test_coverage_property(self) -> None:
        site = MethylationSite("chr1", 100, 3, 15)
        assert site.coverage == 15

    def test_to_dict(self) -> None:
        site = MethylationSite("chr2", 500, 5, 10, strand="-")
        d = site.to_dict()
        assert d["chromosome"] == "chr2"
        assert d["position"] == 500
        assert d["methylation_level"] == pytest.approx(0.5)
        assert d["strand"] == "-"

    def test_strand_default(self) -> None:
        site = MethylationSite("chr1", 100, 5, 10)
        assert site.strand == "+"

    def test_full_methylation(self) -> None:
        site = MethylationSite("chr1", 100, 10, 10)
        assert site.methylation_level == pytest.approx(1.0)

    def test_no_methylation(self) -> None:
        site = MethylationSite("chr1", 100, 0, 10)
        assert site.methylation_level == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# calculate_methylation_statistics
# ---------------------------------------------------------------------------


class TestCalculateMethylationStatistics:
    """Tests for calculate_methylation_statistics."""

    def test_basic_statistics(self) -> None:
        data = {
            "chr1": [
                MethylationSite("chr1", 100, 8, 10),
                MethylationSite("chr1", 200, 3, 10),
            ],
            "chr2": [
                MethylationSite("chr2", 100, 5, 10),
            ],
        }
        stats = calculate_methylation_statistics(data)
        assert stats["total_sites"] == 3
        assert stats["total_chromosomes"] == 2
        assert 0.0 <= stats["mean_methylation"] <= 1.0

    def test_empty_data_returns_empty(self) -> None:
        assert calculate_methylation_statistics({}) == {}

    def test_chromosome_stats_present(self) -> None:
        data = {
            "chr1": [MethylationSite("chr1", 100, 8, 10)],
            "chr2": [MethylationSite("chr2", 200, 2, 10)],
        }
        stats = calculate_methylation_statistics(data)
        assert "chromosome_stats" in stats
        assert "chr1" in stats["chromosome_stats"]
        assert "chr2" in stats["chromosome_stats"]

    def test_coverage_statistics(self) -> None:
        data = {
            "chr1": [
                MethylationSite("chr1", 100, 5, 20),
                MethylationSite("chr1", 200, 3, 10),
            ],
        }
        stats = calculate_methylation_statistics(data)
        assert stats["mean_coverage"] == pytest.approx(15.0)
        assert stats["min_coverage"] == 10
        assert stats["max_coverage"] == 20


# ---------------------------------------------------------------------------
# compute_beta_values and load_cpg_table
# ---------------------------------------------------------------------------


class TestComputeBetaValues:
    """Tests for compute_beta_values."""

    def test_basic_beta_computation(self) -> None:
        df = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1", "chr2"],
                "pos": [100, 200, 300],
                "methylated": [8, 0, 5],
                "unmethylated": [2, 10, 5],
            }
        )
        result = compute_beta_values(df)
        assert "beta" in result.columns
        assert result["beta"].iloc[0] == pytest.approx(0.8)
        assert result["beta"].iloc[1] == pytest.approx(0.0)
        assert result["beta"].iloc[2] == pytest.approx(0.5)

    def test_handles_zero_total(self) -> None:
        df = pd.DataFrame(
            {
                "methylated": [0],
                "unmethylated": [0],
            }
        )
        result = compute_beta_values(df)
        # Should handle division by zero gracefully
        assert result["beta"].iloc[0] == pytest.approx(0.0)

    def test_missing_columns_raises(self) -> None:
        df = pd.DataFrame({"foo": [1], "bar": [2]})
        with pytest.raises(ValueError, match="must contain"):
            compute_beta_values(df)

    def test_does_not_modify_original(self) -> None:
        df = pd.DataFrame(
            {
                "methylated": [5],
                "unmethylated": [5],
            }
        )
        result = compute_beta_values(df)
        assert "beta" not in df.columns  # Original unchanged
        assert "beta" in result.columns

    def test_load_cpg_table(self) -> None:
        repo_root = Path(__file__).resolve().parents[1]
        table_path = repo_root / "tests/data/epigenome/cpg_counts.tsv"
        df = load_cpg_table(table_path)
        assert "chrom" in df.columns
        assert "pos" in df.columns
        assert "methylated" in df.columns
        assert "unmethylated" in df.columns


# ---------------------------------------------------------------------------
# summarize_beta_by_chromosome
# ---------------------------------------------------------------------------


class TestSummarizeBetaByChromosome:
    """Tests for summarize_beta_by_chromosome."""

    def test_basic_summary(self) -> None:
        df = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1", "chr2", "chr2"],
                "beta": [0.8, 0.6, 0.3, 0.5],
            }
        )
        summary = summarize_beta_by_chromosome(df)
        assert "chr1" in summary.index
        assert "chr2" in summary.index
        assert summary.loc["chr1", "mean"] == pytest.approx(0.7)
        assert summary.loc["chr2", "mean"] == pytest.approx(0.4)

    def test_missing_columns_raises(self) -> None:
        df = pd.DataFrame({"foo": [1]})
        with pytest.raises(ValueError, match="must contain"):
            summarize_beta_by_chromosome(df)

    def test_has_expected_agg_columns(self) -> None:
        df = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1", "chr1"],
                "beta": [0.2, 0.5, 0.8],
            }
        )
        summary = summarize_beta_by_chromosome(df)
        assert "mean" in summary.columns
        assert "count" in summary.columns
        assert "min" in summary.columns
        assert "max" in summary.columns


# ---------------------------------------------------------------------------
# ChIPPeak
# ---------------------------------------------------------------------------


class TestChIPPeak:
    """Tests for ChIPPeak data class."""

    def test_basic_creation(self) -> None:
        peak = ChIPPeak("chr1", 1000, 2000, summit=1500, score=100.0)
        assert peak.chromosome == "chr1"
        assert peak.length == 1000
        assert peak.center == 1500

    def test_overlap_detection(self) -> None:
        peak1 = ChIPPeak("chr1", 100, 200)
        peak2 = ChIPPeak("chr1", 150, 250)
        peak3 = ChIPPeak("chr1", 300, 400)
        peak4 = ChIPPeak("chr2", 100, 200)

        assert peak1.overlaps_with(peak2)
        assert not peak1.overlaps_with(peak3)
        assert not peak1.overlaps_with(peak4)  # Different chromosome

    def test_to_dict(self) -> None:
        peak = ChIPPeak("chr1", 100, 500, summit=300, score=50.0, signal_value=25.0, p_value=0.001, q_value=0.01)
        d = peak.to_dict()
        assert d["chromosome"] == "chr1"
        assert d["length"] == 400
        assert d["p_value"] == 0.001

    def test_to_bed_format(self) -> None:
        peak = ChIPPeak("chr1", 100, 500, summit=300, score=50.0)
        bed = peak.to_bed_format()
        assert bed.startswith("chr1\t100\t500")


# ---------------------------------------------------------------------------
# calculate_peak_statistics
# ---------------------------------------------------------------------------


class TestCalculatePeakStatistics:
    """Tests for calculate_peak_statistics."""

    def test_basic_stats(self) -> None:
        peaks = [
            ChIPPeak("chr1", 100, 300, score=50.0, signal_value=10.0),
            ChIPPeak("chr1", 500, 800, score=100.0, signal_value=20.0),
            ChIPPeak("chr2", 200, 400, score=75.0, signal_value=15.0),
        ]
        stats = calculate_peak_statistics(peaks)
        assert stats["total_peaks"] == 3
        assert stats["mean_length"] == pytest.approx(233.33, abs=1.0)
        assert stats["mean_score"] == 75.0

    def test_empty_peaks(self) -> None:
        assert calculate_peak_statistics([]) == {}

    def test_chromosome_distribution(self) -> None:
        peaks = [
            ChIPPeak("chr1", 100, 200, score=1.0),
            ChIPPeak("chr1", 300, 400, score=1.0),
            ChIPPeak("chr2", 100, 200, score=1.0),
        ]
        stats = calculate_peak_statistics(peaks)
        assert stats["chromosome_distribution"]["chr1"] == 2
        assert stats["chromosome_distribution"]["chr2"] == 1


# ---------------------------------------------------------------------------
# filter_peaks_by_score
# ---------------------------------------------------------------------------


class TestFilterPeaksByScore:
    """Tests for filter_peaks_by_score."""

    def test_basic_filtering(self) -> None:
        peaks = [
            ChIPPeak("chr1", 100, 200, score=50.0),
            ChIPPeak("chr1", 300, 400, score=150.0),
            ChIPPeak("chr1", 500, 600, score=25.0),
        ]
        filtered = filter_peaks_by_score(peaks, min_score=100.0)
        assert len(filtered) == 1
        assert filtered[0].score == 150.0

    def test_max_peaks_limit(self) -> None:
        peaks = [ChIPPeak("chr1", i * 200, i * 200 + 100, score=float(i)) for i in range(1, 11)]
        filtered = filter_peaks_by_score(peaks, min_score=0, max_peaks=3)
        assert len(filtered) == 3
        # Should be sorted by score descending
        assert filtered[0].score >= filtered[1].score


# ---------------------------------------------------------------------------
# find_overlapping_peaks
# ---------------------------------------------------------------------------


class TestFindOverlappingPeaks:
    """Tests for find_overlapping_peaks."""

    def test_finds_overlaps(self) -> None:
        peaks1 = [ChIPPeak("chr1", 100, 300)]
        peaks2 = [ChIPPeak("chr1", 200, 400)]
        overlaps = find_overlapping_peaks(peaks1, peaks2)
        assert len(overlaps) == 1

    def test_no_overlap(self) -> None:
        peaks1 = [ChIPPeak("chr1", 100, 200)]
        peaks2 = [ChIPPeak("chr1", 300, 400)]
        overlaps = find_overlapping_peaks(peaks1, peaks2)
        assert len(overlaps) == 0

    def test_different_chromosomes(self) -> None:
        peaks1 = [ChIPPeak("chr1", 100, 200)]
        peaks2 = [ChIPPeak("chr2", 100, 200)]
        overlaps = find_overlapping_peaks(peaks1, peaks2)
        assert len(overlaps) == 0


# ---------------------------------------------------------------------------
# merge_overlapping_peaks
# ---------------------------------------------------------------------------


class TestMergeOverlappingPeaks:
    """Tests for merge_overlapping_peaks (ChIPPeak version)."""

    def test_merges_overlapping(self) -> None:
        peaks = [
            ChIPPeak("chr1", 100, 300, score=50.0, signal_value=10.0),
            ChIPPeak("chr1", 200, 400, score=75.0, signal_value=15.0),
        ]
        merged = merge_overlapping_peaks(peaks)
        assert len(merged) == 1
        assert merged[0].start == 100
        assert merged[0].end == 400

    def test_no_merge_distant_peaks(self) -> None:
        peaks = [
            ChIPPeak("chr1", 100, 200, score=50.0),
            ChIPPeak("chr1", 500, 600, score=50.0),
        ]
        merged = merge_overlapping_peaks(peaks)
        assert len(merged) == 2

    def test_empty_input(self) -> None:
        assert merge_overlapping_peaks([]) == []

    def test_merge_with_distance(self) -> None:
        peaks = [
            ChIPPeak("chr1", 100, 200, score=50.0, signal_value=10.0),
            ChIPPeak("chr1", 210, 300, score=75.0, signal_value=15.0),
        ]
        merged = merge_overlapping_peaks(peaks, max_distance=20)
        assert len(merged) == 1


# ---------------------------------------------------------------------------
# ATACPeak
# ---------------------------------------------------------------------------


class TestATACPeak:
    """Tests for ATACPeak data class."""

    def test_basic_creation(self) -> None:
        peak = ATACPeak("chr1", 1000, 2000, score=50.0, signal_value=25.0)
        assert peak.chromosome == "chr1"
        assert peak.length == 1000
        assert peak.accessibility_score == 25.0

    def test_overlap_detection(self) -> None:
        peak1 = ATACPeak("chr1", 100, 200)
        peak2 = ATACPeak("chr1", 150, 250)
        assert peak1.overlaps_with(peak2)

    def test_to_dict(self) -> None:
        peak = ATACPeak("chr1", 100, 500, score=50.0, signal_value=30.0)
        d = peak.to_dict()
        assert d["accessibility_score"] == 30.0
        assert d["length"] == 400


# ---------------------------------------------------------------------------
# calculate_atac_statistics
# ---------------------------------------------------------------------------


class TestCalculateAtacStatistics:
    """Tests for calculate_atac_statistics."""

    def test_basic_stats(self) -> None:
        peaks = [
            ATACPeak("chr1", 100, 200, score=50.0, signal_value=10.0),
            ATACPeak("chr1", 500, 700, score=100.0, signal_value=20.0),
            ATACPeak("chr2", 200, 350, score=75.0, signal_value=15.0),
        ]
        stats = calculate_atac_statistics(peaks)
        assert stats["total_peaks"] == 3
        assert "mean_length" in stats
        assert "mean_score" in stats

    def test_empty_peaks(self) -> None:
        assert calculate_atac_statistics([]) == {}


# ---------------------------------------------------------------------------
# calculate_atac_specific_metrics
# ---------------------------------------------------------------------------


class TestCalculateAtacSpecificMetrics:
    """Tests for calculate_atac_specific_metrics."""

    def test_nucleosome_fractions(self) -> None:
        peaks = [
            ATACPeak("chr1", 0, 100, signal_value=10.0),  # NFR
            ATACPeak("chr1", 200, 400, signal_value=10.0),  # Mononucleosome
            ATACPeak("chr1", 500, 800, signal_value=10.0),  # Dinucleosome
        ]
        metrics = calculate_atac_specific_metrics(peaks)
        assert "nfr_peak_fraction" in metrics
        assert 0.0 <= metrics["nfr_peak_fraction"] <= 1.0


# ---------------------------------------------------------------------------
# compare_atac_conditions
# ---------------------------------------------------------------------------


class TestCompareAtacConditions:
    """Tests for compare_atac_conditions."""

    def test_identical_conditions(self) -> None:
        peaks = [
            ATACPeak("chr1", 100, 200, signal_value=10.0),
            ATACPeak("chr1", 300, 400, signal_value=15.0),
        ]
        result = compare_atac_conditions(peaks, peaks)
        assert result["condition1_total"] == 2
        assert result["condition2_total"] == 2
        assert result["overlapping_peaks"] >= 1

    def test_no_overlap_conditions(self) -> None:
        c1_peaks = [ATACPeak("chr1", 100, 200, signal_value=10.0)]
        c2_peaks = [ATACPeak("chr1", 500, 600, signal_value=10.0)]
        result = compare_atac_conditions(c1_peaks, c2_peaks)
        assert result["overlapping_peaks"] == 0
        assert result["condition1_only"] == 1

    def test_empty_conditions(self) -> None:
        result = compare_atac_conditions([], [])
        assert result["condition1_total"] == 0
        assert result["condition2_total"] == 0


# ---------------------------------------------------------------------------
# identify_tss_enrichment
# ---------------------------------------------------------------------------


class TestIdentifyTssEnrichment:
    """Tests for identify_tss_enrichment."""

    def test_enrichment_with_overlapping_peaks(self) -> None:
        peaks = [
            ATACPeak("chr1", 900, 1100, signal_value=10.0),
            ATACPeak("chr1", 1900, 2100, signal_value=10.0),
        ]
        tss_positions = {"chr1": [1000, 2000]}
        result = identify_tss_enrichment(peaks, tss_positions, window_size=500)
        assert result["total_tss"] == 2
        assert result["enriched_tss"] == 2

    def test_no_enrichment(self) -> None:
        peaks = [ATACPeak("chr1", 5000, 5100, signal_value=10.0)]
        tss_positions = {"chr1": [100]}
        result = identify_tss_enrichment(peaks, tss_positions, window_size=200)
        assert result["enriched_tss"] == 0


# ---------------------------------------------------------------------------
# calculate_peak_enrichment
# ---------------------------------------------------------------------------


class TestCalculatePeakEnrichment:
    """Tests for calculate_peak_enrichment (ChIPPeak)."""

    def test_basic_enrichment(self) -> None:
        peaks = [
            ChIPPeak("chr1", 100, 300, score=50.0),
            ChIPPeak("chr1", 500, 700, score=75.0),
        ]
        result = calculate_peak_enrichment(peaks, genome_size=10000)
        assert result["total_peaks"] == 2
        assert result["total_peak_bases"] == 400
        assert 0.0 <= result["genome_coverage"] <= 1.0

    def test_empty_peaks(self) -> None:
        assert calculate_peak_enrichment([], genome_size=10000) == {}
