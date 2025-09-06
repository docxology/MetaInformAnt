"""Tests for FASTQ quality control module.

Real implementation testing for FASTQ quality assessment functions.
No mocking used - all tests use real data and computational methods.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pytest

from metainformant.dna.fastq import FastqRecord
from metainformant.quality.fastq import (
    adapter_content,
    analyze_fastq_quality,
    basic_statistics,
    duplication_levels,
    gc_content_distribution,
    n_content_per_position,
    overrepresented_sequences,
    per_base_quality,
    per_sequence_quality,
    quality_score_distribution,
    sequence_length_distribution,
)


class TestBasicStatistics:
    """Test basic FASTQ statistics calculation."""

    def setup_method(self):
        """Setup test FASTQ records with known properties."""
        self.test_reads = [
            FastqRecord("read1", "ATCGATCG", "IIIIIIII"),  # Length 8, high quality
            FastqRecord("read2", "GCGCGCGC", "HHHHHHHH"),  # Length 8, medium quality
            FastqRecord("read3", "ATGCATGCAT", "IIIIIIIIIII"),  # Length 10, high quality
            FastqRecord("read4", "NNNNNNNN", "!!!!!!!!"),  # Length 8, low quality, all N
        ]

    def test_basic_statistics_calculation(self):
        """Test calculation of basic FASTQ statistics."""
        stats = basic_statistics(self.test_reads)

        # Check basic counts and lengths
        assert stats["total_sequences"] == 4
        assert stats["sequence_length_min"] == 8
        assert stats["sequence_length_max"] == 10
        assert stats["sequence_length_mean"] == 8.5  # (8+8+10+8)/4
        assert stats["sequence_length_median"] == 8.0

        # Check GC content calculation
        # read1: ATCGATCG -> A=2, T=2, C=2, G=2 -> GC = 4/8 = 50%
        # read2: GCGCGCGC -> G=4, C=4 -> GC = 8/8 = 100%
        # read3: ATGCATGCAT -> A=3, T=3, G=2, C=2 -> GC = 4/10 = 40%
        # read4: NNNNNNNN -> All N, no ATCG -> GC = 0/0 -> 0%
        expected_gc = (50 + 100 + 40 + 0) / 4  # 47.5%
        assert abs(stats["gc_content_mean"] - expected_gc) < 0.1

        # Check quality scores
        assert "quality_score_mean" in stats
        assert "quality_score_median" in stats
        assert stats["quality_score_min"] >= 0
        assert stats["quality_score_max"] <= 41  # Phred scale

        # Check uniformity
        assert stats["uniform_length"] is False  # Different lengths

    def test_basic_statistics_uniform_length(self):
        """Test statistics with uniform-length reads."""
        uniform_reads = [
            FastqRecord("read1", "ATCG", "IIII"),
            FastqRecord("read2", "GCTA", "HHHH"),
            FastqRecord("read3", "TTTT", "JJJJ"),
        ]

        stats = basic_statistics(uniform_reads)
        assert stats["uniform_length"] is True
        assert stats["sequence_length_min"] == 4
        assert stats["sequence_length_max"] == 4
        assert stats["sequence_length_mean"] == 4.0

    def test_basic_statistics_empty_input(self):
        """Test handling of empty input."""
        stats = basic_statistics([])
        assert stats == {}


class TestPerBaseQuality:
    """Test per-base quality analysis."""

    def test_per_base_quality_calculation(self):
        """Test per-base quality score calculation."""
        reads = [
            FastqRecord("read1", "ATCG", "IIII"),  # Quality 40,40,40,40
            FastqRecord("read2", "GCTA", "HHHH"),  # Quality 39,39,39,39
            FastqRecord("read3", "TTGG", "JJJJ"),  # Quality 41,41,41,41
        ]

        result = per_base_quality(reads)

        assert "per_base_quality" in result
        assert "max_position" in result
        assert result["max_position"] == 4

        per_base_stats = result["per_base_quality"]
        assert len(per_base_stats) == 4

        # Check first position statistics
        pos1_stats = per_base_stats[0]
        assert pos1_stats["position"] == 1
        assert pos1_stats["count"] == 3  # All 3 reads have this position

        # Quality scores: I=40, H=39, J=41, so mean should be 40
        expected_mean = (40 + 39 + 41) / 3
        assert abs(pos1_stats["mean"] - expected_mean) < 0.1

    def test_per_base_quality_variable_lengths(self):
        """Test per-base quality with variable read lengths."""
        reads = [
            FastqRecord("read1", "AT", "II"),  # Length 2
            FastqRecord("read2", "GCT", "HHH"),  # Length 3
            FastqRecord("read3", "ATGC", "JJJJ"),  # Length 4
        ]

        result = per_base_quality(reads)

        assert result["max_position"] == 4
        per_base_stats = result["per_base_quality"]

        # Position 1: all 3 reads
        assert per_base_stats[0]["count"] == 3
        # Position 2: all 3 reads
        assert per_base_stats[1]["count"] == 3
        # Position 3: 2 reads (read2, read3)
        assert per_base_stats[2]["count"] == 2
        # Position 4: 1 read (read3)
        assert per_base_stats[3]["count"] == 1


class TestSequenceLengthDistribution:
    """Test sequence length distribution analysis."""

    def test_length_distribution_calculation(self):
        """Test sequence length distribution calculation."""
        reads = [
            FastqRecord("read1", "ATCG", "IIII"),  # Length 4
            FastqRecord("read2", "ATCGATCG", "IIIIIIII"),  # Length 8
            FastqRecord("read3", "ATCG", "IIII"),  # Length 4 (duplicate)
            FastqRecord("read4", "ATCGATCGATCG", "IIIIIIIIIIII"),  # Length 12
        ]

        result = sequence_length_distribution(reads)

        assert result["min_length"] == 4
        assert result["max_length"] == 12
        assert result["mean_length"] == 7.0  # (4+8+4+12)/4
        assert result["median_length"] == 6.0  # Between 4 and 8
        assert result["mode_length"] == 4  # Most frequent length
        assert result["uniform_length"] is False

        # Check distribution
        length_dist = result["length_distribution"]
        assert length_dist[4] == 2  # Two reads of length 4
        assert length_dist[8] == 1  # One read of length 8
        assert length_dist[12] == 1  # One read of length 12


class TestGCContentDistribution:
    """Test GC content distribution analysis."""

    def test_gc_content_calculation(self):
        """Test GC content calculation."""
        reads = [
            FastqRecord("read1", "ATCG", "IIII"),  # 50% GC (2/4)
            FastqRecord("read2", "AAAA", "IIII"),  # 0% GC (0/4)
            FastqRecord("read3", "GCGC", "IIII"),  # 100% GC (4/4)
            FastqRecord("read4", "ATNNNN", "IIIIII"),  # 0% GC (0/2 excluding N)
        ]

        result = gc_content_distribution(reads)

        # Expected GC contents: 50%, 0%, 100%, 0%
        expected_mean = (50 + 0 + 100 + 0) / 4  # 37.5%
        assert abs(result["mean_gc"] - expected_mean) < 0.1

        assert "gc_contents" in result
        assert "histogram_counts" in result
        assert "gc_distribution" in result

        gc_contents = result["gc_contents"]
        assert len(gc_contents) == 4
        assert 50.0 in gc_contents
        assert 0.0 in gc_contents
        assert 100.0 in gc_contents


class TestAdapterContent:
    """Test adapter contamination detection."""

    def test_adapter_detection_basic(self):
        """Test basic adapter detection."""
        reads = [
            FastqRecord("read1", "ATCGATCGATCG", "IIIIIIIIIIII"),  # No adapter
            FastqRecord("read2", "AGATCGGAAGAGATCG", "IIIIIIIIIIIIIIII"),  # Illumina adapter
            FastqRecord("read3", "ATCGAGATCGGAAGAG", "IIIIIIIIIIIIIIII"),  # Adapter at end
        ]

        result = adapter_content(reads)

        assert "Illumina_Universal" in result
        illumina_result = result["Illumina_Universal"]

        # Should detect adapter in read2 and read3
        assert illumina_result["total_contaminated_reads"] >= 1
        assert illumina_result["contamination_percentage"] > 0
        assert "position_contamination" in illumina_result

    def test_custom_adapters(self):
        """Test custom adapter sequences."""
        reads = [
            FastqRecord("read1", "ATCGCUSTOMADAPTER", "IIIIIIIIIIIIIIIII"),
        ]

        custom_adapters = {"Custom": "CUSTOMADAPTER"}
        result = adapter_content(reads, adapters=custom_adapters)

        assert "Custom" in result
        assert result["Custom"]["total_contaminated_reads"] == 1
        assert result["Custom"]["contamination_percentage"] == 100.0


class TestOverrepresentedSequences:
    """Test overrepresented sequence detection."""

    def test_overrepresented_detection(self):
        """Test detection of overrepresented sequences."""
        # Create reads with some duplicates
        reads = []

        # Add 5 copies of the same sequence (overrepresented)
        for i in range(5):
            reads.append(FastqRecord(f"dup_{i}", "ATCGATCGATCG", "IIIIIIIIIIII"))

        # Add 10 unique sequences
        for i in range(10):
            unique_seq = "A" * (8 + i % 4)  # Varying lengths
            reads.append(FastqRecord(f"unique_{i}", unique_seq, "I" * len(unique_seq)))

        result = overrepresented_sequences(reads, min_count=3, min_percentage=10.0)

        assert "overrepresented_sequences" in result
        overrep_seqs = result["overrepresented_sequences"]

        # Should detect the duplicated sequence
        assert len(overrep_seqs) >= 1

        # Check that ATCGATCGATCG is detected as overrepresented
        found_overrep = False
        for seq, data in overrep_seqs.items():
            if seq == "ATCGATCGATCG":
                found_overrep = True
                assert data["count"] == 5
                assert abs(data["percentage"] - (5 / 15 * 100)) < 0.1  # 33.33%
                break

        assert found_overrep, "Expected overrepresented sequence not found"

    def test_sequence_source_identification(self):
        """Test identification of sequence sources."""
        from metainformant.quality.fastq import _identify_sequence_source

        # Test various sequence types
        assert _identify_sequence_source("AAAAAAAAAA") == "Poly-A tail"
        assert _identify_sequence_source("TTTTTTTTTT") == "Poly-T sequence"
        assert _identify_sequence_source("AGATCGGAAGAG") == "Illumina adapter"
        assert _identify_sequence_source("ATCGATCGATCG") == "Unknown"


class TestDuplicationLevels:
    """Test sequence duplication analysis."""

    def test_duplication_analysis(self):
        """Test sequence duplication level analysis."""
        reads = [
            # Unique sequences (appear once)
            FastqRecord("unique1", "ATCG", "IIII"),
            FastqRecord("unique2", "GCTA", "IIII"),
            # Duplicated sequence (appears twice)
            FastqRecord("dup1", "TTTT", "IIII"),
            FastqRecord("dup2", "TTTT", "IIII"),
            # Highly duplicated (appears 3 times)
            FastqRecord("high1", "AAAA", "IIII"),
            FastqRecord("high2", "AAAA", "IIII"),
            FastqRecord("high3", "AAAA", "IIII"),
        ]

        result = duplication_levels(reads)

        assert result["total_sequences"] == 7
        assert result["unique_sequences"] == 4  # 4 distinct sequences

        # Check duplication levels
        dup_levels = result["duplication_levels"]

        # Should have: 2 sequences appearing 1x, 1 sequence appearing 2x, 1 sequence appearing 3x
        assert 1 in dup_levels  # Unique sequences
        assert 2 in dup_levels  # Duplicated sequences
        assert 3 in dup_levels  # Highly duplicated sequences

        assert dup_levels[1]["sequences_with_this_dup_level"] == 2  # ATCG, GCTA
        assert dup_levels[2]["sequences_with_this_dup_level"] == 1  # TTTT
        assert dup_levels[3]["sequences_with_this_dup_level"] == 1  # AAAA


class TestQualityScoreDistribution:
    """Test quality score distribution analysis."""

    def test_quality_distribution_calculation(self):
        """Test quality score distribution calculation."""
        reads = [
            FastqRecord("read1", "ATCG", "!!!!"),  # Quality 0,0,0,0
            FastqRecord("read2", "GCTA", "IIII"),  # Quality 40,40,40,40
            FastqRecord("read3", "TTGG", "@@@@"),  # Quality 31,31,31,31
        ]

        result = quality_score_distribution(reads)

        assert "quality_scores" in result
        assert "histogram_counts" in result
        assert "quality_distribution" in result

        quality_scores = result["quality_scores"]
        assert len(quality_scores) == 12  # 4 bases * 3 reads

        # Check that we have the expected quality scores
        assert quality_scores.count(0) == 4  # From read1
        assert quality_scores.count(40) == 4  # From read2
        assert quality_scores.count(31) == 4  # From read3

        expected_mean = (0 * 4 + 40 * 4 + 31 * 4) / 12  # 23.67
        assert abs(result["mean_quality"] - expected_mean) < 0.1


class TestIntegratedAnalysis:
    """Test integrated FASTQ quality analysis."""

    def test_analyze_fastq_quality_with_file(self, tmp_path):
        """Test complete FASTQ analysis with actual file."""
        # Create test FASTQ file
        fastq_content = """@read1
ATCGATCGATCGATCG
+
IIIIIIIIIIIIIIII
@read2
GCGCGCGCGCGCGCGC
+
HHHHHHHHHHHHHHHH
@read3
ATGCATGCATGCATGC
+
JJJJJJJJJJJJJJJJ
@read4
TTTTTTTTTTTTTTTT
+
!!!!!!!!!!!!!!!!
"""

        fastq_path = tmp_path / "test.fastq"
        fastq_path.write_text(fastq_content)

        # Run complete analysis
        results = analyze_fastq_quality(fastq_path)

        # Check that all expected analyses are present
        expected_analyses = [
            "basic_stats",
            "per_base_quality",
            "per_sequence_quality",
            "sequence_lengths",
            "gc_content",
            "adapter_content",
            "overrepresented_seqs",
            "duplication_levels",
            "n_content",
            "quality_scores",
            "metadata",
        ]

        for analysis in expected_analyses:
            assert analysis in results, f"Missing analysis: {analysis}"

        # Check metadata
        metadata = results["metadata"]
        assert metadata["filename"] == "test.fastq"
        assert metadata["reads_analyzed"] == 4
        assert metadata["subsampled"] is False

        # Basic validation of results
        basic_stats = results["basic_stats"]
        assert basic_stats["total_sequences"] == 4
        assert basic_stats["uniform_length"] is True  # All length 16

    def test_analyze_fastq_quality_with_subsampling(self, tmp_path):
        """Test FASTQ analysis with subsampling."""
        # Create FASTQ with more reads than subsample limit
        fastq_content = ""
        for i in range(10):
            fastq_content += f"""@read{i}
ATCGATCGATCGATCG
+
IIIIIIIIIIIIIIII
"""

        fastq_path = tmp_path / "test_large.fastq"
        fastq_path.write_text(fastq_content)

        # Run analysis with subsampling
        results = analyze_fastq_quality(fastq_path, subsample=5, seed=42)

        metadata = results["metadata"]
        assert metadata["total_reads_in_file"] == 10
        assert metadata["reads_analyzed"] == 5
        assert metadata["subsampled"] is True


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_fastq_file(self, tmp_path):
        """Test handling of empty FASTQ file."""
        empty_path = tmp_path / "empty.fastq"
        empty_path.write_text("")

        with pytest.raises(ValueError, match="No reads found"):
            analyze_fastq_quality(empty_path)

    def test_nonexistent_file(self):
        """Test handling of nonexistent file."""
        with pytest.raises(FileNotFoundError):
            analyze_fastq_quality("nonexistent.fastq")

    def test_single_read_analysis(self):
        """Test analysis with single read."""
        reads = [FastqRecord("single", "ATCGATCG", "IIIIIIII")]

        stats = basic_statistics(reads)
        assert stats["total_sequences"] == 1
        assert stats["uniform_length"] is True

        length_dist = sequence_length_distribution(reads)
        assert length_dist["min_length"] == 8
        assert length_dist["max_length"] == 8
        assert length_dist["mean_length"] == 8.0

    def test_extreme_quality_scores(self):
        """Test handling of extreme quality scores."""
        reads = [
            FastqRecord("low", "ATCG", "!!!!"),  # Quality 0
            FastqRecord("high", "GCTA", "~~~~"),  # Quality 93 (very high)
        ]

        result = quality_score_distribution(reads)

        assert result["quality_score_min"] == 0
        assert result["quality_score_max"] == 93

        # Should handle without errors
        assert "mean_quality" in result
        assert "median_quality" in result
