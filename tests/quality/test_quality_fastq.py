"""Tests for FASTQ quality control module.

Real implementation testing for FASTQ quality assessment functions.
No mocking used - all tests use real data and computational methods.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.quality.io.fastq import (
    FastqRecord,
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


def _make_record(name: str, sequence: str, quality: str) -> FastqRecord:
    """Helper to create a FastqRecord with proper FASTQ formatting.

    Adds the required '@' prefix to the header and '+' quality header.
    """
    return FastqRecord(f"@{name}", sequence, "+", quality)


class TestBasicStatistics:
    """Test basic FASTQ statistics calculation."""

    def setup_method(self):
        """Setup test FASTQ records with known properties."""
        self.test_reads = [
            _make_record("read1", "ATCGATCG", "IIIIIIII"),  # Length 8, high quality
            _make_record("read2", "GCGCGCGC", "HHHHHHHH"),  # Length 8, medium quality
            _make_record("read3", "ATGCATGCAT", "IIIIIIIIII"),  # Length 10, high quality
            _make_record("read4", "NNNNNNNN", "!!!!!!!!"),  # Length 8, low quality, all N
        ]

    def test_basic_statistics_calculation(self):
        """Test calculation of basic FASTQ statistics."""
        stats = basic_statistics(self.test_reads)

        # Check basic counts and lengths
        assert stats["total_reads"] == 4
        assert stats["min_length"] == 8
        assert stats["max_length"] == 10
        assert stats["mean_length"] == 8.5  # (8+8+10+8)/4
        assert stats["median_length"] == 8.0

        # Check GC content calculation
        # read1: ATCGATCG -> A=2, T=2, C=2, G=2 -> GC = 4/8 = 50%
        # read2: GCGCGCGC -> G=4, C=4 -> GC = 8/8 = 100%
        # read3: ATGCATGCAT -> A=3, T=3, G=2, C=2 -> GC = 4/10 = 40%
        # read4: NNNNNNNN -> All N, no ATCG -> GC = 0/8 = 0%
        expected_gc = (50 + 100 + 40 + 0) / 4  # 47.5%
        assert abs(stats["mean_gc"] - expected_gc) < 0.1

        # Check quality scores
        assert "mean_quality" in stats
        assert "median_quality" in stats
        assert stats["min_quality"] >= 0
        assert stats["max_quality"] <= 41  # Phred scale

    def test_basic_statistics_uniform_length(self):
        """Test statistics with uniform-length reads."""
        uniform_reads = [
            _make_record("read1", "ATCG", "IIII"),
            _make_record("read2", "GCTA", "HHHH"),
            _make_record("read3", "TTTT", "JJJJ"),
        ]

        stats = basic_statistics(uniform_reads)
        # Source does not return 'uniform_length'; verify uniform by checking min == max
        assert stats["min_length"] == 4
        assert stats["max_length"] == 4
        assert stats["mean_length"] == 4.0

    def test_basic_statistics_empty_input(self):
        """Test handling of empty input."""
        stats = basic_statistics([])
        assert stats == {}


class TestPerBaseQuality:
    """Test per-base quality analysis."""

    def test_per_base_quality_calculation(self):
        """Test per-base quality score calculation."""
        reads = [
            _make_record("read1", "ATCG", "IIII"),  # Quality 40,40,40,40
            _make_record("read2", "GCTA", "HHHH"),  # Quality 39,39,39,39
            _make_record("read3", "TTGG", "JJJJ"),  # Quality 41,41,41,41
        ]

        result = per_base_quality(reads)

        assert "positions" in result
        positions = result["positions"]
        assert len(positions) == 4

        # Check first position statistics
        pos1_stats = positions[0]
        assert pos1_stats["position"] == 1
        assert pos1_stats["count"] == 3  # All 3 reads have this position

        # Quality scores: I=40, H=39, J=41, so mean should be 40
        expected_mean = (40 + 39 + 41) / 3
        assert abs(pos1_stats["mean"] - expected_mean) < 0.1

    def test_per_base_quality_variable_lengths(self):
        """Test per-base quality with variable read lengths."""
        reads = [
            _make_record("read1", "AT", "II"),  # Length 2
            _make_record("read2", "GCT", "HHH"),  # Length 3
            _make_record("read3", "ATGC", "JJJJ"),  # Length 4
        ]

        result = per_base_quality(reads)

        positions = result["positions"]
        max_position = len(positions)
        assert max_position == 4

        # Position 1: all 3 reads
        assert positions[0]["count"] == 3
        # Position 2: all 3 reads
        assert positions[1]["count"] == 3
        # Position 3: 2 reads (read2, read3)
        assert positions[2]["count"] == 2
        # Position 4: 1 read (read3)
        assert positions[3]["count"] == 1


class TestSequenceLengthDistribution:
    """Test sequence length distribution analysis."""

    def test_length_distribution_calculation(self):
        """Test sequence length distribution calculation."""
        reads = [
            _make_record("read1", "ATCG", "IIII"),  # Length 4
            _make_record("read2", "ATCGATCG", "IIIIIIII"),  # Length 8
            _make_record("read3", "ATCG", "IIII"),  # Length 4 (duplicate)
            _make_record("read4", "ATCGATCGATCG", "IIIIIIIIIIII"),  # Length 12
        ]

        result = sequence_length_distribution(reads)

        # Source returns {"distribution": [{"length": N, "count": N, "percentage": N}, ...]}
        assert "distribution" in result
        distribution = result["distribution"]

        # Build a lookup by length for assertions
        dist_by_length = {entry["length"]: entry for entry in distribution}

        assert 4 in dist_by_length
        assert dist_by_length[4]["count"] == 2  # Two reads of length 4
        assert 8 in dist_by_length
        assert dist_by_length[8]["count"] == 1  # One read of length 8
        assert 12 in dist_by_length
        assert dist_by_length[12]["count"] == 1  # One read of length 12

        # Verify min/max from distribution
        lengths = [entry["length"] for entry in distribution]
        assert min(lengths) == 4
        assert max(lengths) == 12

        # Verify percentages sum to ~100%
        total_pct = sum(entry["percentage"] for entry in distribution)
        assert abs(total_pct - 100.0) < 0.1


class TestGCContentDistribution:
    """Test GC content distribution analysis."""

    def test_gc_content_calculation(self):
        """Test GC content calculation."""
        reads = [
            _make_record("read1", "ATCG", "IIII"),  # 50% GC (2/4)
            _make_record("read2", "AAAA", "IIII"),  # 0% GC (0/4)
            _make_record("read3", "GCGA", "IIII"),  # 75% GC (3/4)
            _make_record("read4", "ATNNNN", "IIIIII"),  # 0% GC (0/6 total chars)
        ]

        result = gc_content_distribution(reads)

        # Source returns {"bins": [{"bin_start": N, "bin_end": N, "count": N, "percentage": N}, ...]}
        assert "bins" in result
        bins = result["bins"]

        # Should have bins covering 0%, 50%, and 75% GC values
        assert len(bins) >= 1

        # Total count across bins should equal number of reads
        total_count = sum(b["count"] for b in bins)
        assert total_count == 4

        # Each bin should have bin_start, bin_end, count, percentage
        for b in bins:
            assert "bin_start" in b
            assert "bin_end" in b
            assert "count" in b
            assert "percentage" in b
            assert b["count"] > 0


class TestAdapterContent:
    """Test adapter contamination detection."""

    def test_adapter_detection_basic(self):
        """Test basic adapter detection with default adapters."""
        reads = [
            _make_record("read1", "ATCGATCGATCG", "IIIIIIIIIIII"),  # No adapter
            _make_record("read2", "AGATCGGAAGAGATCG", "IIIIIIIIIIIIIIII"),  # Illumina adapter
            _make_record("read3", "ATCGAGATCGGAAGAG", "IIIIIIIIIIIIIIII"),  # Adapter at end
        ]

        result = adapter_content(reads)

        # Source returns {"adapters": {"adapter_1": {"adapter_sequence": "...", "positions": [...]}, ...}}
        assert "adapters" in result
        adapters_result = result["adapters"]

        # Should have adapter entries (default adapters are a list)
        assert len(adapters_result) >= 1

        # Each adapter entry should have adapter_sequence and positions
        for adapter_key, adapter_data in adapters_result.items():
            assert "adapter_sequence" in adapter_data
            assert "positions" in adapter_data
            assert len(adapter_data["positions"]) > 0

    def test_custom_adapters(self):
        """Test custom adapter sequences (passed as list)."""
        reads = [
            _make_record("read1", "ATCGCUSTOMADAPTER", "IIIIIIIIIIIIIIIII"),
        ]

        # Source accepts adapters as List[str], not Dict[str, str]
        custom_adapters = ["CUSTOMADAPTER"]
        result = adapter_content(reads, adapters=custom_adapters)

        assert "adapters" in result
        adapters_result = result["adapters"]
        assert len(adapters_result) >= 1


class TestOverrepresentedSequences:
    """Test overrepresented sequence detection."""

    def test_overrepresented_detection(self):
        """Test detection of overrepresented sequences."""
        # Create reads with some duplicates
        reads = []

        # Add 5 copies of the same sequence (overrepresented) - must be >= min_length (20)
        overrep_seq = "ATCGATCGATCGATCGATCG"  # Exactly 20 chars
        for i in range(5):
            reads.append(_make_record(f"dup_{i}", overrep_seq, "I" * len(overrep_seq)))

        # Add 10 unique sequences of length >= 20
        for i in range(10):
            unique_seq = "A" * 20 + chr(ord("A") + i % 4)  # 21 chars, varying last char
            reads.append(_make_record(f"unique_{i}", unique_seq, "I" * len(unique_seq)))

        # Source signature: overrepresented_sequences(records, min_length=20)
        result = overrepresented_sequences(reads)

        # Source returns {"overrepresented": [{"sequence": "...", "count": N, "percentage": N}, ...]}
        assert "overrepresented" in result
        overrep_list = result["overrepresented"]

        # The duplicated sequence should be detected
        found_overrep = False
        for entry in overrep_list:
            if entry["sequence"] == overrep_seq:
                found_overrep = True
                assert entry["count"] == 5
                expected_pct = (5 / 15) * 100  # 33.33%
                assert abs(entry["percentage"] - expected_pct) < 0.1
                break

        assert found_overrep, "Expected overrepresented sequence not found"


class TestDuplicationLevels:
    """Test sequence duplication analysis."""

    def test_duplication_analysis(self):
        """Test sequence duplication level analysis."""
        reads = [
            # Unique sequences (appear once)
            _make_record("unique1", "ATCG", "IIII"),
            _make_record("unique2", "GCTA", "IIII"),
            # Duplicated sequence (appears twice)
            _make_record("dup1", "TTTT", "IIII"),
            _make_record("dup2", "TTTT", "IIII"),
            # Highly duplicated (appears 3 times)
            _make_record("high1", "AAAA", "IIII"),
            _make_record("high2", "AAAA", "IIII"),
            _make_record("high3", "AAAA", "IIII"),
        ]

        result = duplication_levels(reads)

        assert result["total_sequences"] == 7
        assert result["unique_sequences"] == 4  # 4 distinct sequences

        # Source returns duplication_bins with string keys
        dup_bins = result["duplication_bins"]

        # Should have: 2 sequences appearing 1x, 1 sequence appearing 2x, 1 sequence appearing 3x
        assert "1" in dup_bins  # Unique sequences
        assert "2" in dup_bins  # Duplicated sequences
        assert "3" in dup_bins  # Highly duplicated sequences

        assert dup_bins["1"] == 2  # ATCG, GCTA
        assert dup_bins["2"] == 1  # TTTT
        assert dup_bins["3"] == 1  # AAAA

        # Check duplication rate
        assert "duplication_rate" in result
        # 4 unique out of 7 total -> duplication_rate = (1 - 4/7) * 100
        expected_rate = (1 - 4 / 7) * 100
        assert abs(result["duplication_rate"] - expected_rate) < 0.1


class TestQualityScoreDistribution:
    """Test quality score distribution analysis."""

    def test_quality_distribution_calculation(self):
        """Test quality score distribution calculation."""
        reads = [
            _make_record("read1", "ATCG", "!!!!"),  # Quality 0,0,0,0
            _make_record("read2", "GCTA", "IIII"),  # Quality 40,40,40,40
            _make_record("read3", "TTGG", "@@@@"),  # Quality 31,31,31,31
        ]

        result = quality_score_distribution(reads)

        # Source returns {"distribution": [{"quality_score": N, "count": N, "percentage": N}, ...]}
        assert "distribution" in result
        distribution = result["distribution"]

        # Build lookup by quality score
        dist_by_score = {entry["quality_score"]: entry for entry in distribution}

        # Check expected quality scores are present
        assert 0 in dist_by_score
        assert dist_by_score[0]["count"] == 4  # From read1 (4 bases)
        assert 40 in dist_by_score
        assert dist_by_score[40]["count"] == 4  # From read2 (4 bases)
        assert 31 in dist_by_score
        assert dist_by_score[31]["count"] == 4  # From read3 (4 bases)

        # Total count should be 12 (4 bases * 3 reads)
        total_count = sum(entry["count"] for entry in distribution)
        assert total_count == 12


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

        # Check that all expected analyses are present (actual return keys)
        expected_analyses = [
            "basic_statistics",
            "per_base_quality",
            "per_sequence_quality",
            "sequence_length_distribution",
            "gc_content_distribution",
            "adapter_content",
            "overrepresented_sequences",
            "duplication_levels",
            "n_content_per_position",
            "quality_score_distribution",
        ]

        for analysis in expected_analyses:
            assert analysis in results, f"Missing analysis: {analysis}"

        # Basic validation of results
        basic_stats = results["basic_statistics"]
        assert basic_stats["total_reads"] == 4
        # All reads are length 16, verify min == max
        assert basic_stats["min_length"] == 16
        assert basic_stats["max_length"] == 16

    def test_analyze_fastq_quality_with_n_reads(self, tmp_path):
        """Test FASTQ analysis with n_reads limit."""
        # Create FASTQ with more reads than n_reads limit
        fastq_content = ""
        for i in range(10):
            fastq_content += f"@read{i}\nATCGATCGATCGATCG\n+\nIIIIIIIIIIIIIIII\n"

        fastq_path = tmp_path / "test_large.fastq"
        fastq_path.write_text(fastq_content)

        # Run analysis with n_reads limit (source uses n_reads, not subsample)
        results = analyze_fastq_quality(fastq_path, n_reads=5)

        # Should only analyze 5 reads
        basic_stats = results["basic_statistics"]
        assert basic_stats["total_reads"] == 5


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_fastq_file(self, tmp_path):
        """Test handling of empty FASTQ file."""
        empty_path = tmp_path / "empty.fastq"
        empty_path.write_text("")

        # Source raises ValueError("No reads found in FASTQ file") for empty files
        with pytest.raises((ValueError, Exception), match="No reads found|empty|no read"):
            analyze_fastq_quality(empty_path)

    def test_nonexistent_file(self):
        """Test handling of nonexistent file."""
        # Source uses validate_path_exists which raises ValidationError, not FileNotFoundError
        with pytest.raises(Exception):
            analyze_fastq_quality("nonexistent.fastq")

    def test_single_read_analysis(self):
        """Test analysis with single read."""
        reads = [_make_record("single", "ATCGATCG", "IIIIIIII")]

        stats = basic_statistics(reads)
        assert stats["total_reads"] == 1
        assert stats["min_length"] == 8
        assert stats["max_length"] == 8

        length_result = sequence_length_distribution(reads)
        distribution = length_result["distribution"]
        assert len(distribution) == 1
        assert distribution[0]["length"] == 8
        assert distribution[0]["count"] == 1

    def test_extreme_quality_scores(self):
        """Test handling of extreme quality scores."""
        reads = [
            _make_record("low", "ATCG", "!!!!"),  # Quality 0
            _make_record("high", "GCTA", "~~~~"),  # Quality 93 (very high)
        ]

        result = quality_score_distribution(reads)

        # Source returns {"distribution": [{"quality_score": N, "count": N, "percentage": N}, ...]}
        assert "distribution" in result
        distribution = result["distribution"]

        dist_by_score = {entry["quality_score"]: entry for entry in distribution}

        assert 0 in dist_by_score
        assert 93 in dist_by_score

        # Should handle without errors
        assert len(distribution) >= 2
