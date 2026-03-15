"""Tests for GWAS summary statistics output module."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from metainformant.gwas.analysis.summary_stats import (
    _compute_lambda_gc,
    create_results_summary,
    write_significant_hits,
    write_summary_statistics,
)


@pytest.fixture
def sample_results() -> list:
    """Generate sample association results."""
    return [
        {"beta": 0.5, "se": 0.1, "p_value": 1e-6, "n_samples": 100, "maf": 0.3},
        {"beta": 0.1, "se": 0.05, "p_value": 0.5, "n_samples": 100, "maf": 0.2},
        {"beta": -0.3, "se": 0.08, "p_value": 1e-10, "n_samples": 100, "maf": 0.15},
        {"beta": 0.0, "se": 0.2, "p_value": 0.99, "n_samples": 100, "maf": 0.4},
    ]


@pytest.fixture
def sample_variant_info() -> list:
    """Generate sample variant info."""
    return [
        {"chrom": "chr1", "pos": 1000, "id": "rs1", "ref": "A", "alt": ["G"]},
        {"chrom": "chr1", "pos": 2000, "id": "rs2", "ref": "T", "alt": ["C"]},
        {"chrom": "chr2", "pos": 3000, "id": "rs3", "ref": "G", "alt": ["A"]},
        {"chrom": "chr2", "pos": 4000, "id": "rs4", "ref": "C", "alt": ["T"]},
    ]


class TestWriteSummaryStatistics:
    """Tests for summary statistics TSV output."""

    def test_write_basic(self, tmp_path: Path, sample_results: list, sample_variant_info: list) -> None:
        """Write and verify basic summary statistics."""
        output = tmp_path / "stats.tsv"
        result_path = write_summary_statistics(sample_results, sample_variant_info, output)

        assert result_path.exists()
        lines = result_path.read_text().strip().split("\n")
        assert len(lines) == 5  # header + 4 variants

        # Check header
        header = lines[0].split("\t")
        assert "CHR" in header
        assert "POS" in header
        assert "P" in header
        assert "BETA" in header

    def test_write_creates_directories(self, tmp_path: Path, sample_results: list, sample_variant_info: list) -> None:
        """Should create parent directories if needed."""
        output = tmp_path / "nested" / "dir" / "stats.tsv"
        write_summary_statistics(sample_results, sample_variant_info, output)
        assert output.exists()

    def test_mismatched_lengths(self, tmp_path: Path) -> None:
        """Mismatched results/variant_info lengths should raise ValueError."""
        with pytest.raises(ValueError):
            write_summary_statistics(
                [{"beta": 0.1}],
                [{"chrom": "1"}, {"chrom": "2"}],
                tmp_path / "bad.tsv",
            )

    def test_alt_as_list(self, tmp_path: Path) -> None:
        """Alt alleles as list should be joined with comma."""
        results = [{"beta": 0.1, "se": 0.05, "p_value": 0.5, "n_samples": 50, "maf": 0.2}]
        variants = [{"chrom": "1", "pos": 100, "id": "rs1", "ref": "A", "alt": ["G", "T"]}]
        output = tmp_path / "multi_alt.tsv"
        write_summary_statistics(results, variants, output)

        content = output.read_text()
        assert "G,T" in content


class TestWriteSignificantHits:
    """Tests for significant hits output."""

    def test_filter_significant(self, tmp_path: Path, sample_results: list, sample_variant_info: list) -> None:
        """Only significant hits should be written."""
        output = tmp_path / "sig.tsv"
        write_significant_hits(sample_results, sample_variant_info, output, threshold=1e-5)

        lines = output.read_text().strip().split("\n")
        # header + 1 hit (p=1e-10 is below 1e-5, but p=1e-6 is also below)
        assert len(lines) == 3  # header + 2 significant

    def test_no_significant(self, tmp_path: Path) -> None:
        """If nothing is significant, only header should be written."""
        results = [{"beta": 0.1, "se": 0.1, "p_value": 0.5, "n_samples": 50, "maf": 0.2}]
        variants = [{"chrom": "1", "pos": 100, "id": "rs1", "ref": "A", "alt": "G"}]
        output = tmp_path / "no_sig.tsv"
        write_significant_hits(results, variants, output, threshold=5e-8)

        lines = output.read_text().strip().split("\n")
        assert len(lines) == 1  # Just header

    def test_sorted_by_pvalue(self, tmp_path: Path, sample_results: list, sample_variant_info: list) -> None:
        """Significant hits should be sorted by p-value ascending."""
        output = tmp_path / "sorted.tsv"
        write_significant_hits(sample_results, sample_variant_info, output, threshold=0.01)

        lines = output.read_text().strip().split("\n")
        # Skip header, parse p-values
        p_values = []
        for line in lines[1:]:
            cols = line.split("\t")
            p_values.append(float(cols[7]))  # P column

        # Verify sorted ascending
        assert p_values == sorted(p_values)


class TestCreateResultsSummary:
    """Tests for JSON results summary."""

    def test_basic_summary(self, tmp_path: Path, sample_results: list) -> None:
        """Create and verify JSON summary."""
        output = tmp_path / "summary.json"
        summary = create_results_summary(sample_results, output)

        assert output.exists()
        assert summary["n_variants_tested"] == 4
        assert "lambda_gc" in summary
        assert "n_significant_genome_wide" in summary
        assert "top_hits" in summary

        # Verify JSON is valid
        with open(output) as f:
            loaded = json.load(f)
        assert loaded["n_variants_tested"] == 4

    def test_lambda_gc_reasonable(self, tmp_path: Path) -> None:
        """Lambda GC should be near 1.0 for random p-values."""
        import random

        random.seed(42)
        results = [{"p_value": random.random()} for _ in range(1000)]
        output = tmp_path / "summary.json"
        summary = create_results_summary(results, output)

        # For uniform p-values, lambda should be close to 1
        assert 0.5 < summary["lambda_gc"] < 2.0

    def test_top_hits_limited(self, tmp_path: Path) -> None:
        """Top hits should be limited to 20."""
        results = [{"p_value": 1.0 / (i + 1), "beta": 0.1} for i in range(100)]
        output = tmp_path / "summary.json"
        summary = create_results_summary(results, output)

        assert len(summary["top_hits"]) <= 20


class TestComputeLambdaGC:
    """Tests for genomic inflation factor calculation."""

    def test_empty(self) -> None:
        """Empty p-values should return 1.0."""
        assert _compute_lambda_gc([]) == 1.0

    def test_all_significant(self) -> None:
        """Very small p-values should give lambda >> 1."""
        p_values = [1e-10] * 100
        lambda_gc = _compute_lambda_gc(p_values)
        assert lambda_gc > 5  # Strong inflation

    def test_no_inflation(self) -> None:
        """Uniform p-values should give lambda near 1."""
        import numpy as np

        np.random.seed(42)
        p_values = list(np.random.uniform(0, 1, 10000))
        lambda_gc = _compute_lambda_gc(p_values)
        assert 0.8 < lambda_gc < 1.2
