"""Tests for longread assembly summary functions not covered by test_longread.py.

Tests generate_assembly_summary, generate_sv_summary, and related summary functions.
All tests use real implementations -- NO MOCKING.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from metainformant.longread.utils.summary import (
    RunSummary,
    build_run_summary,
    compare_run_summaries,
    export_run_summary,
    generate_assembly_summary,
    generate_sv_summary,
)

# ---------------------------------------------------------------------------
# generate_assembly_summary
# ---------------------------------------------------------------------------


class TestGenerateAssemblySummary:
    """Tests for generate_assembly_summary function."""

    def test_assembly_summary_basic(self) -> None:
        """generate_assembly_summary computes stats from contig lengths."""
        contig_lengths = [500000, 300000, 200000, 100000, 50000]
        summary = generate_assembly_summary(contig_lengths)

        assert summary["total_contigs"] == 5
        assert summary["total_bases"] == 1150000
        assert summary["largest_contig"] == 500000
        assert summary["smallest_contig"] == 50000
        assert summary["mean_length"] == pytest.approx(230000.0)
        assert summary["n50"] > 0
        assert summary["n90"] > 0
        assert summary["status"] == "complete"

    def test_assembly_summary_empty(self) -> None:
        """generate_assembly_summary returns no_assembly status for empty input."""
        summary = generate_assembly_summary([])
        assert summary["total_contigs"] == 0
        assert summary["total_bases"] == 0
        assert summary["status"] == "no_assembly"

    def test_assembly_summary_single_contig(self) -> None:
        """generate_assembly_summary works with a single contig."""
        summary = generate_assembly_summary([1000000])
        assert summary["total_contigs"] == 1
        assert summary["total_bases"] == 1000000
        assert summary["largest_contig"] == 1000000
        assert summary["smallest_contig"] == 1000000
        assert summary["n50"] == 1000000

    def test_assembly_summary_with_metadata(self) -> None:
        """generate_assembly_summary includes reads_used and polish info."""
        contig_lengths = [100000, 200000, 300000]
        summary = generate_assembly_summary(
            contig_lengths,
            num_reads_used=5000,
            polish_iterations=3,
            coverage=30.5,
        )
        assert summary["reads_used"] == 5000
        assert summary["polish_iterations"] == 3
        assert summary["mean_coverage"] == 30.5

    def test_assembly_summary_sorted_descending(self) -> None:
        """generate_assembly_summary sorts contigs descending for N50 calculation."""
        # Provide unsorted input -- function should handle it
        contig_lengths = [1000, 50000, 10000, 30000]
        summary = generate_assembly_summary(contig_lengths)
        assert summary["largest_contig"] == 50000
        assert summary["smallest_contig"] == 1000
        # N50: total=91000, threshold=45500
        # Sorted desc: 50000 cumulative 50000 >= 45500 -> N50=50000
        assert summary["n50"] == 50000

    def test_assembly_summary_gc_content_is_none(self) -> None:
        """generate_assembly_summary sets gc_content to None (computed separately)."""
        summary = generate_assembly_summary([10000, 20000])
        assert summary["gc_content"] is None


# ---------------------------------------------------------------------------
# generate_sv_summary
# ---------------------------------------------------------------------------


class TestGenerateSvSummary:
    """Tests for generate_sv_summary function."""

    def test_sv_summary_basic(self) -> None:
        """generate_sv_summary computes correct counts and size distributions."""
        variants = [
            {"sv_type": "DEL", "size": 500},
            {"sv_type": "DEL", "size": 1000},
            {"sv_type": "DEL", "size": 800},
            {"sv_type": "INS", "size": 200},
            {"sv_type": "INS", "size": 350},
            {"sv_type": "INV", "size": 50000},
        ]
        summary = generate_sv_summary(variants)

        assert summary["total_variants"] == 6
        assert summary["by_type"]["DEL"] == 3
        assert summary["by_type"]["INS"] == 2
        assert summary["by_type"]["INV"] == 1
        assert summary["status"] == "complete"

        # Size distribution should be present
        size_dist = summary["size_distribution"]
        assert size_dist["min_size"] == 200
        assert size_dist["max_size"] == 50000

    def test_sv_summary_empty(self) -> None:
        """generate_sv_summary returns no_variants for empty input."""
        summary = generate_sv_summary([])
        assert summary["total_variants"] == 0
        assert summary["status"] == "no_variants"

    def test_sv_summary_phased_count(self) -> None:
        """generate_sv_summary counts phased variants."""
        variants = [
            {"sv_type": "DEL", "size": 500, "haplotype": 1},
            {"sv_type": "DEL", "size": 800, "haplotype": 2},
            {"sv_type": "INS", "size": 300, "haplotype": 0},  # unphased
        ]
        summary = generate_sv_summary(variants)
        assert summary["phased_count"] == 2

    def test_sv_summary_single_variant(self) -> None:
        """generate_sv_summary works with a single variant."""
        variants = [{"sv_type": "INV", "size": 10000}]
        summary = generate_sv_summary(variants)
        assert summary["total_variants"] == 1
        assert summary["by_type"]["INV"] == 1
        assert summary["size_distribution"]["min_size"] == 10000
        assert summary["size_distribution"]["max_size"] == 10000

    def test_sv_summary_missing_size(self) -> None:
        """generate_sv_summary handles variants with missing size field."""
        variants = [
            {"sv_type": "BND"},
            {"sv_type": "DEL", "size": 500},
        ]
        summary = generate_sv_summary(variants)
        assert summary["total_variants"] == 2
        # Only 1 variant has a valid size
        size_dist = summary["size_distribution"]
        assert size_dist["min_size"] == 500
        assert size_dist["max_size"] == 500


# ---------------------------------------------------------------------------
# RunSummary dataclass
# ---------------------------------------------------------------------------


class TestRunSummary:
    """Tests for RunSummary dataclass."""

    def test_run_summary_default_values(self) -> None:
        """RunSummary has correct defaults."""
        summary = RunSummary()
        assert summary.run_id == ""
        assert summary.sample_name == ""
        assert summary.platform == "ont"
        assert summary.timestamp == ""
        assert summary.input_stats == {}
        assert summary.qc_stats == {}
        assert summary.assembly_stats == {}
        assert summary.methylation_stats == {}
        assert summary.sv_stats == {}
        assert summary.phasing_stats == {}
        assert summary.pipeline_steps == []
        assert summary.warnings == []

    def test_run_summary_custom_values(self) -> None:
        """RunSummary stores custom values correctly."""
        summary = RunSummary(
            run_id="lr_12345_sample1",
            sample_name="sample1",
            platform="pacbio",
            qc_stats={"total_reads": 10000, "n50": 8000},
            warnings=["Low coverage in region X"],
        )
        assert summary.run_id == "lr_12345_sample1"
        assert summary.platform == "pacbio"
        assert summary.qc_stats["total_reads"] == 10000
        assert len(summary.warnings) == 1


# ---------------------------------------------------------------------------
# build_run_summary
# ---------------------------------------------------------------------------


class TestBuildRunSummary:
    """Tests for build_run_summary function."""

    def test_build_run_summary_minimal(self) -> None:
        """build_run_summary creates a RunSummary with minimal input."""
        summary = build_run_summary("test_sample")
        assert isinstance(summary, RunSummary)
        assert summary.sample_name == "test_sample"
        assert summary.platform == "ont"
        assert summary.run_id.startswith("lr_")
        assert "test_sample" in summary.run_id
        assert summary.timestamp != ""

    def test_build_run_summary_with_all_stats(self) -> None:
        """build_run_summary populates all stats sections."""
        qc = {"total_reads": 5000, "n50": 10000}
        assembly = {"total_contigs": 3, "n50": 500000}
        meth = {"total_sites": 100000, "methylated_sites": 70000}
        sv = {"total_variants": 50}
        phasing = {"num_blocks": 10}

        summary = build_run_summary(
            "full_sample",
            platform="pacbio",
            qc_stats=qc,
            assembly_stats=assembly,
            methylation_stats=meth,
            sv_stats=sv,
            phasing_stats=phasing,
        )

        assert summary.platform == "pacbio"
        assert summary.qc_stats == qc
        assert summary.assembly_stats == assembly
        assert summary.methylation_stats == meth
        assert summary.sv_stats == sv
        assert summary.phasing_stats == phasing


# ---------------------------------------------------------------------------
# export_run_summary
# ---------------------------------------------------------------------------


class TestExportRunSummary:
    """Tests for export_run_summary function."""

    def test_export_run_summary_json(self, tmp_path: Path) -> None:
        """export_run_summary writes valid JSON file."""
        import json

        summary = build_run_summary(
            "json_test",
            qc_stats={"total_reads": 1000, "n50": 5000, "mean_length": 4500.0},
        )
        output_file = tmp_path / "summary.json"
        result = export_run_summary(summary, output_file, format="json")

        assert result == output_file
        assert output_file.exists()
        data = json.loads(output_file.read_text())
        assert data["sample_name"] == "json_test"
        assert data["qc_stats"]["total_reads"] == 1000

    def test_export_run_summary_text(self, tmp_path: Path) -> None:
        """export_run_summary writes human-readable text file."""
        summary = build_run_summary(
            "text_test",
            qc_stats={"total_reads": 2000, "total_bases": 10000000, "n50": 8000, "mean_length": 5000.0},
            assembly_stats={"total_contigs": 5, "total_bases": 4000000, "n50": 900000, "largest_contig": 1500000},
        )
        output_file = tmp_path / "summary.txt"
        result = export_run_summary(summary, output_file, format="text")

        assert result == output_file
        assert output_file.exists()
        content = output_file.read_text()
        assert "LONG-READ SEQUENCING RUN SUMMARY" in content
        assert "text_test" in content
        assert "Quality Control" in content
        assert "Assembly" in content

    def test_export_run_summary_text_with_methylation(self, tmp_path: Path) -> None:
        """export_run_summary text output includes methylation section."""
        summary = build_run_summary(
            "meth_test",
            methylation_stats={
                "modification_type": "5mC",
                "total_sites": 50000,
                "methylated_sites": 35000,
                "global_methylation_rate": 0.7,
            },
        )
        output_file = tmp_path / "meth_summary.txt"
        export_run_summary(summary, output_file, format="text")

        content = output_file.read_text()
        assert "Methylation" in content
        assert "5mC" in content

    def test_export_run_summary_text_with_sv(self, tmp_path: Path) -> None:
        """export_run_summary text output includes SV section."""
        summary = build_run_summary(
            "sv_test",
            sv_stats={"total_variants": 42, "by_type": {"DEL": 20, "INS": 15, "INV": 7}},
        )
        output_file = tmp_path / "sv_summary.txt"
        export_run_summary(summary, output_file, format="text")

        content = output_file.read_text()
        assert "Structural Variants" in content

    def test_export_run_summary_text_with_warnings(self, tmp_path: Path) -> None:
        """export_run_summary text output includes warnings section."""
        summary = build_run_summary("warn_test")
        summary.warnings = ["Low coverage detected", "Adapter contamination high"]

        output_file = tmp_path / "warn_summary.txt"
        export_run_summary(summary, output_file, format="text")

        content = output_file.read_text()
        assert "Warnings" in content
        assert "Low coverage detected" in content

    def test_export_run_summary_invalid_format(self, tmp_path: Path) -> None:
        """export_run_summary raises ValueError for unsupported format."""
        summary = build_run_summary("test")
        with pytest.raises(ValueError, match="Unsupported format"):
            export_run_summary(summary, tmp_path / "out.xyz", format="csv")

    def test_export_run_summary_creates_parent_dirs(self, tmp_path: Path) -> None:
        """export_run_summary creates parent directories if needed."""
        summary = build_run_summary("dir_test")
        output_file = tmp_path / "nested" / "dir" / "summary.json"
        result = export_run_summary(summary, output_file, format="json")

        assert result == output_file
        assert output_file.exists()


# ---------------------------------------------------------------------------
# compare_run_summaries
# ---------------------------------------------------------------------------


class TestCompareRunSummaries:
    """Tests for compare_run_summaries function."""

    def test_compare_empty(self) -> None:
        """compare_run_summaries returns minimal result for empty list."""
        result = compare_run_summaries([])
        assert result["samples"] == 0

    def test_compare_single_summary(self) -> None:
        """compare_run_summaries works with a single summary."""
        summary = RunSummary(
            sample_name="sample1",
            platform="ont",
            qc_stats={"n50": 8000, "total_bases": 5000000, "total_reads": 1000, "mean_quality": 18.5},
        )
        result = compare_run_summaries([summary])

        assert result["samples"] == 1
        assert result["platforms"] == ["ont"]
        assert "sample1" in result["per_sample"]
        assert result["per_sample"]["sample1"]["n50"] == 8000
        assert result["aggregate_n50"]["min"] == 8000
        assert result["aggregate_n50"]["max"] == 8000

    def test_compare_multiple_summaries(self) -> None:
        """compare_run_summaries computes aggregates across multiple samples."""
        summaries = [
            RunSummary(
                sample_name="sample_a",
                platform="ont",
                qc_stats={"n50": 5000, "total_bases": 3000000, "total_reads": 500, "mean_quality": 15.0},
            ),
            RunSummary(
                sample_name="sample_b",
                platform="ont",
                qc_stats={"n50": 10000, "total_bases": 8000000, "total_reads": 1000, "mean_quality": 20.0},
            ),
            RunSummary(
                sample_name="sample_c",
                platform="pacbio",
                qc_stats={"n50": 15000, "total_bases": 12000000, "total_reads": 800, "mean_quality": 25.0},
            ),
        ]
        result = compare_run_summaries(summaries)

        assert result["samples"] == 3
        assert set(result["platforms"]) == {"ont", "pacbio"}

        # Aggregate N50
        assert result["aggregate_n50"]["min"] == 5000
        assert result["aggregate_n50"]["max"] == 15000
        assert result["aggregate_n50"]["mean"] == round((5000 + 10000 + 15000) / 3)

        # Aggregate bases
        assert result["aggregate_bases"]["total"] == 23000000

        # Aggregate quality
        assert result["aggregate_quality"]["min"] == 15.0
        assert result["aggregate_quality"]["max"] == 25.0

    def test_compare_summaries_missing_qc(self) -> None:
        """compare_run_summaries handles summaries without QC stats."""
        summaries = [
            RunSummary(sample_name="no_qc", platform="ont"),
            RunSummary(
                sample_name="with_qc",
                platform="ont",
                qc_stats={"n50": 5000, "total_bases": 1000000, "total_reads": 200},
            ),
        ]
        result = compare_run_summaries(summaries)

        assert result["samples"] == 2
        assert "no_qc" in result["per_sample"]
        assert "with_qc" in result["per_sample"]
        # Only one sample contributes to aggregates
        assert result["aggregate_n50"]["min"] == 5000
        assert result["aggregate_n50"]["max"] == 5000
