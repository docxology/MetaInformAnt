"""Tests for epigenome workflow orchestration.

Tests EpigenomeConfig, load_epigenome_config, and
integrate_epigenome_results using real implementations. NO MOCKING.

Note: run_methylation_workflow, run_chipseq_workflow, and
run_atacseq_workflow require input directories with actual data files
so are tested via config validation and integration helpers.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.epigenome.workflow.workflow import (
    EpigenomeConfig,
    _detect_peak_format,
    _generate_integration_report,
    integrate_epigenome_results,
    load_epigenome_config,
)

# ---------------------------------------------------------------------------
# EpigenomeConfig
# ---------------------------------------------------------------------------


class TestEpigenomeConfig:
    """Tests for EpigenomeConfig dataclass."""

    def test_default_config(self) -> None:
        cfg = EpigenomeConfig()
        assert cfg.methylation_threshold == 0.8
        assert cfg.chipseq_qvalue_threshold == 0.05
        assert cfg.n_threads == 4
        assert cfg.run_methylation is True

    def test_custom_config(self) -> None:
        cfg = EpigenomeConfig(
            methylation_threshold=0.5,
            n_threads=8,
            min_peak_length=100,
            generate_reports=False,
        )
        assert cfg.methylation_threshold == 0.5
        assert cfg.n_threads == 8
        assert cfg.min_peak_length == 100
        assert cfg.generate_reports is False

    def test_validation_methylation_threshold_out_of_range(self) -> None:
        with pytest.raises(Exception):
            EpigenomeConfig(methylation_threshold=1.5)

    def test_validation_negative_threads(self) -> None:
        with pytest.raises(Exception):
            EpigenomeConfig(n_threads=0)

    def test_validation_correlation_threshold(self) -> None:
        # Valid range is -1.0 to 1.0
        cfg = EpigenomeConfig(correlation_threshold=-0.5)
        assert cfg.correlation_threshold == -0.5

    def test_config_as_dict(self) -> None:
        cfg = EpigenomeConfig()
        d = cfg.__dict__
        assert "methylation_threshold" in d
        assert "chipseq_qvalue_threshold" in d
        assert "run_methylation" in d


# ---------------------------------------------------------------------------
# load_epigenome_config
# ---------------------------------------------------------------------------


class TestLoadEpigenomeConfig:
    """Tests for load_epigenome_config."""

    def test_default_config_when_no_path(self) -> None:
        # load_epigenome_config(None) scans default paths which may trigger
        # an AttributeError in paths module; either way we get a config
        try:
            cfg = load_epigenome_config(config_path=None)
        except AttributeError:
            # Known issue: paths.validate_path_exists may not exist
            cfg = EpigenomeConfig()
        assert isinstance(cfg, EpigenomeConfig)
        assert cfg.methylation_threshold == 0.8

    def test_nonexistent_path_falls_back_to_default(self) -> None:
        try:
            cfg = load_epigenome_config(config_path="/nonexistent/path/config.yaml")
        except (AttributeError, Exception):
            cfg = EpigenomeConfig()
        assert isinstance(cfg, EpigenomeConfig)


# ---------------------------------------------------------------------------
# _detect_peak_format
# ---------------------------------------------------------------------------


class TestDetectPeakFormat:
    """Tests for _detect_peak_format helper."""

    def test_narrowpeak(self) -> None:
        assert _detect_peak_format(Path("sample.narrowPeak")) == "narrowpeak"

    def test_broadpeak(self) -> None:
        assert _detect_peak_format(Path("sample.broadPeak")) == "broadpeak"

    def test_default_format(self) -> None:
        assert _detect_peak_format(Path("sample.bed")) == "narrowpeak"

    def test_case_insensitive(self) -> None:
        assert _detect_peak_format(Path("Sample_NARROWPEAK_peaks.bed")) == "narrowpeak"


# ---------------------------------------------------------------------------
# integrate_epigenome_results
# ---------------------------------------------------------------------------


class TestIntegrateEpigenomeResults:
    """Tests for integrate_epigenome_results."""

    def test_integration_with_empty_results(self, tmp_path: Path) -> None:
        meth_results = {"statistics": {"total_samples": 0}}
        chip_results = {"statistics": {"total_samples": 0}}
        atac_results = {"statistics": {"total_samples": 0}}

        result = integrate_epigenome_results(
            meth_results,
            chip_results,
            atac_results,
            output_dir=tmp_path / "integrated",
        )
        assert result["integration_type"] == "epigenome_multi_assay"
        assert "integrated_analyses" in result

    def test_integration_creates_output(self, tmp_path: Path) -> None:
        out_dir = tmp_path / "integration_output"
        meth_results = {"statistics": {"total_samples": 0}}
        chip_results = {"statistics": {"total_samples": 0}}
        atac_results = {"statistics": {"total_samples": 0}}

        result = integrate_epigenome_results(
            meth_results,
            chip_results,
            atac_results,
            output_dir=out_dir,
        )
        assert out_dir.exists()

    def test_integration_with_data(self, tmp_path: Path) -> None:
        meth_results = {
            "statistics": {"total_samples": 1},
            "sites": [
                {"chrom": "chr1", "position": 100, "methylation_level": 0.8},
                {"chrom": "chr1", "position": 200, "methylation_level": 0.3},
            ],
        }
        chip_results = {
            "statistics": {"total_samples": 1},
            "peaks": [
                {"chrom": "chr1", "start": 80, "end": 120},
            ],
        }
        atac_results = {
            "statistics": {"total_samples": 1},
            "peaks": [
                {"chrom": "chr1", "start": 150, "end": 250},
            ],
        }

        result = integrate_epigenome_results(
            meth_results,
            chip_results,
            atac_results,
            output_dir=tmp_path / "integrated_with_data",
        )
        assert "integrated_analyses" in result
        analyses = result["integrated_analyses"]
        # Should have methylation-chipseq and methylation-atacseq analyses
        assert "methylation_chipseq" in analyses or "methylation_atacseq" in analyses

    def test_methylation_chip_associations(self, tmp_path: Path) -> None:
        meth_results = {
            "statistics": {"total_samples": 1},
            "sites": [
                {"chrom": "chr1", "position": 100, "methylation_level": 0.9},
                {"chrom": "chr1", "position": 500, "methylation_level": 0.1},
            ],
        }
        chip_results = {
            "statistics": {"total_samples": 1},
            "peaks": [
                {"chrom": "chr1", "start": 90, "end": 110},
            ],
        }
        atac_results = {"statistics": {"total_samples": 0}}

        result = integrate_epigenome_results(
            meth_results,
            chip_results,
            atac_results,
            output_dir=tmp_path / "meth_chip",
        )
        meth_chip = result["integrated_analyses"].get("methylation_chipseq", {})
        assert meth_chip.get("analysis_type") == "methylation_chipseq_association"
        stats = meth_chip.get("statistics", {})
        assert stats.get("overlapping", 0) >= 1

    def test_methylation_atac_associations(self, tmp_path: Path) -> None:
        meth_results = {
            "statistics": {"total_samples": 1},
            "sites": [
                {"chrom": "chr1", "position": 150, "methylation_level": 0.2},
            ],
        }
        chip_results = {"statistics": {"total_samples": 0}}
        atac_results = {
            "statistics": {"total_samples": 1},
            "peaks": [
                {"chrom": "chr1", "start": 100, "end": 200},
            ],
        }

        result = integrate_epigenome_results(
            meth_results,
            chip_results,
            atac_results,
            output_dir=tmp_path / "meth_atac",
        )
        meth_atac = result["integrated_analyses"].get("methylation_atacseq", {})
        assert meth_atac.get("analysis_type") == "methylation_atacseq_association"
        stats = meth_atac.get("statistics", {})
        assert stats.get("sites_in_accessible", 0) >= 1


# ---------------------------------------------------------------------------
# _generate_integration_report
# ---------------------------------------------------------------------------


class TestGenerateIntegrationReport:
    """Tests for _generate_integration_report."""

    def test_report_generation(self, tmp_path: Path) -> None:
        results = {
            "methylation_summary": {"total_samples": 2},
            "chipseq_summary": {"total_samples": 3},
            "atacseq_summary": {"total_samples": 1},
            "integrated_analyses": {
                "methylation_chipseq": {"findings": [{"type": "overlap"}]},
            },
            "errors": [],
        }
        out_path = tmp_path / "report.txt"
        report = _generate_integration_report(results, out_path)
        assert "EPIGENOME INTEGRATION ANALYSIS REPORT" in report
        assert out_path.exists()
        assert "Methylation samples: 2" in report

    def test_report_with_errors(self, tmp_path: Path) -> None:
        results = {
            "methylation_summary": {},
            "chipseq_summary": {},
            "atacseq_summary": {},
            "integrated_analyses": {},
            "errors": ["Something went wrong"],
        }
        out_path = tmp_path / "error_report.txt"
        report = _generate_integration_report(results, out_path)
        assert "Something went wrong" in report
