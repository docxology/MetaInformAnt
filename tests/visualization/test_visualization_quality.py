"""Tests for quality control visualization functions."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from metainformant.core.utils.errors import ValidationError
from metainformant.visualization.analysis.quality import (
    plot_adapter_content,
    plot_batch_effects_qc,
    plot_coverage_uniformity,
    plot_data_integrity_metrics,
    plot_error_profiles,
    plot_gc_distribution,
    plot_kmer_profiles,
    plot_length_distribution,
    plot_multiomics_quality_overview,
    plot_overrepresented_sequences,
    plot_per_base_quality_boxplot,
    plot_protein_structure_quality,
    plot_quality_metrics,
    plot_sequence_duplication_levels,
    plot_singlecell_qc_metrics,
    plot_vcf_quality_metrics,
)


class TestPlotQualityMetrics:
    """Test plot_quality_metrics function."""

    def test_basic_quality_metrics_plot(self):
        """Test basic quality metrics plot creation."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # Create sample QC data
        qc_data = {
            "per_base_quality": {
                "positions": list(range(100)),
                "mean_qualities": np.random.uniform(20, 40, 100),
                "read_counts_at_position": np.random.randint(1000, 2000, 100),
            },
            "gc_content_distribution": {
                "bins": np.linspace(0, 100, 21),
                "counts": np.random.poisson(50, 20),
                "mean_gc_content": 45.0,
                "median_gc_content": 44.5,
            },
            "sequence_length_distribution": {"lengths": list(range(50, 151)), "counts": np.random.poisson(100, 101)},
            "basic_statistics": {
                "num_reads": 10000,
                "total_bases": 1500000,
                "min_length": 50,
                "max_length": 150,
                "mean_length": 150.0,
            },
        }

        ax = plot_quality_metrics(qc_data)
        assert ax is not None
        plt.close("all")

    def test_quality_metrics_partial_data(self):
        """Test quality metrics plot with partial data."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # Only some QC metrics
        qc_data = {
            "per_base_quality": {"positions": list(range(50)), "mean_qualities": np.random.uniform(25, 35, 50)},
            "basic_statistics": {"num_reads": 5000, "total_bases": 750000, "mean_length": 150.0},
        }

        ax = plot_quality_metrics(qc_data)
        assert ax is not None
        plt.close("all")

    def test_quality_metrics_empty_data(self):
        """Test quality metrics plot with empty data."""
        qc_data = {}

        with pytest.raises(ValueError, match="cannot be empty"):
            plot_quality_metrics(qc_data)


class TestPlotAdapterContent:
    """Test plot_adapter_content function."""

    def test_basic_adapter_content_plot(self):
        """Test basic adapter content plot creation."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # Create sample adapter data
        adapter_data = {
            "Illumina_Universal_Adapter": [0.5, 1.2, 0.8, 2.1, 1.5],
            "Illumina_Small_RNA_3p_Adapter": [0.2, 0.5, 0.3, 0.8, 0.4],
            "TruSeq_Adapter_Index_1": [0.1, 0.3, 0.2, 0.4, 0.2],
        }

        ax = plot_adapter_content(adapter_data)
        assert ax is not None
        plt.close("all")

    def test_adapter_content_with_output_path(self, tmp_path: Path):
        """Test adapter content plot with output path."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        adapter_data = {"Adapter1": [1.0, 1.5], "Adapter2": [0.5, 0.8]}
        output_path = tmp_path / "adapter_content.png"

        ax = plot_adapter_content(adapter_data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close("all")

    def test_adapter_content_empty_data(self):
        """Test adapter content plot with empty data."""
        adapter_data = {}

        with pytest.raises(ValueError, match="cannot be empty"):
            plot_adapter_content(adapter_data)


class TestPlotGcDistribution:
    """Test plot_gc_distribution function."""

    def test_basic_gc_distribution_plot(self):
        """Test basic GC distribution plot creation."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # Create sample GC content data
        gc_data = np.random.normal(45, 10, 1000)
        gc_data = np.clip(gc_data, 0, 100)  # Ensure valid range

        ax = plot_gc_distribution(gc_data)
        assert ax is not None
        plt.close("all")

    def test_gc_distribution_with_output_path(self, tmp_path: Path):
        """Test GC distribution plot with output path."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        gc_data = [40, 45, 50, 35, 55, 42, 48, 38, 52, 46]
        output_path = tmp_path / "gc_dist.png"

        ax = plot_gc_distribution(gc_data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close("all")

    def test_gc_distribution_empty_data(self):
        """Test GC distribution plot with empty data."""
        gc_data = []

        with pytest.raises(ValueError, match="cannot be empty"):
            plot_gc_distribution(gc_data)


class TestPlotLengthDistribution:
    """Test plot_length_distribution function."""

    def test_basic_length_distribution_plot(self):
        """Test basic length distribution plot creation."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # Create sample length data
        length_data = np.random.normal(150, 20, 1000).astype(int)
        length_data = np.clip(length_data, 50, 250)  # Ensure reasonable range

        ax = plot_length_distribution(length_data)
        assert ax is not None
        plt.close("all")

    def test_length_distribution_with_output_path(self, tmp_path: Path):
        """Test length distribution plot with output path."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        length_data = [100, 150, 120, 140, 160, 130, 145, 155, 125, 135]
        output_path = tmp_path / "length_dist.png"

        ax = plot_length_distribution(length_data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close("all")

    def test_length_distribution_empty_data(self):
        """Test length distribution plot with empty data."""
        length_data = []

        with pytest.raises(ValueError, match="cannot be empty"):
            plot_length_distribution(length_data)


class TestPlotPerBaseQualityBoxplot:
    """Test plot_per_base_quality_boxplot function."""

    def test_basic_boxplot_mixed_keys(self):
        """Test boxplot with range-string and integer position keys."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        rng = np.random.default_rng(42)
        per_base_qualities = {
            "1-10": rng.uniform(20, 40, 50),
            "11-20": rng.uniform(25, 40, 50),
            21: rng.uniform(25, 38, 50),
        }

        ax = plot_per_base_quality_boxplot(per_base_qualities)
        assert ax is not None
        assert len(ax.lines) > 0  # boxplot whiskers/medians drawn
        assert ax.get_title() == "Per-Base Quality Scores"
        plt.close("all")

    def test_boxplot_with_output_path(self, tmp_path: Path):
        """Test boxplot saved to output path."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        rng = np.random.default_rng(7)
        per_base_qualities = {"1-5": rng.uniform(30, 40, 30), "6-10": rng.uniform(28, 38, 30)}
        output_path = tmp_path / "per_base_quality.png"

        ax = plot_per_base_quality_boxplot(per_base_qualities, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close("all")


class TestPlotSequenceDuplicationLevels:
    """Test plot_sequence_duplication_levels function."""

    def test_basic_duplication_levels(self):
        """Test duplication level bar chart creation."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        duplication_levels = {"1": 60.0, "2": 25.0, "3": 10.0, "4": 5.0}

        ax = plot_sequence_duplication_levels(duplication_levels)
        assert ax is not None
        assert len(ax.patches) == 4
        assert ax.get_title() == "Sequence Duplication Levels"
        plt.close("all")

    def test_duplication_levels_with_output_path(self, tmp_path: Path):
        """Test duplication levels plot saved to output path."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        duplication_levels = {"1": 70.0, "2": 20.0, "3": 10.0}
        output_path = tmp_path / "duplication.png"

        ax = plot_sequence_duplication_levels(duplication_levels, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close("all")


class TestPlotOverrepresentedSequences:
    """Test plot_overrepresented_sequences function."""

    def test_basic_overrepresented_sequences(self, tmp_path: Path):
        """Test horizontal bar chart for overrepresented sequences."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        seqs = [
            {"sequence": "ACGTACGTACGTACGTACGTACGT", "count": 1200, "percentage": 2.5},
            {"sequence": "TTGGCATTGCATGGACTTGA", "count": 800, "percentage": 1.8},
        ]
        output_path = tmp_path / "overrepresented.png"

        ax = plot_overrepresented_sequences(seqs, output_path=output_path)
        assert ax is not None
        assert len(ax.patches) == 2
        assert output_path.exists()
        plt.close("all")

    def test_empty_list_renders_placeholder(self):
        """Test empty sequence list renders a placeholder without raising."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        ax = plot_overrepresented_sequences([])
        assert ax is not None
        assert ax.get_title() == "Overrepresented Sequences"
        plt.close("all")


class TestPlotKmerProfiles:
    """Test plot_kmer_profiles function."""

    def test_basic_kmer_profiles(self):
        """Test top-N k-mer bar chart creation."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        rng = np.random.default_rng(3)
        kmer_counts = {
            f"{'ACGT'[i % 4]}{'TGCA'[i % 4]}{'GCAT'[i % 4]}{i:02d}": int(rng.integers(10, 1000)) for i in range(30)
        }

        ax = plot_kmer_profiles(kmer_counts, top_n=10)
        assert ax is not None
        assert len(ax.patches) == 10
        assert "Top 10" in ax.get_title()
        plt.close("all")

    def test_kmer_profiles_with_output_path(self, tmp_path: Path):
        """Test k-mer profile plot saved to output path."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        kmer_counts = {"AAAA": 500, "AAAC": 400, "AAGG": 300}
        output_path = tmp_path / "kmers.png"

        ax = plot_kmer_profiles(kmer_counts, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close("all")


class TestPlotVcfQualityMetrics:
    """Test plot_vcf_quality_metrics function."""

    def test_full_vcf_metrics(self, tmp_path: Path):
        """Test all six VCF QC panels render and save."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        rng = np.random.default_rng(11)
        vcf_qc_data = {
            "qual_distribution": {"qualities": [10, 20, 30, 40], "counts": [5, 12, 8, 2]},
            "depth_distribution": {"depths": [5, 10, 20, 30], "counts": [8, 12, 6, 1]},
            "allele_frequencies": rng.uniform(0.0, 0.5, 100),
            "variant_types": {"SNP": 80, "INDEL": 20},
            "titv_by_qual": {"qualities": [20, 30, 40], "titv_ratios": [1.9, 2.1, 2.2]},
            "summary_stats": {"n_variants": 1000, "mean_depth": 25.0},
        }
        output_path = tmp_path / "vcf_qc.png"

        ax = plot_vcf_quality_metrics(vcf_qc_data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close("all")

    def test_partial_vcf_metrics(self):
        """Test VCF QC plot with only some sections present."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        ax = plot_vcf_quality_metrics({"allele_frequencies": [0.1, 0.2, 0.3]})
        assert ax is not None
        plt.close("all")


class TestPlotSinglecellQcMetrics:
    """Test plot_singlecell_qc_metrics function."""

    def test_basic_singlecell_metrics(self, tmp_path: Path):
        """Test single-cell QC histograms and correlation panel."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        rng = np.random.default_rng(5)
        qc_metrics = {
            "n_counts": rng.exponential(5000, 200),
            "n_genes": rng.exponential(3000, 200),
            "percent_mito": rng.uniform(0, 20, 200),
            "doublet_score": rng.uniform(0, 1, 200),
        }
        output_path = tmp_path / "singlecell_qc.png"

        ax = plot_singlecell_qc_metrics(qc_metrics, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close("all")

    def test_minimal_singlecell_metrics(self):
        """Test single-cell QC plot with a single metric (no correlation panel)."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        rng = np.random.default_rng(6)
        ax = plot_singlecell_qc_metrics({"n_counts": rng.exponential(5000, 100)})
        assert ax is not None
        plt.close("all")


class TestPlotProteinStructureQuality:
    """Test plot_protein_structure_quality function."""

    def test_full_structure_metrics(self, tmp_path: Path):
        """Test protein structure panels with dict-valued scores."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        rng = np.random.default_rng(13)
        structure_quality = {
            "b_factors": rng.normal(30, 5, 100),
            "ramachandran_stats": {"favored": 0.9, "allowed": 0.08, "outliers": 0.02},
            "clash_score": {"severe": 2, "mild": 10},
            "overall_quality": {"ramachandran": 0.95, "clashscore": 0.8},
        }
        output_path = tmp_path / "protein_qc.png"

        ax = plot_protein_structure_quality(structure_quality, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close("all")

    def test_scalar_scores(self):
        """Test protein structure panels with scalar clash and quality scores."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        structure_quality = {"clash_score": 4.5, "overall_quality": 0.88}

        ax = plot_protein_structure_quality(structure_quality)
        assert ax is not None
        plt.close("all")


class TestPlotMultiomicsQualityOverview:
    """Test plot_multiomics_quality_overview function."""

    def test_basic_overview_with_output_path(self, tmp_path: Path):
        """Test radar overview across omics layers, including an empty report."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        quality_reports = {
            "Transcriptomics": {"overall_quality": 0.92},
            "Proteomics": {"quality_score": 0.85},
            "Metabolomics": {},
        }
        output_path = tmp_path / "multiomics_qc.png"

        ax = plot_multiomics_quality_overview(quality_reports, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close("all")


class TestPlotCoverageUniformity:
    """Test plot_coverage_uniformity function."""

    def test_basic_coverage_with_output_path(self, tmp_path: Path):
        """Test coverage plot with default x axis."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        rng = np.random.default_rng(17)
        coverage = rng.poisson(30, 200).astype(float)
        output_path = tmp_path / "coverage.png"

        ax = plot_coverage_uniformity(coverage, output_path=output_path)
        assert ax is not None
        assert ax.get_xlabel() == "Position"
        assert output_path.exists()
        plt.close("all")

    def test_coverage_with_genomic_positions(self):
        """Test coverage plot with explicit positions labels the x axis."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        rng = np.random.default_rng(19)
        coverage = rng.poisson(25, 100).astype(float)
        positions = np.arange(1, 101)

        ax = plot_coverage_uniformity(coverage, positions=positions)
        assert ax is not None
        assert ax.get_xlabel() == "Genomic Position"
        plt.close("all")

    def test_rejects_non_array_coverage(self):
        """Test that a plain list is rejected by type validation."""
        with pytest.raises(ValidationError):
            plot_coverage_uniformity([10.0, 20.0, 30.0])


class TestPlotErrorProfiles:
    """Test plot_error_profiles function."""

    def test_basic_error_profiles_with_output_path(self, tmp_path: Path):
        """Test error profile lines render on a log axis and save."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        rng = np.random.default_rng(23)
        profiles = {
            "mismatch": np.clip(rng.exponential(0.001, 100), 1e-6, None),
            "indel": np.clip(rng.exponential(0.0005, 100), 1e-6, None),
        }
        output_path = tmp_path / "error_profiles.png"

        ax = plot_error_profiles(profiles, output_path=output_path)
        assert ax is not None
        assert ax.get_yscale() == "log"
        assert len(ax.lines) == 2
        assert output_path.exists()
        plt.close("all")


class TestPlotBatchEffectsQc:
    """Test plot_batch_effects_qc function."""

    def test_full_batch_panels(self, tmp_path: Path):
        """Test all six batch-effect panels render and save."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        rng = np.random.default_rng(29)
        batch_qc_data = {
            "batch_sizes": {"batchA": 50, "batchB": 60},
            "batch_pca": {
                "pc1": rng.normal(0, 1, 20),
                "pc2": rng.normal(0, 1, 20),
                "batches": ["batchA"] * 10 + ["batchB"] * 10,
            },
            "silhouette_scores": {"batchA": 0.45, "batchB": 0.52},
            "batch_de_stats": {"batches": ["batchA_vs_batchB"], "n_de_genes": [120]},
            "batch_variance": {"PC1": 0.5, "PC2": 0.3},
            "correction_metrics": {"kBET": 0.1},
        }
        output_path = tmp_path / "batch_qc.png"

        ax = plot_batch_effects_qc(batch_qc_data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close("all")


class TestPlotDataIntegrityMetrics:
    """Test plot_data_integrity_metrics function."""

    def test_numeric_metrics_are_plotted(self):
        """Test that only numeric integrity metrics become bars."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        integrity_metrics = {
            "missing_rate": 0.05,
            "error_rate": 0.005,
            "completeness": 0.95,
            "valid_records": 0.99,
            "dataset_name": "should_be_ignored",
        }

        ax = plot_data_integrity_metrics(integrity_metrics)
        assert ax is not None
        assert len(ax.patches) == 4  # non-numeric value skipped
        plt.close("all")

    def test_integrity_metrics_with_output_path(self, tmp_path: Path):
        """Test integrity metrics plot saved to output path."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        integrity_metrics = {"missing_rate": 0.02, "completeness": 0.98}
        output_path = tmp_path / "integrity.png"

        ax = plot_data_integrity_metrics(integrity_metrics, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close("all")


class TestQualityModuleAggregation:
    """Test the quality module's aggregation of the focused submodules."""

    def test_quality_module_reexports_sibling_implementations(self):
        """Test that quality re-exports the focused submodule functions."""
        from metainformant.visualization.analysis import (
            quality,
            quality_assessment,
            quality_omics,
            quality_sequencing,
        )

        assert quality.plot_quality_metrics is quality_sequencing.plot_quality_metrics
        assert quality.plot_gc_distribution is quality_sequencing.plot_gc_distribution
        assert quality.plot_kmer_profiles is quality_sequencing.plot_kmer_profiles
        assert quality.plot_vcf_quality_metrics is quality_omics.plot_vcf_quality_metrics
        assert quality.plot_singlecell_qc_metrics is quality_omics.plot_singlecell_qc_metrics
        assert quality.plot_coverage_uniformity is quality_assessment.plot_coverage_uniformity
        assert quality.plot_batch_effects_qc is quality_assessment.plot_batch_effects_qc

    def test_gc_distribution_explicit_bins(self):
        """Regression: explicit bins kwarg must not collide with the internal default."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        ax = plot_gc_distribution([40, 45, 50, 35, 55], bins=5)
        assert ax is not None
        assert len(ax.patches) == 5
        plt.close("all")

    def test_length_distribution_explicit_bins(self):
        """Regression: explicit bins kwarg must not collide with the internal default."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        ax = plot_length_distribution([100, 120, 140, 160, 180], bins=5)
        assert ax is not None
        assert len(ax.patches) == 5
        plt.close("all")
