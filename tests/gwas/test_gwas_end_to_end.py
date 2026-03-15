"""End-to-end integration tests for the GWAS pipeline.

Tests the complete workflow from VCF input through association testing,
summary statistics output, and annotation.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from metainformant.gwas.analysis.association import association_test_linear
from metainformant.gwas.analysis.quality import parse_vcf_full
from metainformant.gwas.data.genome import normalize_chromosome_name
from metainformant.gwas.workflow.workflow import _normalize_config, load_gwas_config, run_gwas
from tests.fixtures.gwas.generate_test_data import generate_complete_test_dataset


@pytest.fixture
def test_dataset(tmp_path: Path) -> dict:
    """Generate a complete test dataset."""
    return generate_complete_test_dataset(
        output_dir=tmp_path / "data",
        n_variants_per_chrom=10,
        n_samples=30,
        n_chroms=3,
        n_causal=2,
        effect_size=5.0,
        seed=42,
    )


class TestFullLinearGWAS:
    """Test complete linear model GWAS pipeline."""

    def test_full_linear_pipeline(self, tmp_path: Path, test_dataset: dict) -> None:
        """Run full linear GWAS and verify all steps complete."""
        config = {
            "model": "linear",
            "trait": "varroa_resistance",
        }
        results_dir = tmp_path / "results"

        result = run_gwas(
            test_dataset["vcf_path"],
            test_dataset["phenotype_path"],
            config,
            output_dir=results_dir,
        )

        assert result["status"] == "success"
        assert "parse_vcf" in result["steps_completed"]
        assert "qc_filters" in result["steps_completed"]
        assert "population_structure" in result["steps_completed"]
        assert "load_phenotypes" in result["steps_completed"]
        assert "association_testing" in result["steps_completed"]
        assert "multiple_testing_correction" in result["steps_completed"]

    def test_association_results_not_hardcoded(self, tmp_path: Path, test_dataset: dict) -> None:
        """Verify association results are NOT the old hardcoded values."""
        config = {"model": "linear", "trait": "varroa_resistance"}
        result = run_gwas(
            test_dataset["vcf_path"],
            test_dataset["phenotype_path"],
            config,
            output_dir=tmp_path / "results",
        )

        assoc_results = result["results"].get("association_results", [])
        assert len(assoc_results) > 0

        # The old hardcoded values were: beta=0.1, se=0.05, t_stat=2.0, p_value=0.05, r_squared=0.5
        for r in assoc_results:
            # At least some results should differ from the hardcoded values
            if r.get("status") == "success":
                # Not all should be exactly the same
                assert not (
                    abs(r.get("beta", 0) - 0.1) < 1e-10
                    and abs(r.get("se", 0) - 0.05) < 1e-10
                    and abs(r.get("p_value", 0) - 0.05) < 1e-10
                ), "Association results appear to be hardcoded!"

    def test_p_values_vary_across_variants(self, tmp_path: Path, test_dataset: dict) -> None:
        """Different variants should produce different p-values (not all identical)."""
        config = {"model": "linear", "trait": "varroa_resistance"}
        result = run_gwas(
            test_dataset["vcf_path"],
            test_dataset["phenotype_path"],
            config,
            output_dir=tmp_path / "results",
        )

        p_values = [
            r["p_value"] for r in result["results"].get("association_results", []) if r.get("status") == "success"
        ]
        assert len(p_values) > 1
        # Not all p-values should be identical
        assert len(set(round(p, 8) for p in p_values)) > 1


class TestFullMixedModelGWAS:
    """Test complete mixed model GWAS pipeline."""

    def test_full_mixed_pipeline(self, tmp_path: Path, test_dataset: dict) -> None:
        """Run full mixed model GWAS."""
        config = {
            "model": "mixed",
            "trait": "varroa_resistance",
        }
        results_dir = tmp_path / "results"

        result = run_gwas(
            test_dataset["vcf_path"],
            test_dataset["phenotype_path"],
            config,
            output_dir=results_dir,
        )

        assert result["status"] == "success"
        assert "association_testing" in result["steps_completed"]

        assoc_results = result["results"].get("association_results", [])
        assert len(assoc_results) > 0

        # Mixed model results should have test_type
        for r in assoc_results:
            assert r.get("test_type") == "mixed"


class TestChromosomeHandling:
    """Test chromosome name normalization in the pipeline."""

    def test_ncbi_accession_names(self, tmp_path: Path) -> None:
        """NCBI accession names should be normalized."""
        dataset = generate_complete_test_dataset(
            output_dir=tmp_path / "ncbi_data",
            n_variants_per_chrom=5,
            n_samples=20,
            n_chroms=3,
            use_ncbi_names=True,
            seed=123,
        )

        vcf_data = parse_vcf_full(dataset["vcf_path"])
        variants = vcf_data.get("variants", [])
        assert len(variants) > 0

        # Verify NCBI names are present in variants
        chroms = set(v["chrom"] for v in variants)
        has_ncbi = any(c.startswith("NC_") for c in chroms)
        assert has_ncbi

        # Verify normalization works for all
        for chrom in chroms:
            num = normalize_chromosome_name(chrom)
            assert num is not None
            assert 1 <= num <= 16


class TestSummaryStatsWritten:
    """Test that summary statistics files are produced."""

    def test_summary_stats_output(self, tmp_path: Path, test_dataset: dict) -> None:
        """Summary statistics should be written to output directory."""
        config = {
            "model": "linear",
            "trait": "varroa_resistance",
        }
        results_dir = tmp_path / "results"

        result = run_gwas(
            test_dataset["vcf_path"],
            test_dataset["phenotype_path"],
            config,
            output_dir=results_dir,
        )

        assert result["status"] == "success"

        # Check summary stats file was written
        stats_path = results_dir / "summary_statistics.tsv"
        assert stats_path.exists()
        lines = stats_path.read_text().strip().split("\n")
        assert len(lines) > 1  # header + at least one variant

        # Check significant hits file
        sig_path = results_dir / "significant_hits.tsv"
        assert sig_path.exists()

        # Check JSON summary
        summary_path = results_dir / "results_summary.json"
        assert summary_path.exists()


class TestConfigDrivenWorkflow:
    """Test YAML config-driven workflow."""

    def test_nested_config_works(self, tmp_path: Path, test_dataset: dict) -> None:
        """Nested YAML config should be normalized and work."""
        nested_config = {
            "work_dir": str(tmp_path / "work"),
            "variants": {
                "vcf_files": [test_dataset["vcf_path"]],
            },
            "samples": {
                "phenotype_file": test_dataset["phenotype_path"],
            },
            "qc": {
                "min_maf": 0.01,
                "max_missing": 0.1,
                "hwe_pval": 1e-6,
            },
            "association": {
                "model": "linear",
                "trait": "varroa_resistance",
            },
        }

        # Normalize should extract flat keys
        flat = _normalize_config(nested_config)
        assert flat["vcf_path"] == test_dataset["vcf_path"]
        assert flat["phenotype_path"] == test_dataset["phenotype_path"]
        assert flat["model"] == "linear"
        assert flat["trait"] == "varroa_resistance"
        assert "quality_control" in flat

    def test_flat_config_unchanged(self) -> None:
        """Flat config should pass through normalization unchanged."""
        flat_config = {
            "work_dir": "/tmp/work",
            "vcf_path": "/tmp/test.vcf",
            "phenotype_path": "/tmp/pheno.tsv",
            "model": "linear",
            "trait": "height",
        }
        normalized = _normalize_config(flat_config)
        assert normalized["vcf_path"] == "/tmp/test.vcf"
        assert normalized["model"] == "linear"


class TestLDPruningIntegration:
    """Test LD pruning integration in workflow."""

    def test_ld_pruning_enabled(self, tmp_path: Path, test_dataset: dict) -> None:
        """LD pruning should run when enabled in config."""
        config = {
            "model": "linear",
            "trait": "varroa_resistance",
            "ld_pruning": {
                "enabled": True,
                "window_size": 10,
                "step_size": 5,
                "r2_threshold": 0.2,
            },
        }

        result = run_gwas(
            test_dataset["vcf_path"],
            test_dataset["phenotype_path"],
            config,
            output_dir=tmp_path / "results",
        )

        assert result["status"] == "success"
        assert "ld_pruning" in result["steps_completed"]
        ld_info = result["results"].get("ld_pruning", {})
        assert ld_info.get("variants_before", 0) > 0


class TestMultiCovariateRegression:
    """Test that multi-covariate regression produces real results."""

    def test_covariates_change_results(self) -> None:
        """Adding covariates should change regression results."""
        np.random.seed(42)
        n = 50
        genotypes = list(np.random.choice([0, 1, 2], size=n))
        phenotypes = [float(g * 0.5 + np.random.randn()) for g in genotypes]

        # Without covariates
        result_no_cov = association_test_linear(genotypes, phenotypes)

        # With covariates
        cov1 = [float(np.random.randn()) for _ in range(n)]
        cov2 = [float(np.random.randn()) for _ in range(n)]
        result_with_cov = association_test_linear(genotypes, phenotypes, covariates=[cov1, cov2])

        assert result_no_cov["status"] == "success"
        assert result_with_cov["status"] == "success"

        # Results should differ (not hardcoded)
        assert result_no_cov["beta"] != result_with_cov["beta"] or result_no_cov["se"] != result_with_cov["se"]

        # Both should have valid (non-zero) standard errors
        assert result_no_cov["se"] > 0
        assert result_with_cov["se"] > 0

    def test_covariates_not_hardcoded(self) -> None:
        """Multi-covariate regression should NOT return hardcoded 0.1/0.05/2.0/0.05/0.5."""
        genotypes = [0, 1, 2, 0, 1, 2, 0, 1, 2, 0]
        phenotypes = [10.0, 11.0, 12.0, 10.0, 11.0, 12.0, 10.0, 11.0, 12.0, 10.0]
        covariates = [[30, 40, 50, 30, 40, 50, 30, 40, 50, 30]]

        result = association_test_linear(genotypes, phenotypes, covariates=covariates)
        assert result["status"] == "success"

        # The old hardcoded return was (0.1, 0.05, 2.0, 0.05, 0.5)
        assert not (
            abs(result["beta"] - 0.1) < 1e-10 and abs(result["se"] - 0.05) < 1e-10
        ), "Multi-covariate regression still returns hardcoded values!"
