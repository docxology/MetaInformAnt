"""Tests for GWAS phenotype visualization module.

Tests cover phenotype distribution, correlation matrix, genotype-phenotype
boxplots, top hits auto-generation, and phenotype-PCA correlation panels.
All tests use synthetic data and verify output files exist with success status.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pytest

from metainformant.gwas.visualization.interactive.phenotype import (
    genotype_phenotype_boxplot,
    phenotype_correlation_matrix,
    phenotype_distribution,
    phenotype_pca_correlation,
    top_hits_genotype_phenotype,
)

# ---------------------------------------------------------------------------
# Fixtures for synthetic data
# ---------------------------------------------------------------------------


@pytest.fixture
def rng() -> np.random.Generator:
    """Seeded random number generator for reproducibility."""
    return np.random.default_rng(12345)


@pytest.fixture
def phenotype_values(rng: np.random.Generator) -> List[float]:
    """200 normally distributed phenotype values."""
    return rng.normal(loc=50.0, scale=10.0, size=200).tolist()


@pytest.fixture
def sample_ids() -> List[str]:
    """200 sample identifiers."""
    return [f"sample_{i}" for i in range(200)]


@pytest.fixture
def population_map(sample_ids: List[str]) -> Dict[str, str]:
    """Map first 100 samples to PopA and the rest to PopB."""
    mapping: Dict[str, str] = {}
    for i, sid in enumerate(sample_ids):
        mapping[sid] = "PopA" if i < 100 else "PopB"
    return mapping


@pytest.fixture
def three_traits(rng: np.random.Generator) -> Dict[str, List[float]]:
    """Three correlated traits of length 150."""
    base = rng.normal(size=150)
    return {
        "Height": (base * 5.0 + 170.0).tolist(),
        "Weight": (base * 8.0 + rng.normal(size=150) * 3.0 + 75.0).tolist(),
        "BMI": (base * 1.5 + rng.normal(size=150) * 1.0 + 24.0).tolist(),
    }


@pytest.fixture
def genotype_data(rng: np.random.Generator) -> tuple[List[int], List[float]]:
    """Genotypes (0,1,2) and phenotypes for 300 samples."""
    genotypes = rng.choice([0, 1, 2], size=300, p=[0.5, 0.35, 0.15]).tolist()
    phenotypes = []
    for g in genotypes:
        phenotypes.append(float(rng.normal(loc=50.0 + g * 3.0, scale=8.0)))
    return genotypes, phenotypes


@pytest.fixture
def assoc_results_10(rng: np.random.Generator) -> tuple[List[Dict[str, Any]], List[List[int]], List[float]]:
    """10 synthetic association results with genotype matrix and phenotypes."""
    n_variants = 10
    n_samples = 100
    phenotypes = rng.normal(loc=50.0, scale=10.0, size=n_samples).tolist()

    assoc_results: List[Dict[str, Any]] = []
    genotype_matrix: List[List[int]] = []

    for i in range(n_variants):
        p_val = float(rng.uniform(1e-10, 0.5))
        assoc_results.append(
            {
                "variant_id": f"rs{100000 + i}",
                "p_value": p_val,
                "beta": float(rng.normal(0, 0.5)),
                "chromosome": "1",
                "position": 1000000 + i * 50000,
            }
        )
        genos = rng.choice([0, 1, 2], size=n_samples, p=[0.5, 0.35, 0.15]).tolist()
        genotype_matrix.append(genos)

    return assoc_results, genotype_matrix, phenotypes


@pytest.fixture
def pca_data(rng: np.random.Generator) -> tuple[List[List[float]], List[float]]:
    """PCA components (5 PCs) and phenotype for 120 samples."""
    n_samples = 120
    pcs = rng.normal(size=(n_samples, 5)).tolist()
    phenotypes = rng.normal(loc=100, scale=15, size=n_samples).tolist()
    return pcs, phenotypes


# ---------------------------------------------------------------------------
# Test: phenotype_distribution
# ---------------------------------------------------------------------------


class TestPhenotypeDistribution:
    """Tests for phenotype_distribution function."""

    def test_basic_distribution(self, tmp_path: Path, phenotype_values: List[float]) -> None:
        """Test basic phenotype distribution without population stratification."""
        output_file = tmp_path / "pheno_dist.png"
        result = phenotype_distribution(
            phenotype_values=phenotype_values,
            trait_name="Body Weight",
            output_file=output_file,
        )

        assert result["status"] == "success"
        assert output_file.exists()
        assert result["output_path"] == str(output_file)
        assert result["n_samples"] == 200
        assert "mean" in result
        assert "median" in result
        assert "std" in result

    def test_distribution_by_population(
        self,
        tmp_path: Path,
        phenotype_values: List[float],
        sample_ids: List[str],
        population_map: Dict[str, str],
    ) -> None:
        """Test phenotype distribution stratified by population."""
        output_file = tmp_path / "pheno_dist_pop.png"
        result = phenotype_distribution(
            phenotype_values=phenotype_values,
            trait_name="Body Weight",
            output_file=output_file,
            by_population=population_map,
            sample_ids=sample_ids,
        )

        assert result["status"] == "success"
        assert output_file.exists()
        assert result["n_samples"] == 200

    def test_empty_values(self, tmp_path: Path) -> None:
        """Test with empty phenotype values returns failed status."""
        result = phenotype_distribution(
            phenotype_values=[],
            output_file=tmp_path / "empty.png",
        )
        assert result["status"] == "failed"

    def test_no_output_file(self, phenotype_values: List[float]) -> None:
        """Test without output file still computes statistics."""
        result = phenotype_distribution(phenotype_values=phenotype_values)
        assert result["status"] == "success"
        assert result["output_path"] is None


# ---------------------------------------------------------------------------
# Test: phenotype_correlation_matrix
# ---------------------------------------------------------------------------


class TestPhenotypeCorrelationMatrix:
    """Tests for phenotype_correlation_matrix function."""

    def test_three_trait_correlation(self, tmp_path: Path, three_traits: Dict[str, List[float]]) -> None:
        """Test correlation matrix with three correlated traits."""
        output_file = tmp_path / "corr_matrix.png"
        result = phenotype_correlation_matrix(
            phenotypes_dict=three_traits,
            output_file=output_file,
        )

        assert result["status"] == "success"
        assert output_file.exists()
        assert result["n_traits"] == 3
        assert set(result["trait_names"]) == {"Height", "Weight", "BMI"}
        assert "correlation_matrix" in result

        corr = result["correlation_matrix"]
        # Diagonal should be 1.0
        for trait in ["Height", "Weight", "BMI"]:
            assert abs(corr[trait][trait] - 1.0) < 1e-6

        # Height and Weight should be positively correlated (same base signal)
        assert corr["Height"]["Weight"] > 0.3

    def test_single_trait_rejected(self, tmp_path: Path) -> None:
        """Test that a single trait is rejected."""
        result = phenotype_correlation_matrix(
            phenotypes_dict={"OnlyTrait": [1.0, 2.0, 3.0]},
            output_file=tmp_path / "single.png",
        )
        assert result["status"] == "failed"

    def test_unequal_lengths_rejected(self, tmp_path: Path) -> None:
        """Test that unequal-length trait lists are rejected."""
        result = phenotype_correlation_matrix(
            phenotypes_dict={"A": [1.0, 2.0], "B": [1.0, 2.0, 3.0]},
            output_file=tmp_path / "unequal.png",
        )
        assert result["status"] == "failed"


# ---------------------------------------------------------------------------
# Test: genotype_phenotype_boxplot
# ---------------------------------------------------------------------------


class TestGenotypePhenotypeBoxplot:
    """Tests for genotype_phenotype_boxplot function."""

    def test_three_genotype_classes(
        self,
        tmp_path: Path,
        genotype_data: tuple[List[int], List[float]],
    ) -> None:
        """Test boxplot with all three genotype classes present."""
        genotypes, phenotypes = genotype_data
        output_file = tmp_path / "gp_boxplot.png"
        result = genotype_phenotype_boxplot(
            genotypes=genotypes,
            phenotypes=phenotypes,
            variant_id="rs12345",
            output_file=output_file,
        )

        assert result["status"] == "success"
        assert output_file.exists()
        assert result["variant_id"] == "rs12345"
        assert result["n_samples"] == 300
        assert result["n_groups"] == 3
        assert "group_stats" in result
        assert result["f_statistic"] is not None
        assert result["p_value"] is not None
        # With 3.0 effect per allele, F-stat should be substantial
        assert result["f_statistic"] > 1.0

    def test_length_mismatch(self, tmp_path: Path) -> None:
        """Test that mismatched lengths return failed status."""
        result = genotype_phenotype_boxplot(
            genotypes=[0, 1, 2],
            phenotypes=[1.0, 2.0],
            output_file=tmp_path / "mismatch.png",
        )
        assert result["status"] == "failed"

    def test_empty_data(self, tmp_path: Path) -> None:
        """Test that empty data returns failed status."""
        result = genotype_phenotype_boxplot(
            genotypes=[],
            phenotypes=[],
            output_file=tmp_path / "empty.png",
        )
        assert result["status"] == "failed"


# ---------------------------------------------------------------------------
# Test: top_hits_genotype_phenotype
# ---------------------------------------------------------------------------


class TestTopHitsGenotypePhenotype:
    """Tests for top_hits_genotype_phenotype function."""

    def test_top_5_from_10(
        self,
        tmp_path: Path,
        assoc_results_10: tuple[List[Dict[str, Any]], List[List[int]], List[float]],
    ) -> None:
        """Test generating top 5 boxplots from 10 association results."""
        assoc_results, genotype_matrix, phenotypes = assoc_results_10
        output_dir = tmp_path / "top_hits"

        result = top_hits_genotype_phenotype(
            assoc_results=assoc_results,
            genotype_matrix=genotype_matrix,
            phenotypes=phenotypes,
            output_dir=output_dir,
            n_top=5,
        )

        assert result["status"] == "success"
        assert result["n_plotted"] == 5
        assert len(result["generated_files"]) == 5

        for fpath in result["generated_files"]:
            assert Path(fpath).exists()

        # Verify summaries are sorted by p-value (rank order)
        for i, summary in enumerate(result["variant_summaries"]):
            assert summary["rank"] == i + 1

    def test_top_more_than_available(
        self,
        tmp_path: Path,
        assoc_results_10: tuple[List[Dict[str, Any]], List[List[int]], List[float]],
    ) -> None:
        """Test requesting more top hits than available variants."""
        assoc_results, genotype_matrix, phenotypes = assoc_results_10
        output_dir = tmp_path / "top_hits_all"

        result = top_hits_genotype_phenotype(
            assoc_results=assoc_results,
            genotype_matrix=genotype_matrix,
            phenotypes=phenotypes,
            output_dir=output_dir,
            n_top=20,
        )

        assert result["status"] == "success"
        assert result["n_plotted"] == 10

    def test_empty_results(self, tmp_path: Path) -> None:
        """Test with empty association results."""
        result = top_hits_genotype_phenotype(
            assoc_results=[],
            genotype_matrix=[],
            phenotypes=[],
            output_dir=tmp_path / "empty",
        )
        assert result["status"] == "failed"


# ---------------------------------------------------------------------------
# Test: phenotype_pca_correlation
# ---------------------------------------------------------------------------


class TestPhenotypePCACorrelation:
    """Tests for phenotype_pca_correlation function."""

    def test_pca_correlation_panel(
        self,
        tmp_path: Path,
        pca_data: tuple[List[List[float]], List[float]],
    ) -> None:
        """Test PCA-phenotype correlation with 5 PCs and 120 samples."""
        pcs, phenotypes = pca_data
        output_file = tmp_path / "pca_pheno.png"

        result = phenotype_pca_correlation(
            pcs=pcs,
            phenotype_values=phenotypes,
            output_file=output_file,
            trait_name="Blood Pressure",
        )

        assert result["status"] == "success"
        assert output_file.exists()
        assert result["n_samples"] == 120
        assert "r_pc1" in result
        assert "r_pc2" in result
        assert -1.0 <= result["r_pc1"] <= 1.0
        assert -1.0 <= result["r_pc2"] <= 1.0

    def test_insufficient_components(self, tmp_path: Path) -> None:
        """Test with only 1 PC fails gracefully."""
        result = phenotype_pca_correlation(
            pcs=[[1.0], [2.0], [3.0]],
            phenotype_values=[10.0, 20.0, 30.0],
            output_file=tmp_path / "bad_pcs.png",
        )
        assert result["status"] == "failed"

    def test_length_mismatch(self, tmp_path: Path) -> None:
        """Test mismatched lengths between PCs and phenotype."""
        result = phenotype_pca_correlation(
            pcs=[[1.0, 2.0], [3.0, 4.0]],
            phenotype_values=[10.0, 20.0, 30.0],
            output_file=tmp_path / "mismatch.png",
        )
        assert result["status"] == "failed"

    def test_no_output_file(
        self,
        pca_data: tuple[List[List[float]], List[float]],
    ) -> None:
        """Test without output file still computes correlations."""
        pcs, phenotypes = pca_data
        result = phenotype_pca_correlation(pcs=pcs, phenotype_values=phenotypes)
        assert result["status"] == "success"
        assert result["output_path"] is None
