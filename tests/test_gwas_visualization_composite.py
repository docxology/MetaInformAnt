"""Tests for GWAS composite multi-panel visualization functions."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pytest

from metainformant.gwas.visualization.interactive.composite import (
    gwas_summary_panel,
    population_structure_panel,
    top_hit_detail_panel,
)

# ---------------------------------------------------------------------------
# Helpers for synthetic data generation
# ---------------------------------------------------------------------------


def _make_assoc_results(n: int = 100) -> List[Dict[str, Any]]:
    """Generate synthetic association results across several chromosomes."""
    rng = np.random.default_rng(42)
    results: List[Dict[str, Any]] = []
    for i in range(n):
        chrom = str((i % 22) + 1)
        pos = int(rng.integers(1_000_000, 250_000_000))
        pval = float(rng.uniform(1e-10, 1.0))
        beta = float(rng.normal(0, 0.5))
        se = float(abs(rng.normal(0.1, 0.05)))
        results.append(
            {
                "variant_id": f"rs{100000 + i}",
                "chrom": chrom,
                "pos": pos,
                "p_value": pval,
                "beta": beta,
                "se": se,
            }
        )
    return results


def _make_pca_data(n_samples: int = 30, n_pcs: int = 5) -> Dict[str, Any]:
    """Generate synthetic PCA data."""
    rng = np.random.default_rng(123)
    pcs = rng.standard_normal((n_samples, n_pcs))
    explained = np.sort(rng.uniform(0.02, 0.3, size=n_pcs))[::-1]
    explained = explained / explained.sum()
    sample_ids = [f"sample_{i}" for i in range(n_samples)]
    return {
        "pcs": pcs,
        "explained_variance_ratio": explained.tolist(),
        "sample_ids": sample_ids,
    }


def _make_kinship(n: int = 30) -> np.ndarray:
    """Generate a symmetric positive-semi-definite kinship matrix."""
    rng = np.random.default_rng(99)
    a = rng.standard_normal((n, n))
    km = a @ a.T / n
    np.fill_diagonal(km, 1.0)
    return km


def _make_metadata(sample_ids: List[str], populations: List[str]) -> Dict[str, Dict]:
    """Assign populations to sample ids round-robin."""
    meta: Dict[str, Dict] = {}
    for i, sid in enumerate(sample_ids):
        meta[sid] = {"population": populations[i % len(populations)]}
    return meta


# ---------------------------------------------------------------------------
# Tests: gwas_summary_panel
# ---------------------------------------------------------------------------


class TestGwasSummaryPanel:
    """Tests for gwas_summary_panel."""

    def test_full_panel(self, tmp_path: Path) -> None:
        """Full 2x2 panel with all four subplots populated."""
        results = _make_assoc_results(100)
        pca = _make_pca_data(30)
        kinship = _make_kinship(30)
        output = tmp_path / "summary_full.png"

        res = gwas_summary_panel(
            assoc_results=results,
            pca_data=pca,
            kinship_matrix=kinship,
            output_file=output,
        )

        assert res["status"] == "success"
        assert res["panels_generated"] == 4
        assert output.exists()
        assert output.stat().st_size > 0

    def test_results_only(self, tmp_path: Path) -> None:
        """Panel with only association results -- PCA and kinship show placeholder text."""
        results = _make_assoc_results(50)
        output = tmp_path / "summary_results_only.png"

        res = gwas_summary_panel(
            assoc_results=results,
            output_file=output,
        )

        assert res["status"] == "success"
        # Manhattan + QQ generated, PCA and Kinship skipped
        assert res["panels_generated"] == 2
        assert output.exists()
        assert output.stat().st_size > 0


# ---------------------------------------------------------------------------
# Tests: population_structure_panel
# ---------------------------------------------------------------------------


class TestPopulationStructurePanel:
    """Tests for population_structure_panel."""

    def test_full_panel(self, tmp_path: Path) -> None:
        """Full population structure panel with metadata and three populations."""
        n_samples = 30
        pca = _make_pca_data(n_samples, n_pcs=5)
        kinship = _make_kinship(n_samples)
        metadata = _make_metadata(
            pca["sample_ids"],
            populations=["PopA", "PopB", "PopC"],
        )
        output = tmp_path / "pop_structure.png"

        res = population_structure_panel(
            pca_data=pca,
            kinship_matrix=kinship,
            metadata=metadata,
            output_file=output,
        )

        assert res["status"] == "success"
        assert res["panels_generated"] == 4
        assert output.exists()
        assert output.stat().st_size > 0


# ---------------------------------------------------------------------------
# Tests: top_hit_detail_panel
# ---------------------------------------------------------------------------


class TestTopHitDetailPanel:
    """Tests for top_hit_detail_panel."""

    def test_with_genotypes(self, tmp_path: Path) -> None:
        """Detail panel with genotypes and phenotypes."""
        rng = np.random.default_rng(7)
        n_variants = 50
        n_samples = 80

        results = _make_assoc_results(n_variants)
        # Ensure several variants share the same chromosome / nearby positions
        # so the regional plot has content.
        base_chrom = results[0]["chrom"]
        base_pos = results[0]["pos"]
        for i in range(min(10, n_variants)):
            results[i]["chrom"] = base_chrom
            results[i]["pos"] = base_pos + int(rng.integers(-200_000, 200_000))

        genotypes = [rng.integers(0, 3, size=n_samples).tolist() for _ in range(n_variants)]
        phenotypes = rng.normal(10, 2, size=n_samples).tolist()

        output = tmp_path / "top_hit_detail.png"

        res = top_hit_detail_panel(
            assoc_results=results,
            genotypes=genotypes,
            phenotypes=phenotypes,
            variant_index=0,
            output_file=output,
        )

        assert res["status"] == "success"
        assert res["panels_generated"] == 3
        assert output.exists()
        assert output.stat().st_size > 0

    def test_without_genotypes(self, tmp_path: Path) -> None:
        """Detail panel without genotypes -- boxplot shows placeholder."""
        results = _make_assoc_results(20)
        output = tmp_path / "top_hit_no_geno.png"

        res = top_hit_detail_panel(
            assoc_results=results,
            output_file=output,
        )

        assert res["status"] == "success"
        # Regional plot + info panel generated, boxplot skipped
        assert res["panels_generated"] == 2
        assert output.exists()
        assert output.stat().st_size > 0
