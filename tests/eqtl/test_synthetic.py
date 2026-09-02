"""Tests for metainformant.eqtl.synthetic (zero-mocks, real frames)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from metainformant.eqtl.synthetic import (
    create_synthetic_data,
    create_synthetic_genotypes,
    parse_gene_positions,
)


class TestCreateSyntheticData:
    def test_shapes_and_determinism(self):
        expr, geno, gene_pos, var_pos = create_synthetic_data(
            n_genes=10, n_variants=50, n_samples=20, seed=42
        )
        assert expr.shape == (10, 20)
        assert geno.shape == (50, 20)
        assert len(gene_pos) == 10
        assert len(var_pos) == 50
        # Determinism: same seed -> identical matrices
        expr2, geno2, _, _ = create_synthetic_data(
            n_genes=10, n_variants=50, n_samples=20, seed=42
        )
        pd.testing.assert_frame_equal(expr, expr2)
        pd.testing.assert_frame_equal(geno, geno2)

    def test_columns_align_across_matrices(self):
        expr, geno, _, _ = create_synthetic_data(
            n_genes=6, n_variants=30, n_samples=15, seed=7
        )
        assert list(expr.columns) == list(geno.columns)
        assert all(col.startswith("SRR") for col in expr.columns)

    def test_genotype_dosages_valid(self):
        _, geno, _, _ = create_synthetic_data(
            n_genes=6, n_variants=30, n_samples=40, seed=3
        )
        assert set(np.unique(geno.values)).issubset({0, 1, 2})

    def test_true_eqtl_effects_present(self):
        """With 30 genes and seed 42, some genes must carry injected effects."""
        expr, geno, _, _ = create_synthetic_data(
            n_genes=30, n_variants=150, n_samples=60, seed=42
        )
        # Correlate each gene against its first cis variant
        n_correlated = 0
        for g in range(30):
            r = np.corrcoef(expr.iloc[g].values, geno.iloc[g * 5].values)[0, 1]
            if abs(r) > 0.3:
                n_correlated += 1
        assert n_correlated >= 3

    def test_invalid_variant_count_raises(self):
        with pytest.raises(ValueError, match="n_variants"):
            create_synthetic_data(n_genes=10, n_variants=10, n_samples=20)

    def test_invalid_sample_count_raises(self):
        with pytest.raises(ValueError, match="n_samples"):
            create_synthetic_data(n_genes=2, n_variants=10, n_samples=1)


class TestParseGenePositions:
    def test_nc_accession_parsed(self):
        df = parse_gene_positions(["lcl|NC_037638.1_mrna_XM_623972.6_1"])
        assert df.iloc[0]["chrom"].startswith("NC_")

    def test_unknown_target_gets_unknown_chrom(self):
        df = parse_gene_positions(["weird_id"])
        assert df.iloc[0]["chrom"] in {"unknown", "weird"}

    def test_positions_monotonic(self):
        df = parse_gene_positions([f"g{i}_extra_parts" for i in range(5)])
        pos = df["tss_position"].tolist()
        assert pos == sorted(pos)


class TestCreateSyntheticGenotypes:
    def test_matches_gene_positions(self):
        gene_pos = pd.DataFrame(
            {
                "gene_id": ["LOC1", "LOC2"],
                "chrom": ["1", "2"],
                "tss_position": [1_000_000, 2_000_000],
            }
        )
        geno, var_pos = create_synthetic_genotypes(
            [f"S{i}" for i in range(12)], gene_pos, variants_per_gene=3, seed=1
        )
        assert geno.shape == (6, 12)
        assert len(var_pos) == 6
        assert set(np.unique(geno.values)).issubset({0, 1, 2})
        # Positions anchored to TSS +/- 5000
        assert var_pos.iloc[0]["position"] == 1_000_000 - 5000

    def test_determinism(self):
        gene_pos = pd.DataFrame(
            {"gene_id": ["LOC1"], "chrom": ["1"], "tss_position": [1_000_000]}
        )
        g1, _ = create_synthetic_genotypes(["A", "B"], gene_pos, seed=5)
        g2, _ = create_synthetic_genotypes(["A", "B"], gene_pos, seed=5)
        pd.testing.assert_frame_equal(g1, g2)
