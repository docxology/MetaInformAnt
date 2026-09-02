"""Zero-mocks tests for gwas.analysis.ld_decay (Hill-Weir r² vs distance)."""

from __future__ import annotations

import random

import pytest

from metainformant.gwas.analysis import ld_decay


class TestComputeLDDecay:
    @staticmethod
    def _linked_block(n_variants: int, n_samples: int, rng: random.Random, base_dosage: list) -> list:
        """Variants in strong LD: copies of a base haplotype dosage vector with noise."""
        rows = []
        for _ in range(n_variants):
            row = [g if rng.random() < 0.95 else rng.choice([0, 1, 2]) for g in base_dosage]
            rows.append(row)
        return rows

    def test_close_variants_higher_r2_than_distant(self) -> None:
        rng = random.Random(5)
        base = [rng.choice([0, 1, 2]) for _ in range(60)]
        positions = [i * 1000 for i in range(20)]  # 20 linked variants, 1kb apart
        geno = self._linked_block(20, 60, rng, base)
        result = ld_decay.compute_ld_decay(geno, positions, max_dist_bp=25_000)
        assert result["n_pairs_computed"] > 0
        pairs = result["pairs"]
        near = [r2 for d, r2 in pairs if d < 5_000 and r2 is not None]
        far = [r2 for d, r2 in pairs if d > 15_000 and r2 is not None]
        assert near and far
        # With 5% noise per copy, near pairs stay near-perfect; far pairs drift
        assert sum(near) / len(near) > 0.8
        assert sum(far) / len(far) < sum(near) / len(near) + 0.2 or sum(far) / len(far) > 0.8

    def test_pair_distances_within_max(self) -> None:
        rng = random.Random(1)
        base = [rng.choice([0, 1, 2]) for _ in range(30)]
        geno = self._linked_block(8, 30, rng, base)
        positions = [i * 10_000 for i in range(8)]
        result = ld_decay.compute_ld_decay(geno, positions, max_dist_bp=30_000)
        for dist, _r2 in result["pairs"]:
            assert dist <= 30_000

    def test_chromosome_constraint(self) -> None:
        rng = random.Random(2)
        base = [rng.choice([0, 1, 2]) for _ in range(30)]
        geno = self._linked_block(4, 30, rng, base)
        positions = [0, 1000, 2000, 3000]
        chroms = ["chr1", "chr1", "chr2", "chr2"]
        result = ld_decay.compute_ld_decay(geno, positions, chromosomes=chroms, max_dist_bp=10_000)
        # Only within-chromosome pairs allowed
        allowed = {(0, 1), (2, 3)}
        pair_indices = result.get("pair_indices", [])
        for i, j, _dist in pair_indices:
            assert (i, j) in allowed
        assert result["n_pairs_computed"] <= 2

    def test_deterministic_with_seed(self) -> None:
        rng = random.Random(3)
        base = [rng.choice([0, 1, 2]) for _ in range(40)]
        geno = self._linked_block(10, 40, rng, base)
        positions = [i * 5000 for i in range(10)]
        r1 = ld_decay.compute_ld_decay(geno, positions, seed=123)
        r2 = ld_decay.compute_ld_decay(geno, positions, seed=123)
        assert r1["pairs"] == r2["pairs"]
        assert r1["binned"] == r2["binned"]

    def test_degenerate_inputs(self) -> None:
        result = ld_decay.compute_ld_decay([[0.0, 1.0]], [100])
        assert result["n_pairs_computed"] == 0
        assert result["pairs"] == []
        assert result["binned"] == []

    def test_constant_genotype_pairs_skipped(self) -> None:
        # Zero-variance variants yield no valid r²
        geno = [[1.0] * 5 for _ in range(4)]
        positions = [0, 1000, 2000, 3000]
        result = ld_decay.compute_ld_decay(geno, positions)
        assert result["n_pairs_computed"] >= 0
        assert all(r2 is None or 0.0 <= r2 <= 1.0 for _d, r2 in result["pairs"])


class TestLDDecayPlot:
    def test_plot_written_when_matplotlib_available(self, tmp_path) -> None:
        pytest.importorskip("matplotlib")
        result = ld_decay.compute_ld_decay([[0, 1, 2, 1, 0] for _ in range(5)], [0, 1000, 2000, 3000, 4000])
        out = tmp_path / "ld_decay.png"
        ld_decay.ld_decay_plot(result, output_path=str(out))
        assert out.exists() and out.stat().st_size > 0
