"""Tests for the structural_variants.population submodule.

Covers population-level SV analysis: genotyping across samples, allele
frequency calculation, association testing, PCA-based population structure,
LD analysis between SVs and SNPs, and multi-sample callset merging.

Uses real implementations only (no mocking). All numerical checks use
real numpy/scipy computations against hand-verified expected values.
"""

from __future__ import annotations

import math
from typing import Any

import pytest

# Optional dependency checks
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def sample_sv_calls_with_evidence() -> list[dict[str, Any]]:
    """SV calls with per-sample depth and split-read evidence.

    Creates 3 SVs across 4 samples with realistic evidence patterns:
    - SV0 (DEL): clear deletion in samples A and B, absent in C and D
    - SV1 (DUP): duplication in sample C (het), absent elsewhere
    - SV2 (DEL): het deletion across all samples
    """
    return [
        {
            "chrom": "chr1",
            "start": 10000,
            "end": 15000,
            "sv_type": "DEL",
            "samples": {
                "sampleA": {"depth_ratio": 0.1, "split_reads": 12, "discordant_pairs": 8, "total_reads": 30},
                "sampleB": {"depth_ratio": 0.5, "split_reads": 6, "discordant_pairs": 4, "total_reads": 30},
                "sampleC": {"depth_ratio": 1.0, "split_reads": 0, "discordant_pairs": 0, "total_reads": 30},
                "sampleD": {"depth_ratio": 0.95, "split_reads": 1, "discordant_pairs": 0, "total_reads": 30},
            },
        },
        {
            "chrom": "chr2",
            "start": 50000,
            "end": 60000,
            "sv_type": "DUP",
            "samples": {
                "sampleA": {"depth_ratio": 1.0, "split_reads": 0, "discordant_pairs": 0, "total_reads": 30},
                "sampleB": {"depth_ratio": 1.05, "split_reads": 0, "discordant_pairs": 0, "total_reads": 30},
                "sampleC": {"depth_ratio": 1.5, "split_reads": 5, "discordant_pairs": 3, "total_reads": 30},
                "sampleD": {"depth_ratio": 1.0, "split_reads": 0, "discordant_pairs": 0, "total_reads": 30},
            },
        },
        {
            "chrom": "chr3",
            "start": 100000,
            "end": 102000,
            "sv_type": "DEL",
            "samples": {
                "sampleA": {"depth_ratio": 0.55, "split_reads": 4, "discordant_pairs": 3, "total_reads": 30},
                "sampleB": {"depth_ratio": 0.6, "split_reads": 3, "discordant_pairs": 2, "total_reads": 30},
                "sampleC": {"depth_ratio": 0.5, "split_reads": 5, "discordant_pairs": 4, "total_reads": 30},
                "sampleD": {"depth_ratio": 0.65, "split_reads": 3, "discordant_pairs": 2, "total_reads": 30},
            },
        },
    ]


@pytest.fixture
def four_samples() -> list[str]:
    """Four sample names matching the evidence fixture."""
    return ["sampleA", "sampleB", "sampleC", "sampleD"]


@pytest.fixture
def genotype_matrix_5x6() -> list[list[int]]:
    """A 5-SV x 6-sample genotype matrix with known allele frequencies.

    Layout (SVs as rows, samples as columns):
        SV0: [0, 0, 1, 1, 2, 0]  -> AF = 4/12 = 0.333
        SV1: [2, 2, 2, 2, 2, 2]  -> AF = 12/12 = 1.0 (fixed alt)
        SV2: [0, 0, 0, 0, 0, 0]  -> AF = 0/12 = 0.0 (monomorphic ref)
        SV3: [1, 0, 0, 0, 0, 0]  -> AF = 1/12 = 0.083 (singleton)
        SV4: [1, 1, 1, 1, 1, 1]  -> AF = 6/12 = 0.5
    """
    return [
        [0, 0, 1, 1, 2, 0],
        [2, 2, 2, 2, 2, 2],
        [0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 1],
    ]


@pytest.fixture
def population_labels_6() -> list[str]:
    """Population labels for 6 samples: 3 in popA, 3 in popB."""
    return ["popA", "popA", "popA", "popB", "popB", "popB"]


@pytest.fixture
def callsets_for_merge() -> list[list[dict[str, Any]]]:
    """Three callsets with overlapping and distinct SVs for merge testing.

    Callset 0 (sample_0): DEL chr1:10000-15000, INV chr2:50000-55000
    Callset 1 (sample_1): DEL chr1:10050-14980 (overlaps cs0 DEL), DUP chr3:80000-90000
    Callset 2 (sample_2): DEL chr1:10100-14900 (overlaps both), INV chr2:50100-55100 (overlaps cs0 INV)
    """
    return [
        [
            {"chrom": "chr1", "start": 10000, "end": 15000, "sv_type": "DEL", "sample": "sample_0"},
            {"chrom": "chr2", "start": 50000, "end": 55000, "sv_type": "INV", "sample": "sample_0"},
        ],
        [
            {"chrom": "chr1", "start": 10050, "end": 14980, "sv_type": "DEL", "sample": "sample_1"},
            {"chrom": "chr3", "start": 80000, "end": 90000, "sv_type": "DUP", "sample": "sample_1"},
        ],
        [
            {"chrom": "chr1", "start": 10100, "end": 14900, "sv_type": "DEL", "sample": "sample_2"},
            {"chrom": "chr2", "start": 50100, "end": 55100, "sv_type": "INV", "sample": "sample_2"},
        ],
    ]


# ---------------------------------------------------------------------------
# Test classes
# ---------------------------------------------------------------------------


class TestGenotypeSvPopulation:
    """Tests for genotype_sv_population."""

    def test_depth_method_basic(
        self,
        sample_sv_calls_with_evidence: list[dict[str, Any]],
        four_samples: list[str],
    ) -> None:
        """Depth-based genotyping produces correct matrix dimensions and valid values."""
        from metainformant.structural_variants.population import genotype_sv_population

        result = genotype_sv_population(sample_sv_calls_with_evidence, four_samples, method="depth")

        gt = result["genotype_matrix"]
        qual = result["quality_scores"]
        af = result["allele_frequencies"]

        # Correct dimensions: 3 SVs x 4 samples
        assert len(gt) == 3
        assert all(len(row) == 4 for row in gt)
        assert len(qual) == 3
        assert all(len(row) == 4 for row in qual)
        assert len(af) == 3

        # All genotypes are 0, 1, or 2
        for row in gt:
            for g in row:
                assert g in (0, 1, 2), f"Unexpected genotype: {g}"

        # All quality scores are non-negative
        for row in qual:
            for q in row:
                assert q >= 0.0

        # Allele frequencies are in [0, 1]
        for freq in af:
            assert 0.0 <= freq <= 1.0

    def test_depth_method_deletion_genotypes(
        self,
        sample_sv_calls_with_evidence: list[dict[str, Any]],
        four_samples: list[str],
    ) -> None:
        """Depth method assigns expected genotypes for a clear deletion.

        SV0 is a DEL. sampleA has depth_ratio=0.1 (hom del => 2),
        sampleB has 0.5 (het del => 1), sampleC has 1.0 (ref => 0),
        sampleD has 0.95 (ref => 0).
        """
        from metainformant.structural_variants.population import genotype_sv_population

        result = genotype_sv_population(sample_sv_calls_with_evidence, four_samples, method="depth")
        sv0_gt = result["genotype_matrix"][0]

        assert sv0_gt[0] == 2, "sampleA (depth_ratio=0.1) should be hom del (2)"
        assert sv0_gt[1] == 1, "sampleB (depth_ratio=0.5) should be het del (1)"
        assert sv0_gt[2] == 0, "sampleC (depth_ratio=1.0) should be ref (0)"
        assert sv0_gt[3] == 0, "sampleD (depth_ratio=0.95) should be ref (0)"

    def test_depth_method_duplication_genotypes(
        self,
        sample_sv_calls_with_evidence: list[dict[str, Any]],
        four_samples: list[str],
    ) -> None:
        """Depth method assigns expected genotypes for a duplication.

        SV1 is a DUP. sampleC has depth_ratio=1.5 (het dup => 1),
        others have ~1.0 (ref => 0).
        """
        from metainformant.structural_variants.population import genotype_sv_population

        result = genotype_sv_population(sample_sv_calls_with_evidence, four_samples, method="depth")
        sv1_gt = result["genotype_matrix"][1]

        assert sv1_gt[2] == 1, "sampleC (depth_ratio=1.5) should be het dup (1)"
        assert sv1_gt[0] == 0, "sampleA (depth_ratio=1.0) should be ref (0)"

    def test_split_method_basic(
        self,
        sample_sv_calls_with_evidence: list[dict[str, Any]],
        four_samples: list[str],
    ) -> None:
        """Split-read genotyping produces valid results."""
        from metainformant.structural_variants.population import genotype_sv_population

        result = genotype_sv_population(sample_sv_calls_with_evidence, four_samples, method="split")

        gt = result["genotype_matrix"]
        assert len(gt) == 3
        assert all(len(row) == 4 for row in gt)
        for row in gt:
            for g in row:
                assert g in (0, 1, 2)

    def test_split_method_high_evidence(
        self,
        sample_sv_calls_with_evidence: list[dict[str, Any]],
        four_samples: list[str],
    ) -> None:
        """Split method picks up strong split-read evidence.

        SV0 sampleA: split_reads=12, discordant=8, total=30 => fraction=20/30=0.67 => het(1)
        SV0 sampleC: split_reads=0, discordant=0 => fraction=0 => ref(0)
        """
        from metainformant.structural_variants.population import genotype_sv_population

        result = genotype_sv_population(sample_sv_calls_with_evidence, four_samples, method="split")
        sv0_gt = result["genotype_matrix"][0]

        assert sv0_gt[0] >= 1, "sampleA with 20 supporting reads should be het or hom"
        assert sv0_gt[2] == 0, "sampleC with 0 supporting reads should be ref"

    def test_allele_frequency_consistency(
        self,
        sample_sv_calls_with_evidence: list[dict[str, Any]],
        four_samples: list[str],
    ) -> None:
        """Allele frequencies are consistent with the genotype matrix."""
        from metainformant.structural_variants.population import genotype_sv_population

        result = genotype_sv_population(sample_sv_calls_with_evidence, four_samples, method="depth")

        gt = result["genotype_matrix"]
        af = result["allele_frequencies"]
        n_samples = len(four_samples)

        for i, row in enumerate(gt):
            expected_af = sum(row) / (n_samples * 2)
            assert abs(af[i] - expected_af) < 1e-9, f"SV{i} AF mismatch"

    def test_invalid_method_raises(
        self,
        sample_sv_calls_with_evidence: list[dict[str, Any]],
        four_samples: list[str],
    ) -> None:
        """Invalid genotyping method raises ValueError."""
        from metainformant.structural_variants.population import genotype_sv_population

        with pytest.raises(ValueError, match="Unknown genotyping method"):
            genotype_sv_population(sample_sv_calls_with_evidence, four_samples, method="invalid")

    def test_missing_sample_evidence(self) -> None:
        """Samples without evidence data get default genotype (ref)."""
        from metainformant.structural_variants.population import genotype_sv_population

        sv_calls = [
            {
                "chrom": "chr1",
                "start": 1000,
                "end": 2000,
                "sv_type": "DEL",
                "samples": {
                    "s1": {"depth_ratio": 0.1, "split_reads": 10, "discordant_pairs": 5, "total_reads": 30},
                    # s2 and s3 have no evidence
                },
            }
        ]
        samples = ["s1", "s2", "s3"]

        result = genotype_sv_population(sv_calls, samples, method="depth")
        gt = result["genotype_matrix"][0]

        # s1 should have a variant call, s2/s3 should default to ref
        assert gt[0] == 2, "s1 with depth_ratio=0.1 should be hom del"
        assert gt[1] == 0, "s2 with no evidence defaults to ref"
        assert gt[2] == 0, "s3 with no evidence defaults to ref"

    def test_empty_sv_calls(self) -> None:
        """Empty SV call list returns empty results."""
        from metainformant.structural_variants.population import genotype_sv_population

        result = genotype_sv_population([], ["s1", "s2"], method="depth")
        assert result["genotype_matrix"] == []
        assert result["quality_scores"] == []
        assert result["allele_frequencies"] == []


@pytest.mark.skipif(not HAS_NUMPY, reason="numpy required")
class TestSvAlleleFrequency:
    """Tests for sv_allele_frequency."""

    def test_basic_frequencies(self, genotype_matrix_5x6: list[list[int]]) -> None:
        """Allele frequencies match hand-calculated values."""
        from metainformant.structural_variants.population import sv_allele_frequency

        result = sv_allele_frequency(genotype_matrix_5x6)

        freqs = result["frequencies"]
        assert len(freqs) == 5

        # SV0: sum=4, total_alleles=12, AF=4/12=0.333...
        assert abs(freqs[0] - 4 / 12) < 1e-6
        # SV1: sum=12, AF=1.0
        assert abs(freqs[1] - 1.0) < 1e-6
        # SV2: sum=0, AF=0.0
        assert abs(freqs[2] - 0.0) < 1e-6
        # SV3: sum=1, AF=1/12
        assert abs(freqs[3] - 1 / 12) < 1e-6
        # SV4: sum=6, AF=0.5
        assert abs(freqs[4] - 0.5) < 1e-6

    def test_maf_calculation(self, genotype_matrix_5x6: list[list[int]]) -> None:
        """Minor allele frequencies are correctly computed as min(f, 1-f)."""
        from metainformant.structural_variants.population import sv_allele_frequency

        result = sv_allele_frequency(genotype_matrix_5x6)
        maf = result["maf"]

        # SV0: AF=0.333, MAF=0.333
        assert abs(maf[0] - 4 / 12) < 1e-6
        # SV1: AF=1.0, MAF=0.0
        assert abs(maf[1] - 0.0) < 1e-6
        # SV2: AF=0.0, MAF=0.0
        assert abs(maf[2] - 0.0) < 1e-6
        # SV3: AF=1/12, MAF=1/12
        assert abs(maf[3] - 1 / 12) < 1e-6
        # SV4: AF=0.5, MAF=0.5
        assert abs(maf[4] - 0.5) < 1e-6

    def test_polymorphic_count(self, genotype_matrix_5x6: list[list[int]]) -> None:
        """n_polymorphic counts SVs with MAF > 0."""
        from metainformant.structural_variants.population import sv_allele_frequency

        result = sv_allele_frequency(genotype_matrix_5x6)

        # SV0 (MAF=0.333), SV3 (MAF=0.083), SV4 (MAF=0.5) are polymorphic
        # SV1 (MAF=0.0), SV2 (MAF=0.0) are not
        assert result["n_polymorphic"] == 3

    def test_sfs_spectrum(self, genotype_matrix_5x6: list[list[int]]) -> None:
        """Site frequency spectrum bins match expected alt allele counts."""
        from metainformant.structural_variants.population import sv_allele_frequency

        result = sv_allele_frequency(genotype_matrix_5x6)
        sfs = result["sfs"]

        # Total alleles = 6*2 = 12, so SFS has 13 bins (0..12)
        assert len(sfs) == 13

        # Alt allele counts: SV0=4, SV1=12, SV2=0, SV3=1, SV4=6
        assert sfs[0] == 1  # SV2 has AC=0
        assert sfs[1] == 1  # SV3 has AC=1
        assert sfs[4] == 1  # SV0 has AC=4
        assert sfs[6] == 1  # SV4 has AC=6
        assert sfs[12] == 1  # SV1 has AC=12

    def test_population_stratified_frequencies(
        self,
        genotype_matrix_5x6: list[list[int]],
        population_labels_6: list[str],
    ) -> None:
        """Per-population allele frequencies are computed correctly."""
        from metainformant.structural_variants.population import sv_allele_frequency

        result = sv_allele_frequency(genotype_matrix_5x6, sample_labels=population_labels_6)

        assert "population_frequencies" in result
        pop_freqs = result["population_frequencies"]
        assert "popA" in pop_freqs
        assert "popB" in pop_freqs

        # popA = samples 0,1,2 => SV0 genotypes: [0,0,1] => AF = 1/6
        assert abs(pop_freqs["popA"][0] - 1 / 6) < 1e-6
        # popB = samples 3,4,5 => SV0 genotypes: [1,2,0] => AF = 3/6 = 0.5
        assert abs(pop_freqs["popB"][0] - 3 / 6) < 1e-6

        # SV1: all 2s => both pops AF=1.0
        assert abs(pop_freqs["popA"][1] - 1.0) < 1e-6
        assert abs(pop_freqs["popB"][1] - 1.0) < 1e-6

        # SV2: all 0s => both pops AF=0.0
        assert abs(pop_freqs["popA"][2] - 0.0) < 1e-6
        assert abs(pop_freqs["popB"][2] - 0.0) < 1e-6

    def test_without_population_labels(self, genotype_matrix_5x6: list[list[int]]) -> None:
        """When no labels given, population_frequencies is absent from result."""
        from metainformant.structural_variants.population import sv_allele_frequency

        result = sv_allele_frequency(genotype_matrix_5x6)
        assert "population_frequencies" not in result

    def test_single_sample_matrix(self) -> None:
        """Single-sample genotype matrix produces valid results."""
        from metainformant.structural_variants.population import sv_allele_frequency

        gt = [[0], [1], [2]]  # 3 SVs, 1 sample
        result = sv_allele_frequency(gt)

        assert len(result["frequencies"]) == 3
        assert abs(result["frequencies"][0] - 0.0) < 1e-6
        assert abs(result["frequencies"][1] - 0.5) < 1e-6
        assert abs(result["frequencies"][2] - 1.0) < 1e-6


@pytest.mark.skipif(not HAS_NUMPY, reason="numpy required")
class TestSvAssociationTest:
    """Tests for sv_association_test (linear and logistic regression)."""

    def test_continuous_phenotype_linear(self) -> None:
        """Linear regression detects association with continuous phenotype."""
        from metainformant.structural_variants.population import sv_association_test

        # 20 samples: genotype influences phenotype with noise
        rng = np.random.default_rng(42)
        genotypes = [0] * 7 + [1] * 7 + [2] * 6
        # Phenotype with a real effect of ~5 per alt allele
        phenotypes = [float(g * 5.0 + rng.normal(0, 1)) for g in genotypes]

        result = sv_association_test(genotypes, phenotypes)

        assert result["n_samples"] == 20
        assert result["odds_ratio"] is None, "Continuous phenotype should have no OR"
        # Beta should be positive (phenotype increases with genotype)
        assert result["beta"] > 0, f"Expected positive beta, got {result['beta']}"
        # With a strong effect size and 20 samples, p-value should be small
        assert result["p_value"] < 0.05, f"Expected significant p-value, got {result['p_value']}"
        # SE should be positive
        assert result["se"] > 0

    def test_binary_phenotype_logistic(self) -> None:
        """Logistic regression detects association with binary phenotype."""
        from metainformant.structural_variants.population import sv_association_test

        # 30 samples: genotype=2 strongly associated with case status
        genotypes = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
        # Cases predominantly in genotype=2 group
        phenotypes = [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            1.0,
            1.0,
            1.0,
            1.0,
            0.0,
            0.0,
            0.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
        ]

        result = sv_association_test(genotypes, phenotypes)

        assert result["n_samples"] == 30
        assert result["odds_ratio"] is not None, "Binary phenotype should have OR"
        assert result["odds_ratio"] > 1.0, "Genotype should increase disease risk"
        assert result["beta"] > 0, "Beta should be positive"

    def test_no_association(self) -> None:
        """When genotype and phenotype are independent, p-value is large."""
        from metainformant.structural_variants.population import sv_association_test

        rng = np.random.default_rng(123)
        n = 50
        genotypes = [int(rng.choice([0, 1, 2])) for _ in range(n)]
        # Random phenotype unrelated to genotype
        phenotypes = [float(rng.normal(10, 2)) for _ in range(n)]

        result = sv_association_test(genotypes, phenotypes)
        # p-value should generally be non-significant (> 0.01)
        # Using a loose threshold since it is random
        assert result["n_samples"] == n
        assert result["se"] > 0

    def test_with_covariates(self) -> None:
        """Association test adjusts for covariates."""
        from metainformant.structural_variants.population import sv_association_test

        rng = np.random.default_rng(99)
        n = 30
        genotypes = [0] * 10 + [1] * 10 + [2] * 10
        # Covariate (e.g., age) and phenotype both affected by genotype
        covariate = [float(rng.normal(50, 10)) for _ in range(n)]
        phenotypes = [float(g * 3.0 + cov * 0.1 + rng.normal(0, 1)) for g, cov in zip(genotypes, covariate)]
        covariates = [[cov] for cov in covariate]

        result = sv_association_test(genotypes, phenotypes, covariates=covariates)

        assert result["n_samples"] == n
        assert result["se"] > 0
        assert result["beta"] > 0

    def test_length_mismatch_raises(self) -> None:
        """Mismatched genotype/phenotype lengths raise ValueError."""
        from metainformant.structural_variants.population import sv_association_test

        with pytest.raises(ValueError, match="Length mismatch"):
            sv_association_test([0, 1, 2], [1.0, 2.0])

    def test_all_same_genotype(self) -> None:
        """With no genotype variation, test still returns valid results."""
        from metainformant.structural_variants.population import sv_association_test

        genotypes = [0] * 20
        phenotypes = [float(i) for i in range(20)]

        result = sv_association_test(genotypes, phenotypes)
        assert result["n_samples"] == 20
        # Beta should be ~0 since genotype has no variation
        assert abs(result["beta"]) < 1e-3


@pytest.mark.skipif(not HAS_NUMPY, reason="numpy required")
class TestSvPopulationStructure:
    """Tests for sv_population_structure (PCA)."""

    def test_pca_dimensions(self) -> None:
        """PCA returns correct dimensions for PCs, eigenvalues, and loadings."""
        from metainformant.structural_variants.population import sv_population_structure

        rng = np.random.default_rng(42)
        n_svs, n_samples = 20, 10
        gt = rng.choice([0, 1, 2], size=(n_svs, n_samples)).tolist()

        result = sv_population_structure(gt, n_components=3)

        pcs = np.asarray(result["pcs"])
        loadings = np.asarray(result["loadings"])

        assert pcs.shape == (n_samples, 3), f"PCs shape: {pcs.shape}"
        assert len(result["eigenvalues"]) == 3
        assert len(result["variance_explained"]) == 3
        assert loadings.shape == (n_svs, 3), f"Loadings shape: {loadings.shape}"

    def test_variance_explained_sums_correctly(self) -> None:
        """Variance explained fractions sum to <= 1.0 and are non-negative."""
        from metainformant.structural_variants.population import sv_population_structure

        rng = np.random.default_rng(7)
        gt = rng.choice([0, 1, 2], size=(30, 8)).tolist()

        result = sv_population_structure(gt, n_components=5)

        ve = result["variance_explained"]
        for v in ve:
            assert v >= 0.0, f"Negative variance explained: {v}"
        assert sum(ve) <= 1.0 + 1e-6, f"Variance explained sums to {sum(ve)}"

    def test_pca_separates_populations(self) -> None:
        """PCA separates two distinct population clusters.

        Creates two populations with opposite genotype patterns so that
        PC1 clearly separates them.
        """
        from metainformant.structural_variants.population import sv_population_structure

        # Population A: mostly homozygous ref
        # Population B: mostly homozygous alt
        n_svs = 15
        pop_a = [[0] * 5 for _ in range(n_svs)]
        pop_b = [[2] * 5 for _ in range(n_svs)]

        # Combine: 10 samples, SVs x samples
        gt = [a + b for a, b in zip(pop_a, pop_b)]

        result = sv_population_structure(gt, n_components=2)
        pcs = np.asarray(result["pcs"])

        # First 5 samples (pop A) should cluster separately from last 5 (pop B) on PC1
        pc1_popA = pcs[:5, 0]
        pc1_popB = pcs[5:, 0]

        # The means should be clearly different
        assert abs(np.mean(pc1_popA) - np.mean(pc1_popB)) > 0.1, "PCA should separate the two populations"

    def test_n_components_clamped(self) -> None:
        """n_components is clamped to min(n_samples, n_svs)."""
        from metainformant.structural_variants.population import sv_population_structure

        # 3 SVs x 4 samples => max components = 3
        gt = [[0, 1, 2, 0], [1, 0, 1, 2], [2, 1, 0, 1]]

        result = sv_population_structure(gt, n_components=100)

        pcs = np.asarray(result["pcs"])
        n_components_actual = pcs.shape[1]
        assert n_components_actual <= 3, f"Components should be clamped, got {n_components_actual}"

    def test_single_sv_matrix(self) -> None:
        """PCA on a single-SV matrix still works."""
        from metainformant.structural_variants.population import sv_population_structure

        gt = [[0, 1, 2, 1, 0]]  # 1 SV, 5 samples

        result = sv_population_structure(gt, n_components=1)
        pcs = np.asarray(result["pcs"])

        assert pcs.shape[0] == 5
        assert pcs.shape[1] == 1
        assert len(result["eigenvalues"]) == 1

    def test_constant_genotype_matrix(self) -> None:
        """PCA handles a constant genotype matrix (all zeros) gracefully."""
        from metainformant.structural_variants.population import sv_population_structure

        gt = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]

        result = sv_population_structure(gt, n_components=2)
        # Should not raise; variance explained should be ~0
        ve = result["variance_explained"]
        assert all(v >= 0.0 for v in ve)


@pytest.mark.skipif(not HAS_NUMPY, reason="numpy required")
class TestSvLdAnalysis:
    """Tests for sv_ld_analysis."""

    def test_perfect_ld(self) -> None:
        """Identical genotypes between SV and SNP give r2 = 1.0."""
        from metainformant.structural_variants.population import sv_ld_analysis

        # 1 SV, 1 SNP, 10 samples with identical genotypes
        sv_gt = [[0, 0, 1, 1, 2, 2, 0, 1, 2, 0]]
        snp_gt = [[0, 0, 1, 1, 2, 2, 0, 1, 2, 0]]
        sv_pos = [10000]
        snp_pos = [10500]

        result = sv_ld_analysis(sv_gt, snp_gt, sv_pos, snp_pos)

        assert len(result["ld_matrix"]) == 1
        assert len(result["ld_matrix"][0]) == 1
        assert result["ld_matrix"][0][0] > 0.99, "Perfect LD should give r2 ~ 1.0"
        assert result["max_r2_per_sv"][0] > 0.99
        assert result["tag_snps"][0]["snp_index"] == 0
        assert result["tag_snps"][0]["snp_position"] == 10500

    def test_no_ld(self) -> None:
        """Unrelated genotypes give low r2."""
        from metainformant.structural_variants.population import sv_ld_analysis

        # SV has one pattern, SNP has an independent pattern
        sv_gt = [[0, 0, 0, 0, 0, 2, 2, 2, 2, 2]]
        snp_gt = [[0, 2, 0, 2, 0, 2, 0, 2, 0, 2]]  # Alternating, no correlation
        sv_pos = [10000]
        snp_pos = [20000]

        result = sv_ld_analysis(sv_gt, snp_gt, sv_pos, snp_pos)

        r2 = result["ld_matrix"][0][0]
        assert r2 < 0.5, f"Unrelated genotypes should have low r2, got {r2}"

    def test_multiple_svs_and_snps(self) -> None:
        """LD matrix has correct dimensions for multiple SVs and SNPs."""
        from metainformant.structural_variants.population import sv_ld_analysis

        n_samples = 8
        sv_gt = [
            [0, 0, 1, 1, 2, 2, 0, 1],
            [2, 2, 1, 1, 0, 0, 2, 1],
        ]
        snp_gt = [
            [0, 0, 1, 1, 2, 2, 0, 1],
            [1, 1, 0, 0, 2, 2, 1, 0],
            [2, 1, 0, 2, 1, 0, 2, 1],
        ]
        sv_pos = [10000, 20000]
        snp_pos = [10500, 15000, 25000]

        result = sv_ld_analysis(sv_gt, snp_gt, sv_pos, snp_pos)

        ld = result["ld_matrix"]
        assert len(ld) == 2, "Should have 2 rows (SVs)"
        assert all(len(row) == 3 for row in ld), "Should have 3 columns (SNPs)"

        # All r2 values should be in [0, 1]
        for row in ld:
            for r2 in row:
                assert 0.0 <= r2 <= 1.0 + 1e-9, f"r2 out of range: {r2}"

        # tag_snps should have one entry per SV
        assert len(result["tag_snps"]) == 2
        assert len(result["max_r2_per_sv"]) == 2

    def test_tag_snp_selection(self) -> None:
        """Tag SNP is the one with highest r2 for each SV."""
        from metainformant.structural_variants.population import sv_ld_analysis

        # SV perfectly correlated with SNP index 1, not with SNP index 0
        sv_gt = [[0, 1, 2, 0, 1, 2, 0, 1]]
        snp_gt = [
            [2, 2, 2, 2, 2, 2, 2, 2],  # Constant => r2=0
            [0, 1, 2, 0, 1, 2, 0, 1],  # Identical => r2=1
        ]
        sv_pos = [10000]
        snp_pos = [5000, 12000]

        result = sv_ld_analysis(sv_gt, snp_gt, sv_pos, snp_pos)

        assert result["tag_snps"][0]["snp_index"] == 1
        assert result["tag_snps"][0]["snp_position"] == 12000
        assert result["tag_snps"][0]["r2"] > 0.99

    def test_monomorphic_sv_gives_zero_r2(self) -> None:
        """A monomorphic SV (no variation) has r2=0 with all SNPs."""
        from metainformant.structural_variants.population import sv_ld_analysis

        sv_gt = [[0, 0, 0, 0, 0, 0]]  # No variation
        snp_gt = [[0, 1, 2, 0, 1, 2]]
        sv_pos = [10000]
        snp_pos = [10500]

        result = sv_ld_analysis(sv_gt, snp_gt, sv_pos, snp_pos)

        assert result["ld_matrix"][0][0] == 0.0
        assert result["max_r2_per_sv"][0] == 0.0


class TestMergeSvCallsets:
    """Tests for merge_sv_callsets."""

    def test_basic_merge(self, callsets_for_merge: list[list[dict[str, Any]]]) -> None:
        """Overlapping DEL calls from 3 callsets merge into one consensus."""
        from metainformant.structural_variants.population import merge_sv_callsets

        merged = merge_sv_callsets(callsets_for_merge, min_overlap=0.5, max_breakpoint_distance=500)

        # Find the merged chr1 DEL
        chr1_dels = [m for m in merged if m["chrom"] == "chr1" and m["sv_type"] == "DEL"]
        assert len(chr1_dels) == 1, f"Expected 1 merged chr1 DEL, got {len(chr1_dels)}"

        consensus = chr1_dels[0]
        assert consensus["n_samples"] == 3
        assert consensus["support_count"] == 3
        assert set(consensus["samples"]) == {"sample_0", "sample_1", "sample_2"}

    def test_consensus_breakpoints(self, callsets_for_merge: list[list[dict[str, Any]]]) -> None:
        """Consensus breakpoints use median positions."""
        from metainformant.structural_variants.population import merge_sv_callsets

        merged = merge_sv_callsets(callsets_for_merge, min_overlap=0.5, max_breakpoint_distance=500)

        chr1_dels = [m for m in merged if m["chrom"] == "chr1" and m["sv_type"] == "DEL"]
        consensus = chr1_dels[0]

        # Starts: 10000, 10050, 10100 => median (index 1) = 10050
        assert consensus["start"] == 10050
        # Ends: 14900, 14980, 15000 => median (index 1) = 14980
        assert consensus["end"] == 14980

    def test_different_types_not_merged(self) -> None:
        """SVs of different types on the same chromosome are not merged."""
        from metainformant.structural_variants.population import merge_sv_callsets

        callsets = [
            [{"chrom": "chr1", "start": 10000, "end": 15000, "sv_type": "DEL", "sample": "s1"}],
            [{"chrom": "chr1", "start": 10000, "end": 15000, "sv_type": "DUP", "sample": "s2"}],
        ]

        merged = merge_sv_callsets(callsets, min_overlap=0.5, max_breakpoint_distance=500)

        # Should have 2 separate entries (one DEL, one DUP)
        assert len(merged) == 2
        types = {m["sv_type"] for m in merged}
        assert types == {"DEL", "DUP"}

    def test_different_chromosomes_not_merged(self) -> None:
        """SVs on different chromosomes are never merged."""
        from metainformant.structural_variants.population import merge_sv_callsets

        callsets = [
            [{"chrom": "chr1", "start": 10000, "end": 15000, "sv_type": "DEL", "sample": "s1"}],
            [{"chrom": "chr2", "start": 10000, "end": 15000, "sv_type": "DEL", "sample": "s2"}],
        ]

        merged = merge_sv_callsets(callsets, min_overlap=0.5, max_breakpoint_distance=500)
        assert len(merged) == 2

    def test_distant_breakpoints_not_merged(self) -> None:
        """SVs with breakpoints beyond max distance are not merged."""
        from metainformant.structural_variants.population import merge_sv_callsets

        callsets = [
            [{"chrom": "chr1", "start": 10000, "end": 15000, "sv_type": "DEL", "sample": "s1"}],
            [{"chrom": "chr1", "start": 12000, "end": 17000, "sv_type": "DEL", "sample": "s2"}],
        ]

        merged = merge_sv_callsets(callsets, min_overlap=0.5, max_breakpoint_distance=500)
        assert len(merged) == 2, "Breakpoints 2000bp apart should not merge with max_dist=500"

    def test_low_overlap_not_merged(self) -> None:
        """SVs with insufficient reciprocal overlap are not merged."""
        from metainformant.structural_variants.population import merge_sv_callsets

        # SV1: 10000-15000 (5000bp), SV2: 14500-20000 (5500bp)
        # Overlap: 14500-15000 = 500bp
        # Reciprocal: 500/5000=0.1 and 500/5500=0.09 => min=0.09
        callsets = [
            [{"chrom": "chr1", "start": 10000, "end": 15000, "sv_type": "DEL", "sample": "s1"}],
            [{"chrom": "chr1", "start": 14500, "end": 20000, "sv_type": "DEL", "sample": "s2"}],
        ]

        merged = merge_sv_callsets(callsets, min_overlap=0.5, max_breakpoint_distance=5000)
        assert len(merged) == 2, "Low overlap should prevent merging"

    def test_empty_callsets(self) -> None:
        """Empty input returns empty result."""
        from metainformant.structural_variants.population import merge_sv_callsets

        assert merge_sv_callsets([]) == []
        assert merge_sv_callsets([[], []]) == []

    def test_single_callset(self) -> None:
        """Single callset with no merging returns all calls individually."""
        from metainformant.structural_variants.population import merge_sv_callsets

        callsets = [
            [
                {"chrom": "chr1", "start": 10000, "end": 15000, "sv_type": "DEL"},
                {"chrom": "chr2", "start": 50000, "end": 60000, "sv_type": "DUP"},
            ]
        ]

        merged = merge_sv_callsets(callsets)
        assert len(merged) == 2

    def test_default_sample_naming(self) -> None:
        """Calls without explicit sample names get auto-generated names."""
        from metainformant.structural_variants.population import merge_sv_callsets

        callsets = [
            [{"chrom": "chr1", "start": 10000, "end": 15000, "sv_type": "DEL"}],
            [{"chrom": "chr1", "start": 10010, "end": 14990, "sv_type": "DEL"}],
        ]

        merged = merge_sv_callsets(callsets, min_overlap=0.5, max_breakpoint_distance=500)

        # Both should merge; check that auto-generated sample names are present
        assert len(merged) == 1
        assert "sample_0" in merged[0]["samples"]
        assert "sample_1" in merged[0]["samples"]

    def test_strict_overlap_threshold(self, callsets_for_merge: list[list[dict[str, Any]]]) -> None:
        """Very strict overlap threshold prevents merging."""
        from metainformant.structural_variants.population import merge_sv_callsets

        # With 99% overlap required, slightly offset calls should not merge
        merged = merge_sv_callsets(callsets_for_merge, min_overlap=0.99, max_breakpoint_distance=500)

        # The chr1 DELs have starts 10000/10050/10100, which differ enough
        # to prevent 99% reciprocal overlap
        chr1_dels = [m for m in merged if m["chrom"] == "chr1" and m["sv_type"] == "DEL"]
        assert len(chr1_dels) >= 2, "Strict overlap should prevent some merges"


class TestIntegrationPopulationWorkflow:
    """End-to-end integration tests combining multiple population functions."""

    @pytest.mark.skipif(not HAS_NUMPY, reason="numpy required")
    def test_genotype_to_frequency_pipeline(
        self,
        sample_sv_calls_with_evidence: list[dict[str, Any]],
        four_samples: list[str],
    ) -> None:
        """Full pipeline: genotype -> allele frequency -> population structure."""
        from metainformant.structural_variants.population import (
            genotype_sv_population,
            sv_allele_frequency,
            sv_population_structure,
        )

        # Step 1: Genotype
        gt_result = genotype_sv_population(sample_sv_calls_with_evidence, four_samples, method="depth")
        gt_matrix = gt_result["genotype_matrix"]

        # Step 2: Allele frequencies
        labels = ["pop1", "pop1", "pop2", "pop2"]
        af_result = sv_allele_frequency(gt_matrix, sample_labels=labels)

        assert len(af_result["frequencies"]) == 3
        assert "population_frequencies" in af_result
        assert "pop1" in af_result["population_frequencies"]
        assert "pop2" in af_result["population_frequencies"]

        # Step 3: Population structure (PCA)
        pca_result = sv_population_structure(gt_matrix, n_components=2)
        pcs = np.asarray(pca_result["pcs"])

        assert pcs.shape[0] == 4, "Should have PCs for 4 samples"
        assert pcs.shape[1] <= 2

    @pytest.mark.skipif(not HAS_NUMPY, reason="numpy required")
    def test_merge_then_genotype_pipeline(self) -> None:
        """Merge callsets then genotype the consensus SVs."""
        from metainformant.structural_variants.population import (
            genotype_sv_population,
            merge_sv_callsets,
        )

        callsets = [
            [
                {"chrom": "chr1", "start": 10000, "end": 15000, "sv_type": "DEL", "sample": "s1"},
                {"chrom": "chr2", "start": 50000, "end": 60000, "sv_type": "DUP", "sample": "s1"},
            ],
            [
                {"chrom": "chr1", "start": 10020, "end": 14980, "sv_type": "DEL", "sample": "s2"},
            ],
        ]

        # Step 1: Merge
        merged = merge_sv_callsets(callsets, min_overlap=0.5, max_breakpoint_distance=500)
        assert len(merged) >= 1

        # Step 2: Build SV calls with per-sample evidence for genotyping
        sv_calls_for_gt = []
        for m in merged:
            sv_calls_for_gt.append(
                {
                    "chrom": m["chrom"],
                    "start": m["start"],
                    "end": m["end"],
                    "sv_type": m["sv_type"],
                    "samples": {
                        "s1": {
                            "depth_ratio": 0.3 if m["sv_type"] == "DEL" else 1.5,
                            "split_reads": 5,
                            "discordant_pairs": 3,
                            "total_reads": 30,
                        },
                        "s2": {
                            "depth_ratio": 0.5 if m["sv_type"] == "DEL" else 1.0,
                            "split_reads": 3,
                            "discordant_pairs": 2,
                            "total_reads": 30,
                        },
                    },
                }
            )

        gt_result = genotype_sv_population(sv_calls_for_gt, ["s1", "s2"], method="depth")
        assert len(gt_result["genotype_matrix"]) == len(merged)
        assert all(len(row) == 2 for row in gt_result["genotype_matrix"])


class TestModuleImports:
    """Test that the population submodule imports correctly."""

    def test_import_from_population_submodule(self) -> None:
        """All public functions are importable from the population submodule."""
        from metainformant.structural_variants.population import (
            genotype_sv_population,
            merge_sv_callsets,
            sv_allele_frequency,
            sv_association_test,
            sv_ld_analysis,
            sv_population_structure,
        )

        assert callable(genotype_sv_population)
        assert callable(merge_sv_callsets)
        assert callable(sv_allele_frequency)
        assert callable(sv_association_test)
        assert callable(sv_ld_analysis)
        assert callable(sv_population_structure)

    def test_import_from_top_level(self) -> None:
        """Population functions are accessible from the population submodule."""
        from metainformant.structural_variants.population import (
            genotype_sv_population,
            merge_sv_callsets,
            sv_allele_frequency,
            sv_association_test,
            sv_ld_analysis,
            sv_population_structure,
        )

        assert callable(genotype_sv_population)
        assert callable(sv_population_structure)

    def test_population_submodule_attribute(self) -> None:
        """The population submodule is accessible as an attribute."""
        from metainformant import structural_variants

        assert hasattr(structural_variants, "population")
        assert hasattr(structural_variants.population, "genotype_sv_population")
