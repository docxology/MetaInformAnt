"""Zero-mocks tests for gwas.analysis.hwe (exactness-checked chi-square HWE QC)."""

from __future__ import annotations

import math

import pytest

from metainformant.gwas.analysis import hwe


class TestHweChi2Test:
    def test_perfect_hwe_genotype_counts(self) -> None:
        # p=0.5, n=100: expected 25/50/25; observed identical → chi2=0, p=1
        chi2, p = hwe.hwe_chi2_test(25, 50, 25)
        assert chi2 == pytest.approx(0.0, abs=1e-9)
        assert p == pytest.approx(1.0)

    def test_extreme_heterozygote_excess(self) -> None:
        # All heterozygotes: massive deviation from HWE
        chi2, p = hwe.hwe_chi2_test(0, 100, 0)
        assert chi2 > 50.0
        assert p < 1e-10

    def test_heterozygote_deficit(self) -> None:
        # 45/10/45 vs HWE expected ~27/45/27 for p≈0.5
        chi2, p = hwe.hwe_chi2_test(45, 10, 45)
        assert chi2 > 30.0
        assert p < 1e-6

    def test_monomorphic_returns_neutral(self) -> None:
        assert hwe.hwe_chi2_test(100, 0, 0) == (0.0, 1.0)
        assert hwe.hwe_chi2_test(0, 0, 0) == (0.0, 1.0)

    def test_chi2_statistic_closed_form(self) -> None:
        # Hand-computed: n=50, obs 20/20/10.
        n = 50
        p = (2 * 20 + 20) / (2 * n)  # allele freq = 0.6
        exp_het = 2 * p * (1 - p) * n  # 24
        exp_hom_ref = p**2 * n  # 18
        exp_hom_alt = (1 - p) ** 2 * n  # 8
        expected_chi2 = (
            (20 - exp_hom_ref) ** 2 / exp_hom_ref
            + (20 - exp_het) ** 2 / exp_het
            + (10 - exp_hom_alt) ** 2 / exp_hom_alt
        )
        chi2, _ = hwe.hwe_chi2_test(20, 20, 10)
        assert abs(chi2 - expected_chi2) < 1e-9

    def test_p_value_consistency_with_chi2(self) -> None:
        # p-value must equal the chi2(1) survival function of the statistic
        chi2, p = hwe.hwe_chi2_test(10, 80, 10)
        # chi2(1) sf = erfc(sqrt(chi2/2))
        expected_p = math.erfc(math.sqrt(chi2 / 2.0))
        assert abs(p - expected_p) < 1e-9


class TestHweFlagVariants:
    @staticmethod
    def _vcf(variant_specs: dict) -> dict:
        """Build a parse_vcf_full-shaped dict (SAMPLE-major genotypes list).

        variant_specs: {variant_index: [dosage per sample]} — converted to
        the sample-major layout parse_vcf_full produces.
        """
        n_variants = (max(variant_specs) + 1) if variant_specs else 0
        n_samples = max((len(v) for v in variant_specs.values()), default=0)
        sample_rows = []
        for s in range(n_samples):
            row = []
            for v in range(n_variants):
                col = variant_specs.get(v, [])
                row.append(col[s] if s < len(col) else 0)
            sample_rows.append(row)
        variants = [{"id": f"var_{i}", "chrom": "NC_037638.1", "pos": 1000 * (i + 1)} for i in range(n_variants)]
        return {
            "variants": variants,
            "genotypes": sample_rows,
            "samples": [f"S{i}" for i in range(n_samples)],
        }

    def test_flags_heterozygote_excess_variant(self) -> None:
        # Variant 0: 0/30/0 → het excess; Variant 1: 10/20/10 HWE
        vcf = self._vcf({0: [1] * 30, 1: [0] * 10 + [1] * 10 + [2] * 10})
        flagged = hwe.hwe_flag_variants(vcf, threshold=1e-6)
        assert len(flagged) == 2
        assert flagged[0]["flagged"] is True
        assert flagged[0]["n_het"] == 30
        assert flagged[1]["flagged"] is False
        assert flagged[1]["n_hom_ref"] == 10

    def test_threshold_sensitivity(self) -> None:
        vcf = self._vcf({0: [1] * 30})
        strict = hwe.hwe_flag_variants(vcf, threshold=1e-300)
        assert strict[0]["flagged"] is False  # p is tiny but not < 1e-300
        loose = hwe.hwe_flag_variants(vcf, threshold=1e-6)
        assert loose[0]["flagged"] is True

    def test_result_fields_present(self) -> None:
        vcf = self._vcf({0: [0, 1, 2, 1, 0]})
        row = hwe.hwe_flag_variants(vcf)[0]
        for key in ("chrom", "pos", "snp_id", "maf", "chi2_hwe", "p_hwe", "flagged", "obs_het", "exp_het"):
            assert key in row
        assert row["obs_het"] == pytest.approx(0.4)  # 2 het of 5 called

    def test_missing_genotypes_ignored(self) -> None:
        # -1 codes missing
        vcf = self._vcf({0: [0, 1, 2, -1, -1]})
        row = hwe.hwe_flag_variants(vcf)[0]
        assert row["n_hom_ref"] == 1
        assert row["n_het"] == 1
        assert row["n_hom_alt"] == 1
        assert row["n_missing"] == 2

    def test_empty_vcf(self) -> None:
        assert hwe.hwe_flag_variants(self._vcf({})) == []
