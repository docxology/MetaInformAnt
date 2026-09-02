"""Zero-mocks tests for gwas.analysis.strain_analysis (Fst, private variants)."""

from __future__ import annotations

import math
import random

import pytest

from metainformant.gwas.analysis import strain_analysis as sa


class TestBuildStrainMap:
    def test_groups_by_leading_letter(self) -> None:
        ids = ["C15ITQ", "C16G", "I12ITW", "M03WORK", "R21G"]
        smap = sa.build_strain_map(ids)
        assert smap == {"C": [0, 1], "I": [2], "M": [3], "R": [4]}


class TestAlleleFrequenciesByStrain:
    def test_known_frequencies(self) -> None:
        ids = ["C1G", "C2G", "M1G", "M2G"]
        geno = [[0, 2, 1, 1], [2, 0, 0, 0]]
        af = sa.compute_allele_frequencies_by_strain(geno, ids)
        # C strain dosages [0,2]: af = 2/(2*2) = 0.5 ; M strain [1,1]: af = 2/4 = 0.5
        assert af["C"][0] == pytest.approx(0.5)
        assert af["M"][0] == pytest.approx(0.5)
        # Variant 2: C [2,0] → 0.5; M [0,0] → 0.0
        assert af["C"][1] == pytest.approx(0.5)
        assert af["M"][1] == pytest.approx(0.0)

    def test_missing_dosages_excluded(self) -> None:
        ids = ["C1G", "M1G"]
        geno = [[-1, 2]]
        af = sa.compute_allele_frequencies_by_strain(geno, ids)
        assert math.isnan(af["C"][0])
        assert af["M"][0] == pytest.approx(1.0)


class TestComputeFst:
    def test_fixed_difference_max_fst(self) -> None:
        # Strain C all alt, strain M all ref → Fst = 1
        ids = ["C1G", "C2G", "M1G", "M2G"]
        geno = [[2, 2, 0, 0]]
        fst = sa.compute_fst_per_variant(geno, ids)
        vals = fst["C_vs_M"]
        assert len(vals) == 1
        assert vals[0] == pytest.approx(1.0)

    def test_identical_frequencies_zero_fst(self) -> None:
        ids = ["C1G", "C2G", "M1G", "M2G"]
        geno = [[1, 1, 1, 1]]
        fst = sa.compute_fst_per_variant(geno, ids)
        assert fst["C_vs_M"][0] == pytest.approx(0.0, abs=1e-9)

    def test_hudson_vs_weir_cockerham_ordering(self) -> None:
        # Moderately differentiated: both estimators positive, same sign
        ids = ["C1G", "C2G", "C3G", "M1G", "M2G", "M3G"]
        geno = [[2, 2, 1, 0, 0, 0]]
        wc = sa.compute_fst_per_variant(geno, ids, method="weir_cockerham")["C_vs_M"][0]
        hud = sa.compute_fst_per_variant(geno, ids, method="hudson")["C_vs_M"][0]
        assert 0.0 < wc <= 1.0
        assert 0.0 < hud <= 1.0

    def test_hudson_closed_form(self) -> None:
        # p1=0.5, p2=0.0, n1=n2=1 → hand-computed Hudson Fst
        fst = sa._hudson_fst(0.5, 0.0, 1, 1)
        p_mean = 0.25
        h_total = 2 * p_mean * (1 - p_mean)
        h_within = (2 * 0.25 + 0.0) / 2
        assert fst == pytest.approx((h_total - h_within) / h_total)

    def test_all_pairs_generated(self) -> None:
        ids = ["C1G", "I1G", "M1G", "R1G"]
        geno = [[1, 1, 1, 1]]
        fst = sa.compute_fst_per_variant(geno, ids)
        assert set(fst.keys()) == {"C_vs_I", "C_vs_M", "C_vs_R", "I_vs_M", "I_vs_R", "M_vs_R"}

    def test_strain_pair_subset(self) -> None:
        ids = ["C1G", "I1G", "M1G"]
        geno = [[1, 1, 1]]
        fst = sa.compute_fst_per_variant(geno, ids, strain_pairs=[("C", "M")])
        assert set(fst.keys()) == {"C_vs_M"}


class TestStrainSpecificVariants:
    def test_private_variant_detected(self) -> None:
        # Variant 0: fixed alt in C, fixed ref in M → strain-differentiating;
        # both strains report it (each carries its own private allele at the
        # locus). Variant 1: polymorphic everywhere → private to nobody.
        ids = ["C1G", "C2G", "M1G", "M2G"]
        geno = [[2, 2, 0, 0], [1, 0, 1, 0]]
        private = sa.strain_specific_variants(geno, ids, maf_threshold=0.05)
        assert 0 in private["C"]
        assert 0 in private["M"]  # mirror image: ref allele private to M
        assert 1 not in private["C"] and 1 not in private["M"]

    def test_asymmetric_three_strain_private(self) -> None:
        # Variant 0: polymorphic in C only (M, R monomorphic ref)
        # Variant 1: fixed alt in R, monomorphic ref elsewhere → R-private
        ids = ["C1G", "C2G", "M1G", "M2G", "R1G", "R2G"]
        geno = [
            [1, 2, 0, 0, 0, 0],
            [0, 0, 0, 0, 2, 2],
        ]
        private = sa.strain_specific_variants(geno, ids, maf_threshold=0.05)
        assert private["C"] == [0]
        assert private["R"] == [1]
        assert private["M"] == []


class TestComputeGlobalFst:
    def test_global_summary_structure(self) -> None:
        ids = ["C1G", "C2G", "M1G", "M2G", "R1G", "R2G"]
        rng = random.Random(19)
        geno = [[rng.choice([0, 1, 2]) for _ in ids] for _ in range(20)]
        summary = sa.compute_global_fst(geno, ids)
        assert "global_mean_fst" in summary
        assert set(summary["pairwise"].keys()) == {"C_vs_M", "C_vs_R", "M_vs_R"}
        for stats in summary["pairwise"].values():
            assert 0.0 <= stats["mean_fst"] <= 1.0
            assert stats["n_variants"] == 20
            assert stats["max_fst"] >= stats["mean_fst"] - 1e-9
