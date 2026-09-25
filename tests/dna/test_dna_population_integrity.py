"""Integrity tests for dna.population statistics.

Pins the hand-computed behavior of the McDonald-Kreitman test, the signed
linkage-disequilibrium coefficient, and Hudson's F_ST estimator. All
implementations are exercised for real (no stubs, no mocks); every expected
value is derived in the comments.
"""

from __future__ import annotations

import pytest

from metainformant.dna.population import analysis, core


class TestLinkageDisequilibriumSign:
    """D = f_AB - p_A p_B with deterministic reference alleles."""

    def test_coupling_is_positive(self):
        # Haplotypes (A,C) x2 and (G,T) x2: the reference alleles (majority,
        # ties broken alphabetically -> A and C) co-occur in half the
        # haplotypes. p_A = p_B = 0.5, f_AB = 0.5 -> D = 0.5 - 0.25 = +0.25.
        assert core.linkage_disequilibrium(
            ["AC", "AC", "GT", "GT"], 0, 1
        ) == pytest.approx(0.25)

    def test_strong_coupling_is_positive(self):
        # Haplotypes (A,T) x3 and (G,C) x1: reference alleles A and T
        # (3/4 each) co-occur in 3/4 of haplotypes.
        # D = 0.75 - 0.75*0.75 = +0.1875.
        assert core.linkage_disequilibrium(
            ["AT", "AT", "AT", "GC"], 0, 1
        ) == pytest.approx(0.1875)

    def test_equal_frequency_repulsion_is_negative(self):
        # Haplotypes (A,T) x2 and (G,C) x2: reference alleles A and C never
        # co-occur -> D = 0 - 0.25 = -0.25. The previous implementation
        # picked reference alleles via set iteration order and could return
        # +0.25 for the same input.
        assert core.linkage_disequilibrium(
            ["AT", "AT", "GC", "GC"], 0, 1
        ) == pytest.approx(-0.25)

    def test_moderate_repulsion_is_negative(self):
        # Haplotypes (A,C) x1, (A,T) x2, (G,T) x1: reference alleles A (3/4)
        # and T (3/4); f_(A,T) = 2/4 -> D = 0.5 - 0.75*0.75 = -0.0625.
        assert core.linkage_disequilibrium(
            ["AC", "AT", "AT", "GT"], 0, 1
        ) == pytest.approx(-0.0625)

    def test_tie_break_is_independent_of_input_order(self):
        # Equal frequencies at both sites: the alphabetically-first allele
        # must be the reference regardless of how the sequences are ordered.
        first = core.linkage_disequilibrium(["AT", "TA", "AT", "TA"], 0, 1)
        second = core.linkage_disequilibrium(["TA", "AT", "TA", "AT"], 0, 1)
        assert first == second == pytest.approx(-0.25)


class TestMcdonaldKreitman:
    """Real 2x2 MK contingency with hand-computed tables."""

    def test_hand_computed_neutral_table(self):
        # 4 ingroup sequences + 1 outgroup, six codon columns:
        #   1 ATG vs ATA  -> nonsynonymous divergence (Dn = 1)  [M -> I]
        #   2 GCA vs GCG  -> synonymous divergence   (Ds = 1)  [A -> A]
        #   3 AAA/AAT mix -> nonsynonymous polymorphism (Pn = 1) [K/N]
        #   4 GGT/GGG mix -> synonymous polymorphism   (Ps = 1) [G/G]
        #   5 TTT == TTT  -> no count
        #   6 CAA == CAA  -> no count
        s1 = "ATGGCATTTAAACAAGGT"
        s2 = "ATGGCATTTAAACAAGGG"
        s3 = "ATGGCATTTAAACAAGGT"
        s4 = "ATGGCATTTAATCAAGGT"
        outgroup = "ATAGCGTTTAAACAAGGT"
        res = analysis.mcdonald_kreitman_contingency([s1, s2, s3, s4, outgroup])
        assert (res["Ps"], res["Pn"], res["Ds"], res["Dn"]) == (1, 1, 1, 1)
        # NI = (Pn/Ps)/(Dn/Ds) = 1 -> alpha = 0, omega = Dn/Ds = 1.
        assert res["neutral_ratio"] == pytest.approx(1.0)
        assert res["alpha"] == pytest.approx(0.0)
        assert res["omega"] == pytest.approx(1.0)
        # Fisher exact two-sided p on [[1, 1], [1, 1]] is 1.0.
        assert res["fisher_p"] == pytest.approx(1.0)

    def test_hand_computed_adaptive_table(self):
        # Six codon columns giving [[Pn, Ps], [Dn, Ds]] = [[1, 2], [2, 1]]:
        #   1 ATG vs ATA  -> nonsynonymous divergence (Dn)   [M -> I]
        #   2 GCA vs GCG  -> synonymous divergence   (Ds)    [A -> A]
        #   3 TTT vs TTA  -> nonsynonymous divergence (Dn)   [F -> L]
        #   4 AAA/AAT mix vs AAA -> nonsynonymous polymorphism (Pn) [K/N]
        #   5 CAA/CAG mix vs CAA -> synonymous polymorphism    (Ps) [Q/Q]
        #   6 GGT/GGG mix vs GGT -> synonymous polymorphism    (Ps) [G/G]
        # NI = (1/2)/(2/1) = 0.25 -> alpha = 0.75; omega = Dn/Ds = 2.0.
        # The margins (3, 3; 3, 3) admit four tables with hypergeometric
        # probabilities 0.05, 0.45, 0.45, 0.05. Every probability is <= the
        # observed 0.45, so the two-sided Fisher p saturates at 1.0.
        seqs = [
            "ATGGCATTTAAACAAGGT",
            "ATGGCATTTAAACAAGGG",
            "ATGGCATTTAAACAGGGT",
            "ATGGCATTTAATCAAGGT",
            "ATAGCGTTAAAACAAGGT",
        ]
        res = analysis.mcdonald_kreitman_contingency(seqs)
        assert (res["Ps"], res["Pn"], res["Ds"], res["Dn"]) == (2, 1, 1, 2)
        assert res["neutral_ratio"] == pytest.approx(0.25)
        assert res["alpha"] == pytest.approx(0.75)
        assert res["omega"] == pytest.approx(2.0)
        assert res["fisher_p"] == pytest.approx(1.0)
        assert analysis.mcdonald_kreitman_test(seqs) == (0.75, 2.0)

    def test_detect_selection_mk_method_uses_real_test(self):
        results = analysis.detect_selection(
            [
                "ATGGCATTTAAACAAGGT",
                "ATGGCATTTAAACAAGGG",
                "ATGGCATTTAAACAGGGT",
                "ATGGCATTTAATCAAGGT",
                "ATAGCGTTAAAACAAGGT",
            ],
            method="mk_test",
        )
        assert results["alpha"] == pytest.approx(0.75)
        assert results["omega"] == pytest.approx(2.0)

    def test_ambiguous_codon_column_is_skipped(self):
        # The single codon column contains N in the ingroup and cannot be
        # classified, so nothing is counted as divergence either.
        res = analysis.mcdonald_kreitman_contingency(["NNG", "NNG", "GCA"])
        assert (res["Ps"], res["Pn"], res["Ds"], res["Dn"]) == (0, 0, 0, 0)
        assert res["fisher_p"] == pytest.approx(1.0)

    def test_misaligned_sequences_raise(self):
        with pytest.raises(ValueError, match="aligned"):
            analysis.mcdonald_kreitman_test(["ATGGCC", "ATGGCCAA", "ATGGCC"])

    def test_too_few_sequences_returns_neutral_defaults(self):
        assert analysis.mcdonald_kreitman_test(["ATGGCC", "ATGGCC"]) == (0.0, 0.0)
        assert analysis.mcdonald_kreitman_test(["ATGGCC"]) == (0.0, 0.0)
        assert analysis.mcdonald_kreitman_test([]) == (0.0, 0.0)

    def test_fisher_exact_fallback_matches_scipy_definition(self):
        # Two-sided p = sum of hypergeometric probabilities <= the observed
        # table's probability (scipy's definition). For margins (2, 3; 2, 3)
        # the probabilities are 0.3 (observed, k=0), 0.6, 0.1 -> p = 0.4.
        p = analysis._fisher_exact_two_sided_2x2(0, 2, 2, 1)
        assert p == pytest.approx(0.4)
        # Degenerate margin: nothing to test.
        assert analysis._fisher_exact_two_sided_2x2(0, 2, 0, 2) == pytest.approx(1.0)


class TestHudsonFst:
    """Hudson's moment estimator, hand-computed on tiny populations."""

    def test_hand_computed_value(self):
        # Sequences are haploid alignments: one allele per sequence per site.
        # pop1 = A,A,G (n1 = 3), pop2 = G,G (n2 = 2); the pooled reference
        # allele is the most common with alphabetical tie-break.
        #   site 0: ref G (pooled G majority); p1 = 1/3, p2 = 1
        #     num = (1/3-1)^2 - (1/3*2/3)/(3-1) - 0 = 4/9 - 1/9 = 1/3
        #     den = (1/3)*(1-1) + 1*(1-1/3)         = 2/3  -> 0.5
        #   site 1: ref C (T,C tie broken alphabetically); p1 = 1/3, p2 = 1
        #     num = 1/3, den = 2/3                  -> 0.5
        # Fst = (1/3 + 1/3) / (2/3 + 2/3) = 0.5.
        assert core.hudson_fst(["AT", "AT", "GC"], ["GC", "GC"]) == pytest.approx(0.5)

    def test_fixed_difference_is_one(self):
        # Single haploid per population: the sampling-corrected estimator is
        # undefined (n < 2 at every site), but the populations are fixed for
        # different alleles at every comparable site -> the documented
        # degenerate fallback returns 1.0.
        assert core.hudson_fst(["AAAA"], ["GGGG"]) == pytest.approx(1.0)

    def test_identically_fixed_is_zero(self):
        assert core.hudson_fst(["AAAA"], ["AAAA"]) == pytest.approx(0.0)

    def test_symmetric_in_population_order(self):
        pop1 = ["AA", "AA"]
        pop2 = ["AG", "GG"]
        assert core.hudson_fst(pop1, pop2) == pytest.approx(core.hudson_fst(pop2, pop1))

    def test_value_differs_from_pooled_heterozygosity_estimator(self):
        # Same data as the hand-computed pin above: pop1 = A,A,G (n1 = 3),
        # pop2 = G,G (n2 = 2), two polymorphic sites.
        #   Hudson (this module's core.hudson_fst): weights site numerators
        #   by the between-population denominator and applies the sampling
        #   correction -> (1/3 + 1/3) / (2/3 + 2/3) = 0.5.
        #   G_ST (analysis.calculate_fst): unweighted mean over polymorphic
        #   sites of (H_T - H_S) / H_T; per site H_T = 12/25, H_S = 2/9
        #   -> 29/54 each -> 29/54.
        # 0.5 != 29/54: the two functions are distinct estimators.
        pop1 = ["AT", "AT", "GC"]
        pop2 = ["GC", "GC"]
        assert core.hudson_fst(pop1, pop2) == pytest.approx(0.5)
        assert analysis.calculate_fst(pop1, pop2) == pytest.approx(29.0 / 54.0)

    def test_rejects_incompatible_input(self):
        with pytest.raises(ValueError):
            core.hudson_fst([], ["AAAA"])
        with pytest.raises(ValueError):
            core.hudson_fst(["AAAA"], ["AAA"])
