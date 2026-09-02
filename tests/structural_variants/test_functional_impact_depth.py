"""Depth tests for structural-variant functional-impact annotation.

Covers dosage-sensitivity classification, TAD disruption per SV type,
aggregate pathogenicity scoring bounds and monotonicity, and the
impact-classification decision table with real annotation dictionaries.
"""

from __future__ import annotations

from metainformant.structural_variants.annotation.functional_impact import (
    _classify_impact,
    _count_coding_genes,
    _find_overlapping_genes,
    assess_dosage_sensitivity,
    predict_tad_disruption,
    score_pathogenicity,
)


class TestAssessDosageSensitivity:
    def test_unknown_gene_defaults(self) -> None:
        ds = assess_dosage_sensitivity("GENE_UNKNOWN", {}, {}, {}, {})
        assert ds.gene_name == "GENE_UNKNOWN"
        assert ds.haploinsufficiency_score == 0.0
        assert ds.triplosensitivity_score == 0.0
        assert ds.pli_score == 0.0
        assert ds.loeuf == 1.0
        assert not ds.is_haploinsufficient
        assert not ds.is_triplosensitive

    def test_hi_score_threshold(self) -> None:
        ds = assess_dosage_sensitivity("G", haploinsufficiency_db={"G": 0.6})
        assert ds.is_haploinsufficient
        below = assess_dosage_sensitivity("G", haploinsufficiency_db={"G": 0.4})
        assert not below.is_haploinsufficient

    def test_pli_threshold(self) -> None:
        ds = assess_dosage_sensitivity("G", pli_db={"G": 0.95})
        assert ds.is_haploinsufficient

    def test_loeuf_threshold(self) -> None:
        ds = assess_dosage_sensitivity("G", loeuf_db={"G": 0.2})
        assert ds.is_haploinsufficient

    def test_triplosensitivity_threshold(self) -> None:
        ds = assess_dosage_sensitivity("G", triplosensitivity_db={"G": 0.7})
        assert ds.is_triplosensitive
        below = assess_dosage_sensitivity("G", triplosensitivity_db={"G": 0.3})
        assert not below.is_triplosensitive

    def test_no_databases(self) -> None:
        ds = assess_dosage_sensitivity("G")
        assert ds.haploinsufficiency_score == 0.0 and not ds.is_haploinsufficient


class TestPredictTadDisruption:
    BOUNDARY = {"chrom": "chr1", "start": 99_000, "end": 101_000}  # midpoint 100_000

    def test_deletion_spanning_boundary(self) -> None:
        pred = predict_tad_disruption(
            {"chrom": "chr1", "start": 50_000, "end": 150_000, "sv_type": "DEL"}, [self.BOUNDARY]
        )
        assert pred.n_boundaries_disrupted == 1
        assert pred.disruption_score > 0.0

    def test_deletion_breakpoint_near_boundary(self) -> None:
        pred = predict_tad_disruption(
            {"chrom": "chr1", "start": 95_000, "end": 96_000, "sv_type": "DEL"}, [self.BOUNDARY]
        )
        assert pred.n_boundaries_disrupted == 1

    def test_deletion_far_away(self) -> None:
        pred = predict_tad_disruption(
            {"chrom": "chr1", "start": 500_000, "end": 510_000, "sv_type": "DEL"}, [self.BOUNDARY]
        )
        assert pred.n_boundaries_disrupted == 0
        assert pred.disruption_score == 0.0

    def test_dup_only_disrupts_when_spanning(self) -> None:
        spanning = predict_tad_disruption(
            {"chrom": "chr1", "start": 90_000, "end": 110_000, "sv_type": "DUP"}, [self.BOUNDARY]
        )
        nearby = predict_tad_disruption(
            {"chrom": "chr1", "start": 95_000, "end": 96_000, "sv_type": "DUP"}, [self.BOUNDARY]
        )
        assert spanning.n_boundaries_disrupted == 1
        assert nearby.n_boundaries_disrupted == 0

    def test_inversion_near_boundary(self) -> None:
        pred = predict_tad_disruption(
            {"chrom": "chr1", "start": 99_995, "end": 200_000, "sv_type": "INV"}, [self.BOUNDARY]
        )
        assert pred.n_boundaries_disrupted == 1

    def test_other_chromosome_ignored(self) -> None:
        pred = predict_tad_disruption(
            {"chrom": "chr2", "start": 99_000, "end": 101_000, "sv_type": "DEL"}, [self.BOUNDARY]
        )
        assert pred.n_boundaries_disrupted == 0

    def test_boundary_genes_collected(self) -> None:
        boundary = dict(self.BOUNDARY, genes=["BRCA1", "TP53"])
        pred = predict_tad_disruption({"chrom": "chr1", "start": 50_000, "end": 150_000, "sv_type": "DEL"}, [boundary])
        assert set(pred.genes_in_affected_tads) == {"BRCA1", "TP53"}

    def test_score_scales_with_sv_type_weight(self) -> None:
        del_pred = predict_tad_disruption(
            {"chrom": "chr1", "start": 50_000, "end": 150_000, "sv_type": "DEL"}, [self.BOUNDARY]
        )
        dup_pred = predict_tad_disruption(
            {"chrom": "chr1", "start": 50_000, "end": 150_000, "sv_type": "DUP"}, [self.BOUNDARY]
        )
        assert del_pred.disruption_score > dup_pred.disruption_score


class TestScorePathogenicity:
    def _base_annotations(self) -> dict:
        return {
            "overlapping_genes": [],
            "dosage_sensitive_genes": [],
            "tad_disrupted": False,
            "impact_level": "MODIFIER",
            "n_coding_genes": 0,
        }

    def test_score_within_bounds_and_deterministic(self) -> None:
        variant = {"chrom": "chr1", "start": 1, "end": 200_000, "sv_type": "DEL"}
        s1 = score_pathogenicity(variant, self._base_annotations())
        s2 = score_pathogenicity(variant, self._base_annotations())
        assert 0.0 <= s1 <= 1.0
        assert s1 == s2

    def test_larger_variants_score_higher(self) -> None:
        small = score_pathogenicity(
            {"chrom": "chr1", "start": 1, "end": 1_000, "sv_type": "DEL"}, self._base_annotations()
        )
        large = score_pathogenicity(
            {"chrom": "chr1", "start": 1, "end": 500_000, "sv_type": "DEL"}, self._base_annotations()
        )
        assert large > small

    def test_dosage_genes_increase_score(self) -> None:
        base = self._base_annotations()
        with_dosage = dict(base, dosage_sensitive_genes=["G1"], n_coding_genes=2)
        s_base = score_pathogenicity({"chrom": "chr1", "start": 1, "end": 100_000, "sv_type": "DEL"}, base)
        s_dosage = score_pathogenicity({"chrom": "chr1", "start": 1, "end": 100_000, "sv_type": "DEL"}, with_dosage)
        assert s_dosage > s_base

    def test_common_frequency_reduces_score(self) -> None:
        base = self._base_annotations()
        absent = base
        rare = dict(base, population_frequency=0.0001)
        common = dict(base, population_frequency=0.4)
        variant = {"chrom": "chr1", "start": 1, "end": 100_000, "sv_type": "DEL"}
        s_absent = score_pathogenicity(variant, absent)
        s_rare = score_pathogenicity(variant, rare)
        s_common = score_pathogenicity(variant, common)
        assert s_absent > s_rare > s_common

    def test_frequency_component_monotone_and_bounded(self) -> None:
        """Regression: the AF term must decrease monotonically with frequency.

        The previous formula had an inverted sign, so rare variants scored
        LOWER than common ones - the opposite of the documented 'rare = more
        pathogenic' behavior.
        """
        base = self._base_annotations()
        variant = {"chrom": "chr1", "start": 1, "end": 100_000, "sv_type": "DEL"}
        freqs = [1e-6, 1e-4, 1e-2, 0.1, 0.4]
        scores = [score_pathogenicity(variant, dict(base, population_frequency=f)) for f in freqs]
        assert scores == sorted(scores, reverse=True)
        assert all(0.0 <= s <= 1.0 for s in scores)

    def test_impact_level_ordering(self) -> None:
        variant = {"chrom": "chr1", "start": 1, "end": 100_000, "sv_type": "DEL"}
        high = score_pathogenicity(variant, dict(self._base_annotations(), impact_level="HIGH"))
        modifier = score_pathogenicity(variant, dict(self._base_annotations(), impact_level="MODIFIER"))
        assert high > modifier


class TestClassifyImpact:
    GENES = [{"name": "A", "is_coding": True}, {"name": "B", "is_coding": False}]

    def test_dosage_loss_high(self) -> None:
        assert _classify_impact("DEL", 10_000, 1, 1, False, self.GENES, ["A"]) == ("HIGH", "dosage_loss")

    def test_dosage_gain_high(self) -> None:
        assert _classify_impact("DUP", 10_000, 1, 1, False, self.GENES, ["A"]) == ("HIGH", "dosage_gain")

    def test_coding_gene_disruption_high(self) -> None:
        assert _classify_impact("DEL", 10_000, 2, 0, False, self.GENES, ["A", "B"]) == ("HIGH", "gene_disruption")

    def test_translocation_gene_fusion(self) -> None:
        assert _classify_impact("TRA", 10_000, 1, 0, False, self.GENES, ["A"]) == ("HIGH", "gene_fusion_candidate")

    def test_dup_without_dosage_moderate(self) -> None:
        assert _classify_impact("DUP", 10_000, 1, 0, False, self.GENES, ["A"]) == ("MODERATE", "gene_duplication")

    def test_inv_moderate(self) -> None:
        assert _classify_impact("INV", 10_000, 1, 0, False, self.GENES, ["A"]) == ("MODERATE", "gene_inversion")

    def test_tad_disruption_moderate(self) -> None:
        assert _classify_impact("INS", 10_000, 0, 0, True, self.GENES, []) == ("MODERATE", "tad_disruption")

    def test_large_intergenic_low(self) -> None:
        assert _classify_impact("DEL", 200_000, 0, 0, False, self.GENES, []) == ("LOW", "large_intergenic")

    def test_intergenic_low(self) -> None:
        assert _classify_impact("DEL", 20_000, 0, 0, False, self.GENES, []) == ("LOW", "intergenic")

    def test_small_intergenic_modifier(self) -> None:
        assert _classify_impact("DEL", 5_000, 0, 0, False, self.GENES, []) == ("MODIFIER", "benign")


class TestGeneHelpers:
    GENES = [
        {"name": "A", "is_coding": True},
        {"name": "B", "is_coding": False},
        {"name": "C", "is_coding": True},
    ]

    def test_find_overlapping_genes(self) -> None:
        found = _find_overlapping_genes(
            "chr1",
            100,
            200,
            [
                {"name": "A", "chrom": "chr1", "start": 150, "end": 300},
                {"name": "B", "chrom": "chr1", "start": 500, "end": 600},
            ],
        )
        assert found == ["A"]

    def test_count_coding_genes(self) -> None:
        assert _count_coding_genes(["A", "B"], self.GENES) == 1
        assert _count_coding_genes(["A", "C"], self.GENES) == 2
        assert _count_coding_genes(["B"], self.GENES) == 0
