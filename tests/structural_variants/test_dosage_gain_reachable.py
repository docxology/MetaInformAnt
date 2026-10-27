"""Regression tests: ``predict_functional_impact`` must use every dosage-sensitivity database.

``assess_dosage_sensitivity`` classifies a gene as triplosensitive from a triplosensitivity database (and
as haploinsufficient from HI, pLI or LOEUF values), but ``predict_functional_impact`` only had a parameter
for the haploinsufficiency database and passed nothing else on. A duplication could therefore never be
reported as ``dosage_gain``: the class was unreachable.
"""

from __future__ import annotations

from metainformant.structural_variants.annotation.functional_impact import predict_functional_impact

GENES = [
    {"chrom": "chr1", "start": 1_000, "end": 5_000, "name": "TS1", "is_coding": True},
    {"chrom": "chr1", "start": 6_000, "end": 9_000, "name": "OTHER", "is_coding": True},
]
DUP_OVER_TS1 = {"chrom": "chr1", "start": 500, "end": 5_500, "sv_type": "DUP"}
DEL_OVER_TS1 = {"chrom": "chr1", "start": 500, "end": 5_500, "sv_type": "DEL"}


class TestTriplosensitivity:
    def test_duplication_over_a_triplosensitive_gene_is_a_dosage_gain(self) -> None:
        impact = predict_functional_impact(DUP_OVER_TS1, GENES, triplosensitivity_db={"TS1": 0.9})

        assert impact.affected_genes == ["TS1"]
        assert impact.dosage_sensitive_genes == ["TS1"]
        assert (impact.impact_level, impact.impact_type) == ("HIGH", "dosage_gain")

    def test_without_the_database_a_duplication_is_only_a_gene_duplication(self) -> None:
        impact = predict_functional_impact(DUP_OVER_TS1, GENES)
        assert impact.dosage_sensitive_genes == []
        assert (impact.impact_level, impact.impact_type) == ("MODERATE", "gene_duplication")

    def test_the_dosage_gene_raises_the_pathogenicity_score(self) -> None:
        with_db = predict_functional_impact(DUP_OVER_TS1, GENES, triplosensitivity_db={"TS1": 0.9})
        without_db = predict_functional_impact(DUP_OVER_TS1, GENES)
        assert with_db.pathogenicity_score > without_db.pathogenicity_score

    def test_a_low_triplosensitivity_score_is_not_a_dosage_gain(self) -> None:
        impact = predict_functional_impact(DUP_OVER_TS1, GENES, triplosensitivity_db={"TS1": 0.3})
        assert impact.dosage_sensitive_genes == []
        assert impact.impact_type == "gene_duplication"

    def test_a_gene_missing_from_the_database_is_not_sensitive(self) -> None:
        impact = predict_functional_impact(DUP_OVER_TS1, GENES, triplosensitivity_db={"OTHER": 0.9})
        assert impact.dosage_sensitive_genes == []

    def test_triplosensitivity_does_not_matter_for_a_deletion(self) -> None:
        impact = predict_functional_impact(DEL_OVER_TS1, GENES, triplosensitivity_db={"TS1": 0.9})
        assert impact.dosage_sensitive_genes == []
        assert impact.impact_type == "gene_disruption"

    def test_haploinsufficiency_does_not_make_a_duplication_a_dosage_gain(self) -> None:
        impact = predict_functional_impact(DUP_OVER_TS1, GENES, haploinsufficiency_db={"TS1": 0.9})
        assert impact.dosage_sensitive_genes == []
        assert impact.impact_type == "gene_duplication"

    def test_precomputed_gene_overlaps_use_the_database_too(self) -> None:
        variant = dict(DUP_OVER_TS1, overlapping_genes=["TS1"])
        impact = predict_functional_impact(variant, [], triplosensitivity_db={"TS1": 0.9})
        assert (impact.impact_level, impact.impact_type) == ("HIGH", "dosage_gain")


class TestOtherDosageDatabases:
    def test_a_constrained_gene_by_pli_makes_a_deletion_a_dosage_loss(self) -> None:
        impact = predict_functional_impact(DEL_OVER_TS1, GENES, pli_db={"TS1": 0.95})
        assert impact.dosage_sensitive_genes == ["TS1"]
        assert (impact.impact_level, impact.impact_type) == ("HIGH", "dosage_loss")

    def test_a_constrained_gene_by_loeuf_makes_a_deletion_a_dosage_loss(self) -> None:
        impact = predict_functional_impact(DEL_OVER_TS1, GENES, loeuf_db={"TS1": 0.2})
        assert (impact.impact_level, impact.impact_type) == ("HIGH", "dosage_loss")

    def test_without_them_the_same_deletion_is_a_plain_gene_disruption(self) -> None:
        impact = predict_functional_impact(DEL_OVER_TS1, GENES, pli_db={"TS1": 0.5}, loeuf_db={"TS1": 0.9})
        assert impact.dosage_sensitive_genes == []
        assert impact.impact_type == "gene_disruption"

    def test_haploinsufficiency_database_still_works(self) -> None:
        impact = predict_functional_impact(DEL_OVER_TS1, GENES, haploinsufficiency_db={"TS1": 0.8})
        assert (impact.impact_level, impact.impact_type) == ("HIGH", "dosage_loss")

    def test_databases_combine(self) -> None:
        impact = predict_functional_impact(
            {"chrom": "chr1", "start": 500, "end": 9_500, "sv_type": "DUP"},
            GENES,
            haploinsufficiency_db={"OTHER": 0.9},
            triplosensitivity_db={"TS1": 0.9},
        )
        assert impact.affected_genes == ["TS1", "OTHER"]
        assert impact.dosage_sensitive_genes == ["TS1"]  # OTHER is haploinsufficient, which a duplication ignores


def test_existing_positional_arguments_keep_their_meaning() -> None:
    """The fourth positional argument is still the TAD boundaries (the new databases come after it)."""
    boundary = {"chrom": "chr1", "start": 2_999, "end": 3_001}
    variant = {"chrom": "chr1", "start": 500, "end": 5_500, "sv_type": "DEL"}
    impact = predict_functional_impact(variant, GENES, {"TS1": 0.8}, [boundary])
    assert impact.tad_disrupted is True
    assert impact.dosage_sensitive_genes == ["TS1"]
