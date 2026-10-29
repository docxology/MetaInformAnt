"""Genes, alleles and diplotypes CPIC has no answer for: they raise or are indeterminate, never normal.

Four rules, each pinned by tests that fail on the code before them:

1. A gene CPIC defines no phenotype for raises ``ValueError`` naming the gene and listing the
   supported genes. It used to be classified with the CYP2D6 score bands. CYP1A2 has activity
   values (it can be scored) but no CPIC phenotype (it cannot be classified).
2. An allele with no activity value is uncharacterized. CPIC calls a diplotype that carries one
   indeterminate, so it is ``MetabolizerPhenotype.INDETERMINATE`` with an activity score of
   ``None`` and the allele named. It used to score 1.0 and read as normal. The ``alleles`` and
   ``metabolism`` APIs say the same thing.
3. ``Diplotype.activity_score`` defaults to ``None`` (not scored), not 0.0, which read as a poor
   metabolizer. Classifying an unscored diplotype scores it from its alleles.
4. CYP1A2's activity values have no CPIC source, which the gene's registry entry states.

The expected genes and values are written from the CPIC documents the package cites, not read
back from the package's own tables.

REAL IMPLEMENTATION -- all tests use real implementations.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Callable
from typing import Any

import pytest

from metainformant.pharmacogenomics.alleles.diplotype import (
    _ACTIVITY_SCORE_TABLES,
    Diplotype,
    calculate_activity_score,
    calculate_activity_score_from_alleles,
    determine_diplotype,
    resolve_ambiguous_diplotypes,
)
from metainformant.pharmacogenomics.alleles.phenotype import (
    _PHENOTYPE_THRESHOLDS,
    MetabolizerPhenotype,
    classify_phenotype,
    get_phenotype_thresholds,
    population_phenotype_frequencies,
    predict_metabolizer_status,
)
from metainformant.pharmacogenomics.metabolism.metabolizer_status import (
    classify_metabolizer,
    compute_activity_score,
    default_allele_function_table,
    dose_adjustment,
    predict_metabolizer_status as predict_metabolizer,
)

PM = MetabolizerPhenotype.POOR
IM = MetabolizerPhenotype.INTERMEDIATE
NM = MetabolizerPhenotype.NORMAL
IND = MetabolizerPhenotype.INDETERMINATE

# The genes CPIC publishes a phenotype definition for (CYP2D6, CYP2C19, CYP2C9, CYP3A5, DPYD,
# TPMT, NUDT15, SLCO1B1 and UGT1A1 guidelines), written out rather than read from the package.
CPIC_PHENOTYPE_GENES = ["CYP2C19", "CYP2C9", "CYP2D6", "CYP3A5", "DPYD", "NUDT15", "SLCO1B1", "TPMT", "UGT1A1"]

# Every function that assigns a phenotype, called for a given gene with inputs that are valid
# for a gene CPIC does define.
PHENOTYPE_CALLS: dict[str, Callable[[str], Any]] = {
    "predict_metabolizer_status": lambda gene: predict_metabolizer_status(1.0, gene),
    "get_phenotype_thresholds": get_phenotype_thresholds,
    "classify_phenotype_string": lambda gene: classify_phenotype("*1/*2", gene),
    "classify_phenotype_object": lambda gene: classify_phenotype(Diplotype("*1", "*2", gene, 2.0), gene),
    "population_phenotype_frequencies": lambda gene: population_phenotype_frequencies({"*1": 0.5, "*2": 0.5}, gene),
    "classify_metabolizer": lambda gene: classify_metabolizer(1.0, gene),
    "facade_predict_metabolizer_status": lambda gene: predict_metabolizer(
        {"diplotype": "*1/*2"}, gene, activity_scores={"*1": 1.0, "*2": 0.0}
    ),
}


class TestGenesWithoutAPhenotypeDefinitionRaise:
    """Rule 1: no other gene's bands stand in."""

    @pytest.mark.parametrize("call", sorted(PHENOTYPE_CALLS))
    def test_unknown_gene_raises_naming_it_and_the_supported_genes(self, call: str) -> None:
        with pytest.raises(ValueError, match="NOTAGENE") as caught:
            PHENOTYPE_CALLS[call]("NOTAGENE")
        message = str(caught.value)
        for supported in CPIC_PHENOTYPE_GENES:
            assert supported in message, f"{supported} missing from: {message}"

    @pytest.mark.parametrize("call", sorted(PHENOTYPE_CALLS))
    def test_the_gene_is_named_in_upper_case_whatever_the_case_given(self, call: str) -> None:
        with pytest.raises(ValueError, match="NOTAGENE"):
            PHENOTYPE_CALLS[call]("notagene")

    def test_the_supported_genes_are_the_nine_cpic_defines(self) -> None:
        from metainformant.pharmacogenomics.alleles import phenotype

        assert phenotype.supported_phenotype_genes() == CPIC_PHENOTYPE_GENES
        assert sorted(_PHENOTYPE_THRESHOLDS) == CPIC_PHENOTYPE_GENES

    def test_score_that_cyp2d6_would_call_intermediate_is_not_classified_for_another_gene(self) -> None:
        # The old fallback returned IM for a score of 1.0 and NM for a score of 2.0 whatever the gene
        for score in (0.0, 1.0, 2.0, 3.0):
            with pytest.raises(ValueError, match="CYP1A2"):
                predict_metabolizer_status(score, "CYP1A2")

    @pytest.mark.parametrize("gene", CPIC_PHENOTYPE_GENES)
    def test_supported_genes_still_classify_with_whitespace_and_lower_case(self, gene: str) -> None:
        assert predict_metabolizer_status(2.0, f"  {gene.lower()} ") == predict_metabolizer_status(2.0, gene)
        assert classify_phenotype("*1/*1", f" {gene.lower()}")["gene"] == gene
        assert get_phenotype_thresholds(f" {gene.lower()}") == get_phenotype_thresholds(gene)

    def test_invalid_score_is_still_reported_as_a_score_error(self) -> None:
        with pytest.raises(ValueError, match="non-negative"):
            predict_metabolizer_status(-1.0, "CYP2D6")


class TestCyp1a2CanBeScoredButNotClassified:
    """Rules 1 and 4: CYP1A2 has activity values and no CPIC phenotype."""

    def test_it_is_the_only_gene_with_values_and_no_phenotype_definition(self) -> None:
        assert set(_ACTIVITY_SCORE_TABLES) - set(_PHENOTYPE_THRESHOLDS) == {"CYP1A2"}
        assert set(_PHENOTYPE_THRESHOLDS) <= set(_ACTIVITY_SCORE_TABLES)

    def test_activity_scoring_still_works(self) -> None:
        assert determine_diplotype("*1A", "*1F", "CYP1A2").activity_score == 2.5
        assert determine_diplotype("*1C", "*1K", "CYP1A2").activity_score == 1.0
        assert calculate_activity_score_from_alleles("*1A", "*1A", "CYP1A2") == 2.0
        assert compute_activity_score("*1A/*1F", default_allele_function_table()["CYP1A2"]) == 2.5
        scored = Diplotype("*1A", "*1C", "CYP1A2")
        assert calculate_activity_score(scored, "CYP1A2") == 1.5
        assert scored.activity_score == 1.5

    @pytest.mark.parametrize("call", sorted(PHENOTYPE_CALLS))
    def test_classification_raises_with_the_reason(self, call: str) -> None:
        with pytest.raises(ValueError, match="CYP1A2") as caught:
            PHENOTYPE_CALLS[call]("CYP1A2")
        message = str(caught.value)
        assert "CPIC defines no phenotype" in message
        assert "can be scored" in message
        for supported in CPIC_PHENOTYPE_GENES:
            assert supported in message


class TestUncharacterizedAllelesAreIndeterminate:
    """Rule 2: an allele with no activity value is never normal function."""

    def test_activity_score_is_none_not_one_per_unknown_allele(self) -> None:
        assert calculate_activity_score_from_alleles("*1", "*999", "CYP2D6") is None  # used to be 2.0
        assert calculate_activity_score_from_alleles("*998", "*999", "CYP2D6") is None  # used to be 2.0
        assert determine_diplotype("*1", "*999", "CYP2D6").activity_score is None

    def test_an_unknown_gene_has_no_alleles_to_score(self) -> None:
        assert calculate_activity_score_from_alleles("*1", "*2", "NOTAGENE") is None  # used to be 2.0

    @pytest.mark.parametrize(
        ("diplotype", "unknown"),
        [
            ("*1/*999", ["*999"]),
            ("*999/*1", ["*999"]),
            ("*998/*999", ["*998", "*999"]),
            ("*999/*999", ["*999"]),
            ("999/*4", ["*999"]),
        ],
    )
    def test_classification_is_indeterminate_naming_the_allele(self, diplotype: str, unknown: list[str]) -> None:
        result = classify_phenotype(diplotype, "CYP2D6")
        assert result["phenotype"] == IND
        assert result["phenotype_abbreviation"] == "IND"
        assert result["activity_score"] is None
        assert result["unknown_alleles"] == unknown
        for name in unknown:
            assert name in result["clinical_significance"]
        assert "indeterminate" in result["clinical_significance"].lower()

    @pytest.mark.parametrize("gene", CPIC_PHENOTYPE_GENES)
    def test_every_supported_gene_calls_an_unknown_allele_indeterminate(self, gene: str) -> None:
        result = classify_phenotype("*1/*999", gene)
        assert result["phenotype"] == IND and result["activity_score"] is None
        assert result["unknown_alleles"] == ["*999"]

    def test_known_alleles_are_not_blamed_and_still_classify(self) -> None:
        result = classify_phenotype("*1/*4", "CYP2D6")
        assert result["phenotype"] == IM and result["activity_score"] == 1.0
        assert result["unknown_alleles"] == []

    def test_unknown_alleles_are_never_scored_as_the_reference(self) -> None:
        # Both used to read as Normal Metabolizer with an activity score of 2.0
        for diplotype in ("*1/*999", "*998/*999"):
            result = classify_phenotype(diplotype, "CYP2D6")
            assert result["phenotype"] != NM
            assert result["activity_score"] != 2.0

    def test_a_diplotype_built_with_an_unknown_allele_classifies_as_indeterminate(self) -> None:
        dip = determine_diplotype("*4", "*999", "CYP2D6")
        assert dip.activity_score is None
        result = classify_phenotype(dip, "CYP2D6")
        assert result["phenotype"] == IND and result["unknown_alleles"] == ["*999"]

    def test_a_scoring_table_that_lists_the_allele_characterizes_it(self) -> None:
        dip = determine_diplotype("*1", "*999", "CYP2D6", scoring_table={"*1": 1.0, "*999": 0.0})
        assert dip.activity_score == 1.0
        assert classify_phenotype(dip, "CYP2D6")["phenotype"] == IM

    def test_a_score_the_caller_supplies_is_classified_not_second_guessed(self) -> None:
        result = classify_phenotype(Diplotype("*1", "*999", "CYP2D6", activity_score=2.0), "CYP2D6")
        assert result["phenotype"] == NM and result["activity_score"] == 2.0
        assert result["unknown_alleles"] == []

    def test_uncharacterized_alleles_helper_lists_each_once_in_order(self) -> None:
        from metainformant.pharmacogenomics.alleles.diplotype import uncharacterized_alleles

        assert uncharacterized_alleles("*1", "*4", "CYP2D6") == []
        assert uncharacterized_alleles("999", "*1", "CYP2D6") == ["*999"]
        assert uncharacterized_alleles("*7777", "*7777", "CYP2D6") == ["*7777"]
        assert uncharacterized_alleles("*2", "*1", "CYP2D6", {"*1": 1.0}) == ["*2"]
        assert uncharacterized_alleles("*1", "*2", "NOTAGENE") == ["*1", "*2"]

    def test_dpyd_non_star_names_are_known_alleles(self) -> None:
        result = classify_phenotype("*1/HapB3", "DPYD")
        assert result["phenotype"] == IM and result["unknown_alleles"] == []
        assert classify_phenotype("*1/HapB4", "DPYD")["phenotype"] == IND


class TestPopulationFrequenciesCountUncharacterizedAllelesAsIndeterminate:
    def test_mass_on_an_unknown_allele_is_indeterminate_not_normal(self) -> None:
        freqs = population_phenotype_frequencies({"*1": 0.5, "*999": 0.5}, "CYP2D6")
        # *1/*1 (0.25) is NM; *1/*999 (0.5) and *999/*999 (0.25) are indeterminate. It used to be {"NM": 1.0}.
        assert freqs == {"NM": 0.25, "IND": 0.75}
        assert sum(freqs.values()) == pytest.approx(1.0)

    def test_without_unknown_alleles_there_is_no_indeterminate_entry(self) -> None:
        freqs = population_phenotype_frequencies({"*1": 0.7, "*4": 0.3}, "CYP2D6")
        assert "IND" not in freqs
        assert sum(freqs.values()) == pytest.approx(1.0)

    def test_names_that_normalize_to_the_same_allele_add_up(self) -> None:
        # "4" and "*4" are one allele: the second used to overwrite the first
        assert population_phenotype_frequencies({"*1": 0.8, "4": 0.1, "*4": 0.1}, "CYP2D6") == (
            population_phenotype_frequencies({"*1": 0.8, "*4": 0.2}, "CYP2D6")
        )

    @pytest.mark.parametrize("bad", [-0.1, math.nan])
    def test_invalid_frequency_is_rejected(self, bad: float) -> None:
        with pytest.raises(ValueError, match="non-negative"):
            population_phenotype_frequencies({"*1": 0.9, "*4": bad}, "CYP2D6")

    @pytest.mark.parametrize("freqs", [{}, {"*1": 0.0, "*4": 0.0}])
    def test_no_positive_frequency_is_an_error_not_an_empty_distribution(self, freqs: dict[str, float]) -> None:
        with pytest.raises(ValueError, match="positive frequency"):
            population_phenotype_frequencies(freqs, "CYP2D6")


class TestMetabolismFacadeAgreesWithTheAllelesPipeline:
    """Rule 2 in the string/dict front end."""

    def test_unknown_allele_is_indeterminate_with_no_score_and_is_named(self) -> None:
        result = predict_metabolizer({"diplotype": "*1/*999"}, "CYP2D6")
        assert result["phenotype"] == "indeterminate"
        assert result["activity_score"] is None
        assert result["unknown_alleles"] == ["*999"]
        assert result["alleles"] == ("*1", "*999")
        assert result["diplotype"] == "*1/*999"
        assert "poor" in result["thresholds_used"]  # the gene's bands are still reported

    def test_alleles_key_form_is_handled_alike(self) -> None:
        result = predict_metabolizer({"alleles": ["*998", "*999"]}, "CYP2D6")
        assert result["phenotype"] == "indeterminate" and result["unknown_alleles"] == ["*998", "*999"]

    def test_both_apis_give_the_same_answer_for_unknown_alleles(self) -> None:
        for gene in CPIC_PHENOTYPE_GENES:
            classified = classify_phenotype("*1/*999", gene)
            predicted = predict_metabolizer({"diplotype": "*1/*999"}, gene)
            assert classified["phenotype"].value.split()[0].lower() == predicted["phenotype"] == "indeterminate"
            assert classified["activity_score"] is None and predicted["activity_score"] is None
            assert classified["unknown_alleles"] == predicted["unknown_alleles"] == ["*999"]

    def test_known_alleles_report_no_unknown_alleles(self) -> None:
        result = predict_metabolizer({"diplotype": "*1/*4"}, "CYP2D6")
        assert result["unknown_alleles"] == [] and result["activity_score"] == 1.0

    def test_a_custom_table_characterizes_the_allele(self) -> None:
        result = predict_metabolizer({"diplotype": "*1/*999"}, "CYP2D6", activity_scores={"*1": 1.0, "*999": 0.0})
        assert result["activity_score"] == 1.0 and result["phenotype"] == "intermediate"
        assert result["unknown_alleles"] == []

    def test_compute_activity_score_stays_strict_about_unknown_alleles(self) -> None:
        # A float-valued score cannot be indeterminate: it raises instead of assuming normal function.
        table = default_allele_function_table()["CYP2D6"]
        with pytest.raises(ValueError, match="'\\*999' not found"):
            compute_activity_score("*1/*999", table)

    def test_a_malformed_diplotype_still_raises(self) -> None:
        with pytest.raises(ValueError, match="Invalid diplotype format"):
            predict_metabolizer({"diplotype": "*1"}, "CYP2D6")

    def test_unknown_gene_without_a_table_names_the_supported_genes(self) -> None:
        with pytest.raises(ValueError, match="No allele function table") as caught:
            predict_metabolizer({"diplotype": "*1/*1"}, "INVALID_GENE")
        assert "INVALID_GENE" in str(caught.value)
        assert "CYP2D6" in str(caught.value) and "CYP1A2" not in str(caught.value)


class TestIndeterminateIsNeverStandardDosing:
    """The indeterminate status must not come out of dose_adjustment as a normal metabolizer."""

    def test_a_facade_result_with_an_unknown_allele_gets_no_dose(self) -> None:
        status = predict_metabolizer({"diplotype": "*1/*999"}, "CYP2D6")["phenotype"]
        result = dose_adjustment(status, "codeine")
        assert result["adjusted_dose_fraction"] is None  # used to be 1.0, "use standard dosing"
        assert "indeterminate" in result["recommendation"].lower()
        assert result["metabolizer_status"] == "indeterminate"
        assert result["guideline_source"] == "None"

    @pytest.mark.parametrize("status", ["indeterminate", "Indeterminate", "IND", IND])
    @pytest.mark.parametrize("drug", ["codeine", "warfarin", "not-a-drug"])
    def test_every_form_of_the_status_and_any_drug(self, status: str | MetabolizerPhenotype, drug: str) -> None:
        result = dose_adjustment(status, drug)
        assert result["adjusted_dose_fraction"] is None
        assert result["evidence_level"] == "D"

    def test_a_guideline_that_defines_indeterminate_is_honoured(self) -> None:
        custom = {
            "exampledrug": {
                "indeterminate": {
                    "recommendation": "Monitor levels",
                    "adjusted_dose_fraction": 0.8,
                    "evidence_level": "C",
                    "guideline_source": "local",
                }
            }
        }
        result = dose_adjustment("Indeterminate", "ExampleDrug", dose_guidelines=custom)
        assert result["adjusted_dose_fraction"] == 0.8 and result["guideline_source"] == "local"

    def test_other_statuses_are_unchanged(self) -> None:
        assert dose_adjustment("normal", "codeine")["adjusted_dose_fraction"] == 1.0
        assert dose_adjustment("poor", "codeine")["adjusted_dose_fraction"] == 0.0
        assert dose_adjustment("normal", "unknown_drug")["adjusted_dose_fraction"] == 1.0


class TestUnscoredDiplotypes:
    """Rule 3: no score is not a score of zero."""

    def test_the_default_activity_score_is_none(self) -> None:
        assert Diplotype("*1", "*4", "CYP2D6").activity_score is None  # used to be 0.0

    def test_an_unscored_diplotype_is_not_a_poor_metabolizer(self) -> None:
        for alleles, expected, score in [(("*1", "*1"), NM, 2.0), (("*1", "*4"), IM, 1.0), (("*4", "*4"), PM, 0.0)]:
            result = classify_phenotype(Diplotype(*alleles, "CYP2D6"), "CYP2D6")
            assert result["phenotype"] == expected, alleles
            assert result["activity_score"] == score, alleles

    def test_the_diplotype_object_is_not_modified_by_classifying_it(self) -> None:
        dip = Diplotype("*1", "*1", "CYP2D6")
        classify_phenotype(dip, "CYP2D6")
        assert dip.activity_score is None

    def test_every_unscored_diplotype_classifies_as_its_string_form(self) -> None:
        checked = 0
        for gene in CPIC_PHENOTYPE_GENES:
            for first, second in itertools.combinations_with_replacement(sorted(_ACTIVITY_SCORE_TABLES[gene]), 2):
                from_object = classify_phenotype(Diplotype(first, second, gene), gene)
                from_string = classify_phenotype(f"{first}/{second}", gene)
                assert from_object["phenotype"] == from_string["phenotype"], f"{gene} {first}/{second}"
                assert from_object["activity_score"] == from_string["activity_score"], f"{gene} {first}/{second}"
                checked += 1
        assert checked > 300

    def test_function_defined_genes_use_the_allele_rules_for_an_unscored_diplotype(self) -> None:
        # *9/*17 sums to 2.0 like *1/*1 but is intermediate by allele function (a stored 0.0 read as PM)
        result = classify_phenotype(Diplotype("*9", "*17", "CYP2C19"), "CYP2C19")
        assert result["phenotype"] == IM and result["activity_score"] == 2.0

    def test_allele_names_set_after_construction_are_normalized_when_scoring(self) -> None:
        dip = Diplotype("*1", "*4", "CYP2D6")
        dip.allele2 = "4"  # not re-normalized by the dataclass
        result = classify_phenotype(dip, "CYP2D6")
        assert result["phenotype"] == IM and result["activity_score"] == 1.0

    def test_an_explicit_score_of_zero_is_a_score(self) -> None:
        result = classify_phenotype(Diplotype("*1", "*1", "CYP2D6", activity_score=0.0), "CYP2D6")
        assert result["phenotype"] == PM and result["activity_score"] == 0.0

    def test_calculate_activity_score_fills_in_the_score(self) -> None:
        dip = Diplotype("*1", "*10", "CYP2D6")
        assert calculate_activity_score(dip, "CYP2D6") == 1.25
        assert dip.activity_score == 1.25
        unknown = Diplotype("*1", "*999", "CYP2D6")
        assert calculate_activity_score(unknown, "CYP2D6") is None
        assert unknown.activity_score is None

    def test_a_diplotype_of_a_gene_cpic_does_not_define_cannot_be_classified(self) -> None:
        with pytest.raises(ValueError, match="NOTAGENE"):
            classify_phenotype(Diplotype("*1", "*1", "NOTAGENE"), "NOTAGENE")

    def test_a_diplotype_of_another_gene_is_rejected(self) -> None:
        with pytest.raises(ValueError, match="belongs to gene CYP2D6, not CYP2C19"):
            classify_phenotype(Diplotype("*1", "*4", "CYP2D6"), "CYP2C19")

    @pytest.mark.parametrize("missing", [None])
    def test_a_missing_score_cannot_be_mapped_to_a_phenotype(self, missing: None) -> None:
        # Diplotype.activity_score fed straight to the score mapping used to raise TypeError
        # (and, before None, mapped 0.0 to Poor Metabolizer)
        with pytest.raises(ValueError, match="non-negative"):
            predict_metabolizer_status(missing, "CYP2D6")  # type: ignore[arg-type]

    def test_ambiguous_diplotypes_rank_unscored_ones_by_the_score_their_alleles_imply(self) -> None:
        poor = Diplotype("*1", "*4", "CYP2D6")  # 1.0
        normal = Diplotype("*1", "*2", "CYP2D6")  # 2.0
        assert resolve_ambiguous_diplotypes([poor, normal]).diplotype_string == "*1/*2"
        assert resolve_ambiguous_diplotypes([normal, poor]).diplotype_string == "*1/*2"

    def test_ambiguous_diplotypes_rank_an_unscorable_one_after_a_scored_one_with_as_many_known_alleles(self) -> None:
        # Each has one known allele: the one with no possible score (None) ranks after the one carrying a score
        unscorable = Diplotype("*1", "*999", "CYP2D6")
        scored = Diplotype("*1", "*998", "CYP2D6", activity_score=1.5)
        assert resolve_ambiguous_diplotypes([unscorable, scored]).diplotype_string == "*1/*998"
        assert resolve_ambiguous_diplotypes([scored, unscorable]).diplotype_string == "*1/*998"

    def test_ambiguous_diplotypes_prefer_two_known_alleles_to_unscorable_ones(self) -> None:
        unknown = Diplotype("*998", "*999", "CYP2D6")
        other_unknown = Diplotype("*996", "*997", "CYP2D6")
        known = Diplotype("*1", "*4", "CYP2D6")
        assert resolve_ambiguous_diplotypes([unknown, other_unknown, known]).diplotype_string == "*1/*4"
