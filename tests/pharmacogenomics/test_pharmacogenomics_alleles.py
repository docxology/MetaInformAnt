"""Tests for pharmacogenomics alleles: star allele calling, diplotype, phenotype.

Tests star allele dataclass, allele definitions loading, variant matching,
novel allele detection, diplotype determination, activity scoring,
ambiguous diplotype resolution, phenotype classification, and population
phenotype frequency estimation.

NO MOCKING -- all tests use real implementations.
"""

from __future__ import annotations

import pytest

from metainformant.pharmacogenomics.alleles.diplotype import (
    Diplotype,
    calculate_activity_score,
    determine_diplotype,
    phased_diplotype,
    resolve_ambiguous_diplotypes,
)
from metainformant.pharmacogenomics.alleles.phenotype import (
    MetabolizerPhenotype,
    classify_phenotype,
    get_phenotype_thresholds,
    population_phenotype_frequencies,
    predict_metabolizer_status,
)
from metainformant.pharmacogenomics.alleles.star_allele import (
    StarAllele,
    call_star_alleles,
    detect_novel_alleles,
    handle_cyp2d6_cnv,
    load_allele_definitions,
    match_allele_definition,
)


class TestStarAllele:
    """Tests for the StarAllele dataclass and its methods."""

    def test_star_allele_creation_basic(self) -> None:
        allele = StarAllele(name="*1", gene="CYP2D6")
        assert allele.name == "*1"
        assert allele.gene == "CYP2D6"
        assert allele.function == "Unknown function"
        assert allele.activity_value == 1.0

    def test_star_allele_auto_prefix(self) -> None:
        allele = StarAllele(name="4", gene="cyp2d6")
        assert allele.name == "*4"
        assert allele.gene == "CYP2D6"

    def test_star_allele_is_reference(self) -> None:
        ref = StarAllele(name="*1", gene="CYP2D6")
        assert ref.is_reference is True
        variant = StarAllele(name="*4", gene="CYP2D6")
        assert variant.is_reference is False

    def test_star_allele_is_no_function(self) -> None:
        no_func = StarAllele(name="*4", gene="CYP2D6", function="No function", activity_value=0.0)
        assert no_func.is_no_function is True
        normal = StarAllele(name="*1", gene="CYP2D6", function="Normal function", activity_value=1.0)
        assert normal.is_no_function is False

    def test_star_allele_matches_variants(self) -> None:
        allele = StarAllele(
            name="*4",
            gene="CYP2D6",
            defining_variants=frozenset({"rs3892097"}),
        )
        assert allele.matches_variants({"rs3892097", "rs16947"}) is True
        assert allele.matches_variants({"rs16947"}) is False

    def test_star_allele_partial_match_score(self) -> None:
        allele = StarAllele(
            name="*3A",
            gene="TPMT",
            defining_variants=frozenset({"rs1800460", "rs1142345"}),
        )
        assert allele.partial_match_score({"rs1800460", "rs1142345"}) == 1.0
        assert allele.partial_match_score({"rs1800460"}) == 0.5
        assert allele.partial_match_score(set()) == 0.0


class TestCallStarAlleles:
    """Tests for star allele calling from observed variants."""

    def test_call_star_alleles_known_variant(self) -> None:
        result = call_star_alleles({"rs3892097"}, "CYP2D6")
        allele_names = [a.name for a in result]
        assert "*4" in allele_names

    def test_call_star_alleles_no_variants_returns_reference(self) -> None:
        result = call_star_alleles(set(), "CYP2D6")
        assert len(result) == 1
        assert result[0].name == "*1"

    def test_call_star_alleles_multiple_variants(self) -> None:
        result = call_star_alleles({"rs3892097", "rs16947"}, "CYP2D6")
        allele_names = [a.name for a in result]
        assert "*2" in allele_names
        assert "*4" in allele_names

    def test_call_star_alleles_cyp2c19(self) -> None:
        result = call_star_alleles({"rs4244285"}, "CYP2C19")
        allele_names = [a.name for a in result]
        assert "*2" in allele_names


class TestLoadAlleleDefinitions:
    """Tests for loading allele definition tables."""

    def test_load_builtin_cyp2d6(self) -> None:
        alleles = load_allele_definitions("CYP2D6")
        assert len(alleles) > 0
        names = [a.name for a in alleles]
        assert "*1" in names
        assert "*4" in names

    def test_load_builtin_cyp2c19(self) -> None:
        alleles = load_allele_definitions("CYP2C19")
        assert len(alleles) > 0
        names = [a.name for a in alleles]
        assert "*17" in names

    def test_load_builtin_case_insensitive(self) -> None:
        alleles = load_allele_definitions("cyp2d6")
        assert len(alleles) > 0

    def test_load_unknown_gene_raises(self) -> None:
        with pytest.raises(ValueError, match="not found"):
            load_allele_definitions("INVALID_GENE")


class TestMatchAlleleDefinition:
    """Tests for matching observed variants against allele definitions."""

    def test_full_match(self) -> None:
        allele_def = StarAllele(name="*4", gene="CYP2D6", defining_variants=frozenset({"rs3892097"}))
        result = match_allele_definition({"rs3892097"}, allele_def)
        assert result["is_match"] is True
        assert result["match_score"] == 1.0
        assert result["allele"] == "*4"

    def test_no_match(self) -> None:
        allele_def = StarAllele(name="*4", gene="CYP2D6", defining_variants=frozenset({"rs3892097"}))
        result = match_allele_definition({"rs16947"}, allele_def)
        assert result["is_match"] is False
        assert "rs3892097" in result["missing_variants"]

    def test_reference_allele_no_defining(self) -> None:
        allele_def = StarAllele(name="*1", gene="CYP2D6", defining_variants=frozenset())
        result = match_allele_definition({"rs16947"}, allele_def)
        assert result["is_match"] is False
        assert result["match_score"] == 0.0


class TestDetectNovelAlleles:
    """Tests for novel allele detection."""

    def test_detect_novel_with_unknown_variant(self) -> None:
        known = load_allele_definitions("CYP2D6")
        result = detect_novel_alleles({"rs3892097", "rs9999999"}, known)
        assert result["has_novel"] is True
        assert "rs9999999" in result["unmatched_variants"]

    def test_no_novel_all_explained(self) -> None:
        known = load_allele_definitions("CYP2D6")
        result = detect_novel_alleles({"rs3892097"}, known)
        assert result["has_novel"] is False

    def test_partial_matches_returned(self) -> None:
        known = load_allele_definitions("TPMT")
        # rs1800460 is part of *3A (needs both rs1800460 + rs1142345)
        # and is also *3B on its own -- it should fully match *3B
        result = detect_novel_alleles({"rs1800460", "rs9999999"}, known)
        assert result["has_novel"] is True


class TestDiplotype:
    """Tests for the Diplotype dataclass."""

    def test_diplotype_creation(self) -> None:
        dip = Diplotype(allele1="*1", allele2="*4", gene="CYP2D6")
        assert dip.gene == "CYP2D6"
        assert dip.allele1 == "*1"
        assert dip.allele2 == "*4"

    def test_diplotype_canonical_ordering(self) -> None:
        dip = Diplotype(allele1="*4", allele2="*1", gene="CYP2D6")
        assert dip.allele1 == "*1"
        assert dip.allele2 == "*4"

    def test_diplotype_string(self) -> None:
        dip = Diplotype(allele1="*1", allele2="*4", gene="CYP2D6")
        assert dip.diplotype_string == "*1/*4"

    def test_diplotype_is_homozygous(self) -> None:
        homo = Diplotype(allele1="*1", allele2="*1", gene="CYP2D6")
        assert homo.is_homozygous is True
        hetero = Diplotype(allele1="*1", allele2="*4", gene="CYP2D6")
        assert hetero.is_homozygous is False

    def test_diplotype_is_reference(self) -> None:
        ref = Diplotype(allele1="*1", allele2="*1", gene="CYP2D6")
        assert ref.is_reference is True


class TestDetermineDiplotype:
    """Tests for diplotype determination from star alleles."""

    def test_determine_diplotype_strings(self) -> None:
        dip = determine_diplotype("*1", "*4", "CYP2D6")
        assert dip.allele1 == "*1"
        assert dip.allele2 == "*4"
        assert dip.activity_score == 1.0  # *1=1.0 + *4=0.0

    def test_determine_diplotype_star_allele_objects(self) -> None:
        a1 = StarAllele(name="*1", gene="CYP2D6", activity_value=1.0)
        a2 = StarAllele(name="*2", gene="CYP2D6", activity_value=1.0)
        dip = determine_diplotype(a1, a2, "CYP2D6")
        assert dip.activity_score == 2.0

    def test_determine_diplotype_cyp2c19_rapid(self) -> None:
        dip = determine_diplotype("*1", "*17", "CYP2C19")
        assert dip.activity_score == 2.5  # *1=1.0 + *17=1.5


class TestCalculateActivityScore:
    """Tests for activity score calculation."""

    def test_calculate_activity_score_normal(self) -> None:
        dip = Diplotype(allele1="*1", allele2="*1", gene="CYP2D6", activity_score=0.0)
        score = calculate_activity_score(dip, "CYP2D6")
        assert score == 2.0
        assert dip.activity_score == 2.0

    def test_calculate_activity_score_poor(self) -> None:
        dip = Diplotype(allele1="*4", allele2="*4", gene="CYP2D6", activity_score=0.0)
        score = calculate_activity_score(dip, "CYP2D6")
        assert score == 0.0


class TestResolveAmbiguousDiplotypes:
    """Tests for resolving ambiguous diplotype calls."""

    def test_resolve_single_diplotype(self) -> None:
        dip = Diplotype(allele1="*1", allele2="*4", gene="CYP2D6", activity_score=1.0)
        result = resolve_ambiguous_diplotypes([dip])
        assert result.diplotype_string == "*1/*4"

    def test_resolve_multiple_diplotypes(self) -> None:
        dip1 = Diplotype(allele1="*1", allele2="*4", gene="CYP2D6", activity_score=1.0)
        dip2 = Diplotype(allele1="*2", allele2="*3", gene="CYP2D6", activity_score=1.0)
        result = resolve_ambiguous_diplotypes([dip1, dip2])
        assert result.confidence == "moderate"

    def test_resolve_empty_raises(self) -> None:
        with pytest.raises(ValueError, match="empty"):
            resolve_ambiguous_diplotypes([])


class TestPhasedDiplotype:
    """Tests for phased diplotype determination."""

    def test_phased_diplotype_basic(self) -> None:
        variants = {"rs3892097", "rs16947"}
        phase_data = {"rs3892097": 0, "rs16947": 1}
        dip = phased_diplotype(variants, phase_data, "CYP2D6")
        assert dip.phased is True
        assert dip.confidence == "high"
        assert dip.gene == "CYP2D6"


class TestMetabolizerPhenotype:
    """Tests for the MetabolizerPhenotype enum."""

    def test_enum_values(self) -> None:
        assert MetabolizerPhenotype.POOR.value == "Poor Metabolizer"
        assert MetabolizerPhenotype.INTERMEDIATE.value == "Intermediate Metabolizer"
        assert MetabolizerPhenotype.NORMAL.value == "Normal Metabolizer"
        assert MetabolizerPhenotype.RAPID.value == "Rapid Metabolizer"
        assert MetabolizerPhenotype.ULTRARAPID.value == "Ultrarapid Metabolizer"
        assert MetabolizerPhenotype.INDETERMINATE.value == "Indeterminate"

    def test_abbreviation(self) -> None:
        assert MetabolizerPhenotype.POOR.abbreviation == "PM"
        assert MetabolizerPhenotype.NORMAL.abbreviation == "NM"
        assert MetabolizerPhenotype.ULTRARAPID.abbreviation == "UM"


class TestClassifyPhenotype:
    """Tests for diplotype-to-phenotype classification."""

    def test_classify_poor_metabolizer(self) -> None:
        result = classify_phenotype("*4/*4", "CYP2D6")
        assert result["phenotype"] == MetabolizerPhenotype.POOR
        assert result["phenotype_abbreviation"] == "PM"
        assert result["activity_score"] == 0.0

    def test_classify_normal_metabolizer(self) -> None:
        result = classify_phenotype("*1/*1", "CYP2D6")
        assert result["phenotype"] == MetabolizerPhenotype.NORMAL
        assert result["activity_score"] == 2.0

    def test_classify_intermediate_metabolizer(self) -> None:
        # CYP2D6 IM threshold is [0.25, 1.0) -- AS=0.5 is intermediate
        result = classify_phenotype("*4/*9", "CYP2D6")
        assert result["phenotype"] == MetabolizerPhenotype.INTERMEDIATE
        assert result["activity_score"] == 0.5

    def test_classify_from_diplotype_object(self) -> None:
        dip = determine_diplotype("*4", "*9", "CYP2D6")
        result = classify_phenotype(dip, "CYP2D6")
        assert result["phenotype"] == MetabolizerPhenotype.INTERMEDIATE

    def test_classify_has_clinical_significance(self) -> None:
        result = classify_phenotype("*4/*4", "CYP2D6")
        assert "clinical_significance" in result
        assert len(result["clinical_significance"]) > 0


class TestPredictMetabolizerStatus:
    """Tests for predict_metabolizer_status (phenotype module)."""

    def test_predict_poor(self) -> None:
        result = predict_metabolizer_status(0.0, "CYP2D6")
        assert result == MetabolizerPhenotype.POOR

    def test_predict_normal(self) -> None:
        result = predict_metabolizer_status(2.0, "CYP2D6")
        assert result == MetabolizerPhenotype.NORMAL

    def test_predict_ultrarapid(self) -> None:
        result = predict_metabolizer_status(3.0, "CYP2D6")
        assert result == MetabolizerPhenotype.ULTRARAPID


class TestGetPhenotypeThresholds:
    """Tests for phenotype threshold retrieval."""

    def test_get_thresholds_cyp2d6(self) -> None:
        thresholds = get_phenotype_thresholds("CYP2D6")
        assert isinstance(thresholds, list)
        assert len(thresholds) > 0
        assert all("lower_bound" in t and "upper_bound" in t for t in thresholds)
        assert all("phenotype" in t and "abbreviation" in t for t in thresholds)

    def test_get_thresholds_unknown_gene(self) -> None:
        thresholds = get_phenotype_thresholds("UNKNOWN_GENE")
        assert isinstance(thresholds, list)
        assert len(thresholds) > 0
        assert "note" in thresholds[0]


class TestPopulationPhenotypeFrequencies:
    """Tests for population phenotype frequency estimation."""

    def test_population_frequencies_cyp2d6(self) -> None:
        allele_freqs = {"*1": 0.70, "*4": 0.20, "*10": 0.10}
        result = population_phenotype_frequencies(allele_freqs, "CYP2D6")
        assert isinstance(result, dict)
        assert len(result) > 0
        total = sum(result.values())
        assert abs(total - 1.0) < 0.01

    def test_population_frequencies_all_reference(self) -> None:
        allele_freqs = {"*1": 1.0}
        result = population_phenotype_frequencies(allele_freqs, "CYP2D6")
        assert "NM" in result
        assert result["NM"] > 0.9


class TestHandleCYP2D6CNV:
    """Tests for CYP2D6 copy number variation handling."""

    def test_homozygous_deletion(self) -> None:
        result = handle_cyp2d6_cnv(set(), 0)
        assert result["cnv_type"] == "deletion"
        assert result["diplotype_string"] == "*5/*5"
        assert result["total_activity_score"] == 0.0

    def test_normal_diploid(self) -> None:
        result = handle_cyp2d6_cnv({"rs3892097"}, 2)
        assert result["cnv_type"] == "normal"
        assert result["copy_number"] == 2

    def test_gene_duplication(self) -> None:
        result = handle_cyp2d6_cnv(set(), 3)
        assert result["cnv_type"] == "duplication"
        assert result["total_activity_score"] > 2.0
