"""Tests for pharmacogenomics annotations: CPIC, PharmGKB, drug labels.

Tests CPIC guideline loading/lookup, dosing recommendations, actionable gene
listing, allele definition parsing, PharmGKB annotation queries, evidence level
interpretation, variant annotations, drug label parsing, biomarker extraction,
and label type classification.

REAL IMPLEMENTATION -- all tests use real implementations.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from metainformant.pharmacogenomics.annotations.cpic import (
    get_dosing_recommendation,
    list_actionable_genes,
    load_cpic_guidelines,
    lookup_drug_gene,
    parse_cpic_allele_definitions,
)
from metainformant.pharmacogenomics.annotations.drug_labels import (
    _contains_pgx_terms,
    classify_label_type,
    extract_biomarker_info,
    parse_drug_label,
    search_labels_by_gene,
)
from metainformant.pharmacogenomics.annotations.pharmgkb import (
    get_evidence_level,
    get_variant_annotations,
    parse_clinical_annotations,
    query_pharmgkb_annotations,
    search_drug_pathways,
)


class TestLoadCPICGuidelines:
    """Tests for CPIC guideline loading."""

    def test_load_builtin_guidelines(self) -> None:
        guidelines = load_cpic_guidelines()
        assert isinstance(guidelines, list)
        assert len(guidelines) > 0

    def test_builtin_guidelines_structure(self) -> None:
        guidelines = load_cpic_guidelines()
        for g in guidelines:
            assert "drug" in g
            assert "gene" in g
            assert "cpic_level" in g
            assert "recommendations" in g

    def test_load_from_json_file(self, tmp_path: "pytest.TempPathFactory") -> None:
        data = [
            {
                "drug": "testdrug",
                "gene": "TESTGENE",
                "cpic_level": "A",
                "guideline_url": "",
                "recommendations": {},
            }
        ]
        filepath = tmp_path / "guidelines.json"
        filepath.write_text(json.dumps(data))
        result = load_cpic_guidelines(filepath=filepath)
        assert len(result) == 1
        assert result[0]["drug"] == "testdrug"


class TestLookupDrugGene:
    """Tests for drug-gene pair lookups."""

    def test_lookup_codeine_cyp2d6(self) -> None:
        result = lookup_drug_gene("codeine", "CYP2D6")
        assert result is not None
        assert result["drug"] == "codeine"
        assert result["gene"] == "CYP2D6"
        assert result["cpic_level"] == "A"

    def test_lookup_clopidogrel_cyp2c19(self) -> None:
        result = lookup_drug_gene("clopidogrel", "CYP2C19")
        assert result is not None
        assert result["gene"] == "CYP2C19"

    def test_lookup_case_insensitive(self) -> None:
        result = lookup_drug_gene("CODEINE", "cyp2d6")
        assert result is not None

    def test_lookup_nonexistent_pair(self) -> None:
        result = lookup_drug_gene("aspirin", "CYP2D6")
        assert result is None

    def test_lookup_warfarin_cyp2c9(self) -> None:
        result = lookup_drug_gene("warfarin", "CYP2C9")
        assert result is not None
        assert result["cpic_level"] == "A"


class TestGetDosingRecommendation:
    """Tests for dosing recommendation retrieval."""

    def test_dosing_codeine_poor_metabolizer(self) -> None:
        result = get_dosing_recommendation("codeine", "Poor Metabolizer")
        assert result is not None
        assert "codeine" in result["drug"]
        assert "recommendation" in result
        assert len(result["recommendation"]) > 0

    def test_dosing_with_abbreviation(self) -> None:
        result = get_dosing_recommendation("codeine", "PM")
        assert result is not None

    def test_dosing_normal_metabolizer(self) -> None:
        result = get_dosing_recommendation("codeine", "Normal Metabolizer")
        assert result is not None
        assert "standard" in result["recommendation"].lower()

    def test_dosing_unknown_drug(self) -> None:
        result = get_dosing_recommendation("nonexistent_drug", "PM")
        assert result is None


class TestListActionableGenes:
    """Tests for listing CPIC actionable gene-drug pairs."""

    def test_list_level_a(self) -> None:
        result = list_actionable_genes(min_level="A")
        assert isinstance(result, list)
        assert len(result) > 0
        for entry in result:
            assert "gene" in entry
            assert "drug" in entry
            assert "cpic_level" in entry

    def test_list_level_b(self) -> None:
        result_a = list_actionable_genes(min_level="A")
        result_b = list_actionable_genes(min_level="B")
        assert len(result_b) >= len(result_a)

    def test_contains_known_pairs(self) -> None:
        result = list_actionable_genes()
        gene_drug_pairs = {(e["gene"], e["drug"]) for e in result}
        assert ("CYP2D6", "codeine") in gene_drug_pairs
        assert ("CYP2C19", "clopidogrel") in gene_drug_pairs


class TestParseCPICAlleleDefinitions:
    """Tests for CPIC allele definition parsing from file."""

    def test_parse_json_file(self, tmp_path: "pytest.TempPathFactory") -> None:
        data = {
            "CYP2D6": {
                "*1": {
                    "defining_variants": [],
                    "function": "Normal function",
                    "activity_value": 1.0,
                },
                "*4": {
                    "defining_variants": ["rs3892097"],
                    "function": "No function",
                    "activity_value": 0.0,
                },
            }
        }
        filepath = tmp_path / "allele_defs.json"
        filepath.write_text(json.dumps(data))
        result = parse_cpic_allele_definitions(filepath)
        assert "CYP2D6" in result
        assert len(result["CYP2D6"]) == 2

    def test_parse_tsv_file(self, tmp_path: "pytest.TempPathFactory") -> None:
        lines = [
            "gene\tallele\trsid\tposition\tref\talt\tfunction\tactivity_value",
            "CYP2D6\t*4\trs3892097\t\t\t\tNo function\t0.0",
            "CYP2C19\t*2\trs4244285\t\t\t\tNo function\t0.0",
        ]
        filepath = tmp_path / "allele_defs.tsv"
        filepath.write_text("\n".join(lines))
        result = parse_cpic_allele_definitions(filepath)
        assert "CYP2D6" in result
        assert "CYP2C19" in result

    def test_parse_missing_file_raises(self) -> None:
        with pytest.raises(FileNotFoundError):
            parse_cpic_allele_definitions("/nonexistent/path.tsv")


class TestQueryPharmGKBAnnotations:
    """Tests for PharmGKB annotation queries."""

    def test_query_by_gene(self) -> None:
        result = query_pharmgkb_annotations(gene="CYP2D6")
        assert isinstance(result, list)
        assert len(result) > 0
        assert all(a["gene"] == "CYP2D6" for a in result)

    def test_query_by_drug(self) -> None:
        result = query_pharmgkb_annotations(drug="codeine")
        assert len(result) > 0
        assert all(a["drug"] == "codeine" for a in result)

    def test_query_by_gene_and_drug(self) -> None:
        result = query_pharmgkb_annotations(gene="CYP2C19", drug="clopidogrel")
        assert len(result) > 0

    def test_query_no_filters_raises(self) -> None:
        with pytest.raises(ValueError, match="At least one"):
            query_pharmgkb_annotations()


class TestGetEvidenceLevel:
    """Tests for evidence level interpretation."""

    def test_level_1a(self) -> None:
        result = get_evidence_level({"evidence_level": "1A"})
        assert result["level"] == "1A"
        assert result["strength"] == "Strong"
        assert result["actionable"] is True

    def test_level_3(self) -> None:
        result = get_evidence_level({"evidence_level": "3"})
        assert result["strength"] == "Weak"
        assert result["actionable"] is False

    def test_unknown_level(self) -> None:
        result = get_evidence_level({"evidence_level": "Z"})
        assert result["strength"] == "Unknown"
        assert result["actionable"] is False


class TestGetVariantAnnotations:
    """Tests for variant-specific annotations."""

    def test_known_variant_rs4244285(self) -> None:
        result = get_variant_annotations("rs4244285")
        assert result is not None
        assert result["gene"] == "CYP2C19"
        assert result["allele"] == "*2"
        assert "allele_frequency" in result

    def test_known_variant_rs3892097(self) -> None:
        result = get_variant_annotations("rs3892097")
        assert result is not None
        assert result["gene"] == "CYP2D6"
        assert result["allele"] == "*4"

    def test_unknown_variant(self) -> None:
        result = get_variant_annotations("rs000000000")
        assert result is None


class TestSearchDrugPathways:
    """Tests for drug pathway queries."""

    def test_search_clopidogrel(self) -> None:
        result = search_drug_pathways("clopidogrel")
        assert result is not None
        assert "CYP2C19" in result["genes"]
        assert result["key_metabolizing_enzyme"] == "CYP2C19"

    def test_search_warfarin(self) -> None:
        result = search_drug_pathways("warfarin")
        assert result is not None
        assert "CYP2C9" in result["genes"]

    def test_search_nonexistent(self) -> None:
        result = search_drug_pathways("nonexistent_drug")
        assert result is None


class TestParseDrugLabel:
    """Tests for FDA drug label parsing."""

    def test_parse_json_label(self, tmp_path: "pytest.TempPathFactory") -> None:
        label_data = {
            "drug": "testdrug",
            "brand_name": "TestBrand",
            "gene_biomarker": "CYP2D6",
            "label_type": "actionable",
            "warnings_and_precautions": "CYP2D6 poor metabolizer status may affect drug efficacy.",
        }
        filepath = tmp_path / "label.json"
        filepath.write_text(json.dumps(label_data))
        result = parse_drug_label(filepath)
        assert result["drug"] == "testdrug"
        assert "biomarker_info" in result
        assert "label_type" in result

    def test_parse_missing_file_raises(self) -> None:
        with pytest.raises(FileNotFoundError):
            parse_drug_label("/nonexistent/path.json")


class TestExtractBiomarkerInfo:
    """Tests for biomarker info extraction."""

    def test_extract_required_testing(self) -> None:
        label = {
            "gene_biomarker": "HLA-B",
            "label_type": "required",
            "sections": {"warnings_and_precautions": "Screen for HLA-B*5701 allele before starting therapy."},
        }
        result = extract_biomarker_info(label)
        assert result["biomarker_gene"] == "HLA-B"
        assert result["testing_required"] is True

    def test_extract_informational(self) -> None:
        label = {
            "gene_biomarker": "CYP2D6",
            "label_type": "informational",
            "sections": {},
        }
        result = extract_biomarker_info(label)
        assert result["testing_required"] is False


class TestClassifyLabelType:
    """Tests for drug label type classification."""

    def test_classify_required(self) -> None:
        label = {
            "label_type": "required",
            "sections": {},
        }
        result = classify_label_type(label)
        assert result["label_type"] == "required"

    def test_classify_actionable(self) -> None:
        label = {
            "label_type": "actionable",
            "sections": {},
        }
        result = classify_label_type(label)
        assert result["label_type"] == "actionable"

    def test_classify_informational(self) -> None:
        label = {
            "label_type": "informational",
            "sections": {},
        }
        result = classify_label_type(label)
        assert result["label_type"] == "informational"


class TestSearchLabelsByGene:
    """Tests for searching drug labels by gene."""

    def test_search_cyp2d6(self) -> None:
        result = search_labels_by_gene("CYP2D6")
        assert isinstance(result, list)
        assert len(result) > 0
        drugs = [r["drug"] for r in result]
        assert "codeine" in drugs

    def test_search_hla_b(self) -> None:
        result = search_labels_by_gene("HLA-B")
        assert len(result) > 0
        drugs = [r["drug"] for r in result]
        assert "abacavir" in drugs

    def test_search_nonexistent_gene(self) -> None:
        result = search_labels_by_gene("NONEXISTENT_GENE")
        assert len(result) == 0

    def test_search_case_insensitive(self) -> None:
        result = search_labels_by_gene("cyp2d6")
        assert len(result) > 0


class TestParseClinicalAnnotations:
    """Tests for PharmGKB clinical annotation parsing from file."""

    def test_parse_tsv_file(self, tmp_path: "pytest.TempPathFactory") -> None:
        lines = [
            "Clinical Annotation ID\tGene\tLevel of Evidence\tDrug(s)\tPhenotype Category\tAnnotation Text",
            "PA001\tCYP2D6\t1A\tcodeine\tEfficacy\tTest annotation text",
            "PA002\tCYP2C19\t1A\tclopidogrel;voriconazole\tEfficacy\tAnother annotation",
        ]
        filepath = tmp_path / "annotations.tsv"
        filepath.write_text("\n".join(lines))
        result = parse_clinical_annotations(filepath)
        assert isinstance(result, list)
        assert len(result) >= 2
        # Second row has two drugs semicolon-separated
        assert len(result) == 3

    def test_parse_missing_file_raises(self) -> None:
        with pytest.raises(FileNotFoundError):
            parse_clinical_annotations("/nonexistent/path.tsv")


class TestDosingRecommendationPhenotypeForms:
    """get_dosing_recommendation takes the forms MetabolizerPhenotype.from_string takes."""

    @pytest.mark.parametrize(
        ("label", "expected"),
        [
            ("poor", "Poor Metabolizer"),
            ("Poor", "Poor Metabolizer"),
            ("PM", "Poor Metabolizer"),
            ("Poor Metabolizer", "Poor Metabolizer"),
            ("ultrarapid", "Ultrarapid Metabolizer"),
            ("ultra-rapid", "Ultrarapid Metabolizer"),
            ("UM", "Ultrarapid Metabolizer"),
            ("Extensive Metabolizer", "Normal Metabolizer"),
            ("EM", "Normal Metabolizer"),
        ],
    )
    def test_codeine_recommendation_for_any_form(self, label: str, expected: str) -> None:
        # The bare tokens ("poor", "ultrarapid") are what metabolism.predict_metabolizer_status returns;
        # they used to match nothing
        result = get_dosing_recommendation("codeine", label)
        assert result is not None
        assert result["phenotype"] == expected

    @pytest.mark.parametrize("label", ["Indeterminate", "superrapid", "", "Likely Poor Metabolizer"])
    def test_a_label_that_names_no_guideline_phenotype_finds_nothing(self, label: str) -> None:
        assert get_dosing_recommendation("codeine", label) is None


class TestListActionableGenesLevels:
    @pytest.mark.parametrize("level", ["Z", "", "A/B", "1", "E"])
    def test_an_unrecognised_minimum_level_raises_instead_of_meaning_a(self, level: str) -> None:
        with pytest.raises(ValueError, match="min_level"):
            list_actionable_genes(min_level=level)

    def test_levels_are_case_insensitive_and_nested(self) -> None:
        assert list_actionable_genes(min_level="b") == list_actionable_genes(min_level="B")
        assert list_actionable_genes(min_level=" a ") == list_actionable_genes(min_level="A")
        by_level = [len(list_actionable_genes(min_level=level)) for level in "ABCD"]
        assert by_level == sorted(by_level)

    def test_levels_below_a_admit_the_a_b_and_lower_entries(self) -> None:
        guidelines = [
            {"gene": "G1", "drug": "d1", "cpic_level": "A"},
            {"gene": "G2", "drug": "d2", "cpic_level": "A/B"},
            {"gene": "G3", "drug": "d3", "cpic_level": "B"},
            {"gene": "G4", "drug": "d4", "cpic_level": "C"},
            {"gene": "G5", "drug": "d5", "cpic_level": "D"},
        ]
        genes = {
            level: {e["gene"] for e in list_actionable_genes(guidelines=guidelines, min_level=level)}
            for level in "ABCD"
        }
        assert genes == {
            "A": {"G1"},
            "B": {"G1", "G2", "G3"},
            "C": {"G1", "G2", "G3", "G4"},
            "D": {"G1", "G2", "G3", "G4", "G5"},
        }


class TestParseCpicAlleleDefinitionsAreStrict:
    """An allele definition without an activity value is not normal function (1.0)."""

    HEADER = "gene\tallele\trsid\tposition\tref\talt\tfunction\tactivity_value"

    def _tsv(self, tmp_path: Path, *rows: str, header: str | None = None) -> Path:
        path = tmp_path / "defs.tsv"
        path.write_text("\n".join([header or self.HEADER, *rows]) + "\n", encoding="utf-8")
        return path

    def test_a_row_without_an_activity_value_raises_naming_the_allele_and_file(self, tmp_path: Path) -> None:
        path = self._tsv(tmp_path, "CYP2D6\t*4\trs3892097\t\t\t\tNo function\t")
        with pytest.raises(ValueError, match=r"CYP2D6 \*4") as caught:
            parse_cpic_allele_definitions(path)
        assert "activity_value" in str(caught.value) and "defs.tsv" in str(caught.value)

    def test_a_file_without_the_activity_column_raises(self, tmp_path: Path) -> None:
        path = self._tsv(tmp_path, "CYP2D6\t*4\trs3892097\tNo function", header="gene\tallele\trsid\tfunction")
        with pytest.raises(ValueError, match="no activity_value"):
            parse_cpic_allele_definitions(path)

    @pytest.mark.parametrize("value", ["fast", "-0.5", "nan", "inf"])
    def test_a_value_that_is_not_a_finite_non_negative_number_raises(self, tmp_path: Path, value: str) -> None:
        path = self._tsv(tmp_path, f"CYP2D6\t*4\trs3892097\t\t\t\tNo function\t{value}")
        with pytest.raises(ValueError, match="activity_value"):
            parse_cpic_allele_definitions(path)

    @pytest.mark.parametrize(
        "row", ["\t*4\trs3892097\t\t\t\tNo function\t0.0", "CYP2D6\t\trs3892097\t\t\t\tNo function\t0.0"]
    )
    def test_a_row_without_a_gene_or_allele_raises_instead_of_filing_it_under_a_blank_gene(
        self, tmp_path: Path, row: str
    ) -> None:
        with pytest.raises(ValueError, match="without a gene or an allele name"):
            parse_cpic_allele_definitions(self._tsv(tmp_path, row))

    def test_complete_rows_parse_with_their_values(self, tmp_path: Path) -> None:
        path = self._tsv(
            tmp_path,
            "CYP2D6\t*4\trs3892097\t\t\t\tNo function\t0.0",
            "CYP2D6\t*10\trs1065852\t\t\t\tDecreased function\t0.25",
            "cyp2c19\t*17\trs12248560\t\t\t\tIncreased function\t1.5",
        )
        parsed = parse_cpic_allele_definitions(path)
        assert {a["allele"]: a["activity_value"] for a in parsed["CYP2D6"]} == {"*4": 0.0, "*10": 0.25}
        assert parsed["CYP2C19"][0]["activity_value"] == 1.5
        assert parsed["CYP2C19"][0]["defining_variants"] == ["rs12248560"]


class TestParseDrugLabelClassifiesFromText:
    """A label file with no ``label_type`` used to be reported "informational" whatever it said."""

    def _label(self, tmp_path: Path, **sections: str) -> Path:
        path = tmp_path / "label.json"
        path.write_text(json.dumps({"drug": "exampledrug", "gene_biomarker": "CYP2C19", **sections}), encoding="utf-8")
        return path

    def test_mandatory_testing_language_makes_the_label_required(self, tmp_path: Path) -> None:
        path = self._label(
            tmp_path,
            warnings_and_precautions="Genotyping is required: CYP2C19 poor metabolizers must be tested.",
        )
        result = parse_drug_label(path)
        assert result["label_type"] == "required"
        assert result["biomarker_info"]["testing_required"] is True

    def test_recommended_testing_language_makes_the_label_recommended(self, tmp_path: Path) -> None:
        path = self._label(
            tmp_path, indications_and_usage="CYP2C19 genotyping should be performed before the first dose."
        )
        result = parse_drug_label(path)
        assert result["label_type"] == "recommended"
        assert result["biomarker_info"]["testing_recommended"] is True

    def test_genotype_based_dosing_language_makes_the_label_actionable(self, tmp_path: Path) -> None:
        path = self._label(
            tmp_path, dosage_and_administration="CYP2C19 poor metabolizer status: consider alternative therapy."
        )
        assert parse_drug_label(path)["label_type"] == "actionable"

    def test_a_label_with_no_testing_language_is_informational(self, tmp_path: Path) -> None:
        path = self._label(tmp_path, clinical_pharmacology="CYP2C19 contributes to the metabolism of the drug.")
        assert parse_drug_label(path)["label_type"] == "informational"

    @pytest.mark.parametrize("label_type", ["required", "recommended", "actionable", "informational"])
    def test_a_label_type_in_the_file_wins_over_the_text(self, tmp_path: Path, label_type: str) -> None:
        path = self._label(
            tmp_path,
            label_type=label_type,
            warnings_and_precautions="Genotyping is required: CYP2C19 poor metabolizers must be tested.",
        )
        assert parse_drug_label(path)["label_type"] == label_type

    def test_an_unrecognised_label_type_in_the_file_is_classified_from_the_text(self, tmp_path: Path) -> None:
        path = self._label(
            tmp_path,
            label_type="mandatory",
            warnings_and_precautions="Genotyping is required: CYP2C19 poor metabolizers must be tested.",
        )
        assert parse_drug_label(path)["label_type"] == "required"

    def test_sections_are_kept_only_when_they_are_about_pharmacogenomics(self, tmp_path: Path) -> None:
        # "General" contains "gene": it used to make this renal-dosing paragraph a PGx section
        path = self._label(
            tmp_path,
            dosage_and_administration="General guidance: use dose reduction in renal impairment.",
            warnings_and_precautions="CYP2C19 poor metabolizers have higher exposure.",
        )
        result = parse_drug_label(path)
        assert set(result["sections"]) == {"warnings_and_precautions"}

    def test_a_non_pgx_dose_reduction_does_not_make_a_label_actionable(self, tmp_path: Path) -> None:
        path = self._label(
            tmp_path,
            dosage_and_administration="General guidance: use dose reduction in renal impairment.",
        )
        result = parse_drug_label(path)
        assert result["sections"] == {}
        assert result["label_type"] == "informational"


class TestPgxTermMatching:
    @pytest.mark.parametrize(
        "text",
        [
            "CYP2D6 poor metabolizers",
            "Carriers of the HLA-B*15:02 allele",
            "genotype-guided dosing",
            "Variants of DPYD",
            "ALLELES differ",
            "the gene encodes an enzyme",
            "biomarker testing",
            "Pharmacogenomic information",
            "Genetic testing is recommended",
            "Genomic variation in the pathway",
            "CYP enzymes",
            "ultrarapid metabolizers",
        ],
    )
    def test_pharmacogenomic_language_is_recognised(self, text: str) -> None:
        assert _contains_pgx_terms(text)

    @pytest.mark.parametrize(
        "text",
        [
            "General precautions apply",
            "Generally well tolerated",
            "Degenerative disease warnings",
            "Cypress oil is not a drug",
            "Take with water",
            "",
        ],
    )
    def test_words_that_only_contain_a_term_are_not(self, text: str) -> None:
        assert not _contains_pgx_terms(text)


class TestBuiltInGuidelineTableStatesItsScope:
    def test_docstring_gives_the_real_number_of_guidelines_and_names_their_drugs(self) -> None:
        # It said the table "covers major CPIC Level A/A-B guidelines"; it holds nine, all Level A
        guidelines = load_cpic_guidelines()
        doc = load_cpic_guidelines.__doc__ or ""
        claimed = re.search(r"holds (\d+) Level A", doc)
        assert claimed is not None and int(claimed.group(1)) == len(guidelines)
        assert {g["cpic_level"] for g in guidelines} == {"A"}
        for guideline in guidelines:
            assert guideline["drug"] in doc
        assert "does not mean CPIC has no guideline" in " ".join(doc.split())
