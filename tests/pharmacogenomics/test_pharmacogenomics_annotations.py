"""Tests for pharmacogenomics annotations: CPIC, PharmGKB, drug labels.

Tests CPIC guideline loading/lookup, dosing recommendations, actionable gene
listing, allele definition parsing, PharmGKB annotation queries, evidence level
interpretation, variant annotations, drug label parsing, biomarker extraction,
and label type classification.

NO MOCKING -- all tests use real implementations.
"""

from __future__ import annotations

import json

import pytest

from metainformant.pharmacogenomics.annotations.cpic import (
    get_dosing_recommendation,
    list_actionable_genes,
    load_cpic_guidelines,
    lookup_drug_gene,
    parse_cpic_allele_definitions,
)
from metainformant.pharmacogenomics.annotations.drug_labels import (
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
