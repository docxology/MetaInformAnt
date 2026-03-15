"""Tests for pharmacogenomics clinical: ACMG classification, drug interactions, reporting.

Tests ACMG variant classification enums and logic, evidence aggregation,
ClinVar queries, gnomAD frequency checks, drug-gene interaction analysis,
interaction severity calculation, contraindication checks, polypharmacy
analysis, clinical report generation, recommendation formatting, summary
table generation, and report export.

NO MOCKING -- all tests use real implementations.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from metainformant.pharmacogenomics.clinical.drug_interaction import (
    DrugRecommendation,
    InteractionSeverity,
    analyze_drug_gene_interactions,
    calculate_interaction_severity,
    check_contraindications,
    polypharmacy_analysis,
    suggest_alternatives,
)
from metainformant.pharmacogenomics.clinical.pathogenicity import (
    ACMGClassification,
    ACMGCriteria,
    aggregate_evidence,
    apply_acmg_criteria,
    check_gnomad_frequency,
    classify_variant_acmg,
    query_clinvar,
)
from metainformant.pharmacogenomics.clinical.reporting import (
    add_disclaimer,
    export_report,
    format_recommendation,
    generate_clinical_report,
    generate_summary_table,
)


class TestACMGClassification:
    """Tests for ACMGClassification enum values."""

    def test_enum_values(self) -> None:
        assert ACMGClassification.PATHOGENIC.value == "Pathogenic"
        assert ACMGClassification.LIKELY_PATHOGENIC.value == "Likely Pathogenic"
        assert ACMGClassification.VUS.value == "Uncertain Significance"
        assert ACMGClassification.LIKELY_BENIGN.value == "Likely Benign"
        assert ACMGClassification.BENIGN.value == "Benign"

    def test_enum_is_string(self) -> None:
        assert isinstance(ACMGClassification.PATHOGENIC, str)
        assert ACMGClassification.PATHOGENIC == "Pathogenic"


class TestACMGCriteria:
    """Tests for ACMGCriteria enum."""

    def test_pathogenic_criteria(self) -> None:
        assert ACMGCriteria.PVS1.value == "PVS1"
        assert ACMGCriteria.PS1.value == "PS1"
        assert ACMGCriteria.PM1.value == "PM1"
        assert ACMGCriteria.PP1.value == "PP1"

    def test_benign_criteria(self) -> None:
        assert ACMGCriteria.BA1.value == "BA1"
        assert ACMGCriteria.BS1.value == "BS1"
        assert ACMGCriteria.BP1.value == "BP1"


class TestClassifyVariantACMG:
    """Tests for ACMG variant classification."""

    def test_classify_frameshift_variant(self) -> None:
        variant = {
            "rsid": "rs_test",
            "gene": "CYP2D6",
            "consequence": "frameshift",
            "allele_frequency": 0.00001,
        }
        result = classify_variant_acmg(variant)
        assert "classification" in result
        assert isinstance(result["classification"], ACMGClassification)
        assert "criteria_met" in result
        assert "PVS1" in result["criteria_met"]

    def test_classify_with_pre_evaluated_evidence(self) -> None:
        variant = {"rsid": "rs_test", "gene": "CYP2D6"}
        evidence = {"PVS1": True, "PS1": True, "PM2": False, "BA1": False}
        result = classify_variant_acmg(variant, evidence=evidence)
        assert result["classification"] == ACMGClassification.PATHOGENIC

    def test_classify_benign_high_frequency(self) -> None:
        variant = {
            "rsid": "rs_common",
            "gene": "CYP2D6",
            "consequence": "missense",
            "allele_frequency": {"Global": 0.10},
        }
        result = classify_variant_acmg(variant)
        assert result["classification"] == ACMGClassification.BENIGN
        assert "BA1" in result["criteria_met"]


class TestApplyACMGCriteria:
    """Tests for evaluating individual ACMG criteria."""

    def test_pvs1_frameshift(self) -> None:
        variant = {"consequence": "frameshift", "allele_frequency": 0.00001}
        criteria = apply_acmg_criteria(variant)
        assert criteria["PVS1"] is True

    def test_ba1_common_variant(self) -> None:
        variant = {"consequence": "missense", "allele_frequency": {"European": 0.10}}
        criteria = apply_acmg_criteria(variant)
        assert criteria["BA1"] is True

    def test_pm2_rare_variant(self) -> None:
        variant = {"consequence": "missense", "allele_frequency": 0.00001}
        criteria = apply_acmg_criteria(variant)
        assert criteria["PM2"] is True

    def test_pp3_computational_evidence(self) -> None:
        variant = {
            "consequence": "missense",
            "allele_frequency": 0.0001,
            "computational_predictions": {
                "CADD": "damaging",
                "REVEL": "damaging",
                "SIFT": "damaging",
            },
        }
        criteria = apply_acmg_criteria(variant)
        assert criteria["PP3"] is True


class TestAggregateEvidence:
    """Tests for combining ACMG criteria into final classification."""

    def test_pathogenic_pvs1_plus_ps1(self) -> None:
        criteria = {"PVS1": True, "PS1": True, "PM2": False, "BA1": False}
        result = aggregate_evidence(criteria)
        assert result == ACMGClassification.PATHOGENIC

    def test_benign_ba1(self) -> None:
        criteria = {"BA1": True, "PVS1": False}
        result = aggregate_evidence(criteria)
        assert result == ACMGClassification.BENIGN

    def test_likely_benign_bs_plus_bp(self) -> None:
        criteria = {"BS1": True, "BP4": True, "PVS1": False}
        result = aggregate_evidence(criteria)
        assert result == ACMGClassification.LIKELY_BENIGN

    def test_vus_default(self) -> None:
        criteria = {"PM2": True, "PVS1": False, "BA1": False}
        result = aggregate_evidence(criteria)
        assert result == ACMGClassification.VUS


class TestQueryClinvar:
    """Tests for ClinVar query."""

    def test_known_variant(self) -> None:
        result = query_clinvar("rs3892097")
        assert result is not None
        assert result["gene"] == "CYP2D6"
        assert "classification" in result

    def test_unknown_variant(self) -> None:
        result = query_clinvar("rs000000000")
        assert result is None


class TestCheckGnomadFrequency:
    """Tests for gnomAD frequency checks."""

    def test_ba1_high_frequency(self) -> None:
        variant = {"allele_frequency": {"European": 0.10, "African": 0.08}}
        result = check_gnomad_frequency(variant)
        assert result["ba1_triggered"] is True
        assert result["max_frequency"] == 0.10

    def test_no_trigger_low_frequency(self) -> None:
        variant = {"allele_frequency": {"European": 0.001}}
        result = check_gnomad_frequency(variant)
        assert result["ba1_triggered"] is False
        assert result["bs1_triggered"] is False

    def test_specific_population(self) -> None:
        variant = {"allele_frequency": {"European": 0.001, "African": 0.06}}
        result = check_gnomad_frequency(variant, population="European")
        assert result["ba1_triggered"] is False


class TestInteractionSeverity:
    """Tests for InteractionSeverity enum."""

    def test_enum_values(self) -> None:
        assert InteractionSeverity.MAJOR.value == "Major"
        assert InteractionSeverity.MODERATE.value == "Moderate"
        assert InteractionSeverity.MINOR.value == "Minor"
        assert InteractionSeverity.NONE.value == "None"

    def test_requires_action(self) -> None:
        assert InteractionSeverity.MAJOR.requires_action is True
        assert InteractionSeverity.MODERATE.requires_action is True
        assert InteractionSeverity.MINOR.requires_action is False
        assert InteractionSeverity.NONE.requires_action is False


class TestDrugRecommendation:
    """Tests for DrugRecommendation dataclass."""

    def test_creation(self) -> None:
        rec = DrugRecommendation(
            drug="codeine",
            gene="CYP2D6",
            phenotype="Poor Metabolizer",
            recommendation="Avoid codeine use.",
            evidence_level="A",
        )
        assert rec.drug == "codeine"
        assert rec.gene == "CYP2D6"
        assert rec.severity == InteractionSeverity.NONE


class TestAnalyzeDrugGeneInteractions:
    """Tests for drug-gene interaction analysis."""

    def test_codeine_poor_metabolizer(self) -> None:
        drugs = ["codeine"]
        genotypes = {
            "CYP2D6": {"phenotype": "Poor Metabolizer", "diplotype": "*4/*4"},
        }
        result = analyze_drug_gene_interactions(drugs, genotypes)
        assert isinstance(result, list)
        assert len(result) > 0
        assert result[0].drug == "codeine"
        assert result[0].severity == InteractionSeverity.MAJOR

    def test_clopidogrel_normal_metabolizer(self) -> None:
        drugs = ["clopidogrel"]
        genotypes = {
            "CYP2C19": {"phenotype": "Normal Metabolizer"},
        }
        result = analyze_drug_gene_interactions(drugs, genotypes)
        assert isinstance(result, list)
        # Normal metabolizer has no severity interaction from contra DB
        for rec in result:
            assert rec.severity != InteractionSeverity.MAJOR

    def test_multiple_drugs(self) -> None:
        drugs = ["codeine", "clopidogrel"]
        genotypes = {
            "CYP2D6": {"phenotype": "Poor Metabolizer"},
            "CYP2C19": {"phenotype": "Poor Metabolizer"},
        }
        result = analyze_drug_gene_interactions(drugs, genotypes)
        assert len(result) >= 2


class TestCalculateInteractionSeverity:
    """Tests for overall interaction severity calculation."""

    def test_no_interactions(self) -> None:
        result = calculate_interaction_severity([])
        assert result["overall_severity"] == "None"
        assert result["risk_level"] == "Low"

    def test_major_interaction(self) -> None:
        rec = DrugRecommendation(
            drug="codeine",
            gene="CYP2D6",
            phenotype="PM",
            recommendation="Avoid",
            severity=InteractionSeverity.MAJOR,
        )
        result = calculate_interaction_severity([rec])
        assert result["overall_severity"] == "Major"
        assert result["risk_level"] == "High"
        assert result["major_count"] == 1

    def test_multiple_major_very_high_risk(self) -> None:
        recs = [
            DrugRecommendation(
                drug="codeine",
                gene="CYP2D6",
                phenotype="PM",
                recommendation="Avoid",
                severity=InteractionSeverity.MAJOR,
            ),
            DrugRecommendation(
                drug="tramadol",
                gene="CYP2D6",
                phenotype="PM",
                recommendation="Avoid",
                severity=InteractionSeverity.MAJOR,
            ),
        ]
        result = calculate_interaction_severity(recs)
        assert result["risk_level"] == "Very High"


class TestCheckContraindications:
    """Tests for contraindication checking."""

    def test_codeine_pm_contraindicated(self) -> None:
        result = check_contraindications("codeine", "Poor Metabolizer")
        assert result["contraindicated"] is True
        assert len(result["gene_interactions"]) > 0

    def test_codeine_nm_not_contraindicated(self) -> None:
        result = check_contraindications("codeine", "Normal Metabolizer")
        assert result["contraindicated"] is False

    def test_unknown_drug(self) -> None:
        result = check_contraindications("aspirin", "Poor Metabolizer")
        assert result["contraindicated"] is False


class TestPolypharmacyAnalysis:
    """Tests for polypharmacy analysis."""

    def test_basic_polypharmacy(self) -> None:
        drugs = ["codeine", "clopidogrel"]
        genotypes = {
            "CYP2D6": {"phenotype": "Poor Metabolizer"},
            "CYP2C19": {"phenotype": "Poor Metabolizer"},
        }
        result = polypharmacy_analysis(drugs, genotypes)
        assert result["total_drugs"] == 2
        assert isinstance(result["interactions"], list)
        assert "severity_summary" in result
        assert "recommendations" in result


class TestSuggestAlternatives:
    """Tests for alternative drug suggestions."""

    def test_codeine_alternatives(self) -> None:
        result = suggest_alternatives("codeine", "Poor Metabolizer")
        assert result["therapeutic_class"] == "Opioid analgesic"
        assert len(result["alternatives"]) > 0

    def test_unknown_drug_alternatives(self) -> None:
        result = suggest_alternatives("unknown_drug", "PM")
        assert len(result["alternatives"]) == 0


class TestGenerateClinicalReport:
    """Tests for clinical report generation."""

    def test_basic_report(self) -> None:
        patient = {"patient_id": "TEST001", "sex": "Male", "ethnicity": "European"}
        genotypes = {
            "CYP2D6": {"diplotype": "*1/*4"},
            "CYP2C19": {"diplotype": "*1/*1"},
        }
        report = generate_clinical_report(patient, genotypes)
        assert "header" in report
        assert "patient_info" in report
        assert "genotype_results" in report
        assert "disclaimer" in report
        assert report["patient_info"]["patient_id"] == "TEST001"

    def test_report_with_drugs(self) -> None:
        patient = {"patient_id": "TEST002"}
        genotypes = {
            "CYP2D6": {"diplotype": "*4/*4", "phenotype": "Poor Metabolizer"},
        }
        drugs = ["codeine"]
        report = generate_clinical_report(patient, genotypes, drugs=drugs)
        assert len(report["drug_recommendations"]) > 0
        assert len(report["genotype_results"]) > 0

    def test_report_header_metadata(self) -> None:
        report = generate_clinical_report({}, {"CYP2D6": {"diplotype": "*1/*1"}})
        assert report["header"]["report_type"] == "Pharmacogenomic Clinical Report"
        assert "generated_at" in report["header"]


class TestFormatRecommendation:
    """Tests for recommendation formatting."""

    def test_format_major_recommendation(self) -> None:
        result = format_recommendation(
            drug="codeine",
            gene="CYP2D6",
            phenotype="Poor Metabolizer",
            guideline={
                "recommendation": "Avoid codeine use.",
                "evidence_level": "A",
                "severity": "Major",
                "source": "CPIC",
            },
        )
        assert result["drug"] == "codeine"
        assert result["urgency"] == "ACTION REQUIRED"
        assert result["severity"] == "Major"

    def test_format_informational(self) -> None:
        result = format_recommendation(
            drug="tamoxifen",
            gene="CYP2D6",
            phenotype="Normal Metabolizer",
            guideline={
                "recommendation": "Standard dosing.",
                "severity": "None",
            },
        )
        assert result["urgency"] == "FOR INFORMATION"


class TestGenerateSummaryTable:
    """Tests for summary table generation."""

    def test_table_from_genotype_results(self) -> None:
        results = [
            {"gene": "CYP2D6", "diplotype": "*1/*4", "phenotype": "Intermediate Metabolizer"},
            {"gene": "CYP2C19", "diplotype": "*1/*1", "phenotype": "Normal Metabolizer"},
        ]
        table = generate_summary_table(results)
        assert len(table) == 2
        assert table[0]["Gene"] == "CYP2D6"
        assert table[1]["Gene"] == "CYP2C19"

    def test_table_empty_input(self) -> None:
        table = generate_summary_table([])
        assert table == []


class TestExportReport:
    """Tests for report export to file."""

    def test_export_text(self, tmp_path: Path) -> None:
        report = generate_clinical_report(
            {"patient_id": "P001"},
            {"CYP2D6": {"diplotype": "*1/*1"}},
        )
        output_file = tmp_path / "report.txt"
        text = export_report(report, format="text", output_path=output_file)
        assert output_file.exists()
        assert "PHARMACOGENOMIC CLINICAL REPORT" in text

    def test_export_json(self, tmp_path: Path) -> None:
        report = generate_clinical_report(
            {"patient_id": "P002"},
            {"CYP2D6": {"diplotype": "*1/*4"}},
        )
        output_file = tmp_path / "report.json"
        text = export_report(report, format="json", output_path=output_file)
        assert output_file.exists()
        parsed = json.loads(text)
        assert "header" in parsed

    def test_export_html(self, tmp_path: Path) -> None:
        report = generate_clinical_report(
            {"patient_id": "P003"},
            {"CYP2D6": {"diplotype": "*4/*4"}},
            drugs=["codeine"],
        )
        output_file = tmp_path / "report.html"
        text = export_report(report, format="html", output_path=output_file)
        assert output_file.exists()
        assert "<html" in text

    def test_export_invalid_format_raises(self) -> None:
        with pytest.raises(ValueError, match="Unsupported"):
            export_report({}, format="pdf")


class TestAddDisclaimer:
    """Tests for adding disclaimers to reports."""

    def test_add_default_disclaimer(self) -> None:
        report: dict = {"header": {}}
        result = add_disclaimer(report)
        assert "disclaimer" in result
        assert "NOT a substitute" in result["disclaimer"]

    def test_add_custom_disclaimer(self) -> None:
        report: dict = {"header": {}}
        result = add_disclaimer(report, custom_disclaimer="Custom text here.")
        assert result["disclaimer"] == "Custom text here."
