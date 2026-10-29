"""Tests for pharmacogenomics clinical: ACMG classification, drug interactions, reporting.

Tests ACMG variant classification enums and logic, evidence aggregation,
ClinVar queries, gnomAD frequency checks, drug-gene interaction analysis,
interaction severity calculation, contraindication checks, polypharmacy
analysis, clinical report generation, recommendation formatting, summary
table generation, and report export.

REAL IMPLEMENTATION -- all tests use real implementations.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from metainformant.pharmacogenomics.alleles.phenotype import MetabolizerPhenotype
from metainformant.pharmacogenomics.clinical.drug_interaction import (
    DrugRecommendation,
    InteractionSeverity,
    _abbreviate_phenotype,
    analyze_drug_gene_interactions,
    calculate_interaction_severity,
    check_contraindications,
    polypharmacy_analysis,
    suggest_alternatives,
)
from metainformant.pharmacogenomics.clinical.pathogenicity import (
    _CLINVAR_DATA,
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
from metainformant.pharmacogenomics.metabolism.metabolizer_status import predict_metabolizer_status


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

    def test_severity_reflects_highest_gene_interaction(self) -> None:
        # warfarin/CYP2C9 PM is major but not contraindicated; the reported
        # severity must use clinical ordering, not string comparison.
        result = check_contraindications("warfarin", "Poor Metabolizer")
        assert result["contraindicated"] is False
        assert result["severity"] == "Major"
        assert "dose modification" in result["overall_recommendation"]

    def test_contraindicated_reports_major_severity(self) -> None:
        result = check_contraindications("voriconazole", "Ultrarapid Metabolizer")
        assert result["contraindicated"] is True
        assert result["severity"] == "Major"
        assert len(result["gene_interactions"]) == 1
        assert result["gene_interactions"][0]["gene"] == "CYP2C19"


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


class TestPhenotypeLabelsInDrugGeneAnalysis:
    """A poor metabolizer must be one whichever way the phenotype is written.

    ``metabolism.predict_metabolizer_status`` returns the bare token ("poor"). That token used to
    become the made-up abbreviation "PO", match nothing, and produce no recommendation at all, or
    "No known contraindication", for a poor metabolizer.
    """

    POOR_FORMS = ["poor", "Poor", " POOR ", "PM", "pm", "Poor Metabolizer", MetabolizerPhenotype.POOR]

    @pytest.mark.parametrize("label", POOR_FORMS, ids=[str(f) for f in POOR_FORMS])
    def test_poor_metabolizer_in_any_form_is_a_major_codeine_interaction(
        self, label: str | MetabolizerPhenotype
    ) -> None:
        recommendations = analyze_drug_gene_interactions(["codeine"], {"CYP2D6": {"phenotype": label}})
        assert len(recommendations) == 1
        assert recommendations[0].severity == InteractionSeverity.MAJOR
        assert recommendations[0].phenotype == "Poor Metabolizer"

    def test_the_token_predict_metabolizer_status_returns_can_be_fed_straight_in(self) -> None:
        status = predict_metabolizer_status({"diplotype": "*4/*4"}, "CYP2D6")["phenotype"]
        assert status == "poor"
        recommendations = analyze_drug_gene_interactions(["codeine"], {"CYP2D6": {"phenotype": status}})
        assert [r.severity for r in recommendations] == [InteractionSeverity.MAJOR]
        assert check_contraindications("codeine", status)["contraindicated"] is True

    @pytest.mark.parametrize("label", POOR_FORMS, ids=[str(f) for f in POOR_FORMS])
    def test_check_contraindications_accepts_the_same_forms(self, label: str | MetabolizerPhenotype) -> None:
        result = check_contraindications("codeine", label)
        assert result["contraindicated"] is True
        assert result["severity"] == "Major"

    @pytest.mark.parametrize("label", ["superrapid", "", "Likely Poor Metabolizer", "xyz"])
    def test_check_contraindications_refuses_a_label_that_is_not_a_phenotype(self, label: str) -> None:
        # It used to answer "No known contraindication for codeine with superrapid."
        with pytest.raises(ValueError, match="Unrecognised metabolizer phenotype"):
            check_contraindications("codeine", label)

    @pytest.mark.parametrize("label", ["Indeterminate", "IND", "indeterminate", MetabolizerPhenotype.INDETERMINATE])
    def test_check_contraindications_refuses_an_indeterminate_phenotype(
        self, label: str | MetabolizerPhenotype
    ) -> None:
        with pytest.raises(ValueError, match="indeterminate phenotype"):
            check_contraindications("codeine", label)

    def test_abbreviations_are_never_made_up_from_the_first_letters(self) -> None:
        # "poor" -> "PO", "Indeterminate" -> "IN" and "Likely Poor Metabolizer" -> "LI" used to be returned
        assert _abbreviate_phenotype("poor") == "PM"
        assert _abbreviate_phenotype("Extensive Metabolizer") == "NM"
        assert _abbreviate_phenotype("Indeterminate") == "IND"
        with pytest.raises(ValueError, match="Unrecognised metabolizer phenotype"):
            _abbreviate_phenotype("Likely Poor Metabolizer")

    @pytest.mark.parametrize("label", [MetabolizerPhenotype.INDETERMINATE, "Likely Poor Metabolizer", "", "xyz"])
    def test_a_gene_without_an_interpretable_phenotype_gets_no_recommendation_and_says_so(
        self, label: str | MetabolizerPhenotype, caplog: pytest.LogCaptureFixture
    ) -> None:
        with caplog.at_level("WARNING"):
            recommendations = analyze_drug_gene_interactions(["codeine"], {"CYP2D6": {"phenotype": label}})
        assert recommendations == []
        warnings = [r.getMessage() for r in caplog.records if r.levelname == "WARNING"]
        assert any("No recommendation for codeine/CYP2D6" in message for message in warnings), warnings

    def test_a_gene_with_no_guideline_for_the_drug_is_skipped_quietly(self, caplog: pytest.LogCaptureFixture) -> None:
        # Nothing to warn about: CYP2C9 has no codeine guideline, whatever its phenotype
        with caplog.at_level("WARNING"):
            assert analyze_drug_gene_interactions(["codeine"], {"CYP2C9": {"phenotype": "xyz"}}) == []
        assert not [r for r in caplog.records if r.levelname == "WARNING"]


class TestFrequencyCriteriaNeedAFrequency:
    """A variant with no allele frequency has not been seen to be absent from controls."""

    def test_a_missing_frequency_meets_none_of_the_frequency_criteria(self) -> None:
        # PM2 used to be met because the missing frequency was read as 0.0
        criteria = apply_acmg_criteria({"consequence": "frameshift"})
        assert criteria["PM2"] is False
        assert criteria["BA1"] is False
        assert criteria["BS1"] is False

    @pytest.mark.parametrize("empty", [None, "", {}])
    def test_so_do_an_empty_or_null_frequency(self, empty: object) -> None:
        criteria = apply_acmg_criteria({"consequence": "frameshift", "allele_frequency": empty})
        assert criteria["PM2"] is False

    def test_a_null_variant_with_no_frequency_data_is_not_upgraded_by_invented_evidence(self) -> None:
        # PVS1 alone is VUS. With a PM2 made from "no data" it used to be Likely Pathogenic.
        result = classify_variant_acmg({"consequence": "frameshift"})
        assert result["criteria_met"] == ["PVS1"]
        assert result["classification"] == ACMGClassification.VUS

    @pytest.mark.parametrize("frequency", [0.0, {"European": 0.0}, {"European": 0.0, "African": 0.00001}])
    def test_a_measured_absence_is_absence(self, frequency: object) -> None:
        criteria = apply_acmg_criteria({"consequence": "frameshift", "allele_frequency": frequency})
        assert criteria["PM2"] is True
        result = classify_variant_acmg({"consequence": "frameshift", "allele_frequency": frequency})
        assert result["classification"] == ACMGClassification.LIKELY_PATHOGENIC

    def test_common_and_intermediate_frequencies_still_meet_the_benign_criteria(self) -> None:
        assert apply_acmg_criteria({"allele_frequency": 0.08})["BA1"] is True
        intermediate = apply_acmg_criteria({"allele_frequency": 0.02})
        assert intermediate["BS1"] is True and intermediate["BA1"] is False and intermediate["PM2"] is False


class TestContradictoryEvidenceIsUncertain:
    """Richards et al. 2015, Table 5: pathogenic and benign criteria that both reach a call give VUS."""

    @pytest.mark.parametrize(
        "criteria",
        [
            {"PVS1": True, "PS1": True, "BA1": True},  # Pathogenic vs Benign
            {"PVS1": True, "PS1": True, "BP4": True, "BP7": True},  # Pathogenic vs Likely Benign
            {"PS1": True, "PS3": True, "BS1": True, "BS3": True},  # Pathogenic (2 strong) vs Benign (2 strong)
            {"PVS1": True, "PM2": True, "BP4": True, "BP7": True},  # Likely Pathogenic vs Likely Benign
            {"PM1": True, "PM2": True, "PM3": True, "BS1": True, "BP4": True},  # LP (3 moderate) vs LB
        ],
        ids=["P-vs-B", "P-vs-LB", "P-vs-B-strong", "LP-vs-LB", "LP-3-moderate-vs-LB"],
    )
    def test_both_sides_met_is_vus(self, criteria: dict[str, bool]) -> None:
        assert aggregate_evidence(criteria) == ACMGClassification.VUS

    def test_one_side_alone_is_unchanged(self) -> None:
        assert aggregate_evidence({"PVS1": True, "PS1": True}) == ACMGClassification.PATHOGENIC
        assert aggregate_evidence({"PVS1": True, "PM2": True}) == ACMGClassification.LIKELY_PATHOGENIC
        assert aggregate_evidence({"BA1": True}) == ACMGClassification.BENIGN
        assert aggregate_evidence({"BS1": True, "BS2": True}) == ACMGClassification.BENIGN
        assert aggregate_evidence({"BS1": True, "BP4": True}) == ACMGClassification.LIKELY_BENIGN

    def test_pathogenic_evidence_short_of_a_call_does_not_cancel_a_benign_call(self) -> None:
        # PM2 alone calls nothing, so it does not contradict two supporting benign criteria
        assert aggregate_evidence({"PM2": True, "BP4": True, "BP7": True}) == ACMGClassification.LIKELY_BENIGN

    def test_classify_variant_acmg_reports_the_contradiction_as_vus(self) -> None:
        result = classify_variant_acmg({}, evidence={"PVS1": True, "PS1": True, "BA1": True})
        assert result["classification"] == ACMGClassification.VUS
        assert result["criteria_met"] == ["BA1", "PS1", "PVS1"]


class TestGnomadFrequencyInputs:
    def test_a_population_the_variant_has_no_frequency_for_raises(self) -> None:
        # It used to be reported as a frequency of 0.0 (no BA1, no BS1), which a typo could produce
        with pytest.raises(ValueError, match="East_Asian") as caught:
            check_gnomad_frequency({"allele_frequency": {"European": 0.2}}, population="East_Asian")
        assert "European" in str(caught.value)  # the populations it does have

    def test_the_population_is_matched_case_insensitively_and_reported_as_stored(self) -> None:
        result = check_gnomad_frequency({"allele_frequency": {"European": 0.2, "African": 0.01}}, population="european")
        assert result["max_frequency"] == 0.2
        assert result["population"] == "European"
        assert result["ba1_triggered"] is True

    def test_a_scalar_frequency_is_the_global_population(self) -> None:
        result = check_gnomad_frequency({"allele_frequency": 0.07})
        assert (result["max_frequency"], result["population"], result["ba1_triggered"]) == (0.07, "Global", True)
        assert check_gnomad_frequency({"allele_frequency": 0.07}, population="Global")["ba1_triggered"] is True

    @pytest.mark.parametrize("variant", [{}, {"allele_frequency": None}, {"allele_frequency": {}}])
    def test_no_frequency_is_reported_as_none_not_zero(self, variant: dict[str, object]) -> None:
        result = check_gnomad_frequency(variant)
        assert result["max_frequency"] is None
        assert result["population"] is None
        assert result["exceeds_threshold"] is False
        assert result["ba1_triggered"] is False and result["bs1_triggered"] is False
        assert result["frequencies"] == {}

    def test_no_frequency_makes_a_named_population_an_error(self) -> None:
        with pytest.raises(ValueError, match="European"):
            check_gnomad_frequency({}, population="European")

    def test_bs1_is_triggered_between_the_threshold_and_ba1(self) -> None:
        result = check_gnomad_frequency({"allele_frequency": {"European": 0.03}}, threshold=0.01)
        assert result["bs1_triggered"] is True and result["ba1_triggered"] is False
        assert result["exceeds_threshold"] is True


class TestQueryClinvarScope:
    """query_clinvar reads a six-record built-in table, which its docstring now says."""

    def test_docstring_states_the_size_of_the_built_in_table(self) -> None:
        doc = " ".join((query_clinvar.__doc__ or "").split())
        claimed = re.search(r"table of (\d+) ClinVar-style records", doc)
        assert claimed is not None and int(claimed.group(1)) == len(_CLINVAR_DATA)
        assert "not in that table" in doc  # None is not "ClinVar has no record"

    def test_module_docstring_does_not_claim_a_local_clinvar_download(self) -> None:
        import metainformant.pharmacogenomics.clinical.pathogenicity as pathogenicity

        doc = " ".join((pathogenicity.__doc__ or "").split())
        assert "not a ClinVar download" in doc

    @pytest.mark.parametrize(
        "identifier", ["rs3892097", " RS3892097 ", "NM_000106.6:c.506-1G>A", "nm_000106.6:c.506-1g>a"]
    )
    def test_lookup_by_rsid_or_hgvs_ignores_case_and_padding(self, identifier: str) -> None:
        record = query_clinvar(identifier)
        assert record is not None and record["gene"] == "CYP2D6"

    @pytest.mark.parametrize("identifier", ["12345", "VCV000012345", "rs000000000"])
    def test_clinvar_variation_ids_and_unlisted_variants_are_not_found(self, identifier: str) -> None:
        assert query_clinvar(identifier) is None

    def test_a_returned_record_is_a_copy(self) -> None:
        record = query_clinvar("rs3892097")
        assert record is not None
        record["classification"] = "Benign"
        again = query_clinvar("rs3892097")
        assert again is not None and again["classification"] == "Pathogenic"
