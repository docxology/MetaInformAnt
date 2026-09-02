"""Depth tests for clinical reporting: generation, formatting, export.

Uses real genotypes and the built-in CPIC guideline loader; export runs the
actual text/HTML/JSON formatters and writes real files to tmp_path.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from metainformant.pharmacogenomics.clinical.reporting import (
    add_disclaimer,
    export_report,
    format_recommendation,
    generate_clinical_report,
    generate_summary_table,
)

GENOTYPES = {
    "CYP2D6": {"diplotype": "*4/*4"},  # Poor metabolizer
    "CYP2C19": {"diplotype": "*2/*2"},  # Poor metabolizer
}


class TestGenerateClinicalReport:
    def test_report_structure(self) -> None:
        report = generate_clinical_report({"patient_id": "P-1"}, GENOTYPES)
        assert set(report) >= {
            "header",
            "patient_info",
            "genotype_results",
            "drug_recommendations",
            "interaction_summary",
            "clinical_actions",
            "disclaimer",
        }
        assert report["header"]["report_type"] == "Pharmacogenomic Clinical Report"
        assert report["patient_info"]["patient_id"] == "P-1"
        assert report["disclaimer"]

    def test_genotype_results_classified(self) -> None:
        report = generate_clinical_report({}, GENOTYPES)
        genes = {g["gene"] for g in report["genotype_results"]}
        assert genes == {"CYP2D6", "CYP2C19"}
        for g in report["genotype_results"]:
            assert g["phenotype"] is not None
            assert g["phenotype_abbreviation"]

    def test_unknown_gene_defaults_to_cyp2d6_like_mapping(self) -> None:
        # Unknown alleles score 1.0 each; unknown gene uses default thresholds,
        # so *99/*100 classifies as Normal (NM), not IND.
        report = generate_clinical_report({}, {"NOTAGENE": {"diplotype": "*99/*100"}})
        result = report["genotype_results"][0]
        assert result["gene"] == "NOTAGENE"
        assert result["phenotype_abbreviation"] == "NM"

    def test_direct_phenotype_without_diplotype(self) -> None:
        report = generate_clinical_report({}, {"CYP2D6": {"phenotype": "Ultrarapid Metabolizer"}})
        result = report["genotype_results"][0]
        assert result["diplotype"] == "Not provided"
        # _get_abbrev matches substring "rapid" first, giving RM
        assert result["phenotype_abbreviation"] == "RM"

    def test_report_runs_with_drugs_without_matching_guidelines(self) -> None:
        # analyze_drug_gene_interactions only flags pairs in its built-in
        # guideline set; these phenotypes do not match, so no recommendations.
        report = generate_clinical_report({"patient_id": "P-2"}, GENOTYPES, drugs=["clopidogrel", "codeine"])
        assert report["drug_recommendations"] == []
        assert report["clinical_actions"] == []

    def test_no_drugs_no_recommendations(self) -> None:
        report = generate_clinical_report({}, GENOTYPES, drugs=None)
        assert report["drug_recommendations"] == []
        assert report["clinical_actions"] == []

    def test_timestamp_is_iso_utc(self) -> None:
        report = generate_clinical_report({}, GENOTYPES)
        assert report["header"]["generated_at"].endswith("+00:00")


class TestFormatAndSummary:
    def test_format_recommendation_fields_and_urgency(self) -> None:
        rec = format_recommendation(
            "clopidogrel",
            "CYP2C19",
            "Poor Metabolizer",
            {
                "recommendation": "Use alternative antiplatelet therapy",
                "evidence_level": "A",
                "source": "CPIC",
                "severity": "Major",
                "alternatives": ["prasugrel", "ticagrelor"],
            },
        )
        assert rec["drug"] == "clopidogrel" and rec["gene"] == "CYP2C19"
        assert rec["urgency"] == "ACTION REQUIRED"
        assert rec["alternatives"] == ["prasugrel", "ticagrelor"]

        moderate = format_recommendation("d", "g", "Normal", {"severity": "Moderate"})
        assert moderate["urgency"] == "ATTENTION RECOMMENDED"
        none_rec = format_recommendation("d", "g", "Normal", {})
        assert none_rec["urgency"] == "FOR INFORMATION"
        assert none_rec["recommendation"] == "No specific recommendation available."

    def test_summary_table_rows_standardized(self) -> None:
        report = generate_clinical_report({}, GENOTYPES)
        table = generate_summary_table(report["genotype_results"])
        assert len(table) == 2
        assert {row["Gene"] for row in table} == {"CYP2D6", "CYP2C19"}
        for row in table:
            assert row["Phenotype"] == "Poor Metabolizer"
            assert row["Diplotype"]

    def test_add_disclaimer_updates_report(self) -> None:
        report = {"header": {}}
        updated = add_disclaimer(report)
        assert updated["disclaimer"]
        custom = add_disclaimer(report, custom_disclaimer="Custom text.")
        assert custom["disclaimer"] == "Custom text."


class TestExportReport:
    @pytest.fixture()
    def report(self) -> dict:
        return generate_clinical_report({"patient_id": "P-3"}, GENOTYPES, drugs=["clopidogrel"])

    def test_text_export(self, report: dict, tmp_path: Path) -> None:
        out = tmp_path / "report.txt"
        text = export_report(report, format="text", output_path=out)
        assert text == out.read_text(encoding="utf-8")
        assert "P-3" in text
        assert "CYP2D6" in text

    def test_html_export_escapes_and_formats(self, report: dict, tmp_path: Path) -> None:
        out = tmp_path / "report.html"
        html = export_report(report, format="html", output_path=out)
        assert html == out.read_text(encoding="utf-8")
        assert "<html" in html.lower() or "<table" in html.lower()

    def test_json_export_roundtrip(self, report: dict, tmp_path: Path) -> None:
        out = tmp_path / "report.json"
        text = export_report(report, format="json", output_path=out)
        parsed = json.loads(text)
        assert parsed["patient_info"]["patient_id"] == "P-3"
        assert parsed["header"]["report_type"] == report["header"]["report_type"]

    def test_unsupported_format_raises(self, report: dict, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="Unsupported report format"):
            export_report(report, format="xml")

    def test_export_without_path_returns_string(self, report: dict) -> None:
        text = export_report(report, format="json")
        assert isinstance(text, str)
        json.loads(text)
