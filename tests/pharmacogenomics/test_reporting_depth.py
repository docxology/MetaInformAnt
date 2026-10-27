"""Depth tests for clinical reporting: generation, formatting, export.

Uses real genotypes and the built-in CPIC guideline loader; export runs the
actual text/HTML/JSON formatters and writes real files to tmp_path.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from metainformant.pharmacogenomics.annotations.cpic import load_cpic_guidelines
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
        # "ultrarapid" contains "rapid"; it must still be reported as UM, not RM
        assert result["phenotype_abbreviation"] == "UM"

    def test_drug_recommendations_follow_the_classified_phenotypes(self) -> None:
        # GENOTYPES give only diplotypes (no "phenotype" key): *4/*4 is a CYP2D6 poor metabolizer and
        # *2/*2 a CYP2C19 poor metabolizer, so both requested drugs carry a major recommendation.
        report = generate_clinical_report({"patient_id": "P-2"}, GENOTYPES, drugs=["clopidogrel", "codeine"])
        recs = {rec["drug"]: rec for rec in report["drug_recommendations"]}
        assert set(recs) == {"clopidogrel", "codeine"}
        assert recs["clopidogrel"]["gene"] == "CYP2C19" and recs["clopidogrel"]["phenotype"] == "Poor Metabolizer"
        assert recs["codeine"]["gene"] == "CYP2D6" and recs["codeine"]["severity"] == "Major"
        assert {action["drug"] for action in report["clinical_actions"]} == {"clopidogrel", "codeine"}

    def test_drugs_none_analyses_every_cpic_drug(self) -> None:
        report = generate_clinical_report({}, GENOTYPES, drugs=None)
        genes = set(GENOTYPES)
        expected = {g["drug"] for g in load_cpic_guidelines() if g["gene"] in genes}
        assert {rec["drug"] for rec in report["drug_recommendations"]} == expected
        assert expected >= {"codeine", "tramadol", "tamoxifen", "clopidogrel", "voriconazole"}
        assert len(report["clinical_actions"]) == len(expected)

    def test_timestamp_is_iso_utc(self) -> None:
        report = generate_clinical_report({}, GENOTYPES)
        assert report["header"]["generated_at"].endswith("+00:00")


class TestDrugRecommendations:
    """generate_clinical_report documents ``drugs=None`` as "all drugs with CPIC guidelines" and a
    genotype entry as a diplotype and/or a phenotype; both must produce recommendations."""

    PM_GENOTYPES = {"CYP2D6": {"diplotype": "*4/*4"}, "CYP2C19": {"diplotype": "*2/*2"}}
    NM_GENOTYPES = {"CYP2D6": {"diplotype": "*1/*1"}, "CYP2C19": {"diplotype": "*1/*1"}}

    def test_diplotype_only_gene_gets_recommendations_for_requested_drugs(self) -> None:
        report = generate_clinical_report({}, {"CYP2D6": {"diplotype": "*4/*4"}}, drugs=["codeine"])
        assert len(report["drug_recommendations"]) == 1
        rec = report["drug_recommendations"][0]
        assert (rec["drug"], rec["gene"], rec["phenotype"]) == ("codeine", "CYP2D6", "Poor Metabolizer")
        assert rec["severity"] == "Major" and rec["urgency"] == "ACTION REQUIRED"
        assert rec["alternatives"]
        assert report["clinical_actions"][0]["priority"] == 1
        assert report["interaction_summary"]["major_count"] == 1

    def test_drugs_none_covers_every_guideline_drug_of_the_genotyped_genes(self) -> None:
        report = generate_clinical_report({}, self.PM_GENOTYPES)
        recs = report["drug_recommendations"]
        assert [(r["drug"], r["gene"]) for r in recs] == [
            ("codeine", "CYP2D6"),
            ("tramadol", "CYP2D6"),
            ("tamoxifen", "CYP2D6"),
            ("clopidogrel", "CYP2C19"),
            ("voriconazole", "CYP2C19"),
        ]
        assert all(r["severity"] == "Major" for r in recs)
        assert [a["priority"] for a in report["clinical_actions"]] == [1, 2, 3, 4, 5]
        assert report["interaction_summary"]["overall_severity"] == "Major"

    def test_drugs_none_matches_naming_the_same_drugs(self) -> None:
        all_drugs = list(dict.fromkeys(g["drug"] for g in load_cpic_guidelines()))
        implicit = generate_clinical_report({}, self.PM_GENOTYPES)
        explicit = generate_clinical_report({}, self.PM_GENOTYPES, drugs=all_drugs)
        assert implicit["drug_recommendations"] == explicit["drug_recommendations"]
        assert implicit["clinical_actions"] == explicit["clinical_actions"]

    def test_empty_drug_list_requests_no_recommendations(self) -> None:
        report = generate_clinical_report({}, self.PM_GENOTYPES, drugs=[])
        assert report["drug_recommendations"] == []
        assert report["clinical_actions"] == []
        assert report["interaction_summary"]["total_interactions"] == 0

    def test_normal_metabolizers_get_informational_recommendations_and_no_actions(self) -> None:
        report = generate_clinical_report({}, self.NM_GENOTYPES)
        recs = report["drug_recommendations"]
        assert {r["drug"] for r in recs} == {"codeine", "tramadol", "tamoxifen", "clopidogrel", "voriconazole"}
        assert all(r["urgency"] == "FOR INFORMATION" and r["severity"] == "None" for r in recs)
        assert report["clinical_actions"] == []

    def test_recommendation_follows_the_corrected_phenotype(self) -> None:
        # *1/*4 is an intermediate metabolizer (it used to be reported as normal), and the codeine
        # recommendation is the intermediate-metabolizer one.
        report = generate_clinical_report({}, {"CYP2D6": {"diplotype": "*1/*4"}}, drugs=["codeine"])
        assert report["genotype_results"][0]["phenotype_abbreviation"] == "IM"
        rec = report["drug_recommendations"][0]
        assert rec["phenotype"] == "Intermediate Metabolizer"
        assert "reduced efficacy" in rec["recommendation"]

    @pytest.mark.parametrize(
        ("label", "abbreviation", "phenotype"),
        [
            ("Poor Metabolizer", "PM", "Poor Metabolizer"),
            ("poor", "PM", "Poor Metabolizer"),
            ("PM", "PM", "Poor Metabolizer"),
            ("Ultrarapid Metabolizer", "UM", "Ultrarapid Metabolizer"),
            ("ultra-rapid", "UM", "Ultrarapid Metabolizer"),
            ("UM", "UM", "Ultrarapid Metabolizer"),
        ],
    )
    def test_directly_supplied_phenotype_in_any_form_is_reported_and_analysed(
        self, label: str, abbreviation: str, phenotype: str
    ) -> None:
        report = generate_clinical_report({}, {"CYP2D6": {"phenotype": label}}, drugs=["codeine"])
        result = report["genotype_results"][0]
        assert result["phenotype"] == label  # shown as supplied
        assert result["phenotype_abbreviation"] == abbreviation
        rec = report["drug_recommendations"][0]
        assert (rec["drug"], rec["phenotype"], rec["severity"]) == ("codeine", phenotype, "Major")

    def test_free_text_phenotype_keeps_its_category_abbreviation_without_a_recommendation(self) -> None:
        report = generate_clinical_report({}, {"CYP2D6": {"phenotype": "Likely Poor Metabolizer"}}, drugs=["codeine"])
        assert report["genotype_results"][0]["phenotype_abbreviation"] == "PM"
        assert report["drug_recommendations"] == []

    def test_diplotype_wins_over_a_contradicting_phenotype_key(self) -> None:
        # The classified diplotype (*4/*4, poor) is what the report shows and analyses.
        genotypes = {"CYP2D6": {"diplotype": "*4/*4", "phenotype": "Normal Metabolizer"}}
        report = generate_clinical_report({}, genotypes, drugs=["codeine"])
        assert report["genotype_results"][0]["phenotype_abbreviation"] == "PM"
        assert report["drug_recommendations"][0]["phenotype"] == "Poor Metabolizer"

    def test_gene_without_a_guideline_adds_no_recommendation(self) -> None:
        report = generate_clinical_report({}, {"CYP3A5": {"diplotype": "*1/*3"}})
        assert report["genotype_results"][0]["phenotype_abbreviation"] == "IM"
        assert report["drug_recommendations"] == []
        assert report["clinical_actions"] == []

    def test_unclassifiable_diplotype_is_indeterminate_and_adds_no_recommendation(self) -> None:
        report = generate_clinical_report({}, {"CYP2D6": {"diplotype": "*4"}})
        assert report["genotype_results"][0]["phenotype_abbreviation"] == "IND"
        assert report["drug_recommendations"] == []

    def test_recommendations_reach_the_exported_text(self) -> None:
        text = export_report(generate_clinical_report({"patient_id": "P-9"}, self.PM_GENOTYPES), format="text")
        assert "DRUG RECOMMENDATIONS" in text
        assert "[ACTION REQUIRED] codeine / CYP2D6 (Poor Metabolizer)" in text
        assert "CLINICAL ACTION ITEMS (PRIORITIZED)" in text


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
