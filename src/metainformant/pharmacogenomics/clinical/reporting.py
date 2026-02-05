"""Clinical pharmacogenomic report generation.

Generates comprehensive clinical reports for pharmacogenomic testing results,
including patient genotype summaries, drug-specific recommendations, metabolizer
status classifications, and clinical disclaimers.

Report formats supported: text (plain text), HTML, and JSON.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from metainformant.core.utils.logging import get_logger

from ..alleles.phenotype import MetabolizerPhenotype, classify_phenotype
from ..annotations.cpic import get_dosing_recommendation
from .drug_interaction import (
    DrugRecommendation,
    InteractionSeverity,
    analyze_drug_gene_interactions,
    calculate_interaction_severity,
)

logger = get_logger(__name__)


_CLINICAL_DISCLAIMER = (
    "DISCLAIMER: This pharmacogenomic report is generated for informational purposes and "
    "is intended to assist healthcare providers in making clinical decisions. It is NOT a "
    "substitute for professional medical judgment. All treatment decisions should be made "
    "by qualified healthcare providers considering the complete clinical picture, including "
    "patient history, current medications, comorbidities, and other relevant factors.\n\n"
    "Pharmacogenomic test results represent one component of personalized medicine. Drug "
    "response is influenced by many factors beyond genetics, including age, weight, organ "
    "function, diet, drug interactions, and adherence. Results should be interpreted in "
    "the context of the individual patient.\n\n"
    "Gene-drug associations and dosing recommendations are based on published CPIC guidelines "
    "and peer-reviewed literature available at the time of report generation. Guidelines are "
    "updated periodically as new evidence emerges.\n\n"
    "This report does not replace the need for therapeutic drug monitoring where indicated."
)


def generate_clinical_report(
    patient_data: dict[str, Any],
    genotypes: dict[str, dict[str, Any]],
    drugs: list[str] | None = None,
) -> dict[str, Any]:
    """Generate a comprehensive clinical pharmacogenomic report.

    Analyzes patient genotype data across all pharmacogenes and generates
    a structured report with phenotype classifications, drug-specific
    recommendations, and clinical action items.

    Args:
        patient_data: Patient demographics dictionary with optional keys:
            - "patient_id": Unique identifier
            - "name": Patient name (optional, may be de-identified)
            - "dob": Date of birth
            - "sex": Biological sex
            - "ethnicity": Self-reported ethnicity
            - "ordering_provider": Ordering provider name
        genotypes: Dictionary mapping gene symbol -> genotype info:
            - "diplotype": Diplotype string (e.g., "*1/*4")
            - "phenotype": Optional pre-determined phenotype
            - "activity_score": Optional pre-computed activity score
        drugs: Optional list of drugs to generate specific recommendations for.
            If None, generates recommendations for all drugs with CPIC guidelines.

    Returns:
        Comprehensive report dictionary with sections:
            - "header": Report metadata
            - "patient_info": Sanitized patient demographics
            - "genotype_results": Per-gene genotype and phenotype results
            - "drug_recommendations": Per-drug clinical recommendations
            - "interaction_summary": Overall interaction risk assessment
            - "clinical_actions": Prioritized action items
            - "disclaimer": Clinical disclaimer text
    """
    report_timestamp = datetime.now(timezone.utc).isoformat()

    # ── Header ────────────────────────────────────────────────────────────
    header = {
        "report_type": "Pharmacogenomic Clinical Report",
        "report_version": "1.0",
        "generated_at": report_timestamp,
        "generator": "METAINFORMANT Pharmacogenomics Module",
    }

    # ── Patient Info (sanitized) ──────────────────────────────────────────
    patient_info = {
        "patient_id": patient_data.get("patient_id", "UNKNOWN"),
        "sex": patient_data.get("sex", "Not specified"),
        "ethnicity": patient_data.get("ethnicity", "Not specified"),
        "ordering_provider": patient_data.get("ordering_provider", "Not specified"),
    }

    # ── Genotype Results ──────────────────────────────────────────────────
    genotype_results: list[dict[str, Any]] = []

    for gene, geno_info in genotypes.items():
        gene_upper = gene.upper()
        diplotype_str = geno_info.get("diplotype", "")

        if diplotype_str:
            try:
                pheno_result = classify_phenotype(diplotype_str, gene_upper)
            except (ValueError, KeyError) as e:
                logger.warning("Could not classify phenotype for %s %s: %s", gene_upper, diplotype_str, e)
                pheno_result = {
                    "gene": gene_upper,
                    "diplotype": diplotype_str,
                    "activity_score": geno_info.get("activity_score", None),
                    "phenotype": geno_info.get("phenotype", MetabolizerPhenotype.INDETERMINATE),
                    "phenotype_abbreviation": "IND",
                    "clinical_significance": "Unable to determine phenotype automatically.",
                }
        else:
            # Use provided phenotype directly
            phenotype = geno_info.get("phenotype", "Indeterminate")
            pheno_result = {
                "gene": gene_upper,
                "diplotype": "Not provided",
                "activity_score": geno_info.get("activity_score", None),
                "phenotype": phenotype,
                "phenotype_abbreviation": _get_abbrev(phenotype),
                "clinical_significance": f"Phenotype provided directly: {phenotype}",
            }

        genotype_results.append(pheno_result)

    # ── Drug Recommendations ──────────────────────────────────────────────
    drug_recommendations: list[dict[str, Any]] = []

    if drugs:
        interactions = analyze_drug_gene_interactions(drugs, genotypes)
        for rec in interactions:
            drug_rec = format_recommendation(
                rec.drug, rec.gene, rec.phenotype,
                {
                    "recommendation": rec.recommendation,
                    "evidence_level": rec.evidence_level,
                    "source": rec.source,
                    "severity": rec.severity.value,
                    "alternatives": rec.alternatives,
                }
            )
            drug_recommendations.append(drug_rec)
    else:
        interactions = []

    # ── Interaction Summary ───────────────────────────────────────────────
    interaction_summary = calculate_interaction_severity(interactions)

    # ── Clinical Actions ──────────────────────────────────────────────────
    clinical_actions: list[dict[str, Any]] = []
    action_priority = 1

    for rec in sorted(interactions, key=lambda r: (
        0 if r.severity == InteractionSeverity.MAJOR else
        1 if r.severity == InteractionSeverity.MODERATE else 2
    )):
        if rec.severity.requires_action:
            clinical_actions.append({
                "priority": action_priority,
                "drug": rec.drug,
                "gene": rec.gene,
                "severity": rec.severity.value,
                "action": rec.recommendation,
                "alternatives": rec.alternatives,
            })
            action_priority += 1

    # ── Assemble Report ───────────────────────────────────────────────────
    report = {
        "header": header,
        "patient_info": patient_info,
        "genotype_results": genotype_results,
        "drug_recommendations": drug_recommendations,
        "interaction_summary": interaction_summary,
        "clinical_actions": clinical_actions,
        "disclaimer": _CLINICAL_DISCLAIMER,
    }

    logger.info(
        "Generated clinical report: %d genes, %d drug recommendations, %d clinical actions",
        len(genotype_results),
        len(drug_recommendations),
        len(clinical_actions),
    )

    return report


def _get_abbrev(phenotype: Any) -> str:
    """Get phenotype abbreviation from various input types."""
    if isinstance(phenotype, MetabolizerPhenotype):
        return phenotype.abbreviation
    s = str(phenotype).lower()
    mapping = {
        "poor": "PM", "intermediate": "IM", "normal": "NM",
        "rapid": "RM", "ultrarapid": "UM",
    }
    for key, abbrev in mapping.items():
        if key in s:
            return abbrev
    return "IND"


def format_recommendation(
    drug: str,
    gene: str,
    phenotype: str,
    guideline: dict[str, Any],
) -> dict[str, Any]:
    """Format a single drug-gene recommendation for the clinical report.

    Args:
        drug: Drug name
        gene: Gene symbol
        phenotype: Metabolizer phenotype string
        guideline: Guideline data with recommendation, evidence_level, etc.

    Returns:
        Formatted recommendation dictionary
    """
    severity = guideline.get("severity", "None")
    recommendation_text = guideline.get("recommendation", "No specific recommendation available.")
    evidence = guideline.get("evidence_level", "")
    source = guideline.get("source", "CPIC")
    alternatives = guideline.get("alternatives", [])

    # Determine urgency indicator
    if severity == "Major":
        urgency = "ACTION REQUIRED"
    elif severity == "Moderate":
        urgency = "ATTENTION RECOMMENDED"
    else:
        urgency = "FOR INFORMATION"

    result = {
        "drug": drug,
        "gene": gene,
        "phenotype": phenotype,
        "recommendation": recommendation_text,
        "evidence_level": evidence,
        "source": source,
        "severity": severity,
        "urgency": urgency,
        "alternatives": alternatives,
    }

    return result


def generate_summary_table(
    results: list[dict[str, Any]],
) -> list[dict[str, str]]:
    """Generate a summary table of all pharmacogenomic findings.

    Creates a condensed tabular view of genotype results and drug
    recommendations suitable for display or embedding in reports.

    Args:
        results: List of genotype result or recommendation dictionaries

    Returns:
        List of row dictionaries with standardized columns:
            - "Gene": Gene symbol
            - "Diplotype": Diplotype string
            - "Phenotype": Metabolizer phenotype
            - "Drug": Drug name (if recommendation)
            - "Recommendation": Brief recommendation
            - "Severity": Interaction severity
    """
    table: list[dict[str, str]] = []

    for entry in results:
        row: dict[str, str] = {}

        # Handle genotype results
        if "gene" in entry:
            row["Gene"] = str(entry.get("gene", ""))
        elif "Gene" in entry:
            row["Gene"] = str(entry.get("Gene", ""))
        else:
            row["Gene"] = ""

        row["Diplotype"] = str(entry.get("diplotype", entry.get("Diplotype", "")))

        phenotype = entry.get("phenotype", entry.get("Phenotype", ""))
        if isinstance(phenotype, MetabolizerPhenotype):
            row["Phenotype"] = phenotype.value
        else:
            row["Phenotype"] = str(phenotype)

        row["Drug"] = str(entry.get("drug", entry.get("Drug", "")))
        row["Recommendation"] = str(entry.get("recommendation", entry.get("Recommendation", "")))[:200]
        row["Severity"] = str(entry.get("severity", entry.get("Severity", "None")))

        table.append(row)

    return table


def export_report(
    report: dict[str, Any],
    format: str = "text",
    output_path: str | Path | None = None,
) -> str:
    """Export a clinical report to the specified format.

    Args:
        report: Report dictionary (from generate_clinical_report)
        format: Output format - "text", "html", or "json"
        output_path: Optional file path to write the report to

    Returns:
        Formatted report string

    Raises:
        ValueError: If format is not supported
    """
    if format == "text":
        output = _format_text_report(report)
    elif format == "html":
        output = _format_html_report(report)
    elif format == "json":
        output = _format_json_report(report)
    else:
        raise ValueError(f"Unsupported report format: '{format}'. Use 'text', 'html', or 'json'.")

    if output_path is not None:
        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(output, encoding="utf-8")
        logger.info("Exported %s report to %s (%d bytes)", format, output_path, len(output))

    return output


def _format_text_report(report: dict[str, Any]) -> str:
    """Format report as plain text."""
    lines: list[str] = []

    # Header
    header = report.get("header", {})
    lines.append("=" * 80)
    lines.append("PHARMACOGENOMIC CLINICAL REPORT")
    lines.append("=" * 80)
    lines.append(f"Generated: {header.get('generated_at', 'N/A')}")
    lines.append(f"Version: {header.get('report_version', '1.0')}")
    lines.append("")

    # Patient Info
    patient = report.get("patient_info", {})
    lines.append("-" * 40)
    lines.append("PATIENT INFORMATION")
    lines.append("-" * 40)
    lines.append(f"Patient ID: {patient.get('patient_id', 'N/A')}")
    lines.append(f"Sex: {patient.get('sex', 'N/A')}")
    lines.append(f"Ethnicity: {patient.get('ethnicity', 'N/A')}")
    lines.append(f"Ordering Provider: {patient.get('ordering_provider', 'N/A')}")
    lines.append("")

    # Genotype Results
    lines.append("-" * 40)
    lines.append("GENOTYPE RESULTS")
    lines.append("-" * 40)
    for result in report.get("genotype_results", []):
        gene = result.get("gene", "")
        diplotype = result.get("diplotype", "")
        phenotype = result.get("phenotype", "")
        if isinstance(phenotype, MetabolizerPhenotype):
            phenotype = phenotype.value
        abbrev = result.get("phenotype_abbreviation", "")
        score = result.get("activity_score", "N/A")
        lines.append(f"  {gene}: {diplotype} | {phenotype} ({abbrev}) | Activity Score: {score}")
        significance = result.get("clinical_significance", "")
        if significance:
            lines.append(f"    -> {significance}")
    lines.append("")

    # Drug Recommendations
    drug_recs = report.get("drug_recommendations", [])
    if drug_recs:
        lines.append("-" * 40)
        lines.append("DRUG RECOMMENDATIONS")
        lines.append("-" * 40)
        for rec in drug_recs:
            urgency = rec.get("urgency", "")
            lines.append(f"  [{urgency}] {rec.get('drug', '')} / {rec.get('gene', '')} ({rec.get('phenotype', '')})")
            lines.append(f"    Recommendation: {rec.get('recommendation', '')}")
            lines.append(f"    Evidence: CPIC Level {rec.get('evidence_level', 'N/A')}")
            alts = rec.get("alternatives", [])
            if alts:
                lines.append(f"    Alternatives: {', '.join(alts)}")
            lines.append("")

    # Clinical Actions
    actions = report.get("clinical_actions", [])
    if actions:
        lines.append("-" * 40)
        lines.append("CLINICAL ACTION ITEMS (PRIORITIZED)")
        lines.append("-" * 40)
        for action in actions:
            lines.append(
                f"  {action['priority']}. [{action['severity']}] {action['drug']}/{action['gene']}: "
                f"{action['action']}"
            )
            alts = action.get("alternatives", [])
            if alts:
                lines.append(f"     Alternatives: {', '.join(alts)}")
        lines.append("")

    # Interaction Summary
    summary = report.get("interaction_summary", {})
    if summary:
        lines.append("-" * 40)
        lines.append("INTERACTION SUMMARY")
        lines.append("-" * 40)
        lines.append(f"  Risk Level: {summary.get('risk_level', 'N/A')}")
        lines.append(f"  Total Interactions: {summary.get('total_interactions', 0)}")
        lines.append(f"  Major: {summary.get('major_count', 0)}")
        lines.append(f"  Moderate: {summary.get('moderate_count', 0)}")
        lines.append(f"  Minor: {summary.get('minor_count', 0)}")
        lines.append(f"  Summary: {summary.get('summary', '')}")
        lines.append("")

    # Disclaimer
    lines.append("=" * 80)
    lines.append("DISCLAIMER")
    lines.append("=" * 80)
    lines.append(report.get("disclaimer", _CLINICAL_DISCLAIMER))
    lines.append("")

    return "\n".join(lines)


def _format_html_report(report: dict[str, Any]) -> str:
    """Format report as HTML."""
    html_parts: list[str] = []

    html_parts.append("<!DOCTYPE html>")
    html_parts.append('<html lang="en">')
    html_parts.append("<head>")
    html_parts.append('<meta charset="utf-8">')
    html_parts.append("<title>Pharmacogenomic Clinical Report</title>")
    html_parts.append("<style>")
    html_parts.append("body { font-family: 'Segoe UI', Arial, sans-serif; margin: 2em; line-height: 1.6; }")
    html_parts.append("h1 { color: #1a1a2e; border-bottom: 2px solid #16213e; padding-bottom: 0.5em; }")
    html_parts.append("h2 { color: #16213e; margin-top: 1.5em; }")
    html_parts.append("table { border-collapse: collapse; width: 100%; margin: 1em 0; }")
    html_parts.append("th, td { border: 1px solid #ccc; padding: 0.5em; text-align: left; }")
    html_parts.append("th { background: #16213e; color: white; }")
    html_parts.append("tr:nth-child(even) { background: #f2f2f2; }")
    html_parts.append(".major { color: #d32f2f; font-weight: bold; }")
    html_parts.append(".moderate { color: #f57c00; font-weight: bold; }")
    html_parts.append(".minor { color: #388e3c; }")
    html_parts.append(".disclaimer { background: #fff3cd; border: 1px solid #ffc107; padding: 1em; "
                       "margin-top: 2em; font-size: 0.9em; }")
    html_parts.append(".action-item { background: #ffebee; border-left: 4px solid #d32f2f; padding: 0.5em 1em; margin: 0.5em 0; }")
    html_parts.append("</style>")
    html_parts.append("</head><body>")

    # Header
    header = report.get("header", {})
    html_parts.append("<h1>Pharmacogenomic Clinical Report</h1>")
    html_parts.append(f"<p>Generated: {_html_escape(header.get('generated_at', 'N/A'))}</p>")

    # Patient Info
    patient = report.get("patient_info", {})
    html_parts.append("<h2>Patient Information</h2>")
    html_parts.append("<table>")
    for key, label in [("patient_id", "Patient ID"), ("sex", "Sex"), ("ethnicity", "Ethnicity"),
                        ("ordering_provider", "Ordering Provider")]:
        html_parts.append(f"<tr><th>{label}</th><td>{_html_escape(str(patient.get(key, 'N/A')))}</td></tr>")
    html_parts.append("</table>")

    # Genotype Results
    html_parts.append("<h2>Genotype Results</h2>")
    html_parts.append("<table>")
    html_parts.append("<tr><th>Gene</th><th>Diplotype</th><th>Phenotype</th><th>Activity Score</th><th>Clinical Significance</th></tr>")
    for result in report.get("genotype_results", []):
        phenotype = result.get("phenotype", "")
        if isinstance(phenotype, MetabolizerPhenotype):
            phenotype = phenotype.value
        html_parts.append(
            f"<tr><td>{_html_escape(str(result.get('gene', '')))}</td>"
            f"<td>{_html_escape(str(result.get('diplotype', '')))}</td>"
            f"<td>{_html_escape(str(phenotype))}</td>"
            f"<td>{result.get('activity_score', 'N/A')}</td>"
            f"<td>{_html_escape(str(result.get('clinical_significance', '')))}</td></tr>"
        )
    html_parts.append("</table>")

    # Drug Recommendations
    drug_recs = report.get("drug_recommendations", [])
    if drug_recs:
        html_parts.append("<h2>Drug Recommendations</h2>")
        html_parts.append("<table>")
        html_parts.append("<tr><th>Drug</th><th>Gene</th><th>Phenotype</th><th>Severity</th>"
                         "<th>Recommendation</th><th>Alternatives</th></tr>")
        for rec in drug_recs:
            severity = rec.get("severity", "None")
            sev_class = severity.lower() if severity.lower() in ("major", "moderate", "minor") else ""
            alts = ", ".join(rec.get("alternatives", []))
            html_parts.append(
                f"<tr><td>{_html_escape(str(rec.get('drug', '')))}</td>"
                f"<td>{_html_escape(str(rec.get('gene', '')))}</td>"
                f"<td>{_html_escape(str(rec.get('phenotype', '')))}</td>"
                f'<td class="{sev_class}">{_html_escape(str(severity))}</td>'
                f"<td>{_html_escape(str(rec.get('recommendation', '')))}</td>"
                f"<td>{_html_escape(alts)}</td></tr>"
            )
        html_parts.append("</table>")

    # Clinical Actions
    actions = report.get("clinical_actions", [])
    if actions:
        html_parts.append("<h2>Clinical Action Items</h2>")
        for action in actions:
            html_parts.append(
                f'<div class="action-item">'
                f'<strong>#{action["priority"]} [{action["severity"]}] '
                f'{_html_escape(action["drug"])}/{_html_escape(action["gene"])}:</strong> '
                f'{_html_escape(action["action"])}</div>'
            )

    # Disclaimer
    html_parts.append(f'<div class="disclaimer"><h3>Disclaimer</h3><p>{_html_escape(report.get("disclaimer", _CLINICAL_DISCLAIMER))}</p></div>')

    html_parts.append("</body></html>")
    return "\n".join(html_parts)


def _format_json_report(report: dict[str, Any]) -> str:
    """Format report as JSON string."""

    def _serialize(obj: Any) -> Any:
        if isinstance(obj, (MetabolizerPhenotype, InteractionSeverity)):
            return obj.value
        if isinstance(obj, DrugRecommendation):
            return {
                "drug": obj.drug,
                "gene": obj.gene,
                "phenotype": obj.phenotype,
                "recommendation": obj.recommendation,
                "evidence_level": obj.evidence_level,
                "source": obj.source,
                "severity": obj.severity.value,
                "alternatives": obj.alternatives,
            }
        if isinstance(obj, datetime):
            return obj.isoformat()
        raise TypeError(f"Object of type {type(obj)} is not JSON serializable")

    return json.dumps(report, indent=2, default=_serialize, ensure_ascii=False)


def _html_escape(text: str) -> str:
    """Basic HTML escaping for report content."""
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&#39;")
    )


def add_disclaimer(
    report: dict[str, Any],
    custom_disclaimer: str | None = None,
) -> dict[str, Any]:
    """Add or update the clinical disclaimer in a report.

    Args:
        report: Report dictionary
        custom_disclaimer: Optional custom disclaimer text.
            If None, uses the standard METAINFORMANT disclaimer.

    Returns:
        Updated report dictionary with disclaimer
    """
    disclaimer_text = custom_disclaimer if custom_disclaimer else _CLINICAL_DISCLAIMER
    report["disclaimer"] = disclaimer_text
    logger.debug("Added disclaimer to report (%d characters)", len(disclaimer_text))
    return report
