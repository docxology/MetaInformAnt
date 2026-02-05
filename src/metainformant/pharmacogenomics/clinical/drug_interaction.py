"""Drug-gene interaction analysis for pharmacogenomics.

Provides multi-drug interaction assessment based on pharmacogenomic genotype data,
contraindication checking, severity scoring, polypharmacy analysis, and alternative
drug suggestions. Integrates CPIC guidelines with drug-gene interaction databases
to provide clinically actionable recommendations.

Interaction severity levels:
    MAJOR: Potentially life-threatening or causing permanent damage
    MODERATE: May result in deterioration of patient's clinical status
    MINOR: Minimally clinically significant; may increase side effect frequency
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any

from metainformant.core.utils.logging import get_logger

from ..annotations.cpic import lookup_drug_gene, get_dosing_recommendation, _CPIC_GUIDELINES
from ..alleles.phenotype import MetabolizerPhenotype

logger = get_logger(__name__)


class InteractionSeverity(str, Enum):
    """Drug-gene interaction severity classification."""

    MAJOR = "Major"
    MODERATE = "Moderate"
    MINOR = "Minor"
    NONE = "None"

    @property
    def requires_action(self) -> bool:
        """Return True if this severity level requires clinical action."""
        return self in (InteractionSeverity.MAJOR, InteractionSeverity.MODERATE)


@dataclass
class DrugRecommendation:
    """Represents a pharmacogenomic drug recommendation.

    Attributes:
        drug: Drug name
        gene: Gene symbol
        phenotype: Metabolizer phenotype
        recommendation: Clinical recommendation text
        evidence_level: CPIC evidence level (A, A/B, B, C, D)
        source: Source of the recommendation (e.g., "CPIC", "PharmGKB", "FDA")
        severity: Interaction severity level
        alternatives: Suggested alternative drugs
    """

    drug: str
    gene: str
    phenotype: str
    recommendation: str
    evidence_level: str = ""
    source: str = "CPIC"
    severity: InteractionSeverity = InteractionSeverity.NONE
    alternatives: list[str] = field(default_factory=list)


# ── Drug-gene interaction database ─────────────────────────────────────────────
# Maps (drug, gene, phenotype_category) -> interaction details

_CONTRAINDICATION_DB: dict[tuple[str, str, str], dict[str, Any]] = {
    ("codeine", "CYP2D6", "PM"): {
        "contraindicated": True,
        "severity": InteractionSeverity.MAJOR,
        "reason": "Insufficient conversion to active metabolite morphine. Lack of therapeutic effect.",
        "alternatives": ["morphine", "oxycodone", "hydromorphone", "acetaminophen", "ibuprofen"],
    },
    ("codeine", "CYP2D6", "UM"): {
        "contraindicated": True,
        "severity": InteractionSeverity.MAJOR,
        "reason": "Excessive conversion to morphine. Risk of respiratory depression and death.",
        "alternatives": ["morphine (with dose reduction)", "acetaminophen", "ibuprofen"],
    },
    ("tramadol", "CYP2D6", "PM"): {
        "contraindicated": True,
        "severity": InteractionSeverity.MAJOR,
        "reason": "Insufficient formation of active metabolite O-desmethyltramadol.",
        "alternatives": ["morphine", "oxycodone", "hydromorphone"],
    },
    ("tramadol", "CYP2D6", "UM"): {
        "contraindicated": True,
        "severity": InteractionSeverity.MAJOR,
        "reason": "Excessive formation of active metabolite. Risk of respiratory depression.",
        "alternatives": ["morphine (with dose adjustment)", "non-opioid analgesics"],
    },
    ("clopidogrel", "CYP2C19", "PM"): {
        "contraindicated": True,
        "severity": InteractionSeverity.MAJOR,
        "reason": "Significantly reduced formation of active metabolite. Inadequate platelet inhibition.",
        "alternatives": ["prasugrel", "ticagrelor"],
    },
    ("clopidogrel", "CYP2C19", "IM"): {
        "contraindicated": False,
        "severity": InteractionSeverity.MODERATE,
        "reason": "Reduced formation of active metabolite. Potentially suboptimal platelet inhibition.",
        "alternatives": ["prasugrel", "ticagrelor"],
    },
    ("voriconazole", "CYP2C19", "PM"): {
        "contraindicated": False,
        "severity": InteractionSeverity.MAJOR,
        "reason": "Higher drug exposure. Increased risk of hepatotoxicity and visual disturbances.",
        "alternatives": ["isavuconazole", "posaconazole"],
    },
    ("voriconazole", "CYP2C19", "UM"): {
        "contraindicated": True,
        "severity": InteractionSeverity.MAJOR,
        "reason": "Subtherapeutic drug levels. High risk of treatment failure.",
        "alternatives": ["isavuconazole", "posaconazole"],
    },
    ("fluorouracil", "DPYD", "PM"): {
        "contraindicated": True,
        "severity": InteractionSeverity.MAJOR,
        "reason": "Complete DPD deficiency. Life-threatening toxicity risk including neutropenia, "
        "mucositis, diarrhea, and neurotoxicity.",
        "alternatives": ["Consult oncologist for alternative regimen"],
    },
    ("fluorouracil", "DPYD", "IM"): {
        "contraindicated": False,
        "severity": InteractionSeverity.MAJOR,
        "reason": "Partial DPD deficiency. Increased risk of severe toxicity.",
        "alternatives": ["fluorouracil at 50% dose with monitoring"],
    },
    ("capecitabine", "DPYD", "PM"): {
        "contraindicated": True,
        "severity": InteractionSeverity.MAJOR,
        "reason": "Complete DPD deficiency. Life-threatening toxicity risk.",
        "alternatives": ["Consult oncologist for alternative regimen"],
    },
    ("azathioprine", "TPMT", "PM"): {
        "contraindicated": False,
        "severity": InteractionSeverity.MAJOR,
        "reason": "Homozygous non-functional TPMT. Very high risk of myelosuppression.",
        "alternatives": ["azathioprine at 10% dose", "mycophenolate mofetil"],
    },
    ("azathioprine", "NUDT15", "PM"): {
        "contraindicated": False,
        "severity": InteractionSeverity.MAJOR,
        "reason": "Homozygous non-functional NUDT15. Very high risk of leukopenia.",
        "alternatives": ["azathioprine at 10% dose", "mycophenolate mofetil"],
    },
    ("simvastatin", "SLCO1B1", "PM"): {
        "contraindicated": False,
        "severity": InteractionSeverity.MAJOR,
        "reason": "Markedly increased simvastatin acid exposure. High myopathy risk.",
        "alternatives": ["rosuvastatin", "pravastatin", "fluvastatin"],
    },
    ("simvastatin", "SLCO1B1", "IM"): {
        "contraindicated": False,
        "severity": InteractionSeverity.MODERATE,
        "reason": "Increased simvastatin acid exposure. Elevated myopathy risk.",
        "alternatives": ["rosuvastatin", "pravastatin", "simvastatin (max 20mg)"],
    },
    ("warfarin", "CYP2C9", "PM"): {
        "contraindicated": False,
        "severity": InteractionSeverity.MAJOR,
        "reason": "Greatly reduced warfarin clearance. High bleeding risk at standard doses.",
        "alternatives": ["warfarin (50-80% dose reduction)", "direct oral anticoagulants"],
    },
    ("warfarin", "CYP2C9", "IM"): {
        "contraindicated": False,
        "severity": InteractionSeverity.MODERATE,
        "reason": "Reduced warfarin clearance. Bleeding risk at standard doses.",
        "alternatives": ["warfarin (20-40% dose reduction)", "direct oral anticoagulants"],
    },
    ("tamoxifen", "CYP2D6", "PM"): {
        "contraindicated": False,
        "severity": InteractionSeverity.MAJOR,
        "reason": "Greatly reduced endoxifen formation. Reduced therapeutic efficacy.",
        "alternatives": ["aromatase inhibitor (if postmenopausal)", "tamoxifen (higher dose with monitoring)"],
    },
}

# ── Alternative drug databases by therapeutic class ────────────────────────────

_ALTERNATIVES_DB: dict[str, dict[str, Any]] = {
    "codeine": {
        "therapeutic_class": "Opioid analgesic",
        "alternatives": {
            "morphine": {"affected_by_cyp2d6": False, "notes": "Direct-acting, not a prodrug"},
            "oxycodone": {"affected_by_cyp2d6": True, "notes": "Partially CYP2D6-dependent"},
            "hydromorphone": {"affected_by_cyp2d6": False, "notes": "Not CYP2D6-dependent"},
            "acetaminophen": {"affected_by_cyp2d6": False, "notes": "Non-opioid alternative"},
            "ibuprofen": {"affected_by_cyp2d6": False, "notes": "NSAID alternative"},
        },
    },
    "clopidogrel": {
        "therapeutic_class": "Antiplatelet agent",
        "alternatives": {
            "prasugrel": {"affected_by_cyp2c19": False, "notes": "Not CYP2C19-dependent activation"},
            "ticagrelor": {"affected_by_cyp2c19": False, "notes": "Direct-acting, no CYP2C19 dependence"},
        },
    },
    "simvastatin": {
        "therapeutic_class": "HMG-CoA reductase inhibitor",
        "alternatives": {
            "rosuvastatin": {"affected_by_slco1b1": True, "notes": "Less OATP1B1-dependent than simvastatin"},
            "pravastatin": {"affected_by_slco1b1": True, "notes": "Less myopathy risk"},
            "fluvastatin": {"affected_by_slco1b1": False, "notes": "Different metabolic pathway"},
            "pitavastatin": {"affected_by_slco1b1": True, "notes": "Lower myopathy risk profile"},
        },
    },
    "warfarin": {
        "therapeutic_class": "Anticoagulant",
        "alternatives": {
            "apixaban": {"pgx_affected": False, "notes": "Direct oral anticoagulant, no PGx dependency"},
            "rivaroxaban": {"pgx_affected": False, "notes": "Direct oral anticoagulant"},
            "dabigatran": {"pgx_affected": False, "notes": "Direct thrombin inhibitor"},
            "edoxaban": {"pgx_affected": False, "notes": "Direct oral anticoagulant"},
        },
    },
    "voriconazole": {
        "therapeutic_class": "Azole antifungal",
        "alternatives": {
            "isavuconazole": {"affected_by_cyp2c19": False, "notes": "Less CYP2C19-dependent"},
            "posaconazole": {"affected_by_cyp2c19": False, "notes": "Not CYP2C19-dependent"},
        },
    },
    "tamoxifen": {
        "therapeutic_class": "Selective estrogen receptor modulator",
        "alternatives": {
            "anastrozole": {"affected_by_cyp2d6": False, "notes": "Aromatase inhibitor (postmenopausal only)"},
            "letrozole": {"affected_by_cyp2d6": False, "notes": "Aromatase inhibitor (postmenopausal only)"},
            "exemestane": {"affected_by_cyp2d6": False, "notes": "Steroidal aromatase inhibitor"},
        },
    },
}


def analyze_drug_gene_interactions(
    drugs: list[str],
    genotype_data: dict[str, dict[str, Any]],
) -> list[DrugRecommendation]:
    """Analyze drug-gene interactions for a set of drugs and genotypes.

    For each drug, checks against all available genotype data to identify
    pharmacogenomic interactions and generate clinical recommendations.

    Args:
        drugs: List of drug names to analyze
        genotype_data: Dictionary mapping gene symbol to genotype information.
            Each gene entry should have:
                - "phenotype": MetabolizerPhenotype or string (e.g., "Poor Metabolizer")
                - "diplotype": Optional diplotype string (e.g., "*1/*4")
                - "activity_score": Optional activity score

    Returns:
        List of DrugRecommendation objects for all identified interactions
    """
    recommendations: list[DrugRecommendation] = []

    for drug in drugs:
        drug_lower = drug.lower().strip()

        # Check each gene in the genotype data against CPIC guidelines
        for gene, geno_info in genotype_data.items():
            gene_upper = gene.upper()
            phenotype = geno_info.get("phenotype", "")

            if isinstance(phenotype, MetabolizerPhenotype):
                phenotype_str = phenotype.value
                phenotype_abbrev = phenotype.abbreviation
            else:
                phenotype_str = str(phenotype)
                phenotype_abbrev = _abbreviate_phenotype(phenotype_str)

            # Look up CPIC guideline
            guideline = lookup_drug_gene(drug_lower, gene_upper)
            if guideline is None:
                continue

            # Get specific recommendation for this phenotype
            dosing = get_dosing_recommendation(drug_lower, phenotype_str)

            # Check contraindication database
            contra_key = (drug_lower, gene_upper, phenotype_abbrev)
            contra_info = _CONTRAINDICATION_DB.get(contra_key)

            severity = InteractionSeverity.NONE
            alternatives: list[str] = []
            rec_text = ""

            if contra_info is not None:
                severity = contra_info["severity"]
                alternatives = contra_info.get("alternatives", [])
                rec_text = contra_info["reason"]
            elif dosing is not None:
                rec_text = dosing.get("recommendation", "")
                # Infer severity from recommendation language
                rec_lower = rec_text.lower()
                if any(w in rec_lower for w in ["avoid", "contraindicated", "do not use"]):
                    severity = InteractionSeverity.MAJOR
                elif any(w in rec_lower for w in ["reduce", "lower dose", "monitor", "consider alternative"]):
                    severity = InteractionSeverity.MODERATE
                elif any(w in rec_lower for w in ["standard of care", "per standard"]):
                    severity = InteractionSeverity.NONE
                else:
                    severity = InteractionSeverity.MINOR

            if rec_text or severity != InteractionSeverity.NONE:
                rec = DrugRecommendation(
                    drug=drug,
                    gene=gene_upper,
                    phenotype=phenotype_str,
                    recommendation=rec_text if rec_text else "No specific recommendation available.",
                    evidence_level=guideline.get("cpic_level", ""),
                    source="CPIC",
                    severity=severity,
                    alternatives=alternatives,
                )
                recommendations.append(rec)

                logger.info(
                    "Drug-gene interaction: %s/%s (%s) -> Severity: %s",
                    drug,
                    gene_upper,
                    phenotype_str,
                    severity.value,
                )

    logger.info(
        "Analyzed %d drugs against %d genes: %d interactions found",
        len(drugs),
        len(genotype_data),
        len(recommendations),
    )

    return recommendations


def _abbreviate_phenotype(phenotype: str) -> str:
    """Convert phenotype string to abbreviation.

    Args:
        phenotype: Full phenotype name

    Returns:
        Abbreviation (PM, IM, NM, RM, UM)
    """
    mapping = {
        "poor metabolizer": "PM",
        "intermediate metabolizer": "IM",
        "normal metabolizer": "NM",
        "extensive metabolizer": "NM",
        "rapid metabolizer": "RM",
        "ultrarapid metabolizer": "UM",
    }
    return mapping.get(phenotype.lower().strip(), phenotype.upper()[:2])


def check_contraindications(
    drug: str,
    phenotype: str | MetabolizerPhenotype,
) -> dict[str, Any]:
    """Check if a drug is contraindicated for a given metabolizer phenotype.

    Searches across all gene-drug contraindication records for the specified drug.

    Args:
        drug: Drug name (case-insensitive)
        phenotype: Metabolizer phenotype (string or enum)

    Returns:
        Dictionary with:
            - "contraindicated": Whether drug is contraindicated
            - "severity": Interaction severity
            - "gene_interactions": List of specific gene-level contraindication details
            - "overall_recommendation": Summary recommendation
    """
    drug_lower = drug.lower().strip()

    if isinstance(phenotype, MetabolizerPhenotype):
        phenotype_abbrev = phenotype.abbreviation
        phenotype_str = phenotype.value
    else:
        phenotype_abbrev = _abbreviate_phenotype(str(phenotype))
        phenotype_str = str(phenotype)

    gene_interactions: list[dict[str, Any]] = []
    any_contraindicated = False
    max_severity = InteractionSeverity.NONE

    for (d, g, p), info in _CONTRAINDICATION_DB.items():
        if d == drug_lower and p == phenotype_abbrev:
            gene_interactions.append(
                {
                    "gene": g,
                    "contraindicated": info["contraindicated"],
                    "severity": info["severity"].value,
                    "reason": info["reason"],
                    "alternatives": info.get("alternatives", []),
                }
            )
            if info["contraindicated"]:
                any_contraindicated = True
            if info["severity"].value > max_severity.value:
                max_severity = info["severity"]

    if any_contraindicated:
        overall = f"{drug} is contraindicated for {phenotype_str}. Use alternative therapy."
    elif max_severity == InteractionSeverity.MAJOR:
        overall = f"{drug} requires significant dose modification for {phenotype_str}."
    elif max_severity == InteractionSeverity.MODERATE:
        overall = f"{drug} may require dose adjustment for {phenotype_str}. Monitor closely."
    else:
        overall = f"No known contraindication for {drug} with {phenotype_str}."

    return {
        "contraindicated": any_contraindicated,
        "severity": max_severity.value,
        "gene_interactions": gene_interactions,
        "overall_recommendation": overall,
    }


def calculate_interaction_severity(
    interactions: list[DrugRecommendation],
) -> dict[str, Any]:
    """Calculate overall severity for a set of drug-gene interactions.

    Assesses the combined risk from multiple pharmacogenomic interactions,
    providing an overall severity classification and risk summary.

    Args:
        interactions: List of DrugRecommendation objects

    Returns:
        Dictionary with:
            - "overall_severity": Highest severity level across all interactions
            - "major_count": Number of major interactions
            - "moderate_count": Number of moderate interactions
            - "minor_count": Number of minor interactions
            - "total_interactions": Total interaction count
            - "risk_level": Qualitative risk assessment ("High", "Moderate", "Low")
            - "summary": Text summary of findings
    """
    if not interactions:
        return {
            "overall_severity": InteractionSeverity.NONE.value,
            "major_count": 0,
            "moderate_count": 0,
            "minor_count": 0,
            "total_interactions": 0,
            "risk_level": "Low",
            "summary": "No pharmacogenomic interactions identified.",
        }

    major = sum(1 for i in interactions if i.severity == InteractionSeverity.MAJOR)
    moderate = sum(1 for i in interactions if i.severity == InteractionSeverity.MODERATE)
    minor = sum(1 for i in interactions if i.severity == InteractionSeverity.MINOR)

    if major > 0:
        overall = InteractionSeverity.MAJOR
        risk = "High"
    elif moderate > 0:
        overall = InteractionSeverity.MODERATE
        risk = "Moderate"
    elif minor > 0:
        overall = InteractionSeverity.MINOR
        risk = "Low"
    else:
        overall = InteractionSeverity.NONE
        risk = "Low"

    # Enhanced risk for multiple major interactions
    if major >= 2:
        risk = "Very High"

    summary_parts: list[str] = []
    if major:
        summary_parts.append(f"{major} major interaction(s)")
    if moderate:
        summary_parts.append(f"{moderate} moderate interaction(s)")
    if minor:
        summary_parts.append(f"{minor} minor interaction(s)")

    summary = f"Identified {', '.join(summary_parts)}." if summary_parts else "No significant interactions."

    return {
        "overall_severity": overall.value,
        "major_count": major,
        "moderate_count": moderate,
        "minor_count": minor,
        "total_interactions": len(interactions),
        "risk_level": risk,
        "summary": summary,
    }


def polypharmacy_analysis(
    drug_list: list[str],
    genotype_data: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    """Perform comprehensive polypharmacy pharmacogenomic assessment.

    Analyzes all drugs in a patient's medication list against their genotype
    data to identify all PGx interactions, drug-drug-gene overlaps, and
    cumulative risk.

    Args:
        drug_list: List of all drugs the patient is taking
        genotype_data: Dictionary mapping gene -> genotype info (see analyze_drug_gene_interactions)

    Returns:
        Dictionary with:
            - "total_drugs": Number of drugs analyzed
            - "drugs_with_interactions": Drugs that have PGx interactions
            - "interactions": Full list of DrugRecommendation objects
            - "severity_summary": Interaction severity breakdown
            - "gene_overlap": Genes affected by multiple drugs
            - "high_risk_combinations": Drug combinations with compounding PGx risk
            - "recommendations": Prioritized clinical recommendations
    """
    all_interactions = analyze_drug_gene_interactions(drug_list, genotype_data)
    severity_info = calculate_interaction_severity(all_interactions)

    # Identify drugs with interactions
    drugs_with_ix = set()
    gene_drug_map: dict[str, list[str]] = {}

    for rec in all_interactions:
        if rec.severity != InteractionSeverity.NONE:
            drugs_with_ix.add(rec.drug)

            if rec.gene not in gene_drug_map:
                gene_drug_map[rec.gene] = []
            gene_drug_map[rec.gene].append(rec.drug)

    # Identify gene overlaps (multiple drugs affected by same gene)
    gene_overlap: dict[str, list[str]] = {gene: drugs for gene, drugs in gene_drug_map.items() if len(drugs) > 1}

    # Identify high-risk combinations
    high_risk: list[dict[str, Any]] = []
    for gene, drugs in gene_overlap.items():
        major_drugs = [
            rec.drug for rec in all_interactions if rec.gene == gene and rec.severity == InteractionSeverity.MAJOR
        ]
        if len(major_drugs) >= 2:
            high_risk.append(
                {
                    "gene": gene,
                    "drugs": major_drugs,
                    "risk": "Multiple drugs with major PGx interactions via the same gene pathway",
                }
            )

    # Generate prioritized recommendations
    sorted_interactions = sorted(
        all_interactions,
        key=lambda r: (
            0
            if r.severity == InteractionSeverity.MAJOR
            else (
                1 if r.severity == InteractionSeverity.MODERATE else 2 if r.severity == InteractionSeverity.MINOR else 3
            )
        ),
    )

    prioritized_recs: list[str] = []
    for rec in sorted_interactions:
        if rec.severity.requires_action:
            action = f"[{rec.severity.value.upper()}] {rec.drug}/{rec.gene}: {rec.recommendation}"
            if rec.alternatives:
                action += f" Alternatives: {', '.join(rec.alternatives)}"
            prioritized_recs.append(action)

    result = {
        "total_drugs": len(drug_list),
        "drugs_with_interactions": sorted(drugs_with_ix),
        "interactions": all_interactions,
        "severity_summary": severity_info,
        "gene_overlap": gene_overlap,
        "high_risk_combinations": high_risk,
        "recommendations": prioritized_recs,
    }

    logger.info(
        "Polypharmacy analysis: %d drugs, %d with PGx interactions, %d gene overlaps",
        len(drug_list),
        len(drugs_with_ix),
        len(gene_overlap),
    )

    return result


def suggest_alternatives(
    drug: str,
    phenotype: str | MetabolizerPhenotype,
    alternatives_db: dict[str, dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Suggest alternative drugs based on pharmacogenomic phenotype.

    Queries the alternatives database for drugs in the same therapeutic class
    that are less affected by the patient's pharmacogenomic profile.

    Args:
        drug: Drug name (case-insensitive)
        phenotype: Metabolizer phenotype
        alternatives_db: Optional custom alternatives database.
            If None, uses built-in data.

    Returns:
        Dictionary with:
            - "drug": Original drug name
            - "therapeutic_class": Drug's therapeutic class
            - "alternatives": List of alternative drug options with details
            - "phenotype": Patient's phenotype
            - "rationale": Why alternatives are suggested
    """
    if alternatives_db is None:
        alternatives_db = _ALTERNATIVES_DB

    drug_lower = drug.lower().strip()
    if isinstance(phenotype, MetabolizerPhenotype):
        phenotype_str = phenotype.value
    else:
        phenotype_str = str(phenotype)

    drug_info = alternatives_db.get(drug_lower)
    if drug_info is None:
        return {
            "drug": drug,
            "therapeutic_class": "Unknown",
            "alternatives": [],
            "phenotype": phenotype_str,
            "rationale": f"No alternative drug information available for {drug}.",
        }

    alternatives_list: list[dict[str, Any]] = []
    for alt_name, alt_info in drug_info.get("alternatives", {}).items():
        alternatives_list.append(
            {
                "drug": alt_name,
                "notes": alt_info.get("notes", ""),
                "pgx_considerations": {k: v for k, v in alt_info.items() if k not in ("notes",)},
            }
        )

    result = {
        "drug": drug,
        "therapeutic_class": drug_info.get("therapeutic_class", "Unknown"),
        "alternatives": alternatives_list,
        "phenotype": phenotype_str,
        "rationale": f"Patient has {phenotype_str} status affecting {drug} metabolism. "
        f"Alternatives from the same therapeutic class ({drug_info.get('therapeutic_class', 'unknown')}) "
        f"with reduced or no pharmacogenomic interaction risk are listed.",
    }

    logger.info(
        "Suggested %d alternatives for %s (%s)",
        len(alternatives_list),
        drug,
        phenotype_str,
    )

    return result
