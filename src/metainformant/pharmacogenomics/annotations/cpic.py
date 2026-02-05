"""CPIC (Clinical Pharmacogenetics Implementation Consortium) guideline support.

Provides parsing and querying of CPIC guideline data, including drug-gene pair
lookups, dosing recommendations based on metabolizer phenotype, actionable gene
lists, and allele definition table parsing.

CPIC levels:
    A: Strong genetic evidence + strong prescribing action recommended
    A/B: Strong evidence, moderate action
    B: Moderate genetic evidence + prescribing action recommended
    C: Moderate evidence, optional action
    D: Weak evidence or limited data
"""

from __future__ import annotations

import csv
import os
from pathlib import Path
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


# ── Built-in CPIC guideline data ──────────────────────────────────────────────
# Comprehensive drug-gene interaction data derived from CPIC guidelines.
# Each entry contains the CPIC level, affected gene, drug, phenotype-specific
# recommendations, and guideline reference.

_CPIC_GUIDELINES: list[dict[str, Any]] = [
    # CYP2D6 guidelines
    {
        "drug": "codeine",
        "gene": "CYP2D6",
        "cpic_level": "A",
        "guideline_url": "https://cpicpgx.org/guidelines/guideline-for-codeine-and-cyp2d6/",
        "recommendations": {
            "Poor Metabolizer": {
                "recommendation": "Avoid codeine use. Use alternative analgesics not metabolized by CYP2D6.",
                "classification": "Strong",
                "implications": "Greatly reduced morphine formation. Insufficient pain relief.",
            },
            "Intermediate Metabolizer": {
                "recommendation": "Use codeine per standard of care. Monitor for reduced efficacy.",
                "classification": "Moderate",
                "implications": "Reduced morphine formation. May have reduced analgesia.",
            },
            "Normal Metabolizer": {
                "recommendation": "Use codeine per standard of care.",
                "classification": "Strong",
                "implications": "Normal morphine formation following codeine administration.",
            },
            "Ultrarapid Metabolizer": {
                "recommendation": "Avoid codeine use. Use alternative analgesics not metabolized by CYP2D6.",
                "classification": "Strong",
                "implications": "Greatly increased morphine formation. Risk of toxicity and respiratory depression.",
            },
        },
    },
    {
        "drug": "tramadol",
        "gene": "CYP2D6",
        "cpic_level": "A",
        "guideline_url": "https://cpicpgx.org/guidelines/cpic-guideline-for-codeine-and-cyp2d6/",
        "recommendations": {
            "Poor Metabolizer": {
                "recommendation": "Avoid tramadol. Use alternative analgesics.",
                "classification": "Strong",
                "implications": "Reduced formation of active metabolite O-desmethyltramadol.",
            },
            "Intermediate Metabolizer": {
                "recommendation": "Use tramadol per standard of care. Be alert for reduced efficacy.",
                "classification": "Moderate",
                "implications": "Reduced formation of O-desmethyltramadol.",
            },
            "Normal Metabolizer": {
                "recommendation": "Use tramadol per standard of care.",
                "classification": "Strong",
                "implications": "Normal metabolism.",
            },
            "Ultrarapid Metabolizer": {
                "recommendation": "Avoid tramadol. Use alternative analgesics.",
                "classification": "Strong",
                "implications": "Increased formation of O-desmethyltramadol. Risk of toxicity.",
            },
        },
    },
    {
        "drug": "tamoxifen",
        "gene": "CYP2D6",
        "cpic_level": "A",
        "guideline_url": "https://cpicpgx.org/guidelines/cpic-guideline-for-tamoxifen-based-on-cyp2d6-genotype/",
        "recommendations": {
            "Poor Metabolizer": {
                "recommendation": "Avoid tamoxifen. Use aromatase inhibitor if postmenopausal.",
                "classification": "Strong",
                "implications": "Greatly reduced endoxifen concentrations. Reduced efficacy expected.",
            },
            "Intermediate Metabolizer": {
                "recommendation": "Consider alternative endocrine therapy or higher dose tamoxifen.",
                "classification": "Moderate",
                "implications": "Reduced endoxifen concentrations.",
            },
            "Normal Metabolizer": {
                "recommendation": "Use tamoxifen per standard of care.",
                "classification": "Strong",
                "implications": "Normal endoxifen formation.",
            },
            "Ultrarapid Metabolizer": {
                "recommendation": "Use tamoxifen per standard of care.",
                "classification": "Strong",
                "implications": "Normal to increased endoxifen formation.",
            },
        },
    },
    # CYP2C19 guidelines
    {
        "drug": "clopidogrel",
        "gene": "CYP2C19",
        "cpic_level": "A",
        "guideline_url": "https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/",
        "recommendations": {
            "Poor Metabolizer": {
                "recommendation": "Use alternative antiplatelet therapy (prasugrel or ticagrelor) if no contraindication.",
                "classification": "Strong",
                "implications": "Significantly reduced platelet inhibition. Increased risk of cardiovascular events.",
            },
            "Intermediate Metabolizer": {
                "recommendation": "Use alternative antiplatelet therapy (prasugrel or ticagrelor) if no contraindication.",
                "classification": "Moderate",
                "implications": "Reduced platelet inhibition.",
            },
            "Normal Metabolizer": {
                "recommendation": "Use clopidogrel per standard of care.",
                "classification": "Strong",
                "implications": "Normal clopidogrel metabolism and platelet inhibition.",
            },
            "Rapid Metabolizer": {
                "recommendation": "Use clopidogrel per standard of care.",
                "classification": "Strong",
                "implications": "Normal to increased active metabolite formation.",
            },
            "Ultrarapid Metabolizer": {
                "recommendation": "Use clopidogrel per standard of care.",
                "classification": "Strong",
                "implications": "Increased active metabolite formation.",
            },
        },
    },
    {
        "drug": "voriconazole",
        "gene": "CYP2C19",
        "cpic_level": "A",
        "guideline_url": "https://cpicpgx.org/guidelines/guideline-for-voriconazole-and-cyp2c19/",
        "recommendations": {
            "Poor Metabolizer": {
                "recommendation": "Use alternative antifungal agent not metabolized by CYP2C19.",
                "classification": "Strong",
                "implications": "Higher voriconazole exposure. Increased risk of adverse events.",
            },
            "Intermediate Metabolizer": {
                "recommendation": "Use voriconazole at standard dose. Monitor drug levels.",
                "classification": "Moderate",
                "implications": "Moderately higher voriconazole exposure.",
            },
            "Normal Metabolizer": {
                "recommendation": "Use voriconazole per standard of care.",
                "classification": "Strong",
                "implications": "Normal voriconazole metabolism.",
            },
            "Rapid Metabolizer": {
                "recommendation": "Use alternative antifungal or increase dose with therapeutic monitoring.",
                "classification": "Moderate",
                "implications": "Lower voriconazole exposure. Risk of subtherapeutic levels.",
            },
            "Ultrarapid Metabolizer": {
                "recommendation": "Use alternative antifungal agent not metabolized by CYP2C19.",
                "classification": "Strong",
                "implications": "Much lower voriconazole exposure. High risk of treatment failure.",
            },
        },
    },
    # CYP2C9 + VKORC1 guidelines (warfarin)
    {
        "drug": "warfarin",
        "gene": "CYP2C9",
        "cpic_level": "A",
        "guideline_url": "https://cpicpgx.org/guidelines/guideline-for-warfarin-and-cyp2c9-and-vkorc1/",
        "recommendations": {
            "Poor Metabolizer": {
                "recommendation": "Reduce warfarin dose significantly (typically 50-80% reduction). "
                "Consider pharmacogenomic dosing algorithm.",
                "classification": "Strong",
                "implications": "Greatly reduced warfarin clearance. High bleeding risk with standard dose.",
            },
            "Intermediate Metabolizer": {
                "recommendation": "Reduce warfarin dose (typically 20-40% reduction). "
                "Consider pharmacogenomic dosing algorithm.",
                "classification": "Moderate",
                "implications": "Reduced warfarin clearance.",
            },
            "Normal Metabolizer": {
                "recommendation": "Use warfarin per standard of care. Consider pharmacogenomic dosing algorithm.",
                "classification": "Strong",
                "implications": "Normal warfarin metabolism.",
            },
        },
    },
    # DPYD guidelines
    {
        "drug": "fluorouracil",
        "gene": "DPYD",
        "cpic_level": "A",
        "guideline_url": "https://cpicpgx.org/guidelines/guideline-for-fluoropyrimidines-and-dpyd/",
        "recommendations": {
            "Poor Metabolizer": {
                "recommendation": "Avoid fluorouracil, capecitabine, or tegafur. "
                "If unavoidable, use strongly reduced dose (< 25% of standard).",
                "classification": "Strong",
                "implications": "Complete DPD deficiency. Very high risk of severe, potentially fatal toxicity.",
            },
            "Intermediate Metabolizer": {
                "recommendation": "Reduce starting dose by 50%. Titrate based on toxicity and/or drug levels.",
                "classification": "Strong",
                "implications": "Partial DPD deficiency. Increased risk of severe toxicity.",
            },
            "Normal Metabolizer": {
                "recommendation": "Use fluoropyrimidines per standard of care.",
                "classification": "Strong",
                "implications": "Normal DPD activity.",
            },
        },
    },
    # TPMT / NUDT15 guidelines
    {
        "drug": "azathioprine",
        "gene": "TPMT",
        "cpic_level": "A",
        "guideline_url": "https://cpicpgx.org/guidelines/guideline-for-thiopurines-and-tpmt-and-nudt15/",
        "recommendations": {
            "Poor Metabolizer": {
                "recommendation": "Consider alternative agent or reduce dose to 10% of standard. "
                "Allow 4-6 weeks to reach steady state before dose adjustment.",
                "classification": "Strong",
                "implications": "Homozygous for non-functional TPMT alleles. Very high risk of myelosuppression.",
            },
            "Intermediate Metabolizer": {
                "recommendation": "Reduce starting dose to 30-80% of standard. " "Titrate based on tolerance.",
                "classification": "Strong",
                "implications": "Heterozygous. Moderate risk of myelosuppression.",
            },
            "Normal Metabolizer": {
                "recommendation": "Use azathioprine per standard of care. Start at full dose.",
                "classification": "Strong",
                "implications": "Normal TPMT activity. Low risk of myelosuppression.",
            },
        },
    },
    # SLCO1B1 guidelines
    {
        "drug": "simvastatin",
        "gene": "SLCO1B1",
        "cpic_level": "A",
        "guideline_url": "https://cpicpgx.org/guidelines/guideline-for-statins/",
        "recommendations": {
            "Poor Metabolizer": {
                "recommendation": "Use alternative statin (rosuvastatin, pravastatin, fluvastatin). "
                "If simvastatin required, use lower dose (max 20mg).",
                "classification": "Strong",
                "implications": "Increased simvastatin acid exposure. High risk of myopathy.",
            },
            "Intermediate Metabolizer": {
                "recommendation": "Use lower dose simvastatin (max 20mg) or consider alternative statin.",
                "classification": "Moderate",
                "implications": "Moderately increased simvastatin acid exposure.",
            },
            "Normal Metabolizer": {
                "recommendation": "Use simvastatin per standard of care.",
                "classification": "Strong",
                "implications": "Normal OATP1B1 transporter function.",
            },
        },
    },
]


def load_cpic_guidelines(
    filepath: str | Path | None = None,
) -> list[dict[str, Any]]:
    """Load CPIC guideline data.

    Loads from either an external file (TSV or JSON) or the built-in guideline
    table. The built-in table covers major CPIC Level A/A-B guidelines.

    Args:
        filepath: Optional path to external CPIC guideline file.
            Supports JSON and TSV formats. If None, uses built-in data.

    Returns:
        List of guideline dictionaries with keys: drug, gene, cpic_level,
        guideline_url, recommendations
    """
    if filepath is not None:
        return _load_guidelines_from_file(filepath)

    logger.info("Loaded %d built-in CPIC guidelines", len(_CPIC_GUIDELINES))
    return list(_CPIC_GUIDELINES)


def _load_guidelines_from_file(filepath: str | Path) -> list[dict[str, Any]]:
    """Load CPIC guidelines from an external file.

    Args:
        filepath: Path to the guideline file (JSON or TSV)

    Returns:
        List of parsed guideline dictionaries
    """
    from metainformant.core import io

    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"CPIC guideline file not found: {filepath}")

    if path.suffix in (".json", ".gz"):
        data = io.load_json(path)
        if isinstance(data, list):
            return data
        elif isinstance(data, dict):
            return data.get("guidelines", [data])
        return [data]

    # TSV format
    guidelines: list[dict[str, Any]] = []
    with open(path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            guidelines.append(
                {
                    "drug": row.get("drug", "").lower(),
                    "gene": row.get("gene", "").upper(),
                    "cpic_level": row.get("cpic_level", row.get("level", "")),
                    "guideline_url": row.get("guideline_url", row.get("url", "")),
                    "recommendations": {},
                }
            )

    logger.info("Loaded %d CPIC guidelines from %s", len(guidelines), filepath)
    return guidelines


def lookup_drug_gene(
    drug: str,
    gene: str,
    guidelines: list[dict[str, Any]] | None = None,
) -> dict[str, Any] | None:
    """Look up CPIC guideline for a specific drug-gene pair.

    Args:
        drug: Drug name (case-insensitive)
        gene: Gene symbol (case-insensitive)
        guidelines: Optional pre-loaded guidelines. If None, uses built-in data.

    Returns:
        Guideline dictionary if found, None otherwise
    """
    if guidelines is None:
        guidelines = _CPIC_GUIDELINES

    drug_lower = drug.lower().strip()
    gene_upper = gene.upper().strip()

    for guideline in guidelines:
        if guideline["drug"].lower() == drug_lower and guideline["gene"].upper() == gene_upper:
            logger.debug("Found CPIC guideline for %s/%s (Level %s)", drug, gene, guideline["cpic_level"])
            return guideline

    logger.debug("No CPIC guideline found for %s/%s", drug, gene)
    return None


def get_dosing_recommendation(
    drug: str,
    phenotype: str,
    guidelines: list[dict[str, Any]] | None = None,
) -> dict[str, Any] | None:
    """Get dosing recommendation for a drug based on metabolizer phenotype.

    Searches all guidelines for the specified drug and returns the
    recommendation matching the given phenotype. If the drug has guidelines
    for multiple genes, returns the first matching recommendation.

    Args:
        drug: Drug name (case-insensitive)
        phenotype: Metabolizer phenotype string (e.g., "Poor Metabolizer", "PM",
            "Normal Metabolizer", "NM")
        guidelines: Optional pre-loaded guidelines

    Returns:
        Dictionary with recommendation details, or None if not found.
        Keys: drug, gene, cpic_level, phenotype, recommendation, classification, implications
    """
    if guidelines is None:
        guidelines = _CPIC_GUIDELINES

    drug_lower = drug.lower().strip()
    phenotype_normalized = _normalize_phenotype_string(phenotype)

    for guideline in guidelines:
        if guideline["drug"].lower() != drug_lower:
            continue

        recommendations = guideline.get("recommendations", {})
        for pheno_key, rec_data in recommendations.items():
            if _normalize_phenotype_string(pheno_key) == phenotype_normalized:
                result = {
                    "drug": guideline["drug"],
                    "gene": guideline["gene"],
                    "cpic_level": guideline["cpic_level"],
                    "phenotype": pheno_key,
                    "recommendation": rec_data.get("recommendation", ""),
                    "classification": rec_data.get("classification", ""),
                    "implications": rec_data.get("implications", ""),
                    "guideline_url": guideline.get("guideline_url", ""),
                }
                logger.info(
                    "Dosing recommendation for %s (%s): %s",
                    drug,
                    pheno_key,
                    rec_data.get("recommendation", "")[:80],
                )
                return result

    logger.debug("No dosing recommendation found for %s with phenotype %s", drug, phenotype)
    return None


def _normalize_phenotype_string(phenotype: str) -> str:
    """Normalize phenotype string for matching.

    Handles abbreviations (PM, IM, NM, RM, UM) and full names.

    Args:
        phenotype: Raw phenotype string

    Returns:
        Normalized phenotype string
    """
    abbrev_map = {
        "pm": "poor metabolizer",
        "im": "intermediate metabolizer",
        "nm": "normal metabolizer",
        "em": "normal metabolizer",  # "Extensive" is old terminology
        "rm": "rapid metabolizer",
        "um": "ultrarapid metabolizer",
    }

    normalized = phenotype.strip().lower()
    return abbrev_map.get(normalized, normalized)


def list_actionable_genes(
    guidelines: list[dict[str, Any]] | None = None,
    min_level: str = "B",
) -> list[dict[str, str]]:
    """List CPIC Level A and B (actionable) gene-drug pairs.

    Args:
        guidelines: Optional pre-loaded guidelines
        min_level: Minimum CPIC level to include ("A" for Level A only,
            "B" for Level A and B)

    Returns:
        List of dictionaries with keys: gene, drug, cpic_level
    """
    if guidelines is None:
        guidelines = _CPIC_GUIDELINES

    accepted_levels = {"A"}
    if min_level.upper() in ("B", "C", "D"):
        accepted_levels.add("A/B")
        accepted_levels.add("B")
    if min_level.upper() in ("C", "D"):
        accepted_levels.add("C")
    if min_level.upper() == "D":
        accepted_levels.add("D")

    actionable: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()

    for guideline in guidelines:
        level = guideline.get("cpic_level", "")
        if level not in accepted_levels:
            continue

        key = (guideline["gene"].upper(), guideline["drug"].lower())
        if key in seen:
            continue
        seen.add(key)

        actionable.append(
            {
                "gene": guideline["gene"].upper(),
                "drug": guideline["drug"],
                "cpic_level": level,
            }
        )

    # Sort by gene, then drug
    actionable.sort(key=lambda x: (x["gene"], x["drug"]))

    logger.info(
        "Found %d actionable gene-drug pairs at CPIC Level %s or higher",
        len(actionable),
        min_level,
    )

    return actionable


def parse_cpic_allele_definitions(
    filepath: str | Path,
) -> dict[str, list[dict[str, Any]]]:
    """Parse CPIC allele definition tables.

    Reads CPIC-formatted allele definition files (TSV format) that define the
    genetic variants comprising each star allele for pharmacogenes.

    Expected TSV columns: gene, allele, rsid, position, ref, alt, function, activity_value

    Args:
        filepath: Path to allele definition TSV/JSON file

    Returns:
        Dictionary mapping gene -> list of allele definition dicts.
        Each allele dict contains: allele, defining_variants, function, activity_value

    Raises:
        FileNotFoundError: If file does not exist
    """
    from metainformant.core import io

    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"Allele definition file not found: {filepath}")

    if path.suffix in (".json", ".gz"):
        data = io.load_json(path)
        if isinstance(data, dict):
            return data
        raise ValueError(f"Expected JSON dict, got {type(data)}")

    # TSV parsing
    result: dict[str, list[dict[str, Any]]] = {}

    with open(path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        allele_variants: dict[tuple[str, str], list[str]] = {}
        allele_meta: dict[tuple[str, str], dict[str, Any]] = {}

        for row in reader:
            gene = row.get("gene", "").upper()
            allele = row.get("allele", "")
            rsid = row.get("rsid", "")
            position = row.get("position", "")
            function = row.get("function", "Unknown function")
            activity = row.get("activity_value", "1.0")

            key = (gene, allele)
            if key not in allele_variants:
                allele_variants[key] = []
                allele_meta[key] = {
                    "allele": allele,
                    "function": function,
                    "activity_value": float(activity) if activity else 1.0,
                }

            # Add variant identifier
            if rsid:
                allele_variants[key].append(rsid)
            elif position:
                ref = row.get("ref", "")
                alt = row.get("alt", "")
                allele_variants[key].append(f"{position}:{ref}:{alt}")

    # Assemble results
    for (gene, allele), variants in allele_variants.items():
        if gene not in result:
            result[gene] = []

        meta = allele_meta[(gene, allele)]
        result[gene].append(
            {
                "allele": allele,
                "defining_variants": variants,
                "function": meta["function"],
                "activity_value": meta["activity_value"],
            }
        )

    total_alleles = sum(len(v) for v in result.values())
    logger.info(
        "Parsed %d allele definitions across %d genes from %s",
        total_alleles,
        len(result),
        filepath,
    )

    return result
