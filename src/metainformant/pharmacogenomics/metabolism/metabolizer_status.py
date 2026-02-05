"""Metabolizer status prediction from pharmacogenomic genotype data.

Provides CPIC-style metabolizer phenotype prediction from diplotype/genotype
information, activity score computation, and evidence-based dose adjustment
recommendations. Implements the Clinical Pharmacogenetics Implementation
Consortium (CPIC) activity score framework for CYP2D6, CYP2C19, and CYP2C9.

Activity Score (AS) Framework:
    Each allele is assigned a function value (0, 0.5, or 1.0). The diplotype
    activity score is the sum of both alleles. Metabolizer phenotype is then
    classified based on gene-specific thresholds.
"""

from __future__ import annotations

import math
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Built-in allele function tables (CPIC-derived)
# ---------------------------------------------------------------------------


def default_allele_function_table() -> dict:
    """Return the built-in allele function assignment table.

    Covers CYP2D6, CYP2C19, and CYP2C9 with CPIC-consensus function
    assignments. Each gene maps to a dict of allele -> function_value.

    Returns:
        Dictionary mapping gene name to a dict of
        ``{allele_name: function_value}``. Function values:
        0.0 = no function, 0.5 = decreased function, 1.0 = normal function,
        1.5 = increased function (CYP2C19 *17).
    """
    return {
        "CYP2D6": {
            "*1": 1.0,  # Normal function
            "*2": 1.0,  # Normal function
            "*3": 0.0,  # No function (frameshift)
            "*4": 0.0,  # No function (splicing defect)
            "*5": 0.0,  # No function (gene deletion)
            "*6": 0.0,  # No function (frameshift)
            "*7": 0.0,  # No function
            "*8": 0.0,  # No function
            "*9": 0.5,  # Decreased function
            "*10": 0.5,  # Decreased function (common in East Asian)
            "*14": 0.5,  # Decreased function
            "*17": 0.5,  # Decreased function
            "*29": 0.5,  # Decreased function
            "*41": 0.5,  # Decreased function (reduced splicing)
            "*1xN": 2.0,  # Increased function (gene duplication, N copies)
            "*2xN": 2.0,  # Increased function (gene duplication)
        },
        "CYP2C19": {
            "*1": 1.0,  # Normal function
            "*2": 0.0,  # No function (splicing defect, 681G>A)
            "*3": 0.0,  # No function (premature stop)
            "*4": 0.0,  # No function
            "*5": 0.0,  # No function
            "*6": 0.0,  # No function
            "*7": 0.0,  # No function
            "*8": 0.0,  # No function
            "*9": 0.5,  # Decreased function
            "*10": 0.5,  # Decreased function
            "*17": 1.5,  # Increased function (promoter variant)
        },
        "CYP2C9": {
            "*1": 1.0,  # Normal function
            "*2": 0.5,  # Decreased function (R144C)
            "*3": 0.0,  # No function (I359L, substantially reduced)
            "*4": 0.5,  # Decreased function
            "*5": 0.0,  # No function
            "*6": 0.0,  # No function
            "*8": 0.5,  # Decreased function
            "*11": 0.5,  # Decreased function
        },
        "CYP3A5": {
            "*1": 1.0,  # Normal function (expresser)
            "*3": 0.0,  # No function (non-expresser, most common)
            "*6": 0.0,  # No function
            "*7": 0.0,  # No function
        },
        "CYP1A2": {
            "*1A": 1.0,  # Normal function
            "*1C": 0.5,  # Decreased function
            "*1F": 1.5,  # Increased inducibility
            "*1K": 0.5,  # Decreased function
        },
        "DPYD": {
            "*1": 1.0,  # Normal function
            "*2A": 0.0,  # No function (IVS14+1G>A)
            "*13": 0.0,  # No function
            "c.2846A>T": 0.5,  # Decreased function
            "HapB3": 0.5,  # Decreased function
        },
        "TPMT": {
            "*1": 1.0,  # Normal function
            "*2": 0.0,  # No function
            "*3A": 0.0,  # No function
            "*3B": 0.0,  # No function
            "*3C": 0.0,  # No function
            "*4": 0.0,  # No function
        },
        "UGT1A1": {
            "*1": 1.0,  # Normal function
            "*6": 0.5,  # Decreased function
            "*28": 0.5,  # Decreased function (7 TA repeats)
            "*36": 1.0,  # Normal function (5 TA repeats)
            "*37": 0.0,  # No function (8 TA repeats)
        },
    }


# ---------------------------------------------------------------------------
# Metabolizer classification thresholds
# ---------------------------------------------------------------------------

_METABOLIZER_THRESHOLDS: dict[str, dict[str, tuple[float, float]]] = {
    "CYP2D6": {
        "poor": (0.0, 0.0),
        "intermediate": (0.25, 1.0),
        "normal": (1.25, 2.25),
        "ultrarapid": (2.5, float("inf")),
    },
    "CYP2C19": {
        "poor": (0.0, 0.0),
        "intermediate": (0.5, 1.0),
        "normal": (1.5, 2.0),
        "rapid": (2.5, 2.5),
        "ultrarapid": (3.0, float("inf")),
    },
    "CYP2C9": {
        "poor": (0.0, 0.0),
        "intermediate": (0.5, 1.0),
        "normal": (1.5, 2.0),
    },
    "DPYD": {
        "poor": (0.0, 0.0),
        "intermediate": (0.5, 1.0),
        "normal": (1.5, 2.0),
    },
    "TPMT": {
        "poor": (0.0, 0.0),
        "intermediate": (0.5, 1.0),
        "normal": (1.5, 2.0),
    },
    "UGT1A1": {
        "poor": (0.0, 0.5),
        "intermediate": (1.0, 1.0),
        "normal": (1.5, 2.0),
    },
}

# Fallback thresholds for genes not in the lookup
_DEFAULT_THRESHOLDS: dict[str, tuple[float, float]] = {
    "poor": (0.0, 0.0),
    "intermediate": (0.5, 1.0),
    "normal": (1.5, 2.0),
}


# ---------------------------------------------------------------------------
# Dose adjustment guidelines
# ---------------------------------------------------------------------------

_DOSE_GUIDELINES: dict[str, dict[str, dict[str, Any]]] = {
    "codeine": {
        "poor": {
            "recommendation": "Avoid codeine; use non-tramadol analgesic",
            "adjusted_dose_fraction": 0.0,
            "evidence_level": "A",
            "guideline_source": "CPIC",
        },
        "intermediate": {
            "recommendation": "Use codeine with caution at lowest effective dose; monitor for efficacy",
            "adjusted_dose_fraction": 0.75,
            "evidence_level": "A",
            "guideline_source": "CPIC",
        },
        "normal": {
            "recommendation": "Use standard dosing",
            "adjusted_dose_fraction": 1.0,
            "evidence_level": "A",
            "guideline_source": "CPIC",
        },
        "ultrarapid": {
            "recommendation": "Avoid codeine due to increased morphine formation and toxicity risk",
            "adjusted_dose_fraction": 0.0,
            "evidence_level": "A",
            "guideline_source": "CPIC",
        },
    },
    "tramadol": {
        "poor": {
            "recommendation": "Avoid tramadol; use non-CYP2D6-dependent analgesic",
            "adjusted_dose_fraction": 0.0,
            "evidence_level": "B",
            "guideline_source": "CPIC",
        },
        "intermediate": {
            "recommendation": "Use tramadol at standard dose; monitor for reduced efficacy",
            "adjusted_dose_fraction": 1.0,
            "evidence_level": "B",
            "guideline_source": "CPIC",
        },
        "normal": {
            "recommendation": "Use standard dosing",
            "adjusted_dose_fraction": 1.0,
            "evidence_level": "A",
            "guideline_source": "CPIC",
        },
        "ultrarapid": {
            "recommendation": "Avoid tramadol due to risk of respiratory depression",
            "adjusted_dose_fraction": 0.0,
            "evidence_level": "B",
            "guideline_source": "CPIC",
        },
    },
    "tamoxifen": {
        "poor": {
            "recommendation": "Consider alternative endocrine therapy (e.g., aromatase inhibitor)",
            "adjusted_dose_fraction": 0.0,
            "evidence_level": "A",
            "guideline_source": "CPIC",
        },
        "intermediate": {
            "recommendation": "Consider alternative or higher-dose tamoxifen (40 mg)",
            "adjusted_dose_fraction": 2.0,
            "evidence_level": "B",
            "guideline_source": "CPIC",
        },
        "normal": {
            "recommendation": "Use standard dosing (20 mg daily)",
            "adjusted_dose_fraction": 1.0,
            "evidence_level": "A",
            "guideline_source": "CPIC",
        },
        "ultrarapid": {
            "recommendation": "Use standard dosing; monitor for increased side effects",
            "adjusted_dose_fraction": 1.0,
            "evidence_level": "B",
            "guideline_source": "CPIC",
        },
    },
    "clopidogrel": {
        "poor": {
            "recommendation": "Use alternative antiplatelet (e.g., prasugrel, ticagrelor)",
            "adjusted_dose_fraction": 0.0,
            "evidence_level": "A",
            "guideline_source": "CPIC",
        },
        "intermediate": {
            "recommendation": "Use alternative antiplatelet if feasible",
            "adjusted_dose_fraction": 0.0,
            "evidence_level": "A",
            "guideline_source": "CPIC",
        },
        "normal": {
            "recommendation": "Use standard dosing",
            "adjusted_dose_fraction": 1.0,
            "evidence_level": "A",
            "guideline_source": "CPIC",
        },
        "rapid": {
            "recommendation": "Use standard dosing",
            "adjusted_dose_fraction": 1.0,
            "evidence_level": "A",
            "guideline_source": "CPIC",
        },
        "ultrarapid": {
            "recommendation": "Use standard dosing",
            "adjusted_dose_fraction": 1.0,
            "evidence_level": "A",
            "guideline_source": "CPIC",
        },
    },
    "warfarin": {
        "poor": {
            "recommendation": "Reduce initial warfarin dose by 50-80%; titrate based on INR",
            "adjusted_dose_fraction": 0.3,
            "evidence_level": "A",
            "guideline_source": "CPIC",
        },
        "intermediate": {
            "recommendation": "Reduce initial warfarin dose by 20-40%; titrate based on INR",
            "adjusted_dose_fraction": 0.7,
            "evidence_level": "A",
            "guideline_source": "CPIC",
        },
        "normal": {
            "recommendation": "Use standard dosing algorithm",
            "adjusted_dose_fraction": 1.0,
            "evidence_level": "A",
            "guideline_source": "CPIC",
        },
    },
    "fluorouracil": {
        "poor": {
            "recommendation": "Avoid fluoropyrimidines; use alternative chemotherapy",
            "adjusted_dose_fraction": 0.0,
            "evidence_level": "A",
            "guideline_source": "CPIC",
        },
        "intermediate": {
            "recommendation": "Reduce starting dose by 25-50%; titrate based on toxicity",
            "adjusted_dose_fraction": 0.5,
            "evidence_level": "A",
            "guideline_source": "CPIC",
        },
        "normal": {
            "recommendation": "Use standard dosing",
            "adjusted_dose_fraction": 1.0,
            "evidence_level": "A",
            "guideline_source": "CPIC",
        },
    },
    "azathioprine": {
        "poor": {
            "recommendation": "Consider alternative agent; if used, reduce dose to 10% and monitor",
            "adjusted_dose_fraction": 0.1,
            "evidence_level": "A",
            "guideline_source": "CPIC",
        },
        "intermediate": {
            "recommendation": "Start at 30-70% of standard dose; titrate based on tolerance",
            "adjusted_dose_fraction": 0.5,
            "evidence_level": "A",
            "guideline_source": "CPIC",
        },
        "normal": {
            "recommendation": "Use standard dosing",
            "adjusted_dose_fraction": 1.0,
            "evidence_level": "A",
            "guideline_source": "CPIC",
        },
    },
    "omeprazole": {
        "poor": {
            "recommendation": "Reduce omeprazole dose by 50%; consider alternative PPI",
            "adjusted_dose_fraction": 0.5,
            "evidence_level": "B",
            "guideline_source": "CPIC",
        },
        "intermediate": {
            "recommendation": "Standard dosing likely adequate; monitor efficacy",
            "adjusted_dose_fraction": 1.0,
            "evidence_level": "B",
            "guideline_source": "CPIC",
        },
        "normal": {
            "recommendation": "Use standard dosing",
            "adjusted_dose_fraction": 1.0,
            "evidence_level": "A",
            "guideline_source": "CPIC",
        },
        "rapid": {
            "recommendation": "Increase dose by 50-100% for H. pylori eradication",
            "adjusted_dose_fraction": 1.5,
            "evidence_level": "B",
            "guideline_source": "CPIC",
        },
        "ultrarapid": {
            "recommendation": "Increase dose by 100-200% or use alternative therapy",
            "adjusted_dose_fraction": 2.0,
            "evidence_level": "B",
            "guideline_source": "CPIC",
        },
    },
}


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------


def compute_activity_score(
    diplotype: str,
    allele_function_table: dict,
) -> float:
    """Compute CPIC-style activity score from a diplotype string.

    Parses a diplotype of the form ``"*1/*4"`` and sums the function values
    of each allele from the provided table.

    Args:
        diplotype: Diplotype string, e.g. ``"*1/*4"``, ``"*1xN/*2"``.
        allele_function_table: Dictionary mapping allele name to numeric
            function value (e.g. ``{"*1": 1.0, "*4": 0.0}``).

    Returns:
        Activity score (sum of both allele function values).

    Raises:
        ValueError: If diplotype format is invalid or alleles not found.
    """
    parts = diplotype.strip().split("/")
    if len(parts) != 2:
        raise ValueError(f"Invalid diplotype format '{diplotype}': expected exactly one '/' separator")

    allele_a = parts[0].strip()
    allele_b = parts[1].strip()

    score_a = allele_function_table.get(allele_a)
    score_b = allele_function_table.get(allele_b)

    if score_a is None:
        raise ValueError(
            f"Allele '{allele_a}' not found in function table. "
            f"Available alleles: {sorted(allele_function_table.keys())}"
        )
    if score_b is None:
        raise ValueError(
            f"Allele '{allele_b}' not found in function table. "
            f"Available alleles: {sorted(allele_function_table.keys())}"
        )

    activity_score = float(score_a) + float(score_b)
    logger.debug(
        "Activity score for %s: %s(%.1f) + %s(%.1f) = %.2f",
        diplotype,
        allele_a,
        score_a,
        allele_b,
        score_b,
        activity_score,
    )
    return activity_score


def classify_metabolizer(activity_score: float, gene: str) -> str:
    """Classify metabolizer phenotype from an activity score.

    Uses gene-specific thresholds from the CPIC framework to map an
    activity score to a phenotype category.

    Args:
        activity_score: Numeric activity score from :func:`compute_activity_score`.
        gene: Gene symbol (e.g. ``"CYP2D6"``).

    Returns:
        Phenotype string: one of ``"poor"``, ``"intermediate"``, ``"normal"``,
        ``"rapid"``, or ``"ultrarapid"``.
    """
    gene_upper = gene.strip().upper()
    thresholds = _METABOLIZER_THRESHOLDS.get(gene_upper, _DEFAULT_THRESHOLDS)

    for phenotype, (low, high) in thresholds.items():
        if low <= activity_score <= high:
            logger.debug(
                "Classified %s AS=%.2f as %s metabolizer",
                gene,
                activity_score,
                phenotype,
            )
            return phenotype

    # If above all ranges, assume ultrarapid; if between ranges, pick nearest
    max_threshold_phenotype = max(thresholds.items(), key=lambda x: x[1][1])[0]
    logger.warning(
        "Activity score %.2f for %s outside defined ranges; defaulting to %s",
        activity_score,
        gene,
        max_threshold_phenotype,
    )
    return max_threshold_phenotype


def predict_metabolizer_status(
    genotype: dict,
    gene: str,
    activity_scores: dict | None = None,
) -> dict:
    """Predict metabolizer phenotype from genotype information.

    Takes a genotype dictionary containing diplotype information and
    computes the metabolizer phenotype using the CPIC activity score
    framework.

    Args:
        genotype: Dictionary with at least a ``"diplotype"`` key containing
            the diplotype string (e.g. ``"*1/*4"``), or a ``"alleles"`` key
            with a list of two allele names.
        gene: Gene symbol (e.g. ``"CYP2D6"``).
        activity_scores: Optional custom allele function table for the gene.
            If ``None``, uses the built-in table.

    Returns:
        Dictionary with keys:
            - phenotype: Metabolizer classification string.
            - activity_score: Computed numeric activity score.
            - alleles: Tuple of two allele names.
            - gene: Gene symbol.
            - diplotype: Diplotype string.
            - thresholds_used: The threshold ranges applied.

    Raises:
        ValueError: If genotype is missing required fields.
    """
    gene_upper = gene.strip().upper()

    # Resolve allele function table
    if activity_scores is not None:
        allele_table = activity_scores
    else:
        full_table = default_allele_function_table()
        allele_table = full_table.get(gene_upper, {})
        if not allele_table:
            raise ValueError(
                f"No allele function table available for gene '{gene}'. "
                f"Available genes: {sorted(full_table.keys())}"
            )

    # Extract diplotype
    if "diplotype" in genotype:
        diplotype = genotype["diplotype"]
    elif "alleles" in genotype:
        alleles = genotype["alleles"]
        if len(alleles) != 2:
            raise ValueError(f"Expected exactly 2 alleles, got {len(alleles)}")
        diplotype = f"{alleles[0]}/{alleles[1]}"
    else:
        raise ValueError("Genotype must contain 'diplotype' or 'alleles' key")

    # Compute activity score and classify
    score = compute_activity_score(diplotype, allele_table)
    phenotype = classify_metabolizer(score, gene_upper)

    # Determine alleles
    allele_parts = diplotype.strip().split("/")
    alleles_tuple = (allele_parts[0].strip(), allele_parts[1].strip())

    thresholds = _METABOLIZER_THRESHOLDS.get(gene_upper, _DEFAULT_THRESHOLDS)

    logger.info(
        "Predicted %s metabolizer status for %s %s: %s (AS=%.2f)",
        gene,
        diplotype,
        alleles_tuple,
        phenotype,
        score,
    )

    return {
        "phenotype": phenotype,
        "activity_score": score,
        "alleles": alleles_tuple,
        "gene": gene_upper,
        "diplotype": diplotype,
        "thresholds_used": {k: {"low": v[0], "high": v[1]} for k, v in thresholds.items()},
    }


def dose_adjustment(
    metabolizer_status: str,
    drug: str,
    dose_guidelines: dict | None = None,
) -> dict:
    """Recommend dose adjustment based on metabolizer status and drug.

    Looks up the drug in the dose guideline table and returns the
    recommendation for the given metabolizer phenotype.

    Args:
        metabolizer_status: Metabolizer phenotype string (e.g. ``"poor"``,
            ``"intermediate"``, ``"normal"``, ``"ultrarapid"``).
        drug: Drug name to look up.
        dose_guidelines: Optional custom dose guidelines. If ``None``, uses
            the built-in CPIC-derived guidelines. Expected format:
            ``{drug: {phenotype: {recommendation, adjusted_dose_fraction, ...}}}``.

    Returns:
        Dictionary with keys:
            - drug: Normalised drug name.
            - metabolizer_status: Input phenotype.
            - recommendation: Clinical recommendation text.
            - adjusted_dose_fraction: Multiplier for standard dose (0.0 = avoid,
              1.0 = standard, 2.0 = double).
            - evidence_level: Evidence grade (``"A"`` through ``"D"``).
            - guideline_source: Source of the guideline.
    """
    if dose_guidelines is None:
        dose_guidelines = _DOSE_GUIDELINES

    drug_lower = drug.strip().lower()
    status_lower = metabolizer_status.strip().lower()

    drug_entry = dose_guidelines.get(drug_lower)
    if drug_entry is None:
        # Check case-insensitive
        for key, val in dose_guidelines.items():
            if key.lower() == drug_lower:
                drug_entry = val
                break

    if drug_entry is None:
        logger.warning("No dose guidelines available for drug '%s'", drug)
        return {
            "drug": drug_lower,
            "metabolizer_status": status_lower,
            "recommendation": (
                f"No specific dose guidelines available for {drug}. " "Use standard dosing and monitor clinically."
            ),
            "adjusted_dose_fraction": 1.0,
            "evidence_level": "D",
            "guideline_source": "None",
        }

    status_entry = drug_entry.get(status_lower)
    if status_entry is None:
        # Try to find closest match
        logger.warning(
            "No guideline for %s metabolizer status '%s'; available: %s",
            drug,
            metabolizer_status,
            sorted(drug_entry.keys()),
        )
        return {
            "drug": drug_lower,
            "metabolizer_status": status_lower,
            "recommendation": (
                f"No specific guideline for {metabolizer_status} metabolizer "
                f"status with {drug}. Use standard dosing with clinical monitoring."
            ),
            "adjusted_dose_fraction": 1.0,
            "evidence_level": "D",
            "guideline_source": "None",
        }

    logger.info(
        "Dose adjustment for %s (%s metabolizer): fraction=%.2f, evidence=%s",
        drug,
        metabolizer_status,
        status_entry["adjusted_dose_fraction"],
        status_entry["evidence_level"],
    )

    return {
        "drug": drug_lower,
        "metabolizer_status": status_lower,
        "recommendation": status_entry["recommendation"],
        "adjusted_dose_fraction": status_entry["adjusted_dose_fraction"],
        "evidence_level": status_entry["evidence_level"],
        "guideline_source": status_entry["guideline_source"],
    }
