"""Metabolizer status prediction from pharmacogenomic genotype data.

Provides CPIC-style metabolizer phenotype prediction from diplotype/genotype
information, activity score computation, and evidence-based dose adjustment
recommendations. Implements the Clinical Pharmacogenetics Implementation
Consortium (CPIC) activity score framework for CYP2D6, CYP2C19, and CYP2C9.

This module is a string/dict front end over the ``alleles`` subpackage and holds
no allele or phenotype data of its own: allele function values come from
``alleles.diplotype`` and phenotype assignment from ``alleles.phenotype``, so
``predict_metabolizer_status`` and ``alleles.phenotype.classify_phenotype`` agree
for every diplotype. That includes the cases without a normal answer: a gene CPIC
defines no phenotype for (CYP1A2 has activity values but no phenotype) raises
``ValueError`` naming the supported genes, and a diplotype with an allele that has no
activity value is ``"indeterminate"`` with an activity score of ``None``, never
normal.

Activity Score (AS) Framework:
    Each allele is assigned a function value (0, 0.5, or 1.0). The diplotype
    activity score is the sum of both alleles. Metabolizer phenotype is then
    classified based on gene-specific thresholds.
"""

from __future__ import annotations

from typing import Any

from metainformant.core.utils.logging import get_logger

from ..alleles.diplotype import _ACTIVITY_SCORE_TABLES
from ..alleles.phenotype import (
    MetabolizerPhenotype,
    _phenotype_for_alleles,
    get_phenotype_thresholds,
    predict_metabolizer_status as _phenotype_from_score,
    supported_phenotype_genes,
)

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Built-in allele function tables (CPIC-derived)
# ---------------------------------------------------------------------------


def default_allele_function_table() -> dict:
    """Return the built-in allele function assignment table.

    A copy of the package's single allele table,
    ``alleles.diplotype._ACTIVITY_SCORE_TABLES`` (CPIC allele functionality tables), so
    the values here can never disagree with ``determine_diplotype`` or
    ``alleles.phenotype.classify_phenotype``. Change allele values there, not here.

    Returns:
        Dictionary mapping gene name to a dict of
        ``{allele_name: function_value}``. Function values:
        0.0 = no function, 0.5 = decreased function, 1.0 = normal function,
        1.5 = increased function (CYP2C19 *17); CYP2D6 *10 is 0.25. Covers
        CYP2D6, CYP2C19, CYP2C9, CYP3A5, CYP1A2, DPYD, TPMT, NUDT15, SLCO1B1
        and UGT1A1. CYP1A2's values are not from CPIC (see the table in
        ``alleles.diplotype``) and CPIC defines no CYP1A2 phenotype, so it can be
        scored but not classified. The copy is independent: changing it does not
        change the built-in table.
    """
    return {gene: dict(table) for gene, table in _ACTIVITY_SCORE_TABLES.items()}


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


def _split_diplotype(diplotype: str) -> tuple[str, str]:
    """Split a ``"*1/*4"`` diplotype string into its two allele names.

    Raises:
        ValueError: If the string does not contain exactly one ``/``.
    """
    parts = diplotype.strip().split("/")
    if len(parts) != 2:
        raise ValueError(f"Invalid diplotype format '{diplotype}': expected exactly one '/' separator")
    return parts[0].strip(), parts[1].strip()


def compute_activity_score(
    diplotype: str,
    allele_function_table: dict,
) -> float:
    """Compute CPIC-style activity score from a diplotype string.

    Parses a diplotype of the form ``"*1/*4"`` and sums the function values
    of each allele from the provided table.

    This is the strict score: an allele with no value in the table is an error, not
    normal function. :func:`predict_metabolizer_status` reports such a diplotype as
    indeterminate instead.

    Args:
        diplotype: Diplotype string, e.g. ``"*1/*4"``, ``"*1xN/*2"``.
        allele_function_table: Dictionary mapping allele name to numeric
            function value (e.g. ``{"*1": 1.0, "*4": 0.0}``).

    Returns:
        Activity score (sum of both allele function values).

    Raises:
        ValueError: If diplotype format is invalid or alleles not found.
    """
    allele_a, allele_b = _split_diplotype(diplotype)

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


def _phenotype_token(phenotype: MetabolizerPhenotype) -> str:
    """Return the lowercase token for a phenotype (``"poor"``, ``"ultrarapid"``, ...)."""
    return phenotype.value.split()[0].lower()


def _status_token(status: str) -> str:
    """Return the lowercase phenotype token for any accepted phenotype label.

    ``"Poor Metabolizer"``, ``"PM"``, ``"poor"`` and ``MetabolizerPhenotype.POOR`` all give
    ``"poor"``. A label that is not a known phenotype is returned stripped and lowercased.
    """
    try:
        return _phenotype_token(MetabolizerPhenotype.from_string(status))
    except ValueError:
        return status.strip().lower()


def classify_metabolizer(activity_score: float, gene: str) -> str:
    """Classify metabolizer phenotype from an activity score.

    Delegates to :func:`metainformant.pharmacogenomics.alleles.phenotype.predict_metabolizer_status`,
    the package's one score-to-phenotype mapping (CPIC definitions, cited there), and returns
    the phenotype as a lowercase token.

    Args:
        activity_score: Numeric activity score from :func:`compute_activity_score`.
        gene: Gene symbol (e.g. ``"CYP2D6"``).

    Returns:
        Phenotype string: one of ``"poor"``, ``"intermediate"``, ``"normal"``,
        ``"rapid"``, or ``"ultrarapid"``.

    Raises:
        ValueError: If the activity score is negative or NaN, or the gene has no CPIC
            phenotype definition (the message names the gene and lists the supported
            genes; CYP1A2 can be scored but not classified). No other gene's bands stand in.
    """
    gene_upper = gene.strip().upper()
    phenotype = _phenotype_from_score(activity_score, gene_upper)
    logger.debug("Classified %s AS=%.2f as %s metabolizer", gene, activity_score, phenotype.value)
    return _phenotype_token(phenotype)


def predict_metabolizer_status(
    genotype: dict,
    gene: str,
    activity_scores: dict | None = None,
) -> dict:
    """Predict metabolizer phenotype from genotype information.

    Takes a genotype dictionary containing diplotype information and
    computes the metabolizer phenotype using the CPIC activity score
    framework. The result is the same phenotype
    :func:`metainformant.pharmacogenomics.alleles.phenotype.classify_phenotype`
    gives for that diplotype (as a lowercase token), including for an allele
    missing from the table: the diplotype is ``"indeterminate"`` with an activity
    score of ``None`` and the allele named in ``unknown_alleles`` (CPIC's term for a
    diplotype with an uncharacterized allele), never normal function.

    Args:
        genotype: Dictionary with at least a ``"diplotype"`` key containing
            the diplotype string (e.g. ``"*1/*4"``), or a ``"alleles"`` key
            with a list of two allele names.
        gene: Gene symbol (e.g. ``"CYP2D6"``).
        activity_scores: Optional custom allele function table for the gene.
            If ``None``, uses the built-in table.

    Returns:
        Dictionary with keys:
            - phenotype: Metabolizer classification string (``"indeterminate"`` when an
              allele has no value in the table).
            - activity_score: Computed numeric activity score, or ``None`` when an allele
              has no value in the table.
            - alleles: Tuple of two allele names.
            - gene: Gene symbol.
            - diplotype: Diplotype string.
            - unknown_alleles: The alleles with no value in the table (empty when none).
            - thresholds_used: The activity-score bands the gene uses, keyed by
              phenotype token: ``{"low", "high", "low_inclusive", "high_inclusive"}``.
              CYP2C19, SLCO1B1 and UGT1A1 are classified from allele function, which
              these bands reproduce for the built-in alleles.

    Raises:
        ValueError: If genotype is missing required fields or the diplotype string is
            malformed, or if the gene has no CPIC phenotype definition (the message names
            the gene and lists the supported genes; CYP1A2 has activity values but no CPIC
            phenotype, so it can be scored with :func:`compute_activity_score` but not
            classified here).
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
                f"Supported genes: {', '.join(supported_phenotype_genes())}"
            )

    # Only a gene CPIC defines a phenotype for is classified: this raises for the rest (CYP1A2
    # has allele values, a custom table may name any gene) rather than borrowing another gene's bands
    thresholds = get_phenotype_thresholds(gene_upper)

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

    # Determine alleles (this also validates the diplotype format)
    alleles_tuple = _split_diplotype(diplotype)
    unknown_alleles = [allele for allele in dict.fromkeys(alleles_tuple) if allele not in allele_table]

    score: float | None
    if unknown_alleles:
        # An allele with no value is uncharacterized: the diplotype is indeterminate, as CPIC says,
        # and is never scored as normal function
        logger.warning(
            "No activity value for allele(s) %s of %s: metabolizer status is indeterminate",
            ", ".join(unknown_alleles),
            gene_upper,
        )
        score = None
        phenotype = _phenotype_token(MetabolizerPhenotype.INDETERMINATE)
    else:
        score = compute_activity_score(diplotype, allele_table)
        # Classify from the two alleles' own values, exactly as alleles.phenotype does
        phenotype = _phenotype_token(
            _phenotype_for_alleles(
                gene_upper,
                float(allele_table[alleles_tuple[0]]),
                float(allele_table[alleles_tuple[1]]),
            )
        )

    logger.info(
        "Predicted %s metabolizer status for %s %s: %s (AS=%s)",
        gene,
        diplotype,
        alleles_tuple,
        phenotype,
        "none" if score is None else f"{score:.2f}",
    )

    return {
        "phenotype": phenotype,
        "activity_score": score,
        "alleles": alleles_tuple,
        "gene": gene_upper,
        "diplotype": diplotype,
        "unknown_alleles": unknown_alleles,
        "thresholds_used": {
            _phenotype_token(band["phenotype"]): {
                "low": band["lower_bound"],
                "high": band["upper_bound"],
                "low_inclusive": band["lower_inclusive"],
                "high_inclusive": band["upper_inclusive"],
            }
            for band in thresholds
        },
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
        metabolizer_status: Metabolizer phenotype in any common form: a token
            (``"poor"``, ``"intermediate"``, ``"normal"``, ``"rapid"``,
            ``"ultrarapid"``), a CPIC name (``"Poor Metabolizer"``), an abbreviation
            (``"PM"``) or a ``MetabolizerPhenotype``, case-insensitively.
        drug: Drug name to look up.
        dose_guidelines: Optional custom dose guidelines. If ``None``, uses
            the built-in CPIC-derived guidelines. Expected format:
            ``{drug: {phenotype: {recommendation, adjusted_dose_fraction, ...}}}``;
            the phenotype keys may be written in any of the forms above.

    Returns:
        Dictionary with keys:
            - drug: Normalised drug name.
            - metabolizer_status: The phenotype token the lookup used (the stripped,
              lowercased input when it is not a known phenotype).
            - recommendation: Clinical recommendation text.
            - adjusted_dose_fraction: Multiplier for standard dose (0.0 = avoid,
              1.0 = standard, 2.0 = double); ``None`` for an ``"indeterminate"`` status
              (see below).
            - evidence_level: Evidence grade (``"A"`` through ``"D"``).
            - guideline_source: Source of the guideline.

    An ``"indeterminate"`` status, which :func:`predict_metabolizer_status` returns for a
    diplotype with an uncharacterized allele, is not normal: unless the guidelines hold an
    entry for it, no dose is recommended (``adjusted_dose_fraction`` is ``None``, not 1.0).
    An unknown drug or an unrecognised status label still returns standard dosing with
    evidence level ``"D"``.
    """
    if dose_guidelines is None:
        dose_guidelines = _DOSE_GUIDELINES

    drug_lower = drug.strip().lower()
    status_lower = _status_token(metabolizer_status)

    drug_entry = dose_guidelines.get(drug_lower)
    if drug_entry is None:
        # Check case-insensitive
        for key, val in dose_guidelines.items():
            if key.lower() == drug_lower:
                drug_entry = val
                break

    status_entry = None
    if drug_entry is not None:
        status_entry = drug_entry.get(status_lower)
        if status_entry is None:
            # Custom guidelines may key phenotypes by full name or abbreviation
            status_entry = next((v for k, v in drug_entry.items() if _status_token(k) == status_lower), None)

    if status_entry is None and status_lower == _phenotype_token(MetabolizerPhenotype.INDETERMINATE):
        logger.warning("No dose recommendation for %s: the metabolizer status is indeterminate", drug)
        return {
            "drug": drug_lower,
            "metabolizer_status": status_lower,
            "recommendation": (
                f"No dose recommendation for {drug}: the metabolizer phenotype is indeterminate (for example a "
                "diplotype with an uncharacterized allele). It is not assumed to be normal; resolve the "
                "phenotype before dosing by genotype."
            ),
            "adjusted_dose_fraction": None,
            "evidence_level": "D",
            "guideline_source": "None",
        }

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

    if status_entry is None:
        # Try to find closest match
        logger.warning(
            "No guideline for %s metabolizer status '%s'; available: %s",
            drug,
            status_lower,
            sorted(drug_entry.keys()),
        )
        return {
            "drug": drug_lower,
            "metabolizer_status": status_lower,
            "recommendation": (
                f"No specific guideline for {status_lower} metabolizer "
                f"status with {drug}. Use standard dosing with clinical monitoring."
            ),
            "adjusted_dose_fraction": 1.0,
            "evidence_level": "D",
            "guideline_source": "None",
        }

    logger.info(
        "Dose adjustment for %s (%s metabolizer): fraction=%.2f, evidence=%s",
        drug,
        status_lower,
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
