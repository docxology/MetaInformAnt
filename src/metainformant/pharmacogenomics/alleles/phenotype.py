"""Metabolizer phenotype prediction from pharmacogene diplotypes.

Maps activity scores to standardized metabolizer phenotype categories per CPIC
guidelines. Supports gene-specific threshold definitions and population-level
phenotype frequency estimation from allele frequency data.

Metabolizer categories:
    - Poor Metabolizer (PM): No or very low enzyme activity
    - Intermediate Metabolizer (IM): Reduced activity between PM and NM
    - Normal Metabolizer (NM): Typical activity level (reference)
    - Rapid Metabolizer (RM): Increased activity (one increased-function allele)
    - Ultrarapid Metabolizer (UM): Very high activity (gene duplication or two increased-function alleles)
"""

from __future__ import annotations

from enum import Enum
from itertools import combinations_with_replacement
from typing import Any

from metainformant.core.utils.logging import get_logger

from .diplotype import _ACTIVITY_SCORE_TABLES, Diplotype, determine_diplotype

logger = get_logger(__name__)


class MetabolizerPhenotype(str, Enum):
    """Standardized metabolizer phenotype categories.

    Values follow CPIC terminology for consistent clinical communication.
    Inherits from str for JSON serialization compatibility.
    """

    POOR = "Poor Metabolizer"
    INTERMEDIATE = "Intermediate Metabolizer"
    NORMAL = "Normal Metabolizer"
    RAPID = "Rapid Metabolizer"
    ULTRARAPID = "Ultrarapid Metabolizer"
    INDETERMINATE = "Indeterminate"

    @property
    def abbreviation(self) -> str:
        """Return standard abbreviation (PM, IM, NM, RM, UM)."""
        _abbrevs = {
            "Poor Metabolizer": "PM",
            "Intermediate Metabolizer": "IM",
            "Normal Metabolizer": "NM",
            "Rapid Metabolizer": "RM",
            "Ultrarapid Metabolizer": "UM",
            "Indeterminate": "IND",
        }
        return _abbrevs.get(self.value, "?")


# ── Gene-specific phenotype threshold definitions ──────────────────────────────
# These define how activity scores map to phenotype categories for each gene.
# Thresholds are [lower_bound_inclusive, upper_bound_exclusive) except for the
# last category which is upper-inclusive.
# Derived from CPIC gene-specific guideline tables.

_PHENOTYPE_THRESHOLDS: dict[str, list[tuple[float, float, MetabolizerPhenotype]]] = {
    "CYP2D6": [
        (0.0, 0.0, MetabolizerPhenotype.POOR),  # AS = 0
        (0.25, 1.0, MetabolizerPhenotype.INTERMEDIATE),  # 0 < AS < 1.25
        (1.0, 2.25, MetabolizerPhenotype.NORMAL),  # 1.25 <= AS <= 2.25
        (2.25, 99.0, MetabolizerPhenotype.ULTRARAPID),  # AS > 2.25
    ],
    "CYP2C19": [
        (0.0, 0.0, MetabolizerPhenotype.POOR),  # AS = 0
        (0.5, 1.0, MetabolizerPhenotype.INTERMEDIATE),  # 0 < AS < 1
        (1.0, 2.0, MetabolizerPhenotype.NORMAL),  # AS = 1 to < 2
        (2.0, 2.5, MetabolizerPhenotype.RAPID),  # AS = 2 (one *17)
        (2.5, 99.0, MetabolizerPhenotype.ULTRARAPID),  # AS >= 2.5 (two *17)
    ],
    "CYP2C9": [
        (0.0, 0.0, MetabolizerPhenotype.POOR),
        (0.5, 1.5, MetabolizerPhenotype.INTERMEDIATE),
        (1.5, 99.0, MetabolizerPhenotype.NORMAL),
    ],
    "CYP3A5": [
        (0.0, 0.0, MetabolizerPhenotype.POOR),
        (0.5, 1.0, MetabolizerPhenotype.INTERMEDIATE),
        (1.0, 99.0, MetabolizerPhenotype.NORMAL),
    ],
    "DPYD": [
        (0.0, 0.0, MetabolizerPhenotype.POOR),
        (0.5, 1.5, MetabolizerPhenotype.INTERMEDIATE),
        (1.5, 2.5, MetabolizerPhenotype.NORMAL),
    ],
    "TPMT": [
        (0.0, 0.0, MetabolizerPhenotype.POOR),
        (0.5, 1.5, MetabolizerPhenotype.INTERMEDIATE),
        (1.5, 99.0, MetabolizerPhenotype.NORMAL),
    ],
    "NUDT15": [
        (0.0, 0.0, MetabolizerPhenotype.POOR),
        (0.5, 1.5, MetabolizerPhenotype.INTERMEDIATE),
        (1.5, 99.0, MetabolizerPhenotype.NORMAL),
    ],
    "SLCO1B1": [
        (0.0, 0.5, MetabolizerPhenotype.POOR),
        (0.5, 1.5, MetabolizerPhenotype.INTERMEDIATE),
        (1.5, 99.0, MetabolizerPhenotype.NORMAL),
    ],
    "UGT1A1": [
        (0.0, 0.5, MetabolizerPhenotype.POOR),
        (0.5, 1.5, MetabolizerPhenotype.INTERMEDIATE),
        (1.5, 99.0, MetabolizerPhenotype.NORMAL),
    ],
}


def predict_metabolizer_status(
    activity_score: float,
    gene: str,
) -> MetabolizerPhenotype:
    """Map an activity score to a metabolizer phenotype category.

    Uses gene-specific threshold tables derived from CPIC guidelines. When a
    gene has no defined thresholds, applies a default mapping:
        AS = 0 -> PM, 0 < AS < 1.25 -> IM, 1.25 <= AS <= 2.25 -> NM, AS > 2.25 -> UM

    Args:
        activity_score: Combined activity score (typically sum of two allele scores)
        gene: Gene symbol (e.g., "CYP2D6")

    Returns:
        MetabolizerPhenotype enum value
    """
    gene_upper = gene.upper()
    thresholds = _PHENOTYPE_THRESHOLDS.get(gene_upper)

    if thresholds is None:
        logger.warning(
            "No phenotype thresholds defined for %s, using default CYP2D6-like mapping",
            gene_upper,
        )
        thresholds = _PHENOTYPE_THRESHOLDS["CYP2D6"]

    # Handle exact zero
    if activity_score == 0.0:
        for low, high, phenotype in thresholds:
            if low == 0.0 and high == 0.0:
                return phenotype
        return MetabolizerPhenotype.POOR

    # Find matching threshold range
    for low, high, phenotype in thresholds:
        if low == 0.0 and high == 0.0:
            continue  # Skip exact-zero entry
        if low <= activity_score < high:
            return phenotype

    # If score exceeds all thresholds, return the last (highest) category
    if thresholds:
        return thresholds[-1][2]

    return MetabolizerPhenotype.INDETERMINATE


def get_phenotype_thresholds(
    gene: str,
) -> list[dict[str, Any]]:
    """Get gene-specific phenotype threshold definitions.

    Args:
        gene: Gene symbol (e.g., "CYP2D6")

    Returns:
        List of threshold dictionaries with keys:
            - "lower_bound": Lower activity score boundary (inclusive)
            - "upper_bound": Upper activity score boundary (exclusive)
            - "phenotype": MetabolizerPhenotype value
            - "abbreviation": Phenotype abbreviation (PM, IM, NM, RM, UM)
    """
    gene_upper = gene.upper()
    thresholds = _PHENOTYPE_THRESHOLDS.get(gene_upper)

    if thresholds is None:
        return [
            {
                "lower_bound": 0.0,
                "upper_bound": 0.0,
                "phenotype": MetabolizerPhenotype.POOR,
                "abbreviation": "PM",
                "note": f"No specific thresholds defined for {gene_upper}, using defaults",
            }
        ]

    result: list[dict[str, Any]] = []
    for low, high, phenotype in thresholds:
        result.append(
            {
                "lower_bound": low,
                "upper_bound": high,
                "phenotype": phenotype,
                "abbreviation": phenotype.abbreviation,
            }
        )

    return result


def classify_phenotype(
    diplotype: Diplotype | str,
    gene: str,
) -> dict[str, Any]:
    """Full diplotype-to-phenotype classification pipeline.

    Takes a diplotype (object or string like "*1/*4") and gene, computes the
    activity score, and maps it to a metabolizer phenotype.

    Args:
        diplotype: Diplotype object or string in format "*X/*Y"
        gene: Gene symbol

    Returns:
        Dictionary with:
            - "gene": Gene symbol
            - "diplotype": Diplotype string
            - "activity_score": Computed activity score
            - "phenotype": MetabolizerPhenotype value
            - "phenotype_abbreviation": Short form (PM, IM, NM, RM, UM)
            - "clinical_significance": Brief clinical note
    """
    gene_upper = gene.upper()

    if isinstance(diplotype, str):
        parts = diplotype.strip().split("/")
        if len(parts) != 2:
            raise ValueError(f"Invalid diplotype string format: '{diplotype}'. Expected '*X/*Y'.")
        allele1 = parts[0].strip()
        allele2 = parts[1].strip()
        dip = determine_diplotype(allele1, allele2, gene_upper)
    else:
        dip = diplotype

    phenotype = predict_metabolizer_status(dip.activity_score, gene_upper)

    # Generate clinical significance note
    clinical_notes = _generate_clinical_note(phenotype, gene_upper)

    result = {
        "gene": gene_upper,
        "diplotype": dip.diplotype_string,
        "activity_score": dip.activity_score,
        "phenotype": phenotype,
        "phenotype_abbreviation": phenotype.abbreviation,
        "clinical_significance": clinical_notes,
    }

    logger.info(
        "Classified %s %s -> AS=%.2f -> %s (%s)",
        gene_upper,
        dip.diplotype_string,
        dip.activity_score,
        phenotype.value,
        phenotype.abbreviation,
    )

    return result


def _generate_clinical_note(phenotype: MetabolizerPhenotype, gene: str) -> str:
    """Generate a brief clinical significance note for a phenotype.

    Args:
        phenotype: Metabolizer phenotype
        gene: Gene symbol

    Returns:
        Clinical significance string
    """
    notes = {
        MetabolizerPhenotype.POOR: (
            f"Poor metabolizer for {gene}. Significantly reduced or absent enzyme activity. "
            "Consider alternative therapy or significant dose reduction per CPIC guidelines."
        ),
        MetabolizerPhenotype.INTERMEDIATE: (
            f"Intermediate metabolizer for {gene}. Reduced enzyme activity. "
            "May require dose adjustment depending on the specific drug."
        ),
        MetabolizerPhenotype.NORMAL: (
            f"Normal metabolizer for {gene}. Expected typical enzyme activity. "
            "Standard dosing recommendations apply."
        ),
        MetabolizerPhenotype.RAPID: (
            f"Rapid metabolizer for {gene}. Increased enzyme activity. "
            "May require dose increase for prodrugs or dose reduction for active drugs."
        ),
        MetabolizerPhenotype.ULTRARAPID: (
            f"Ultrarapid metabolizer for {gene}. Markedly increased enzyme activity. "
            "Consider significant dose adjustment or alternative therapy per CPIC guidelines."
        ),
        MetabolizerPhenotype.INDETERMINATE: (
            f"Indeterminate metabolizer status for {gene}. " "Insufficient data to assign a phenotype category."
        ),
    }
    return notes.get(phenotype, f"Unknown phenotype status for {gene}.")


def population_phenotype_frequencies(
    allele_frequencies: dict[str, float],
    gene: str,
) -> dict[str, float]:
    """Estimate expected phenotype distribution from allele frequencies.

    Calculates the expected frequency of each metabolizer phenotype in a
    population given the allele frequency spectrum, assuming Hardy-Weinberg
    equilibrium. For each possible diplotype (all combinations of alleles),
    the diplotype frequency is 2*p*q for heterozygotes and p^2 for homozygotes.

    Args:
        allele_frequencies: Mapping of allele name -> population frequency.
            Frequencies should sum to approximately 1.0.
        gene: Gene symbol

    Returns:
        Dictionary mapping phenotype abbreviation -> expected frequency.
        Example: {"PM": 0.05, "IM": 0.20, "NM": 0.65, "RM": 0.08, "UM": 0.02}
    """
    gene_upper = gene.upper()

    # Normalize allele names
    alleles: dict[str, float] = {}
    for name, freq in allele_frequencies.items():
        norm_name = name if name.startswith("*") else f"*{name}"
        alleles[norm_name] = freq

    # Verify frequencies sum to ~1.0
    total_freq = sum(alleles.values())
    if abs(total_freq - 1.0) > 0.05:
        logger.warning(
            "Allele frequencies for %s sum to %.4f (expected ~1.0). Normalizing.",
            gene_upper,
            total_freq,
        )
        if total_freq > 0:
            alleles = {k: v / total_freq for k, v in alleles.items()}

    # Calculate diplotype frequencies under HWE
    allele_names = sorted(alleles.keys())
    phenotype_freqs: dict[str, float] = {
        p.abbreviation: 0.0 for p in MetabolizerPhenotype if p != MetabolizerPhenotype.INDETERMINATE
    }

    scoring_table = _ACTIVITY_SCORE_TABLES.get(gene_upper, {})

    for i, a1 in enumerate(allele_names):
        for j, a2 in enumerate(allele_names):
            if j < i:
                continue  # Only count each pair once

            freq1 = alleles[a1]
            freq2 = alleles[a2]

            # HWE: p^2 for homozygous, 2pq for heterozygous
            if i == j:
                diplotype_freq = freq1 * freq2
            else:
                diplotype_freq = 2 * freq1 * freq2

            # Get activity score for this diplotype
            score1 = scoring_table.get(a1, 1.0)
            score2 = scoring_table.get(a2, 1.0)
            activity = score1 + score2

            phenotype = predict_metabolizer_status(activity, gene_upper)
            abbrev = phenotype.abbreviation

            if abbrev in phenotype_freqs:
                phenotype_freqs[abbrev] += diplotype_freq

    # Normalize in case of rounding
    total = sum(phenotype_freqs.values())
    if total > 0 and abs(total - 1.0) > 0.001:
        phenotype_freqs = {k: v / total for k, v in phenotype_freqs.items()}

    # Remove categories with zero frequency
    phenotype_freqs = {k: round(v, 6) for k, v in phenotype_freqs.items() if v > 0.0}

    logger.info(
        "Population phenotype frequencies for %s: %s",
        gene_upper,
        {k: f"{v:.4f}" for k, v in phenotype_freqs.items()},
    )

    return phenotype_freqs
