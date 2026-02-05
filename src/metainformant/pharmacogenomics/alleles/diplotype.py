"""Diplotype determination and activity scoring for pharmacogenes.

Implements diplotype construction from two star alleles, activity score calculation
per CPIC guidelines, ambiguous diplotype resolution, and phased diplotype inference.

A diplotype represents the combination of two haplotypes (one from each chromosome)
for a given gene. The activity score system converts diplotype information into a
numerical value that maps to metabolizer phenotype categories.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from metainformant.core.utils.logging import get_logger

from .star_allele import StarAllele, call_star_alleles, load_allele_definitions

logger = get_logger(__name__)


@dataclass
class Diplotype:
    """Represents a pharmacogene diplotype (pair of star alleles).

    Attributes:
        allele1: First star allele (typically the lower-numbered or reference)
        allele2: Second star allele
        gene: Gene symbol
        activity_score: Combined activity score for the diplotype
        phased: Whether the diplotype assignment is phased (True) or inferred (False)
        confidence: Confidence in the diplotype call ("high", "moderate", "low")
    """

    allele1: str
    allele2: str
    gene: str
    activity_score: float = 0.0
    phased: bool = False
    confidence: str = "high"

    def __post_init__(self) -> None:
        """Normalize diplotype representation (sort alleles alphabetically)."""
        self.gene = self.gene.upper()
        if not self.allele1.startswith("*"):
            self.allele1 = f"*{self.allele1}"
        if not self.allele2.startswith("*"):
            self.allele2 = f"*{self.allele2}"
        # Canonical ordering: lower allele first
        if self.allele1 > self.allele2:
            self.allele1, self.allele2 = self.allele2, self.allele1

    @property
    def diplotype_string(self) -> str:
        """Return formatted diplotype string (e.g., '*1/*2')."""
        return f"{self.allele1}/{self.allele2}"

    @property
    def is_homozygous(self) -> bool:
        """Return True if both alleles are the same."""
        return self.allele1 == self.allele2

    @property
    def is_reference(self) -> bool:
        """Return True if both alleles are reference (*1/*1)."""
        return self.allele1 == "*1" and self.allele2 == "*1"


# ── Activity score tables per CPIC guidelines ──────────────────────────────────
# Maps gene -> allele name -> activity value.
# These are derived from CPIC allele functionality tables.

_ACTIVITY_SCORE_TABLES: dict[str, dict[str, float]] = {
    "CYP2D6": {
        "*1": 1.0,
        "*2": 1.0,
        "*33": 1.0,
        "*35": 1.0,
        "*9": 0.5,
        "*10": 0.25,
        "*17": 0.5,
        "*29": 0.5,
        "*41": 0.5,
        "*3": 0.0,
        "*4": 0.0,
        "*5": 0.0,
        "*6": 0.0,
        "*7": 0.0,
        "*8": 0.0,
        "*11": 0.0,
        "*12": 0.0,
        "*14": 0.0,
        "*15": 0.0,
        "*36": 0.0,
    },
    "CYP2C19": {
        "*1": 1.0,
        "*17": 1.5,
        "*2": 0.0,
        "*3": 0.0,
        "*4": 0.0,
        "*5": 0.0,
        "*6": 0.0,
        "*7": 0.0,
        "*8": 0.0,
    },
    "CYP2C9": {
        "*1": 1.0,
        "*2": 0.5,
        "*11": 0.5,
        "*3": 0.0,
        "*5": 0.0,
        "*6": 0.0,
        "*8": 0.0,
    },
    "CYP3A5": {
        "*1": 1.0,
        "*3": 0.0,
        "*6": 0.0,
        "*7": 0.0,
    },
    "DPYD": {
        "*1": 1.0,
        "c.2846A>T": 0.5,
        "HapB3": 0.5,
        "*2A": 0.0,
        "*13": 0.0,
    },
    "TPMT": {
        "*1": 1.0,
        "*2": 0.0,
        "*3A": 0.0,
        "*3B": 0.0,
        "*3C": 0.0,
    },
    "NUDT15": {
        "*1": 1.0,
        "*2": 0.0,
        "*3": 0.0,
    },
    "SLCO1B1": {
        "*1": 1.0,
        "*1B": 1.0,
        "*5": 0.5,
        "*15": 0.5,
    },
    "UGT1A1": {
        "*1": 1.0,
        "*6": 0.5,
        "*28": 0.5,
        "*37": 0.0,
    },
}


def determine_diplotype(
    allele1: str | StarAllele,
    allele2: str | StarAllele,
    gene: str,
    scoring_table: dict[str, float] | None = None,
) -> Diplotype:
    """Construct a diplotype from two star alleles with activity scoring.

    Creates a Diplotype object from two allele designations, looking up
    activity values from the gene-specific scoring table to compute the
    combined activity score.

    Args:
        allele1: First allele name (str) or StarAllele object
        allele2: Second allele name (str) or StarAllele object
        gene: Gene symbol (e.g., "CYP2D6")
        scoring_table: Optional custom activity score mapping (allele_name -> score).
            If None, uses built-in CPIC-derived tables.

    Returns:
        Diplotype object with computed activity score
    """
    gene_upper = gene.upper()

    # Extract allele names
    name1 = allele1.name if isinstance(allele1, StarAllele) else allele1
    name2 = allele2.name if isinstance(allele2, StarAllele) else allele2

    if not name1.startswith("*"):
        name1 = f"*{name1}"
    if not name2.startswith("*"):
        name2 = f"*{name2}"

    # Calculate activity score
    activity = calculate_activity_score_from_alleles(name1, name2, gene_upper, scoring_table)

    diplotype = Diplotype(
        allele1=name1,
        allele2=name2,
        gene=gene_upper,
        activity_score=activity,
    )

    logger.info(
        "Determined diplotype %s for %s (activity score: %.2f)",
        diplotype.diplotype_string,
        gene_upper,
        activity,
    )

    return diplotype


def calculate_activity_score_from_alleles(
    allele1_name: str,
    allele2_name: str,
    gene: str,
    scoring_table: dict[str, float] | None = None,
) -> float:
    """Calculate the combined activity score for two alleles.

    Args:
        allele1_name: First allele name
        allele2_name: Second allele name
        gene: Gene symbol
        scoring_table: Optional custom scoring table

    Returns:
        Combined activity score (sum of individual allele scores)
    """
    gene_upper = gene.upper()
    table = scoring_table if scoring_table is not None else _ACTIVITY_SCORE_TABLES.get(gene_upper, {})

    # Normalize allele names
    name1 = allele1_name if allele1_name.startswith("*") else f"*{allele1_name}"
    name2 = allele2_name if allele2_name.startswith("*") else f"*{allele2_name}"

    score1 = table.get(name1, 1.0)  # Default to 1.0 (normal function) if unknown
    score2 = table.get(name2, 1.0)

    if name1 not in table:
        logger.warning("Allele %s not found in scoring table for %s, defaulting to 1.0", name1, gene_upper)
    if name2 not in table:
        logger.warning("Allele %s not found in scoring table for %s, defaulting to 1.0", name2, gene_upper)

    return score1 + score2


def calculate_activity_score(
    diplotype: Diplotype,
    gene: str,
    scoring_table: dict[str, float] | None = None,
) -> float:
    """Calculate activity score for an existing Diplotype.

    Args:
        diplotype: Diplotype object
        gene: Gene symbol
        scoring_table: Optional custom scoring table

    Returns:
        Combined activity score, also updates the Diplotype object
    """
    score = calculate_activity_score_from_alleles(diplotype.allele1, diplotype.allele2, gene.upper(), scoring_table)
    diplotype.activity_score = score
    return score


def resolve_ambiguous_diplotypes(
    possible_diplotypes: list[Diplotype],
) -> Diplotype:
    """Resolve ambiguous diplotype calls when phase is unknown.

    When genotype data is unphased, multiple diplotype configurations may be
    possible. This function selects the most clinically relevant assignment
    using the following priority rules (per CPIC recommendations):

    1. Prefer diplotypes with known activity scores for both alleles
    2. Among scored diplotypes, prefer the one with the highest activity score
       (conservative: assume best-case metabolism)
    3. If tied, prefer diplotypes containing more common alleles (*1, *2)

    Args:
        possible_diplotypes: List of possible Diplotype objects

    Returns:
        The selected most likely Diplotype

    Raises:
        ValueError: If possible_diplotypes is empty
    """
    if not possible_diplotypes:
        raise ValueError("Cannot resolve empty list of diplotypes")

    if len(possible_diplotypes) == 1:
        return possible_diplotypes[0]

    # Score each diplotype for ranking
    scored: list[tuple[float, int, Diplotype]] = []
    for dip in possible_diplotypes:
        gene = dip.gene
        table = _ACTIVITY_SCORE_TABLES.get(gene, {})

        # Count how many alleles have known scores
        known_count = 0
        if dip.allele1 in table:
            known_count += 1
        if dip.allele2 in table:
            known_count += 1

        scored.append((dip.activity_score, known_count, dip))

    # Sort: highest known_count first, then highest activity score (conservative)
    scored.sort(key=lambda x: (x[1], x[0]), reverse=True)

    selected = scored[0][2]
    selected.confidence = "moderate" if len(possible_diplotypes) > 1 else "high"

    logger.info(
        "Resolved ambiguous diplotypes (%d candidates) -> %s (confidence: %s)",
        len(possible_diplotypes),
        selected.diplotype_string,
        selected.confidence,
    )

    return selected


def phased_diplotype(
    variants: set[str] | frozenset[str],
    phase_data: dict[str, int],
    gene: str,
    allele_definitions: list[StarAllele] | None = None,
) -> Diplotype:
    """Determine diplotype using phasing information.

    When phase information is available (e.g., from long reads, family data,
    or statistical phasing), this function assigns variants to maternal and
    paternal haplotypes before calling star alleles on each independently.

    Args:
        variants: Set of all observed variant identifiers
        phase_data: Mapping of variant_id -> phase (0 for haplotype A, 1 for haplotype B).
            Variants not in phase_data are treated as unphased.
        gene: Gene symbol (e.g., "CYP2D6")
        allele_definitions: Optional pre-loaded allele definitions

    Returns:
        Diplotype with phased=True and alleles called from each haplotype
    """
    gene_upper = gene.upper()

    if allele_definitions is None:
        allele_definitions = load_allele_definitions(gene_upper)

    # Partition variants by phase
    hap_a_variants: set[str] = set()
    hap_b_variants: set[str] = set()
    unphased: set[str] = set()

    for v in variants:
        phase = phase_data.get(v)
        if phase == 0:
            hap_a_variants.add(v)
        elif phase == 1:
            hap_b_variants.add(v)
        else:
            unphased.add(v)

    if unphased:
        logger.warning(
            "%d variants for %s lack phase information: %s",
            len(unphased),
            gene_upper,
            sorted(unphased),
        )
        # Add unphased variants to both haplotypes (conservative)
        hap_a_variants.update(unphased)
        hap_b_variants.update(unphased)

    # Call star alleles independently on each haplotype
    alleles_a = call_star_alleles(hap_a_variants, gene_upper, allele_definitions)
    alleles_b = call_star_alleles(hap_b_variants, gene_upper, allele_definitions)

    # Take the best (most specific) allele call from each haplotype
    allele1 = alleles_a[0] if alleles_a else StarAllele(name="*1", gene=gene_upper)
    allele2 = alleles_b[0] if alleles_b else StarAllele(name="*1", gene=gene_upper)

    diplotype = determine_diplotype(allele1, allele2, gene_upper)
    diplotype.phased = True
    diplotype.confidence = "high"

    logger.info(
        "Phased diplotype for %s: %s (hap_a=%d variants, hap_b=%d variants, unphased=%d)",
        gene_upper,
        diplotype.diplotype_string,
        len(hap_a_variants - unphased),
        len(hap_b_variants - unphased),
        len(unphased),
    )

    return diplotype
