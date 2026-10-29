"""Diplotype determination and activity scoring for pharmacogenes.

Implements diplotype construction from two star alleles, activity score calculation
per CPIC guidelines, ambiguous diplotype resolution, and phased diplotype inference.

A diplotype represents the combination of two haplotypes (one from each chromosome)
for a given gene. The activity score system converts diplotype information into a
numerical value that maps to metabolizer phenotype categories.

An allele with no activity value in the scoring table is uncharacterized, and it is
never scored as normal function: the diplotype's activity score is ``None`` and
:func:`~metainformant.pharmacogenomics.alleles.phenotype.classify_phenotype` calls the
phenotype indeterminate, as CPIC does.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass

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
        activity_score: Combined activity score for the diplotype, or ``None`` (the default)
            when it has not been scored. ``None`` is not 0.0: it never reads as a poor
            metabolizer. Score it with :func:`calculate_activity_score`, or build the
            diplotype with :func:`determine_diplotype`, which does; ``classify_phenotype``
            computes the score from the alleles of an unscored diplotype. ``None`` also
            results when an allele has no activity value in the scoring table.
        phased: Whether the diplotype assignment is phased (True) or inferred (False)
        confidence: Confidence in the diplotype call ("high", "moderate", "low")
    """

    allele1: str
    allele2: str
    gene: str
    activity_score: float | None = None
    phased: bool = False
    confidence: str = "high"

    def __post_init__(self) -> None:
        """Normalize diplotype representation (sort alleles alphabetically)."""
        self.gene = self.gene.upper()
        self.allele1 = _normalize_allele_name(self.allele1, self.gene)
        self.allele2 = _normalize_allele_name(self.allele2, self.gene)
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
# Maps gene -> allele name -> activity value, from the CPIC allele functionality tables
# (api.cpicpgx.org/v1/allele, checked 2026-10-27): normal function 1.0, decreased
# function 0.5, no function 0.0, increased function CYP2C19 *17 1.5. CYP2D6 carries
# CPIC's own activity values, including 0.25 for *9, *10 and *41. CPIC assigns no
# activity values for CYP2C19, SLCO1B1 or UGT1A1; their entries encode the CPIC function
# class on the same scale, which phenotype.py reads pairwise. CYP1A2 is the one gene in
# the table with no CPIC source at all (see its entry).
#
# An allele that is absent from a gene's table has no activity value. It is never scored
# as normal function: calculate_activity_score_from_alleles() returns None for the
# diplotype, and phenotype.classify_phenotype() reports it as indeterminate, naming the
# allele (CPIC calls a diplotype with an uncharacterized allele indeterminate).
#
# This is the package's single allele-function table. metabolism.metabolizer_status
# derives its table from it, and the star-allele calling definitions in
# star_allele._BUILTIN_ALLELE_DEFINITIONS must carry the same activity values
# (both are enforced by tests). Phenotype assignment from these values lives in
# phenotype.py.

_ACTIVITY_SCORE_TABLES: dict[str, dict[str, float]] = {
    "CYP2D6": {
        "*1": 1.0,
        "*2": 1.0,
        "*33": 1.0,
        "*35": 1.0,
        "*9": 0.25,
        "*10": 0.25,
        "*17": 0.5,
        "*29": 0.5,
        "*41": 0.25,
        "*3": 0.0,
        "*4": 0.0,
        "*5": 0.0,
        "*6": 0.0,
        "*7": 0.0,
        "*8": 0.0,
        "*11": 0.0,
        "*12": 0.0,
        "*14": 0.5,
        "*15": 0.0,
        "*36": 0.0,
        # Symbolic duplication placeholders (xN = N >= 2 copies of the allele). 2.0 is exact
        # for x2 only; handle_cyp2d6_cnv() scores a real copy number.
        "*1xN": 2.0,
        "*2xN": 2.0,
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
        "*9": 0.5,
        "*10": 0.5,
    },
    "CYP2C9": {
        "*1": 1.0,
        "*2": 0.5,
        "*4": 0.5,
        "*11": 0.5,
        "*3": 0.0,
        "*5": 0.5,
        "*6": 0.0,
        "*8": 0.5,
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
        "*4": 0.0,
    },
    "NUDT15": {
        "*1": 1.0,
        "*2": 0.0,
        "*3": 0.0,
    },
    "SLCO1B1": {
        "*1": 1.0,
        "*1B": 1.0,
        "*5": 0.0,
        "*15": 0.0,
    },
    # UGT1A1 and CYP1A2 have activity scores only: no star-allele calling definitions are
    # bundled for them (see star_allele._GENES_WITHOUT_BUILTIN_DEFINITIONS for why), so
    # their diplotypes must be built from allele names called elsewhere.
    # CPIC lists *36 as increased function but gives *1/*36 and *36/*36 the normal
    # phenotype (Gammal et al., Clin Pharmacol Ther 2016, Table 2), so it carries 1.0.
    "UGT1A1": {
        "*1": 1.0,
        "*6": 0.5,
        "*28": 0.5,
        "*36": 1.0,
        "*37": 0.5,
    },
    # CYP1A2 does NOT come from CPIC. CPIC lists no CYP1A2 allele values and defines no CYP1A2
    # phenotype (checked against api.cpicpgx.org by the maintainers on 2026-10-27), so these
    # four alleles can be scored but not classified: phenotype.py raises ValueError for CYP1A2
    # instead of borrowing another gene's bands.
    # Provenance of the values: unattributed. The table first appeared in the metabolism module
    # with the inline labels "Normal function" (*1A), "Decreased function" (*1C, *1K) and
    # "Increased inducibility" (*1F) and no citation, and no source for them is cited anywhere
    # in this repository: not PharmVar's function assignments, not the literature on the
    # inducibility of *1F. They have not been checked against either. 0.5 and 1.5 are only the
    # generic decreased/increased-function values this table uses for the genes CPIC does
    # score, not measured or published CYP1A2 activities. Treat them as unsourced until they
    # are traced to PharmVar or primary literature, and do not read them as CPIC values.
    "CYP1A2": {
        "*1A": 1.0,
        "*1C": 0.5,
        "*1F": 1.5,
        "*1K": 0.5,
    },
}


def _normalize_allele_name(name: str, gene: str, scoring_table: Mapping[str, float] | None = None) -> str:
    """Normalize an allele name for scoring-table lookup.

    Star allele names are prefixed with ``*`` when missing. Names already
    present in the scoring table are returned unchanged, so non-star allele
    designations (e.g. DPYD's "c.2846A>T" and "HapB3") are preserved.
    """
    if name.startswith("*"):
        return name
    table = scoring_table if scoring_table is not None else _ACTIVITY_SCORE_TABLES.get(gene.upper(), {})
    return name if name in table else f"*{name}"


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
        Diplotype object with its computed activity score. The score is ``None`` when
        either allele has no activity value in the scoring table (an uncharacterized
        allele; a warning names it): such a diplotype is indeterminate, not normal.
    """
    gene_upper = gene.upper()

    # Extract allele names
    name1 = allele1.name if isinstance(allele1, StarAllele) else allele1
    name2 = allele2.name if isinstance(allele2, StarAllele) else allele2

    # Calculate activity score
    activity = calculate_activity_score_from_alleles(name1, name2, gene_upper, scoring_table)

    diplotype = Diplotype(
        allele1=name1,
        allele2=name2,
        gene=gene_upper,
        activity_score=activity,
    )

    logger.info(
        "Determined diplotype %s for %s (activity score: %s)",
        diplotype.diplotype_string,
        gene_upper,
        "none (uncharacterized allele)" if activity is None else f"{activity:.2f}",
    )

    return diplotype


def uncharacterized_alleles(
    allele1_name: str,
    allele2_name: str,
    gene: str,
    scoring_table: Mapping[str, float] | None = None,
) -> list[str]:
    """Return the alleles of a diplotype that have no activity value in the scoring table.

    Names are normalized as in :func:`calculate_activity_score_from_alleles`, listed in the
    order given and without repeats (a homozygous diplotype names its allele once).

    Args:
        allele1_name: First allele name
        allele2_name: Second allele name
        gene: Gene symbol
        scoring_table: Optional custom scoring table. If None, uses the built-in table for the
            gene (an empty table for a gene without one, so every allele is uncharacterized).

    Returns:
        The normalized names of the uncharacterized alleles; an empty list when both are scored.
    """
    gene_upper = gene.upper()
    table = scoring_table if scoring_table is not None else _ACTIVITY_SCORE_TABLES.get(gene_upper, {})
    names = (
        _normalize_allele_name(allele1_name, gene_upper, table),
        _normalize_allele_name(allele2_name, gene_upper, table),
    )
    return [name for name in dict.fromkeys(names) if name not in table]


def calculate_activity_score_from_alleles(
    allele1_name: str,
    allele2_name: str,
    gene: str,
    scoring_table: dict[str, float] | None = None,
) -> float | None:
    """Calculate the combined activity score for two alleles.

    Args:
        allele1_name: First allele name
        allele2_name: Second allele name
        gene: Gene symbol
        scoring_table: Optional custom scoring table

    Returns:
        Combined activity score (sum of individual allele scores), or ``None`` if either allele
        has no activity value in the table. An uncharacterized allele is never assumed to be
        normal function (1.0); a warning names it, and the diplotype is indeterminate.
    """
    gene_upper = gene.upper()
    table = scoring_table if scoring_table is not None else _ACTIVITY_SCORE_TABLES.get(gene_upper, {})

    unknown = uncharacterized_alleles(allele1_name, allele2_name, gene_upper, table)
    if unknown:
        logger.warning(
            "No activity value for allele(s) %s in the %s scoring table; the diplotype has no activity "
            "score and its phenotype is indeterminate (pass scoring_table= with a value for the allele to score it)",
            ", ".join(unknown),
            gene_upper,
        )
        return None

    # Normalize allele names (non-star designations already in the table are kept)
    name1 = _normalize_allele_name(allele1_name, gene_upper, table)
    name2 = _normalize_allele_name(allele2_name, gene_upper, table)
    return table[name1] + table[name2]


def calculate_activity_score(
    diplotype: Diplotype,
    gene: str,
    scoring_table: dict[str, float] | None = None,
) -> float | None:
    """Calculate activity score for an existing Diplotype.

    Args:
        diplotype: Diplotype object
        gene: Gene symbol
        scoring_table: Optional custom scoring table

    Returns:
        Combined activity score, also updates the Diplotype object; ``None`` (stored on the
        diplotype too) if an allele has no activity value in the table, see
        :func:`calculate_activity_score_from_alleles`.
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

    A diplotype that has not been scored (``activity_score`` is ``None``) ranks by the score
    its alleles imply; if one of its alleles has no activity value it has no score and ranks
    after every diplotype with as many known alleles.

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

        score = dip.activity_score
        if score is None and known_count == 2:
            score = table[dip.allele1] + table[dip.allele2]
        scored.append((-math.inf if score is None else score, known_count, dip))

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


def _haplotype_allele(called: list[StarAllele], gene: str, haplotype: str) -> tuple[StarAllele, bool]:
    """Return the one allele a haplotype carries, and whether more than one was called for it.

    A haplotype is a single star allele. When the variants on it match several definitions
    (for example a loss-of-function marker beside a normal-function one), the most impaired
    allele is taken: a no-function variant is not offset by a normal-function marker on the same
    haplotype, and the calls come back ordered by name, so taking the first would pick "*2" over
    "*4". The ambiguity is logged and returned so the caller can lower its confidence.

    Args:
        called: Alleles called from the haplotype's variants (:func:`call_star_alleles`)
        gene: Gene symbol
        haplotype: Label of the haplotype, for the log message

    Returns:
        ``(allele, ambiguous)``. An empty ``called`` gives the reference ``*1``.
    """
    if not called:
        return StarAllele(name="*1", gene=gene), False
    if len(called) == 1:
        return called[0], False
    chosen = min(called, key=lambda allele: (allele.activity_value, allele.name))
    logger.warning(
        "Haplotype %s of %s matches %d alleles (%s); calling %s, the most impaired function",
        haplotype,
        gene,
        len(called),
        ", ".join(allele.name for allele in called),
        chosen.name,
    )
    return chosen, True


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

    A haplotype carries one allele. If its variants match more than one definition, the allele
    with the lowest activity value is taken (a loss-of-function variant is not offset by a
    normal-function marker on the same haplotype), a warning names the candidates, and the
    diplotype's confidence is ``"moderate"`` instead of ``"high"``.

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

    # One allele per haplotype: the most impaired one when several definitions match
    allele1, ambiguous_a = _haplotype_allele(alleles_a, gene_upper, "A")
    allele2, ambiguous_b = _haplotype_allele(alleles_b, gene_upper, "B")

    diplotype = determine_diplotype(allele1, allele2, gene_upper)
    diplotype.phased = True
    diplotype.confidence = "moderate" if ambiguous_a or ambiguous_b else "high"

    logger.info(
        "Phased diplotype for %s: %s (hap_a=%d variants, hap_b=%d variants, unphased=%d)",
        gene_upper,
        diplotype.diplotype_string,
        len(hap_a_variants - unphased),
        len(hap_b_variants - unphased),
        len(unphased),
    )

    return diplotype
