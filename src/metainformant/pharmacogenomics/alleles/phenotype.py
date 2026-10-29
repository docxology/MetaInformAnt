"""Metabolizer phenotype prediction from pharmacogene diplotypes.

Maps diplotypes to standardized metabolizer phenotype categories per CPIC
guidelines. CPIC defines the phenotype of CYP2D6, CYP2C9 and DPYD from the
diplotype's total activity score and of the other genes from the function of the
two alleles. Both definitions are implemented: ``_PHENOTYPE_THRESHOLDS`` holds the
score bands of every gene (for TPMT, NUDT15 and CYP3A5 they reproduce the
function-based definition exactly), and ``_FUNCTION_PHENOTYPES`` the pairwise
allele-function rules of CYP2C19, SLCO1B1 and UGT1A1, where a sum of activity
values cannot always tell the allele pairs apart. Also supports population-level
phenotype frequency estimation from allele frequency data.

Only genes that CPIC defines a phenotype for are classified (``supported_phenotype_genes()``).
For any other gene, including CYP1A2 whose alleles have activity values but no CPIC
phenotype, every function here that assigns a phenotype or lists bands raises
``ValueError`` naming the gene and the supported genes: no other gene's bands are borrowed.
A diplotype with an allele that has no activity value is indeterminate (CPIC's term),
never scored as normal function.

Metabolizer categories:
    - Poor Metabolizer (PM): No or very low enzyme activity
    - Intermediate Metabolizer (IM): Reduced activity between PM and NM
    - Normal Metabolizer (NM): Typical activity level (reference)
    - Rapid Metabolizer (RM): Increased activity (one increased-function allele)
    - Ultrarapid Metabolizer (UM): Very high activity (gene duplication or two increased-function alleles)
    - Indeterminate (IND): No phenotype can be assigned (an uncharacterized allele)
"""

from __future__ import annotations

import math
import re
from collections.abc import Sequence
from enum import Enum, IntEnum
from typing import Any

from metainformant.core.utils.logging import get_logger

from .diplotype import (
    _ACTIVITY_SCORE_TABLES,
    Diplotype,
    _normalize_allele_name,
    determine_diplotype,
    uncharacterized_alleles,
)

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

    @classmethod
    def from_string(cls, label: str) -> MetabolizerPhenotype:
        """Parse a phenotype label written in any of the forms used across the package.

        Accepts the full CPIC name ("Poor Metabolizer"), its abbreviation ("PM"), the bare
        lowercase token used by ``metabolism.metabolizer_status`` ("poor", "ultrarapid"),
        variants such as "ultra-rapid", and the pre-2017 CPIC synonym "extensive" for
        normal, all case-insensitively. A MetabolizerPhenotype is returned unchanged.

        Raises:
            ValueError: If the label is not a recognised phenotype.
        """
        if isinstance(label, cls):
            return label
        key = re.sub(r"[\s_-]+", "", str(label).casefold())
        for suffix in ("metabolizer", "metaboliser"):
            if key.endswith(suffix):
                key = key[: -len(suffix)]
                break
        phenotype = _PHENOTYPE_ALIASES.get(key)
        if phenotype is None:
            raise ValueError(f"Unrecognised metabolizer phenotype: {label!r}")
        return phenotype


_PHENOTYPE_ALIASES: dict[str, MetabolizerPhenotype] = {
    "poor": MetabolizerPhenotype.POOR,
    "pm": MetabolizerPhenotype.POOR,
    "intermediate": MetabolizerPhenotype.INTERMEDIATE,
    "im": MetabolizerPhenotype.INTERMEDIATE,
    "normal": MetabolizerPhenotype.NORMAL,
    "nm": MetabolizerPhenotype.NORMAL,
    "extensive": MetabolizerPhenotype.NORMAL,
    "em": MetabolizerPhenotype.NORMAL,
    "rapid": MetabolizerPhenotype.RAPID,
    "rm": MetabolizerPhenotype.RAPID,
    "ultrarapid": MetabolizerPhenotype.ULTRARAPID,
    "um": MetabolizerPhenotype.ULTRARAPID,
    "indeterminate": MetabolizerPhenotype.INDETERMINATE,
    "ind": MetabolizerPhenotype.INDETERMINATE,
}


_PM = MetabolizerPhenotype.POOR
_IM = MetabolizerPhenotype.INTERMEDIATE
_NM = MetabolizerPhenotype.NORMAL
_RM = MetabolizerPhenotype.RAPID
_UM = MetabolizerPhenotype.ULTRARAPID
_INF = math.inf

# ── Gene-specific phenotype bands over the activity score ──────────────────────
# The activity score of a diplotype is the sum of its two allele values in
# diplotype._ACTIVITY_SCORE_TABLES (normal function 1.0, decreased 0.5, no function 0.0,
# CYP2C19 *17 1.5, CYP2D6 *10 0.25). Each gene maps to consecutive bands written as
# (upper bound, upper bound inclusive?, phenotype): a band starts where the previous one
# ends (the first at 0) and the last runs to infinity, so the bands tile every score >= 0
# with no gap and no overlap (a test enforces this). Bounds follow the CPIC definitions
# cited per gene; they are not "round" cut-offs, so do not shift them without the source.

_PHENOTYPE_THRESHOLDS: dict[str, tuple[tuple[float, bool, MetabolizerPhenotype], ...]] = {
    # CPIC/DPWG 2019 consensus (Caudle et al., Clin Transl Sci 2020): PM AS = 0;
    # IM 0 < AS < 1.25 (0.25-1.0); NM 1.25 <= AS <= 2.25; UM AS > 2.25.
    "CYP2D6": (
        (0.0, True, _PM),
        (1.25, False, _IM),
        (2.25, True, _NM),
        (_INF, True, _UM),
    ),
    # CPIC defines CYP2C19 by allele function (Moriyama et al., Clin Pharmacol Ther 2017,
    # Table 1; Lee et al., Clin Pharmacol Ther 2022, Table 1), implemented by
    # _FUNCTION_PHENOTYPES, which classify_phenotype() uses. These bands are the same
    # mapping on the score scale, for callers that only hold a score: two no-function
    # alleles (AS 0) PM; one no-function allele (*1/*2 = 1.0, *2/*17 = 1.5) IM; *1/*1 (2.0)
    # NM; *1/*17 (2.5) RM; *17/*17 (3.0) UM. AS 0.5 (decreased + no function) is CPIC's
    # "likely poor". A sum cannot separate *9/*17 from *1/*1 (both 2.0) or *9/*9 from
    # *1/*2 (both 1.0); only the function rules resolve those.
    "CYP2C19": (
        (0.5, True, _PM),
        (2.0, False, _IM),
        (2.5, False, _NM),
        (3.0, False, _RM),
        (_INF, True, _UM),
    ),
    # CPIC CYP2C9 activity score (Karnes et al., Clin Pharmacol Ther 2021; Theken et al.,
    # Clin Pharmacol Ther 2020): NM AS = 2; IM AS 1-1.5; PM AS 0-0.5.
    "CYP2C9": (
        (0.5, True, _PM),
        (2.0, False, _IM),
        (_INF, True, _NM),
    ),
    # CPIC CYP3A5 / tacrolimus (Birdwell et al., Clin Pharmacol Ther 2015): NM *1/*1 (AS 2);
    # IM one *1 plus one no-function allele, e.g. *1/*3 (AS 1); PM two no-function alleles (AS 0).
    "CYP3A5": (
        (0.0, True, _PM),
        (2.0, False, _IM),
        (_INF, True, _NM),
    ),
    # CPIC DPYD activity score (Amstutz et al., Clin Pharmacol Ther 2018): NM AS = 2;
    # IM AS 1-1.5 (e.g. *1/c.2846A>T = 1.5); PM AS 0-0.5. The fluorouracil entries in
    # annotations.cpic carry the same phenotype labels.
    "DPYD": (
        (0.5, True, _PM),
        (1.5, True, _IM),
        (_INF, True, _NM),
    ),
    # CPIC thiopurines (Relling et al., Clin Pharmacol Ther 2019): NM two normal-function
    # alleles (AS 2); IM one no-function allele (AS 1); PM two no-function alleles (AS 0).
    "TPMT": (
        (0.0, True, _PM),
        (2.0, False, _IM),
        (_INF, True, _NM),
    ),
    "NUDT15": (
        (0.0, True, _PM),
        (2.0, False, _IM),
        (_INF, True, _NM),
    ),
    # CPIC statins (Cooper-DeHoff et al., Clin Pharmacol Ther 2022): normal function (two
    # normal alleles, AS 2) is NM; decreased function (one normal plus one no-function
    # allele, *1/*5 = 1.0) is IM; poor function (two no-function alleles, *5/*5 = 0) is PM.
    # CPIC classes *5 and *15 as no function. The simvastatin entries in annotations.cpic
    # use the same PM/IM/NM labels.
    "SLCO1B1": (
        (0.0, True, _PM),
        (2.0, False, _IM),
        (_INF, True, _NM),
    ),
    # CPIC atazanavir (Gammal et al., Clin Pharmacol Ther 2016, Table 2): NM two normal
    # alleles (AS 2); IM one decreased allele (*1/*28 = 1.5); PM two decreased alleles
    # (*28/*28 = 1.0).
    "UGT1A1": (
        (1.0, True, _PM),
        (2.0, False, _IM),
        (_INF, True, _NM),
    ),
}


def supported_phenotype_genes() -> list[str]:
    """Return the genes that have a CPIC phenotype definition, sorted.

    These are the only genes the functions of this module classify; every other gene raises
    ``ValueError``. A gene can have activity values and still not be listed: CYP1A2 can be
    scored (``diplotype.determine_diplotype``) but CPIC defines no phenotype for it.
    """
    return sorted(_PHENOTYPE_THRESHOLDS)


def _canonical_gene(gene: str) -> str:
    """Return the gene symbol as the tables spell it (stripped, upper case)."""
    return gene.strip().upper()


def _unsupported_gene_error(gene_upper: str) -> ValueError:
    """Build the ValueError for a gene without a CPIC phenotype definition."""
    if gene_upper in _ACTIVITY_SCORE_TABLES:
        detail = " (it has activity values, so it can be scored, but CPIC defines no phenotype for it)"
    else:
        detail = ""
    return ValueError(
        f"No CPIC phenotype definition for gene {gene_upper!r}{detail}. "
        f"Supported genes: {', '.join(supported_phenotype_genes())}"
    )


def _phenotype_bands(gene_upper: str) -> tuple[tuple[float, bool, MetabolizerPhenotype], ...]:
    """Return a gene's score bands.

    Raises:
        ValueError: If CPIC defines no phenotype for the gene. No other gene's bands stand in.
    """
    bands = _PHENOTYPE_THRESHOLDS.get(gene_upper)
    if bands is None:
        raise _unsupported_gene_error(gene_upper)
    return bands


# ── Phenotype from the function of the two alleles ─────────────────────────────
# CYP2C19, SLCO1B1 and UGT1A1 are defined by CPIC as pairs of allele functions, and a sum
# of activity values cannot always tell the pairs apart (see CYP2C19 above), so
# classify_phenotype() translates these genes pairwise. An allele's function class is read
# from its activity value: 0 no function, below 1 decreased, 1 normal, above 1 increased.


class _AlleleFunction(IntEnum):
    """CPIC allele function classes, ordered from least to most active."""

    NO_FUNCTION = 0
    DECREASED = 1
    NORMAL = 2
    INCREASED = 3


_NO, _DEC, _NORM, _INC = (
    _AlleleFunction.NO_FUNCTION,
    _AlleleFunction.DECREASED,
    _AlleleFunction.NORMAL,
    _AlleleFunction.INCREASED,
)

# Genes with reduced-function alleles only (SLCO1B1, UGT1A1): two normal alleles NM; a
# normal allele beside one reduced (decreased or no function) allele IM; two reduced
# alleles PM. Keys are (less active, more active) allele function classes. Pairs with an
# increased-function allele are absent on purpose: there is no category for them here.
_REDUCED_FUNCTION_PHENOTYPES: dict[tuple[_AlleleFunction, _AlleleFunction], MetabolizerPhenotype] = {
    (_NORM, _NORM): _NM,
    (_DEC, _NORM): _IM,
    (_NO, _NORM): _IM,
    (_DEC, _DEC): _PM,
    (_NO, _DEC): _PM,
    (_NO, _NO): _PM,
}

_FUNCTION_PHENOTYPES: dict[str, dict[tuple[_AlleleFunction, _AlleleFunction], MetabolizerPhenotype]] = {
    # Moriyama et al., Clin Pharmacol Ther 2017, Table 1 and Lee et al., Clin Pharmacol Ther
    # 2022, Table 1: UM two increased-function alleles (*17/*17); RM one normal plus one
    # increased (*1/*17); NM two normal (*1/*1); IM one no-function allele with a normal or
    # increased one (*1/*2, *2/*17); PM two no-function alleles (*2/*2, *3/*3). CPIC calls the
    # pairs involving a decreased-function allele (*9, *10) "likely intermediate" (with a
    # normal or increased allele) and "likely poor" (with a decreased or no-function one);
    # MetabolizerPhenotype has no "likely" categories, so they map to IM and PM.
    "CYP2C19": {
        (_INC, _INC): _UM,
        (_NORM, _INC): _RM,
        (_NORM, _NORM): _NM,
        (_NO, _NORM): _IM,
        (_NO, _INC): _IM,
        (_DEC, _NORM): _IM,
        (_DEC, _INC): _IM,
        (_NO, _NO): _PM,
        (_DEC, _DEC): _PM,
        (_NO, _DEC): _PM,
    },
    # Cooper-DeHoff et al., Clin Pharmacol Ther 2022: normal / decreased / poor function.
    "SLCO1B1": _REDUCED_FUNCTION_PHENOTYPES,
    # Gammal et al., Clin Pharmacol Ther 2016, Table 2: *1/*1 normal; *1/*28, *1/*37 intermediate;
    # *28/*28, *28/*37, *37/*37 poor. Going by function class keeps these right whether an
    # allele's activity value is recorded as decreased or as no function.
    "UGT1A1": _REDUCED_FUNCTION_PHENOTYPES,
}


def _function_class(activity: float) -> _AlleleFunction:
    """Return the CPIC function class implied by one allele's activity value."""
    if math.isnan(activity) or activity < 0:
        raise ValueError(f"Allele activity values must be non-negative numbers, got {activity!r}")
    if math.isclose(activity, 0.0, abs_tol=1e-9):
        return _AlleleFunction.NO_FUNCTION
    if math.isclose(activity, 1.0):
        return _AlleleFunction.NORMAL
    return _AlleleFunction.DECREASED if activity < 1.0 else _AlleleFunction.INCREASED


def _band_phenotype(
    bands: tuple[tuple[float, bool, MetabolizerPhenotype], ...],
    activity_score: float,
) -> MetabolizerPhenotype:
    """Return the phenotype of the band containing ``activity_score`` (>= 0)."""
    for upper, inclusive, phenotype in bands[:-1]:
        if activity_score < upper or (inclusive and activity_score == upper):
            return phenotype
    return bands[-1][2]


def _phenotype_for_alleles(gene: str, activity1: float, activity2: float) -> MetabolizerPhenotype:
    """Phenotype of a diplotype from its two alleles' activity values.

    Genes CPIC defines by allele function are translated from the pair of function classes;
    all others from the summed activity score.

    Raises:
        ValueError: If a value is negative or NaN, the gene has no CPIC phenotype definition, or
            the gene has no phenotype for that pair of allele functions (an increased-function
            SLCO1B1 or UGT1A1 allele).
    """
    gene_upper = _canonical_gene(gene)
    rules = _FUNCTION_PHENOTYPES.get(gene_upper)
    if rules is None:
        return predict_metabolizer_status(activity1 + activity2, gene_upper)

    less_active, more_active = sorted((_function_class(activity1), _function_class(activity2)))
    phenotype = rules.get((less_active, more_active))
    if phenotype is None:
        raise ValueError(
            f"No {gene_upper} phenotype is defined for a diplotype of {less_active.name.lower().replace('_', ' ')} "
            f"and {more_active.name.lower().replace('_', ' ')} alleles"
        )
    return phenotype


def predict_metabolizer_status(
    activity_score: float,
    gene: str,
) -> MetabolizerPhenotype:
    """Map an activity score to a metabolizer phenotype category.

    Uses the gene-specific bands in ``_PHENOTYPE_THRESHOLDS`` (CPIC definitions, cited
    there). The bands cover every score >= 0 without gaps. For example, the CYP2D6 bands:
        AS = 0 -> PM, 0 < AS < 1.25 -> IM, 1.25 <= AS <= 2.25 -> NM, AS > 2.25 -> UM

    A gene CPIC defines no phenotype for has no bands, and no other gene's bands stand in
    for it: it raises (see :func:`supported_phenotype_genes`).

    A score is a sum of allele values, so for CYP2C19 it cannot always identify the
    alleles behind it (a decreased-function allele makes *9/*17 sum to the same 2.0 as
    *1/*1). Use :func:`classify_phenotype` when the diplotype is known.

    Args:
        activity_score: Combined activity score (typically sum of two allele scores)
        gene: Gene symbol (e.g., "CYP2D6")

    Returns:
        MetabolizerPhenotype enum value

    Raises:
        ValueError: If the activity score is negative, NaN or None (an unscored diplotype has
            no score to classify), or the gene has no CPIC phenotype definition; the message
            names the gene and lists the supported genes
    """
    if activity_score is None or math.isnan(activity_score) or activity_score < 0:
        raise ValueError(f"Activity score must be a non-negative number, got {activity_score!r}")

    return _band_phenotype(_phenotype_bands(_canonical_gene(gene)), activity_score)


def get_phenotype_thresholds(
    gene: str,
) -> list[dict[str, Any]]:
    """Get gene-specific phenotype threshold definitions.

    These are the activity-score bands :func:`predict_metabolizer_status` applies. For
    CYP2C19, SLCO1B1 and UGT1A1, which :func:`classify_phenotype` translates from allele
    function, they are the equivalent bands on the score scale.

    Args:
        gene: Gene symbol (e.g., "CYP2D6")

    Returns:
        List of threshold dictionaries, lowest activity first, with keys:
            - "lower_bound": Lower activity score boundary
            - "upper_bound": Upper activity score boundary (infinity for the last band)
            - "lower_inclusive": Whether a score equal to the lower bound is in the band
            - "upper_inclusive": Whether a score equal to the upper bound is in the band
            - "phenotype": MetabolizerPhenotype value
            - "abbreviation": Phenotype abbreviation (PM, IM, NM, RM, UM)
            - "gene_specific": Always True. It used to be False for the default (CYP2D6) bands
              returned for a gene with none of its own; such a gene now raises. The key is kept
              so callers written before that still work.

    Raises:
        ValueError: If the gene has no CPIC phenotype definition; the message names the gene
            and lists the supported genes
    """
    bands = _phenotype_bands(_canonical_gene(gene))

    result: list[dict[str, Any]] = []
    lower, lower_inclusive = 0.0, True
    for upper, upper_inclusive, phenotype in bands:
        result.append(
            {
                "lower_bound": lower,
                "upper_bound": upper,
                "lower_inclusive": lower_inclusive,
                "upper_inclusive": upper_inclusive,
                "phenotype": phenotype,
                "abbreviation": phenotype.abbreviation,
                "gene_specific": True,
            }
        )
        lower, lower_inclusive = upper, not upper_inclusive

    return result


def classify_phenotype(
    diplotype: Diplotype | str,
    gene: str,
) -> dict[str, Any]:
    """Full diplotype-to-phenotype classification pipeline.

    Takes a diplotype (object or string like "*1/*4") and gene, computes the
    activity score, and maps it to a metabolizer phenotype. Genes CPIC defines by
    allele function (CYP2C19, SLCO1B1, UGT1A1) are classified from the function of the
    two alleles, the rest from the activity score.

    A ``Diplotype`` whose ``activity_score`` is ``None`` (not scored) is scored here from
    its alleles; it is never read as a score of 0.0. If an allele has no activity value in
    the built-in table (an uncharacterized allele), no score or phenotype can be assigned:
    the result is indeterminate, with ``activity_score`` ``None`` and the allele named, as
    CPIC calls such a diplotype indeterminate. It is never scored as normal function. A
    ``Diplotype`` that carries a score (for example one built with a custom
    ``scoring_table``) is classified from that score.

    Args:
        diplotype: Diplotype object or string in format "*X/*Y"
        gene: Gene symbol

    Returns:
        Dictionary with:
            - "gene": Gene symbol
            - "diplotype": Diplotype string
            - "activity_score": Computed activity score, or None when it cannot be computed
            - "phenotype": MetabolizerPhenotype value
            - "phenotype_abbreviation": Short form (PM, IM, NM, RM, UM, IND)
            - "clinical_significance": Brief clinical note
            - "unknown_alleles": Alleles with no activity value (empty unless the phenotype is
              indeterminate for that reason)

    Raises:
        ValueError: If the gene has no CPIC phenotype definition (the message names the gene and
            lists the supported genes; CYP1A2 can be scored but not classified), the string is not
            ``*X/*Y``, or a ``Diplotype`` object belongs to a different gene.
    """
    gene_upper = _canonical_gene(gene)
    _phenotype_bands(gene_upper)

    if isinstance(diplotype, str):
        parts = diplotype.strip().split("/")
        if len(parts) != 2:
            raise ValueError(f"Invalid diplotype string format: '{diplotype}'. Expected '*X/*Y'.")
        allele1 = parts[0].strip()
        allele2 = parts[1].strip()
        dip = determine_diplotype(allele1, allele2, gene_upper)
    else:
        dip = diplotype
        if dip.gene != gene_upper:
            raise ValueError(f"Diplotype {dip.diplotype_string} belongs to gene {dip.gene}, not {gene_upper}")

    activity_score, phenotype, unknown_alleles = _assess_diplotype(dip, gene_upper)

    # Generate clinical significance note
    clinical_notes = _generate_clinical_note(phenotype, gene_upper, unknown_alleles)

    result = {
        "gene": gene_upper,
        "diplotype": dip.diplotype_string,
        "activity_score": activity_score,
        "phenotype": phenotype,
        "phenotype_abbreviation": phenotype.abbreviation,
        "clinical_significance": clinical_notes,
        "unknown_alleles": unknown_alleles,
    }

    logger.info(
        "Classified %s %s -> AS=%s -> %s (%s)",
        gene_upper,
        dip.diplotype_string,
        "none" if activity_score is None else f"{activity_score:.2f}",
        phenotype.value,
        phenotype.abbreviation,
    )

    return result


def _assess_diplotype(dip: Diplotype, gene: str) -> tuple[float | None, MetabolizerPhenotype, list[str]]:
    """Activity score, phenotype and unscorable alleles of a diplotype of a classifiable gene.

    An unscored diplotype (``activity_score`` None) is scored from its alleles' values in the
    built-in table; with an allele missing from the table it has no score and is indeterminate.
    A scored diplotype uses the two alleles' own activity values (which the genes CPIC defines
    by allele function need) when they reproduce its stored score. One scored with a custom
    scoring table, or containing alleles missing from the built-in table, is classified from its
    stored score alone: the score was supplied, so the alleles are not treated as uncharacterized.

    Returns:
        ``(activity_score, phenotype, unknown_alleles)``
    """
    table = _ACTIVITY_SCORE_TABLES.get(gene, {})
    stored = dip.activity_score
    unknown = uncharacterized_alleles(dip.allele1, dip.allele2, gene)
    name1 = _normalize_allele_name(dip.allele1, gene, table)
    name2 = _normalize_allele_name(dip.allele2, gene, table)

    if stored is None:
        if unknown:
            return None, MetabolizerPhenotype.INDETERMINATE, unknown
        activity1, activity2 = table[name1], table[name2]
        return activity1 + activity2, _phenotype_for_alleles(gene, activity1, activity2), []

    if not unknown:
        activity1, activity2 = table[name1], table[name2]
        if math.isclose(activity1 + activity2, stored, abs_tol=1e-9):
            return stored, _phenotype_for_alleles(gene, activity1, activity2), []
    return stored, predict_metabolizer_status(stored, gene), []


def _generate_clinical_note(
    phenotype: MetabolizerPhenotype,
    gene: str,
    unknown_alleles: Sequence[str] = (),
) -> str:
    """Generate a brief clinical significance note for a phenotype.

    Args:
        phenotype: Metabolizer phenotype
        gene: Gene symbol
        unknown_alleles: Alleles with no activity value, named in an indeterminate note

    Returns:
        Clinical significance string
    """
    if phenotype == MetabolizerPhenotype.INDETERMINATE and unknown_alleles:
        return (
            f"Indeterminate metabolizer status for {gene}. Allele(s) {', '.join(unknown_alleles)} have no "
            "activity value in the CPIC allele table (uncharacterized), so no activity score or phenotype "
            "can be assigned. CPIC calls a diplotype with an uncharacterized allele indeterminate."
        )
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
    Each diplotype is assigned its phenotype exactly as :func:`classify_phenotype`
    does, so a diplotype with an allele that has no activity value (an uncharacterized
    allele) counts as indeterminate ("IND"), never as normal.

    Args:
        allele_frequencies: Mapping of allele name -> population frequency.
            Frequencies should sum to approximately 1.0. Names that normalize to the same
            allele ("4" and "*4") are added together.
        gene: Gene symbol

    Returns:
        Dictionary mapping phenotype abbreviation -> expected frequency, without the
        categories that have none ("IND" appears only when an uncharacterized allele has
        a nonzero frequency).
        Example: {"PM": 0.05, "IM": 0.20, "NM": 0.65, "RM": 0.08, "UM": 0.02}

    Raises:
        ValueError: If the gene has no CPIC phenotype definition (the message names the gene
            and lists the supported genes), a frequency is negative or NaN, or no frequency
            is positive.
    """
    gene_upper = _canonical_gene(gene)
    _phenotype_bands(gene_upper)

    scoring_table = _ACTIVITY_SCORE_TABLES.get(gene_upper, {})

    # Normalize allele names (non-star designations already in the table are kept)
    alleles: dict[str, float] = {}
    for name, freq in allele_frequencies.items():
        if math.isnan(freq) or freq < 0:
            raise ValueError(f"Allele frequency of {name!r} must be a non-negative number, got {freq!r}")
        key = _normalize_allele_name(name, gene_upper, scoring_table)
        alleles[key] = alleles.get(key, 0.0) + freq

    # Verify frequencies sum to ~1.0
    total_freq = sum(alleles.values())
    if total_freq <= 0:
        raise ValueError("allele_frequencies must contain at least one positive frequency")
    if abs(total_freq - 1.0) > 0.05:
        logger.warning(
            "Allele frequencies for %s sum to %.4f (expected ~1.0). Normalizing.",
            gene_upper,
            total_freq,
        )
        alleles = {k: v / total_freq for k, v in alleles.items()}

    # Calculate diplotype frequencies under HWE
    allele_names = sorted(alleles.keys())
    uncharacterized = [name for name in allele_names if name not in scoring_table]
    if uncharacterized:
        logger.warning(
            "Allele(s) %s of %s have no activity value; diplotypes carrying them are counted as indeterminate",
            ", ".join(uncharacterized),
            gene_upper,
        )
    phenotype_freqs: dict[str, float] = {p.abbreviation: 0.0 for p in MetabolizerPhenotype}

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

            # Get the allele activity values for this diplotype; an allele with none makes the
            # diplotype indeterminate (it is not assumed to be normal function)
            score1 = scoring_table.get(a1)
            score2 = scoring_table.get(a2)

            if score1 is None or score2 is None:
                phenotype = MetabolizerPhenotype.INDETERMINATE
            else:
                phenotype = _phenotype_for_alleles(gene_upper, score1, score2)

            phenotype_freqs[phenotype.abbreviation] += diplotype_freq

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
