"""Functional impact prediction for structural variants.

Predicts the functional consequences of structural variants by assessing
gene dosage sensitivity, topologically associating domain (TAD) boundary
disruption, and aggregate pathogenicity scoring.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Sequence

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


@dataclass
class FunctionalImpact:
    """Functional impact assessment for a structural variant.

    Attributes:
        variant_id: Identifier for the variant.
        impact_level: Impact level ('HIGH', 'MODERATE', 'LOW', 'MODIFIER').
        impact_type: Specific impact type (e.g., 'gene_disruption', 'gene_fusion',
            'dosage_change', 'regulatory_disruption', 'tad_disruption').
        affected_genes: List of affected gene names.
        dosage_sensitive_genes: Genes with dosage sensitivity concern.
        tad_disrupted: Whether a TAD boundary is disrupted.
        pathogenicity_score: Aggregate pathogenicity score [0, 1].
        details: Additional details about the impact.
    """

    variant_id: str = ""
    impact_level: str = "MODIFIER"
    impact_type: str = ""
    affected_genes: list[str] = field(default_factory=list)
    dosage_sensitive_genes: list[str] = field(default_factory=list)
    tad_disrupted: bool = False
    pathogenicity_score: float = 0.0
    details: dict[str, Any] = field(default_factory=dict)


@dataclass
class DosageSensitivity:
    """Dosage sensitivity information for a gene.

    Attributes:
        gene_name: Gene name.
        haploinsufficiency_score: Haploinsufficiency score [0, 1].
            Higher means more likely to be haploinsufficient.
        triplosensitivity_score: Triplosensitivity score [0, 1].
            Higher means more sensitive to copy gain.
        pli_score: Probability of loss-of-function intolerance (pLI) [0, 1].
        loeuf: Loss-of-function observed/expected upper fraction (LOEUF).
            Lower values indicate more constraint.
        is_haploinsufficient: Whether gene is classified as haploinsufficient.
        is_triplosensitive: Whether gene is classified as triplosensitive.
    """

    gene_name: str
    haploinsufficiency_score: float = 0.0
    triplosensitivity_score: float = 0.0
    pli_score: float = 0.0
    loeuf: float = 1.0
    is_haploinsufficient: bool = False
    is_triplosensitive: bool = False


@dataclass
class TADPrediction:
    """TAD boundary disruption prediction.

    Attributes:
        variant_chrom: Chromosome.
        variant_start: Variant start position.
        variant_end: Variant end position.
        disrupted_boundaries: List of disrupted TAD boundary positions.
        n_boundaries_disrupted: Number of TAD boundaries disrupted.
        genes_in_affected_tads: Genes within affected TAD domains.
        disruption_score: TAD disruption severity score [0, 1].
    """

    variant_chrom: str = ""
    variant_start: int = 0
    variant_end: int = 0
    disrupted_boundaries: list[int] = field(default_factory=list)
    n_boundaries_disrupted: int = 0
    genes_in_affected_tads: list[str] = field(default_factory=list)
    disruption_score: float = 0.0


def predict_functional_impact(
    variant: dict[str, Any],
    gene_annotations: list[dict[str, Any]],
    haploinsufficiency_db: dict[str, float] | None = None,
    tad_boundaries: list[dict[str, Any]] | None = None,
) -> FunctionalImpact:
    """Predict the functional impact of a structural variant.

    Integrates multiple lines of evidence to classify the functional
    consequence of an SV:

    1. Gene overlap analysis (disrupted vs. duplicated genes)
    2. Dosage sensitivity of affected genes
    3. TAD boundary disruption
    4. Variant type-specific considerations (DEL vs DUP vs INV)
    5. Aggregate pathogenicity scoring

    Args:
        variant: Variant dictionary containing:
            - 'chrom': Chromosome
            - 'start': Start position
            - 'end': End position
            - 'sv_type': Type of SV (DEL, DUP, INV, TRA, INS)
            - 'id': Optional variant identifier
            - 'overlapping_genes': Optional pre-computed gene overlaps
        gene_annotations: List of gene annotation dictionaries with:
            - 'chrom': Chromosome
            - 'start': Gene start
            - 'end': Gene end
            - 'name': Gene name
            - 'is_coding': Whether gene is protein-coding (optional)
        haploinsufficiency_db: Optional dictionary mapping gene names to
            haploinsufficiency scores (0-1).
        tad_boundaries: Optional list of TAD boundary dictionaries with
            'chrom', 'start', 'end' keys.

    Returns:
        FunctionalImpact object with comprehensive impact assessment.
    """
    sv_type = variant.get("sv_type", "UNKNOWN")
    if hasattr(sv_type, "value"):
        sv_type = sv_type.value

    chrom = variant.get("chrom", "")
    start = variant.get("start", 0)
    end = variant.get("end", 0)
    variant_id = variant.get("id", f"{chrom}:{start}-{end}")
    size = abs(end - start)

    # Find overlapping genes
    overlapping_genes = variant.get("overlapping_genes", [])
    if not overlapping_genes:
        overlapping_genes = _find_overlapping_genes(chrom, start, end, gene_annotations)

    # Assess dosage sensitivity
    dosage_genes: list[str] = []
    if haploinsufficiency_db:
        for gene_name in overlapping_genes:
            sensitivity = assess_dosage_sensitivity(gene_name, haploinsufficiency_db)
            if sensitivity.is_haploinsufficient and sv_type in ("DEL", "HOMODEL"):
                dosage_genes.append(gene_name)
            elif sensitivity.is_triplosensitive and sv_type == "DUP":
                dosage_genes.append(gene_name)

    # Assess TAD disruption
    tad_disrupted = False
    tad_pred: TADPrediction | None = None
    if tad_boundaries:
        tad_pred = predict_tad_disruption(variant, tad_boundaries)
        tad_disrupted = tad_pred.n_boundaries_disrupted > 0

    # Determine impact level and type
    impact_level, impact_type = _classify_impact(
        sv_type=sv_type,
        size=size,
        n_genes=len(overlapping_genes),
        n_dosage_genes=len(dosage_genes),
        tad_disrupted=tad_disrupted,
        gene_annotations=gene_annotations,
        overlapping_genes=overlapping_genes,
    )

    # Compute pathogenicity score
    path_score = score_pathogenicity(
        variant,
        {
            "overlapping_genes": overlapping_genes,
            "dosage_sensitive_genes": dosage_genes,
            "tad_disrupted": tad_disrupted,
            "impact_level": impact_level,
            "n_coding_genes": _count_coding_genes(overlapping_genes, gene_annotations),
        },
    )

    details: dict[str, Any] = {
        "sv_type": sv_type,
        "size": size,
        "n_overlapping_genes": len(overlapping_genes),
    }
    if tad_pred:
        details["tad_prediction"] = {
            "n_boundaries": tad_pred.n_boundaries_disrupted,
            "disruption_score": tad_pred.disruption_score,
        }

    return FunctionalImpact(
        variant_id=variant_id,
        impact_level=impact_level,
        impact_type=impact_type,
        affected_genes=overlapping_genes,
        dosage_sensitive_genes=dosage_genes,
        tad_disrupted=tad_disrupted,
        pathogenicity_score=path_score,
        details=details,
    )


def assess_dosage_sensitivity(
    gene: str,
    haploinsufficiency_db: dict[str, float] | None = None,
    triplosensitivity_db: dict[str, float] | None = None,
    pli_db: dict[str, float] | None = None,
    loeuf_db: dict[str, float] | None = None,
) -> DosageSensitivity:
    """Assess dosage sensitivity of a gene.

    Looks up the gene in multiple dosage sensitivity databases and
    produces a composite assessment. A gene is classified as
    haploinsufficient if its HI score exceeds 0.5 or pLI exceeds 0.9.

    Args:
        gene: Gene name to look up.
        haploinsufficiency_db: Dictionary mapping gene names to HI scores (0-1).
        triplosensitivity_db: Dictionary mapping gene names to TS scores (0-1).
        pli_db: Dictionary mapping gene names to pLI scores (0-1).
        loeuf_db: Dictionary mapping gene names to LOEUF values.

    Returns:
        DosageSensitivity object with all available information.
    """
    hi_score = 0.0
    ts_score = 0.0
    pli_score = 0.0
    loeuf = 1.0

    if haploinsufficiency_db:
        hi_score = haploinsufficiency_db.get(gene, 0.0)

    if triplosensitivity_db:
        ts_score = triplosensitivity_db.get(gene, 0.0)

    if pli_db:
        pli_score = pli_db.get(gene, 0.0)

    if loeuf_db:
        loeuf = loeuf_db.get(gene, 1.0)

    # Classification thresholds
    is_hi = hi_score >= 0.5 or pli_score >= 0.9 or loeuf < 0.35
    is_ts = ts_score >= 0.5

    return DosageSensitivity(
        gene_name=gene,
        haploinsufficiency_score=hi_score,
        triplosensitivity_score=ts_score,
        pli_score=pli_score,
        loeuf=loeuf,
        is_haploinsufficient=is_hi,
        is_triplosensitive=is_ts,
    )


def predict_tad_disruption(
    variant: dict[str, Any],
    tad_boundaries: list[dict[str, Any]],
    boundary_window: int = 10_000,
) -> TADPrediction:
    """Predict whether a structural variant disrupts TAD boundaries.

    TAD (Topologically Associating Domain) boundaries are important
    chromatin organization elements. Their disruption can lead to
    ectopic enhancer-promoter contacts and gene misregulation.

    A boundary is considered disrupted if the SV breakpoint falls
    within boundary_window base pairs of the boundary midpoint.

    Args:
        variant: Variant dictionary with 'chrom', 'start', 'end', 'sv_type'.
        tad_boundaries: List of TAD boundary dictionaries with:
            - 'chrom': Chromosome
            - 'start': Boundary start position
            - 'end': Boundary end position
            - 'insulation_score': Optional insulation score
            - 'genes': Optional list of genes in adjacent TADs
        boundary_window: Window size around each boundary for disruption
            detection (default 10kb).

    Returns:
        TADPrediction object with disruption assessment.
    """
    chrom = variant.get("chrom", "")
    start = variant.get("start", 0)
    end = variant.get("end", 0)
    sv_type = variant.get("sv_type", "UNKNOWN")
    if hasattr(sv_type, "value"):
        sv_type = sv_type.value

    disrupted: list[int] = []
    affected_genes: list[str] = []

    for boundary in tad_boundaries:
        b_chrom = boundary.get("chrom", "")
        if b_chrom != chrom:
            continue

        b_start = boundary.get("start", 0)
        b_end = boundary.get("end", 0)
        b_mid = (b_start + b_end) // 2

        # Check if variant spans or is near the boundary
        # Disruption can occur from:
        # 1. Deletion spanning the boundary
        # 2. Inversion with breakpoint near boundary
        # 3. Translocation with breakpoint near boundary

        is_disrupted = False

        if sv_type == "DEL":
            # Deletion spanning the boundary
            if start < b_mid < end:
                is_disrupted = True
            # Deletion breakpoint near boundary
            elif abs(start - b_mid) <= boundary_window or abs(end - b_mid) <= boundary_window:
                is_disrupted = True
        elif sv_type in ("INV", "TRA", "BND"):
            # Breakpoint near boundary
            if abs(start - b_mid) <= boundary_window or abs(end - b_mid) <= boundary_window:
                is_disrupted = True
        elif sv_type == "DUP":
            # Duplication spanning the boundary (tandem dup inserts copy at boundary)
            if start < b_mid < end:
                is_disrupted = True
        else:
            # Generic check
            if start - boundary_window < b_mid < end + boundary_window:
                is_disrupted = True

        if is_disrupted:
            disrupted.append(b_mid)
            # Collect genes from boundary annotation
            b_genes = boundary.get("genes", [])
            affected_genes.extend(b_genes)

    n_disrupted = len(disrupted)

    # Calculate disruption score
    if n_disrupted == 0:
        disruption_score = 0.0
    else:
        # Score based on number of boundaries and SV type severity
        type_weight = {"DEL": 1.0, "INV": 0.9, "TRA": 0.95, "BND": 0.95, "DUP": 0.7, "INS": 0.3}
        weight = type_weight.get(sv_type, 0.5)
        disruption_score = min(1.0, n_disrupted * 0.4 * weight)

    return TADPrediction(
        variant_chrom=chrom,
        variant_start=start,
        variant_end=end,
        disrupted_boundaries=disrupted,
        n_boundaries_disrupted=n_disrupted,
        genes_in_affected_tads=list(set(affected_genes)),
        disruption_score=disruption_score,
    )


def score_pathogenicity(
    variant: dict[str, Any],
    annotations: dict[str, Any],
) -> float:
    """Compute aggregate pathogenicity score for a structural variant.

    Combines multiple evidence lines into a single [0, 1] score using
    a weighted logistic model:

    - Size of the variant (larger = more likely pathogenic)
    - Number of affected coding genes
    - Dosage sensitivity of affected genes
    - TAD boundary disruption
    - Variant type (DEL/DUP more pathogenic than INV on average)
    - Population frequency (if available)

    Args:
        variant: Variant dictionary with 'chrom', 'start', 'end', 'sv_type'.
        annotations: Annotation dictionary containing:
            - 'overlapping_genes': List of overlapping gene names
            - 'dosage_sensitive_genes': List of dosage-sensitive genes
            - 'tad_disrupted': Whether TAD boundary is disrupted
            - 'impact_level': Impact level string
            - 'n_coding_genes': Number of affected coding genes
            - 'population_frequency': Allele frequency in population (optional)

    Returns:
        Pathogenicity score between 0 and 1.
    """
    sv_type = variant.get("sv_type", "UNKNOWN")
    if hasattr(sv_type, "value"):
        sv_type = sv_type.value

    start = variant.get("start", 0)
    end = variant.get("end", 0)
    size = abs(end - start)

    n_genes = len(annotations.get("overlapping_genes", []))
    n_dosage = len(annotations.get("dosage_sensitive_genes", []))
    tad_disrupted = annotations.get("tad_disrupted", False)
    impact_level = annotations.get("impact_level", "MODIFIER")
    n_coding = annotations.get("n_coding_genes", 0)
    pop_freq = annotations.get("population_frequency", None)

    # Component 1: Size (sigmoid, midpoint at 100kb)
    size_score = 1.0 / (1.0 + math.exp(-0.00005 * (size - 100_000)))

    # Component 2: Gene impact
    gene_score = 0.0
    if n_coding > 0:
        gene_score = min(1.0, n_coding * 0.15)
    if n_dosage > 0:
        gene_score = min(1.0, gene_score + n_dosage * 0.3)

    # Component 3: TAD disruption
    tad_score = 0.6 if tad_disrupted else 0.0

    # Component 4: SV type weight
    type_weights = {
        "DEL": 0.7,
        "HOMODEL": 0.9,
        "DUP": 0.5,
        "AMP": 0.6,
        "INV": 0.4,
        "TRA": 0.8,
        "BND": 0.8,
        "INS": 0.3,
        "UNKNOWN": 0.3,
    }
    type_score = type_weights.get(sv_type, 0.3)

    # Component 5: Impact level
    level_weights = {"HIGH": 0.9, "MODERATE": 0.5, "LOW": 0.2, "MODIFIER": 0.1}
    level_score = level_weights.get(impact_level, 0.1)

    # Component 6: Population frequency (rare = more pathogenic)
    freq_score = 1.0
    if pop_freq is not None and pop_freq > 0:
        # Common variants are less likely pathogenic
        freq_score = max(0.0, 1.0 - math.log10(pop_freq + 1e-6) / (-1.0))
        freq_score = min(1.0, max(0.0, freq_score))

    # Weighted combination
    raw_score = (
        0.15 * size_score
        + 0.30 * gene_score
        + 0.10 * tad_score
        + 0.15 * type_score
        + 0.15 * level_score
        + 0.15 * freq_score
    )

    return max(0.0, min(1.0, raw_score))


def _find_overlapping_genes(
    chrom: str,
    start: int,
    end: int,
    gene_annotations: list[dict[str, Any]],
) -> list[str]:
    """Find gene names overlapping a region.

    Args:
        chrom: Chromosome.
        start: Start position.
        end: End position.
        gene_annotations: List of gene annotation dicts.

    Returns:
        List of overlapping gene names.
    """
    genes: list[str] = []
    for gene in gene_annotations:
        g_chrom = gene.get("chrom", "")
        g_start = gene.get("start", 0)
        g_end = gene.get("end", 0)
        g_name = gene.get("name", "")

        if g_chrom == chrom and g_start < end and g_end > start:
            genes.append(g_name)

    return genes


def _count_coding_genes(
    gene_names: list[str],
    gene_annotations: list[dict[str, Any]],
) -> int:
    """Count the number of protein-coding genes in a list.

    Args:
        gene_names: List of gene names.
        gene_annotations: Gene annotation database.

    Returns:
        Number of protein-coding genes.
    """
    coding_set: set[str] = set()
    for gene in gene_annotations:
        if gene.get("is_coding", True):  # Default to coding if not specified
            coding_set.add(gene.get("name", ""))

    return sum(1 for g in gene_names if g in coding_set)


def _classify_impact(
    sv_type: str,
    size: int,
    n_genes: int,
    n_dosage_genes: int,
    tad_disrupted: bool,
    gene_annotations: list[dict[str, Any]],
    overlapping_genes: list[str],
) -> tuple[str, str]:
    """Classify impact level and type.

    Args:
        sv_type: Type of structural variant.
        size: Size in base pairs.
        n_genes: Number of affected genes.
        n_dosage_genes: Number of dosage-sensitive genes.
        tad_disrupted: Whether TAD boundary is disrupted.
        gene_annotations: Gene annotations for context.
        overlapping_genes: List of overlapping gene names.

    Returns:
        Tuple of (impact_level, impact_type).
    """
    # HIGH impact conditions
    if n_dosage_genes > 0:
        if sv_type in ("DEL", "HOMODEL"):
            return "HIGH", "dosage_loss"
        elif sv_type == "DUP":
            return "HIGH", "dosage_gain"

    if sv_type in ("DEL", "HOMODEL") and n_genes > 0:
        n_coding = _count_coding_genes(overlapping_genes, gene_annotations)
        if n_coding > 0:
            return "HIGH", "gene_disruption"

    if sv_type == "TRA" and n_genes > 0:
        return "HIGH", "gene_fusion_candidate"

    # MODERATE impact conditions
    if n_genes > 0:
        if sv_type == "DUP":
            return "MODERATE", "gene_duplication"
        elif sv_type == "INV":
            return "MODERATE", "gene_inversion"
        elif sv_type == "INS":
            return "MODERATE", "gene_insertion"

    if tad_disrupted:
        return "MODERATE", "tad_disruption"

    # LOW impact
    if size > 100_000:
        return "LOW", "large_intergenic"

    if n_genes == 0 and size > 10_000:
        return "LOW", "intergenic"

    return "MODIFIER", "benign"
