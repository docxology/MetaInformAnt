"""Metabolic pathway reconstruction from metagenomic annotations.

Implements pathway mapping, completeness scoring, and cross-sample pathway
comparison. Maps functional annotations to KEGG and MetaCyc pathways to
reconstruct the metabolic potential of metagenomic communities.

Pathway reconstruction approach:
1. Map annotated genes (KOs, ECs) to pathway definitions.
2. Calculate pathway completeness as fraction of required reactions present.
3. Score confidence using coverage redundancy and annotation quality.
4. Compare pathway profiles across samples for differential metabolism.
"""

from __future__ import annotations

import os
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

_ENV_PREFIX = "META_"


@dataclass
class PathwayDefinition:
    """Definition of a metabolic pathway."""

    pathway_id: str
    name: str
    description: str = ""
    database: str = "KEGG"  # 'KEGG' or 'MetaCyc'
    required_kos: list[str] = field(default_factory=list)
    required_ecs: list[str] = field(default_factory=list)
    optional_kos: list[str] = field(default_factory=list)
    # Reactions organized as steps (list of alternatives per step)
    reaction_steps: list[list[str]] = field(default_factory=list)
    category: str = ""
    subcategory: str = ""

    @property
    def total_required(self) -> int:
        """Total number of required reactions/genes."""
        if self.reaction_steps:
            return len(self.reaction_steps)
        return len(self.required_kos) + len(self.required_ecs)


@dataclass
class PathwayResult:
    """Result of pathway reconstruction for one sample."""

    pathway_id: str
    pathway_name: str
    completeness: float  # 0.0 - 1.0
    matched_kos: list[str] = field(default_factory=list)
    matched_ecs: list[str] = field(default_factory=list)
    missing_kos: list[str] = field(default_factory=list)
    missing_ecs: list[str] = field(default_factory=list)
    total_required: int = 0
    total_matched: int = 0
    redundancy: float = 1.0  # Average copies per matched gene
    confidence: float = 0.0


# Built-in reference pathways (subset of KEGG)
_BUILTIN_PATHWAYS: dict[str, PathwayDefinition] = {
    "map00010": PathwayDefinition(
        pathway_id="map00010",
        name="Glycolysis / Gluconeogenesis",
        database="KEGG",
        required_kos=[
            "K00844",
            "K00845",
            "K01810",
            "K00850",
            "K01623",
            "K01624",
            "K01803",
            "K00134",
            "K00927",
            "K01689",
            "K00873",
        ],
        category="Metabolism",
        subcategory="Carbohydrate metabolism",
    ),
    "map00020": PathwayDefinition(
        pathway_id="map00020",
        name="Citrate cycle (TCA cycle)",
        database="KEGG",
        required_kos=[
            "K01647",
            "K01681",
            "K00031",
            "K00030",
            "K00164",
            "K00658",
            "K01902",
            "K01903",
            "K00234",
            "K00235",
            "K01676",
            "K00024",
        ],
        category="Metabolism",
        subcategory="Carbohydrate metabolism",
    ),
    "map00030": PathwayDefinition(
        pathway_id="map00030",
        name="Pentose phosphate pathway",
        database="KEGG",
        required_kos=[
            "K00036",
            "K01057",
            "K00033",
            "K01783",
            "K00615",
            "K00616",
            "K01807",
        ],
        category="Metabolism",
        subcategory="Carbohydrate metabolism",
    ),
    "map00190": PathwayDefinition(
        pathway_id="map00190",
        name="Oxidative phosphorylation",
        database="KEGG",
        required_kos=[
            "K02274",
            "K02275",
            "K02276",
            "K02277",
            "K00330",
            "K00331",
            "K00332",
            "K00333",
            "K00411",
            "K00412",
            "K02256",
            "K02257",
            "K02258",
            "K02259",
        ],
        category="Metabolism",
        subcategory="Energy metabolism",
    ),
    "map00220": PathwayDefinition(
        pathway_id="map00220",
        name="Arginine biosynthesis",
        database="KEGG",
        required_kos=[
            "K01940",
            "K01755",
            "K00611",
            "K01438",
            "K01476",
            "K00145",
            "K00930",
            "K00818",
        ],
        category="Metabolism",
        subcategory="Amino acid metabolism",
    ),
    "map00910": PathwayDefinition(
        pathway_id="map00910",
        name="Nitrogen metabolism",
        database="KEGG",
        required_kos=[
            "K02567",
            "K02568",
            "K00360",
            "K00366",
            "K00362",
            "K00363",
            "K10944",
            "K10945",
            "K10946",
        ],
        category="Metabolism",
        subcategory="Energy metabolism",
    ),
    "map00680": PathwayDefinition(
        pathway_id="map00680",
        name="Methane metabolism",
        database="KEGG",
        required_kos=[
            "K00399",
            "K00401",
            "K00402",
            "K00169",
            "K00170",
            "K00171",
            "K00172",
            "K01895",
        ],
        category="Metabolism",
        subcategory="Energy metabolism",
    ),
    "map00195": PathwayDefinition(
        pathway_id="map00195",
        name="Photosynthesis",
        database="KEGG",
        required_kos=[
            "K02689",
            "K02690",
            "K02691",
            "K02692",
            "K02703",
            "K02706",
            "K02707",
            "K02708",
        ],
        category="Metabolism",
        subcategory="Energy metabolism",
    ),
}


def reconstruct_pathways(
    annotations: dict[str, list[str]],
    database: str = "KEGG",
    pathway_definitions: dict[str, PathwayDefinition] | None = None,
    min_completeness: float = 0.0,
) -> list[PathwayResult]:
    """Reconstruct metabolic pathways from functional annotations.

    Maps annotated KO/EC identifiers to known pathways and calculates
    the completeness (fraction of required reactions present) for each.

    Args:
        annotations: Dict mapping gene IDs to lists of KO/EC annotations.
            Example: {"gene1": ["K00844", "K00845"], "gene2": ["K01810"]}
        database: Pathway database to use ('KEGG' or 'MetaCyc').
        pathway_definitions: Custom pathway definitions. If None, uses built-in set.
        min_completeness: Minimum completeness threshold to include in results.

    Returns:
        List of PathwayResult objects, sorted by completeness descending.
    """
    if pathway_definitions is None:
        pathways = {pid: pdef for pid, pdef in _BUILTIN_PATHWAYS.items() if pdef.database == database}
    else:
        pathways = pathway_definitions

    # Collect all annotated KOs/ECs
    all_kos: set[str] = set()
    all_ecs: set[str] = set()
    ko_counts: dict[str, int] = defaultdict(int)
    ec_counts: dict[str, int] = defaultdict(int)

    for gene_id, gene_annotations in annotations.items():
        for ann in gene_annotations:
            ann_upper = ann.upper().strip()
            if ann_upper.startswith("K") and len(ann_upper) == 6 and ann_upper[1:].isdigit():
                all_kos.add(ann_upper)
                ko_counts[ann_upper] += 1
            elif "." in ann_upper:
                # Likely EC number
                all_ecs.add(ann_upper)
                ec_counts[ann_upper] += 1

    logger.info(
        f"Input annotations: {len(all_kos)} unique KOs, {len(all_ecs)} unique ECs from {len(annotations)} genes"
    )

    results: list[PathwayResult] = []

    for pathway_id, pathway_def in pathways.items():
        matched_kos = [ko for ko in pathway_def.required_kos if ko in all_kos]
        missing_kos = [ko for ko in pathway_def.required_kos if ko not in all_kos]
        matched_ecs = [ec for ec in pathway_def.required_ecs if ec in all_ecs]
        missing_ecs = [ec for ec in pathway_def.required_ecs if ec not in all_ecs]

        total_required = pathway_def.total_required
        total_matched = len(matched_kos) + len(matched_ecs)

        completeness = total_matched / total_required if total_required > 0 else 0.0

        # Calculate redundancy (average copy number of matched genes)
        if matched_kos:
            redundancy = sum(ko_counts[ko] for ko in matched_kos) / len(matched_kos)
        else:
            redundancy = 0.0

        # Confidence score combines completeness and redundancy
        confidence = completeness * min(1.0, redundancy)

        if completeness >= min_completeness:
            results.append(
                PathwayResult(
                    pathway_id=pathway_id,
                    pathway_name=pathway_def.name,
                    completeness=completeness,
                    matched_kos=matched_kos,
                    matched_ecs=matched_ecs,
                    missing_kos=missing_kos,
                    missing_ecs=missing_ecs,
                    total_required=total_required,
                    total_matched=total_matched,
                    redundancy=redundancy,
                    confidence=confidence,
                )
            )

    results.sort(key=lambda r: r.completeness, reverse=True)
    logger.info(
        f"Reconstructed {len(results)} pathways; "
        f"{sum(1 for r in results if r.completeness >= 0.5)} with >=50% completeness"
    )
    return results


def calculate_pathway_completeness(
    pathway: PathwayDefinition,
    annotations: set[str],
) -> float:
    """Calculate completeness of a single pathway.

    Args:
        pathway: Pathway definition with required KOs/ECs.
        annotations: Set of observed KO/EC identifiers.

    Returns:
        Completeness fraction (0.0 to 1.0).
    """
    if pathway.reaction_steps:
        # Step-based: each step is complete if any alternative is present
        steps_complete = 0
        for step_alternatives in pathway.reaction_steps:
            if any(alt in annotations for alt in step_alternatives):
                steps_complete += 1
        return steps_complete / len(pathway.reaction_steps) if pathway.reaction_steps else 0.0

    # Simple KO/EC counting
    total = len(pathway.required_kos) + len(pathway.required_ecs)
    if total == 0:
        return 0.0

    matched = sum(1 for ko in pathway.required_kos if ko in annotations)
    matched += sum(1 for ec in pathway.required_ecs if ec in annotations)

    return matched / total


def compare_pathway_profiles(
    samples: dict[str, dict[str, list[str]]],
    pathway_definitions: dict[str, PathwayDefinition] | None = None,
    database: str = "KEGG",
) -> dict[str, dict[str, float]]:
    """Compare pathway completeness profiles across multiple samples.

    Reconstructs pathways for each sample and builds a comparison matrix
    of completeness scores.

    Args:
        samples: Dict mapping sample IDs to their annotations
            (each sample is {gene_id: [annotations]}).
        pathway_definitions: Custom pathway definitions.
        database: Pathway database to use.

    Returns:
        Dict mapping pathway IDs to {sample_id: completeness} dicts.
    """
    comparison: dict[str, dict[str, float]] = defaultdict(dict)

    for sample_id, sample_annotations in samples.items():
        results = reconstruct_pathways(
            sample_annotations,
            database=database,
            pathway_definitions=pathway_definitions,
            min_completeness=0.0,
        )

        for result in results:
            comparison[result.pathway_id][sample_id] = result.completeness

    logger.info(f"Compared {len(comparison)} pathways across {len(samples)} samples")
    return dict(comparison)


def find_differential_pathways(
    group1_samples: dict[str, dict[str, list[str]]],
    group2_samples: dict[str, dict[str, list[str]]],
    pathway_definitions: dict[str, PathwayDefinition] | None = None,
    min_diff: float = 0.2,
) -> list[dict[str, Any]]:
    """Find pathways differentially present between two sample groups.

    Compares mean pathway completeness between groups and identifies
    pathways with significant differences.

    Args:
        group1_samples: Annotations for group 1 samples.
        group2_samples: Annotations for group 2 samples.
        pathway_definitions: Custom pathway definitions.
        min_diff: Minimum difference in mean completeness to report.

    Returns:
        List of dicts with pathway_id, pathway_name, group1_mean,
        group2_mean, difference, and direction.
    """
    all_samples = {
        **{f"g1_{k}": v for k, v in group1_samples.items()},
        **{f"g2_{k}": v for k, v in group2_samples.items()},
    }

    profiles = compare_pathway_profiles(all_samples, pathway_definitions=pathway_definitions)

    g1_ids = {f"g1_{k}" for k in group1_samples}
    g2_ids = {f"g2_{k}" for k in group2_samples}

    differential: list[dict[str, Any]] = []

    for pathway_id, sample_scores in profiles.items():
        g1_scores = [score for sid, score in sample_scores.items() if sid in g1_ids]
        g2_scores = [score for sid, score in sample_scores.items() if sid in g2_ids]

        g1_mean = sum(g1_scores) / len(g1_scores) if g1_scores else 0.0
        g2_mean = sum(g2_scores) / len(g2_scores) if g2_scores else 0.0
        diff = g1_mean - g2_mean

        if abs(diff) >= min_diff:
            # Get pathway name from definitions or builtins
            pathway_name = pathway_id
            if pathway_definitions and pathway_id in pathway_definitions:
                pathway_name = pathway_definitions[pathway_id].name
            elif pathway_id in _BUILTIN_PATHWAYS:
                pathway_name = _BUILTIN_PATHWAYS[pathway_id].name

            differential.append(
                {
                    "pathway_id": pathway_id,
                    "pathway_name": pathway_name,
                    "group1_mean": round(g1_mean, 4),
                    "group2_mean": round(g2_mean, 4),
                    "difference": round(diff, 4),
                    "abs_difference": round(abs(diff), 4),
                    "direction": "group1_enriched" if diff > 0 else "group2_enriched",
                }
            )

    differential.sort(key=lambda d: d["abs_difference"], reverse=True)
    logger.info(f"Found {len(differential)} differentially present pathways (min_diff={min_diff})")
    return differential


__all__ = [
    "PathwayDefinition",
    "PathwayResult",
    "reconstruct_pathways",
    "calculate_pathway_completeness",
    "compare_pathway_profiles",
    "find_differential_pathways",
]
