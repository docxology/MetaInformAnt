"""Protein function prediction - analysis utilities.

This module provides tools for predicting intrinsically disordered regions,
active sites, and post-translational modification sites. All predictions
use pure Python implementations with biologically grounded algorithms
and pattern matching.

Functions in this module:
    - predict_disordered_regions: Intrinsically disordered region prediction
    - find_active_sites: Active site identification from catalytic motifs
    - predict_post_translational_mods: PTM site prediction from sequence motifs
"""

from __future__ import annotations

import math
import re
from typing import Any, Dict, List, Tuple

from metainformant.core.utils.logging import get_logger
from metainformant.protein.function.prediction_core import _DISORDER_PROPENSITY

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Intrinsically disordered region prediction
# ---------------------------------------------------------------------------


def predict_disordered_regions(
    sequence: str,
    method: str = "iupred_like",
    window: int = 21,
    threshold: float = 0.5,
) -> List[Dict[str, Any]]:
    """Predict intrinsically disordered regions using amino acid propensity.

    Implements a simplified IUPred-like algorithm that scores each residue
    for disorder propensity based on its amino acid type and local sequence
    context. Contiguous regions above the threshold are reported as
    intrinsically disordered regions (IDRs).

    Args:
        sequence: Protein amino acid sequence string.
        method: Prediction method. Currently supports ``"iupred_like"``
            (amino acid propensity with sliding window).
        window: Sliding window size for smoothing (must be odd).
        threshold: Minimum score to classify a region as disordered.

    Returns:
        List of dicts, each with:
            - start: 0-based start position
            - end: 0-based end position (exclusive)
            - length: region length
            - mean_score: average disorder score in the region
            - max_score: maximum disorder score in the region
            - sequence: subsequence of the disordered region

    Raises:
        ValueError: If method is not supported.
    """
    if method not in ("iupred_like",):
        raise ValueError(f"Unsupported disorder prediction method: {method}")

    if not sequence or len(sequence) < window:
        return []

    seq = sequence.upper()
    n = len(seq)

    # Compute raw per-residue disorder propensity
    raw_scores = [_DISORDER_PROPENSITY.get(aa, 0.0) for aa in seq]

    # Smooth with sliding window
    half_w = window // 2
    smoothed: List[float] = []
    for i in range(n):
        start = max(0, i - half_w)
        end = min(n, i + half_w + 1)
        window_vals = raw_scores[start:end]
        smoothed.append(sum(window_vals) / len(window_vals))

    # Normalize to 0-1 range using sigmoid-like transformation
    normalized: List[float] = []
    for s in smoothed:
        # Map propensity scores to probability-like values
        prob = 1.0 / (1.0 + math.exp(-5.0 * s))
        normalized.append(prob)

    # Find contiguous regions above threshold
    regions: List[Dict[str, Any]] = []
    i = 0
    while i < n:
        if normalized[i] >= threshold:
            region_start = i
            region_scores: List[float] = []
            while i < n and normalized[i] >= threshold:
                region_scores.append(normalized[i])
                i += 1
            region_end = i

            # Require minimum length of 5 residues
            if region_end - region_start >= 5:
                regions.append(
                    {
                        "start": region_start,
                        "end": region_end,
                        "length": region_end - region_start,
                        "mean_score": round(sum(region_scores) / len(region_scores), 4),
                        "max_score": round(max(region_scores), 4),
                        "sequence": seq[region_start:region_end],
                    }
                )
        else:
            i += 1

    logger.debug("Predicted %d disordered regions", len(regions))
    return regions


# ---------------------------------------------------------------------------
# Active site prediction
# ---------------------------------------------------------------------------

# Known catalytic motifs within specific domain types
_ACTIVE_SITE_PATTERNS: List[Dict[str, Any]] = [
    {
        "name": "Serine protease catalytic triad",
        "domains": {"Trypsin-like serine protease"},
        "patterns": [
            {"residue": "H", "pattern": r"[GSAC]H[ILVFYWM]", "role": "catalytic histidine"},
            {"residue": "D", "pattern": r"D[ILVFA][ASTG]", "role": "catalytic aspartate"},
            {"residue": "S", "pattern": r"G[DNSG]SG[GAS]", "role": "catalytic serine"},
        ],
    },
    {
        "name": "Protein kinase active site",
        "domains": {"Protein kinase domain"},
        "patterns": [
            {"residue": "K", "pattern": r"V[AI]K", "role": "ATP binding lysine"},
            {"residue": "D", "pattern": r"[HY]RD[ILVMA]", "role": "catalytic aspartate"},
            {"residue": "D", "pattern": r"DFG", "role": "DFG motif aspartate"},
        ],
    },
    {
        "name": "GTPase catalytic site",
        "domains": {"Ras family GTPase"},
        "patterns": [
            {"residue": "T", "pattern": r"DTAGQ", "role": "catalytic threonine"},
            {"residue": "G", "pattern": r"G.{4}GK[ST]", "role": "P-loop glycine"},
        ],
    },
    {
        "name": "Metalloprotease zinc binding",
        "domains": set(),  # generic, no domain restriction
        "patterns": [
            {"residue": "H", "pattern": r"HE..H", "role": "zinc binding histidine"},
            {"residue": "E", "pattern": r"HE..H", "role": "catalytic glutamate"},
        ],
    },
    {
        "name": "Cysteine protease active site",
        "domains": set(),
        "patterns": [
            {"residue": "C", "pattern": r"Q.C[GSAW]", "role": "catalytic cysteine"},
            {"residue": "H", "pattern": r"[ILVFA]H[AG]", "role": "catalytic histidine"},
        ],
    },
]


def find_active_sites(
    sequence: str,
    domains: List[Dict[str, Any]] | None = None,
) -> List[Dict[str, Any]]:
    """Identify potential catalytic residues from conserved motif patterns.

    Scans the protein sequence for known catalytic site motifs. When domain
    annotations are provided, the search is focused to domain-specific
    patterns for higher specificity. Without domain information, all generic
    catalytic motif patterns are tested.

    Args:
        sequence: Protein amino acid sequence string.
        domains: Optional list of domain hit dicts with ``name`` keys.
            Used to restrict pattern search to domain-appropriate motifs.

    Returns:
        List of dicts, each with:
            - site_name: name of the catalytic site type
            - residue: the catalytic amino acid
            - position: 0-based position in the sequence
            - role: functional role of the residue
            - motif_matched: the regex pattern that matched
            - context: surrounding sequence context (5 residues each side)
            - confidence: confidence score (0.0 to 1.0)
    """
    if not sequence:
        return []

    seq = sequence.upper()
    domain_names: set[str] = set()
    if domains:
        domain_names = {d.get("name", "") for d in domains}

    results: List[Dict[str, Any]] = []

    for site_def in _ACTIVE_SITE_PATTERNS:
        # Check if this pattern applies (domain match or generic)
        required_domains = site_def.get("domains", set())
        if required_domains and not (required_domains & domain_names):
            # No matching domain, but still check if no domains provided
            if domains is not None:
                continue

        # Confidence is higher when domain context matches
        base_confidence = 0.8 if (required_domains & domain_names) else 0.4
        if not required_domains:
            base_confidence = 0.5

        for pattern_def in site_def["patterns"]:
            target_residue = pattern_def["residue"]
            pattern = pattern_def["pattern"]
            role = pattern_def["role"]

            try:
                for match in re.finditer(pattern, seq):
                    # Find the target residue within the match
                    matched_str = match.group()
                    for offset, aa in enumerate(matched_str):
                        if aa == target_residue:
                            abs_pos = match.start() + offset
                            # Get context
                            ctx_start = max(0, abs_pos - 5)
                            ctx_end = min(len(seq), abs_pos + 6)
                            context = seq[ctx_start:ctx_end]

                            results.append(
                                {
                                    "site_name": site_def["name"],
                                    "residue": target_residue,
                                    "position": abs_pos,
                                    "role": role,
                                    "motif_matched": pattern,
                                    "context": context,
                                    "confidence": round(base_confidence, 4),
                                }
                            )
                            break  # only first target residue per match

            except re.error as exc:
                logger.warning("Regex error in active site pattern: %s", exc)

    # Deduplicate by position
    seen_positions: set[int] = set()
    deduped: List[Dict[str, Any]] = []
    for result in sorted(results, key=lambda r: (-r["confidence"], r["position"])):
        if result["position"] not in seen_positions:
            deduped.append(result)
            seen_positions.add(result["position"])

    deduped.sort(key=lambda r: r["position"])

    logger.debug("Found %d potential active site residues", len(deduped))
    return deduped


# ---------------------------------------------------------------------------
# Post-translational modification site prediction
# ---------------------------------------------------------------------------

# PTM motif definitions
_PTM_PATTERNS: List[Dict[str, Any]] = [
    # Phosphorylation
    {
        "type": "phosphorylation",
        "subtype": "Ser/Thr kinase (basophilic)",
        "pattern": r"[RK]{2,3}.{0,2}[ST]",
        "target_residue": "ST",
        "description": "Basophilic kinase substrate motif (PKA/PKC-like)",
    },
    {
        "type": "phosphorylation",
        "subtype": "Ser/Thr kinase (proline-directed)",
        "pattern": r"[ST]P",
        "target_residue": "ST",
        "description": "Proline-directed kinase substrate (CDK/MAPK-like)",
    },
    {
        "type": "phosphorylation",
        "subtype": "Tyr kinase",
        "pattern": r"[EDQN].{0,2}Y.{0,2}[ILVFM]",
        "target_residue": "Y",
        "description": "Tyrosine kinase substrate motif",
    },
    {
        "type": "phosphorylation",
        "subtype": "CK2 substrate",
        "pattern": r"[ST].{2}[DE]",
        "target_residue": "ST",
        "description": "Casein kinase 2 substrate motif",
    },
    # N-linked glycosylation
    {
        "type": "glycosylation",
        "subtype": "N-linked",
        "pattern": r"N[^P][ST][^P]",
        "target_residue": "N",
        "description": "N-X-S/T sequon for N-linked glycosylation",
    },
    # O-linked glycosylation (simplified)
    {
        "type": "glycosylation",
        "subtype": "O-linked",
        "pattern": r"[ST]P.P",
        "target_residue": "ST",
        "description": "O-linked glycosylation motif (mucin-type)",
    },
    # Ubiquitination
    {
        "type": "ubiquitination",
        "subtype": "ubiquitin conjugation",
        "pattern": r"[IVLMA].K.[DERED]",
        "target_residue": "K",
        "description": "Ubiquitination site (lysine in degradation context)",
    },
    # SUMOylation
    {
        "type": "sumoylation",
        "subtype": "SUMO conjugation",
        "pattern": r"[IVLMF]K.E",
        "target_residue": "K",
        "description": "SUMOylation consensus motif (psi-K-X-E)",
    },
    # Acetylation
    {
        "type": "acetylation",
        "subtype": "lysine acetylation",
        "pattern": r"[GAS]K[RKPS]",
        "target_residue": "K",
        "description": "Lysine acetylation motif",
    },
    # Myristoylation (N-terminal)
    {
        "type": "myristoylation",
        "subtype": "N-myristoylation",
        "pattern": r"^MG[^EDRKHPFYW].[STAGCN][^P]",
        "target_residue": "G",
        "description": "N-myristoylation motif (position 2 glycine)",
    },
    # GPI anchor signal
    {
        "type": "gpi_anchor",
        "subtype": "GPI attachment",
        "pattern": r"[GANDSTC][GASV][GANDSTC].{4,12}[LIVMFWG]{5,}$",
        "target_residue": "GANDSTC",
        "description": "GPI anchor signal (C-terminal)",
    },
]


def predict_post_translational_mods(sequence: str) -> List[Dict[str, Any]]:
    """Predict post-translational modification sites from sequence motifs.

    Scans for phosphorylation sites (Ser/Thr/Tyr with kinase-specific
    motifs), N-linked glycosylation (N-X-S/T sequon), O-linked
    glycosylation, ubiquitination, SUMOylation, acetylation,
    myristoylation, and GPI anchor signals.

    Args:
        sequence: Protein amino acid sequence string.

    Returns:
        List of dicts sorted by position, each with:
            - type: PTM type (phosphorylation, glycosylation, etc.)
            - subtype: specific PTM subtype
            - position: 0-based position of the modified residue
            - residue: the modified amino acid
            - motif: the matched motif sequence
            - description: human-readable description
            - confidence: prediction confidence (0.0 to 1.0)
    """
    if not sequence:
        return []

    seq = sequence.upper()
    results: List[Dict[str, Any]] = []

    for ptm_def in _PTM_PATTERNS:
        try:
            for match in re.finditer(ptm_def["pattern"], seq):
                matched_str = match.group()
                target_chars = set(ptm_def["target_residue"])

                # Find the target residue within the match
                for offset, aa in enumerate(matched_str):
                    if aa in target_chars:
                        abs_pos = match.start() + offset

                        # Confidence based on motif specificity
                        confidence = 0.6
                        motif_len = len(matched_str)
                        if motif_len >= 5:
                            confidence = 0.75
                        if motif_len >= 7:
                            confidence = 0.85

                        # Surface accessibility heuristic: charged/polar
                        # neighbors increase confidence
                        window_start = max(0, abs_pos - 3)
                        window_end = min(len(seq), abs_pos + 4)
                        local = seq[window_start:window_end]
                        polar_count = sum(1 for c in local if c in "STNQKRDEH")
                        if polar_count >= 3:
                            confidence = min(confidence + 0.1, 1.0)

                        results.append(
                            {
                                "type": ptm_def["type"],
                                "subtype": ptm_def["subtype"],
                                "position": abs_pos,
                                "residue": aa,
                                "motif": matched_str,
                                "description": ptm_def["description"],
                                "confidence": round(confidence, 4),
                            }
                        )
                        break  # only first target per match

        except re.error as exc:
            logger.warning("Regex error in PTM pattern: %s", exc)

    # Sort by position and remove duplicate positions for same PTM type
    results.sort(key=lambda r: (r["position"], r["type"]))

    seen: set[Tuple[int, str]] = set()
    deduped: List[Dict[str, Any]] = []
    for r in results:
        key = (r["position"], r["type"])
        if key not in seen:
            deduped.append(r)
            seen.add(key)

    logger.debug("Predicted %d PTM sites", len(deduped))
    return deduped
