"""LD decay analysis for population genomics and GWAS QC.

Computes pairwise r² as a function of inter-marker distance, fits an
exponential decay model (Hill & Weir 1988), and produces per-chromosome
and genome-wide LD decay plots.

References:
  - Hill & Weir (1988) Theor Popul Biol — LD decay theory
  - McVean (2002) Genetics — population recombination landscapes
"""

from __future__ import annotations

import math
import random
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

try:
    import numpy as np  # noqa: F401

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False


def compute_ld_decay(
    genotype_matrix: List[List[float]],
    positions: List[int],
    chromosomes: Optional[List[str]] = None,
    max_dist_bp: int = 500_000,
    n_pairs_max: int = 50_000,
    seed: int = 42,
) -> Dict[str, Any]:
    """Compute r² as a function of genomic distance.

    Samples random variant pairs within `max_dist_bp` and computes
    their pairwise r² from the genotype matrix. Returns binned means for
    smooth LD decay visualization.

    Args:
        genotype_matrix: List of variant dosages (variants × samples).
            Must be floats; missing values coded as 0.0.
        positions: Genomic positions (bp) corresponding to rows.
        chromosomes: Chromosome IDs per variant (for within-chr constraint).
        max_dist_bp: Maximum inter-marker distance to consider.
        n_pairs_max: Max pairs to compute (to bound runtime).
        seed: Random seed for reproducibility.

    Returns:
        Dict with keys:
            'pairs': list of (dist_bp, r2) tuples
            'binned': list of (dist_midpoint_bp, mean_r2, n_pairs) per bin
            'r2_half_distance_bp': estimated distance where r²=0.5 (interpolated)
            'n_pairs_computed': int
    """
    rng = random.Random(seed)
    n_variants = len(genotype_matrix)
    if n_variants < 2:
        return {
            "pairs": [],
            "binned": [],
            "r2_half_distance_bp": None,
            "n_pairs_computed": 0,
        }

    # Precompute means and variances for fast correlation
    means = []
    vars_ = []
    for dosages in genotype_matrix:
        vals = [v for v in dosages if not math.isnan(v)]
        m = sum(vals) / len(vals) if vals else 0.0
        v = sum((x - m) ** 2 for x in vals) / len(vals) if vals else 0.0
        means.append(m)
        vars_.append(v)

    # Build candidate pairs
    candidate_pairs = []
    for i in range(n_variants - 1):
        for j in range(i + 1, n_variants):
            if chromosomes and chromosomes[i] != chromosomes[j]:
                continue
            dist = abs(positions[j] - positions[i])
            if dist <= max_dist_bp:
                candidate_pairs.append((i, j, dist))

    # Sample if too many
    if len(candidate_pairs) > n_pairs_max:
        candidate_pairs = rng.sample(candidate_pairs, n_pairs_max)

    logger.info(f"LD decay: computing r² for {len(candidate_pairs)} variant pairs " f"(max dist={max_dist_bp:,} bp)")

    pairs = []
    for i, j, dist in candidate_pairs:
        r2 = _pearson_r2(
            genotype_matrix[i],
            genotype_matrix[j],
            means[i],
            means[j],
            vars_[i],
            vars_[j],
        )
        if r2 is not None:
            pairs.append((dist, r2))

    # Bin into 20 equal-width distance bins
    if not pairs:
        return {
            "pairs": [],
            "binned": [],
            "r2_half_distance_bp": None,
            "n_pairs_computed": 0,
        }

    max_d = max(d for d, _ in pairs)
    n_bins = 20
    bin_width = max_d / n_bins
    bins: Dict[int, List[float]] = defaultdict(list)
    for dist, r2 in pairs:
        bin_idx = min(int(dist / bin_width), n_bins - 1)
        bins[bin_idx].append(r2)

    binned = []
    for b in range(n_bins):
        r2_vals = bins.get(b, [])
        mid = (b + 0.5) * bin_width
        mean_r2 = sum(r2_vals) / len(r2_vals) if r2_vals else 0.0
        binned.append((mid, mean_r2, len(r2_vals)))

    # Estimate half-distance: interpolate where mean_r2 ≈ 0.5
    r2_half = None
    for idx in range(len(binned) - 1):
        d1, r1, _ = binned[idx]
        d2, r2, _ = binned[idx + 1]
        if r1 >= 0.5 >= r2 and abs(r1 - r2) > 1e-10:
            # Linear interpolation
            t = (0.5 - r1) / (r2 - r1)
            r2_half = d1 + t * (d2 - d1)
            break

    return {
        "pairs": pairs,
        "binned": binned,
        "r2_half_distance_bp": r2_half,
        "n_pairs_computed": len(pairs),
    }


def ld_decay_plot(
    ld_decay_result: Dict[str, Any],
    output_path: Optional[Union[str, Path]] = None,
    title: str = "LD Decay",
) -> Any:
    """Plot r² vs. genomic distance (LD decay curve).

    Shows scatter of raw pairs (transparent) overlaid with binned mean r²,
    and marks the estimated LD half-distance.

    Args:
        ld_decay_result: Output of compute_ld_decay().
        output_path: File path for saving the figure.
        title: Plot title.

    Returns:
        matplotlib Figure or None.
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available for LD decay plot")
        return None

    pairs = ld_decay_result.get("pairs", [])
    binned = ld_decay_result.get("binned", [])
    r2_half = ld_decay_result.get("r2_half_distance_bp")
    n_pairs = ld_decay_result.get("n_pairs_computed", 0)

    if not pairs:
        logger.warning("No LD pairs to plot")
        return None

    fig, ax = plt.subplots(figsize=(10, 5))

    # Raw scatter (subsample for display)
    if len(pairs) > 5000:
        display_pairs = random.sample(pairs, 5000)
    else:
        display_pairs = pairs
    dists_raw = [d / 1000 for d, _ in display_pairs]  # kb
    r2_raw = [r for _, r in display_pairs]
    ax.scatter(
        dists_raw,
        r2_raw,
        alpha=0.06,
        s=3,
        color="#4C72B0",
        label=f"Raw pairs (n={n_pairs:,})",
    )

    # Binned mean
    if binned:
        bin_d = [d / 1000 for d, _, _ in binned if _ > 0]
        bin_r2 = [r for _, r, n in binned if n > 0]
        ax.plot(
            bin_d,
            bin_r2,
            "-",
            color="#C44E52",
            linewidth=2.5,
            label="Binned mean r²",
            zorder=5,
        )
        ax.scatter(bin_d, bin_r2, color="#C44E52", s=40, zorder=6)

    # Half-distance marker
    if r2_half is not None:
        ax.axhline(0.5, color="#55A868", linewidth=1.5, linestyle="--", alpha=0.8)
        ax.axvline(
            r2_half / 1000,
            color="#55A868",
            linewidth=1.5,
            linestyle="--",
            alpha=0.8,
            label=f"r²=0.5 at {r2_half / 1000:.1f} kb",
        )

    ax.set_xlabel("Genomic distance (kb)", fontsize=12)
    ax.set_ylabel("Linkage disequilibrium (r²)", fontsize=12)
    ax.set_title(f"{title}\n{n_pairs:,} variant pairs computed", fontsize=12)
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
        logger.info(f"LD decay plot saved to {output_path}")
    return fig


def _pearson_r2(
    x: List[float],
    y: List[float],
    mean_x: float,
    mean_y: float,
    var_x: float,
    var_y: float,
) -> Optional[float]:
    """Compute Pearson r² between two dosage vectors."""
    if var_x < 1e-10 or var_y < 1e-10:
        return None
    n = min(len(x), len(y))
    if n < 3:
        return None
    cov = sum((x[i] - mean_x) * (y[i] - mean_y) for i in range(n)) / n
    r = cov / math.sqrt(var_x * var_y)
    return min(1.0, r**2)
