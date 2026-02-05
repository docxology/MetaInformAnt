"""Fine-mapping visualization functions for GWAS.

This module provides plots for fine-mapping analysis including credible sets,
posterior inclusion probabilities (PIPs), and conditional analysis.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core import logging

logger = logging.get_logger(__name__)

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    plt = None

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None


def compute_credible_set(
    assoc_results: List[Dict[str, Any]],
    credible_level: float = 0.95,
) -> Dict[str, Any]:
    """Compute approximate Bayes factor credible set from association results.

    Uses approximate Bayes factors derived from p-values to compute posterior
    inclusion probabilities (PIPs) and identify the credible set of variants.

    Args:
        assoc_results: List of dicts, each with at least a "p_value" key.
        credible_level: Cumulative PIP threshold for the credible set (default 0.95).

    Returns:
        Dict with status, PIPs, credible set indices, size, and cumulative probability.
    """
    if not HAS_NUMPY:
        logger.warning("numpy not available for credible set computation")
        return {"status": "skipped", "reason": "numpy not available"}

    if not assoc_results:
        logger.error("No association results provided")
        return {"status": "failed", "reason": "No association results provided"}

    # Extract p-values
    p_values = []
    for result in assoc_results:
        p = result.get("p_value")
        if p is None:
            p = result.get("P", 1.0)
        p_values.append(float(p))

    p_values_arr = np.array(p_values, dtype=np.float64)

    # Clamp p-values to avoid log(0)
    p_values_arr = np.clip(p_values_arr, 1e-300, 1.0)

    # Compute approximate chi-squared statistic: chi2 ~ -2 * log(p)
    chi2_stats = -2.0 * np.log(p_values_arr)

    # Compute approximate Bayes factors: ABF_i = exp(0.5 * chi2_i)
    # Larger chi2 (smaller p-value) => larger ABF => more evidence for association
    # Using log-space for numerical stability
    log_abf = 0.5 * chi2_stats

    # Normalize to avoid overflow: subtract max before exponentiation
    log_abf_shifted = log_abf - np.max(log_abf)
    abf = np.exp(log_abf_shifted)

    # Compute posterior inclusion probabilities (PIPs)
    total_abf = np.sum(abf)
    if total_abf == 0:
        logger.error("Total ABF is zero; cannot compute PIPs")
        return {"status": "failed", "reason": "Total ABF is zero"}

    pips = abf / total_abf

    # Sort by PIP descending and accumulate until credible_level
    sorted_indices = np.argsort(-pips)
    cumulative_pip = 0.0
    credible_set_indices: List[int] = []

    for idx in sorted_indices:
        credible_set_indices.append(int(idx))
        cumulative_pip += float(pips[idx])
        if cumulative_pip >= credible_level:
            break

    return {
        "status": "success",
        "pips": pips.tolist(),
        "credible_set_indices": credible_set_indices,
        "credible_set_size": len(credible_set_indices),
        "cumulative_probability": float(cumulative_pip),
    }


def credible_set_plot(
    assoc_results: List[Dict[str, Any]],
    output_file: Optional[Union[str, Path]] = None,
    credible_level: float = 0.95,
    ld_matrix: Optional[Any] = None,
    title: str = "Fine-Mapping Credible Set",
) -> Dict[str, Any]:
    """Create a credible set plot showing posterior inclusion probabilities.

    Scatter plot with variant position on x-axis and PIP on y-axis. Variants
    in the credible set are highlighted. Optionally colored by LD r-squared
    with the lead variant.

    Args:
        assoc_results: List of dicts with "p_value" and optionally "position" keys.
        output_file: Optional output file path.
        credible_level: Cumulative PIP threshold for the credible set.
        ld_matrix: Optional LD matrix (numpy array, n_variants x n_variants).
            Used to color variants by LD r-squared with the lead variant.
        title: Plot title.

    Returns:
        Dict with status, output_path, and credible_set_size.
    """
    if not HAS_MATPLOTLIB or not HAS_NUMPY:
        logger.warning("matplotlib/numpy not available for credible set plot")
        return {"status": "skipped", "reason": "matplotlib or numpy not available", "output_path": None}

    if not assoc_results:
        logger.error("No association results provided")
        return {"status": "failed", "reason": "No association results", "output_path": None}

    # Compute credible set
    cs_result = compute_credible_set(assoc_results, credible_level=credible_level)
    if cs_result["status"] != "success":
        return {
            "status": cs_result["status"],
            "reason": cs_result.get("reason", "Credible set computation failed"),
            "output_path": None,
        }

    pips = np.array(cs_result["pips"])
    credible_indices = set(cs_result["credible_set_indices"])

    # Extract positions (use variant index if no position provided)
    positions = []
    for i, result in enumerate(assoc_results):
        pos = result.get("position", result.get("POS", result.get("BP", i)))
        positions.append(float(pos))
    positions_arr = np.array(positions)

    # Determine the lead variant (highest PIP)
    lead_idx = int(np.argmax(pips))

    # Determine PIP threshold at the boundary of the credible set
    if credible_indices:
        pip_threshold = min(pips[idx] for idx in credible_indices)
    else:
        pip_threshold = 0.0

    fig, ax = plt.subplots(figsize=(10, 6))

    if ld_matrix is not None:
        # Color by LD r-squared with lead variant
        ld_arr = np.array(ld_matrix)
        ld_with_lead = ld_arr[lead_idx, :]

        # Scatter all points colored by LD
        scatter = ax.scatter(
            positions_arr,
            pips,
            c=ld_with_lead,
            cmap="coolwarm",
            vmin=0,
            vmax=1,
            s=50,
            alpha=0.8,
            edgecolors="black",
            linewidths=0.5,
            zorder=2,
        )
        cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
        cbar.set_label("LD r$^2$ with lead variant")

        # Highlight credible set with larger markers
        cs_mask = np.array([i in credible_indices for i in range(len(pips))])
        ax.scatter(
            positions_arr[cs_mask],
            pips[cs_mask],
            facecolors="none",
            edgecolors="red",
            s=120,
            linewidths=2,
            zorder=3,
            label=f"Credible set (n={cs_result['credible_set_size']})",
        )
    else:
        # Color by credible set membership
        colors = ["red" if i in credible_indices else "gray" for i in range(len(pips))]
        ax.scatter(
            positions_arr,
            pips,
            c=colors,
            s=50,
            alpha=0.8,
            edgecolors="black",
            linewidths=0.5,
            zorder=2,
        )
        # Legend entries
        ax.scatter([], [], c="red", s=50, label=f"In credible set (n={cs_result['credible_set_size']})")
        ax.scatter([], [], c="gray", s=50, label="Outside credible set")

    # Draw horizontal dashed line at PIP threshold
    ax.axhline(
        y=pip_threshold,
        color="blue",
        linestyle="--",
        alpha=0.6,
        linewidth=1,
        label=f"PIP threshold ({pip_threshold:.3f})",
    )

    ax.set_xlabel("Variant Position", fontsize=12)
    ax.set_ylabel("Posterior Inclusion Probability (PIP)", fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.set_ylim(-0.02, max(pips) * 1.1 + 0.02)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    output_path_str: Optional[str] = None
    if output_file:
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        output_path_str = str(output_file)
        logger.info(f"Saved credible set plot to {output_file}")

    plt.close(fig)

    return {
        "status": "success",
        "output_path": output_path_str,
        "credible_set_size": cs_result["credible_set_size"],
    }


def conditional_analysis_plot(
    assoc_results_list: List[List[Dict[str, Any]]],
    output_file: Optional[Union[str, Path]] = None,
    labels: Optional[List[str]] = None,
    title: str = "Conditional Analysis",
) -> Dict[str, Any]:
    """Create an overlay plot of multiple rounds of conditional association analysis.

    Each round is plotted as a separate trace with its own color/marker.

    Args:
        assoc_results_list: List of association result lists (one per round).
        output_file: Optional output file path.
        labels: Optional list of labels for each round
            (default: "Round 1", "Round 2", etc.).
        title: Plot title.

    Returns:
        Dict with status, output_path, and n_rounds.
    """
    if not HAS_MATPLOTLIB or not HAS_NUMPY:
        logger.warning("matplotlib/numpy not available for conditional analysis plot")
        return {"status": "skipped", "reason": "matplotlib or numpy not available", "output_path": None}

    if not assoc_results_list:
        logger.error("No association result rounds provided")
        return {"status": "failed", "reason": "No association results", "output_path": None}

    n_rounds = len(assoc_results_list)

    if labels is None:
        labels = [f"Round {i + 1}" for i in range(n_rounds)]

    # Color and marker cycles
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f"]
    markers = ["o", "s", "^", "D", "v", "<", ">", "p"]

    fig, ax = plt.subplots(figsize=(12, 6))

    for round_idx, (results, label) in enumerate(zip(assoc_results_list, labels)):
        if not results:
            continue

        # Extract positions and p-values
        positions = []
        neg_log_p = []
        for i, result in enumerate(results):
            pos = result.get("position", result.get("POS", result.get("BP", i)))
            p = result.get("p_value", result.get("P", 1.0))
            p = max(float(p), 1e-300)
            positions.append(float(pos))
            neg_log_p.append(-math.log10(p))

        color = colors[round_idx % len(colors)]
        marker = markers[round_idx % len(markers)]

        ax.scatter(
            positions,
            neg_log_p,
            c=color,
            marker=marker,
            s=30,
            alpha=0.7,
            label=label,
            zorder=2 + round_idx,
        )

    # Genome-wide significance threshold
    ax.axhline(
        y=-math.log10(5e-8),
        color="red",
        linestyle="--",
        alpha=0.7,
        linewidth=1,
        label="Genome-wide significance",
    )

    ax.set_xlabel("Variant Position", fontsize=12)
    ax.set_ylabel("-log$_{10}$(p-value)", fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    output_path_str: Optional[str] = None
    if output_file:
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        output_path_str = str(output_file)
        logger.info(f"Saved conditional analysis plot to {output_file}")

    plt.close(fig)

    return {
        "status": "success",
        "output_path": output_path_str,
        "n_rounds": n_rounds,
    }


def pip_vs_ld_plot(
    pips: List[float],
    ld_with_lead: List[float],
    output_file: Optional[Union[str, Path]] = None,
    variant_ids: Optional[List[str]] = None,
    title: str = "PIP vs LD with Lead Variant",
) -> Dict[str, Any]:
    """Create a scatter plot of PIP vs LD r-squared with the lead variant.

    Highlights variants in the credible set (PIP > 0.5 by default) and draws
    quadrant lines at PIP=0.5 and LD=0.5.

    Args:
        pips: List of posterior inclusion probabilities.
        ld_with_lead: List of LD r-squared values with the lead variant.
        output_file: Optional output file path.
        variant_ids: Optional list of variant IDs for annotation.
        title: Plot title.

    Returns:
        Dict with status, output_path, and n_variants.
    """
    if not HAS_MATPLOTLIB or not HAS_NUMPY:
        logger.warning("matplotlib/numpy not available for PIP vs LD plot")
        return {"status": "skipped", "reason": "matplotlib or numpy not available", "output_path": None}

    if not pips or not ld_with_lead:
        logger.error("PIP or LD data not provided")
        return {"status": "failed", "reason": "Missing PIP or LD data", "output_path": None}

    if len(pips) != len(ld_with_lead):
        logger.error("PIP and LD arrays must have the same length")
        return {"status": "failed", "reason": "Length mismatch between PIP and LD arrays", "output_path": None}

    pips_arr = np.array(pips, dtype=np.float64)
    ld_arr = np.array(ld_with_lead, dtype=np.float64)
    n_variants = len(pips_arr)

    fig, ax = plt.subplots(figsize=(8, 8))

    # Color by PIP threshold (high PIP = in credible set)
    high_pip_mask = pips_arr >= 0.5
    low_pip_mask = ~high_pip_mask

    ax.scatter(
        ld_arr[low_pip_mask],
        pips_arr[low_pip_mask],
        c="gray",
        s=40,
        alpha=0.6,
        label="PIP < 0.5",
        zorder=2,
    )
    ax.scatter(
        ld_arr[high_pip_mask],
        pips_arr[high_pip_mask],
        c="red",
        s=60,
        alpha=0.8,
        edgecolors="black",
        linewidths=0.5,
        label="PIP >= 0.5",
        zorder=3,
    )

    # Quadrant lines
    ax.axhline(y=0.5, color="blue", linestyle="--", alpha=0.5, linewidth=1)
    ax.axvline(x=0.5, color="blue", linestyle="--", alpha=0.5, linewidth=1)

    # Annotate top-PIP variants
    if variant_ids is not None and len(variant_ids) == n_variants:
        # Annotate the top 5 PIP variants
        top_indices = np.argsort(-pips_arr)[:5]
        for idx in top_indices:
            if pips_arr[idx] > 0.1:  # Only annotate meaningful PIPs
                ax.annotate(
                    variant_ids[idx],
                    (ld_arr[idx], pips_arr[idx]),
                    textcoords="offset points",
                    xytext=(5, 5),
                    fontsize=8,
                    alpha=0.9,
                )

    ax.set_xlabel("LD r$^2$ with Lead Variant", fontsize=12)
    ax.set_ylabel("Posterior Inclusion Probability (PIP)", fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.legend(fontsize=10, loc="upper left")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    output_path_str: Optional[str] = None
    if output_file:
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        output_path_str = str(output_file)
        logger.info(f"Saved PIP vs LD plot to {output_file}")

    plt.close(fig)

    return {
        "status": "success",
        "output_path": output_path_str,
        "n_variants": n_variants,
    }
