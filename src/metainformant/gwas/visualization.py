"""Visualization for GWAS results (Manhattan plots, Q-Q plots, etc.)."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use("Agg", force=True)

from ..core.io import ensure_directory
from ..core.logging import get_logger

logger = get_logger(__name__)

# Try importing pandas for easier data handling
try:
    import pandas as pd

    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False
    pd = None


def manhattan_plot(
    results: list[dict[str, Any]] | Path,
    output_path: str | Path,
    *,
    significance_threshold: float | None = None,
    suggestiveness_threshold: float | None = None,
    chrom_colors: list[str] | None = None,
    title: str | None = None,
) -> dict[str, Any]:
    """Generate Manhattan plot for GWAS results.

    Args:
        results: List of result dictionaries or path to results file (TSV)
        output_path: Path to save plot
        significance_threshold: P-value threshold for genome-wide significance (default: 5e-8)
        suggestiveness_threshold: P-value threshold for suggestiveness (default: 1e-5)
        chrom_colors: List of colors to alternate for chromosomes
        title: Plot title

    Returns:
        Dictionary with plot metadata
    """
    logger.info(f"manhattan_plot: Generating Manhattan plot")

    if significance_threshold is None:
        significance_threshold = 5e-8
    if suggestiveness_threshold is None:
        suggestiveness_threshold = 1e-5

    # Load results if path provided
    if isinstance(results, Path) or (isinstance(results, str) and Path(results).exists()):
        if PANDAS_AVAILABLE:
            df = pd.read_csv(results, sep="\t")
            results_list = df.to_dict("records")
        else:
            # Simple TSV parsing using read_tsv for list of lists
            from ..core.io import read_tsv

            data_list = read_tsv(results)
            if len(data_list) < 2:
                return {"status": "failed", "error": "Invalid results file"}
            # First row is header
            header = data_list[0]
            results_list = []
            for row in data_list[1:]:
                if len(row) >= len(header):
                    results_list.append({header[i]: row[i] for i in range(len(header))})
    else:
        results_list = results

    if not results_list:
        return {"status": "failed", "error": "No results to plot"}

    # Extract data
    chroms: list[str] = []
    positions: list[int] = []
    pvalues: list[float] = []

    for result in results_list:
        chrom = str(result.get("CHROM", ""))
        pos = result.get("POS", 0)
        pval = result.get("p_value", 1.0)

        try:
            pos_int = int(pos)
            pval_float = float(pval)
            if pval_float > 0:
                chroms.append(chrom)
                positions.append(pos_int)
                pvalues.append(pval_float)
        except (ValueError, TypeError):
            continue

    if not pvalues:
        return {"status": "failed", "error": "No valid p-values found"}

    # Prepare chromosome positions for plotting
    unique_chroms = sorted(set(chroms), key=lambda x: (int(x[3:]) if x.startswith("chr") and x[3:].isdigit() else 999, x))
    chrom_positions: dict[str, list[int]] = {chrom: [] for chrom in unique_chroms}
    chrom_pvalues: dict[str, list[float]] = {chrom: [] for chrom in unique_chroms}

    for chrom, pos, pval in zip(chroms, positions, pvalues):
        if chrom in chrom_positions:
            chrom_positions[chrom].append(pos)
            chrom_pvalues[chrom].append(pval)

    # Calculate cumulative positions for x-axis
    cumulative_pos = 0
    chrom_cumulative: dict[str, int] = {}
    chrom_spans: dict[str, tuple[int, int]] = {}

    for chrom in unique_chroms:
        chrom_cumulative[chrom] = cumulative_pos
        if chrom_positions[chrom]:
            min_pos = min(chrom_positions[chrom])
            max_pos = max(chrom_positions[chrom])
            span = max_pos - min_pos
            chrom_spans[chrom] = (cumulative_pos, cumulative_pos + span)
            cumulative_pos += span + 10000000  # Gap between chromosomes

    # Create plot
    fig, ax = plt.subplots(figsize=(14, 6))

    # Colors for chromosomes
    if chrom_colors is None:
        chrom_colors = ["#2E86AB", "#A23B72", "#F18F01", "#C73E1D"]
    color_cycle = [chrom_colors[i % len(chrom_colors)] for i in range(len(unique_chroms))]

    # Plot points
    for chrom_idx, chrom in enumerate(unique_chroms):
        if chrom not in chrom_positions:
            continue

        positions_chrom = chrom_positions[chrom]
        pvalues_chrom = chrom_pvalues[chrom]
        x_pos = [chrom_cumulative[chrom] + pos - min(positions_chrom) for pos in positions_chrom]
        y_neg_log_p = [-math.log10(p) if p > 0 else 0 for p in pvalues_chrom]

        ax.scatter(
            x_pos,
            y_neg_log_p,
            c=color_cycle[chrom_idx],
            s=20,
            alpha=0.6,
            label=chrom if chrom_idx < 20 else None,  # Limit legend entries
        )

    # Add significance lines
    sig_line = -math.log10(significance_threshold)
    sugg_line = -math.log10(suggestiveness_threshold)

    ax.axhline(y=sig_line, color="r", linestyle="--", linewidth=1, label=f"Genome-wide significance ({significance_threshold:.0e})")
    ax.axhline(y=sugg_line, color="orange", linestyle="--", linewidth=1, label=f"Suggestive ({suggestiveness_threshold:.0e})")

    # Labels and title
    ax.set_xlabel("Chromosome", fontsize=12)
    ax.set_ylabel("-log10(p-value)", fontsize=12)
    if title:
        ax.set_title(title, fontsize=14, fontweight="bold")
    else:
        ax.set_title("Manhattan Plot", fontsize=14, fontweight="bold")

    # X-axis ticks at chromosome centers
    chrom_centers = [(chrom_spans[chrom][0] + chrom_spans[chrom][1]) / 2 for chrom in unique_chroms if chrom in chrom_spans]
    ax.set_xticks(chrom_centers)
    ax.set_xticklabels(unique_chroms, rotation=45, ha="right", fontsize=8)

    ax.grid(True, alpha=0.3, axis="y")
    ax.legend(loc="upper right", fontsize=8, ncol=2)

    plt.tight_layout()

    # Save plot
    output_path_obj = Path(output_path)
    ensure_directory(output_path_obj.parent)
    plt.savefig(output_path_obj, dpi=300, bbox_inches="tight")
    plt.close()

    logger.info(f"manhattan_plot: Saved Manhattan plot to {output_path_obj}")

    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "num_variants": len(pvalues),
        "num_significant": sum(1 for p in pvalues if p < significance_threshold),
    }


def qq_plot(
    pvalues: list[float] | Path,
    output_path: str | Path,
    *,
    title: str | None = None,
) -> dict[str, Any]:
    """Generate Q-Q plot for p-value distribution.

    Args:
        pvalues: List of p-values or path to results file with p_value column
        output_path: Path to save plot
        title: Plot title

    Returns:
        Dictionary with plot metadata
    """
    logger.info(f"qq_plot: Generating Q-Q plot")

    # Load p-values if path provided
    if isinstance(pvalues, Path) or (isinstance(pvalues, str) and Path(pvalues).exists()):
        if PANDAS_AVAILABLE:
            df = pd.read_csv(pvalues, sep="\t")
            pvalues_list = df["p_value"].tolist()
        else:
            from ..core.io import read_tsv

            data_list = read_tsv(pvalues)
            if len(data_list) < 2:
                return {"status": "failed", "error": "Invalid results file"}
            
            # First row is header
            header = data_list[0]
            pval_idx = header.index("p_value") if "p_value" in header else None
            if pval_idx is None:
                return {"status": "failed", "error": "p_value column not found"}
            
            pvalues_list = []
            for row in data_list[1:]:
                if len(row) > pval_idx:
                    try:
                        pval = float(row[pval_idx])
                        if 0 < pval <= 1:
                            pvalues_list.append(pval)
                    except (ValueError, TypeError):
                        continue
    else:
        pvalues_list = [float(p) for p in pvalues if isinstance(p, (int, float)) and 0 < p <= 1]

    if not pvalues_list:
        return {"status": "failed", "error": "No valid p-values found"}

    # Remove invalid p-values
    pvalues_arr = np.array(pvalues_list)
    pvalues_arr = pvalues_arr[(pvalues_arr > 0) & (pvalues_arr <= 1)]

    if len(pvalues_arr) == 0:
        return {"status": "failed", "error": "No valid p-values after filtering"}

    # Sort p-values
    observed_pvals = np.sort(pvalues_arr)

    # Expected p-values under null (uniform distribution)
    n = len(observed_pvals)
    expected_pvals = np.linspace(1.0 / n, 1.0, n)

    # Convert to -log10 scale
    observed_neg_log = -np.log10(observed_pvals)
    expected_neg_log = -np.log10(expected_pvals)

    # Calculate lambda_GC (genomic inflation factor)
    median_obs = np.median(observed_neg_log)
    median_exp = np.median(expected_neg_log)
    lambda_gc = median_obs / median_exp if median_exp > 0 else 1.0

    # Create plot
    fig, ax = plt.subplots(figsize=(8, 8))

    # Plot points
    ax.scatter(expected_neg_log, observed_neg_log, s=20, alpha=0.6, c="#2E86AB")

    # Plot diagonal line (expected under null)
    max_val = max(np.max(expected_neg_log), np.max(observed_neg_log))
    ax.plot([0, max_val], [0, max_val], "r--", linewidth=2, label="Expected (null)")

    # Labels and title
    ax.set_xlabel("Expected -log10(p-value)", fontsize=12)
    ax.set_ylabel("Observed -log10(p-value)", fontsize=12)
    if title:
        ax.set_title(f"{title}\nλ_GC = {lambda_gc:.4f}", fontsize=14, fontweight="bold")
    else:
        ax.set_title(f"Q-Q Plot\nλ_GC = {lambda_gc:.4f}", fontsize=14, fontweight="bold")

    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", fontsize=10)

    plt.tight_layout()

    # Save plot
    output_path_obj = Path(output_path)
    ensure_directory(output_path_obj.parent)
    plt.savefig(output_path_obj, dpi=300, bbox_inches="tight")
    plt.close()

    logger.info(f"qq_plot: Saved Q-Q plot to {output_path_obj}")

    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "lambda_gc": float(lambda_gc),
        "num_pvalues": len(pvalues_arr),
    }


def regional_plot(
    results: list[dict[str, Any]] | Path,
    region: str,
    output_path: str | Path,
    *,
    window: int = 500000,
    title: str | None = None,
) -> dict[str, Any]:
    """Generate regional association plot for a specific genomic region.

    Args:
        results: List of result dictionaries or path to results file
        region: Genomic region (format: "chr:start-end" or "chr:center")
        output_path: Path to save plot
        window: Window size around center if single position given (in bp)
        title: Plot title

    Returns:
        Dictionary with plot metadata
    """
    logger.info(f"regional_plot: Generating regional plot for {region}")

    # Parse region
    if ":" in region:
        parts = region.split(":")
        chrom = parts[0]
        if "-" in parts[1]:
            start_end = parts[1].split("-")
            start = int(start_end[0])
            end = int(start_end[1])
        else:
            center = int(parts[1])
            start = center - window
            end = center + window
    else:
        return {"status": "failed", "error": "Invalid region format"}

    # Load results if path provided
    if isinstance(results, Path) or (isinstance(results, str) and Path(results).exists()):
        if PANDAS_AVAILABLE:
            df = pd.read_csv(results, sep="\t")
            results_list = df.to_dict("records")
        else:
            from ..core.io import read_tsv

            data_list = read_tsv(results)
            if len(data_list) < 2:
                return {"status": "failed", "error": "Invalid results file"}
            # First row is header
            header = data_list[0]
            results_list = []
            for row in data_list[1:]:
                if len(row) >= len(header):
                    results_list.append({header[i]: row[i] for i in range(len(header))})
    else:
        results_list = results

    # Filter results by region
    regional_results = []
    for result in results_list:
        result_chrom = str(result.get("CHROM", ""))
        result_pos = result.get("POS", 0)
        try:
            result_pos_int = int(result_pos)
            if result_chrom == chrom and start <= result_pos_int <= end:
                regional_results.append(result)
        except (ValueError, TypeError):
            continue

    if not regional_results:
        return {"status": "failed", "error": f"No variants found in region {region}"}

    # Extract data
    positions = []
    pvalues = []
    for result in regional_results:
        try:
            pos = int(result.get("POS", 0))
            pval = float(result.get("p_value", 1.0))
            if pval > 0:
                positions.append(pos)
                pvalues.append(-math.log10(pval))
        except (ValueError, TypeError):
            continue

    if not positions:
        return {"status": "failed", "error": "No valid positions"}

    # Create plot
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.scatter(positions, pvalues, s=30, alpha=0.7, c="#2E86AB")

    # Significance line
    sig_line = -math.log10(5e-8)
    ax.axhline(y=sig_line, color="r", linestyle="--", linewidth=1, label="Genome-wide significance")

    # Labels and title
    ax.set_xlabel(f"Position on {chrom}", fontsize=12)
    ax.set_ylabel("-log10(p-value)", fontsize=12)
    if title:
        ax.set_title(title, fontsize=14, fontweight="bold")
    else:
        ax.set_title(f"Regional Association Plot: {region}", fontsize=14, fontweight="bold")

    ax.grid(True, alpha=0.3, axis="y")
    ax.legend(loc="upper right", fontsize=10)

    plt.tight_layout()

    # Save plot
    output_path_obj = Path(output_path)
    ensure_directory(output_path_obj.parent)
    plt.savefig(output_path_obj, dpi=300, bbox_inches="tight")
    plt.close()

    logger.info(f"regional_plot: Saved regional plot to {output_path_obj}")

    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "region": region,
        "num_variants": len(positions),
    }

