"""Statistical visualization functions for GWAS.

This module provides plots for statistical analysis and quality control.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

from metainformant.core import logging

logger = logging.get_logger(__name__)


def lambda_gc_plot(
    results: List[Any],
    output_file: Optional[str | Path] = None,
    title: str = "Genomic Control Lambda Distribution",
) -> dict[str, Any]:
    """Create a plot showing the distribution of genomic control lambda values.

    Args:
        results: List of GWAS result dicts (with 'p_value' and 'CHROM' keys) or list of float lambda values
        output_file: Optional output file path
        title: Plot title

    Returns:
        Dictionary with plot status and metadata

    Example:
        >>> lambdas = [1.05, 1.12, 1.08, 1.15]
        >>> plot = lambda_gc_plot(lambdas)
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for lambda GC plot")
        return {"status": "failed", "error": "matplotlib not available"}

    if not results:
        logger.error("No data provided")
        return {"status": "failed", "error": "No data provided"}

    try:
        # Accept list-of-dicts: compute per-chromosome lambda GC from p-values
        if isinstance(results[0], dict):
            from collections import defaultdict

            chrom_pvals: dict[str, list[float]] = defaultdict(list)
            for r in results:
                p = r.get("p_value")
                chrom = r.get("CHROM", "unknown")
                if p is not None:
                    chrom_pvals[chrom].append(p)

            lambda_by_chrom: dict[str, float] = {}
            lambda_values: list[float] = []
            for chrom, pvals in chrom_pvals.items():
                if len(pvals) >= 2:
                    obs = sorted([-math.log10(max(p, 1e-300)) for p in pvals])
                    n = len(obs)
                    exp = sorted([-math.log10((i + 1) / (n + 1)) for i in range(n)])
                    med_obs = float(np.median(obs))
                    med_exp = float(np.median(exp))
                    lam = med_obs / med_exp if med_exp != 0 else 1.0
                    lambda_by_chrom[chrom] = lam
                    lambda_values.append(lam)

            if not lambda_values:
                return {"status": "failed", "error": "Not enough data to compute lambda GC"}

            # Also compute overall lambda GC
            all_pvals = []
            for pvals in chrom_pvals.values():
                all_pvals.extend(pvals)
            obs_all = sorted([-math.log10(max(p, 1e-300)) for p in all_pvals])
            n_all = len(obs_all)
            exp_all = sorted([-math.log10((i + 1) / (n_all + 1)) for i in range(n_all)])
            overall_lambda = (
                float(np.median(obs_all)) / float(np.median(exp_all)) if float(np.median(exp_all)) != 0 else 1.0
            )
        else:
            lambda_values = list(results)
            lambda_by_chrom = {}
            overall_lambda = float(np.mean(lambda_values))

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        # Plot 1: Histogram of lambda values
        ax1.hist(lambda_values, bins=max(1, min(20, len(lambda_values))), alpha=0.7, color="skyblue", edgecolor="black")
        ax1.axvline(
            x=np.mean(lambda_values),
            color="red",
            linestyle="--",
            linewidth=2,
            label=f"Mean: {np.mean(lambda_values):.3f}",
        )
        ax1.axvline(x=1.0, color="green", linestyle="-", alpha=0.7, label="Expected (1.0)")
        ax1.set_xlabel("Genomic Control Lambda", fontsize=12)
        ax1.set_ylabel("Frequency", fontsize=12)
        ax1.set_title("Lambda GC Distribution", fontsize=14)
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Plot 2: Q-Q plot against expected distribution
        expected_lambdas = np.random.normal(1.0, 0.1, len(lambda_values))
        expected_lambdas = np.sort(expected_lambdas)
        observed_lambdas = np.sort(lambda_values)

        ax2.scatter(expected_lambdas, observed_lambdas, alpha=0.7, color="orange")
        ax2.plot(
            [min(expected_lambdas), max(expected_lambdas)],
            [min(expected_lambdas), max(expected_lambdas)],
            "k--",
            alpha=0.7,
            label="Expected",
        )

        ax2.set_xlabel("Expected Lambda GC", fontsize=12)
        ax2.set_ylabel("Observed Lambda GC", fontsize=12)
        ax2.set_title("Q-Q Plot", fontsize=14)
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        # Overall title
        fig.suptitle(title, fontsize=16, y=0.98)

        # Add summary statistics as text
        stats_text = f"""Summary Statistics:
Mean: {np.mean(lambda_values):.3f}
Median: {np.median(lambda_values):.3f}
SD: {np.std(lambda_values):.3f}
Range: {np.min(lambda_values):.3f} - {np.max(lambda_values):.3f}
Inflated: {np.mean(lambda_values) > 1.1}"""

        fig.text(
            0.02,
            0.02,
            stats_text,
            fontsize=10,
            verticalalignment="bottom",
            fontfamily="monospace",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
        )

        plt.tight_layout()

        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches="tight")
            logger.info(f"Saved lambda GC plot to {output_file}")

        result: dict[str, Any] = {
            "status": "success",
            "overall_lambda_gc": overall_lambda,
            "output_path": str(output_file) if output_file else None,
        }
        if lambda_by_chrom:
            result["lambda_gc_by_chrom"] = lambda_by_chrom
        return result

    except Exception as e:
        return {"status": "failed", "error": str(e)}


def power_plot(
    sample_sizes: List[int],
    effect_sizes: List[float],
    alpha: float = 0.05,
    output_file: Optional[str | Path] = None,
    output_path: Optional[str | Path] = None,
    title: str = "Statistical Power Analysis",
) -> dict[str, Any]:
    """Create a power analysis plot showing statistical power vs sample size and effect size.

    Args:
        sample_sizes: List of sample sizes to evaluate
        effect_sizes: List of effect sizes to evaluate
        alpha: Significance level
        output_file: Optional output file path
        output_path: Alias for output_file
        title: Plot title

    Returns:
        Dictionary with plot status and metadata
    """
    # Allow output_path as alias for output_file
    if output_path is not None and output_file is None:
        output_file = output_path

    try:
        import matplotlib.pyplot as plt
        import numpy as np
        from scipy import stats
    except ImportError:
        logger.warning("matplotlib or scipy not available for power plot")
        return {"status": "failed", "error": "matplotlib or scipy not available"}

    if not sample_sizes or not effect_sizes:
        logger.error("No sample sizes or effect sizes provided")
        return {"status": "failed", "error": "No sample sizes or effect sizes provided"}

    # Create meshgrid for power calculation
    n_grid, effect_grid = np.meshgrid(sample_sizes, effect_sizes)
    power_grid = np.zeros_like(n_grid, dtype=float)

    # Calculate power for each combination
    for i in range(len(effect_sizes)):
        for j in range(len(sample_sizes)):
            n = sample_sizes[j]
            effect = effect_sizes[i]

            # For simplicity, assume two-sample t-test power
            # Power = 1 - β, where β is probability of Type II error
            try:
                # Calculate non-centrality parameter
                ncp = effect * np.sqrt(n / 2)  # Approximation for equal sample sizes
                # Critical value for alpha
                t_crit = stats.t.ppf(1 - alpha / 2, df=2 * n - 2)
                # Power calculation
                power = 1 - stats.nct.cdf(t_crit, df=2 * n - 2, nc=ncp)
                power_grid[i, j] = power
            except (ValueError, ZeroDivisionError):
                power_grid[i, j] = np.nan

    # Create plot
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot power surface
    cs = ax.contourf(n_grid, effect_grid, power_grid, levels=np.linspace(0, 1, 11), cmap="RdYlBu_r", alpha=0.8)
    plt.colorbar(cs, ax=ax, label="Statistical Power")

    # Add contour lines
    contours = ax.contour(n_grid, effect_grid, power_grid, levels=[0.8, 0.9], colors="black", linewidths=2)
    ax.clabel(contours, inline=True, fontsize=10)

    ax.set_xlabel("Sample Size", fontsize=12)
    ax.set_ylabel("Effect Size", fontsize=12)
    ax.set_title(title, fontsize=14, pad=20)
    ax.grid(True, alpha=0.3)

    # Add reference lines
    ax.axhline(y=0.2, color="red", linestyle="--", alpha=0.7, label="Small effect (d=0.2)")
    ax.axhline(y=0.5, color="orange", linestyle="--", alpha=0.7, label="Medium effect (d=0.5)")
    ax.axhline(y=0.8, color="green", linestyle="--", alpha=0.7, label="Large effect (d=0.8)")

    ax.legend()

    # Add text annotations for power thresholds
    ax.text(
        0.02,
        0.98,
        "Power Contours:\n0.8 (solid)\n0.9 (dashed)",
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved power analysis plot to {output_file}")

    return {
        "status": "success",
        "n_sample_sizes": len(sample_sizes),
        "n_effect_sizes": len(effect_sizes),
        "alpha": alpha,
        "output_path": str(output_file) if output_file else None,
    }


def qq_plot(
    results: list[Any], output_path: Optional[str | Path] = None, figsize: tuple[int, int] = (8, 8)
) -> dict[str, Any]:
    """Create a Q-Q plot for p-values.

    Args:
        results: List of p-values (floats) or list of result dicts with 'p_value' key
        output_path: Path to save the plot
        figsize: Figure size

    Returns:
        Dictionary with plot status and metadata
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np

        HAS_MATPLOTLIB = True
    except ImportError:
        HAS_MATPLOTLIB = False

    if not HAS_MATPLOTLIB:
        return {"status": "failed", "error": "matplotlib not available"}

    try:
        # Accept list-of-dicts: extract p_value from each dict
        if results and isinstance(results[0], dict):
            p_values: list[float] = [r["p_value"] for r in results if r.get("p_value") is not None]
        else:
            p_values = list(results)

        # Remove NA values
        p_values = [p for p in p_values if p is not None and not math.isnan(p)]

        if not p_values:
            return {"status": "failed", "error": "No valid p-values provided"}

        # Sort p-values
        p_values = sorted(p_values)
        n = len(p_values)

        # Expected p-values under null hypothesis
        expected = [-math.log10((i + 1) / (n + 1)) for i in range(n)]

        # Observed p-values
        observed = [-math.log10(p) for p in p_values]

        # Calculate lambda GC (genomic control inflation factor)
        median_observed = np.median(observed)
        median_expected = np.median(expected)
        lambda_gc = median_observed / median_expected if median_expected != 0 else 1.0

        # Create plot
        fig, ax = plt.subplots(figsize=figsize)

        ax.scatter(expected, observed, alpha=0.6, s=20, color="blue")

        # Add diagonal line
        max_val = max(max(expected), max(observed))
        ax.plot([0, max_val], [0, max_val], "r--", alpha=0.7, label="Expected")

        # Add confidence interval (simplified)
        ci_lower = [e - 0.5 for e in expected]
        ci_upper = [e + 0.5 for e in expected]
        ax.fill_between(expected, ci_lower, ci_upper, alpha=0.1, color="gray", label="95% CI")

        ax.set_xlabel("Expected -log₁₀(p)")
        ax.set_ylabel("Observed -log₁₀(p)")
        ax.set_title(f"Q-Q Plot\nλ = {lambda_gc:.3f}")
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches="tight")

        return {
            "status": "success",
            "n_p_values": len(p_values),
            "lambda_gc": lambda_gc,
            "median_expected": median_expected,
            "median_observed": median_observed,
            "output_path": str(output_path) if output_path else None,
        }

    except Exception as e:
        return {"status": "failed", "error": str(e)}


def qq_plot_stratified(
    results: list[Any],
    output_path: Optional[str | Path] = None,
    strata: Optional[list[str]] = None,
    maf_bins: Optional[list[tuple[float, float]]] = None,
    figsize: tuple[int, int] = (15, 10),
) -> dict[str, Any]:
    """Create stratified Q-Q plots for p-values by different groups.

    Args:
        results: List of result dicts with 'p_value' (and optionally 'MAF') keys, or list of p-values
        output_path: Path to save the plot
        strata: List of stratum labels (same length as results, used when results are floats)
        maf_bins: Optional list of (min, max) MAF bin tuples for stratification
        figsize: Figure size

    Returns:
        Dictionary with plot status and metadata
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np

        HAS_MATPLOTLIB = True
    except ImportError:
        HAS_MATPLOTLIB = False

    if not HAS_MATPLOTLIB:
        return {"status": "failed", "error": "matplotlib not available"}

    try:
        # Accept list-of-dicts with MAF bins
        if results and isinstance(results[0], dict) and maf_bins is not None:
            stratum_data: dict[str, list[float]] = {}
            for r in results:
                p = r.get("p_value")
                maf = r.get("MAF")
                if p is None or maf is None:
                    continue
                for lo, hi in maf_bins:
                    label = f"MAF [{lo:.2f}, {hi:.2f})"
                    if lo <= maf < hi:
                        if label not in stratum_data:
                            stratum_data[label] = []
                        stratum_data[label].append(p)
                        break
            if not stratum_data:
                return {"status": "failed", "error": "No valid data for MAF bins"}
        elif results and isinstance(results[0], dict):
            # List-of-dicts without MAF bins: single stratum
            p_values_list = [r["p_value"] for r in results if r.get("p_value") is not None]
            stratum_data = {"all": p_values_list}
        elif strata is not None:
            p_values_list = list(results)
            if len(p_values_list) != len(strata):
                return {"status": "failed", "error": "p_values and strata must have same length"}
            stratum_data = {}
            for p_val, stratum in zip(p_values_list, strata):
                if stratum not in stratum_data:
                    stratum_data[stratum] = []
                stratum_data[stratum].append(p_val)
        else:
            return {"status": "failed", "error": "Must provide strata or maf_bins for stratification"}

        if not stratum_data:
            return {"status": "failed", "error": "No valid data"}

        n_strata = len(stratum_data)
        n_cols = min(3, n_strata)
        n_rows = (n_strata + n_cols - 1) // n_cols

        fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, squeeze=False)
        axes = axes.flatten()

        strata_results: dict[str, Any] = {}

        for i, (stratum, stratum_p_values) in enumerate(stratum_data.items()):
            if i >= len(axes):
                break

            ax = axes[i]

            # Remove NA values
            stratum_p_values = [p for p in stratum_p_values if p is not None and not math.isnan(p)]

            if not stratum_p_values:
                ax.text(0.5, 0.5, f"No valid p-values\nfor {stratum}", ha="center", va="center", transform=ax.transAxes)
                continue

            # Sort p-values
            stratum_p_values = sorted(stratum_p_values)
            n = len(stratum_p_values)

            # Expected p-values under null hypothesis
            expected = [-math.log10((j + 1) / (n + 1)) for j in range(n)]

            # Observed p-values
            observed = [-math.log10(p) for p in stratum_p_values]

            # Calculate lambda GC
            median_observed = np.median(observed)
            median_expected = np.median(expected)
            lambda_gc = median_observed / median_expected if median_expected != 0 else 1.0

            # Plot
            ax.scatter(expected, observed, alpha=0.6, s=20, color=f"C{i}")
            ax.plot([0, max(expected + observed)], [0, max(expected + observed)], "k--", alpha=0.7)

            ax.set_xlabel("Expected -log₁₀(p)")
            ax.set_ylabel("Observed -log₁₀(p)")
            ax.set_title(f"{stratum}\nλ = {lambda_gc:.3f}")
            ax.grid(True, alpha=0.3)

            strata_results[stratum] = {
                "n_p_values": n,
                "lambda_gc": lambda_gc,
                "median_expected": median_expected,
                "median_observed": median_observed,
            }

        # Hide unused subplots
        for i in range(len(stratum_data), len(axes)):
            axes[i].set_visible(False)

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches="tight")

        return {
            "status": "success",
            "n_strata": len(stratum_data),
            "strata_results": strata_results,
            "output_path": str(output_path) if output_path else None,
        }

    except Exception as e:
        return {"status": "failed", "error": str(e)}


def volcano_plot(
    results: list[Any],
    output_path: Optional[str | Path] = None,
    significance_threshold: float = 5e-8,
    effect_size_threshold: float = 0.5,
    figsize: tuple[int, int] = (10, 8),
) -> dict[str, Any]:
    """Create a volcano plot for GWAS results.

    Args:
        results: List of result dicts with 'p_value' and 'beta' keys, or list of effect sizes
            (if list of floats, a second positional arg for p_values is expected via output_path)
        output_path: Path to save the plot
        significance_threshold: P-value threshold for significance
        effect_size_threshold: Effect size threshold for highlighting
        figsize: Figure size

    Returns:
        Dictionary with plot status and metadata
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np

        HAS_MATPLOTLIB = True
    except ImportError:
        HAS_MATPLOTLIB = False

    if not HAS_MATPLOTLIB:
        return {"status": "failed", "error": "matplotlib not available"}

    try:
        # Accept list-of-dicts: extract effect_sizes and p_values
        if results and isinstance(results[0], dict):
            effect_sizes: list[float] = [r.get("beta", 0.0) for r in results]
            p_values: list[float] = [r["p_value"] for r in results]
        else:
            effect_sizes = list(results)
            p_values = list(output_path) if output_path is not None else []
            output_path = None

        if len(effect_sizes) != len(p_values):
            return {"status": "failed", "error": "effect_sizes and p_values must have same length"}

        # Convert p-values to -log10
        neg_log_p = [-math.log10(max(p, 1e-300)) for p in p_values]

        # Create plot
        fig, ax = plt.subplots(figsize=figsize)

        # Color points based on significance and effect size
        colors = []
        for effect, p_val in zip(effect_sizes, p_values):
            if p_val < significance_threshold and abs(effect) > effect_size_threshold:
                colors.append("red")  # Significant and large effect
            elif p_val < significance_threshold:
                colors.append("orange")  # Significant but small effect
            elif abs(effect) > effect_size_threshold:
                colors.append("blue")  # Large effect but not significant
            else:
                colors.append("gray")  # Not significant, small effect

        ax.scatter(effect_sizes, neg_log_p, c=colors, alpha=0.6, s=20)

        # Add significance line
        sig_threshold = -math.log10(significance_threshold)
        ax.axhline(
            y=sig_threshold,
            color="red",
            linestyle="--",
            alpha=0.7,
            label=f"Significance threshold ({significance_threshold})",
        )

        # Add effect size lines
        ax.axvline(
            x=effect_size_threshold,
            color="blue",
            linestyle="--",
            alpha=0.7,
            label=f"Effect size threshold ({effect_size_threshold})",
        )
        ax.axvline(x=-effect_size_threshold, color="blue", linestyle="--", alpha=0.7)

        ax.set_xlabel("Effect Size (Beta)")
        ax.set_ylabel("-log₁₀(p-value)")
        ax.set_title("Volcano Plot")
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Add statistics text
        n_significant = sum(1 for p in p_values if p < significance_threshold)
        n_large_effect = sum(1 for effect in effect_sizes if abs(effect) > effect_size_threshold)
        n_both = sum(
            1
            for effect, p in zip(effect_sizes, p_values)
            if abs(effect) > effect_size_threshold and p < significance_threshold
        )

        stats_text = f"Total variants: {len(effect_sizes)}\nSignificant: {n_significant}\nLarge effect: {n_large_effect}\nBoth: {n_both}"
        ax.text(
            0.02,
            0.98,
            stats_text,
            transform=ax.transAxes,
            verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
        )

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches="tight")

        return {
            "status": "success",
            "n_variants": len(effect_sizes),
            "n_significant": n_significant,
            "n_large_effect": n_large_effect,
            "n_significant_large_effect": n_both,
            "significance_threshold": significance_threshold,
            "effect_size_threshold": effect_size_threshold,
            "output_path": str(output_path) if output_path else None,
        }

    except Exception as e:
        return {"status": "failed", "error": str(e)}
