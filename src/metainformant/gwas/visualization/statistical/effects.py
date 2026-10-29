"""GWAS effect size visualization utilities.

This module provides tools for visualizing and analyzing effect sizes
from GWAS results, including effect size distributions, comparisons,
and regional effect plots.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, Mapping, Optional

import numpy as np

from metainformant.core.utils import logging
from metainformant.gwas.analysis.correction import lambda_gc_from_p_values
from metainformant.gwas.visualization.utils import (
    chromosome_label,
    chromosome_sort_key,
    cumulative_chromosome_offsets,
    get_beta,
    get_chromosome,
    get_position,
    record_placement,
)

logger = logging.get_logger(__name__)


def effect_size_plot(results: Any, output_file: Optional[str | Path] = None) -> Optional[Any]:
    """Create effect size distribution and analysis plot.

    Args:
        results: GWAS results DataFrame or dictionary
        output_file: Optional output file path

    Returns:
        Plot object if matplotlib available, None otherwise

    Example:
        >>> plot = effect_size_plot(gwas_results)
    """
    try:
        import matplotlib.pyplot as plt
        import pandas as pd
        import seaborn as sns  # noqa: F401
    except ImportError:
        logger.warning("matplotlib, seaborn, or pandas not available for plotting")
        return None

    # Convert to DataFrame if needed
    if isinstance(results, (dict, list)):
        results = pd.DataFrame(results)

    if not hasattr(results, "columns") or "BETA" not in results.columns:
        logger.warning("Results must contain BETA column for effect size plotting")
        return None

    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle("GWAS Effect Size Analysis", fontsize=16)

    # 1. Effect size distribution
    ax1 = axes[0, 0]
    beta_values = results["BETA"].dropna()

    # Plot histogram
    ax1.hist(beta_values, bins=50, alpha=0.7, color="skyblue", edgecolor="black")
    ax1.axvline(x=0, color="red", linestyle="--", alpha=0.7, label="No effect")
    ax1.set_title("Effect Size Distribution")
    ax1.set_xlabel("Effect Size (Beta)")
    ax1.set_ylabel("Frequency")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # 2. Effect size by chromosome
    ax2 = axes[0, 1]
    if "CHR" in results.columns:
        # Group by chromosome
        chrom_effects = []
        chrom_labels = []

        for chrom in sorted(results["CHR"].unique()):
            chrom_data = results[results["CHR"] == chrom]["BETA"].dropna()
            if len(chrom_data) > 0:
                chrom_effects.append(chrom_data.values)
                chrom_labels.append(str(chrom))

        if chrom_effects:
            ax2.boxplot(chrom_effects, tick_labels=chrom_labels)
            ax2.set_title("Effect Size by Chromosome")
            ax2.set_xlabel("Chromosome")
            ax2.set_ylabel("Effect Size (Beta)")
            ax2.grid(True, alpha=0.3)

    # 3. Effect size vs p-value
    ax3 = axes[1, 0]
    if "P" in results.columns:
        p_values = results["P"].dropna()
        beta_subset = results["BETA"].dropna()

        # Align the data
        min_len = min(len(p_values), len(beta_subset))
        p_subset = p_values[:min_len]
        beta_aligned = beta_subset[:min_len]

        scatter = ax3.scatter(
            -np.log10(p_subset),
            beta_aligned,
            alpha=0.6,
            s=2,
            c=beta_aligned,
            cmap="RdYlBu_r",
        )
        ax3.set_title("Effect Size vs Significance")
        ax3.set_xlabel("-log10(P-value)")
        ax3.set_ylabel("Effect Size (Beta)")
        plt.colorbar(scatter, ax=ax3, label="Effect Size")

    # 4. Effect size summary statistics
    ax4 = axes[1, 1]
    ax4.axis("off")

    # Calculate statistics
    stats = {
        "Mean Effect": f"{beta_values.mean():.4f}",
        "Median Effect": f"{beta_values.median():.4f}",
        "Std Effect": f"{beta_values.std():.4f}",
        "Min Effect": f"{beta_values.min():.4f}",
        "Max Effect": f"{beta_values.max():.4f}",
        "Large Effects (>0.5)": int((abs(beta_values) > 0.5).sum()),
        "Total SNPs": len(beta_values),
    }

    # Create table
    table_data = [[key, value] for key, value in stats.items()]
    table = ax4.table(
        cellText=table_data,
        colLabels=["Statistic", "Value"],
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.2)
    ax4.set_title("Effect Size Summary", pad=20)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved effect size plot to {output_file}")

    return plt.gcf()


def compare_effect_sizes(study_results: Dict[str, Any], output_file: Optional[str | Path] = None) -> Optional[Any]:
    """Compare effect sizes across different studies or methods.

    Args:
        study_results: Dictionary mapping study names to results
        output_file: Optional output file path

    Returns:
        Plot object if matplotlib available, None otherwise

    Example:
        >>> studies = {"Study1": results1, "Study2": results2}
        >>> plot = compare_effect_sizes(studies)
    """
    try:
        import matplotlib.pyplot as plt
        import pandas as pd
        import seaborn as sns  # noqa: F401
    except ImportError:
        logger.warning("matplotlib, seaborn, or pandas not available for plotting")
        return None

    if not study_results:
        return None

    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle("Effect Size Comparison", fontsize=16)

    # Prepare data
    study_data = {}
    for study_name, results in study_results.items():
        if isinstance(results, dict):
            results = pd.DataFrame(results)

        if hasattr(results, "columns") and "BETA" in results.columns:
            study_data[study_name] = results["BETA"].dropna().values

    if not study_data:
        return None

    study_names = list(study_data.keys())

    # 1. Effect size distributions
    ax1 = axes[0, 0]
    for study_name, beta_values in study_data.items():
        ax1.hist(beta_values, bins=30, alpha=0.5, label=study_name, density=True)

    ax1.set_title("Effect Size Distributions")
    ax1.set_xlabel("Effect Size (Beta)")
    ax1.set_ylabel("Density")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # 2. Box plot comparison
    ax2 = axes[0, 1]
    beta_lists = list(study_data.values())
    ax2.boxplot(beta_lists, tick_labels=study_names)
    ax2.set_title("Effect Size Comparison")
    ax2.set_ylabel("Effect Size (Beta)")
    ax2.grid(True, alpha=0.3)

    # 3. Effect size correlation (if multiple studies)
    ax3 = axes[1, 0]
    if len(study_names) >= 2:
        # Compare first two studies
        study1_betas = study_data[study_names[0]]
        study2_betas = study_data[study_names[1]]

        # Align by SNP if possible (simplified)
        min_len = min(len(study1_betas), len(study2_betas))
        study1_subset = study1_betas[:min_len]
        study2_subset = study2_betas[:min_len]

        ax3.scatter(study1_subset, study2_subset, alpha=0.6, s=2)
        ax3.set_title(f"Effect Size Correlation: {study_names[0]} vs {study_names[1]}")
        ax3.set_xlabel(f"{study_names[0]} Effect Size")
        ax3.set_ylabel(f"{study_names[1]} Effect Size")

        # Add correlation line
        if len(study1_subset) > 1:
            corr = np.corrcoef(study1_subset, study2_subset)[0, 1]
            ax3.text(
                0.05,
                0.95,
                f"r = {corr:.3f}",
                transform=ax3.transAxes,
                fontsize=12,
                verticalalignment="top",
                bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
            )

        # Add diagonal line
        min_val = min(min(study1_subset), min(study2_subset))
        max_val = max(max(study1_subset), max(study2_subset))
        ax3.plot([min_val, max_val], [min_val, max_val], "r--", alpha=0.7)

    # 4. Effect size statistics comparison
    ax4 = axes[1, 1]
    ax4.axis("off")

    # Calculate statistics for each study
    stats_data = []
    for study_name, beta_values in study_data.items():
        stats_data.append([study_name, ".4f", ".4f", (abs(beta_values) > 0.5).sum(), len(beta_values)])

    if stats_data:
        table = ax4.table(
            cellText=stats_data,
            colLabels=[
                "Study",
                "Mean |Effect|",
                "Max |Effect|",
                "Large Effects",
                "Total SNPs",
            ],
            loc="center",
            cellLoc="center",
        )
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1, 1.2)

    ax4.set_title("Effect Size Statistics", pad=20)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved effect size comparison plot to {output_file}")

    return plt.gcf()


def effect_size_manhattan(
    results: Any,
    significance_threshold: float = 5e-8,
    output_file: Optional[str | Path] = None,
) -> Optional[Any]:
    """Create Manhattan plot colored by effect size.

    Args:
        results: GWAS results DataFrame or dictionary
        significance_threshold: P-value threshold for significance
        output_file: Optional output file path

    Returns:
        Plot object if matplotlib available, None otherwise

    Example:
        >>> plot = effect_size_manhattan(gwas_results)
    """
    try:
        import matplotlib.pyplot as plt
        import pandas as pd
    except ImportError:
        logger.warning("matplotlib or pandas not available for plotting")
        return None

    # Convert to DataFrame if needed
    if isinstance(results, dict):
        results = pd.DataFrame(results)

    required_cols = ["CHR", "BP", "P", "BETA"]
    if not all(col in results.columns for col in required_cols):
        logger.warning(f"Results must contain columns: {required_cols}")
        return None

    fig, ax = plt.subplots(1, 1, figsize=(15, 8))

    # Prepare data
    df = results.dropna(subset=required_cols).copy()
    df["-logP"] = -np.log10(df["P"].clip(lower=1e-300))

    # Chromosomes in natural order, laid end to end as in the other Manhattan plots: a variant's x is its
    # chromosome's offset plus its genomic position (not its row number within the chromosome)
    chrom_label = df["CHR"].map({value: chromosome_label(value) for value in df["CHR"].unique()})
    chromosomes = sorted(set(chrom_label), key=chromosome_sort_key)
    chrom_offsets = cumulative_chromosome_offsets(df["BP"].groupby(chrom_label).max().to_dict(), chromosomes)
    df["pos"] = df["BP"] + chrom_label.map(chrom_offsets)

    # Color by effect size
    effect_sizes = df["BETA"].abs()
    max_effect = effect_sizes.max()

    # Plot points colored by effect size
    scatter = ax.scatter(
        df["pos"],
        df["-logP"],
        c=effect_sizes,
        cmap="RdYlBu_r",
        alpha=0.8,
        s=3,
        vmin=0,
        vmax=max_effect,
    )

    # Add significance threshold line
    ax.axhline(
        y=-np.log10(significance_threshold),
        color="red",
        linestyle="--",
        alpha=0.7,
        label=f"Significance ({significance_threshold:.0e})",
    )

    # Add chromosome labels at the centre of each chromosome's variants
    for chrom in chromosomes:
        chrom_x = df.loc[chrom_label == chrom, "pos"]
        ax.text(
            (chrom_x.min() + chrom_x.max()) / 2,
            ax.get_ylim()[0] - 0.5,
            chrom,
            ha="center",
            va="top",
            fontsize=10,
        )

    ax.set_title("Manhattan Plot (Colored by Effect Size)", fontsize=14)
    ax.set_xlabel("Chromosome", fontsize=12)
    ax.set_ylabel("-log10(P-value)", fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label("Absolute Effect Size", fontsize=12)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved effect size Manhattan plot to {output_file}")

    return plt.gcf()


def effect_size_qq(results: Any, output_file: Optional[str | Path] = None) -> Optional[Any]:
    """Create Q-Q plot with effect size information.

    Args:
        results: GWAS results DataFrame or dictionary
        output_file: Optional output file path

    Returns:
        Plot object if matplotlib available, None otherwise

    Example:
        >>> plot = effect_size_qq(gwas_results)
    """
    try:
        import matplotlib.pyplot as plt
        import pandas as pd
    except ImportError:
        logger.warning("matplotlib or pandas not available for plotting")
        return None

    # Convert to DataFrame if needed
    if isinstance(results, dict):
        results = pd.DataFrame(results)

    if "P" not in results.columns or "BETA" not in results.columns:
        logger.warning("Results must contain P and BETA columns")
        return None

    fig, ax = plt.subplots(1, 1, figsize=(10, 8))

    # Prepare Q-Q data
    p_values = results["P"].dropna().values
    p_values = np.sort(p_values)

    # Expected p-values under null
    n = len(p_values)
    expected = np.arange(1, n + 1) / (n + 1)

    # Convert to -log10
    observed_log = -np.log10(p_values.clip(lower=1e-300))
    expected_log = -np.log10(expected)

    # Color by effect size
    beta_values = results["BETA"].dropna().values
    beta_sorted = beta_values[np.argsort(results["P"].dropna().values)]
    effect_sizes = np.abs(beta_sorted)

    # Plot
    scatter = ax.scatter(expected_log, observed_log, c=effect_sizes, cmap="RdYlBu_r", alpha=0.7, s=4)

    # Add diagonal line
    ax.plot(
        [0, max(expected_log.max(), observed_log.max())],
        [0, max(expected_log.max(), observed_log.max())],
        "k--",
        alpha=0.7,
        label="Expected",
    )

    ax.set_title("Q-Q Plot (Colored by Effect Size)", fontsize=14)
    ax.set_xlabel("Expected -log10(P)", fontsize=12)
    ax.set_ylabel("Observed -log10(P)", fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label("Absolute Effect Size", fontsize=12)

    # Add inflation factor annotation
    lambda_gc = lambda_gc_from_p_values(p_values.tolist()) or 1.0
    ax.text(
        0.05,
        0.95,
        f"λ = {lambda_gc:.3f}",
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
    )

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved effect size Q-Q plot to {output_file}")

    return plt.gcf()


def effect_size_forest_plot(
    effect_data: Any,
    output_file: Optional[str | Path] = None,
    title: str = "Effect Size Forest Plot",
    top_n: Optional[int] = None,
) -> dict[str, Any]:
    """Create a forest plot showing effect sizes and confidence intervals.

    Args:
        effect_data: Either a dictionary containing effect size data with keys
                    ('variants', 'effects', 'lower_ci', 'upper_ci'), or a list
                    of dicts with keys like CHROM, POS, BETA, SE, p_value.
        output_file: Optional output file path
        title: Plot title
        top_n: Optional limit on number of variants to show (sorted by significance)

    Returns:
        Dictionary with plot status and metadata

    Example:
        >>> data = {
        ...     'variants': ['rs123', 'rs456', 'rs789'],
        ...     'effects': [0.1, 0.05, 0.15],
        ...     'lower_ci': [0.05, 0.02, 0.08],
        ...     'upper_ci': [0.15, 0.08, 0.22]
        ... }
        >>> result = effect_size_forest_plot(data)
    """
    try:
        import matplotlib.patches as patches  # noqa: F401
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for forest plot")
        return {"status": "failed", "error": "matplotlib not available"}

    try:
        # Handle list-of-dicts input (GWAS results format)
        if isinstance(effect_data, list):
            if not effect_data:
                return {"status": "failed", "error": "No effect data provided"}

            # Sort by p-value (most significant first) and apply top_n
            sorted_data = sorted(
                effect_data,
                key=lambda x: x.get("p_value", x.get("P", x.get("pval", 1.0))),
            )
            if top_n is not None:
                sorted_data = sorted_data[:top_n]

            variants = []
            effects = []
            lower_ci = []
            upper_ci = []

            for item in sorted_data:
                chrom = get_chromosome(item, "?")
                pos = item.get("POS", item.get("pos", 0))
                beta = item.get("BETA", item.get("beta", 0.0))
                se = item.get("SE", item.get("se", abs(beta) * 0.3))

                variant_name = f"chr{chrom}:{pos}"
                variants.append(variant_name)
                effects.append(float(beta))
                lower_ci.append(float(beta) - 1.96 * float(se))
                upper_ci.append(float(beta) + 1.96 * float(se))
        else:
            # Legacy dict-of-lists format
            required_keys = ["variants", "effects", "lower_ci", "upper_ci"]
            missing_keys = [key for key in required_keys if key not in effect_data]
            if missing_keys:
                return {
                    "status": "failed",
                    "error": f"Missing required keys: {missing_keys}",
                }

            variants = effect_data["variants"]
            effects = effect_data["effects"]
            lower_ci = effect_data["lower_ci"]
            upper_ci = effect_data["upper_ci"]

            if not (len(variants) == len(effects) == len(lower_ci) == len(upper_ci)):
                return {
                    "status": "failed",
                    "error": "All effect data arrays must have the same length",
                }

            if top_n is not None:
                variants = variants[:top_n]
                effects = effects[:top_n]
                lower_ci = lower_ci[:top_n]
                upper_ci = upper_ci[:top_n]

        if not variants:
            return {"status": "failed", "error": "No variants to plot"}

        fig, ax = plt.subplots(1, 1, figsize=(10, max(6, len(variants) * 0.3)))

        # Plot effect sizes and confidence intervals
        y_positions = range(len(variants))

        # Plot confidence intervals
        for i, (effect, lower, upper, variant) in enumerate(zip(effects, lower_ci, upper_ci, variants)):
            # Confidence interval line
            ax.plot([lower, upper], [i, i], "k-", linewidth=2, alpha=0.8)

            # Effect size point
            ax.plot(effect, i, "ro", markersize=8, alpha=0.8)

            # Add variant label
            ax.text(-0.02, i, variant, ha="right", va="center", fontsize=10)

        # Add vertical line at zero effect
        ax.axvline(x=0, color="gray", linestyle="--", alpha=0.7)

        # Add vertical line at reference effect (e.g., 1 for OR, 0 for beta)
        if all(isinstance(e, (int, float)) and e > 0 for e in effects):
            # Likely odds ratios - add reference line at 1
            ax.axvline(x=1, color="blue", linestyle="-", alpha=0.5, label="Reference (OR=1)")
            ax.legend()

        # Format plot
        ax.set_xlabel("Effect Size", fontsize=12)
        ax.set_title(title, fontsize=14, pad=20)
        ax.set_yticks(y_positions)
        ax.set_yticklabels([""] * len(variants))  # Hide y-tick labels (we use text labels)
        ax.grid(True, alpha=0.3, axis="x")

        # Add summary statistics
        mean_effect = np.mean(effects)
        ax.text(
            0.02,
            0.98,
            f"Mean Effect: {mean_effect:.3f}",
            transform=ax.transAxes,
            fontsize=11,
            verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
        )

        plt.tight_layout()

        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches="tight")
            logger.info(f"Saved forest plot to {output_file}")

        return {
            "status": "success",
            "n_variants": len(variants),
            "output_path": str(output_file) if output_file else None,
        }

    except Exception as e:
        logger.error(f"Error creating forest plot: {e}")
        return {"status": "failed", "error": str(e)}


def effect_direction_plot(
    results: list[dict],
    output_path: Optional[str | Path] = None,
    figsize: tuple[int, int] = (10, 6),
) -> dict[str, Any]:
    """Create an effect direction plot showing beta coefficients by chromosome.

    Chromosome, position and effect are read under every spelling of
    :mod:`metainformant.gwas.visualization.utils` (``CHROM``/``chrom``/..., ``POS``/``pos``/``BP``/...,
    ``BETA``/``beta``). A result without an effect estimate (NaN, None or ``"NA"``) is not drawn, rather than
    drawn as an effect of 0, and ``n_without_effect`` reports how many were left out. A result without a
    chromosome or position fails the plot with an error that names the accepted keys.

    Args:
        results: GWAS results, each with a chromosome, a position and a beta
        output_path: Path to save the plot
        figsize: Figure size

    Returns:
        Dictionary with plot status and metadata
    """
    try:
        import matplotlib.pyplot as plt

        HAS_MATPLOTLIB = True
    except ImportError:
        HAS_MATPLOTLIB = False

    if not HAS_MATPLOTLIB:
        return {"status": "failed", "error": "matplotlib not available"}

    try:
        if not results:
            return {"status": "failed", "error": "No results provided"}

        # Extract data
        chromosomes: list[str] = []
        positions: list[float] = []
        betas: list[float] = []
        n_without_effect = 0

        for index, result in enumerate(results):
            chrom, pos = record_placement(result, index, plot="effect_direction_plot")
            beta = get_beta(result)
            if math.isnan(beta):
                n_without_effect += 1
                continue

            chromosomes.append(chrom)
            positions.append(pos)
            betas.append(beta)

        if n_without_effect:
            logger.warning(
                f"{n_without_effect} of {len(results)} results have no effect estimate and are not drawn "
                "in the effect direction plot"
            )
        if not betas:
            return {"status": "failed", "error": "No beta values found"}

        # Chromosomes in natural order (1, 2, ..., 10, ..., X), laid end to end: a variant's x position is
        # its chromosome's offset plus its position on that chromosome
        unique_chroms = sorted(set(chromosomes), key=chromosome_sort_key)
        chrom_extent: dict[str, float] = {}
        for chrom, pos in zip(chromosomes, positions):
            chrom_extent[chrom] = max(chrom_extent.get(chrom, pos), pos)
        chrom_offsets = cumulative_chromosome_offsets(chrom_extent, unique_chroms)
        x_positions = [chrom_offsets[chrom] + pos for chrom, pos in zip(chromosomes, positions)]

        # Create plot
        fig, ax = plt.subplots(figsize=figsize)

        # Color points by effect direction
        colors = ["red" if beta < 0 else "blue" for beta in betas]
        ax.scatter(x_positions, betas, c=colors, alpha=0.6, s=20)

        # Add zero line
        ax.axhline(y=0, color="black", linestyle="-", alpha=0.5)

        # Add chromosome labels at the mean position of each chromosome's variants
        chrom_centers = []
        for chrom in unique_chroms:
            chrom_x = [x for x, c in zip(x_positions, chromosomes) if c == chrom]
            if chrom_x:
                chrom_centers.append((chrom, sum(chrom_x) / len(chrom_x)))

        for chrom, center in chrom_centers:
            ax.text(center, min(betas) - 0.1, str(chrom), ha="center", va="top", fontsize=8)

        ax.set_xlabel("Genomic Position")
        ax.set_ylabel("Effect Size (Beta)")
        ax.set_title("Effect Direction Plot")
        ax.grid(True, alpha=0.3)

        # Add legend
        import matplotlib.patches as mpatches

        red_patch = mpatches.Patch(color="red", label="Negative effect", alpha=0.6)
        blue_patch = mpatches.Patch(color="blue", label="Positive effect", alpha=0.6)
        ax.legend(handles=[red_patch, blue_patch])

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches="tight")

        return {
            "status": "success",
            "n_variants": len(betas),
            "n_chromosomes": len(unique_chroms),
            "positive_effects": sum(1 for b in betas if b > 0),
            "negative_effects": sum(1 for b in betas if b < 0),
            "n_without_effect": n_without_effect,
            "output_path": str(output_path) if output_path else None,
        }

    except Exception as e:
        return {"status": "failed", "error": str(e)}


#: Variant-identifier keys searched when annotations are given as a mapping (the keys ``allelic_series_plot`` reads).
_VARIANT_ID_KEYS = ("ID", "rsid", "variant_id")


def _clean_category(value: Any) -> Optional[str]:
    """A supplied functional category as text, or None when it is absent (None, NaN or blank)."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    text = str(value).strip()
    return text or None


def _normalise_annotation_keys(annotations: Mapping[Any, Any]) -> Dict[Any, Any]:
    """Copy of ``annotations`` whose ``(chromosome, position)`` keys are normalised (``1``, ``1.0``, ``"1"`` agree)."""
    normalised: Dict[Any, Any] = {}
    for key, category in annotations.items():
        if isinstance(key, tuple) and len(key) == 2:
            normalised[(chromosome_label(key[0]), int(key[1]))] = category
        else:
            normalised[key] = category
    return normalised


def _variant_category(variant: Mapping[str, Any], annotation_key: str, annotations: Mapping[Any, Any]) -> Optional[str]:
    """The functional category the caller supplied for ``variant``, or None. It is never guessed.

    The variant's own ``annotation_key`` field wins; otherwise ``annotations`` (keys already normalised) is
    searched by the variant's id (``ID``, ``rsid`` or ``variant_id``) and then by ``(chromosome, position)``.
    """
    category = _clean_category(variant.get(annotation_key))
    if category is not None:
        return category
    for key in _VARIANT_ID_KEYS:
        identifier = variant.get(key)
        if identifier is not None and identifier in annotations:
            return _clean_category(annotations[identifier])
    chrom = get_chromosome(variant)
    position = get_position(variant)
    if chrom is not None and position is not None:
        try:
            location = (chromosome_label(chrom), int(position))
        except (TypeError, ValueError):  # e.g. a NaN position: the variant simply has no (chromosome, position) key
            return None
        return _clean_category(annotations.get(location))
    return None


def functional_enrichment_plot(
    enrichment_results: Any,
    output_path: Optional[str | Path] = None,
    top_n: int = 20,
    figsize: tuple[int, int] = (12, 8),
    significance_threshold: Optional[float] = None,
    annotations: Optional[Mapping[Any, Any]] = None,
    annotation_key: str = "annotation",
) -> dict[str, Any]:
    """Create a functional enrichment plot for GWAS results.

    Args:
        enrichment_results: Either a dictionary mapping functions to enrichment
            statistics, or a list of dicts with GWAS variant results (CHROM,
            POS, p_value) whose functional categories you supply: as a field on each
            variant (``annotation_key``) or through ``annotations``. The function never
            assigns categories itself; the bars count the significant variants in each
            category you gave.
        output_path: Path to save the plot
        top_n: Number of top enriched functions to show
        figsize: Figure size
        significance_threshold: P-value threshold for significance (default 5e-8).
            Must be between 0 and 1 (exclusive).
        annotations: Optional mapping from a variant to its category, for variants that carry no
            ``annotation_key`` field. Keys are a variant id (the record's ``ID``, ``rsid`` or ``variant_id``)
            or a ``(chromosome, position)`` tuple (``1``, ``1.0`` and ``"1"`` name the same chromosome).
        annotation_key: Name of the per-variant category field (default ``"annotation"``). A blank, None or
            NaN value counts as no annotation.

    Returns:
        Dictionary with plot status and metadata. For a list of variants: ``category_counts``
        (significant variants per supplied category), ``categories_found``,
        ``total_significant_variants`` and ``n_unannotated`` (significant variants with no category, which
        are not plotted). When no significant variant has a category the status is ``"skipped"`` and no
        figure is drawn.
    """
    # Validate significance_threshold if provided
    if significance_threshold is not None:
        if significance_threshold <= 0 or significance_threshold >= 1:
            raise ValueError("significance_threshold must be between 0 and 1 (exclusive)")

    threshold = significance_threshold if significance_threshold is not None else 5e-8

    try:
        import matplotlib.pyplot as plt

        HAS_MATPLOTLIB = True
    except ImportError:
        HAS_MATPLOTLIB = False

    if not HAS_MATPLOTLIB:
        return {"status": "failed", "error": "matplotlib not available"}

    try:
        if not enrichment_results:
            return {"status": "failed", "error": "No results provided"}

        # Handle list-of-dicts input (GWAS variant results)
        if isinstance(enrichment_results, list):
            # Filter significant variants
            significant_variants = [
                v for v in enrichment_results if v.get("p_value", v.get("P", v.get("pval", 1.0))) < threshold
            ]

            if not significant_variants:
                return {
                    "status": "failed",
                    "error": "No significant variants found at threshold {:.1e}".format(threshold),
                }

            # Count the categories the caller supplied. Nothing is inferred from the variant's position or
            # anything else: a variant without an annotation is reported as such, not given a category.
            lookup = _normalise_annotation_keys(annotations) if annotations is not None else {}
            category_counts: Dict[str, int] = {}
            for variant in significant_variants:
                category = _variant_category(variant, annotation_key, lookup)
                if category is not None:
                    category_counts[category] = category_counts.get(category, 0) + 1
            n_unannotated = len(significant_variants) - sum(category_counts.values())

            if not category_counts:
                return {
                    "status": "skipped",
                    "message": (
                        f"No functional annotations were supplied for the {len(significant_variants)} significant "
                        f"variants, so no category plot was drawn: give each variant a '{annotation_key}' field "
                        "or pass annotations={variant id or (chromosome, position): category}"
                    ),
                    "total_significant_variants": len(significant_variants),
                    "n_unannotated": n_unannotated,
                }

            # Create plot
            fig, ax = plt.subplots(figsize=figsize)

            sorted_cats = sorted(category_counts.items(), key=lambda x: x[1], reverse=True)[:top_n]
            cat_names = [c[0] for c in sorted_cats]
            cat_values = [c[1] for c in sorted_cats]

            bars = ax.barh(range(len(cat_names)), cat_values, color="skyblue", edgecolor="black")
            ax.set_yticks(range(len(cat_names)))
            ax.set_yticklabels(cat_names)
            ax.set_xlabel("Number of Significant Variants")
            ax.set_title("Functional Annotation of Significant GWAS Variants")
            ax.grid(True, alpha=0.3)
            if n_unannotated:
                ax.text(
                    0.99,
                    0.02,
                    f"{n_unannotated} significant variant(s) without an annotation are not shown",
                    transform=ax.transAxes,
                    ha="right",
                    va="bottom",
                    fontsize=9,
                )

            plt.tight_layout()

            if output_path:
                plt.savefig(output_path, dpi=300, bbox_inches="tight")

            return {
                "status": "success",
                "total_significant_variants": len(significant_variants),
                "categories_found": len(category_counts),
                "category_counts": category_counts,
                "n_unannotated": n_unannotated,
                "output_path": str(output_path) if output_path else None,
            }

        # Legacy dict-of-dicts format (enrichment results)
        if not enrichment_results:
            return {"status": "failed", "error": "No enrichment results provided"}

        # Sort by enrichment ratio and get top N
        sorted_results = sorted(
            [(func, stats) for func, stats in enrichment_results.items() if stats.get("enrichment_ratio", 0) > 1],
            key=lambda x: x[1].get("enrichment_ratio", 0),
            reverse=True,
        )[:top_n]

        if not sorted_results:
            return {"status": "failed", "error": "No significantly enriched functions"}

        functions = [func for func, _ in sorted_results]
        enrichment_ratios = [stats.get("enrichment_ratio", 1) for _, stats in sorted_results]
        p_values = [stats.get("p_value", 1) for _, stats in sorted_results]

        # Create plot
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

        # Plot enrichment ratios
        bars = ax1.barh(range(len(functions)), enrichment_ratios, color="skyblue")
        ax1.set_yticks(range(len(functions)))
        ax1.set_yticklabels([f[:30] + "..." if len(f) > 30 else f for f in functions])
        ax1.set_xlabel("Enrichment Ratio")
        ax1.set_title("Functional Enrichment Ratios")
        ax1.grid(True, alpha=0.3)

        # Add enrichment ratio values on bars
        for i, (bar, ratio) in enumerate(zip(bars, enrichment_ratios)):
            ax1.text(
                bar.get_width() + 0.1,
                bar.get_y() + bar.get_height() / 2,
                f"{ratio:.1f}",
                ha="left",
                va="center",
                fontsize=8,
            )

        # Plot -log10(p-values)
        import math

        neg_log_p = [-math.log10(max(p, 1e-10)) for p in p_values]
        ax2.scatter(enrichment_ratios, neg_log_p, s=50, alpha=0.7, color="red")

        # Add function labels
        for i, func in enumerate(functions):
            short_name = func[:20] + "..." if len(func) > 20 else func
            ax2.annotate(
                short_name,
                (enrichment_ratios[i], neg_log_p[i]),
                xytext=(5, 5),
                textcoords="offset points",
                fontsize=6,
            )

        ax2.set_xlabel("Enrichment Ratio")
        ax2.set_ylabel("-log10(p-value)")
        ax2.set_title("Enrichment vs Significance")
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches="tight")

        return {
            "status": "success",
            "n_functions": len(functions),
            "top_functions": functions[:5],  # Return top 5 for summary
            "enrichment_range": f"{min(enrichment_ratios):.1f}-{max(enrichment_ratios):.1f}",
            "output_path": str(output_path) if output_path else None,
        }

    except Exception as e:
        return {"status": "failed", "error": str(e)}


def allelic_series_plot(
    results: list[dict],
    output_path: Optional[str | Path] = None,
    gene_region: Optional[tuple[str, int, int]] = None,
) -> dict[str, Any]:
    """Create an allelic series plot showing effect sizes across alleles at a locus.

    This visualizes multiple variants within a gene/region, ordered by position,
    with effect sizes shown as a bar chart to identify allelic series patterns.

    Args:
        results: GWAS results with CHROM, POS, BETA columns
        output_path: Path to save the plot
        gene_region: Optional tuple of (chrom, start, end) for gene region

    Returns:
        Dictionary with plot status and metadata
    """
    try:
        import matplotlib.pyplot as plt

        HAS_MATPLOTLIB = True
    except ImportError:
        HAS_MATPLOTLIB = False

    if not HAS_MATPLOTLIB:
        return {"status": "failed", "error": "matplotlib not available"}

    if not results:
        return {"status": "failed", "error": "No results provided"}

    try:
        # Check if any results have effect size data
        has_beta = any(
            result.get("BETA") is not None or result.get("beta") is not None or result.get("effect_size") is not None
            for result in results
        )
        if not has_beta:
            return {
                "status": "skipped",
                "message": "No effect size data (BETA/beta/effect_size) available in results",
            }

        # Extract and filter variants
        variants = []
        for result in results:
            chrom = get_chromosome(result)
            pos = result.get("POS", result.get("pos", result.get("position")))
            beta = result.get("BETA", result.get("beta", result.get("effect_size")))
            rsid = result.get("ID", result.get("rsid", result.get("variant_id", f"var_{pos}")))
            p_value = result.get("P", result.get("p_value", result.get("pval", 1.0)))

            if beta is None or pos is None:
                continue

            # Filter by gene region if specified
            if gene_region:
                region_chrom, region_start, region_end = gene_region
                if str(chrom) != str(region_chrom):
                    continue
                if not (region_start <= pos <= region_end):
                    continue

            variants.append(
                {
                    "chrom": chrom,
                    "pos": pos,
                    "beta": float(beta),
                    "rsid": rsid,
                    "p_value": float(p_value),
                }
            )

        if not variants:
            return {
                "status": "failed",
                "error": "No variants found in specified region",
            }

        # Sort by position
        variants.sort(key=lambda x: x["pos"])

        # Create the plot
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
        fig.suptitle("Allelic Series Plot", fontsize=14)

        [v["pos"] for v in variants]
        betas = [v["beta"] for v in variants]
        p_values = [v["p_value"] for v in variants]
        rsids = [v["rsid"] for v in variants]

        # Top panel: Effect sizes
        colors = ["red" if b < 0 else "blue" for b in betas]
        ax1.bar(range(len(variants)), betas, color=colors, alpha=0.7, edgecolor="black")

        ax1.axhline(y=0, color="black", linestyle="-", alpha=0.5)
        ax1.set_ylabel("Effect Size (Beta)", fontsize=12)
        ax1.set_title("Effect Sizes by Variant", fontsize=12)
        ax1.grid(True, alpha=0.3, axis="y")

        # Add legend
        import matplotlib.patches as mpatches

        red_patch = mpatches.Patch(color="red", label="Negative effect", alpha=0.7)
        blue_patch = mpatches.Patch(color="blue", label="Positive effect", alpha=0.7)
        ax1.legend(handles=[red_patch, blue_patch], loc="upper right")

        # Bottom panel: -log10(p-value)
        neg_log_p = [-np.log10(max(p, 1e-300)) for p in p_values]
        ax2.bar(range(len(variants)), neg_log_p, color="green", alpha=0.7, edgecolor="black")
        ax2.axhline(y=-np.log10(5e-8), color="red", linestyle="--", alpha=0.7, label="p = 5e-8")
        ax2.set_ylabel("-log₁₀(p-value)", fontsize=12)
        ax2.set_xlabel("Variant", fontsize=12)
        ax2.set_title("Significance by Variant", fontsize=12)
        ax2.grid(True, alpha=0.3, axis="y")
        ax2.legend()

        # Set x-axis labels
        ax2.set_xticks(range(len(variants)))
        ax2.set_xticklabels(rsids, rotation=45, ha="right", fontsize=8)

        # Add summary statistics
        positive_effects = sum(1 for b in betas if b > 0)
        negative_effects = sum(1 for b in betas if b < 0)
        significant = sum(1 for p in p_values if p < 5e-8)

        stats_text = (
            f"Variants: {len(variants)} | Positive: {positive_effects} | "
            f"Negative: {negative_effects} | Significant: {significant}"
        )
        fig.text(
            0.5,
            0.02,
            stats_text,
            ha="center",
            fontsize=10,
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
        )

        plt.tight_layout()
        plt.subplots_adjust(bottom=0.15)

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches="tight")
            logger.info(f"Saved allelic series plot to {output_path}")

        return {
            "status": "success",
            "n_variants": len(variants),
            "positive_effects": positive_effects,
            "negative_effects": negative_effects,
            "significant_variants": significant,
            "gene_region": gene_region,
            "output_path": str(output_path) if output_path else None,
        }

    except Exception as e:
        logger.error(f"Error creating allelic series plot: {e}")
        return {"status": "failed", "error": str(e)}
