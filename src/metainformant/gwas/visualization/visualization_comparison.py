"""GWAS comparison visualization utilities.

This module provides tools for comparing GWAS results across different
studies, populations, or analysis methods.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from metainformant.core import logging

logger = logging.get_logger(__name__)


def compare_gwas_studies(
    study_results: Dict[str, Dict[str, Any]], output_file: Optional[str | Path] = None
) -> Optional[Any]:
    """Create comparison plots across multiple GWAS studies.

    Args:
        study_results: Dictionary mapping study names to GWAS results
        output_file: Optional output file path

    Returns:
        Plot object if matplotlib available, None otherwise

    Example:
        >>> studies = {
        ...     "Study1": {"results": df1, "significance_threshold": 5e-8},
        ...     "Study2": {"results": df2, "significance_threshold": 5e-8}
        ... }
        >>> plot = compare_gwas_studies(studies)
    """
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
        import pandas as pd
    except ImportError:
        logger.warning("matplotlib, seaborn, or pandas not available for plotting")
        return None

    if not study_results:
        return None

    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle("GWAS Study Comparison", fontsize=16)

    # Prepare data
    study_names = list(study_results.keys())
    colors = sns.color_palette("husl", len(study_names))

    # 1. Manhattan plot comparison
    ax1 = axes[0, 0]
    for i, (study_name, study_data) in enumerate(study_results.items()):
        results = study_data.get("results")
        if isinstance(results, dict):
            # Convert dict to DataFrame if needed
            import pandas as pd

            results = pd.DataFrame(results)

        if hasattr(results, "columns") and "CHR" in results.columns and "P" in results.columns:
            plot_manhattan_overlay(ax1, results, study_name, colors[i])

    ax1.set_title("Manhattan Plot Comparison")
    ax1.set_xlabel("Chromosome")
    ax1.set_ylabel("-log10(P-value)")

    # 2. QQ plot comparison
    ax2 = axes[0, 1]
    for i, (study_name, study_data) in enumerate(study_results.items()):
        results = study_data.get("results")
        if isinstance(results, dict):
            import pandas as pd

            results = pd.DataFrame(results)

        if hasattr(results, "columns") and "P" in results.columns:
            p_values = results["P"].dropna().values
            plot_qq_overlay(ax2, p_values, study_name, colors[i])

    ax2.set_title("Q-Q Plot Comparison")
    ax2.plot([0, 1], [0, 1], "k--", alpha=0.5, label="Expected")

    # 3. Effect size distribution
    ax3 = axes[1, 0]
    effect_sizes = []
    labels = []

    for study_name, study_data in study_results.items():
        results = study_data.get("results")
        if isinstance(results, dict):
            import pandas as pd

            results = pd.DataFrame(results)

        if hasattr(results, "columns") and "BETA" in results.columns:
            beta_values = results["BETA"].dropna().values
            effect_sizes.append(beta_values)
            labels.append(study_name)

    if effect_sizes:
        ax3.hist(effect_sizes, bins=50, alpha=0.7, label=labels)
        ax3.set_title("Effect Size Distribution")
        ax3.set_xlabel("Effect Size (Beta)")
        ax3.set_ylabel("Frequency")
        ax3.legend()

    # 4. Summary statistics table
    ax4 = axes[1, 1]
    ax4.axis("off")

    # Create summary table
    summary_data = []
    for study_name, study_data in study_results.items():
        results = study_data.get("results")
        if isinstance(results, dict):
            import pandas as pd

            results = pd.DataFrame(results)

        n_snps = len(results) if hasattr(results, "__len__") else 0
        sig_threshold = study_data.get("significance_threshold", 5e-8)

        if hasattr(results, "columns") and "P" in results.columns:
            n_sig = (results["P"] < sig_threshold).sum()
        else:
            n_sig = 0

        summary_data.append([study_name, n_snps, n_sig])

    if summary_data:
        table = ax4.table(
            cellText=summary_data, colLabels=["Study", "Total SNPs", "Significant SNPs"], loc="center", cellLoc="center"
        )
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 1.5)

    ax4.set_title("Study Summary", pad=20)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved GWAS comparison plot to {output_file}")

    return plt.gcf()


def plot_manhattan_overlay(ax: Any, results: Any, label: str, color: str) -> None:
    """Plot Manhattan plot for one study on given axis."""
    try:
        import pandas as pd
    except ImportError:
        return

    if not hasattr(results, "columns") or "CHR" not in results.columns or "P" not in results.columns:
        return

    # Prepare data
    df = results.copy()
    df["-logP"] = -np.log10(df["P"].clip(lower=1e-300))  # Avoid log(0)

    # Create chromosome positions
    chromosomes = sorted(df["CHR"].unique())
    chrom_starts = {}
    current_pos = 0

    for chrom in chromosomes:
        chrom_starts[chrom] = current_pos
        chrom_size = len(df[df["CHR"] == chrom])
        current_pos += chrom_size + 1000  # Gap between chromosomes

    df["pos"] = df.apply(lambda row: chrom_starts[row["CHR"]] + row.name, axis=1)

    # Plot
    ax.scatter(df["pos"], df["-logP"], alpha=0.6, s=2, color=color, label=label)

    # Add chromosome labels
    for chrom in chromosomes:
        chrom_center = chrom_starts[chrom] + len(df[df["CHR"] == chrom]) / 2
        ax.text(chrom_center, ax.get_ylim()[0] - 0.5, str(chrom), ha="center", va="top", fontsize=8)

    ax.legend()


def plot_qq_overlay(ax: Any, p_values: np.ndarray, label: str, color: str) -> None:
    """Plot Q-Q plot for one study on given axis."""
    # Remove NA values and sort
    p_values = p_values[~np.isnan(p_values)]
    p_values = np.sort(p_values)

    # Expected p-values under null
    n = len(p_values)
    expected = np.arange(1, n + 1) / (n + 1)

    # Convert to -log10
    observed_log = -np.log10(p_values.clip(lower=1e-300))
    expected_log = -np.log10(expected)

    # Plot
    ax.scatter(expected_log, observed_log, alpha=0.6, s=2, color=color, label=label)


def compare_populations(
    pop_data: Dict[str, Any], trait_name: str = "Trait", output_file: Optional[str | Path] = None
) -> Optional[Any]:
    """Compare GWAS results across different populations.

    Args:
        pop_data: Dictionary mapping population names to GWAS results
        trait_name: Name of the trait being analyzed
        output_file: Optional output file path

    Returns:
        Plot object if matplotlib available, None otherwise

    Example:
        >>> populations = {
        ...     "EUR": {"results": eur_df, "color": "blue"},
        ...     "AFR": {"results": afr_df, "color": "red"}
        ... }
        >>> plot = compare_populations(populations, "Height")
    """
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        logger.warning("matplotlib or seaborn not available for plotting")
        return None

    if not pop_data:
        return None

    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle(f"{trait_name} - Population Comparison", fontsize=16)

    # 1. Effect size comparison
    ax1 = axes[0, 0]
    pop_names = []
    effect_sizes = []

    for pop_name, pop_info in pop_data.items():
        results = pop_info.get("results")
        if hasattr(results, "columns") and "BETA" in results.columns:
            beta_values = results["BETA"].dropna().values
            effect_sizes.append(beta_values)
            pop_names.append(pop_name)

    if effect_sizes:
        ax1.boxplot(effect_sizes, labels=pop_names)
        ax1.set_title("Effect Size Distribution by Population")
        ax1.set_ylabel("Effect Size (Beta)")
        ax1.grid(True, alpha=0.3)

    # 2. Significant SNP overlap
    ax2 = axes[0, 1]
    sig_snps = {}

    threshold = 5e-8  # Common GWAS threshold

    for pop_name, pop_info in pop_data.items():
        results = pop_info.get("results")
        if hasattr(results, "columns") and "P" in results.columns:
            sig_mask = results["P"] < threshold
            sig_snps[pop_name] = set(results[sig_mask].index if hasattr(results, "index") else [])

    if len(sig_snps) >= 2:
        # Create overlap matrix
        pops = list(sig_snps.keys())
        overlap_matrix = np.zeros((len(pops), len(pops)))

        for i, pop1 in enumerate(pops):
            for j, pop2 in enumerate(pops):
                if i == j:
                    overlap_matrix[i, j] = len(sig_snps[pop1])
                else:
                    overlap = len(sig_snps[pop1] & sig_snps[pop2])
                    overlap_matrix[i, j] = overlap

        sns.heatmap(overlap_matrix, annot=True, fmt=".0f", cmap="Blues", xticklabels=pops, yticklabels=pops, ax=ax2)
        ax2.set_title("Significant SNP Overlap")

    # 3. Population-specific enrichment
    ax3 = axes[1, 0]
    enrichment_data = calculate_population_enrichment(pop_data)

    if enrichment_data:
        pops = list(enrichment_data.keys())
        enrichments = [enrichment_data[pop]["enrichment"] for pop in pops]
        p_values = [enrichment_data[pop]["p_value"] for pop in pops]

        bars = ax3.bar(pops, enrichments, alpha=0.7)
        ax3.set_title("Population-Specific Enrichment")
        ax3.set_ylabel("Enrichment Score")
        ax3.set_yscale("log")

        # Add p-value annotations
        for bar, p_val in zip(bars, p_values):
            height = bar.get_height()
            ax3.text(bar.get_x() + bar.get_width() / 2, height, f"p={p_val:.2e}", ha="center", va="bottom", fontsize=8)

    # 4. Genetic architecture comparison
    ax4 = axes[1, 1]
    arch_data = analyze_genetic_architecture(pop_data)

    if arch_data:
        pops = list(arch_data.keys())
        h2_estimates = [arch_data[pop].get("heritability", 0) for pop in pops]
        polygenicity = [arch_data[pop].get("polygenicity", 0) for pop in pops]

        x = np.arange(len(pops))
        width = 0.35

        ax4.bar(x - width / 2, h2_estimates, width, label="Heritability", alpha=0.7)
        ax4.bar(x + width / 2, polygenicity, width, label="Polygenicity", alpha=0.7)

        ax4.set_title("Genetic Architecture")
        ax4.set_xticks(x)
        ax4.set_xticklabels(pops)
        ax4.legend()
        ax4.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved population comparison plot to {output_file}")

    return plt.gcf()


def calculate_population_enrichment(pop_data: Dict[str, Any]) -> Dict[str, Dict[str, float]]:
    """Calculate population-specific enrichment scores.

    Args:
        pop_data: Population GWAS data

    Returns:
        Enrichment analysis results
    """
    # Simplified enrichment calculation
    enrichment_results = {}

    # Get all unique SNPs across populations
    all_snps = set()
    for pop_info in pop_data.values():
        results = pop_info.get("results")
        if hasattr(results, "index"):
            all_snps.update(results.index)

    for pop_name, pop_info in pop_data.items():
        results = pop_info.get("results")
        if not hasattr(results, "columns") or "P" not in results.columns:
            continue

        # Count significant SNPs in this population
        sig_snps = (results["P"] < 5e-8).sum()
        total_snps = len(results)

        # Expected proportion
        expected_prop = len([p for p in pop_data.keys() if p != pop_name]) / len(pop_data)

        # Enrichment calculation (simplified)
        observed_prop = sig_snps / total_snps if total_snps > 0 else 0

        if expected_prop > 0:
            enrichment = observed_prop / expected_prop
        else:
            enrichment = 1.0

        # Chi-square test (simplified p-value calculation)
        from scipy import stats

        try:
            contingency = [
                [sig_snps, total_snps - sig_snps],
                [int(expected_prop * total_snps), int((1 - expected_prop) * total_snps)],
            ]
            chi2, p_value = stats.chi2_contingency(contingency)[:2]
        except ImportError:
            p_value = 0.5  # Conservative estimate

        enrichment_results[pop_name] = {
            "enrichment": enrichment,
            "p_value": p_value,
            "observed_sig": sig_snps,
            "expected_sig": int(expected_prop * total_snps),
        }

    return enrichment_results


def analyze_genetic_architecture(pop_data: Dict[str, Any]) -> Dict[str, Dict[str, float]]:
    """Analyze genetic architecture differences across populations.

    Args:
        pop_data: Population GWAS data

    Returns:
        Genetic architecture analysis
    """
    # Simplified genetic architecture analysis
    architecture = {}

    for pop_name, pop_info in pop_data.items():
        results = pop_info.get("results")
        if not hasattr(results, "columns") or "P" in results.columns:
            continue

        # Estimate heritability (simplified)
        p_values = results["P"].dropna().values
        lambda_gc = np.median(p_values) / 0.456  # Approximate genomic control
        h2_estimate = min(lambda_gc / (lambda_gc + 1), 1.0)  # Simplified

        # Estimate polygenicity (simplified)
        sig_snps = (p_values < 5e-8).sum()
        polygenicity = sig_snps / len(p_values)

        architecture[pop_name] = {"heritability": h2_estimate, "polygenicity": polygenicity, "lambda_gc": lambda_gc}

    return architecture


def compare_analysis_methods(method_results: Dict[str, Any], output_file: Optional[str | Path] = None) -> Optional[Any]:
    """Compare different GWAS analysis methods.

    Args:
        method_results: Dictionary mapping method names to results
        output_file: Optional output file path

    Returns:
        Plot object if matplotlib available, None otherwise

    Example:
        >>> methods = {
        ...     "Linear": {"results": linear_df},
        ...     "Logistic": {"results": logistic_df},
        ...     "Mixed Model": {"results": mixed_df}
        ... }
        >>> plot = compare_analysis_methods(methods)
    """
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        logger.warning("matplotlib or seaborn not available for plotting")
        return None

    if not method_results:
        return None

    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle("GWAS Method Comparison", fontsize=16)

    # 1. P-value distribution comparison
    ax1 = axes[0, 0]
    for method_name, method_data in method_results.items():
        results = method_data.get("results")
        if hasattr(results, "columns") and "P" in results.columns:
            p_values = results["P"].dropna().values
            ax1.hist(p_values, bins=50, alpha=0.5, label=method_name, density=True)

    ax1.set_title("P-value Distribution")
    ax1.set_xlabel("P-value")
    ax1.set_ylabel("Density")
    ax1.legend()
    ax1.set_yscale("log")

    # 2. Q-Q plot comparison
    ax2 = axes[0, 1]
    for method_name, method_data in method_results.items():
        results = method_data.get("results")
        if hasattr(results, "columns") and "P" in results.columns:
            p_values = results["P"].dropna().values
            plot_qq_overlay(ax2, p_values, method_name, None)

    ax2.set_title("Q-Q Plot Comparison")
    ax2.plot([0, 1], [0, 1], "k--", alpha=0.5, label="Expected")
    ax2.legend()

    # 3. Power comparison
    ax3 = axes[1, 0]
    power_data = calculate_method_power(method_results)

    if power_data:
        methods = list(power_data.keys())
        powers = [power_data[method]["power"] for method in methods]

        bars = ax3.bar(methods, powers, alpha=0.7)
        ax3.set_title("Method Power Comparison")
        ax3.set_ylabel("Statistical Power")
        ax3.set_ylim(0, 1)

        # Add value labels
        for bar, power in zip(bars, powers):
            ax3.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01, ".3f", ha="center", va="bottom")

    # 4. Computational performance
    ax4 = axes[1, 1]
    perf_data = method_results  # Assume timing data is included

    methods = list(perf_data.keys())
    times = [perf_data[method].get("runtime_minutes", 0) for method in methods]

    if any(times):
        bars = ax4.bar(methods, times, alpha=0.7, color="orange")
        ax4.set_title("Computational Performance")
        ax4.set_ylabel("Runtime (minutes)")
        ax4.set_yscale("log")

        # Add value labels
        for bar, time_val in zip(bars, times):
            if time_val > 0:
                ax4.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01, ".1f", ha="center", va="bottom")

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved method comparison plot to {output_file}")

    return plt.gcf()


def calculate_method_power(method_results: Dict[str, Any]) -> Dict[str, Dict[str, float]]:
    """Calculate statistical power for different methods.

    Args:
        method_results: Method comparison data

    Returns:
        Power analysis results
    """
    power_results = {}

    for method_name, method_data in method_results.items():
        results = method_data.get("results")
        if not hasattr(results, "columns") or "P" in results.columns:
            continue

        # Assume some SNPs are truly associated (simplified)
        # In practice, would need ground truth data
        p_values = results["P"].dropna().values

        # Estimate power as proportion of significant findings
        # (This is a very simplified calculation)
        sig_snps = (p_values < 5e-8).sum()
        estimated_power = min(sig_snps / len(p_values), 1.0)

        power_results[method_name] = {
            "power": estimated_power,
            "significant_snps": sig_snps,
            "total_snps": len(p_values),
        }

    return power_results


def miami_plot(
    study1_results: Any,
    study2_results: Any,
    study1_name: str = "Study 1",
    study2_name: str = "Study 2",
    output_file: Optional[str | Path] = None,
    significance_threshold: float = 5e-8,
) -> Optional[Any]:
    """Create a Miami plot comparing two GWAS studies.

    Args:
        study1_results: Results from first GWAS study (DataFrame)
        study2_results: Results from second GWAS study (DataFrame)
        study1_name: Name for first study
        study2_name: Name for second study
        output_file: Optional output file path
        significance_threshold: P-value threshold for significance

    Returns:
        Plot object if matplotlib available, None otherwise

    Example:
        >>> plot = miami_plot(df1, df2, "European", "African")
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
    except ImportError:
        logger.warning("matplotlib not available for Miami plot")
        return None

    # Validate input data
    if not hasattr(study1_results, "columns") or not hasattr(study2_results, "columns"):
        logger.error("Input data must be DataFrames with appropriate columns")
        return None

    required_cols = ["CHR", "BP", "P"]
    for df, name in [(study1_results, study1_name), (study2_results, study2_name)]:
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            logger.error(f"Missing required columns in {name}: {missing_cols}")
            return None

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    # Color scheme for chromosomes
    chrom_colors = ["#1f77b4", "#ff7f0e"]  # Blue, Orange
    chroms = sorted(study1_results["CHR"].unique())

    # Plot study 1 (top panel)
    plot_manhattan_with_colors(ax1, study1_results, chrom_colors, significance_threshold, chroms)
    ax1.set_title(f"{study1_name} GWAS Results", fontsize=14)
    ax1.set_ylabel("-log₁₀(P-value)", fontsize=12)

    # Plot study 2 (bottom panel, inverted)
    plot_manhattan_with_colors(ax2, study2_results, chrom_colors, significance_threshold, chroms)
    ax2.invert_yaxis()  # Flip the y-axis
    ax2.set_title(f"{study2_name} GWAS Results", fontsize=14)
    ax2.set_ylabel("-log₁₀(P-value)", fontsize=12)
    ax2.set_xlabel("Chromosome Position", fontsize=12)

    # Add chromosome labels
    chrom_positions = []
    chrom_labels = []
    current_pos = 0

    for chrom in chroms:
        chrom_data = study1_results[study1_results["CHR"] == chrom]
        if not chrom_data.empty:
            chrom_start = current_pos
            chrom_end = current_pos + chrom_data["BP"].max()
            chrom_center = (chrom_start + chrom_end) / 2
            chrom_positions.append(chrom_center)
            chrom_labels.append(str(chrom))
            current_pos = chrom_end

    plt.xticks(chrom_positions, chrom_labels)

    # Overall title
    fig.suptitle(f"Miami Plot: {study1_name} vs {study2_name}", fontsize=16, y=0.95)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved Miami plot to {output_file}")

    return plt.gcf()


def plot_manhattan_with_colors(ax, results_df, colors, threshold, chroms):
    """Helper function to plot Manhattan plot with alternating colors."""
    current_pos = 0
    color_idx = 0

    for chrom in chroms:
        chrom_data = results_df[results_df["CHR"] == chrom].copy()
        if chrom_data.empty:
            continue

        # Adjust positions for chromosome
        chrom_data["adjusted_bp"] = chrom_data["BP"] + current_pos

        # Plot points
        ax.scatter(
            chrom_data["adjusted_bp"], -np.log10(chrom_data["P"]), c=colors[color_idx % len(colors)], s=2, alpha=0.7
        )

        # Add significance line
        max_pos = chrom_data["adjusted_bp"].max()
        ax.axhline(y=-np.log10(threshold), color="red", linestyle="--", alpha=0.7)

        current_pos = max_pos
        color_idx += 1


def multi_trait_manhattan(
    trait_results: Dict[str, Dict[str, Any]],
    output_file: Optional[str | Path] = None,
    significance_threshold: float = 5e-8,
    title: str = "Multi-Trait Manhattan Plot",
) -> Optional[Any]:
    """Create Manhattan plots for multiple traits side-by-side.

    Args:
        trait_results: Dictionary mapping trait names to GWAS results
        output_file: Optional output file path
        significance_threshold: P-value threshold for significance line
        title: Plot title

    Returns:
        Plot object if matplotlib available, None otherwise
    """
    try:
        import matplotlib.pyplot as plt
        import pandas as pd
    except ImportError:
        logger.warning("matplotlib or pandas not available for multi-trait Manhattan plot")
        return None

    if not trait_results:
        logger.error("No trait results provided")
        return None

    n_traits = len(trait_results)
    if n_traits > 6:
        logger.warning(f"Many traits ({n_traits}), plot may be crowded")

    # Create subplots
    n_cols = min(3, n_traits)
    n_rows = (n_traits + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 4 * n_rows), squeeze=False, sharey=True)
    fig.suptitle(title, fontsize=16, y=0.95)

    # Colors for chromosomes
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]

    trait_names = list(trait_results.keys())

    for i, trait_name in enumerate(trait_names):
        row = i // n_cols
        col = i % n_cols
        ax = axes[row, col]

        results = trait_results[trait_name]

        # Convert to DataFrame if needed
        if isinstance(results, dict):
            if "results" in results:
                df = results["results"]
            else:
                # Assume it's a dict with required columns
                df = pd.DataFrame(results)
        elif hasattr(results, "to_dict"):  # DataFrame-like
            df = results
        else:
            logger.error(f"Unsupported results format for trait {trait_name}")
            continue

        if df.empty:
            ax.text(0.5, 0.5, f"No data for\n{trait_name}", ha="center", va="center", transform=ax.transAxes)
            continue

        # Ensure required columns exist
        required_cols = ["CHR", "BP", "P"]
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            logger.error(f"Missing required columns {missing_cols} for trait {trait_name}")
            continue

        # Get unique chromosomes and sort
        chroms = sorted(df["CHR"].unique(), key=lambda x: int(x) if str(x).isdigit() else float("inf"))

        # Plot Manhattan for this trait
        plot_manhattan_with_colors(ax, df, colors, significance_threshold, chroms)

        ax.set_title(f"{trait_name}", fontsize=12, pad=10)
        if col == 0:  # Leftmost column
            ax.set_ylabel("-log₁₀(p-value)", fontsize=10)
        ax.set_xlabel("Chromosome", fontsize=10)

        # Set chromosome labels at chromosome centers
        current_pos = 0
        chrom_centers = []
        chrom_names = []

        for chrom in chroms:
            chrom_data = df[df["CHR"] == chrom]
            if not chrom_data.empty:
                chrom_start = current_pos
                chrom_end = current_pos + chrom_data["BP"].max()
                chrom_centers.append((chrom_start + chrom_end) / 2)
                chrom_names.append(str(chrom))
                current_pos = chrom_end

        if chrom_centers:
            ax.set_xticks(chrom_centers)
            ax.set_xticklabels(chrom_names, rotation=45, ha="right")

    # Hide unused subplots
    for i in range(n_traits, n_rows * n_cols):
        row = i // n_cols
        col = i % n_cols
        axes[row, col].set_visible(False)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved multi-trait Manhattan plot to {output_file}")

    return fig


def cross_cohort_forest(
    cohort_results: Dict[str, Dict[str, Any]],
    trait_name: str,
    output_file: Optional[str | Path] = None,
    title: str = "Cross-Cohort Forest Plot",
) -> Optional[Any]:
    """Create a forest plot comparing effect sizes across cohorts.

    Args:
        cohort_results: Dictionary mapping cohort names to GWAS results
        trait_name: Name of the trait being analyzed
        output_file: Optional output file path
        title: Plot title

    Returns:
        Plot object if matplotlib available, None otherwise
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        logger.warning("matplotlib not available for cross-cohort forest plot")
        return None

    if not cohort_results:
        logger.error("No cohort results provided")
        return None

    fig, ax = plt.subplots(figsize=(10, len(cohort_results) * 0.5 + 2))

    cohorts = list(cohort_results.keys())
    effect_sizes = []
    errors = []
    valid_cohorts = []

    # Extract effect sizes and confidence intervals
    for cohort in cohorts:
        results = cohort_results[cohort]
        if isinstance(results, dict) and "effect_size" in results:
            effect_sizes.append(results["effect_size"])
            # Calculate error bars (assume 95% CI if available)
            if "ci_lower" in results and "ci_upper" in results:
                error_lower = results["effect_size"] - results["ci_lower"]
                error_upper = results["ci_upper"] - results["effect_size"]
                errors.append([error_lower, error_upper])
            else:
                # Default error bars
                errors.append([0.1, 0.1])
            valid_cohorts.append(cohort)
        else:
            logger.warning(f"No effect size found for cohort {cohort}")

    if not valid_cohorts:
        logger.error("No valid effect sizes found")
        return None

    # Plot forest plot
    y_positions = range(len(valid_cohorts))

    # Plot error bars
    for i, (cohort, effect_size, error) in enumerate(zip(valid_cohorts, effect_sizes, errors)):
        ax.errorbar(effect_size, i, xerr=np.array(error).reshape(2, 1), fmt="o", color="blue", capsize=3, markersize=6)

    # Add vertical line at zero effect
    ax.axvline(x=0, color="red", linestyle="--", alpha=0.7, label="No effect")

    # Add vertical line at overall effect (simple average)
    overall_effect = np.mean(effect_sizes)
    ax.axvline(x=overall_effect, color="green", linestyle="-", alpha=0.7, label=f"Overall effect: {overall_effect:.3f}")

    ax.set_yticks(y_positions)
    ax.set_yticklabels(valid_cohorts)
    ax.set_xlabel("Effect Size", fontsize=12)
    ax.set_title(f"{title} - {trait_name}", fontsize=14, pad=20)
    ax.grid(True, alpha=0.3, axis="x")
    ax.legend()

    # Add effect size values as text
    for i, (effect_size, error) in enumerate(zip(effect_sizes, errors)):
        ax.text(effect_size + max(error) + 0.01, i, f"{effect_size:.3f}", va="center", fontsize=9)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved cross-cohort forest plot to {output_file}")

    return fig


def concordance_plot(
    study1_results: list[dict],
    study2_results: list[dict],
    output_path: Optional[str | Path] = None,
    figsize: tuple[int, int] = (8, 8),
) -> dict[str, Any]:
    """Create a concordance plot comparing two GWAS studies.

    Args:
        study1_results: Results from first study
        study2_results: Results from second study
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
        return {"status": "error", "message": "matplotlib not available"}

    try:
        # Extract p-values from both studies
        p_values1 = []
        p_values2 = []

        # Create SNP ID to p-value mapping for both studies
        snp_to_p1 = {}
        snp_to_p2 = {}

        for result in study1_results:
            snp_id = result.get("ID", f"{result.get('CHROM')}:{result.get('POS')}")
            p_val = result.get("p_value", result.get("P", 1.0))
            snp_to_p1[snp_id] = p_val

        for result in study2_results:
            snp_id = result.get("ID", f"{result.get('CHROM')}:{result.get('POS')}")
            p_val = result.get("p_value", result.get("P", 1.0))
            snp_to_p2[snp_id] = p_val

        # Find overlapping SNPs
        common_snps = set(snp_to_p1.keys()) & set(snp_to_p2.keys())

        if not common_snps:
            return {"status": "error", "message": "No overlapping SNPs found"}

        # Get p-values for common SNPs
        for snp in common_snps:
            p_values1.append(-math.log10(max(snp_to_p1[snp], 1e-300)))
            p_values2.append(-math.log10(max(snp_to_p2[snp], 1e-300)))

        # Create plot
        fig, ax = plt.subplots(figsize=figsize)

        ax.scatter(p_values1, p_values2, alpha=0.6, s=20, color="blue")

        # Add diagonal line
        max_val = max(max(p_values1), max(p_values2))
        ax.plot([0, max_val], [0, max_val], "r--", alpha=0.7, label="Perfect concordance")

        ax.set_xlabel("Study 1 -log₁₀(p)")
        ax.set_ylabel("Study 2 -log₁₀(p)")
        ax.set_title(f"GWAS Concordance Plot\n{n:,} overlapping SNPs")
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches="tight")

        return {
            "status": "success",
            "n_overlapping_snps": len(common_snps),
            "study1_snps": len(snp_to_p1),
            "study2_snps": len(snp_to_p2),
            "output_path": str(output_path) if output_path else None,
        }

    except Exception as e:
        return {"status": "error", "message": str(e)}
