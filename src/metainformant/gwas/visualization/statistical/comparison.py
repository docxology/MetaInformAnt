"""GWAS comparison visualization utilities.

This module provides tools for comparing GWAS results across different
studies, populations, or analysis methods.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

import numpy as np

from metainformant.core.utils import logging
from metainformant.gwas.analysis.correction import lambda_gc_from_p_values

logger = logging.get_logger(__name__)


def _chi_square_p_value_2x2(table: Sequence[Sequence[float]]) -> float:
    """P-value of the chi-square test of independence for a 2x2 table, with Yates' continuity correction.

    Equal to ``scipy.stats.chi2_contingency(table)[1]`` (1 degree of freedom, ``correction=True``) but needs
    no SciPy: with 1 degree of freedom the chi-square survival function is ``erfc(sqrt(x / 2))``.

    Returns:
        The p-value, or NaN when the test is undefined (a zero row or column total leaves a zero expected
        count). It is never replaced by a made-up value.
    """
    (a, b), (c, d) = table
    row_totals = (a + b, c + d)
    col_totals = (a + c, b + d)
    grand_total = a + b + c + d
    if 0 in row_totals or 0 in col_totals:
        return float("nan")

    chi2 = 0.0
    for i, observed_row in enumerate(((a, b), (c, d))):
        for j, observed in enumerate(observed_row):
            expected = row_totals[i] * col_totals[j] / grand_total
            difference = expected - observed
            # Yates: move each observed count 0.5 towards its expected count (never past it)
            corrected = observed + math.copysign(min(0.5, abs(difference)), difference) if difference else observed
            chi2 += (corrected - expected) ** 2 / expected
    return math.erfc(math.sqrt(chi2 / 2.0))


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
        import pandas as pd
        import seaborn as sns
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
        if results is None:
            continue
        if isinstance(results, dict):
            # Convert dict to DataFrame if needed
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
        if results is None:
            continue
        if isinstance(results, dict):
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
        if results is None:
            continue
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
        if results is None:
            continue
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
            cellText=summary_data,
            colLabels=["Study", "Total SNPs", "Significant SNPs"],
            loc="center",
            cellLoc="center",
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
        import pandas as pd  # noqa: F401
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
        ax.text(
            chrom_center,
            ax.get_ylim()[0] - 0.5,
            str(chrom),
            ha="center",
            va="top",
            fontsize=8,
        )

    ax.legend()


def plot_qq_overlay(ax: Any, p_values: np.ndarray, label: str, color: str | None) -> None:
    """Plot Q-Q plot for one study on given axis."""
    # Remove NA values and sort
    p_values = p_values[~np.isnan(p_values)]
    p_values = np.sort(p_values)

    # Expected p-values under null
    n = len(p_values)
    expected = np.arange(1, n + 1) / (n + 1)

    # Convert to -log10
    observed_log = -np.log10(np.clip(p_values, 1e-300, None))
    expected_log = -np.log10(expected)

    # Plot
    ax.scatter(expected_log, observed_log, alpha=0.6, s=2, color=color, label=label)


def compare_populations(
    pop_data: Dict[str, Any],
    trait_name: str = "Trait",
    output_file: Optional[str | Path] = None,
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
        ax1.boxplot(effect_sizes, tick_labels=pop_names)
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

        sns.heatmap(
            overlap_matrix,
            annot=True,
            fmt=".0f",
            cmap="Blues",
            xticklabels=pops,
            yticklabels=pops,
            ax=ax2,
        )
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

        # Add p-value annotations (a test that is undefined for a population is marked, not given a value)
        for bar, p_val in zip(bars, p_values):
            height = bar.get_height()
            ax3.text(
                bar.get_x() + bar.get_width() / 2,
                height,
                f"p={p_val:.2e}" if math.isfinite(p_val) else "p=n/a",
                ha="center",
                va="bottom",
                fontsize=8,
            )

    # 4. Genome-wide significant fraction per population (the one architecture summary the p-values support)
    ax4 = axes[1, 1]
    arch_data = analyze_genetic_architecture(pop_data)

    if arch_data:
        pops = list(arch_data.keys())
        significant_fraction = [arch_data[pop]["polygenicity"] for pop in pops]

        ax4.bar(pops, significant_fraction, alpha=0.7)
        ax4.set_title("Genome-wide Significant Fraction")
        ax4.set_ylabel("Fraction of SNPs genome-wide significant")
        ax4.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved population comparison plot to {output_file}")

    return plt.gcf()


def calculate_population_enrichment(
    pop_data: Dict[str, Any],
) -> Dict[str, Dict[str, float]]:
    """Calculate population-specific enrichment scores.

    For each population whose ``results`` has a ``P`` column, ``observed_sig`` is its number of genome-wide
    significant SNPs (P < 5e-8) among ``len(results)`` SNPs. ``enrichment`` is the significant fraction divided
    by ``expected_prop``, the share of the *other* populations (0.5 for two populations); ``expected_sig`` is
    ``int(expected_prop * len(results))``. ``p_value`` is the chi-square test (1 degree of freedom, Yates'
    continuity correction; no SciPy needed) of the observed significant / non-significant counts against the
    counts expected at that share. The expected share is a modelling assumption of this summary, not data.

    ``p_value`` is NaN where the test is undefined (for example a single population, whose expected share is
    0, with no significant SNP). It is never replaced by a made-up value.

    Args:
        pop_data: Population GWAS data

    Returns:
        Enrichment analysis results
    """
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

        # Enrichment: significant fraction relative to the expected share
        observed_prop = sig_snps / total_snps if total_snps > 0 else 0

        if expected_prop > 0:
            enrichment = observed_prop / expected_prop
        else:
            enrichment = 1.0

        # Chi-square test of the observed counts against the counts expected at that share
        p_value = _chi_square_p_value_2x2(
            [
                [sig_snps, total_snps - sig_snps],
                [
                    int(expected_prop * total_snps),
                    int((1 - expected_prop) * total_snps),
                ],
            ]
        )

        enrichment_results[pop_name] = {
            "enrichment": enrichment,
            "p_value": p_value,
            "observed_sig": sig_snps,
            "expected_sig": int(expected_prop * total_snps),
        }

    return enrichment_results


def analyze_genetic_architecture(
    pop_data: Dict[str, Any],
) -> Dict[str, Dict[str, float]]:
    """Summarise the association signal of each population from its p-values.

    For every population whose ``results`` has a ``P`` column the result holds:

    * ``lambda_gc``: the genomic-control inflation factor of its p-values (1.0 when it cannot be computed);
    * ``polygenicity``: the fraction of tested SNPs reaching genome-wide significance (P < 5e-8). Despite the
      key name this measures *detected* signal, which grows with sample size, not the number of causal variants.

    There is deliberately no heritability. SNP heritability cannot be estimated from p-values alone (it needs
    effect sizes and sample sizes, or LD scores). An earlier version reported lambda / (lambda + 1) under that
    name; that is a fixed function of the inflation factor and carries no genetic information.

    Args:
        pop_data: Population GWAS data

    Returns:
        Per-population ``{"polygenicity": ..., "lambda_gc": ...}``
    """
    architecture = {}

    for pop_name, pop_info in pop_data.items():
        results = pop_info.get("results")
        if not hasattr(results, "columns") or "P" not in results.columns:
            continue

        p_values = results["P"].dropna().values
        lambda_gc = lambda_gc_from_p_values(p_values.tolist()) or 1.0

        # Fraction of tested SNPs that reach genome-wide significance
        sig_snps = (p_values < 5e-8).sum()
        polygenicity = sig_snps / len(p_values)

        architecture[pop_name] = {
            "polygenicity": polygenicity,
            "lambda_gc": lambda_gc,
        }

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
        import seaborn as sns  # noqa: F401
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
        # Only a comparison against known causal SNPs is power; otherwise the bars say what they are
        is_power = all(power_data[method]["basis"] == "causal_snps" for method in methods)

        bars = ax3.bar(methods, powers, alpha=0.7)
        ax3.set_title("Method Power Comparison" if is_power else "Genome-wide Significant Fraction")
        ax3.set_ylabel("Statistical Power" if is_power else "Fraction of SNPs significant")
        ax3.set_ylim(0, 1)

        # Add value labels
        for bar, power in zip(bars, powers):
            ax3.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.01,
                f"{power:.3f}",
                ha="center",
                va="bottom",
            )

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
                ax4.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.01,
                    f"{time_val:.1f}",
                    ha="center",
                    va="bottom",
                )

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved method comparison plot to {output_file}")

    return plt.gcf()


def calculate_method_power(
    method_results: Dict[str, Any],
    significance_threshold: float = 5e-8,
) -> Dict[str, Dict[str, Any]]:
    """Summarise what each method calls significant, and the power that implies when truth is known.

    Power is the probability of calling a SNP significant *given that it is truly associated*, so it can only
    be measured against known truth. Each method therefore gets one of two honestly-labelled values under
    the ``"power"`` key, told apart by ``"basis"``:

    * ``"causal_snps"`` -- the method's entry carries ``"causal_snps"``, the identifiers of the truly
      associated SNPs (matched against ``results.index``, as in a simulation study). ``power`` is the
      fraction of those that were tested and called significant; ``"causal_snps"`` and
      ``"causal_detected"`` in the output give the counts. Causal SNPs that were not tested are not counted.
    * ``"significant_fraction"`` -- no usable ground truth. ``power`` is only the fraction of *all* tested
      SNPs reaching the threshold, a discovery rate that stays small whenever most SNPs are null. It is
      reported under ``"power"`` so existing callers keep working, but it is not statistical power.

    Methods whose ``"results"`` is not a table with a ``P`` column, or that have no valid p-value, are omitted.

    Args:
        method_results: Mapping of method name to ``{"results": DataFrame with a P column,
            "causal_snps": optional collection of truly associated SNP identifiers, ...}``
        significance_threshold: P-value below which a SNP counts as called (default genome-wide 5e-8)

    Returns:
        Per method: ``power``, ``basis``, ``significant_snps`` and ``total_snps`` (counts over all tested
        SNPs), plus ``causal_snps`` / ``causal_detected`` when ground truth was used.
    """
    power_results: Dict[str, Dict[str, Any]] = {}

    for method_name, method_data in method_results.items():
        results = method_data.get("results")
        if not hasattr(results, "columns") or "P" not in results.columns:
            continue

        p_values = results["P"].dropna()
        total_snps = len(p_values)
        if total_snps == 0:
            continue

        sig_snps = int((p_values < significance_threshold).sum())
        entry: Dict[str, Any] = {
            "power": sig_snps / total_snps,
            "basis": "significant_fraction",
            "significant_snps": sig_snps,
            "total_snps": total_snps,
        }

        causal = method_data.get("causal_snps")
        if causal is not None:
            causal_ids = [causal] if isinstance(causal, str) else list(causal)
            tested_causal = results.index.isin(causal_ids)
            n_causal = int(tested_causal.sum())
            if n_causal:
                detected = int((results.loc[tested_causal, "P"] < significance_threshold).sum())
                entry.update(
                    {
                        "power": detected / n_causal,
                        "basis": "causal_snps",
                        "causal_snps": n_causal,
                        "causal_detected": detected,
                    }
                )
            else:
                logger.warning(
                    f"None of the {len(causal_ids)} causal SNPs given for method {method_name!r} were tested; "
                    "reporting the significant fraction instead of power"
                )

        power_results[method_name] = entry

    return power_results


def miami_plot(
    trait1_results: Any,
    trait2_results: Any,
    output_file: Optional[str | Path] = None,
    trait1_name: str = "Trait 1",
    trait2_name: str = "Trait 2",
    significance_threshold: float = 5e-8,
    **kwargs: Any,
) -> Dict[str, Any]:
    """Create a Miami plot comparing two GWAS studies.

    Args:
        trait1_results: Results from first GWAS study (DataFrame or list of dicts)
        trait2_results: Results from second GWAS study (DataFrame or list of dicts)
        output_file: Optional output file path
        trait1_name: Name for first trait
        trait2_name: Name for second trait
        significance_threshold: P-value threshold for significance
        **kwargs: Additional keyword arguments

    Returns:
        Dictionary with status and metadata
    """
    try:
        import matplotlib.pyplot as plt
        import pandas as pd
    except ImportError:
        return {"status": "failed", "error": "matplotlib or pandas not available"}

    try:
        # Convert list-of-dicts to DataFrame if needed
        if isinstance(trait1_results, list):
            trait1_df = pd.DataFrame(trait1_results)
        else:
            trait1_df = trait1_results

        if isinstance(trait2_results, list):
            trait2_df = pd.DataFrame(trait2_results)
        else:
            trait2_df = trait2_results

        # Validate non-empty
        if trait1_df.empty or trait2_df.empty:
            return {
                "status": "failed",
                "error": "One or both trait result sets are empty",
            }

        # Normalize column names: accept CHROM->CHR, POS->BP, p_value->P
        for df in [trait1_df, trait2_df]:
            if "CHROM" in df.columns and "CHR" not in df.columns:
                df.rename(columns={"CHROM": "CHR"}, inplace=True)
            if "POS" in df.columns and "BP" not in df.columns:
                df.rename(columns={"POS": "BP"}, inplace=True)
            if "p_value" in df.columns and "P" not in df.columns:
                df.rename(columns={"p_value": "P"}, inplace=True)

        required_cols = ["CHR", "BP", "P"]
        for df, name in [(trait1_df, trait1_name), (trait2_df, trait2_name)]:
            missing_cols = [col for col in required_cols if col not in df.columns]
            if missing_cols:
                return {
                    "status": "failed",
                    "error": f"Missing required columns in {name}: {missing_cols}",
                }

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

        # Color scheme for chromosomes
        chrom_colors = ["#1f77b4", "#ff7f0e"]
        chroms = sorted(trait1_df["CHR"].unique())

        # Plot trait 1 (top panel)
        plot_manhattan_with_colors(ax1, trait1_df, chrom_colors, significance_threshold, chroms)
        ax1.set_title(f"{trait1_name} GWAS Results", fontsize=14)
        ax1.set_ylabel("-log10(P-value)", fontsize=12)

        # Plot trait 2 (bottom panel, inverted)
        plot_manhattan_with_colors(ax2, trait2_df, chrom_colors, significance_threshold, chroms)
        ax2.invert_yaxis()
        ax2.set_title(f"{trait2_name} GWAS Results", fontsize=14)
        ax2.set_ylabel("-log10(P-value)", fontsize=12)
        ax2.set_xlabel("Chromosome Position", fontsize=12)

        fig.suptitle(f"Miami Plot: {trait1_name} vs {trait2_name}", fontsize=16, y=0.95)
        plt.tight_layout()

        if output_file:
            plt.savefig(str(output_file), dpi=300, bbox_inches="tight")
            logger.info(f"Saved Miami plot to {output_file}")

        plt.close(fig)

        return {
            "status": "success",
            "trait1_variants": len(trait1_df),
            "trait2_variants": len(trait2_df),
            "output_file": str(output_file) if output_file else None,
        }

    except Exception as e:
        return {"status": "failed", "error": str(e)}


def plot_manhattan_with_colors(
    ax: Any,
    results_df: Any,
    colors: Sequence[Any],
    threshold: float,
    chroms: Sequence[Any],
) -> None:
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
            chrom_data["adjusted_bp"],
            -np.log10(chrom_data["P"]),
            c=colors[color_idx % len(colors)],
            s=2,
            alpha=0.7,
        )

        # Add significance line
        max_pos = chrom_data["adjusted_bp"].max()
        ax.axhline(y=-np.log10(threshold), color="red", linestyle="--", alpha=0.7)

        current_pos = max_pos
        color_idx += 1


def multi_trait_manhattan(
    trait_results: Dict[str, Any],
    output_file: Optional[str | Path] = None,
    significance_threshold: float = 5e-8,
    title: str = "Multi-Trait Manhattan Plot",
    **kwargs: Any,
) -> Dict[str, Any]:
    """Create Manhattan plots for multiple traits side-by-side.

    Args:
        trait_results: Dictionary mapping trait names to GWAS results (DataFrame or list of dicts)
        output_file: Optional output file path
        significance_threshold: P-value threshold for significance line
        title: Plot title
        **kwargs: Additional keyword arguments

    Returns:
        Dictionary with status and metadata
    """
    try:
        import matplotlib.pyplot as plt
        import pandas as pd
    except ImportError:
        return {"status": "failed", "error": "matplotlib or pandas not available"}

    try:
        if not trait_results:
            return {"status": "failed", "error": "No trait results provided"}

        n_traits = len(trait_results)
        if n_traits < 2:
            return {
                "status": "failed",
                "error": f"Need at least 2 traits for multi-trait plot, got {n_traits}",
            }

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
            if isinstance(results, list):
                df = pd.DataFrame(results)
            elif isinstance(results, dict):
                if "results" in results:
                    df = results["results"]
                else:
                    df = pd.DataFrame(results)
            elif hasattr(results, "to_dict"):  # DataFrame-like
                df = results
            else:
                logger.error(f"Unsupported results format for trait {trait_name}")
                continue

            if df.empty:
                ax.text(
                    0.5,
                    0.5,
                    f"No data for\n{trait_name}",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
                continue

            # Normalize column names: accept CHROM->CHR, POS->BP, p_value->P
            if "CHROM" in df.columns and "CHR" not in df.columns:
                df = df.rename(columns={"CHROM": "CHR"})
            if "POS" in df.columns and "BP" not in df.columns:
                df = df.rename(columns={"POS": "BP"})
            if "p_value" in df.columns and "P" not in df.columns:
                df = df.rename(columns={"p_value": "P"})

            # Ensure required columns exist
            required_cols = ["CHR", "BP", "P"]
            missing = [c for c in required_cols if c not in df.columns]
            if missing:
                logger.error(f"Missing required columns {missing} for trait {trait_name}")
                continue

            # Get unique chromosomes and sort
            chroms = sorted(
                df["CHR"].unique(),
                key=lambda x: int(x) if str(x).isdigit() else float("inf"),
            )

            # Plot Manhattan for this trait
            plot_manhattan_with_colors(ax, df, colors, significance_threshold, chroms)

            ax.set_title(f"{trait_name}", fontsize=12, pad=10)
            if col == 0:
                ax.set_ylabel("-log10(p-value)", fontsize=10)
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
            plt.savefig(str(output_file), dpi=300, bbox_inches="tight")
            logger.info(f"Saved multi-trait Manhattan plot to {output_file}")

        plt.close(fig)

        return {
            "status": "success",
            "num_traits": n_traits,
            "output_file": str(output_file) if output_file else None,
        }

    except Exception as e:
        return {"status": "failed", "error": str(e)}


def cross_cohort_forest(
    cohort_results: Dict[str, Any],
    output_file: Optional[str | Path] = None,
    variant_id: Optional[str] = None,
    trait_name: str = "",
    title: str = "Cross-Cohort Forest Plot",
    **kwargs: Any,
) -> Dict[str, Any]:
    """Create a forest plot comparing effect sizes across cohorts.

    Args:
        cohort_results: Dictionary mapping cohort names to GWAS results
        output_file: Optional output file path
        variant_id: Optional variant identifier for the forest plot
        trait_name: Name of the trait being analyzed
        title: Plot title
        **kwargs: Additional keyword arguments

    Returns:
        Dictionary with status and metadata. Returns skipped status when
        cohort data contains raw variant lists rather than pre-computed
        effect sizes, as full meta-analysis is not yet implemented.
    """
    if not cohort_results:
        return {"status": "failed", "error": "No cohort results provided"}

    # Check if cohort data contains pre-computed effect sizes or raw variant lists
    # If raw variant lists, return skipped since meta-analysis is not yet implemented
    has_precomputed = False
    for cohort_name, results in cohort_results.items():
        if isinstance(results, dict) and "effect_size" in results:
            has_precomputed = True
            break
        if isinstance(results, list):
            # Raw variant list -- meta-analysis not implemented
            return {
                "status": "skipped",
                "message": "Meta-analysis forest plot not yet implemented for raw variant data",
                "variant_id": variant_id,
                "num_cohorts": len(cohort_results),
            }

    if not has_precomputed:
        return {
            "status": "skipped",
            "message": "No pre-computed effect sizes found in cohort data",
            "variant_id": variant_id,
            "num_cohorts": len(cohort_results),
        }

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return {"status": "failed", "error": "matplotlib not available"}

    try:
        fig, ax = plt.subplots(figsize=(10, len(cohort_results) * 0.5 + 2))

        cohorts = list(cohort_results.keys())
        effect_sizes = []
        errors = []
        valid_cohorts = []

        for cohort in cohorts:
            results = cohort_results[cohort]
            if isinstance(results, dict) and "effect_size" in results:
                effect_sizes.append(results["effect_size"])
                if "ci_lower" in results and "ci_upper" in results:
                    error_lower = results["effect_size"] - results["ci_lower"]
                    error_upper = results["ci_upper"] - results["effect_size"]
                    errors.append([error_lower, error_upper])
                else:
                    errors.append([0.1, 0.1])
                valid_cohorts.append(cohort)

        if not valid_cohorts:
            plt.close(fig)
            return {"status": "failed", "error": "No valid effect sizes found"}

        y_positions = range(len(valid_cohorts))
        for i, (cohort, effect_size, error) in enumerate(zip(valid_cohorts, effect_sizes, errors)):
            ax.errorbar(
                effect_size,
                i,
                xerr=np.array(error).reshape(2, 1),
                fmt="o",
                color="blue",
                capsize=3,
                markersize=6,
            )

        ax.axvline(x=0, color="red", linestyle="--", alpha=0.7, label="No effect")
        overall_effect = np.mean(effect_sizes)
        ax.axvline(
            x=overall_effect,
            color="green",
            linestyle="-",
            alpha=0.7,
            label=f"Overall effect: {overall_effect:.3f}",
        )

        ax.set_yticks(list(y_positions))
        ax.set_yticklabels(valid_cohorts)
        ax.set_xlabel("Effect Size", fontsize=12)
        plot_title = f"{title} - {trait_name}" if trait_name else title
        if variant_id:
            plot_title += f" ({variant_id})"
        ax.set_title(plot_title, fontsize=14, pad=20)
        ax.grid(True, alpha=0.3, axis="x")
        ax.legend()

        for i, (effect_size, error) in enumerate(zip(effect_sizes, errors)):
            ax.text(
                effect_size + max(error) + 0.01,
                i,
                f"{effect_size:.3f}",
                va="center",
                fontsize=9,
            )

        plt.tight_layout()

        if output_file:
            plt.savefig(str(output_file), dpi=300, bbox_inches="tight")
            logger.info(f"Saved cross-cohort forest plot to {output_file}")

        plt.close(fig)

        return {
            "status": "success",
            "num_cohorts": len(valid_cohorts),
            "overall_effect": float(overall_effect),
            "variant_id": variant_id,
            "output_file": str(output_file) if output_file else None,
        }

    except Exception as e:
        return {"status": "failed", "error": str(e)}


def concordance_plot(
    study1_results: list[dict],
    study2_results: list[dict],
    output_path: Optional[str | Path] = None,
    figsize: tuple[int, int] = (8, 8),
    title: Optional[str] = None,
    **kwargs: Any,
) -> dict[str, Any]:
    """Create a concordance plot comparing two GWAS studies.

    Args:
        study1_results: Results from first study
        study2_results: Results from second study
        output_path: Path to save the plot
        figsize: Figure size
        title: Optional plot title
        **kwargs: Additional keyword arguments

    Returns:
        Dictionary with plot status and metadata
    """
    import math

    try:
        import matplotlib.pyplot as plt

        HAS_MATPLOTLIB = True
    except ImportError:
        HAS_MATPLOTLIB = False

    if not HAS_MATPLOTLIB:
        return {"status": "failed", "error": "matplotlib not available"}

    try:
        # Create SNP ID to data mapping for both studies
        snp_to_data1: Dict[str, Dict[str, Any]] = {}
        snp_to_data2: Dict[str, Dict[str, Any]] = {}

        for result in study1_results:
            snp_id = result.get("ID", f"{result.get('CHROM')}:{result.get('POS')}")
            snp_to_data1[snp_id] = result

        for result in study2_results:
            snp_id = result.get("ID", f"{result.get('CHROM')}:{result.get('POS')}")
            snp_to_data2[snp_id] = result

        # Find overlapping SNPs
        common_snps = set(snp_to_data1.keys()) & set(snp_to_data2.keys())

        if not common_snps:
            return {"status": "failed", "error": "No overlapping SNPs found"}

        # Extract p-values and betas for common SNPs
        p_values1 = []
        p_values2 = []
        betas1 = []
        betas2 = []

        for snp in sorted(common_snps):
            d1 = snp_to_data1[snp]
            d2 = snp_to_data2[snp]

            p1 = d1.get("p_value", d1.get("P", 1.0))
            p2 = d2.get("p_value", d2.get("P", 1.0))
            p_values1.append(-math.log10(max(p1, 1e-300)))
            p_values2.append(-math.log10(max(p2, 1e-300)))

            b1 = d1.get("beta", d1.get("BETA"))
            b2 = d2.get("beta", d2.get("BETA"))
            if b1 is not None and b2 is not None:
                betas1.append(b1)
                betas2.append(b2)

        n = len(common_snps)

        # Calculate correlation (use betas if available, else p-values)
        correlation = None
        if len(betas1) >= 2:
            arr1 = np.array(betas1, dtype=float)
            arr2 = np.array(betas2, dtype=float)
            if np.std(arr1) > 0 and np.std(arr2) > 0:
                correlation = float(np.corrcoef(arr1, arr2)[0, 1])
        if correlation is None and len(p_values1) >= 2:
            arr1 = np.array(p_values1, dtype=float)
            arr2 = np.array(p_values2, dtype=float)
            if np.std(arr1) > 0 and np.std(arr2) > 0:
                correlation = float(np.corrcoef(arr1, arr2)[0, 1])

        # Create plot
        fig, ax = plt.subplots(figsize=figsize)

        ax.scatter(p_values1, p_values2, alpha=0.6, s=20, color="blue")

        # Add diagonal line
        max_val = max(max(p_values1), max(p_values2))
        ax.plot([0, max_val], [0, max_val], "r--", alpha=0.7, label="Perfect concordance")

        ax.set_xlabel("Study 1 -log10(p)")
        ax.set_ylabel("Study 2 -log10(p)")
        plot_title = title if title else f"GWAS Concordance Plot\n{n:,} overlapping SNPs"
        ax.set_title(plot_title)
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        if output_path:
            plt.savefig(str(output_path), dpi=300, bbox_inches="tight")

        plt.close(fig)

        concordance_result: Dict[str, Any] = {
            "status": "success",
            "n_overlapping_snps": n,
            "study1_snps": len(snp_to_data1),
            "study2_snps": len(snp_to_data2),
            "output_path": str(output_path) if output_path else None,
        }
        if correlation is not None:
            concordance_result["correlation"] = correlation

        return concordance_result

    except Exception as e:
        return {"status": "failed", "error": str(e)}
