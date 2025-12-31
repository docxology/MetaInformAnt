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


def compare_gwas_studies(study_results: Dict[str, Dict[str, Any]],
                        output_file: Optional[str | Path] = None) -> Optional[Any]:
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
    fig.suptitle('GWAS Study Comparison', fontsize=16)

    # Prepare data
    study_names = list(study_results.keys())
    colors = sns.color_palette("husl", len(study_names))

    # 1. Manhattan plot comparison
    ax1 = axes[0, 0]
    for i, (study_name, study_data) in enumerate(study_results.items()):
        results = study_data.get('results')
        if isinstance(results, dict):
            # Convert dict to DataFrame if needed
            import pandas as pd
            results = pd.DataFrame(results)

        if hasattr(results, 'columns') and 'CHR' in results.columns and 'P' in results.columns:
            plot_manhattan_overlay(ax1, results, study_name, colors[i])

    ax1.set_title('Manhattan Plot Comparison')
    ax1.set_xlabel('Chromosome')
    ax1.set_ylabel('-log10(P-value)')

    # 2. QQ plot comparison
    ax2 = axes[0, 1]
    for i, (study_name, study_data) in enumerate(study_results.items()):
        results = study_data.get('results')
        if isinstance(results, dict):
            import pandas as pd
            results = pd.DataFrame(results)

        if hasattr(results, 'columns') and 'P' in results.columns:
            p_values = results['P'].dropna().values
            plot_qq_overlay(ax2, p_values, study_name, colors[i])

    ax2.set_title('Q-Q Plot Comparison')
    ax2.plot([0, 1], [0, 1], 'k--', alpha=0.5, label='Expected')

    # 3. Effect size distribution
    ax3 = axes[1, 0]
    effect_sizes = []
    labels = []

    for study_name, study_data in study_results.items():
        results = study_data.get('results')
        if isinstance(results, dict):
            import pandas as pd
            results = pd.DataFrame(results)

        if hasattr(results, 'columns') and 'BETA' in results.columns:
            beta_values = results['BETA'].dropna().values
            effect_sizes.append(beta_values)
            labels.append(study_name)

    if effect_sizes:
        ax3.hist(effect_sizes, bins=50, alpha=0.7, label=labels)
        ax3.set_title('Effect Size Distribution')
        ax3.set_xlabel('Effect Size (Beta)')
        ax3.set_ylabel('Frequency')
        ax3.legend()

    # 4. Summary statistics table
    ax4 = axes[1, 1]
    ax4.axis('off')

    # Create summary table
    summary_data = []
    for study_name, study_data in study_results.items():
        results = study_data.get('results')
        if isinstance(results, dict):
            import pandas as pd
            results = pd.DataFrame(results)

        n_snps = len(results) if hasattr(results, '__len__') else 0
        sig_threshold = study_data.get('significance_threshold', 5e-8)

        if hasattr(results, 'columns') and 'P' in results.columns:
            n_sig = (results['P'] < sig_threshold).sum()
        else:
            n_sig = 0

        summary_data.append([study_name, n_snps, n_sig])

    if summary_data:
        table = ax4.table(cellText=summary_data,
                         colLabels=['Study', 'Total SNPs', 'Significant SNPs'],
                         loc='center', cellLoc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 1.5)

    ax4.set_title('Study Summary', pad=20)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved GWAS comparison plot to {output_file}")

    return plt.gcf()


def plot_manhattan_overlay(ax: Any, results: Any, label: str, color: str) -> None:
    """Plot Manhattan plot for one study on given axis."""
    try:
        import pandas as pd
    except ImportError:
        return

    if not hasattr(results, 'columns') or 'CHR' not in results.columns or 'P' not in results.columns:
        return

    # Prepare data
    df = results.copy()
    df['-logP'] = -np.log10(df['P'].clip(lower=1e-300))  # Avoid log(0)

    # Create chromosome positions
    chromosomes = sorted(df['CHR'].unique())
    chrom_starts = {}
    current_pos = 0

    for chrom in chromosomes:
        chrom_starts[chrom] = current_pos
        chrom_size = len(df[df['CHR'] == chrom])
        current_pos += chrom_size + 1000  # Gap between chromosomes

    df['pos'] = df.apply(lambda row: chrom_starts[row['CHR']] + row.name, axis=1)

    # Plot
    ax.scatter(df['pos'], df['-logP'], alpha=0.6, s=2, color=color, label=label)

    # Add chromosome labels
    for chrom in chromosomes:
        chrom_center = chrom_starts[chrom] + len(df[df['CHR'] == chrom]) / 2
        ax.text(chrom_center, ax.get_ylim()[0] - 0.5, str(chrom),
               ha='center', va='top', fontsize=8)

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


def compare_populations(pop_data: Dict[str, Any], trait_name: str = "Trait",
                       output_file: Optional[str | Path] = None) -> Optional[Any]:
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
    fig.suptitle(f'{trait_name} - Population Comparison', fontsize=16)

    # 1. Effect size comparison
    ax1 = axes[0, 0]
    pop_names = []
    effect_sizes = []

    for pop_name, pop_info in pop_data.items():
        results = pop_info.get('results')
        if hasattr(results, 'columns') and 'BETA' in results.columns:
            beta_values = results['BETA'].dropna().values
            effect_sizes.append(beta_values)
            pop_names.append(pop_name)

    if effect_sizes:
        ax1.boxplot(effect_sizes, labels=pop_names)
        ax1.set_title('Effect Size Distribution by Population')
        ax1.set_ylabel('Effect Size (Beta)')
        ax1.grid(True, alpha=0.3)

    # 2. Significant SNP overlap
    ax2 = axes[0, 1]
    sig_snps = {}

    threshold = 5e-8  # Common GWAS threshold

    for pop_name, pop_info in pop_data.items():
        results = pop_info.get('results')
        if hasattr(results, 'columns') and 'P' in results.columns:
            sig_mask = results['P'] < threshold
            sig_snps[pop_name] = set(results[sig_mask].index if hasattr(results, 'index') else [])

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

        sns.heatmap(overlap_matrix, annot=True, fmt='.0f', cmap='Blues',
                   xticklabels=pops, yticklabels=pops, ax=ax2)
        ax2.set_title('Significant SNP Overlap')

    # 3. Population-specific enrichment
    ax3 = axes[1, 0]
    enrichment_data = calculate_population_enrichment(pop_data)

    if enrichment_data:
        pops = list(enrichment_data.keys())
        enrichments = [enrichment_data[pop]['enrichment'] for pop in pops]
        p_values = [enrichment_data[pop]['p_value'] for pop in pops]

        bars = ax3.bar(pops, enrichments, alpha=0.7)
        ax3.set_title('Population-Specific Enrichment')
        ax3.set_ylabel('Enrichment Score')
        ax3.set_yscale('log')

        # Add p-value annotations
        for bar, p_val in zip(bars, p_values):
            height = bar.get_height()
            ax3.text(bar.get_x() + bar.get_width()/2, height,
                    f'p={p_val:.2e}', ha='center', va='bottom', fontsize=8)

    # 4. Genetic architecture comparison
    ax4 = axes[1, 1]
    arch_data = analyze_genetic_architecture(pop_data)

    if arch_data:
        pops = list(arch_data.keys())
        h2_estimates = [arch_data[pop].get('heritability', 0) for pop in pops]
        polygenicity = [arch_data[pop].get('polygenicity', 0) for pop in pops]

        x = np.arange(len(pops))
        width = 0.35

        ax4.bar(x - width/2, h2_estimates, width, label='Heritability', alpha=0.7)
        ax4.bar(x + width/2, polygenicity, width, label='Polygenicity', alpha=0.7)

        ax4.set_title('Genetic Architecture')
        ax4.set_xticks(x)
        ax4.set_xticklabels(pops)
        ax4.legend()
        ax4.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
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
        results = pop_info.get('results')
        if hasattr(results, 'index'):
            all_snps.update(results.index)

    for pop_name, pop_info in pop_data.items():
        results = pop_info.get('results')
        if not hasattr(results, 'columns') or 'P' not in results.columns:
            continue

        # Count significant SNPs in this population
        sig_snps = (results['P'] < 5e-8).sum()
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
            contingency = [[sig_snps, total_snps - sig_snps],
                          [int(expected_prop * total_snps), int((1 - expected_prop) * total_snps)]]
            chi2, p_value = stats.chi2_contingency(contingency)[:2]
        except ImportError:
            p_value = 0.5  # Conservative estimate

        enrichment_results[pop_name] = {
            'enrichment': enrichment,
            'p_value': p_value,
            'observed_sig': sig_snps,
            'expected_sig': int(expected_prop * total_snps)
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
        results = pop_info.get('results')
        if not hasattr(results, 'columns') or 'P' in results.columns:
            continue

        # Estimate heritability (simplified)
        p_values = results['P'].dropna().values
        lambda_gc = np.median(p_values) / 0.456  # Approximate genomic control
        h2_estimate = min(lambda_gc / (lambda_gc + 1), 1.0)  # Simplified

        # Estimate polygenicity (simplified)
        sig_snps = (p_values < 5e-8).sum()
        polygenicity = sig_snps / len(p_values)

        architecture[pop_name] = {
            'heritability': h2_estimate,
            'polygenicity': polygenicity,
            'lambda_gc': lambda_gc
        }

    return architecture


def compare_analysis_methods(method_results: Dict[str, Any],
                           output_file: Optional[str | Path] = None) -> Optional[Any]:
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
    fig.suptitle('GWAS Method Comparison', fontsize=16)

    # 1. P-value distribution comparison
    ax1 = axes[0, 0]
    for method_name, method_data in method_results.items():
        results = method_data.get('results')
        if hasattr(results, 'columns') and 'P' in results.columns:
            p_values = results['P'].dropna().values
            ax1.hist(p_values, bins=50, alpha=0.5, label=method_name, density=True)

    ax1.set_title('P-value Distribution')
    ax1.set_xlabel('P-value')
    ax1.set_ylabel('Density')
    ax1.legend()
    ax1.set_yscale('log')

    # 2. Q-Q plot comparison
    ax2 = axes[0, 1]
    for method_name, method_data in method_results.items():
        results = method_data.get('results')
        if hasattr(results, 'columns') and 'P' in results.columns:
            p_values = results['P'].dropna().values
            plot_qq_overlay(ax2, p_values, method_name, None)

    ax2.set_title('Q-Q Plot Comparison')
    ax2.plot([0, 1], [0, 1], 'k--', alpha=0.5, label='Expected')
    ax2.legend()

    # 3. Power comparison
    ax3 = axes[1, 0]
    power_data = calculate_method_power(method_results)

    if power_data:
        methods = list(power_data.keys())
        powers = [power_data[method]['power'] for method in methods]

        bars = ax3.bar(methods, powers, alpha=0.7)
        ax3.set_title('Method Power Comparison')
        ax3.set_ylabel('Statistical Power')
        ax3.set_ylim(0, 1)

        # Add value labels
        for bar, power in zip(bars, powers):
            ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                    '.3f', ha='center', va='bottom')

    # 4. Computational performance
    ax4 = axes[1, 1]
    perf_data = method_results  # Assume timing data is included

    methods = list(perf_data.keys())
    times = [perf_data[method].get('runtime_minutes', 0) for method in methods]

    if any(times):
        bars = ax4.bar(methods, times, alpha=0.7, color='orange')
        ax4.set_title('Computational Performance')
        ax4.set_ylabel('Runtime (minutes)')
        ax4.set_yscale('log')

        # Add value labels
        for bar, time_val in zip(bars, times):
            if time_val > 0:
                ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                        '.1f', ha='center', va='bottom')

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
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
        results = method_data.get('results')
        if not hasattr(results, 'columns') or 'P' in results.columns:
            continue

        # Assume some SNPs are truly associated (simplified)
        # In practice, would need ground truth data
        p_values = results['P'].dropna().values

        # Estimate power as proportion of significant findings
        # (This is a very simplified calculation)
        sig_snps = (p_values < 5e-8).sum()
        estimated_power = min(sig_snps / len(p_values), 1.0)

        power_results[method_name] = {
            'power': estimated_power,
            'significant_snps': sig_snps,
            'total_snps': len(p_values)
        }

    return power_results


