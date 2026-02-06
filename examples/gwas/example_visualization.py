#!/usr/bin/env python3
"""GWAS visualization example.

This example demonstrates genome-wide association study (GWAS) visualization techniques using METAINFORMANT's plotting capabilities.

Usage:
    python examples/gwas/example_visualization.py

Expected output:
    output/examples/gwas/gwas_plots.json
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from metainformant.core import io
from metainformant.gwas.visualization import manhattan_plot, qq_plot


def main():
    """Demonstrate GWAS visualization techniques."""
    # Setup output directory
    output_dir = Path("output/examples/gwas")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== METAINFORMANT GWAS Visualization Example ===")

    # 1. Create simulated GWAS results
    print("\n1. Creating simulated GWAS results...")

    np.random.seed(42)

    n_snps = 5000  # Simulate 5000 SNPs across genome

    # Create SNP positions across 5 chromosomes
    chromosomes = []
    positions = []
    current_pos = 0

    # Fixed chromosome sizes to ensure total equals n_snps
    chr_sizes = [1000, 1000, 1000, 1000, 1000]  # 5000 total SNPs

    for chr_num, chr_size in enumerate(chr_sizes, 1):
        chromosomes.extend([f"chr{chr_num}"] * chr_size)
        positions.extend(range(current_pos, current_pos + chr_size))
        current_pos += chr_size + 100  # Gap between chromosomes

    # Generate p-values: most null, some associated
    p_values = np.random.uniform(0, 1, n_snps)

    # Add some significant associations
    significant_indices = np.random.choice(n_snps, size=50, replace=False)
    p_values[significant_indices] = np.random.uniform(1e-15, 1e-8, len(significant_indices))

    # Add genome-wide significant hits
    gw_sig_indices = significant_indices[:5]
    p_values[gw_sig_indices] = np.random.uniform(1e-20, 1e-15, len(gw_sig_indices))

    # Create effect sizes
    effect_sizes = np.random.normal(0, 0.1, n_snps)
    effect_sizes[significant_indices] += np.random.normal(0, 0.3, len(significant_indices))

    # Create DataFrame for visualization
    gwas_results = pd.DataFrame(
        {
            "CHR": chromosomes,
            "BP": positions,
            "P": p_values,
            "BETA": effect_sizes,
            "SNP": [f"rs{i}" for i in range(n_snps)],
        }
    )

    print(f"✓ Created GWAS results for {n_snps} SNPs across {len(set(chromosomes))} chromosomes")
    print(f"  Genome-wide significant SNPs: {len(gw_sig_indices)}")
    print(f"  Suggestive SNPs: {len(significant_indices) - len(gw_sig_indices)}")

    # 2. Manhattan plot
    print("\n2. Generating Manhattan plot...")

    try:
        manhattan_fig = manhattan_plot(
            gwas_results,
            significance_threshold=5e-8,  # Genome-wide significance
            suggestive_threshold=1e-5,  # Suggestive significance
            title="GWAS Manhattan Plot - Simulated Data",
        )

        manhattan_path = output_dir / "manhattan_plot.png"
        manhattan_fig.savefig(manhattan_path, dpi=300, bbox_inches="tight")
        manhattan_fig.savefig(output_dir / "manhattan_plot.pdf", bbox_inches="tight")
        plt.close(manhattan_fig)

        print(f"✓ Manhattan plot saved to: {manhattan_path}")

        manhattan_info = {
            "file": str(manhattan_path.relative_to(output_dir)),
            "significance_threshold": 5e-8,
            "suggestive_threshold": 1e-5,
            "snps_plotted": len(gwas_results),
            "chromosomes": sorted(list(set(gwas_results["CHR"]))),
            "genome_wide_significant": len(gwas_results[gwas_results["P"] < 5e-8]),
            "suggestive_significant": len(gwas_results[(gwas_results["P"] < 1e-5) & (gwas_results["P"] >= 5e-8)]),
        }

    except Exception as e:
        print(f"Manhattan plot failed: {e}")
        manhattan_info = {"error": str(e)}

    # 3. Q-Q plot
    print("\n3. Generating Q-Q plot...")

    try:
        qq_fig = qq_plot(gwas_results["P"].values, title="GWAS Q-Q Plot - Simulated Data")

        qq_path = output_dir / "qq_plot.png"
        qq_fig.savefig(qq_path, dpi=300, bbox_inches="tight")
        qq_fig.savefig(output_dir / "qq_plot.pdf", bbox_inches="tight")
        plt.close(qq_fig)

        print(f"✓ Q-Q plot saved to: {qq_path}")

        # Calculate lambda GC (genomic control inflation factor)
        # Simplified calculation
        observed_median = np.median(-np.log10(gwas_results["P"]))
        expected_median = np.median(-np.log10(np.random.uniform(0, 1, len(gwas_results))))
        lambda_gc = observed_median / expected_median

        qq_info = {
            "file": str(qq_path.relative_to(output_dir)),
            "lambda_gc": lambda_gc,
            "interpretation": "λ > 1 suggests test statistic inflation (population structure, etc.)",
            "snps_analyzed": len(gwas_results),
            "inflation_status": "minimal" if lambda_gc < 1.05 else "moderate" if lambda_gc < 1.1 else "severe",
        }

        print(f"  Genomic control λ: {lambda_gc:.3f} ({qq_info['inflation_status']} inflation)")

    except Exception as e:
        print(f"Q-Q plot failed: {e}")
        qq_info = {"error": str(e)}

    # 4. Effect size analysis
    print("\n4. Analyzing effect sizes...")

    effect_analysis = {
        "summary_stats": {
            "mean_beta": gwas_results["BETA"].mean(),
            "median_beta": gwas_results["BETA"].median(),
            "std_beta": gwas_results["BETA"].std(),
            "min_beta": gwas_results["BETA"].min(),
            "max_beta": gwas_results["BETA"].max(),
        },
        "significant_effects": {
            "gw_sig_mean_beta": gwas_results[gwas_results["P"] < 5e-8]["BETA"].mean(),
            "gw_sig_median_beta": gwas_results[gwas_results["P"] < 5e-8]["BETA"].median(),
            "gw_sig_count": len(gwas_results[gwas_results["P"] < 5e-8]),
        },
    }

    print("  Effect size distribution:")
    print(f"    Mean β: {effect_analysis['summary_stats']['mean_beta']:.3f}")
    print(f"    GW-significant mean β: {effect_analysis['significant_effects']['gw_sig_mean_beta']:.3f}")

    # 5. Chromosome-wise statistics
    print("\n5. Calculating chromosome-wise statistics...")

    chr_stats = []

    for chr_name in sorted(set(gwas_results["CHR"])):
        chr_data = gwas_results[gwas_results["CHR"] == chr_name]
        chr_stat = {
            "chromosome": chr_name,
            "n_snps": len(chr_data),
            "min_p": chr_data["P"].min(),
            "gw_significant": len(chr_data[chr_data["P"] < 5e-8]),
            "suggestive": len(chr_data[(chr_data["P"] < 1e-5) & (chr_data["P"] >= 5e-8)]),
            "mean_beta": chr_data["BETA"].mean(),
            "median_p": chr_data["P"].median(),
        }
        chr_stats.append(chr_stat)

        if chr_stat["gw_significant"] > 0:
            print(f"  {chr_name}: {chr_stat['n_snps']} SNPs, {chr_stat['gw_significant']} GW-significant")

    # 6. Visualization metadata
    print("\n6. Creating visualization metadata...")

    visualization_metadata = {
        "plot_files": {
            "manhattan_plot": "manhattan_plot.png",
            "qq_plot": "qq_plot.png",
            "vector_formats": ["manhattan_plot.pdf", "qq_plot.pdf"],
        },
        "plot_characteristics": {
            "manhattan": {
                "type": "Genome-wide association results",
                "x_axis": "Chromosomal position",
                "y_axis": "-log10(p-value)",
                "color_scheme": "Alternating chromosome colors",
                "thresholds": {"genome_wide": 5e-8, "suggestive": 1e-5},
            },
            "qq": {
                "type": "Expected vs observed p-value distribution",
                "x_axis": "Expected -log10(p)",
                "y_axis": "Observed -log10(p)",
                "diagonal": "Null expectation line",
                "lambda_gc": "Genomic control inflation factor",
            },
        },
    }

    # 7. Create comprehensive visualization results
    print("\n7. Creating comprehensive visualization analysis...")

    visualization_results = {
        "example": "gwas_visualization",
        "domain": "gwas",
        "results": {
            "timestamp": "2024-12-26T10:00:00Z",
            "dataset_summary": {
                "total_snps": len(gwas_results),
                "chromosomes": sorted(list(set(gwas_results["CHR"]))),
                "genome_wide_significant": len(gwas_results[gwas_results["P"] < 5e-8]),
                "suggestive_significant": len(gwas_results[(gwas_results["P"] < 1e-5) & (gwas_results["P"] >= 5e-8)]),
                "lambda_gc": qq_info.get("lambda_gc", "N/A"),
            },
            "plots_generated": {"manhattan_plot": manhattan_info, "qq_plot": qq_info},
            "effect_size_analysis": effect_analysis,
            "chromosome_statistics": chr_stats,
            "visualization_metadata": visualization_metadata,
            "gwas_visualization_best_practices": [
                "Use Manhattan plots to visualize genome-wide associations",
                "Include Q-Q plots to check for test statistic inflation",
                "Set appropriate significance thresholds",
                "Use consistent color schemes across plots",
                "Save both raster (PNG) and vector (PDF) formats",
                "Include lambda GC values to assess inflation",
                "Consider regional plots for significant loci",
            ],
            "interpretation_guidelines": [
                "Manhattan plot peaks indicate associated regions",
                "Q-Q plot deviation from diagonal suggests inflation",
                "λ > 1.05 may indicate population structure or other artifacts",
                "Effect sizes help prioritize variants for follow-up",
                "Chromosome-specific patterns may indicate technical artifacts",
            ],
        },
    }

    results_file = output_dir / "gwas_plots.json"
    io.dump_json(visualization_results, results_file, indent=2)

    print(f"✓ Comprehensive GWAS visualization analysis saved to: {results_file}")

    print("\n=== GWAS Visualization Example Complete ===")
    print("This example demonstrated METAINFORMANT's GWAS visualization capabilities:")
    print("- Manhattan plot generation with significance thresholds")
    print("- Q-Q plot for distribution checking and inflation assessment")
    print("- Effect size analysis and chromosome-wise statistics")
    print("- Best practices for GWAS data visualization")

    print(f"\nAll outputs saved to: {output_dir}")


if __name__ == "__main__":
    main()
