#!/usr/bin/env python3
"""Visualization functions for population genetics analysis results.

This module contains functions for generating comprehensive visualizations
of population genetics analysis results including diversity comparisons,
neutrality tests, population structure, and demographic model comparisons.
"""

from __future__ import annotations

import random
from pathlib import Path
from typing import Any

from metainformant.core.logging import setup_logger
from metainformant.dna.population_viz import (
    plot_allele_frequency_spectrum,
    plot_bootstrap_distribution,
    plot_demographic_comparison,
    plot_diversity_comparison,
    plot_fst_comparison,
    plot_fst_matrix,
    plot_hardy_weinberg_test,
    plot_heterozygosity_distribution,
    plot_kinship_matrix,
    plot_linkage_disequilibrium_decay,
    plot_neutrality_test_suite,
    plot_neutrality_test_summary,
    plot_outlier_detection,
    plot_pca_results,
    plot_pi_vs_theta,
    plot_permutation_test,
    plot_statistic_correlation_matrix,
    plot_statistic_distribution,
    plot_site_frequency_spectrum,
    plot_summary_statistics_grid,
    plot_tajimas_d_comparison,
)
from metainformant.simulation.popgen import generate_site_frequency_spectrum


def generate_visualizations(
    analysis_results: dict[str, Any],
    output_dir: Path,
) -> None:
    """Generate comprehensive visualizations for analysis results.

    Args:
        analysis_results: Analysis results dictionary
        output_dir: Output directory for plots
    """
    logger = setup_logger("metainformant.popgen.visualization")
    plots_dir = output_dir / "plots"
    from metainformant.core.io import ensure_directory
    ensure_directory(str(plots_dir))

    scenarios = analysis_results["scenario_analyses"]
    comp = analysis_results["comparative_analysis"]
    demo = analysis_results.get("demographic_model_comparisons", {})

    # 1. Diversity comparison
    logger.info("Generating diversity comparison plot")
    plot_diversity_comparison(
        comp["diversity_comparison"],
        output_path=plots_dir / "diversity_comparison.png",
        title="Nucleotide Diversity (π) Across Scenarios",
    )

    # 2. Tajima's D comparison
    logger.info("Generating Tajima's D comparison plot")
    plot_tajimas_d_comparison(
        comp["tajimas_d_comparison"],
        output_path=plots_dir / "tajimas_d_comparison.png",
        title="Tajima's D Across Scenarios",
    )

    # 3. Fst comparison
    logger.info("Generating Fst comparison plot")
    plot_fst_comparison(
        comp["fst_comparison"],
        output_path=plots_dir / "fst_comparison.png",
        title="Fst Comparison Between Populations",
    )

    # 4. Neutrality test summary
    logger.info("Generating neutrality test summary plot")
    neutrality_data = {}
    for scenario_name in ["neutral", "high_diversity", "low_diversity", "bottleneck", "expansion"]:
        if scenario_name in scenarios:
            scenario_data = scenarios[scenario_name]
            if "neutrality_tests" in scenario_data:
                neutrality_data[scenario_name.replace("_", " ").title()] = scenario_data["neutrality_tests"]

    if neutrality_data:
        plot_neutrality_test_summary(
            neutrality_data,
            output_path=plots_dir / "neutrality_test_summary.png",
            title="Neutrality Test Summary",
        )

    # 5. PCA plots (if available)
    if "large_genotypes" in scenarios and scenarios["large_genotypes"].get("pca", {}).get("status") == "success":
        logger.info("Generating PCA plots")
        pca_data = scenarios["large_genotypes"]["pca"]
        if "full_result" in pca_data:
            plot_pca_results(
                pca_data["full_result"],
                output_path=plots_dir / "pca_analysis.png",
                title="Principal Component Analysis (Large Genotype Matrix)",
                n_components=10,
            )

    # 6. Kinship matrix (if available)
    if "large_genotypes" in scenarios and scenarios["large_genotypes"].get("kinship", {}).get("status") == "success":
        logger.info("Generating kinship matrix plot")
        kinship_data = scenarios["large_genotypes"]["kinship"]
        if "full_result" in kinship_data:
            plot_kinship_matrix(
                kinship_data["full_result"],
                output_path=plots_dir / "kinship_matrix.png",
                title="Kinship Matrix (VanRaden Method)",
                max_samples=100,
            )

    # 7. Site frequency spectrum (generate example)
    logger.info("Generating site frequency spectrum example")

    example_sfs = generate_site_frequency_spectrum(
        sample_size=30,
        n_sites=100,
        theta=0.01,
        folded=True,
        rng=random.Random(42),
    )
    plot_site_frequency_spectrum(
        example_sfs,
        output_path=plots_dir / "site_frequency_spectrum_example.png",
        title="Site Frequency Spectrum (Example)",
    )

    # 8. Demographic model comparison
    if demo:
        logger.info("Generating demographic model comparison plot")
        plot_demographic_comparison(
            demo,
            output_path=plots_dir / "demographic_model_comparison.png",
            title="Demographic Model Comparison",
        )

    # 9. Summary statistics grid
    logger.info("Generating summary statistics grid")
    summary_stats_data = {}
    for scenario_name in ["neutral", "high_diversity", "low_diversity", "bottleneck", "expansion"]:
        if scenario_name in scenarios:
            scenario_data = scenarios[scenario_name]
            if "summary_statistics" in scenario_data:
                summary_stats_data[scenario_name.replace("_", " ").title()] = scenario_data["summary_statistics"]

    if summary_stats_data:
        plot_summary_statistics_grid(
            summary_stats_data,
            output_path=plots_dir / "summary_statistics_grid.png",
            title="Summary Statistics Grid",
        )

    # 10. Linkage disequilibrium decay (if available)
    if "linkage_disequilibrium" in scenarios:
        logger.info("Generating LD decay plot")
        ld_data = scenarios["linkage_disequilibrium"]
        if "r_squared_values" in ld_data and ld_data["r_squared_values"]:
            plot_linkage_disequilibrium_decay(
                ld_data["r_squared_values"],
                distances=list(range(len(ld_data["r_squared_values"]))),
                output_path=plots_dir / "linkage_disequilibrium_decay.png",
                title="Linkage Disequilibrium Decay",
            )

    # Additional statistical test visualizations
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(exist_ok=True)

    # Collect neutrality test results for comprehensive suite plot
    neutrality_test_results = {}
    for scenario_name, analysis in analysis_results["scenario_analyses"].items():
        scenario_tests = {}
        # Collect all neutrality test statistics for this scenario
        if "neutrality_tests" in analysis:
            scenario_tests.update(analysis["neutrality_tests"])
        # Add Fu & Li tests if available
        if "fu_and_li_d_star" in analysis:
            scenario_tests["fu_and_li_d_star"] = analysis["fu_and_li_d_star"]
        if "fu_and_li_f_star" in analysis:
            scenario_tests["fu_and_li_f_star"] = analysis["fu_and_li_f_star"]
        if "fay_wu_h" in analysis:
            scenario_tests["fay_wu_h"] = analysis["fay_wu_h"]

        if scenario_tests:
            neutrality_test_results[scenario_name] = scenario_tests

    if neutrality_test_results:
        logger.info("Generating comprehensive neutrality test suite plot")
        plot_neutrality_test_suite(
            neutrality_test_results,
            output_path=plots_dir / "comprehensive_neutrality_test_suite.png",
            title="Comprehensive Neutrality Test Suite - All Scenarios",
        )

    # Generate correlation matrix for statistics
    logger.info("Generating statistic correlation matrix")
    all_stats = {}
    for scenario_name, analysis in analysis_results["scenario_analyses"].items():
        if "summary_statistics" in analysis:
            stats = analysis["summary_statistics"]
            for stat_name, stat_value in stats.items():
                if isinstance(stat_value, (int, float)):
                    if stat_name not in all_stats:
                        all_stats[stat_name] = []
                    all_stats[stat_name].append(stat_value)

    if len(all_stats) > 1:
        plot_statistic_correlation_matrix(
            all_stats,
            output_path=plots_dir / "statistic_correlation_matrix.png",
        )

    # Generate π vs θ plot
    logger.info("Generating π vs θ comparison plot")
    pi_values = []
    theta_values = []
    for scenario_name, analysis in analysis_results["scenario_analyses"].items():
        if "summary_statistics" in analysis:
            stats = analysis["summary_statistics"]
            if "nucleotide_diversity" in stats and "wattersons_theta" in stats:
                pi_values.append(stats["nucleotide_diversity"])
                theta_values.append(stats["wattersons_theta"])

    if len(pi_values) > 1 and len(theta_values) > 1:
        plot_pi_vs_theta(
            pi_values,
            theta_values,
            output_path=plots_dir / "pi_vs_theta_comparison.png",
        )

    # Generate Fst matrix if multiple populations
    logger.info("Generating Fst matrix")
    fst_matrix = {}
    for scenario_name, analysis in analysis_results["scenario_analyses"].items():
        if "fst" in analysis:
            fst_matrix[scenario_name] = {"comparison": analysis["fst"]}

    if len(fst_matrix) > 1:
        plot_fst_matrix(
            fst_matrix,
            output_path=plots_dir / "fst_matrix.png",
        )

    logger.info(f"All visualizations saved to {plots_dir}")




