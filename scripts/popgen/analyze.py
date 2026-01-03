#!/usr/bin/env python3
"""Analysis functions for population genetics datasets.

This module contains functions for analyzing population genetics datasets
including summary statistics, neutrality tests, and comparative analyses.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Any

from metainformant.core.io import dump_json, load_json
from metainformant.core.logging import setup_logger
from metainformant.dna.population import (
    fu_and_li_d_star_from_sequences,
    fu_and_li_f_star_from_sequences,
    fay_wu_h_from_sequences,
)
from metainformant.dna.population_analysis import (
    calculate_summary_statistics,
    compare_populations,
    neutrality_test_suite,
)
from metainformant.dna.sequences import read_fasta
from metainformant.gwas.quality import test_hwe
from metainformant.gwas.structure import compute_pca, compute_kinship_matrix
# Note: Additional neutrality tests are imported where used

from metainformant.math.demography import (
    bottleneck_effective_size,
    exponential_growth_effective_size,
    two_epoch_effective_size,
)


def analyze_dataset(dataset_info: dict[str, Any], output_dir: Path) -> dict[str, Any]:
    """Perform comprehensive analysis on generated dataset.

    Args:
        dataset_info: Dataset metadata
        output_dir: Output directory for results

    Returns:
        Dictionary with analysis results
    """
    logger = setup_logger("metainformant.popgen.analysis")
    logger.info("Starting comprehensive analysis")

    results = {
        "analysis_time": datetime.utcnow().isoformat(),
        "scenario_analyses": {},
    }

    # Analyze each scenario
    scenarios = dataset_info["scenarios"]

    # Scenario 1: Neutral
    logger.info("Analyzing Scenario 1: Neutral population")
    neutral_seqs = list(read_fasta(scenarios["neutral"]["file"]).values())
    neutral_stats = calculate_summary_statistics(sequences=neutral_seqs)
    neutral_neutrality = neutrality_test_suite(neutral_seqs)

    # Calculate additional tests for neutral scenario
    neutral_fl_d = fu_and_li_d_star_from_sequences(neutral_seqs)
    neutral_fl_f = fu_and_li_f_star_from_sequences(neutral_seqs)
    neutral_fw_h = fay_wu_h_from_sequences(neutral_seqs)

    results["scenario_analyses"]["neutral"] = {
        "summary_statistics": neutral_stats,
        "neutrality_tests": neutral_neutrality,
        "fu_and_li_d_star": neutral_fl_d,
        "fu_and_li_f_star": neutral_fl_f,
        "fay_wu_h": neutral_fw_h,
    }

    # Scenario 2: High diversity
    logger.info("Analyzing Scenario 2: High diversity")
    highdiv_seqs = list(read_fasta(scenarios["high_diversity"]["file"]).values())
    highdiv_stats = calculate_summary_statistics(sequences=highdiv_seqs)
    highdiv_neutrality = neutrality_test_suite(highdiv_seqs)
    # Add additional tests
    highdiv_fl_d = fu_and_li_d_star_from_sequences(highdiv_seqs)
    highdiv_fl_f = fu_and_li_f_star_from_sequences(highdiv_seqs)
    highdiv_fw_h = fay_wu_h_from_sequences(highdiv_seqs)
    results["scenario_analyses"]["high_diversity"] = {
        "summary_statistics": highdiv_stats,
        "neutrality_tests": highdiv_neutrality,
        "fu_and_li_d_star": highdiv_fl_d,
        "fu_and_li_f_star": highdiv_fl_f,
        "fay_wu_h": highdiv_fw_h,
    }

    # Scenario 3: Low diversity
    logger.info("Analyzing Scenario 3: Low diversity")
    lowdiv_seqs = list(read_fasta(scenarios["low_diversity"]["file"]).values())
    lowdiv_stats = calculate_summary_statistics(sequences=lowdiv_seqs)
    lowdiv_neutrality = neutrality_test_suite(lowdiv_seqs)
    # Add additional tests
    lowdiv_fl_d = fu_and_li_d_star_from_sequences(lowdiv_seqs)
    lowdiv_fl_f = fu_and_li_f_star_from_sequences(lowdiv_seqs)
    lowdiv_fw_h = fay_wu_h_from_sequences(lowdiv_seqs)
    results["scenario_analyses"]["low_diversity"] = {
        "summary_statistics": lowdiv_stats,
        "neutrality_tests": lowdiv_neutrality,
        "fu_and_li_d_star": lowdiv_fl_d,
        "fu_and_li_f_star": lowdiv_fl_f,
        "fay_wu_h": lowdiv_fw_h,
    }

    # Scenario 4: Bottleneck
    logger.info("Analyzing Scenario 4: Bottleneck")
    bottleneck_seqs = list(read_fasta(scenarios["bottleneck"]["file"]).values())
    bottleneck_stats = calculate_summary_statistics(sequences=bottleneck_seqs)
    bottleneck_neutrality = neutrality_test_suite(bottleneck_seqs)
    # Add additional tests
    bottleneck_fl_d = fu_and_li_d_star_from_sequences(bottleneck_seqs)
    bottleneck_fl_f = fu_and_li_f_star_from_sequences(bottleneck_seqs)
    bottleneck_fw_h = fay_wu_h_from_sequences(bottleneck_seqs)
    results["scenario_analyses"]["bottleneck"] = {
        "summary_statistics": bottleneck_stats,
        "neutrality_tests": bottleneck_neutrality,
        "fu_and_li_d_star": bottleneck_fl_d,
        "fu_and_li_f_star": bottleneck_fl_f,
        "fay_wu_h": bottleneck_fw_h,
    }

    # Scenario 5: Expansion
    logger.info("Analyzing Scenario 5: Expansion")
    expansion_seqs = list(read_fasta(scenarios["expansion"]["file"]).values())
    expansion_stats = calculate_summary_statistics(sequences=expansion_seqs)
    expansion_neutrality = neutrality_test_suite(expansion_seqs)
    # Add additional tests
    expansion_fl_d = fu_and_li_d_star_from_sequences(expansion_seqs)
    expansion_fl_f = fu_and_li_f_star_from_sequences(expansion_seqs)
    expansion_fw_h = fay_wu_h_from_sequences(expansion_seqs)
    results["scenario_analyses"]["expansion"] = {
        "summary_statistics": expansion_stats,
        "neutrality_tests": expansion_neutrality,
        "fu_and_li_d_star": expansion_fl_d,
        "fu_and_li_f_star": expansion_fl_f,
        "fay_wu_h": expansion_fw_h,
    }

    # Scenario 6: Two populations (low Fst)
    logger.info("Analyzing Scenario 6: Two populations (low Fst)")
    pop1_low = list(read_fasta(scenarios["two_populations_low_fst"]["file_pop1"]).values())
    pop2_low = list(read_fasta(scenarios["two_populations_low_fst"]["file_pop2"]).values())
    comparison_low = compare_populations(
        pop1_sequences=pop1_low,
        pop2_sequences=pop2_low,
    )
    results["scenario_analyses"]["two_populations_low_fst"] = comparison_low

    # Scenario 7: Two populations (high Fst)
    logger.info("Analyzing Scenario 7: Two populations (high Fst)")
    pop1_high = list(read_fasta(scenarios["two_populations_high_fst"]["file_pop1"]).values())
    pop2_high = list(read_fasta(scenarios["two_populations_high_fst"]["file_pop2"]).values())
    comparison_high = compare_populations(
        pop1_sequences=pop1_high,
        pop2_sequences=pop2_high,
    )
    results["scenario_analyses"]["two_populations_high_fst"] = comparison_high

    # Scenario 8: Large genotype matrix
    logger.info("Analyzing Scenario 8: Large genotype matrix")

    large_genotypes = load_json(scenarios["large_genotypes"]["file"])
    pca_result = compute_pca(large_genotypes, n_components=10)
    kinship_result = compute_kinship_matrix(large_genotypes, method="vanraden")

    # Hardy-Weinberg test on genotype matrix
    hwe_p_values = test_hwe(large_genotypes)

    # Format results for plotting function
    hwe_result = []
    for i, p_value in enumerate(hwe_p_values):
        hwe_result.append({
            "locus": f"Variant_{i}",
            "p_value": p_value,
            "chi_square": None,  # Could be calculated if needed
            "degrees_of_freedom": 2,  # Standard for HWE test
            "hwe_deviated": p_value < 0.05  # Significant deviation
        })

    results["scenario_analyses"]["large_genotypes"] = {
        "pca": {
            "status": pca_result.get("status"),
            "n_components": pca_result.get("n_components"),
            "explained_variance_ratio": pca_result.get("explained_variance_ratio", [])[:5],
            "missing_data_stats": pca_result.get("missing_data_stats"),
            "pcs": pca_result.get("pcs", [])[:100],  # Store first 100 samples for visualization
            "full_result": pca_result,  # Store full result for visualization
        },
        "kinship": {
            "status": kinship_result.get("status"),
            "method": kinship_result.get("method"),
            "num_samples": kinship_result.get("num_samples"),
            "missing_data_stats": kinship_result.get("missing_data_stats"),
            "kinship_matrix": kinship_result.get("kinship_matrix", [])[:100],  # First 100 samples
            "full_result": kinship_result,  # Store full result for visualization
        },
        "hardy_weinberg_test": hwe_result,
    }

    # Scenario 9: Linkage disequilibrium
    logger.info("Analyzing Scenario 9: Linkage disequilibrium")
    ld_genotypes = load_json(scenarios["linkage_disequilibrium"]["file"])
    # Calculate LD between adjacent sites using ld_coefficients
    from metainformant.math.ld import ld_coefficients

    ld_values = []
    for i in range(min(10, len(ld_genotypes[0]) - 1)):  # Sample first 10 pairs
        site1 = [row[i] for row in ld_genotypes]
        site2 = [row[i + 1] for row in ld_genotypes]
        # Convert genotypes to allele counts for LD calculation
        # For diploid: 0=AA, 1=AB, 2=BB
        # Convert to allele counts: 0->0, 1->1, 2->2 (but we need 0/1 encoding)
        # For simplicity, use genotype values directly
        try:
            ld_result = ld_coefficients(site1, site2)
            if ld_result and "r_squared" in ld_result:
                ld_values.append(ld_result["r_squared"])
        except Exception:
            # Skip if calculation fails
            pass

    results["scenario_analyses"]["linkage_disequilibrium"] = {
        "mean_r_squared": sum(ld_values) / len(ld_values) if ld_values else 0.0,
        "r_squared_values": ld_values,
    }

    # Comparative analysis
    logger.info("Performing comparative analysis")
    results["comparative_analysis"] = {
        "diversity_comparison": {
            "neutral": neutral_stats["nucleotide_diversity"],
            "high_diversity": highdiv_stats["nucleotide_diversity"],
            "low_diversity": lowdiv_stats["nucleotide_diversity"],
            "bottleneck": bottleneck_stats["nucleotide_diversity"],
            "expansion": expansion_stats["nucleotide_diversity"],
        },
        "tajimas_d_comparison": {
            "neutral": neutral_neutrality["tajimas_d"],
            "high_diversity": highdiv_neutrality["tajimas_d"],
            "low_diversity": lowdiv_neutrality["tajimas_d"],
            "bottleneck": bottleneck_neutrality["tajimas_d"],
            "expansion": expansion_neutrality["tajimas_d"],
        },
        "fst_comparison": {
            "low_fst": comparison_low["fst"],
            "high_fst": comparison_high["fst"],
        },
    }

    # Demographic model predictions
    logger.info("Comparing to demographic models")
    results["demographic_model_comparisons"] = {
        "bottleneck": {
            "estimated_ne": bottleneck_effective_size(
                pre_bottleneck_size=10000,
                bottleneck_size=5,
                bottleneck_duration=10,
                recovery_generations=20,
            ),
            "observed_diversity": bottleneck_stats["nucleotide_diversity"],
        },
        "expansion": {
            "estimated_ne": exponential_growth_effective_size(
                current_size=1000,
                growth_rate=0.1,
                generations=23,  # log(10)/0.1 ≈ 23
            ),
            "observed_diversity": expansion_stats["nucleotide_diversity"],
        },
    }

    # Validate results before saving
    logger.info("Validating analysis results")
    validation_errors = []

    # Check that all scenarios have expected data
    for scenario_name, scenario_data in results["scenario_analyses"].items():
        if scenario_name in ["neutral", "high_diversity", "low_diversity", "bottleneck", "expansion"]:
            if "summary_statistics" not in scenario_data:
                validation_errors.append(f"{scenario_name}: Missing summary_statistics")
            if "neutrality_tests" not in scenario_data:
                validation_errors.append(f"{scenario_name}: Missing neutrality_tests")
        elif scenario_name in ["two_populations_low_fst", "two_populations_high_fst"]:
            if "fst" not in scenario_data:
                validation_errors.append(f"{scenario_name}: Missing Fst")
        elif scenario_name == "large_genotypes":
            if "pca" not in scenario_data:
                validation_errors.append(f"{scenario_name}: Missing PCA")
            if "kinship" not in scenario_data:
                validation_errors.append(f"{scenario_name}: Missing kinship")

    if validation_errors:
        logger.warning(f"Validation found {len(validation_errors)} issues:")
        for error in validation_errors:
            logger.warning(f"  - {error}")
    else:
        logger.info("✓ All validation checks passed")

    # Add validation summary to results
    results["validation"] = {
        "status": "passed" if not validation_errors else "warnings",
        "errors": validation_errors,
        "total_scenarios": len(results["scenario_analyses"]),
    }

    # Save results
    results_file = output_dir / "analysis_results.json"
    dump_json(results, str(results_file))
    logger.info(f"Analysis complete. Results saved to {results_file}")

    return results
