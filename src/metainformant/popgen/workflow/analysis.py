"""Population genetics analysis pipeline methods.

Extracted from ``scripts/popgen/analysis.py`` (scripts were importing the
defunct ``metainformant.dna.population_analysis`` and
``metainformant.simulation.popgen`` module paths; both are re-pointed at
their real homes: ``metainformant.dna.population`` and
``metainformant.simulation.models.popgen``).

All statistics are descriptive. Inferential cross-species claims remain
gated behind the evidence-manifest freeze (see repo RNA boundary policy).
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from metainformant.core.io import dump_json, load_json
from metainformant.core.utils.logging import setup_logger
from metainformant.dna.population.analysis import (
    calculate_summary_statistics,
    compare_populations,
    neutrality_test_suite,
)
from metainformant.dna.population.core import (
    fu_and_li_d_star_from_sequences,
    fu_and_li_f_star_from_sequences,
    fay_wu_h_from_sequences,
)
from metainformant.dna.sequence.core import read_fasta
from metainformant.gwas.analysis.quality import test_hwe
from metainformant.gwas.analysis.structure import compute_kinship_matrix, compute_pca
from metainformant.math.population_genetics.demography import (
    bottleneck_effective_size,
    exponential_growth_effective_size,
)
from metainformant.simulation.models.popgen import (
    generate_genotype_matrix,
    generate_linkage_disequilibrium_data,
    generate_population_sequences,
    generate_site_frequency_spectrum,
    generate_two_populations,
    simulate_bottleneck_population,
    simulate_population_expansion,
)

__all__ = [
    "generate_genotype_matrix",
    "generate_linkage_disequilibrium_data",
    "generate_population_sequences",
    "generate_site_frequency_spectrum",
    "generate_two_populations",
    "simulate_bottleneck_population",
    "simulate_population_expansion",
]

logger = setup_logger("metainformant.popgen.workflow.analysis")

_SEQUENCE_SCENARIOS = (
    "neutral",
    "high_diversity",
    "low_diversity",
    "bottleneck",
    "expansion",
)


def _utc_now_iso() -> str:
    """Return a timezone-aware UTC timestamp (naive utcnow is deprecated)."""
    return datetime.now(timezone.utc).isoformat()


def summarize_scenario(sequences: list[str], label: str) -> dict[str, Any]:
    """Run the full descriptive statistics suite for one sequence set.

    Args:
        sequences: DNA sequence strings from one population/scenario.
        label: Scenario name used in log lines.

    Returns:
        Dict with ``summary_statistics``, ``neutrality_tests``,
        ``fu_and_li_d_star``, ``fu_and_li_f_star``, and ``fay_wu_h``.
    """
    logger.info(f"Analyzing scenario: {label}")
    stats = calculate_summary_statistics(sequences=sequences)
    neutrality = neutrality_test_suite(sequences)
    return {
        "summary_statistics": stats,
        "neutrality_tests": neutrality,
        "fu_and_li_d_star": fu_and_li_d_star_from_sequences(sequences),
        "fu_and_li_f_star": fu_and_li_f_star_from_sequences(sequences),
        "fay_wu_h": fay_wu_h_from_sequences(sequences),
    }


def sequence_scenario_suite(fasta_path: str | Path, label: str) -> dict[str, Any]:
    """Load a FASTA file and run :func:`summarize_scenario` on its records.

    Args:
        fasta_path: Path to a FASTA file of one population's sequences.
        label: Scenario name used in log lines.

    Returns:
        The scenario suite dict from :func:`summarize_scenario`.
    """
    sequences = list(read_fasta(str(fasta_path)).values())
    return summarize_scenario(sequences, label)


def compare_two_population_sequences(
    pop1_fasta: str | Path,
    pop2_fasta: str | Path,
) -> dict[str, Any]:
    """Compare two populations from FASTA files (descriptive Fst-based).

    Computes Hudson's Fst directly from the sequences and adds the
    per-population summary statistics.

    Args:
        pop1_fasta: FASTA file for population 1.
        pop2_fasta: FASTA file for population 2.

    Returns:
        Dict with ``fst``, per-population ``summary_statistics``, and the
        statistic comparison from
        :func:`metainformant.dna.population.analysis.compare_populations`.
    """
    from metainformant.dna.population.analysis import calculate_fst

    pop1 = list(read_fasta(str(pop1_fasta)).values())
    pop2 = list(read_fasta(str(pop2_fasta)).values())

    fst = calculate_fst(pop1, pop2)
    pop1_stats = calculate_summary_statistics(sequences=pop1)
    pop2_stats = calculate_summary_statistics(sequences=pop2)
    comparison = compare_populations(pop1_data=pop1_stats, pop2_data=pop2_stats)

    result: dict[str, Any] = {
        "fst": fst,
        "pop1_stats": pop1_stats,
        "pop2_stats": pop2_stats,
    }
    result.update(comparison)
    return result


def genotype_structure_analysis(
    genotype_matrix: list[list[int]],
    n_components: int = 10,
    kinship_method: str = "vanraden",
) -> dict[str, Any]:
    """Run PCA, kinship, and HWE analyses on a genotype matrix.

    Args:
        genotype_matrix: Individuals x sites dosage matrix.
        n_components: Number of PCA components.
        kinship_method: Kinship estimator name (default VanRaden).

    Returns:
        Dict with ``pca``, ``kinship``, and ``hardy_weinberg_test`` keys.
    """
    pca_result = compute_pca(genotype_matrix, n_components=n_components)
    kinship_result = compute_kinship_matrix(genotype_matrix, method=kinship_method)
    hwe_p_values = test_hwe(genotype_matrix)

    hwe_result = [
        {
            "locus": f"Variant_{i}",
            "p_value": p_value,
            "chi_square": None,
            "degrees_of_freedom": 2,
            "hwe_deviated": p_value < 0.05,
        }
        for i, p_value in enumerate(hwe_p_values)
    ]

    return {
        "pca": {
            "status": pca_result.get("status"),
            "n_components": pca_result.get("n_components"),
            "explained_variance_ratio": pca_result.get("explained_variance_ratio", [])[:5],
            "missing_data_stats": pca_result.get("missing_data_stats"),
            "pcs": pca_result.get("pcs", [])[:100],
            "full_result": pca_result,
        },
        "kinship": {
            "status": kinship_result.get("status"),
            "method": kinship_result.get("method"),
            "num_samples": kinship_result.get("num_samples"),
            "missing_data_stats": kinship_result.get("missing_data_stats"),
            "kinship_matrix": kinship_result.get("kinship_matrix", [])[:100],
            "full_result": kinship_result,
        },
        "hardy_weinberg_test": hwe_result,
    }


def ld_summary(genotype_matrix: list[list[int]], max_pairs: int = 10) -> dict[str, Any]:
    """Summarize linkage disequilibrium between adjacent genotype sites.

    Uses pairwise allele-frequency LD (r-squared) computed directly from
    the dosage columns, which is the format produced by
    :func:`metainformant.simulation.models.popgen.generate_linkage_disequilibrium_data`.

    Args:
        genotype_matrix: Individuals x sites dosage matrix (0/1/2).
        max_pairs: Maximum adjacent site pairs to test.

    Returns:
        Dict with ``mean_r_squared`` and per-pair ``r_squared_values``.
    """
    import numpy as np

    ld_values: list[float] = []
    n_sites = len(genotype_matrix[0]) if genotype_matrix else 0
    for i in range(min(max_pairs, n_sites - 1)):
        site1 = np.asarray([row[i] for row in genotype_matrix], dtype=float)
        site2 = np.asarray([row[i + 1] for row in genotype_matrix], dtype=float)
        # Skip degenerate pairs (invariant sites -> undefined r)
        if site1.std() == 0 or site2.std() == 0:
            continue
        # Allele frequencies of the "1" allele
        p1 = site1.mean() / 2.0
        p2 = site2.mean() / 2.0
        if p1 in (0.0, 1.0) or p2 in (0.0, 1.0):
            continue
        # Diploid dosage correlation coefficient equals r (the LD correlation)
        r = float(np.corrcoef(site1, site2)[0, 1])
        ld_values.append(r * r)
    return {
        "mean_r_squared": sum(ld_values) / len(ld_values) if ld_values else 0.0,
        "r_squared_values": ld_values,
    }


def demographic_model_comparisons(
    bottleneck_diversity: float,
    expansion_diversity: float,
    pre_bottleneck_size: float = 10000,
    bottleneck_size: float = 5,
    post_bottleneck_size: float = 1000,
    recovery_generations: int = 20,
    initial_size: float = 1000,
    growth_rate: float = 0.1,
    expansion_generations: int = 23,
) -> dict[str, Any]:
    """Compare observed diversities to simple demographic model predictions.

    Args:
        bottleneck_diversity: Observed nucleotide diversity post-bottleneck.
        expansion_diversity: Observed nucleotide diversity post-expansion.
        pre_bottleneck_size: Pre-bottleneck effective size.
        bottleneck_size: Bottleneck effective size.
        post_bottleneck_size: Effective size after recovery.
        recovery_generations: Recovery period in generations.
        initial_size: Effective size at expansion onset.
        growth_rate: Exponential growth rate.
        expansion_generations: Generations since expansion onset.

    Returns:
        Dict with per-model estimated Ne and observed diversity.
    """
    return {
        "bottleneck": {
            "estimated_ne": bottleneck_effective_size(
                pre_bottleneck_size=pre_bottleneck_size,
                bottleneck_size=bottleneck_size,
                post_bottleneck_size=post_bottleneck_size,
                recovery_generations=recovery_generations,
            ),
            "observed_diversity": bottleneck_diversity,
        },
        "expansion": {
            "estimated_ne": exponential_growth_effective_size(
                initial_size=initial_size,
                growth_rate=growth_rate,
                generations=expansion_generations,
            ),
            "observed_diversity": expansion_diversity,
        },
    }


def _validate_results(results: dict[str, Any]) -> list[str]:
    """Return a list of validation issues found in an analysis result dict."""
    validation_errors: list[str] = []
    for scenario_name, scenario_data in results["scenario_analyses"].items():
        if scenario_name in _SEQUENCE_SCENARIOS:
            if "summary_statistics" not in scenario_data:
                validation_errors.append(f"{scenario_name}: Missing summary_statistics")
            if "neutrality_tests" not in scenario_data:
                validation_errors.append(f"{scenario_name}: Missing neutrality_tests")
        elif scenario_name.startswith("two_populations_"):
            if "fst" not in scenario_data:
                validation_errors.append(f"{scenario_name}: Missing Fst")
        elif scenario_name == "large_genotypes":
            if "pca" not in scenario_data:
                validation_errors.append(f"{scenario_name}: Missing PCA")
            if "kinship" not in scenario_data:
                validation_errors.append(f"{scenario_name}: Missing kinship")
    return validation_errors


def analyze_dataset(dataset_info: dict[str, Any], output_dir: Path) -> dict[str, Any]:
    """Perform comprehensive analysis on a generated dataset.

    Runs every scenario declared in ``dataset_info["scenarios"]`` and
    assembles comparative + demographic summaries.

    Args:
        dataset_info: Dataset metadata (from generate_dataset.py workflow).
        output_dir: Output directory for results JSON.

    Returns:
        Dictionary with per-scenario analyses, comparisons, validation.
    """
    logger.info("Starting comprehensive analysis")
    results: dict[str, Any] = {
        "analysis_time": _utc_now_iso(),
        "scenario_analyses": {},
    }

    scenarios = dataset_info["scenarios"]

    for scenario_name in _SEQUENCE_SCENARIOS:
        if scenario_name in scenarios:
            results["scenario_analyses"][scenario_name] = sequence_scenario_suite(
                scenarios[scenario_name]["file"], scenario_name
            )

    for scenario_name in ("two_populations_low_fst", "two_populations_high_fst"):
        if scenario_name in scenarios:
            logger.info(f"Comparing populations for {scenario_name}")
            results["scenario_analyses"][scenario_name] = compare_two_population_sequences(
                scenarios[scenario_name]["file_pop1"],
                scenarios[scenario_name]["file_pop2"],
            )

    if "large_genotypes" in scenarios:
        logger.info("Analyzing large genotype matrix")
        large_genotypes = load_json(scenarios["large_genotypes"]["file"])
        results["scenario_analyses"]["large_genotypes"] = genotype_structure_analysis(large_genotypes)

    if "linkage_disequilibrium" in scenarios:
        logger.info("Analyzing linkage disequilibrium")
        ld_genotypes = load_json(scenarios["linkage_disequilibrium"]["file"])
        results["scenario_analyses"]["linkage_disequilibrium"] = ld_summary(ld_genotypes)

    # Comparative analysis
    logger.info("Performing comparative analysis")
    seq_analyses = {
        name: results["scenario_analyses"][name] for name in _SEQUENCE_SCENARIOS if name in results["scenario_analyses"]
    }
    results["comparative_analysis"] = {
        "diversity_comparison": {
            name: data["summary_statistics"]["nucleotide_diversity"] for name, data in seq_analyses.items()
        },
        "tajimas_d_comparison": {name: data["neutrality_tests"]["tajima_d"] for name, data in seq_analyses.items()},
        "fst_comparison": {
            name.replace("two_populations_", ""): results["scenario_analyses"][name]["fst"]
            for name in ("two_populations_low_fst", "two_populations_high_fst")
            if name in results["scenario_analyses"]
        },
    }

    # Demographic model predictions
    logger.info("Comparing to demographic models")
    bottleneck_stats = seq_analyses.get("bottleneck", {}).get("summary_statistics", {})
    expansion_stats = seq_analyses.get("expansion", {}).get("summary_statistics", {})
    if bottleneck_stats and expansion_stats:
        results["demographic_model_comparisons"] = demographic_model_comparisons(
            bottleneck_diversity=bottleneck_stats["nucleotide_diversity"],
            expansion_diversity=expansion_stats["nucleotide_diversity"],
        )

    # Validate results
    logger.info("Validating analysis results")
    validation_errors = _validate_results(results)
    if validation_errors:
        logger.warning(f"Validation found {len(validation_errors)} issues:")
        for error in validation_errors:
            logger.warning(f"  - {error}")
    else:
        logger.info("All validation checks passed")

    results["validation"] = {
        "status": "passed" if not validation_errors else "warnings",
        "errors": validation_errors,
        "total_scenarios": len(results["scenario_analyses"]),
    }

    # Save results
    results_file = Path(output_dir) / "analysis_results.json"
    dump_json(results, str(results_file))
    logger.info(f"Analysis complete. Results saved to {results_file}")

    return results
