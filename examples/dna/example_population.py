#!/usr/bin/env python3
"""DNA population genetics example.

This example demonstrates population genetics statistics using METAINFORMANT's population genetics toolkit, including nucleotide diversity, Tajima's D, and other neutrality tests.

Usage:
    python examples/dna/example_population.py

Output:
    output/examples/dna/population_stats.json
"""

from __future__ import annotations

from pathlib import Path
from metainformant.core import io
from metainformant.dna import population

def main():
    """Demonstrate DNA population genetics analysis."""
    # Setup output directory
    output_dir = Path("output/examples/dna")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== METAINFORMANT DNA Population Genetics Example ===")

    # 1. Create simulated population data
    print("\n1. Creating simulated population data...")

    # Population 1: Low diversity (mostly similar sequences)
    population_1 = [
        "ATCGATCGATCGATCGATCG",  # Reference
        "ATCGATCGATCGATCGATCG",  # Identical
        "ATCGATCGATCGATCGATCA",  # One SNP
        "ATCGATCGATCGATCGATCG",  # Identical
        "ATCGATCGATCGATCGATCG",  # Identical
    ]

    # Population 2: High diversity (many differences)
    population_2 = [
        "ATCGATCGATCGATCGATCG",  # Reference
        "GTCGATCGATCGATCGATCG",  # SNP at position 1
        "ATCGATCGATCGATCGTTTG",  # SNPs at positions 15-19
        "ATCGATCGATCGATCGAAAA",  # SNPs at positions 15-19
        "ATCGTTTTTTTTTTTTTTTT",  # Many SNPs
    ]

    # Combined dataset for neutrality tests
    combined_population = [
        "ATCGATCGATCGATCGATCG",
        "ATCGATCGATCGATCGATCA",
        "ATCGATCGATCGATCGATCG",
        "ATCGATCGATCGATCGTTTG",
        "ATCGATCGATCGATCGATCG",
        "GTCGATCGATCGATCGATCG",
        "ATCGATCGATCGATCGAAAA",
        "ATCGATCGATCGATCGATCG",
        "ATCGTTTTTTTTTTTTTTTT",
        "ATCGATCGATCGATCGATCG",
    ]

    populations = {
        "low_diversity": population_1,
        "high_diversity": population_2,
        "combined_neutrality_test": combined_population
    }

    print(f"✓ Created {len(populations)} populations for analysis")
    print(f"  Low diversity: {len(population_1)} sequences")
    print(f"  High diversity: {len(population_2)} sequences")
    print(f"  Combined: {len(combined_population)} sequences")

    # 2. Basic diversity statistics
    print("\n2. Calculating basic diversity statistics...")

    diversity_stats = {}

    for pop_name, sequences in populations.items():
        if pop_name != "combined_neutrality_test":  # Skip combined for basic stats
            pi = population.nucleotide_diversity(sequences)
            theta = population.wattersons_theta(sequences)
            s = population.segregating_sites(sequences)

            diversity_stats[pop_name] = {
                "nucleotide_diversity_pi": pi,
                "wattersons_theta": theta,
                "segregating_sites": s,
                "pi_theta_ratio": pi / theta if theta > 0 else 0
            }

            print(f"  {pop_name}: π={pi:.4f}, θ={theta:.4f}, S={s}")

    # 3. Neutrality tests
    print("\n3. Performing neutrality tests...")

    neutrality_tests = {}

    # Tajima's D
    tajima_d = population.tajimas_d(combined_population)
    neutrality_tests["tajimas_d"] = {
        "statistic": tajima_d,
        "interpretation": "Negative: purifying selection or population expansion; Positive: balancing selection or bottleneck",
        "significance": "Typically significant if |D| > 2"
    }

    # Fu and Li's D*
    fu_li_d_star = population.fu_and_li_d_star_from_sequences(combined_population)
    neutrality_tests["fu_li_d_star"] = {
        "statistic": fu_li_d_star,
        "interpretation": "Tests for recent population growth or selection",
        "significance": "Negative values suggest population expansion"
    }

    # Fu and Li's F*
    fu_li_f_star = population.fu_and_li_f_star_from_sequences(combined_population)
    neutrality_tests["fu_li_f_star"] = {
        "statistic": fu_li_f_star,
        "interpretation": "More powerful for detecting population growth",
        "significance": "Negative values suggest population expansion"
    }

    # Fay and Wu's H
    fay_wu_h = population.fay_wu_h_from_sequences(combined_population)
    neutrality_tests["fay_wu_h"] = {
        "statistic": fay_wu_h,
        "interpretation": "Detects recent positive selection",
        "significance": "Negative values suggest recent selective sweeps"
    }

    print("  Neutrality test results:")
    for test_name, result in neutrality_tests.items():
        stat = result["statistic"]
        interpretation = "positive selection" if stat > 0 else "negative selection/growth" if stat < 0 else "neutral evolution"
        print(f"    {test_name}: {stat:.4f} ({interpretation})")

    # 4. Population differentiation (Fst)
    print("\n4. Calculating population differentiation (Fst)...")

    fst_result = population.hudson_fst(population_1, population_2)
    fst_analysis = {
        "fst_value": fst_result,
        "interpretation": {
            "0.0-0.05": "Little genetic differentiation",
            "0.05-0.15": "Moderate genetic differentiation",
            "0.15-0.25": "Great genetic differentiation",
            "0.25+": "Very great genetic differentiation"
        },
        "population_sizes": {
            "population_1": len(population_1),
            "population_2": len(population_2)
        }
    }

    if fst_result < 0.05:
        differentiation_level = "Little differentiation"
    elif fst_result < 0.15:
        differentiation_level = "Moderate differentiation"
    elif fst_result < 0.25:
        differentiation_level = "Great differentiation"
    else:
        differentiation_level = "Very great differentiation"

    print(f"  Fst between populations: {fst_result:.4f} ({differentiation_level})")

    # 5. Genotype matrix analysis
    print("\n5. Analyzing genotype matrices...")

    # Create simple genotype matrices (biallelic SNPs)
    # 0 = homozygous reference, 1 = heterozygous/alternative

    genotype_matrix_1 = [
        [0, 0, 1, 0, 0],  # SNP 1
        [0, 0, 0, 1, 0],  # SNP 2
        [1, 0, 0, 0, 0],  # SNP 3
    ]

    genotype_matrix_2 = [
        [1, 1, 1, 1, 1],  # SNP 1 - different alleles (all heterozygous)
        [1, 1, 1, 0, 1],  # SNP 2
        [0, 1, 1, 1, 1],  # SNP 3
    ]

    genotype_matrices = {
        "population_1": genotype_matrix_1,
        "population_2": genotype_matrix_2
    }

    genotype_analysis = {}

    for pop_name, matrix in genotype_matrices.items():
        allele_freqs = population.allele_frequencies(matrix)

        # Calculate heterozygosity per SNP (biallelic: 0=homo ref, 1=hetero, 2=homo alt)
        het_per_snp = []
        for snp_genotypes in matrix:
            # For biallelic SNPs, heterozygosity is count of 1's divided by total
            het_count = sum(1 for g in snp_genotypes if g == 1)
            het_per_snp.append(het_count / len(snp_genotypes) if snp_genotypes else 0)

        avg_heterozygosity = sum(het_per_snp) / len(het_per_snp) if het_per_snp else 0

        genotype_analysis[pop_name] = {
            "allele_frequencies": allele_freqs,
            "observed_heterozygosity": avg_heterozygosity,
            "heterozygosity_per_snp": het_per_snp,
            "num_snps": len(matrix),
            "num_individuals": len(matrix[0]) if matrix else 0
        }

        print(f"  {pop_name}: {len(allele_freqs)} SNPs, heterozygosity = {avg_heterozygosity:.4f}")

    # 6. Comparative analysis
    print("\n6. Comparative population analysis...")

    comparison = {
        "diversity_comparison": {
            "low_vs_high_diversity": {
                "pi_ratio": diversity_stats["high_diversity"]["nucleotide_diversity_pi"] / diversity_stats["low_diversity"]["nucleotide_diversity_pi"],
                "theta_ratio": diversity_stats["high_diversity"]["wattersons_theta"] / diversity_stats["low_diversity"]["wattersons_theta"],
                "segregating_sites_ratio": diversity_stats["high_diversity"]["segregating_sites"] / diversity_stats["low_diversity"]["segregating_sites"]
            }
        },
        "neutrality_test_summary": {
            "tests_performed": len(neutrality_tests),
            "significant_deviations": sum(1 for test in neutrality_tests.values() if abs(test["statistic"]) > 2),
            "directional_bias": {
                "negative": sum(1 for test in neutrality_tests.values() if test["statistic"] < -2),
                "positive": sum(1 for test in neutrality_tests.values() if test["statistic"] > 2)
            }
        }
    }

    print("  Diversity ratios (high/low):")
    for metric, ratio in comparison["diversity_comparison"]["low_vs_high_diversity"].items():
        print(f"    {metric}: {ratio:.2f}x")

    # 7. Create comprehensive population genetics results
    print("\n7. Creating comprehensive population genetics results...")

    population_results = {
        "dna_population_genetics_analysis": {
            "timestamp": "2024-12-26T10:00:00Z",
            "population_data": populations,
            "analysis_methods": {
                "nucleotide_diversity": {
                    "description": "Average pairwise differences between sequences (π)",
                    "formula": "π = (sum of pairwise differences) / (n choose 2)",
                    "interpretation": "Higher values indicate more genetic diversity"
                },
                "wattersons_theta": {
                    "description": "Expected nucleotide diversity under neutral model",
                    "formula": "θ = S / sum(1/i for i in 1 to n-1)",
                    "interpretation": "Expected diversity based on segregating sites"
                },
                "tajimas_d": {
                    "description": "Test for departure from neutral evolution",
                    "formula": "D = (π - θ) / sqrt(Var(π - θ))",
                    "interpretation": "Negative: purifying selection; Positive: balancing selection"
                },
                "fst": {
                    "description": "Fixation index measuring population differentiation",
                    "formula": "Fst = (Ht - Hs) / Ht",
                    "interpretation": "0: no differentiation; 1: complete differentiation"
                }
            },
            "results": {
                "diversity_statistics": diversity_stats,
                "neutrality_tests": neutrality_tests,
                "population_differentiation": fst_analysis,
                "genotype_analysis": genotype_analysis,
                "comparative_analysis": comparison
            },
            "summary_metrics": {
                "populations_analyzed": len(populations),
                "total_sequences": sum(len(seqs) for seqs in populations.values()),
                "neutrality_tests_performed": len(neutrality_tests),
                "fst_calculated": True,
                "genotype_matrices_analyzed": len(genotype_matrices),
                "diversity_range": {
                    "min_pi": min(stats["nucleotide_diversity_pi"] for stats in diversity_stats.values()),
                    "max_pi": max(stats["nucleotide_diversity_pi"] for stats in diversity_stats.values())
                }
            },
            "functions_demonstrated": [
                "nucleotide_diversity", "wattersons_theta", "tajimas_d",
                "hudson_fst", "fu_and_li_d_star_from_sequences",
                "fu_and_li_f_star_from_sequences", "fay_wu_h_from_sequences",
                "segregating_sites", "allele_frequencies", "observed_heterozygosity"
            ],
            "key_insights": [
                "High nucleotide diversity indicates recent population expansion or balancing selection",
                "Tajima's D < 0 suggests purifying selection or population growth",
                "Fst values near 0 indicate well-mixed populations",
                "Heterozygosity measures genetic variation within populations",
                "Multiple neutrality tests provide robust assessment of evolutionary forces"
            ]
        }
    }

    results_file = output_dir / "population_stats.json"
    io.dump_json(population_results, results_file, indent=2)

    print(f"✓ Comprehensive population genetics analysis saved to: {results_file}")

    print("\n=== DNA Population Genetics Example Complete ===")
    print("This example demonstrated METAINFORMANT's population genetics capabilities:")
    print("- Nucleotide diversity and Waterson's theta calculations")
    print("- Neutrality tests (Tajima's D, Fu & Li's D*/F*, Fay & Wu's H)")
    print("- Population differentiation (Fst) analysis")
    print("- Genotype matrix analysis and heterozygosity calculations")
    print("- Comparative analysis of population genetic parameters")

    print(f"\nAll outputs saved to: {output_dir}")

if __name__ == "__main__":
    main()
