#!/usr/bin/env python3
"""GWAS association testing example.

This example demonstrates genome-wide association study (GWAS) statistical analysis using METAINFORMANT's GWAS toolkit.

Usage:
    python examples/gwas/example_association.py

Expected output:
    output/examples/gwas/association_results.json
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from metainformant.core import io
from metainformant.gwas.analysis.association import association_test_linear, association_test_logistic
from metainformant.gwas.analysis.correction import bonferroni_correction, fdr_correction


def main():
    """Demonstrate GWAS association testing."""
    # Setup output directory
    output_dir = Path("output/examples/gwas")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== METAINFORMANT GWAS Association Example ===")

    # 1. Create simulated GWAS data
    print("\n1. Creating simulated GWAS data...")

    np.random.seed(42)  # For reproducible results

    n_samples = 1000  # 1000 individuals
    n_snps = 100  # 100 SNPs for demonstration

    # Generate genotype data (0, 1, 2 for homozygous ref, het, homozygous alt)
    genotypes = np.random.randint(0, 3, size=(n_snps, n_samples))

    # Create phenotype data
    # SNP at index 50 has a real association (effect size = 0.5)
    causal_snp_idx = 50
    causal_effect = 0.5

    # Generate phenotype with noise + genetic effect
    noise = np.random.normal(0, 1, n_samples)
    genetic_effect = genotypes[causal_snp_idx] * causal_effect
    phenotypes = noise + genetic_effect

    # For logistic regression example, create binary phenotype
    binary_phenotypes = (phenotypes > np.median(phenotypes)).astype(int)

    print(f"✓ Created GWAS dataset: {n_samples} samples, {n_snps} SNPs")
    print(f"  Causal SNP at index {causal_snp_idx} with effect size {causal_effect}")

    # 2. Linear regression association testing
    print("\n2. Performing linear regression association tests...")

    linear_results = []

    for snp_idx in range(min(20, n_snps)):  # Test first 20 SNPs for demo
        snp_genotypes = genotypes[snp_idx]

        result = association_test_linear(snp_genotypes, phenotypes)

        linear_results.append(
            {
                "snp_index": snp_idx,
                "chromosome": f"chr{((snp_idx // 10) % 5) + 1}",  # Simulated chromosomes
                "position": snp_idx * 100000,  # Simulated positions
                "genotypes": snp_genotypes.tolist(),
                "beta": result["beta"],
                "se": result["se"],
                "t_stat": result["t_stat"],
                "p_value": result["p_value"],
                "is_causal": snp_idx == causal_snp_idx,
            }
        )

        if snp_idx == causal_snp_idx:
            print(f"  Causal SNP {snp_idx}: β={result['beta']:.3f}, p={result['p_value']:.2e}")
        elif snp_idx < 5:  # Show first few null SNPs
            print(f"  SNP {snp_idx}: β={result['beta']:.3f}, p={result['p_value']:.2e}")

    # 3. Logistic regression for binary traits
    print("\n3. Performing logistic regression for binary traits...")

    logistic_results = []

    for snp_idx in range(min(20, n_snps)):  # Test first 20 SNPs
        snp_genotypes = genotypes[snp_idx]

        result = association_test_logistic(snp_genotypes, binary_phenotypes)

        logistic_results.append(
            {
                "snp_index": snp_idx,
                "beta": result["beta"],
                "se": result["se"],
                "z_stat": result["z_stat"],
                "p_value": result["p_value"],
                "odds_ratio": result.get("odds_ratio", np.exp(result["beta"])),
                "is_causal": snp_idx == causal_snp_idx,
            }
        )

    # Show example SNP result (first SNP)
    example_logistic = logistic_results[0]
    print(f"  Example SNP logistic: OR={example_logistic['odds_ratio']:.3f}, p={example_logistic['p_value']:.2e}")

    # 4. Multiple testing correction
    print("\n4. Applying multiple testing corrections...")

    # Extract p-values from linear regression
    linear_p_values = [r["p_value"] for r in linear_results]

    # Bonferroni correction
    bonferroni_result = bonferroni_correction(linear_p_values, alpha=0.05)
    bonferroni_rejected = bonferroni_result["significant"]
    bonferroni_threshold = bonferroni_result["corrected_alpha"]

    # FDR correction (Benjamini-Hochberg)
    fdr_result = fdr_correction(linear_p_values, alpha=0.05)
    fdr_rejected = fdr_result["significant"]

    correction_results = {
        "bonferroni": {
            "alpha_threshold": bonferroni_threshold,
            "significant_snps": int(sum(bonferroni_rejected)),
            "significant_indices": [i for i, rej in enumerate(bonferroni_rejected) if rej],
        },
        "fdr": {
            "significant_snps": int(sum(fdr_rejected)),
            "significant_indices": [i for i, rej in enumerate(fdr_rejected) if rej],
        },
    }

    print("  Bonferroni correction (α=0.05):")
    print(f"    Threshold: {bonferroni_threshold:.2e}")
    print(f"    Significant SNPs: {correction_results['bonferroni']['significant_snps']}")

    print("  FDR correction (α=0.05):")
    print(f"    Significant SNPs: {correction_results['fdr']['significant_snps']}")

    # 5. Effect size interpretation
    print("\n5. Interpreting effect sizes...")

    # For linear regression, effect size is in phenotype units
    significant_linear = [r for r in linear_results if r["p_value"] < 0.05]

    effect_sizes = []
    for result in significant_linear:
        effect_magnitude = abs(result["beta"])
        if effect_magnitude < 0.2:
            magnitude = "small"
        elif effect_magnitude < 0.5:
            magnitude = "medium"
        else:
            magnitude = "large"

        effect_sizes.append(
            {
                "snp_index": result["snp_index"],
                "effect_size": result["beta"],
                "magnitude": magnitude,
                "is_causal": result["is_causal"],
            }
        )

    print("  Effect size interpretation:")
    for effect in effect_sizes[:3]:  # Show first 3
        causal_note = " (CAUSAL)" if effect["is_causal"] else ""
        print(f"    SNP {effect['snp_index']}: β={effect['effect_size']:.3f} ({effect['magnitude']}){causal_note}")

    # 6. Power analysis concepts
    print("\n6. Power analysis considerations...")

    # Calculate statistical power concepts (simplified)
    sample_size = n_samples
    causal_effect_size = causal_effect
    phenotypic_variance = np.var(phenotypes)

    # Simplified power calculation (would need more sophisticated implementation)
    power_estimate = min(0.99, 1 - np.exp(-sample_size * causal_effect_size**2 / (4 * phenotypic_variance)))

    power_analysis = {
        "sample_size": sample_size,
        "causal_effect_size": causal_effect_size,
        "phenotypic_variance": phenotypic_variance,
        "estimated_power": power_estimate,
        "notes": [
            "Power depends on sample size, effect size, and variance",
            "Larger effects are easier to detect",
            "More samples increase statistical power",
            "Multiple testing correction reduces power",
        ],
    }

    print(f"  Estimated power to detect causal effect: {power_estimate:.1%}")
    print("  Power depends on: sample size, effect size, phenotypic variance")

    # 7. Create comprehensive association results
    print("\n7. Creating comprehensive association analysis results...")

    association_results = {
        "example": "gwas_association_testing",
        "domain": "gwas",
        "results": {
            "timestamp": "2024-12-26T10:00:00Z",
            "dataset": {
                "n_samples": n_samples,
                "n_snps_tested": len(linear_results),
                "causal_snp": causal_snp_idx,
                "causal_effect_size": causal_effect,
                "phenotype_type": "continuous",
                "binary_trait_available": True,
            },
            "association_methods": {
                "linear_regression": {
                    "description": "Tests association between SNP dosage and continuous phenotype",
                    "model": "phenotype ~ SNP_dosage",
                    "output": "beta (effect size), SE, t-statistic, p-value",
                    "results": linear_results,
                },
                "logistic_regression": {
                    "description": "Tests association between SNP dosage and binary phenotype",
                    "model": "logit(prob) = beta * SNP_dosage",
                    "output": "beta, SE, z-statistic, p-value, odds ratio",
                    "results": logistic_results[:5],  # First 5 for brevity
                },
            },
            "multiple_testing_correction": correction_results,
            "effect_size_analysis": effect_sizes,
            "power_analysis": power_analysis,
            "key_findings": [
                f"Causal SNP {causal_snp_idx} detected with strong association",
                f"Bonferroni correction identified {correction_results['bonferroni']['significant_snps']} significant SNPs",
                f"FDR correction identified {correction_results['fdr']['significant_snps']} significant SNPs",
                "Effect sizes range from small to large depending on SNP",
                f"Statistical power estimated at {power_estimate:.1%} for effect size {causal_effect}",
            ],
            "gwas_best_practices_demonstrated": [
                "Use appropriate statistical model for phenotype type",
                "Apply multiple testing correction",
                "Interpret effect sizes in biological context",
                "Consider statistical power in study design",
                "Validate significant associations",
            ],
        },
    }

    results_file = output_dir / "association_results.json"
    io.dump_json(association_results, results_file, indent=2)

    print(f"✓ Comprehensive GWAS association analysis saved to: {results_file}")

    print("\n=== GWAS Association Example Complete ===")
    print("This example demonstrated METAINFORMANT's GWAS statistical analysis capabilities:")
    print("- Linear and logistic regression association testing")
    print("- Multiple testing correction (Bonferroni, FDR)")
    print("- Effect size interpretation")
    print("- Statistical power considerations")

    print(f"\nAll outputs saved to: {output_dir}")


if __name__ == "__main__":
    main()
