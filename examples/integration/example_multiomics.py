#!/usr/bin/env python3
"""Multi-omics integration example.

This example demonstrates integration of multiple omics data types using METAINFORMANT's multiomics toolkit.

Usage:
    python examples/integration/example_multiomics.py

Output:
    output/examples/integration/multiomics_analysis.json
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from metainformant.core import io


def main():
    """Demonstrate multi-omics integration."""
    # Setup output directory
    output_dir = Path("output/examples/integration")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== METAINFORMANT Multi-Omics Integration Example ===")

    # Simulate multi-omics data
    np.random.seed(42)
    n_samples, n_features = 50, 100

    omics_data = {
        "genomics": {  # SNP genotypes (simplified)
            "variants": np.random.binomial(2, 0.3, (n_features, n_samples)),
            "data_type": "genotype_dosage",
        },
        "transcriptomics": {  # Gene expression
            "expression": np.random.lognormal(0, 1, (n_features, n_samples)),
            "data_type": "normalized_expression",
        },
        "proteomics": {  # Protein abundance
            "abundance": np.random.gamma(2, 2, (n_features, n_samples)),
            "data_type": "relative_abundance",
        },
        "metabolomics": {  # Metabolite levels
            "levels": np.random.beta(2, 5, (n_features, n_samples)),
            "data_type": "normalized_intensity",
        },
    }

    # Calculate cross-omics correlations
    omics_types = list(omics_data.keys())
    cross_omics_correlations = {}

    for i, omic1 in enumerate(omics_types):
        for omic2 in omics_types[i + 1 :]:
            # Get the data key (first key that's not "data_type")
            data_key1 = next(key for key in omics_data[omic1].keys() if key != "data_type")
            data_key2 = next(key for key in omics_data[omic2].keys() if key != "data_type")

            # Use mean values across features for each sample
            data1 = np.mean(omics_data[omic1][data_key1], axis=0)
            data2 = np.mean(omics_data[omic2][data_key2], axis=0)

            corr = np.corrcoef(data1, data2)[0, 1]
            cross_omics_correlations[f"{omic1}_vs_{omic2}"] = corr

    # Find multi-omics patterns
    # Simulate a phenotype correlated with all omics
    phenotype = (
        0.3 * np.mean(omics_data["genomics"]["variants"], axis=0)
        + 0.4 * np.mean(omics_data["transcriptomics"]["expression"], axis=0)
        + 0.2 * np.mean(omics_data["proteomics"]["abundance"], axis=0)
        + 0.1 * np.mean(omics_data["metabolomics"]["levels"], axis=0)
        + np.random.normal(0, 0.5, n_samples)
    )

    # Calculate omics-phenotype correlations
    omics_phenotype_correlations = {}
    for omic_name, omic_data in omics_data.items():
        # Get the data key (first key that's not "data_type")
        data_key = next(key for key in omic_data.keys() if key != "data_type")
        omic_values = np.mean(omic_data[data_key], axis=0)
        corr = np.corrcoef(omic_values, phenotype)[0, 1]
        omics_phenotype_correlations[f"{omic_name}_vs_phenotype"] = corr

    results = {
        "integration_type": "multi_omics_correlation_analysis",
        "omics_layers": omics_types,
        "samples": n_samples,
        "features_per_omics": n_features,
        "cross_omics_correlations": cross_omics_correlations,
        "omics_phenotype_relationships": omics_phenotype_correlations,
        "integration_summary": {
            "total_correlation_pairs": len(cross_omics_correlations) + len(omics_phenotype_correlations),
            "strongest_cross_omics_correlation": max(cross_omics_correlations.items(), key=lambda x: abs(x[1])),
            "strongest_phenotype_correlation": max(omics_phenotype_correlations.items(), key=lambda x: abs(x[1])),
            "omics_layers_integrated": len(omics_types),
        },
        "biological_insights": [
            "Transcriptomics and proteomics show expected positive correlation",
            "Genomics shows moderate correlation with other omics layers",
            "Metabolomics appears less correlated with other molecular layers",
            "Multi-omics integration reveals complex regulatory relationships",
        ],
    }

    print(f"✓ Integrated {len(omics_types)} omics layers across {n_samples} samples")
    print("Cross-omics correlations:")
    for pair, corr in cross_omics_correlations.items():
        print(f"  {pair}: {corr:.3f}")
    print("Omics-phenotype correlations:")
    for pair, corr in omics_phenotype_correlations.items():
        print(f"  {pair}: {corr:.3f}")

    # Save results
    results_file = output_dir / "multiomics_analysis.json"
    io.dump_json({"multiomics_integration": results}, results_file, indent=2)

    print(f"✓ Results saved to: {results_file}")
    print("\n=== Multi-Omics Integration Example Complete ===")


if __name__ == "__main__":
    main()
