#!/usr/bin/env python3
"""DNA-RNA integration example.

This example demonstrates integration of DNA sequence features with RNA expression analysis using METAINFORMANT.

Usage:
    python examples/integration/example_dna_rna.py

Output:
    output/examples/integration/dna_rna_integration.json
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from metainformant.core import io


def main():
    """Demonstrate DNA-RNA integration."""
    # Setup output directory
    output_dir = Path("output/examples/integration")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== METAINFORMANT DNA-RNA Integration Example ===")

    # Simulate DNA sequence features and RNA expression
    np.random.seed(42)
    n_genes = 100

    # DNA features
    dna_features = {
        "gene_lengths": np.random.normal(2000, 500, n_genes),
        "gc_content": np.random.beta(2, 2, n_genes),  # GC content 0-1
        "intron_count": np.random.poisson(3, n_genes),
        "tss_motifs": np.random.binomial(1, 0.3, n_genes),  # Presence of TATA box
    }

    # RNA expression (correlated with DNA features)
    expression_levels = (
        0.5 * dna_features["gc_content"] * 1000  # Higher GC -> higher expression
        + 0.2 * (dna_features["tss_motifs"] * 500)  # TATA box -> higher expression
        + np.random.lognormal(0, 0.5, n_genes) * 100  # Biological noise
    )

    # Calculate correlations between DNA features and RNA expression
    correlations = {}
    for feature_name, feature_values in dna_features.items():
        corr = np.corrcoef(feature_values, expression_levels)[0, 1]
        correlations[f"{feature_name}_vs_expression"] = corr

    # Find genes with high expression and specific DNA features
    high_expression_threshold = np.percentile(expression_levels, 75)
    high_expression_genes = expression_levels > high_expression_threshold

    dna_features_high_expr = {}
    for feature_name, feature_values in dna_features.items():
        high_expr_values = feature_values[high_expression_genes]
        dna_features_high_expr[feature_name] = {
            "mean": float(np.mean(high_expr_values)),
            "median": float(np.median(high_expr_values)),
            "vs_all_genes": float(np.mean(high_expr_values) / np.mean(feature_values)),
        }

    results = {
        "integration_type": "dna_features_to_rna_expression",
        "genes_analyzed": n_genes,
        "dna_features_integrated": list(dna_features.keys()),
        "correlations_found": correlations,
        "top_correlation": max(correlations.items(), key=lambda x: abs(x[1])),
        "high_expression_analysis": {
            "threshold": float(high_expression_threshold),
            "genes_above_threshold": int(np.sum(high_expression_genes)),
            "dna_features_in_high_expr_genes": dna_features_high_expr,
        },
        "integration_insights": [
            "GC content shows moderate positive correlation with expression",
            "TATA box presence correlates with higher expression levels",
            "High-expressing genes tend to have different DNA features",
            "Integration reveals regulatory relationships across omics layers",
        ],
    }

    print(f"✓ Integrated DNA features with RNA expression for {n_genes} genes")
    print(f"Top correlation: {results['top_correlation'][0]} = {results['top_correlation'][1]:.3f}")
    print(f"High-expressing genes: {results['high_expression_analysis']['genes_above_threshold']}")

    # Save results
    results_file = output_dir / "dna_rna_integration.json"
    io.dump_json({"dna_rna_integration": results}, results_file, indent=2)

    print(f"✓ Results saved to: {results_file}")
    print("\n=== DNA-RNA Integration Example Complete ===")


if __name__ == "__main__":
    main()
