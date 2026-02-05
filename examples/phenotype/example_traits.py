#!/usr/bin/env python3
"""Phenotype trait analysis example.

This example demonstrates phenotypic trait analysis using METAINFORMANT's phenotype toolkit.

Usage:
    python examples/phenotype/example_traits.py

Output:
    output/examples/phenotype/trait_analysis.json
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from metainformant.core import io


def main():
    """Demonstrate phenotype trait analysis."""
    # Setup output directory
    output_dir = Path("output/examples/phenotype")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== METAINFORMANT Phenotype Example ===")

    # Create simulated trait data
    np.random.seed(42)
    n_individuals = 100

    traits = {
        "body_length": np.random.normal(10, 2, n_individuals),
        "wing_span": np.random.normal(8, 1.5, n_individuals),
        "colony_size": np.random.normal(500, 100, n_individuals),
    }

    # Calculate correlations
    correlations = {}
    trait_names = list(traits.keys())

    for i, trait1 in enumerate(trait_names):
        for trait2 in trait_names[i + 1 :]:
            corr = np.corrcoef(traits[trait1], traits[trait2])[0, 1]
            correlations[f"{trait1}_vs_{trait2}"] = corr

    print(f"✓ Analyzed {len(traits)} traits for {n_individuals} individuals")
    print("Trait correlations:")
    for pair, corr in correlations.items():
        print(f"  {pair}: {corr:.3f}")

    # Save results
    results_file = output_dir / "trait_analysis.json"
    io.dump_json(
        {
            "phenotype_analysis": {
                "traits_analyzed": trait_names,
                "individuals": n_individuals,
                "correlations": correlations,
            }
        },
        results_file,
        indent=2,
    )

    print(f"✓ Results saved to: {results_file}")
    print("\n=== Phenotype Example Complete ===")


if __name__ == "__main__":
    main()
