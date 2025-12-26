#!/usr/bin/env python3
"""Ecology community analysis example.

This example demonstrates ecological community analysis using METAINFORMANT's ecology toolkit.

Usage:
    python examples/ecology/example_community.py

Output:
    output/examples/ecology/community_analysis.json
"""

from __future__ import annotations

import numpy as np
from pathlib import Path
from metainformant.core import io
from metainformant.ecology.community import community_metrics

def main():
    """Demonstrate ecological community analysis."""
    # Setup output directory
    output_dir = Path("output/examples/ecology")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== METAINFORMANT Ecology Example ===")

    # Create simulated community data
    species = ["Species_A", "Species_B", "Species_C", "Species_D", "Species_E"]
    n_sites = 5

    # Abundance matrix (species x sites)
    abundance_matrix = np.random.poisson(lam=[20, 15, 10, 8, 5], size=(len(species), n_sites))

    print(f"✓ Created community data: {len(species)} species across {n_sites} sites")

    # Calculate diversity for each site
    diversity_indices = []
    for site_abundances in abundance_matrix.T:  # Each column is a site
        metrics = community_metrics(site_abundances)
        diversity_indices.append(metrics["shannon"])

    results = {
        "species": species,
        "sites": n_sites,
        "abundance_matrix": abundance_matrix.tolist(),
        "diversity_indices": diversity_indices,
        "mean_diversity": float(np.mean(diversity_indices))
    }

    print(f"Mean Shannon diversity: {results['mean_diversity']:.3f}")

    # Save results
    results_file = output_dir / "community_analysis.json"
    io.dump_json({
        "ecology_analysis": results
    }, results_file, indent=2)

    print(f"✓ Results saved to: {results_file}")
    print("\n=== Ecology Example Complete ===")

if __name__ == "__main__":
    main()
