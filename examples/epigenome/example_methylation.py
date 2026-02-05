#!/usr/bin/env python3
"""Epigenome methylation analysis example.

This example demonstrates DNA methylation analysis using METAINFORMANT's epigenome toolkit.

Usage:
    python examples/epigenome/example_methylation.py

Output:
    output/examples/epigenome/methylation_analysis.json
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from metainformant.core import io


def main():
    """Demonstrate methylation analysis."""
    # Setup output directory
    output_dir = Path("output/examples/epigenome")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== METAINFORMANT Epigenome Example ===")

    # Simulate methylation data
    np.random.seed(42)
    n_sites, n_samples = 1000, 20

    # Generate methylation levels (beta values 0-1)
    methylation_matrix = np.random.beta(2, 5, (n_sites, n_samples))

    # Calculate methylation statistics
    mean_methylation = np.mean(methylation_matrix, axis=1)
    hyper_methylated = np.sum(methylation_matrix > 0.8, axis=1)  # Sites with high methylation
    hypo_methylated = np.sum(methylation_matrix < 0.2, axis=1)  # Sites with low methylation

    # Differential methylation (compare first 10 vs last 10 samples)
    group1_mean = np.mean(methylation_matrix[:, :10], axis=1)
    group2_mean = np.mean(methylation_matrix[:, 10:], axis=1)
    methylation_diff = group2_mean - group1_mean

    # Find differentially methylated sites (|diff| > 0.2)
    dm_sites = np.abs(methylation_diff) > 0.2
    n_dm_sites = np.sum(dm_sites)

    results = {
        "sites_analyzed": n_sites,
        "samples_analyzed": n_samples,
        "global_methylation": {
            "mean": float(np.mean(methylation_matrix)),
            "median": float(np.median(methylation_matrix)),
            "range": [float(np.min(methylation_matrix)), float(np.max(methylation_matrix))],
        },
        "site_classification": {
            "hyper_methylated_sites": int(np.sum(np.mean(methylation_matrix, axis=1) > 0.8)),
            "hypo_methylated_sites": int(np.sum(np.mean(methylation_matrix, axis=1) < 0.2)),
            "intermediate_sites": int(
                np.sum((np.mean(methylation_matrix, axis=1) >= 0.2) & (np.mean(methylation_matrix, axis=1) <= 0.8))
            ),
        },
        "differential_methylation": {
            "comparison": "samples_1-10_vs_11-20",
            "dm_sites_found": int(n_dm_sites),
            "dm_percentage": float(n_dm_sites / n_sites * 100),
            "max_difference": float(np.max(np.abs(methylation_diff))),
        },
    }

    print(f"✓ Analyzed {n_sites} CpG sites across {n_samples} samples")
    print(
        f"Methylation range: {results['global_methylation']['range'][0]:.3f} - {results['global_methylation']['range'][1]:.3f}"
    )
    print(
        f"Differentially methylated sites: {results['differential_methylation']['dm_sites_found']} ({results['differential_methylation']['dm_percentage']:.1f}%)"
    )

    # Save results
    results_file = output_dir / "methylation_analysis.json"
    io.dump_json({"epigenome_analysis": results}, results_file, indent=2)

    print(f"✓ Results saved to: {results_file}")
    print("\n=== Epigenome Example Complete ===")


if __name__ == "__main__":
    main()
