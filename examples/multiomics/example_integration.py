#!/usr/bin/env python3
"""Multi-omics integration example.

This example demonstrates multi-omic data integration using METAINFORMANT's multiomics toolkit.

Usage:
    python examples/multiomics/example_integration.py

Output:
    output/examples/multiomics/integration_analysis.json
"""

from __future__ import annotations

import numpy as np
from pathlib import Path
from metainformant.core import io

def main():
    """Demonstrate multi-omics integration."""
    # Setup output directory
    output_dir = Path("output/examples/multiomics")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== METAINFORMANT Multi-Omics Example ===")

    # Simulated multi-omics data
    np.random.seed(42)
    n_samples = 50

    omics_data = {
        "dna_methylation": np.random.beta(2, 5, n_samples),  # Methylation levels
        "rna_expression": np.random.lognormal(0, 1, n_samples),  # Expression values
        "protein_abundance": np.random.gamma(2, 2, n_samples)  # Protein levels
    }

    # Calculate correlations between omics layers
    correlations = {}
    omics_types = list(omics_data.keys())

    for i, omic1 in enumerate(omics_types):
        for omic2 in omics_types[i+1:]:
            corr = np.corrcoef(omics_data[omic1], omics_data[omic2])[0, 1]
            correlations[f"{omic1}_vs_{omic2}"] = corr

    results = {
        "omics_layers": omics_types,
        "samples": n_samples,
        "correlations": correlations,
        "integration_summary": {
            "layers_integrated": len(omics_types),
            "correlation_pairs": len(correlations),
            "strongest_correlation": max(correlations.items(), key=lambda x: abs(x[1]))
        }
    }

    print(f"✓ Integrated {len(omics_types)} omics layers")
    print("Strongest correlation:", results["integration_summary"]["strongest_correlation"])

    # Save results
    results_file = output_dir / "integration_analysis.json"
    io.dump_json({
        "multiomics_integration": results
    }, results_file, indent=2)

    print(f"✓ Results saved to: {results_file}")
    print("\n=== Multi-Omics Example Complete ===")

if __name__ == "__main__":
    main()
