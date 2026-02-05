#!/usr/bin/env python3
"""Biological simulation example.

This example demonstrates synthetic data generation using METAINFORMANT's simulation toolkit.

Usage:
    python examples/simulation/example_simulation.py

Output:
    output/examples/simulation/simulation_results.json
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from metainformant.core import io


def main():
    """Demonstrate biological simulation."""
    # Setup output directory
    output_dir = Path("output/examples/simulation")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== METAINFORMANT Simulation Example ===")

    # Simulate gene expression evolution
    np.random.seed(42)
    n_generations, n_genes = 10, 50

    # Initial expression levels
    expression = np.random.lognormal(0, 1, n_genes)

    # Simulate evolution over generations
    expression_trajectory = [expression.copy()]

    for gen in range(1, n_generations):
        # Add random variation (mutation) and selection pressure
        mutation = np.random.normal(0, 0.1, n_genes)
        selection = -0.05 * (expression - np.mean(expression))  # Favor average expression
        expression = expression + mutation + selection
        expression = np.maximum(expression, 0.1)  # Floor at minimum expression
        expression_trajectory.append(expression.copy())

    results = {
        "simulation_type": "gene_expression_evolution",
        "generations": n_generations,
        "genes": n_genes,
        "initial_mean_expression": float(np.mean(expression_trajectory[0])),
        "final_mean_expression": float(np.mean(expression_trajectory[-1])),
        "expression_stability": float(np.std(expression_trajectory[-1]) / np.mean(expression_trajectory[-1])),
        "trajectory_summary": {
            "start": [float(x) for x in expression_trajectory[0][:5]],  # First 5 genes
            "end": [float(x) for x in expression_trajectory[-1][:5]],
        },
    }

    print(f"✓ Simulated {n_generations} generations of evolution")
    print(f"Expression stability: {results['expression_stability']:.3f}")

    # Save results
    results_file = output_dir / "simulation_results.json"
    io.dump_json({"biological_simulation": results}, results_file, indent=2)

    print(f"✓ Results saved to: {results_file}")
    print("\n=== Simulation Example Complete ===")


if __name__ == "__main__":
    main()
