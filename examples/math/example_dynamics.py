#!/usr/bin/env python3
"""Mathematical biology population dynamics example.

This example demonstrates population dynamics modeling using METAINFORMANT's math toolkit.

Usage:
    python examples/math/example_dynamics.py

Output:
    output/examples/math/dynamics_model.json
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from metainformant.core import io


def logistic_growth(N0, r, K, t_max):
    """Simple logistic growth model."""
    t = np.arange(t_max)
    N = K / (1 + (K - N0) / N0 * np.exp(-r * t))
    return t.tolist(), N.tolist()


def main():
    """Demonstrate population dynamics modeling."""
    # Setup output directory
    output_dir = Path("output/examples/math")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== METAINFORMANT Math Biology Example ===")

    # Logistic growth parameters
    N0 = 10  # Initial population
    r = 0.5  # Growth rate
    K = 100  # Carrying capacity
    t_max = 20  # Time steps

    t, N = logistic_growth(N0, r, K, t_max)

    results = {
        "model": "logistic_growth",
        "parameters": {"N0": N0, "r": r, "K": K, "t_max": t_max},
        "time_series": {"time": t, "population": N},
        "final_population": N[-1],
        "carrying_capacity_reached": N[-1] >= K * 0.95,
    }

    print(f"✓ Modeled population growth: {N0} → {results['final_population']:.1f}")

    # Save results
    results_file = output_dir / "dynamics_model.json"
    io.dump_json({"math_biology_modeling": results}, results_file, indent=2)

    print(f"✓ Results saved to: {results_file}")
    print("\n=== Math Biology Example Complete ===")


if __name__ == "__main__":
    main()
