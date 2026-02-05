#!/usr/bin/env python3
"""Visualization example.

This example demonstrates plotting capabilities using METAINFORMANT's visualization toolkit.

Usage:
    python examples/visualization/example_plots.py

Output:
    output/examples/visualization/plots_demo.json
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from metainformant.core import io
from metainformant.visualization.plots import lineplot


def main():
    """Demonstrate visualization capabilities."""
    # Setup output directory
    output_dir = Path("output/examples/visualization")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== METAINFORMANT Visualization Example ===")

    # Generate sample data
    x = np.linspace(0, 10, 50)
    y1 = np.sin(x) + np.random.normal(0, 0.1, len(x))
    y2 = np.cos(x) + np.random.normal(0, 0.1, len(x))

    # Create line plot
    try:
        ax = lineplot(x, y1, label="Signal 1")
        ax.plot(x, y2, label="Signal 2", linestyle="--")
        ax.set_xlabel("Time")
        ax.set_ylabel("Amplitude")
        ax.set_title("Biological Signal Comparison")
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Save plot
        plot_path = output_dir / "biological_signals.png"
        ax.figure.savefig(plot_path, dpi=300, bbox_inches="tight")

        print(f"✓ Created visualization: {plot_path}")

        results = {
            "plot_type": "line_plot",
            "data_points": len(x),
            "signals": ["signal_1", "signal_2"],
            "output_format": "PNG",
            "resolution": "300 DPI",
        }

    except Exception as e:
        print(f"Plotting failed: {e}")
        results = {"error": str(e)}

    # Save results
    results_file = output_dir / "plots_demo.json"
    io.dump_json({"visualization_demo": results}, results_file, indent=2)

    print(f"✓ Results saved to: {results_file}")
    print("\n=== Visualization Example Complete ===")


if __name__ == "__main__":
    main()
