#!/usr/bin/env python3
"""Visualization data simulation script.

This script generates synthetic data for various plot types including
time-series, multi-dimensional, and statistical test data.

Usage:
    python3 scripts/simulation/simulate_visualization.py --type timeseries --n-points 100
    python3 scripts/simulation/simulate_visualization.py --type multidim --n-samples 50 --n-features 10
    python3 scripts/simulation/simulate_visualization.py --type statistical --n-samples 100
"""

import argparse
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core import io, logging, paths, validation

logger = logging.get_logger(__name__)


def simulate_timeseries(
    output_dir: Path,
    n_points: int,
    n_series: int,
    noise_level: float,
    seed: int,
) -> dict:
    """Simulate time-series data."""
    logger.info(f"Generating time-series data: {n_points} points, {n_series} series")
    rng = random.Random(seed)
    np.random.seed(seed)
    
    time_points = np.linspace(0, 10, n_points)
    data = {"time": time_points}
    
    for i in range(n_series):
        # Generate different patterns
        if i % 3 == 0:
            # Sine wave
            values = np.sin(2 * np.pi * time_points) + np.random.normal(0, noise_level, n_points)
        elif i % 3 == 1:
            # Exponential growth
            values = np.exp(0.2 * time_points) + np.random.normal(0, noise_level, n_points)
        else:
            # Linear trend
            values = 0.5 * time_points + np.random.normal(0, noise_level, n_points)
        
        data[f"series_{i:03d}"] = values
    
    df = pd.DataFrame(data)
    csv_file = output_dir / "timeseries_data.csv"
    df.to_csv(csv_file, index=False)
    
    logger.info(f"Time-series data saved to {csv_file}")
    
    return {
        "type": "timeseries",
        "n_points": n_points,
        "n_series": n_series,
        "output_file": str(csv_file),
    }


def simulate_multidim(
    output_dir: Path,
    n_samples: int,
    n_features: int,
    n_clusters: int,
    seed: int,
) -> dict:
    """Simulate multi-dimensional data with clusters."""
    logger.info(f"Generating multi-dimensional data: {n_samples} samples, {n_features} features")
    rng = random.Random(seed)
    np.random.seed(seed)
    
    # Generate cluster centers
    cluster_centers = np.random.randn(n_clusters, n_features) * 5
    
    # Generate data points around centers
    data = []
    labels = []
    
    samples_per_cluster = n_samples // n_clusters
    
    for cluster_idx in range(n_clusters):
        center = cluster_centers[cluster_idx]
        for _ in range(samples_per_cluster):
            point = center + np.random.randn(n_features)
            data.append(point)
            labels.append(cluster_idx)
    
    # Add remaining samples
    for _ in range(n_samples - len(data)):
        cluster_idx = rng.randint(0, n_clusters - 1)
        center = cluster_centers[cluster_idx]
        point = center + np.random.randn(n_features)
        data.append(point)
        labels.append(cluster_idx)
    
    # Create DataFrame
    feature_names = [f"feature_{i:03d}" for i in range(n_features)]
    df = pd.DataFrame(data, columns=feature_names)
    df.insert(0, "sample_id", [f"sample_{i:04d}" for i in range(n_samples)])
    df["cluster"] = labels
    
    csv_file = output_dir / "multidim_data.csv"
    df.to_csv(csv_file, index=False)
    
    logger.info(f"Multi-dimensional data saved to {csv_file}")
    
    return {
        "type": "multidim",
        "n_samples": n_samples,
        "n_features": n_features,
        "n_clusters": n_clusters,
        "output_file": str(csv_file),
    }


def simulate_statistical(
    output_dir: Path,
    n_samples: int,
    n_groups: int,
    seed: int,
) -> dict:
    """Simulate statistical test data."""
    logger.info(f"Generating statistical test data: {n_samples} samples, {n_groups} groups")
    rng = random.Random(seed)
    np.random.seed(seed)
    
    samples_per_group = n_samples // n_groups
    
    data = []
    for group_idx in range(n_groups):
        # Different means per group
        mean = group_idx * 2.0
        std = 1.0
        
        for _ in range(samples_per_group):
            value = np.random.normal(mean, std)
            data.append({
                "group": f"group_{group_idx:02d}",
                "value": value,
            })
    
    # Add remaining samples
    for _ in range(n_samples - len(data)):
        group_idx = rng.randint(0, n_groups - 1)
        mean = group_idx * 2.0
        value = np.random.normal(mean, 1.0)
        data.append({
            "group": f"group_{group_idx:02d}",
            "value": value,
        })
    
    df = pd.DataFrame(data)
    csv_file = output_dir / "statistical_data.csv"
    df.to_csv(csv_file, index=False)
    
    logger.info(f"Statistical data saved to {csv_file}")
    
    return {
        "type": "statistical",
        "n_samples": n_samples,
        "n_groups": n_groups,
        "output_file": str(csv_file),
    }


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Visualization data simulation for testing and validation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Simulate time-series data
  %(prog)s --type timeseries --n-points 100 --n-series 5

  # Simulate multi-dimensional data
  %(prog)s --type multidim --n-samples 50 --n-features 10 --n-clusters 3

  # Simulate statistical test data
  %(prog)s --type statistical --n-samples 100 --n-groups 4
        """,
    )
    
    parser.add_argument(
        "--type",
        required=True,
        choices=["timeseries", "multidim", "statistical"],
        help="Type of simulation to run",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/simulation/visualization"),
        help="Output directory (default: output/simulation/visualization)",
    )
    parser.add_argument("--n-points", type=int, default=100, help="Number of time points (timeseries type)")
    parser.add_argument("--n-series", type=int, default=3, help="Number of series (timeseries type)")
    parser.add_argument("--noise-level", type=float, default=0.1, help="Noise level (timeseries type)")
    parser.add_argument("--n-samples", type=int, default=50, help="Number of samples (multidim/statistical types)")
    parser.add_argument("--n-features", type=int, default=10, help="Number of features (multidim type)")
    parser.add_argument("--n-clusters", type=int, default=3, help="Number of clusters (multidim type)")
    parser.add_argument("--n-groups", type=int, default=3, help="Number of groups (statistical type)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")
    
    args = parser.parse_args()
    
    if args.verbose:
        import logging as std_logging
        logger.setLevel(std_logging.DEBUG)
    
    output_dir = paths.ensure_directory(args.output)
    
    try:
        if args.type == "timeseries":
            results = simulate_timeseries(
                output_dir, args.n_points, args.n_series, args.noise_level, args.seed
            )
        elif args.type == "multidim":
            results = simulate_multidim(
                output_dir, args.n_samples, args.n_features, args.n_clusters, args.seed
            )
        elif args.type == "statistical":
            results = simulate_statistical(output_dir, args.n_samples, args.n_groups, args.seed)
        
        # Save summary
        summary_file = output_dir / "simulation_summary.json"
        io.dump_json(results, summary_file, indent=2)
        logger.info(f"Simulation complete. Summary saved to {summary_file}")
        
        return 0
    except Exception as e:
        logger.error(f"Simulation failed: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())

