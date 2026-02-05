#!/usr/bin/env python3
"""Phenotype simulation script.

This script generates synthetic phenotypic trait data including morphological
and behavioral traits with correlations.

Usage:
    python3 scripts/simulation/simulate_phenotype.py --type morphological --n-samples 100
    python3 scripts/simulation/simulate_phenotype.py --type behavioral --n-samples 200
    python3 scripts/simulation/simulate_phenotype.py --type correlation --n-samples 150 --n-traits 10
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


def simulate_morphological(
    output_dir: Path,
    n_samples: int,
    seed: int,
) -> dict:
    """Simulate morphological traits."""
    logger.info(f"Generating morphological traits: {n_samples} samples")
    rng = random.Random(seed)
    np.random.seed(seed)

    # Common morphological traits
    traits = {
        "body_length": np.random.normal(10.0, 2.0, n_samples),
        "head_width": np.random.normal(2.5, 0.5, n_samples),
        "antenna_length": np.random.normal(5.0, 1.0, n_samples),
        "leg_length": np.random.normal(8.0, 1.5, n_samples),
        "wing_length": np.random.normal(12.0, 2.5, n_samples),
    }

    # Ensure positive values
    for trait_name, values in traits.items():
        traits[trait_name] = np.maximum(values, 0.1)

    # Create DataFrame
    df = pd.DataFrame(traits)
    df.insert(0, "sample_id", [f"sample_{i:04d}" for i in range(n_samples)])

    csv_file = output_dir / "morphological_traits.csv"
    df.to_csv(csv_file, index=False)

    logger.info(f"Morphological traits saved to {csv_file}")

    return {
        "type": "morphological",
        "n_samples": n_samples,
        "n_traits": len(traits),
        "output_file": str(csv_file),
    }


def simulate_behavioral(
    output_dir: Path,
    n_samples: int,
    seed: int,
) -> dict:
    """Simulate behavioral traits."""
    logger.info(f"Generating behavioral traits: {n_samples} samples")
    rng = random.Random(seed)
    np.random.seed(seed)

    # Common behavioral traits
    traits = {
        "foraging_time": np.random.normal(120.0, 30.0, n_samples),  # minutes
        "nest_construction_time": np.random.normal(60.0, 15.0, n_samples),
        "aggression_score": np.random.normal(5.0, 1.5, n_samples),
        "social_interactions": np.random.poisson(20, n_samples),
        "exploration_distance": np.random.normal(50.0, 15.0, n_samples),  # meters
    }

    # Ensure reasonable values
    traits["foraging_time"] = np.maximum(traits["foraging_time"], 0)
    traits["nest_construction_time"] = np.maximum(traits["nest_construction_time"], 0)
    traits["aggression_score"] = np.maximum(traits["aggression_score"], 0)
    traits["social_interactions"] = np.maximum(traits["social_interactions"], 0)
    traits["exploration_distance"] = np.maximum(traits["exploration_distance"], 0)

    # Create DataFrame
    df = pd.DataFrame(traits)
    df.insert(0, "sample_id", [f"sample_{i:04d}" for i in range(n_samples)])

    csv_file = output_dir / "behavioral_traits.csv"
    df.to_csv(csv_file, index=False)

    logger.info(f"Behavioral traits saved to {csv_file}")

    return {
        "type": "behavioral",
        "n_samples": n_samples,
        "n_traits": len(traits),
        "output_file": str(csv_file),
    }


def simulate_correlation(
    output_dir: Path,
    n_samples: int,
    n_traits: int,
    correlation_strength: float,
    seed: int,
) -> dict:
    """Simulate correlated traits."""
    logger.info(f"Generating correlated traits: {n_samples} samples, {n_traits} traits")
    rng = random.Random(seed)
    np.random.seed(seed)

    # Generate correlated data using multivariate normal
    mean = np.zeros(n_traits)
    # Create correlation matrix
    corr_matrix = np.eye(n_traits)
    for i in range(n_traits):
        for j in range(i + 1, n_traits):
            if rng.random() < 0.3:  # 30% chance of correlation
                corr = rng.uniform(-correlation_strength, correlation_strength)
                corr_matrix[i, j] = corr
                corr_matrix[j, i] = corr

    # Convert to covariance matrix
    std = np.ones(n_traits)
    cov_matrix = np.outer(std, std) * corr_matrix

    # Generate data
    data = np.random.multivariate_normal(mean, cov_matrix, n_samples)

    # Create DataFrame
    trait_names = [f"trait_{i:03d}" for i in range(n_traits)]
    df = pd.DataFrame(data, columns=trait_names)
    df.insert(0, "sample_id", [f"sample_{i:04d}" for i in range(n_samples)])

    csv_file = output_dir / "correlated_traits.csv"
    df.to_csv(csv_file, index=False)

    # Save correlation matrix
    corr_df = pd.DataFrame(corr_matrix, index=trait_names, columns=trait_names)
    corr_file = output_dir / "correlation_matrix.csv"
    corr_df.to_csv(corr_file)

    logger.info(f"Correlated traits saved to {csv_file}")

    return {
        "type": "correlation",
        "n_samples": n_samples,
        "n_traits": n_traits,
        "output_file": str(csv_file),
        "correlation_file": str(corr_file),
    }


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Phenotype simulation for testing and validation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Simulate morphological traits
  %(prog)s --type morphological --n-samples 100

  # Simulate behavioral traits
  %(prog)s --type behavioral --n-samples 200

  # Simulate correlated traits
  %(prog)s --type correlation --n-samples 150 --n-traits 10
        """,
    )

    parser.add_argument(
        "--type",
        required=True,
        choices=["morphological", "behavioral", "correlation"],
        help="Type of simulation to run",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/simulation/phenotype"),
        help="Output directory (default: output/simulation/phenotype)",
    )
    parser.add_argument("--n-samples", type=int, default=100, help="Number of samples")
    parser.add_argument("--n-traits", type=int, default=10, help="Number of traits (correlation type)")
    parser.add_argument(
        "--correlation-strength", type=float, default=0.5, help="Correlation strength (correlation type)"
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")

    args = parser.parse_args()

    if args.verbose:
        import logging as std_logging

        logger.setLevel(std_logging.DEBUG)

    output_dir = paths.ensure_directory(args.output)

    try:
        if args.type == "morphological":
            results = simulate_morphological(output_dir, args.n_samples, args.seed)
        elif args.type == "behavioral":
            results = simulate_behavioral(output_dir, args.n_samples, args.seed)
        elif args.type == "correlation":
            results = simulate_correlation(
                output_dir, args.n_samples, args.n_traits, args.correlation_strength, args.seed
            )

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
