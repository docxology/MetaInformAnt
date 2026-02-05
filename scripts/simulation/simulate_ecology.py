#!/usr/bin/env python3
"""Ecology simulation script.

This script generates synthetic ecological data including species abundance matrices,
community composition, and environmental metadata.

Usage:
    python3 scripts/simulation/simulate_ecology.py --type abundance --n-species 50 --n-samples 20
    python3 scripts/simulation/simulate_ecology.py --type community --n-species 100 --n-samples 30
    python3 scripts/simulation/simulate_ecology.py --type environmental --n-samples 25
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


def simulate_abundance(
    output_dir: Path,
    n_species: int,
    n_samples: int,
    seed: int,
) -> dict:
    """Simulate species abundance matrix."""
    logger.info(f"Generating abundance matrix: {n_species} species x {n_samples} samples")
    rng = random.Random(seed)
    np.random.seed(seed)

    # Generate abundance data (log-normal distribution)
    abundances = np.random.lognormal(mean=2.0, sigma=1.5, size=(n_species, n_samples))

    # Create DataFrame
    species_ids = [f"species_{i:04d}" for i in range(n_species)]
    sample_ids = [f"sample_{i:03d}" for i in range(n_samples)]
    df = pd.DataFrame(abundances, index=species_ids, columns=sample_ids)

    csv_file = output_dir / "species_abundance.csv"
    df.to_csv(csv_file)

    logger.info(f"Abundance matrix saved to {csv_file}")

    return {
        "type": "abundance",
        "n_species": n_species,
        "n_samples": n_samples,
        "output_file": str(csv_file),
    }


def simulate_community(
    output_dir: Path,
    n_species: int,
    n_samples: int,
    diversity_index: float,
    seed: int,
) -> dict:
    """Simulate community composition with specified diversity."""
    logger.info(f"Generating community composition: {n_species} species")
    rng = random.Random(seed)
    np.random.seed(seed)

    # Generate communities with controlled diversity
    communities = []
    for sample_idx in range(n_samples):
        # Generate relative abundances
        abundances = np.random.exponential(scale=1.0, size=n_species)
        abundances = abundances / abundances.sum()  # Normalize

        # Adjust for diversity
        if diversity_index < 1.0:
            # Lower diversity: concentrate abundance in fewer species
            n_dominant = max(1, int(n_species * diversity_index))
            dominant_indices = rng.sample(range(n_species), n_dominant)
            new_abundances = np.zeros(n_species)
            for idx in dominant_indices:
                new_abundances[idx] = abundances[idx]
            new_abundances = new_abundances / new_abundances.sum()
            abundances = new_abundances

        communities.append(
            {
                "sample_id": f"sample_{sample_idx:03d}",
                "abundances": abundances.tolist(),
            }
        )

    # Create abundance matrix
    abundance_matrix = np.array([c["abundances"] for c in communities]).T
    species_ids = [f"species_{i:04d}" for i in range(n_species)]
    sample_ids = [f"sample_{i:03d}" for i in range(n_samples)]
    df = pd.DataFrame(abundance_matrix, index=species_ids, columns=sample_ids)

    csv_file = output_dir / "community_abundance.csv"
    df.to_csv(csv_file)

    logger.info(f"Community composition saved to {csv_file}")

    return {
        "type": "community",
        "n_species": n_species,
        "n_samples": n_samples,
        "diversity_index": diversity_index,
        "output_file": str(csv_file),
    }


def simulate_environmental(
    output_dir: Path,
    n_samples: int,
    seed: int,
) -> dict:
    """Simulate environmental metadata."""
    logger.info(f"Generating environmental metadata: {n_samples} samples")
    rng = random.Random(seed)
    np.random.seed(seed)

    # Common environmental variables
    metadata = {
        "sample_id": [f"sample_{i:03d}" for i in range(n_samples)],
        "temperature": np.random.normal(20.0, 5.0, n_samples),
        "humidity": np.random.normal(60.0, 15.0, n_samples),
        "precipitation": np.random.exponential(scale=50.0, size=n_samples),
        "elevation": np.random.normal(500.0, 200.0, n_samples),
        "latitude": np.random.uniform(-90.0, 90.0, n_samples),
        "longitude": np.random.uniform(-180.0, 180.0, n_samples),
        "soil_ph": np.random.normal(6.5, 1.0, n_samples),
    }

    # Ensure reasonable ranges
    metadata["humidity"] = np.clip(metadata["humidity"], 0, 100)
    metadata["precipitation"] = np.maximum(metadata["precipitation"], 0)
    metadata["elevation"] = np.maximum(metadata["elevation"], 0)
    metadata["soil_ph"] = np.clip(metadata["soil_ph"], 3.0, 10.0)

    # Create DataFrame
    df = pd.DataFrame(metadata)

    csv_file = output_dir / "environmental_metadata.csv"
    df.to_csv(csv_file, index=False)

    logger.info(f"Environmental metadata saved to {csv_file}")

    return {
        "type": "environmental",
        "n_samples": n_samples,
        "n_variables": len(metadata) - 1,  # Exclude sample_id
        "output_file": str(csv_file),
    }


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Ecology simulation for testing and validation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Simulate species abundance
  %(prog)s --type abundance --n-species 50 --n-samples 20

  # Simulate community composition
  %(prog)s --type community --n-species 100 --n-samples 30 --diversity 0.7

  # Simulate environmental metadata
  %(prog)s --type environmental --n-samples 25
        """,
    )

    parser.add_argument(
        "--type",
        required=True,
        choices=["abundance", "community", "environmental"],
        help="Type of simulation to run",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/simulation/ecology"),
        help="Output directory (default: output/simulation/ecology)",
    )
    parser.add_argument("--n-species", type=int, default=50, help="Number of species (abundance/community types)")
    parser.add_argument("--n-samples", type=int, default=20, help="Number of samples")
    parser.add_argument("--diversity", type=float, default=0.7, help="Diversity index (community type, 0-1)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")

    args = parser.parse_args()

    if args.verbose:
        import logging as std_logging

        logger.setLevel(std_logging.DEBUG)

    output_dir = paths.ensure_directory(args.output)

    try:
        if args.type == "abundance":
            results = simulate_abundance(output_dir, args.n_species, args.n_samples, args.seed)
        elif args.type == "community":
            results = simulate_community(output_dir, args.n_species, args.n_samples, args.diversity, args.seed)
        elif args.type == "environmental":
            results = simulate_environmental(output_dir, args.n_samples, args.seed)

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
