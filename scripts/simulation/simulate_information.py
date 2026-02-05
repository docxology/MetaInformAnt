#!/usr/bin/env python3
"""Information theory simulation script.

This script generates synthetic data for information-theoretic analysis
including sequence information content and mutual information test data.

Usage:
    python3 scripts/simulation/simulate_information.py --type sequences --n-sequences 100 --length 1000
    python3 scripts/simulation/simulate_information.py --type mutual --n-samples 200
"""

import argparse
import random
import sys
from pathlib import Path

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core import io, logging, paths, validation
from metainformant.simulation.sequences import generate_random_dna

logger = logging.get_logger(__name__)


def simulate_sequences(
    output_dir: Path,
    n_sequences: int,
    sequence_length: int,
    gc_content: float,
    complexity_level: str,
    seed: int,
) -> dict:
    """Simulate sequences with varying information content."""
    logger.info(f"Generating sequences: {n_sequences} sequences, length {sequence_length}")
    rng = random.Random(seed)

    sequences = []

    for seq_idx in range(n_sequences):
        if complexity_level == "low":
            # Low complexity: repetitive sequences
            base = rng.choice("ACGT")
            seq = base * sequence_length
        elif complexity_level == "high":
            # High complexity: random sequence
            seq = generate_random_dna(sequence_length, gc_content=gc_content, rng=rng)
        else:  # medium
            # Medium complexity: some repetition
            seq = ""
            for _ in range(sequence_length // 10):
                base = rng.choice("ACGT")
                seq += base * 10
            # Fill remainder
            seq += generate_random_dna(sequence_length - len(seq), gc_content=gc_content, rng=rng)

        sequences.append(
            {
                "sequence_id": f"seq_{seq_idx:04d}",
                "sequence": seq,
                "length": len(seq),
                "complexity": complexity_level,
            }
        )

    # Save as JSON
    json_file = output_dir / "information_sequences.json"
    io.dump_json(sequences, json_file, indent=2)

    # Also save as FASTA
    from metainformant.dna.sequences import write_fasta

    fasta_dict = {s["sequence_id"]: s["sequence"] for s in sequences}
    fasta_file = output_dir / "information_sequences.fasta"
    write_fasta(fasta_dict, str(fasta_file))

    logger.info(f"Information sequences saved to {json_file}")

    return {
        "type": "sequences",
        "n_sequences": n_sequences,
        "length": sequence_length,
        "complexity_level": complexity_level,
        "json_file": str(json_file),
        "fasta_file": str(fasta_file),
    }


def simulate_mutual(
    output_dir: Path,
    n_samples: int,
    correlation_strength: float,
    seed: int,
) -> dict:
    """Simulate data for mutual information analysis."""
    logger.info(f"Generating mutual information test data: {n_samples} samples")
    rng = random.Random(seed)

    import numpy as np

    np.random.seed(seed)

    # Generate correlated variables
    # X: random variable
    X = np.random.normal(0, 1, n_samples)

    # Y: correlated with X
    Y = correlation_strength * X + np.sqrt(1 - correlation_strength**2) * np.random.normal(0, 1, n_samples)

    # Z: independent of X
    Z = np.random.normal(0, 1, n_samples)

    # Create DataFrame
    import pandas as pd

    df = pd.DataFrame(
        {
            "sample_id": [f"sample_{i:04d}" for i in range(n_samples)],
            "X": X,
            "Y": Y,
            "Z": Z,
        }
    )

    csv_file = output_dir / "mutual_information_data.csv"
    df.to_csv(csv_file, index=False)

    # Calculate theoretical MI (for validation)
    # MI(X;Y) ≈ -0.5 * log(1 - correlation^2) for Gaussian
    theoretical_mi = -0.5 * np.log(1 - correlation_strength**2) if correlation_strength < 1.0 else float("inf")

    metadata = {
        "n_samples": n_samples,
        "correlation_strength": correlation_strength,
        "theoretical_mi_xy": float(theoretical_mi),
        "theoretical_mi_xz": 0.0,  # Independent
    }

    metadata_file = output_dir / "mi_metadata.json"
    io.dump_json(metadata, metadata_file, indent=2)

    logger.info(f"Mutual information data saved to {csv_file}")

    return {
        "type": "mutual",
        "n_samples": n_samples,
        "correlation_strength": correlation_strength,
        "output_file": str(csv_file),
        "metadata_file": str(metadata_file),
    }


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Information theory simulation for testing and validation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Simulate sequences with varying complexity
  %(prog)s --type sequences --n-sequences 100 --length 1000 --complexity high

  # Simulate mutual information test data
  %(prog)s --type mutual --n-samples 200 --correlation 0.7
        """,
    )

    parser.add_argument(
        "--type",
        required=True,
        choices=["sequences", "mutual"],
        help="Type of simulation to run",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/simulation/information"),
        help="Output directory (default: output/simulation/information)",
    )
    parser.add_argument("--n-sequences", type=int, default=100, help="Number of sequences (sequences type)")
    parser.add_argument("--length", type=int, default=1000, help="Sequence length (sequences type)")
    parser.add_argument("--gc-content", type=float, default=0.5, help="GC content (sequences type)")
    parser.add_argument(
        "--complexity",
        type=str,
        choices=["low", "medium", "high"],
        default="high",
        help="Sequence complexity level (sequences type)",
    )
    parser.add_argument("--n-samples", type=int, default=200, help="Number of samples (mutual type)")
    parser.add_argument("--correlation", type=float, default=0.7, help="Correlation strength (mutual type, 0-1)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")

    args = parser.parse_args()

    if args.verbose:
        import logging as std_logging

        logger.setLevel(std_logging.DEBUG)

    output_dir = paths.ensure_directory(args.output)

    try:
        if args.type == "sequences":
            results = simulate_sequences(
                output_dir, args.n_sequences, args.length, args.gc_content, args.complexity, args.seed
            )
        elif args.type == "mutual":
            results = simulate_mutual(output_dir, args.n_samples, args.correlation, args.seed)

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
