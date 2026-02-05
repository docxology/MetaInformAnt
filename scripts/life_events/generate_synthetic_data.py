#!/usr/bin/env python3
"""Generate synthetic life event sequences.

This script generates realistic synthetic event sequences for testing and demos.

Usage:
    python3 scripts/life_events/generate_synthetic_data.py --n-sequences 100 --output output/life_events/synthetic.json
    python3 scripts/life_events/generate_synthetic_data.py --realistic --seasonal-patterns --n-sequences 50
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core import io, paths
from metainformant.core.utils.logging import setup_logger
from metainformant.life_events import (
    generate_realistic_life_events,
    generate_synthetic_life_events,
)

logger = setup_logger(__name__, level="INFO")


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate synthetic life event sequences",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--n-sequences",
        type=int,
        default=100,
        help="Number of sequences to generate (default: 100)",
    )
    parser.add_argument(
        "--min-events",
        type=int,
        default=5,
        help="Minimum events per sequence (default: 5)",
    )
    parser.add_argument(
        "--max-events",
        type=int,
        default=30,
        help="Maximum events per sequence (default: 30)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/life_events/synthetic_sequences.json"),
        help="Output file path (default: output/life_events/synthetic_sequences.json)",
    )
    parser.add_argument(
        "--generate-outcomes",
        action="store_true",
        help="Generate outcome labels",
    )
    parser.add_argument(
        "--outcome-relationship",
        type=str,
        default="complex",
        choices=["random", "health_focused", "education_focused", "complex"],
        help="Outcome relationship pattern (default: complex)",
    )
    parser.add_argument(
        "--realistic",
        action="store_true",
        help="Use realistic generation with dependencies and patterns",
    )
    parser.add_argument(
        "--seasonal-patterns",
        action="store_true",
        help="Apply seasonal patterns (requires --realistic)",
    )
    parser.add_argument(
        "--temporal-noise",
        type=float,
        default=0.0,
        help="Level of temporal noise (0.0-1.0, default: 0.0)",
    )
    parser.add_argument(
        "--missing-data",
        type=float,
        default=0.0,
        help="Probability of missing events (0.0-1.0, default: 0.0)",
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=42,
        help="Random seed (default: 42)",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose logging",
    )
    return parser.parse_args()


def main():
    """Main function."""
    args = parse_args()

    if args.verbose:
        logger.setLevel("DEBUG")

    logger.info(f"Generating {args.n_sequences} synthetic event sequences...")

    # Ensure output directory exists
    paths.ensure_directory(args.output.parent)

    # Generate sequences
    if args.realistic:
        sequences, outcomes = generate_realistic_life_events(
            n_sequences=args.n_sequences,
            min_events_per_sequence=args.min_events,
            max_events_per_sequence=args.max_events,
            generate_outcomes=args.generate_outcomes,
            outcome_relationship=args.outcome_relationship,
            seasonal_patterns=args.seasonal_patterns,
            temporal_noise=args.temporal_noise,
            missing_data_probability=args.missing_data,
            random_state=args.random_state,
        )
    else:
        sequences, outcomes = generate_synthetic_life_events(
            n_sequences=args.n_sequences,
            min_events_per_sequence=args.min_events,
            max_events_per_sequence=args.max_events,
            generate_outcomes=args.generate_outcomes,
            outcome_relationship=args.outcome_relationship,
            random_state=args.random_state,
        )

    # Save sequences
    sequences_data = [seq.to_dict() for seq in sequences]
    io.dump_json(sequences_data, args.output, indent=2)
    logger.info(f"✅ Saved {len(sequences)} sequences to {args.output}")

    # Save outcomes if generated
    if outcomes is not None:
        outcomes_file = args.output.parent / f"{args.output.stem}_outcomes.json"
        io.dump_json(
            {
                "outcomes": outcomes.tolist() if hasattr(outcomes, "tolist") else list(outcomes),
                "n_sequences": len(outcomes),
            },
            outcomes_file,
            indent=2,
        )
        logger.info(f"✅ Saved outcomes to {outcomes_file}")

    logger.info("✅ Synthetic data generation complete!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
