#!/usr/bin/env python3
"""Generate synthetic life event sequences.

This script generates realistic synthetic event sequences for testing and demos.

Usage:
    python3 scripts/life_events/generate_synthetic_data.py --n-sequences 100 --output output/life_events/synthetic.json
    python3 scripts/life_events/generate_synthetic_data.py --realistic --seasonal-patterns --n-sequences 50
"""

import argparse
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core import io  # noqa: E402
from metainformant.core.io import paths  # noqa: E402
from metainformant.core.utils.logging import setup_logger  # noqa: E402
from metainformant.life_events.core.utils import add_temporal_noise, generate_synthetic_life_events  # noqa: E402

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
        default="random",
        choices=["random", "event_based", "temporal"],
        help=(
            "How outcomes relate to the sequences: random (coin flip), event_based (1 if the sequence has a "
            "diagnosis, divorce or financial difficulty), temporal (1 if it has more than 12 events) "
            "(default: random)"
        ),
    )
    parser.add_argument(
        "--realistic",
        action="store_true",
        help="Use realistic generation with dependencies and patterns",
    )
    parser.add_argument(
        "--seasonal-patterns",
        action="store_true",
        help=(
            "Apply seasonal patterns (requires --realistic; not implemented by metainformant.life_events, "
            "so it only logs a warning)"
        ),
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

    # Generate sequences (generate_synthetic_life_events already builds domain-aware event chains)
    sequences, outcomes = generate_synthetic_life_events(
        n_sequences=args.n_sequences,
        min_events_per_sequence=args.min_events,
        max_events_per_sequence=args.max_events,
        generate_outcomes=args.generate_outcomes,
        outcome_relationship=args.outcome_relationship,
        random_state=args.random_state,
    )

    if args.realistic:
        # Realistic mode perturbs the generated timelines per sequence with the library's add_temporal_noise.
        if args.seasonal_patterns:
            logger.warning("--seasonal-patterns is not implemented by metainformant.life_events; ignoring it")
        if args.temporal_noise > 0 or args.missing_data > 0:
            rng = random.Random(args.random_state)
            sequences = [
                add_temporal_noise(
                    seq,
                    noise_level=args.temporal_noise,
                    rng=rng,
                    missing_probability=args.missing_data,
                )
                for seq in sequences
            ]

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
