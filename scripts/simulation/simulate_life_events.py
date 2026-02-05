#!/usr/bin/env python3
"""Life events simulation script.

This script generates synthetic life event sequences and life course patterns
for testing life course analysis workflows.

Usage:
    python3 scripts/simulation/simulate_life_events.py --type sequences --n-sequences 100
    python3 scripts/simulation/simulate_life_events.py --type patterns --n-sequences 200
"""

import argparse
import random
import sys
from datetime import datetime, timedelta
from pathlib import Path

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core import io, logging, paths, validation

logger = logging.get_logger(__name__)


def simulate_sequences(
    output_dir: Path,
    n_sequences: int,
    min_events: int,
    max_events: int,
    seed: int,
) -> dict:
    """Simulate life event sequences."""
    logger.info(f"Generating life event sequences: {n_sequences} sequences")
    rng = random.Random(seed)

    # Event types by domain
    event_types = {
        "education": ["enrollment", "graduation", "degree", "certification"],
        "occupation": ["job_start", "job_change", "promotion", "retirement"],
        "health": ["diagnosis", "treatment", "recovery", "checkup"],
        "family": ["marriage", "divorce", "birth", "adoption"],
        "residence": ["move", "purchase", "rental"],
    }

    sequences = []
    start_date = datetime(2000, 1, 1)

    for seq_idx in range(n_sequences):
        person_id = f"person_{seq_idx:04d}"
        n_events = rng.randint(min_events, max_events)

        events = []
        current_date = start_date

        for event_idx in range(n_events):
            # Select domain and event type
            domain = rng.choice(list(event_types.keys()))
            event_type = rng.choice(event_types[domain])

            # Advance date
            days_advance = rng.randint(30, 365)
            current_date = current_date + timedelta(days=days_advance)

            events.append(
                {
                    "domain": domain,
                    "event_type": event_type,
                    "date": current_date.isoformat(),
                    "metadata": {
                        "sequence_position": event_idx,
                        "age": (current_date - start_date).days // 365,
                    },
                }
            )

        sequences.append(
            {
                "person_id": person_id,
                "events": events,
                "n_events": len(events),
            }
        )

    # Save as JSON
    json_file = output_dir / "life_event_sequences.json"
    io.dump_json(sequences, json_file, indent=2)

    logger.info(f"Life event sequences saved to {json_file}")

    return {
        "type": "sequences",
        "n_sequences": n_sequences,
        "output_file": str(json_file),
    }


def simulate_patterns(
    output_dir: Path,
    n_sequences: int,
    n_patterns: int,
    seed: int,
) -> dict:
    """Simulate life course patterns."""
    logger.info(f"Generating life course patterns: {n_sequences} sequences, {n_patterns} patterns")
    rng = random.Random(seed)

    # Define common life course patterns
    patterns = [
        {
            "name": "traditional_education",
            "events": [
                ("education", "enrollment"),
                ("education", "graduation"),
                ("occupation", "job_start"),
            ],
        },
        {
            "name": "career_focused",
            "events": [
                ("occupation", "job_start"),
                ("occupation", "promotion"),
                ("occupation", "job_change"),
            ],
        },
        {
            "name": "family_oriented",
            "events": [
                ("family", "marriage"),
                ("family", "birth"),
                ("residence", "purchase"),
            ],
        },
    ]

    # Extend patterns if needed
    while len(patterns) < n_patterns:
        pattern_name = f"pattern_{len(patterns)}"
        n_pattern_events = rng.randint(2, 5)
        pattern_events = []

        domains = ["education", "occupation", "health", "family", "residence"]
        event_types_map = {
            "education": ["enrollment", "graduation"],
            "occupation": ["job_start", "promotion"],
            "health": ["diagnosis", "treatment"],
            "family": ["marriage", "birth"],
            "residence": ["move", "purchase"],
        }

        for _ in range(n_pattern_events):
            domain = rng.choice(domains)
            event_type = rng.choice(event_types_map[domain])
            pattern_events.append((domain, event_type))

        patterns.append({"name": pattern_name, "events": pattern_events})

    # Generate sequences following patterns
    sequences = []
    start_date = datetime(2000, 1, 1)

    for seq_idx in range(n_sequences):
        person_id = f"person_{seq_idx:04d}"
        pattern = rng.choice(patterns)

        events = []
        current_date = start_date

        for domain, event_type in pattern["events"]:
            days_advance = rng.randint(180, 730)
            current_date = current_date + timedelta(days=days_advance)

            events.append(
                {
                    "domain": domain,
                    "event_type": event_type,
                    "date": current_date.isoformat(),
                    "pattern": pattern["name"],
                }
            )

        sequences.append(
            {
                "person_id": person_id,
                "pattern": pattern["name"],
                "events": events,
                "n_events": len(events),
            }
        )

    # Save sequences
    json_file = output_dir / "pattern_sequences.json"
    io.dump_json(sequences, json_file, indent=2)

    # Save pattern definitions
    patterns_file = output_dir / "patterns.json"
    io.dump_json(patterns, patterns_file, indent=2)

    logger.info(f"Life course patterns saved to {json_file}")

    return {
        "type": "patterns",
        "n_sequences": n_sequences,
        "n_patterns": len(patterns),
        "output_file": str(json_file),
        "patterns_file": str(patterns_file),
    }


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Life events simulation for testing and validation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Simulate event sequences
  %(prog)s --type sequences --n-sequences 100 --min-events 5 --max-events 20

  # Simulate life course patterns
  %(prog)s --type patterns --n-sequences 200 --n-patterns 5
        """,
    )

    parser.add_argument(
        "--type",
        required=True,
        choices=["sequences", "patterns"],
        help="Type of simulation to run",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/simulation/life_events"),
        help="Output directory (default: output/simulation/life_events)",
    )
    parser.add_argument("--n-sequences", type=int, default=100, help="Number of sequences")
    parser.add_argument("--min-events", type=int, default=5, help="Minimum events per sequence (sequences type)")
    parser.add_argument("--max-events", type=int, default=20, help="Maximum events per sequence (sequences type)")
    parser.add_argument("--n-patterns", type=int, default=3, help="Number of patterns (patterns type)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")

    args = parser.parse_args()

    if args.verbose:
        import logging as std_logging

        logger.setLevel(std_logging.DEBUG)

    output_dir = paths.ensure_directory(args.output)

    try:
        if args.type == "sequences":
            results = simulate_sequences(output_dir, args.n_sequences, args.min_events, args.max_events, args.seed)
        elif args.type == "patterns":
            results = simulate_patterns(output_dir, args.n_sequences, args.n_patterns, args.seed)

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
