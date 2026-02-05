#!/usr/bin/env python3
"""Compare event patterns between two population groups.

This script compares event sequences from two groups and generates comparison visualizations.

Usage:
    python3 scripts/life_events/compare_groups.py --group1 data/group1.json --group2 data/group2.json --output output/life_events/comparison/
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core import io, paths
from metainformant.core.utils.logging import setup_logger
from metainformant.life_events import (
    compare_populations,
    load_sequences_from_json,
    plot_population_comparison,
)

logger = setup_logger(__name__, level="INFO")


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Compare event patterns between two groups",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--group1",
        type=Path,
        required=True,
        help="Sequences file for group 1 (JSON format)",
    )
    parser.add_argument(
        "--group2",
        type=Path,
        required=True,
        help="Sequences file for group 2 (JSON format)",
    )
    parser.add_argument(
        "--group1-label",
        type=str,
        default="Group 1",
        help="Label for group 1 (default: Group 1)",
    )
    parser.add_argument(
        "--group2-label",
        type=str,
        default="Group 2",
        help="Label for group 2 (default: Group 2)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/life_events/comparison"),
        help="Output directory (default: output/life_events/comparison)",
    )
    return parser.parse_args()


def main():
    """Main function."""
    args = parse_args()

    if not args.group1.exists():
        logger.error(f"Group 1 file not found: {args.group1}")
        return 1

    if not args.group2.exists():
        logger.error(f"Group 2 file not found: {args.group2}")
        return 1

    logger.info(f"Loading group 1 sequences from {args.group1}")
    sequences_group1 = load_sequences_from_json(args.group1)
    logger.info(f"✅ Loaded {len(sequences_group1)} sequences for group 1")

    logger.info(f"Loading group 2 sequences from {args.group2}")
    sequences_group2 = load_sequences_from_json(args.group2)
    logger.info(f"✅ Loaded {len(sequences_group2)} sequences for group 2")

    paths.ensure_directory(args.output)

    logger.info("Comparing populations...")
    comparison = compare_populations(
        sequences_group1,
        sequences_group2,
        output_dir=args.output,
    )

    logger.info("Generating comparison visualization...")
    plot_population_comparison(
        sequences_group1,
        sequences_group2,
        group1_label=args.group1_label,
        group2_label=args.group2_label,
        output_path=args.output / "population_comparison.png",
    )

    logger.info(f"✅ Comparison complete. Results saved to {args.output}")
    logger.info(f"   Common event types: {len(comparison['comparison']['common_event_types'])}")
    logger.info(f"   Unique to group 1: {len(comparison['comparison']['unique_to_group1'])}")
    logger.info(f"   Unique to group 2: {len(comparison['comparison']['unique_to_group2'])}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
