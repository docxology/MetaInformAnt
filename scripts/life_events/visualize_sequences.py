#!/usr/bin/env python3
"""Visualize event sequences.

This script generates various visualizations for event sequences.

Usage:
    python3 scripts/life_events/visualize_sequences.py --input data/sequences.json --output output/life_events/visualizations/
    python3 scripts/life_events/visualize_sequences.py --input data/sequences.json --all-visualizations
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core import paths
from metainformant.core.utils.logging import setup_logger
from metainformant.life_events import (
    load_sequences_from_json,
    plot_domain_distribution,
    plot_event_timeline,
    plot_sequence_length_distribution,
    plot_temporal_density,
)

logger = setup_logger(__name__, level="INFO")


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Visualize event sequences",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Input sequences file (JSON format)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/life_events/visualizations"),
        help="Output directory (default: output/life_events/visualizations)",
    )
    parser.add_argument(
        "--all-visualizations",
        action="store_true",
        help="Generate all visualization types",
    )
    parser.add_argument(
        "--timeline",
        action="store_true",
        help="Generate timeline visualizations",
    )
    parser.add_argument(
        "--domain-distribution",
        action="store_true",
        help="Generate domain distribution plot",
    )
    parser.add_argument(
        "--temporal-density",
        action="store_true",
        help="Generate temporal density plot",
    )
    parser.add_argument(
        "--sequence-lengths",
        action="store_true",
        help="Generate sequence length distribution",
    )
    return parser.parse_args()


def main():
    """Main function."""
    args = parse_args()
    
    if not args.input.exists():
        logger.error(f"Input file not found: {args.input}")
        return 1
    
    logger.info(f"Loading sequences from {args.input}")
    sequences = load_sequences_from_json(args.input)
    logger.info(f"✅ Loaded {len(sequences)} sequences")
    
    paths.ensure_directory(args.output)
    
    generate_all = args.all_visualizations or (
        not args.timeline
        and not args.domain_distribution
        and not args.temporal_density
        and not args.sequence_lengths
    )
    
    if generate_all or args.timeline:
        logger.info("Generating timeline visualizations...")
        for i in range(min(5, len(sequences))):
            plot_event_timeline(
                sequences[i],
                output_path=args.output / f"timeline_sequence_{i+1}.png",
            )
        logger.info("✅ Timeline visualizations generated")
    
    if generate_all or args.domain_distribution:
        logger.info("Generating domain distribution plot...")
        plot_domain_distribution(
            sequences,
            output_path=args.output / "domain_distribution.png",
        )
        logger.info("✅ Domain distribution plot generated")
    
    if generate_all or args.temporal_density:
        logger.info("Generating temporal density plot...")
        plot_temporal_density(
            sequences,
            output_path=args.output / "temporal_density.png",
        )
        logger.info("✅ Temporal density plot generated")
    
    if generate_all or args.sequence_lengths:
        logger.info("Generating sequence length distribution...")
        plot_sequence_length_distribution(
            sequences,
            output_path=args.output / "sequence_length_distribution.png",
        )
        logger.info("✅ Sequence length distribution generated")
    
    logger.info(f"✅ All visualizations saved to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

