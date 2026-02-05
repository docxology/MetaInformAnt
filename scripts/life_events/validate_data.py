#!/usr/bin/env python3
"""Validate event sequence data.

This script validates event sequences for common issues and generates validation reports.

Usage:
    python3 scripts/life_events/validate_data.py --input data/sequences.json --output output/life_events/validation_report.json
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core import io, paths
from metainformant.core.utils.logging import setup_logger
from metainformant.life_events import (
    get_event_statistics,
    load_sequences_from_json,
    validate_sequence,
)

logger = setup_logger(__name__, level="INFO")


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate event sequence data",
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
        default=Path("output/life_events/validation_report.json"),
        help="Output validation report file (default: output/life_events/validation_report.json)",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Treat warnings as errors",
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

    logger.info("Validating sequences...")

    validation_results = []
    invalid_count = 0
    warning_count = 0

    for seq in sequences:
        is_valid, errors = validate_sequence(seq)
        warnings = [e for e in errors if e.startswith("warning")]
        errors_only = [e for e in errors if not e.startswith("warning")]

        result = {
            "person_id": seq.person_id,
            "is_valid": is_valid and (not args.strict or len(warnings) == 0),
            "errors": errors_only,
            "warnings": warnings,
        }
        validation_results.append(result)

        if not result["is_valid"]:
            invalid_count += 1
        if warnings:
            warning_count += len(warnings)

    # Get statistics
    stats = get_event_statistics(sequences)

    report = {
        "total_sequences": len(sequences),
        "valid_sequences": len(sequences) - invalid_count,
        "invalid_sequences": invalid_count,
        "total_warnings": warning_count,
        "statistics": stats,
        "validation_results": validation_results,
    }

    paths.ensure_directory(args.output.parent)
    io.dump_json(report, args.output, indent=2)

    logger.info(f"✅ Validation complete. Report saved to {args.output}")
    logger.info(f"   Valid sequences: {report['valid_sequences']}/{report['total_sequences']}")
    logger.info(f"   Invalid sequences: {report['invalid_sequences']}")
    logger.info(f"   Total warnings: {report['total_warnings']}")

    if invalid_count > 0:
        logger.warning(f"⚠️  Found {invalid_count} invalid sequences")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
