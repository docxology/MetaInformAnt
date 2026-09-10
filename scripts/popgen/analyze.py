#!/usr/bin/env python3
"""Thin CLI wrapper for the population genetics dataset workflow.

Business logic lives in ``metainformant.popgen.workflow`` (dataset generation,
analysis, reporting, visualization); this script only bootstraps the src tree
and delegates.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core.io import ensure_directory, load_json
from metainformant.core.utils.logging import setup_logger
from metainformant.popgen.workflow import (
    analyze_dataset,
    generate_comprehensive_dataset,
    generate_summary_report,
    generate_visualizations,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Population genetics dataset generation and analysis")
    parser.add_argument(
        "--output-dir",
        type=str,
        default="output/popgen",
        help="Output directory (default: output/popgen)",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    parser.add_argument(
        "--n-sequences",
        type=int,
        default=50,
        help="Number of sequences per scenario (default: 50)",
    )
    parser.add_argument("--sequence-length", type=int, default=5000, help="Sequence length (default: 5000)")
    parser.add_argument("--skip-generation", action="store_true", help="Skip generation, only analyze existing dataset")

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    ensure_directory(str(output_dir))

    logger = setup_logger("metainformant.popgen.comprehensive")

    # Generate dataset
    if not args.skip_generation:
        logger.info("Step 1: Generating comprehensive dataset")
        dataset_info = generate_comprehensive_dataset(
            output_dir,
            seed=args.seed,
            n_sequences_per_scenario=args.n_sequences,
            sequence_length=args.sequence_length,
        )
    else:
        logger.info("Skipping generation, loading existing dataset info")
        dataset_info = load_json(str(output_dir / "dataset_info.json"))

    # Analyze dataset (business logic in metainformant.popgen)
    logger.info("Step 2: Analyzing dataset")
    analysis_results = analyze_dataset(dataset_info, output_dir)

    # Generate report
    logger.info("Step 3: Generating summary report")
    generate_summary_report(dataset_info, analysis_results, output_dir)

    # Generate visualizations
    logger.info("Step 4: Generating visualizations")
    generate_visualizations(analysis_results, output_dir)

    logger.info("Comprehensive analysis complete!")
    logger.info(f"Results saved to: {output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
