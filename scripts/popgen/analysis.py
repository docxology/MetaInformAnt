#!/usr/bin/env python3
"""Comprehensive population genetics synthetic data generation and analysis.

This script generates large synthetic datasets with multiple demographic scenarios
and performs comprehensive population genetics analysis, demonstrating the full
capabilities of the METAINFORMANT population genetics modules.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Add project to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root / "src"))

from metainformant.core.io import ensure_directory
from metainformant.core.utils.logging import setup_logger

from generate_dataset import generate_comprehensive_dataset
from analyze import analyze_dataset
from visualize import generate_visualizations
from report import generate_summary_report


def main():
    """Main execution function."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Population genetics dataset generation and analysis"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="output/popgen",
        help="Output directory (default: output/popgen)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed (default: 42)",
    )
    parser.add_argument(
        "--n-sequences",
        type=int,
        default=50,
        help="Number of sequences per scenario (default: 50)",
    )
    parser.add_argument(
        "--sequence-length",
        type=int,
        default=5000,
        help="Sequence length (default: 5000)",
    )
    parser.add_argument(
        "--skip-generation",
        action="store_true",
        help="Skip generation, only analyze existing dataset",
    )
    
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
        from metainformant.core.io import load_json
        dataset_info = load_json(str(output_dir / "dataset_info.json"))
    
    # Analyze dataset
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


if __name__ == "__main__":
    main()