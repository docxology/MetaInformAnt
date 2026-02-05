#!/usr/bin/env python3
"""Phenotype analysis workflow orchestrator.

This script provides comprehensive orchestration for phenotype analysis workflows,
including trait data loading, statistics calculation, and correlation analysis.

Usage:
    python3 scripts/phenotype/run_phenotype_analysis.py --input data/phenotype/traits.json --output output/phenotype/results
    python3 scripts/phenotype/run_phenotype_analysis.py --input traits.tsv --analyze-correlations --analyze-statistics
    python3 scripts/phenotype/run_phenotype_analysis.py --help
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Any

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core.io import dump_json, ensure_directory, load_json, read_csv
from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Phenotype analysis workflow orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Load AntWiki JSON data
  %(prog)s --input antwiki_species.json --output output/phenotype/antwiki

  # Analyze CSV/TSV phenotype data
  %(prog)s --input traits.csv --analyze-statistics --analyze-correlations

  # Basic statistics only
  %(prog)s --input traits.tsv --analyze-statistics
        """,
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Input phenotype data (JSON, CSV, or TSV format)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/phenotype"),
        help="Output directory (default: output/phenotype)",
    )
    parser.add_argument(
        "--analyze-statistics",
        action="store_true",
        help="Calculate basic statistics for numeric traits",
    )
    parser.add_argument(
        "--analyze-correlations",
        action="store_true",
        help="Calculate correlations between traits",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose logging",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without executing",
    )
    return parser.parse_args()


def load_phenotype_data(input_path: Path) -> dict[str, Any]:
    """Load phenotype data from various formats.

    Args:
        input_path: Path to input file

    Returns:
        Dictionary with loaded data and metadata
    """
    logger.info(f"Loading phenotype data from {input_path}")

    if input_path.suffix == ".json":
        # Try AntWiki format first
        try:
            from metainformant.phenotype.antwiki import load_antwiki_json

            data = load_antwiki_json(input_path)
            return {
                "format": "antwiki_json",
                "data": data,
                "num_entries": len(data),
            }
        except Exception:
            # Fall back to generic JSON
            data = load_json(input_path)
            return {
                "format": "json",
                "data": data,
            }
    elif input_path.suffix in [".csv", ".tsv"]:
        # Load CSV/TSV
        try:
            df = read_csv(input_path, sep="\t" if input_path.suffix == ".tsv" else ",")
            return {
                "format": "tabular",
                "data": df.to_dict("records"),
                "columns": list(df.columns),
                "num_rows": len(df),
            }
        except Exception as e:
            logger.error(f"Failed to load CSV/TSV: {e}")
            raise
    else:
        raise ValueError(f"Unsupported file format: {input_path.suffix}")


def analyze_statistics(data: dict[str, Any], output_dir: Path) -> dict[str, Any]:
    """Calculate basic statistics for numeric traits.

    Args:
        data: Loaded phenotype data
        output_dir: Output directory for results

    Returns:
        Dictionary with statistics results
    """
    logger.info("Calculating trait statistics...")

    results = {
        "statistics": {},
        "numeric_traits": [],
    }

    if data["format"] == "antwiki_json":
        # Analyze AntWiki measurements
        all_measurements = {}
        for entry in data["data"]:
            measurements = entry.get("measurements", {})
            for trait, values in measurements.items():
                if trait not in all_measurements:
                    all_measurements[trait] = []

                # Handle both single values and lists
                if isinstance(values, list):
                    all_measurements[trait].extend(values)
                elif isinstance(values, (int, float)):
                    all_measurements[trait].append(values)

        # Calculate statistics for each measurement type
        for trait, values in all_measurements.items():
            if values and all(isinstance(v, (int, float)) for v in values):
                numeric_values = [float(v) for v in values]
                results["statistics"][trait] = {
                    "mean": sum(numeric_values) / len(numeric_values),
                    "min": min(numeric_values),
                    "max": max(numeric_values),
                    "count": len(numeric_values),
                }
                results["numeric_traits"].append(trait)

    elif data["format"] == "tabular":
        # Analyze tabular data
        import pandas as pd

        df = pd.DataFrame(data["data"])
        numeric_cols = df.select_dtypes(include=["number"]).columns

        for col in numeric_cols:
            results["statistics"][col] = {
                "mean": float(df[col].mean()),
                "std": float(df[col].std()),
                "min": float(df[col].min()),
                "max": float(df[col].max()),
                "count": int(df[col].count()),
            }
            results["numeric_traits"].append(col)

    # Save results
    output_file = output_dir / "statistics_analysis.json"
    dump_json(results, output_file)
    logger.info(f"Statistics analysis saved to {output_file}")

    return results


def analyze_correlations(data: dict[str, Any], output_dir: Path) -> dict[str, Any]:
    """Calculate correlations between numeric traits.

    Args:
        data: Loaded phenotype data
        output_dir: Output directory for results

    Returns:
        Dictionary with correlation results
    """
    logger.info("Calculating trait correlations...")

    results = {
        "correlations": {},
    }

    try:
        import numpy as np
        import pandas as pd
    except ImportError:
        logger.warning("pandas/numpy required for correlation analysis")
        results["error"] = "pandas/numpy not available"
        return results

    if data["format"] == "antwiki_json":
        # Convert AntWiki data to dataframe for correlation
        records = []
        for entry in data["data"]:
            record = {"species": entry.get("species", "unknown")}
            measurements = entry.get("measurements", {})
            for trait, values in measurements.items():
                # Use mean if list, otherwise use value directly
                if isinstance(values, list) and values:
                    record[trait] = sum(v for v in values if isinstance(v, (int, float))) / len(
                        [v for v in values if isinstance(v, (int, float))]
                    )
                elif isinstance(values, (int, float)):
                    record[trait] = values
            records.append(record)

        if records:
            df = pd.DataFrame(records)
            numeric_cols = df.select_dtypes(include=["number"]).columns

            if len(numeric_cols) >= 2:
                corr_matrix = df[numeric_cols].corr()
                results["correlations"] = corr_matrix.to_dict()
                logger.info(f"Calculated correlations for {len(numeric_cols)} traits")
            else:
                results["note"] = "Need at least 2 numeric traits for correlation"
        else:
            results["note"] = "No valid records for correlation analysis"

    elif data["format"] == "tabular":
        df = pd.DataFrame(data["data"])
        numeric_cols = df.select_dtypes(include=["number"]).columns

        if len(numeric_cols) >= 2:
            corr_matrix = df[numeric_cols].corr()
            results["correlations"] = corr_matrix.to_dict()
            logger.info(f"Calculated correlations for {len(numeric_cols)} traits")
        else:
            results["note"] = "Need at least 2 numeric traits for correlation"

    # Save results
    output_file = output_dir / "correlation_analysis.json"
    dump_json(results, output_file)
    logger.info(f"Correlation analysis saved to {output_file}")

    return results


def run_workflow(args):
    """Execute phenotype analysis workflow.

    Args:
        args: Parsed command-line arguments
    """
    logger.info("Starting phenotype analysis workflow")
    logger.info(f"Input: {args.input}")
    logger.info(f"Output: {args.output}")

    if not args.input.exists():
        raise FileNotFoundError(f"Input file not found: {args.input}")

    if args.dry_run:
        logger.info("DRY RUN - no changes will be made")
        logger.info(f"Would analyze: {args.input}")
        logger.info(f"Would write to: {args.output}")
        if args.analyze_statistics:
            logger.info("Would calculate statistics")
        if args.analyze_correlations:
            logger.info("Would calculate correlations")
        return

    # Ensure output directory exists
    output_dir = ensure_directory(args.output)
    logger.info(f"Output directory: {output_dir}")

    # Load phenotype data
    try:
        data = load_phenotype_data(args.input)
        logger.info(
            f"Loaded data: {data.get('format', 'unknown')} format, {data.get('num_entries', data.get('num_rows', 'unknown'))} entries"
        )
    except Exception as e:
        logger.error(f"Failed to load phenotype data: {e}")
        raise

    # Run analyses
    workflow_results = {
        "input_file": str(args.input),
        "output_dir": str(output_dir),
        "data_format": data.get("format", "unknown"),
        "analyses": {},
    }

    # Statistics analysis (always run if no specific analysis requested, or if explicitly requested)
    if args.analyze_statistics or not args.analyze_correlations:
        try:
            stats_results = analyze_statistics(data, output_dir)
            workflow_results["analyses"]["statistics"] = stats_results
        except Exception as e:
            logger.error(f"Statistics analysis failed: {e}", exc_info=True)
            workflow_results["analyses"]["statistics"] = {"error": str(e)}

    # Correlation analysis
    if args.analyze_correlations:
        try:
            corr_results = analyze_correlations(data, output_dir)
            workflow_results["analyses"]["correlations"] = corr_results
        except Exception as e:
            logger.error(f"Correlation analysis failed: {e}", exc_info=True)
            workflow_results["analyses"]["correlations"] = {"error": str(e)}

    # Save workflow summary
    summary_file = output_dir / "workflow_summary.json"
    dump_json(workflow_results, summary_file, indent=2)
    logger.info(f"Workflow summary saved to {summary_file}")

    logger.info("Workflow complete")


def main():
    """Main entry point."""
    args = parse_args()

    # Setup logging level
    if args.verbose:
        logger.setLevel(logging.DEBUG)

    try:
        run_workflow(args)
        return 0
    except Exception as e:
        logger.error(f"Workflow failed: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
