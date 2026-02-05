#!/usr/bin/env python3
"""Visualization workflow orchestrator.

This script provides comprehensive orchestration for visualization workflows,
including plot generation from data files, heatmaps, line plots, and animations.

Usage:
    python3 scripts/visualization/run_visualization.py --input data.csv --plot-type lineplot --output output/visualization/plots
    python3 scripts/visualization/run_visualization.py --input matrix.csv --plot-type heatmap --output output/visualization
    python3 scripts/visualization/run_visualization.py --help
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Any

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core.io import ensure_directory, dump_json, read_csv, load_json
from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Visualization workflow orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Line plot from CSV data
  %(prog)s --input data.csv --plot-type lineplot --output output/visualization/lineplot

  # Heatmap from matrix
  %(prog)s --input matrix.csv --plot-type heatmap --output output/visualization/heatmap

  # Animated time series
  %(prog)s --input timeseries.json --plot-type animation --output output/visualization/animation
        """
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Input data file (CSV, TSV, or JSON)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/visualization"),
        help="Output directory (default: output/visualization)",
    )
    parser.add_argument(
        "--plot-type",
        type=str,
        required=True,
        choices=["lineplot", "heatmap", "animation", "histogram"],
        help="Type of plot to generate",
    )
    parser.add_argument(
        "--x-column",
        type=str,
        help="Column name for x-axis (CSV input)",
    )
    parser.add_argument(
        "--y-column",
        type=str,
        help="Column name for y-axis (CSV input)",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Output resolution DPI (default: 300)",
    )
    parser.add_argument(
        "--format",
        type=str,
        default="png",
        choices=["png", "pdf", "svg"],
        help="Output format (default: png)",
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


def run_workflow(args):
    """Execute visualization workflow."""
    logger.info("Starting visualization workflow")
    logger.info(f"Input: {args.input}")
    logger.info(f"Plot type: {args.plot_type}")
    logger.info(f"Output: {args.output}")
    
    if not args.input.exists():
        raise FileNotFoundError(f"Input file not found: {args.input}")
    
    if args.dry_run:
        logger.info("DRY RUN - no changes will be made")
        return
    
    output_dir = ensure_directory(args.output)
    logger.info(f"Output directory: {output_dir}")
    
    # Load data
    logger.info(f"Loading data from {args.input}")
    try:
        if args.input.suffix == ".json":
            data = load_json(args.input)
        else:
            df = read_csv(args.input)
            data = df
    except Exception as e:
        logger.error(f"Failed to load data: {e}")
        raise
    
    workflow_results = {
        "input_file": str(args.input),
        "plot_type": args.plot_type,
        "output_dir": str(output_dir),
        "plots_generated": [],
    }
    
    # Generate plots
    try:
        if args.plot_type == "lineplot":
            logger.info("Generating line plot...")
            from metainformant.visualization import lineplot
            
            if isinstance(data, dict):
                # JSON data
                y_data = data.get("values", list(data.values())[0] if data else [])
                x_data = data.get("x", None)
            else:
                # DataFrame
                if args.y_column and args.y_column in data.columns:
                    y_data = data[args.y_column].values
                else:
                    y_data = data.iloc[:, -1].values  # Use last column
                
                if args.x_column and args.x_column in data.columns:
                    x_data = data[args.x_column].values
                else:
                    x_data = None
            
            ax = lineplot(x_data, y_data)
            ax.set_xlabel(args.x_column or "Index")
            ax.set_ylabel(args.y_column or "Value")
            ax.set_title(f"Line Plot: {args.input.stem}")
            
            output_file = output_dir / f"{args.input.stem}_lineplot.{args.format}"
            ax.figure.savefig(output_file, dpi=args.dpi, format=args.format)
            logger.info(f"Line plot saved to {output_file}")
            workflow_results["plots_generated"].append(str(output_file))
        
        elif args.plot_type == "heatmap":
            logger.info("Generating heatmap...")
            from metainformant.visualization import heatmap
            
            if isinstance(data, dict):
                # Try to extract matrix from JSON
                import numpy as np
                matrix = np.array(data.get("matrix", list(data.values())[0] if data else []))
            else:
                # DataFrame
                matrix = data.values
            
            ax = heatmap(matrix, cmap="viridis")
            ax.set_title(f"Heatmap: {args.input.stem}")
            
            output_file = output_dir / f"{args.input.stem}_heatmap.{args.format}"
            ax.figure.savefig(output_file, dpi=args.dpi, format=args.format)
            logger.info(f"Heatmap saved to {output_file}")
            workflow_results["plots_generated"].append(str(output_file))
        
        elif args.plot_type == "animation":
            logger.info("Generating animation...")
            from metainformant.visualization import animate_time_series
            
            if isinstance(data, dict):
                series = data.get("series", list(data.values())[0] if data else [])
            else:
                series = data.iloc[:, -1].values  # Use last column
            
            fig, anim = animate_time_series(series, interval_ms=100)
            
            output_file = output_dir / f"{args.input.stem}_animation.gif"
            try:
                from matplotlib.animation import PillowWriter
                writer = PillowWriter(fps=10)
                anim.save(output_file, writer=writer)
                logger.info(f"Animation saved to {output_file}")
                workflow_results["plots_generated"].append(str(output_file))
            except Exception as e:
                logger.warning(f"Could not save animation: {e}")
        
        elif args.plot_type == "histogram":
            logger.info("Generating histogram...")
            import matplotlib.pyplot as plt
            
            if isinstance(data, dict):
                values = data.get("values", list(data.values())[0] if data else [])
            else:
                values = data.iloc[:, -1].values  # Use last column
            
            fig, ax = plt.subplots()
            ax.hist(values, bins=30, edgecolor='black')
            ax.set_xlabel("Value")
            ax.set_ylabel("Frequency")
            ax.set_title(f"Histogram: {args.input.stem}")
            
            output_file = output_dir / f"{args.input.stem}_histogram.{args.format}"
            fig.savefig(output_file, dpi=args.dpi, format=args.format)
            logger.info(f"Histogram saved to {output_file}")
            workflow_results["plots_generated"].append(str(output_file))
        
    except Exception as e:
        logger.error(f"Plot generation failed: {e}", exc_info=True)
        workflow_results["error"] = str(e)
    
    # Save summary
    summary_file = output_dir / "workflow_summary.json"
    dump_json(workflow_results, summary_file, indent=2)
    logger.info(f"Workflow summary saved to {summary_file}")
    
    logger.info("Workflow complete")


def main():
    """Main entry point."""
    args = parse_args()
    
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
