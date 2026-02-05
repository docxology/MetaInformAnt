#!/usr/bin/env python3
"""Epigenome analysis workflow orchestrator.

This script provides comprehensive orchestration for epigenome analysis workflows,
including methylation analysis, chromatin track processing, and beta value computation.

Usage:
    python3 scripts/epigenome/run_epigenome_analysis.py --methylation methylation.tsv --output output/epigenome/results
    python3 scripts/epigenome/run_epigenome_analysis.py --bedgraph track.bedgraph --output output/epigenome/tracks
    python3 scripts/epigenome/run_epigenome_analysis.py --help
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Any

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core.io import ensure_directory, dump_json
from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Epigenome analysis workflow orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Methylation analysis
  %(prog)s --methylation cpg_table.tsv --compute-beta --summarize-by-chromosome

  # Chromatin track analysis
  %(prog)s --bedgraph track.bedgraph --output output/epigenome/tracks

  # Full methylation pipeline
  %(prog)s --methylation cpg.tsv --compute-beta --summarize-by-chromosome --output output/epigenome/methylation
        """
    )
    parser.add_argument(
        "--methylation",
        type=Path,
        help="Input CpG methylation table (TSV format)",
    )
    parser.add_argument(
        "--bedgraph",
        type=Path,
        help="Input bedGraph track file",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/epigenome"),
        help="Output directory (default: output/epigenome)",
    )
    parser.add_argument(
        "--compute-beta",
        action="store_true",
        help="Compute beta values from methylation counts",
    )
    parser.add_argument(
        "--summarize-by-chromosome",
        action="store_true",
        help="Summarize methylation by chromosome",
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
    """Execute epigenome analysis workflow."""
    logger.info("Starting epigenome analysis workflow")
    logger.info(f"Output: {args.output}")
    
    if args.dry_run:
        logger.info("DRY RUN - no changes will be made")
        if args.methylation:
            logger.info(f"Would analyze methylation: {args.methylation}")
        if args.bedgraph:
            logger.info(f"Would analyze bedGraph: {args.bedgraph}")
        return
    
    output_dir = ensure_directory(args.output)
    logger.info(f"Output directory: {output_dir}")
    
    workflow_results = {
        "output_dir": str(output_dir),
        "analyses": {},
    }
    
    # Methylation analysis
    if args.methylation:
        if not args.methylation.exists():
            raise FileNotFoundError(f"Methylation file not found: {args.methylation}")
        
        try:
            logger.info(f"Loading methylation data from {args.methylation}")
            from metainformant.epigenome import load_cpg_table, compute_beta_values, summarize_beta_by_chromosome
            
            cpg_df = load_cpg_table(args.methylation)
            logger.info(f"Loaded {len(cpg_df)} CpG sites")
            workflow_results["input_file"] = str(args.methylation)
            
            # Compute beta values
            if args.compute_beta:
                logger.info("Computing beta values...")
                cpg_df = compute_beta_values(cpg_df)
                
                # Save beta values
                beta_file = output_dir / "beta_values.tsv"
                cpg_df.to_csv(beta_file, sep="\t", index=False)
                logger.info(f"Beta values saved to {beta_file}")
                
                workflow_results["analyses"]["beta_values"] = {
                    "n_sites": len(cpg_df),
                    "output_file": str(beta_file),
                }
            
            # Summarize by chromosome
            if args.summarize_by_chromosome:
                logger.info("Summarizing by chromosome...")
                if args.compute_beta and 'beta' in cpg_df.columns:
                    summary_df = summarize_beta_by_chromosome(cpg_df)
                else:
                    # Compute beta first if not done
                    cpg_df = compute_beta_values(cpg_df)
                    summary_df = summarize_beta_by_chromosome(cpg_df)
                
                summary_file = output_dir / "chromosome_summary.tsv"
                summary_df.to_csv(summary_file, sep="\t", index=False)
                logger.info(f"Chromosome summary saved to {summary_file}")
                
                workflow_results["analyses"]["chromosome_summary"] = {
                    "n_chromosomes": len(summary_df),
                    "output_file": str(summary_file),
                }
        except Exception as e:
            logger.error(f"Methylation analysis failed: {e}", exc_info=True)
            workflow_results["analyses"]["methylation"] = {"error": str(e)}
    
    # BedGraph track analysis
    if args.bedgraph:
        if not args.bedgraph.exists():
            raise FileNotFoundError(f"BedGraph file not found: {args.bedgraph}")
        
        try:
            logger.info(f"Loading bedGraph track from {args.bedgraph}")
            from metainformant.epigenome import read_bedgraph
            
            track_df = read_bedgraph(args.bedgraph)
            logger.info(f"Loaded {len(track_df)} track regions")
            
            # Save summary
            track_summary = {
                "n_regions": len(track_df),
                "chromosomes": track_df['chromosome'].unique().tolist() if 'chromosome' in track_df.columns else [],
                "mean_value": float(track_df['value'].mean()) if 'value' in track_df.columns else None,
            }
            
            output_file = output_dir / "bedgraph_analysis.json"
            dump_json(track_summary, output_file)
            workflow_results["analyses"]["bedgraph"] = track_summary
            logger.info(f"BedGraph analysis saved to {output_file}")
        except Exception as e:
            logger.error(f"BedGraph analysis failed: {e}", exc_info=True)
            workflow_results["analyses"]["bedgraph"] = {"error": str(e)}
    
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
