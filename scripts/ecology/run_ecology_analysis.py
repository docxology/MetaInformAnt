#!/usr/bin/env python3
"""Ecology analysis workflow orchestrator.

This script provides comprehensive orchestration for ecology analysis workflows,
including diversity indices, community composition, and beta diversity analysis.

Usage:
    python3 scripts/ecology/run_ecology_analysis.py --input abundance.tsv --output output/ecology/results
    python3 scripts/ecology/run_ecology_analysis.py --input abundance.tsv --diversity --beta-diversity
    python3 scripts/ecology/run_ecology_analysis.py --help
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Any

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core.io import ensure_directory, dump_json, read_csv
from metainformant.core.logging import get_logger

logger = get_logger(__name__)


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Ecology analysis workflow orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Diversity analysis
  %(prog)s --input abundance.tsv --diversity --output output/ecology/diversity

  # Full community analysis
  %(prog)s --input abundance.tsv --diversity --beta-diversity --rarefaction

  # Beta diversity only
  %(prog)s --input abundance.tsv --beta-diversity
        """
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Input abundance table (CSV/TSV format, species x sites)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/ecology"),
        help="Output directory (default: output/ecology)",
    )
    parser.add_argument(
        "--diversity",
        action="store_true",
        help="Calculate diversity indices (Shannon, Simpson, etc.)",
    )
    parser.add_argument(
        "--beta-diversity",
        action="store_true",
        help="Calculate beta diversity metrics",
    )
    parser.add_argument(
        "--rarefaction",
        action="store_true",
        help="Generate rarefaction curves",
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
    """Execute ecology analysis workflow."""
    logger.info("Starting ecology analysis workflow")
    logger.info(f"Input: {args.input}")
    logger.info(f"Output: {args.output}")
    
    if not args.input.exists():
        raise FileNotFoundError(f"Input file not found: {args.input}")
    
    if args.dry_run:
        logger.info("DRY RUN - no changes will be made")
        return
    
    output_dir = ensure_directory(args.output)
    logger.info(f"Output directory: {output_dir}")
    
    # Load abundance data
    logger.info(f"Loading abundance data from {args.input}")
    try:
        df = read_csv(args.input, sep="\t" if args.input.suffix == ".tsv" else ",")
        logger.info(f"Loaded abundance table: {df.shape[0]} species, {df.shape[1]} sites")
    except Exception as e:
        logger.error(f"Failed to load data: {e}")
        raise
    
    workflow_results = {
        "input_file": str(args.input),
        "output_dir": str(output_dir),
        "data_shape": list(df.shape),
        "analyses": {},
    }
    
    # Convert to site abundances (each column is a site)
    site_abundances = []
    site_names = []
    for col in df.columns:
        abundances = df[col].values.tolist()
        site_abundances.append(abundances)
        site_names.append(col)
    
    # Diversity analysis
    if args.diversity:
        try:
            logger.info("Calculating diversity indices...")
            from metainformant.ecology.community import (
                shannon_diversity,
                simpson_diversity,
                species_richness,
                pielou_evenness,
                chao1_estimator,
            )
            
            diversity_results = {}
            for i, (site_name, abundances) in enumerate(zip(site_names, site_abundances)):
                # Convert to float and filter zeros
                abundances_float = [float(a) for a in abundances if float(a) > 0]
                abundances_int = [int(a) for a in abundances if float(a) > 0]
                
                diversity_results[site_name] = {
                    "shannon_diversity": float(shannon_diversity(abundances_float)),
                    "simpson_diversity": float(simpson_diversity(abundances_float)),
                    "species_richness": int(species_richness(abundances_float)),
                    "pielou_evenness": float(pielou_evenness(abundances_float)),
                }
                
                if abundances_int:
                    diversity_results[site_name]["chao1_estimator"] = float(chao1_estimator(abundances_int))
            
            output_file = output_dir / "diversity_analysis.json"
            dump_json(diversity_results, output_file)
            workflow_results["analyses"]["diversity"] = diversity_results
            logger.info(f"Diversity analysis saved to {output_file}")
        except Exception as e:
            logger.error(f"Diversity analysis failed: {e}", exc_info=True)
            workflow_results["analyses"]["diversity"] = {"error": str(e)}
    
    # Beta diversity analysis
    if args.beta_diversity:
        try:
            logger.info("Calculating beta diversity...")
            from metainformant.ecology.community import (
                bray_curtis_dissimilarity,
                jaccard_similarity,
                beta_diversity_partitioning,
            )
            
            # Calculate pairwise dissimilarities
            dissimilarities = {}
            for i, site1_name in enumerate(site_names):
                for j, site2_name in enumerate(site_names[i+1:], start=i+1):
                    site1_ab = [float(a) for a in site_abundances[i]]
                    site2_ab = [float(a) for a in site_abundances[j]]
                    
                    pair_key = f"{site1_name}_vs_{site2_name}"
                    dissimilarities[pair_key] = {
                        "bray_curtis": float(bray_curtis_dissimilarity(site1_ab, site2_ab)),
                        "jaccard": float(jaccard_similarity(site1_ab, site2_ab)),
                    }
            
            # Beta diversity partitioning
            beta_partition = beta_diversity_partitioning(site_abundances)
            
            beta_results = {
                "pairwise_dissimilarities": dissimilarities,
                "beta_diversity_partitioning": beta_partition,
            }
            
            output_file = output_dir / "beta_diversity_analysis.json"
            dump_json(beta_results, output_file)
            workflow_results["analyses"]["beta_diversity"] = beta_results
            logger.info(f"Beta diversity analysis saved to {output_file}")
        except Exception as e:
            logger.error(f"Beta diversity analysis failed: {e}", exc_info=True)
            workflow_results["analyses"]["beta_diversity"] = {"error": str(e)}
    
    # Rarefaction curves
    if args.rarefaction:
        try:
            logger.info("Generating rarefaction curves...")
            from metainformant.ecology.community import rarefaction_curve
            
            rarefaction_results = {}
            for site_name, abundances in zip(site_names, site_abundances):
                abundances_int = [int(a) for a in abundances if float(a) > 0]
                if abundances_int:
                    curve = rarefaction_curve(abundances_int)
                    rarefaction_results[site_name] = curve[:100]  # Limit output size
            
            output_file = output_dir / "rarefaction_curves.json"
            dump_json(rarefaction_results, output_file)
            workflow_results["analyses"]["rarefaction"] = {
                "n_sites": len(rarefaction_results),
                "output_file": str(output_file),
            }
            logger.info(f"Rarefaction curves saved to {output_file}")
        except Exception as e:
            logger.error(f"Rarefaction analysis failed: {e}", exc_info=True)
            workflow_results["analyses"]["rarefaction"] = {"error": str(e)}
    
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
