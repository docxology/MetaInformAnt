#!/usr/bin/env python3
"""Ecology analysis workflow orchestrator.

This script provides comprehensive orchestration for ecology analysis workflows,
including diversity indices, community composition, beta diversity, ordination,
indicator species analysis, functional ecology, and macroecology.

Usage:
    python3 scripts/ecology/run_ecology_analysis.py --input abundance.tsv --output output/ecology/results
    python3 scripts/ecology/run_ecology_analysis.py --input abundance.tsv --diversity --beta-diversity
    python3 scripts/ecology/run_ecology_analysis.py --input abundance.tsv --all
    python3 scripts/ecology/run_ecology_analysis.py --help
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Any

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core.io import dump_json, ensure_directory, read_csv
from metainformant.core.utils.logging import get_logger

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
  %(prog)s --input abundance.tsv --all

  # Beta diversity with ordination
  %(prog)s --input abundance.tsv --beta-diversity --ordination

  # Indicator species analysis with groups
  %(prog)s --input abundance.tsv --indicators --groups groups.csv

  # Macroecological analysis
  %(prog)s --input abundance.tsv --macroecology
        """,
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
        "--groups",
        type=Path,
        default=None,
        help="Optional group labels file (CSV with site,group columns) for indicator/ANOSIM analyses",
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
        "--ordination",
        action="store_true",
        help="Run ordination analysis (PCoA)",
    )
    parser.add_argument(
        "--indicators",
        action="store_true",
        help="Run indicator species analysis (requires --groups)",
    )
    parser.add_argument(
        "--macroecology",
        action="store_true",
        help="Fit species abundance distribution models",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Run all available analyses",
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

    # If --all, enable everything
    if args.all:
        args.diversity = True
        args.beta_diversity = True
        args.rarefaction = True
        args.ordination = True
        args.macroecology = True
        if args.groups:
            args.indicators = True

    if args.dry_run:
        logger.info("DRY RUN - no changes will be made")
        analyses = []
        if args.diversity:
            analyses.append("diversity")
        if args.beta_diversity:
            analyses.append("beta_diversity")
        if args.rarefaction:
            analyses.append("rarefaction")
        if args.ordination:
            analyses.append("ordination")
        if args.indicators:
            analyses.append("indicators")
        if args.macroecology:
            analyses.append("macroecology")
        logger.info(f"Would run: {', '.join(analyses) or 'none specified'}")
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

    # Load group labels if provided
    group_labels = None
    if args.groups and args.groups.exists():
        try:
            groups_df = read_csv(args.groups)
            group_labels = groups_df.iloc[:, 1].tolist()
            logger.info(f"Loaded group labels: {len(set(group_labels))} groups")
        except Exception as e:
            logger.warning(f"Could not load group labels: {e}")

    workflow_results: dict[str, Any] = {
        "input_file": str(args.input),
        "output_dir": str(output_dir),
        "data_shape": list(df.shape),
        "analyses": {},
    }

    # Convert to site abundances (each column is a site)
    site_abundances: list[list[float]] = []
    site_names: list[str] = []
    for col in df.columns:
        abundances = df[col].values.tolist()
        site_abundances.append(abundances)
        site_names.append(col)

    # Diversity analysis
    if args.diversity:
        try:
            logger.info("Calculating diversity indices...")
            from metainformant.ecology.analysis.community import (
                chao1_estimator,
                pielou_evenness,
                shannon_diversity,
                simpson_diversity,
                species_richness,
            )

            diversity_results = {}
            for site_name, abundances in zip(site_names, site_abundances):
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
            from metainformant.ecology.analysis.community import (
                alpha_beta_gamma_diversity,
                beta_diversity,
                community_similarity_matrix,
            )

            # Calculate pairwise dissimilarities
            dissimilarities = {}
            for i, site1_name in enumerate(site_names):
                for j, site2_name in enumerate(site_names[i + 1 :], start=i + 1):
                    site1_ab = [float(a) for a in site_abundances[i]]
                    site2_ab = [float(a) for a in site_abundances[j]]

                    pair_key = f"{site1_name}_vs_{site2_name}"
                    dissimilarities[pair_key] = {
                        "bray_curtis": float(beta_diversity(site1_ab, site2_ab, "bray_curtis")),
                        "jaccard": float(beta_diversity(site1_ab, site2_ab, "jaccard")),
                        "sorensen": float(beta_diversity(site1_ab, site2_ab, "sorensen")),
                    }

            # Alpha-beta-gamma partitioning
            abg = alpha_beta_gamma_diversity(site_abundances)

            beta_results = {
                "pairwise_dissimilarities": dissimilarities,
                "alpha_beta_gamma": abg,
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
            from metainformant.ecology.analysis.community import rarefaction_curve

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

    # Ordination analysis
    if args.ordination:
        try:
            logger.info("Running ordination analysis...")
            from metainformant.ecology.analysis.ordination import distance_matrix, pcoa

            # Build distance matrix from site abundances
            dist = distance_matrix(site_abundances, method="bray_curtis")

            # Run PCoA
            pcoa_result = pcoa(dist, n_components=min(3, len(site_abundances) - 1))

            ordination_results = {
                "method": "PCoA",
                "coordinates": pcoa_result["coordinates"],
                "eigenvalues": pcoa_result["eigenvalues"],
                "variance_explained": pcoa_result["variance_explained"],
            }

            output_file = output_dir / "ordination_analysis.json"
            dump_json(ordination_results, output_file)
            workflow_results["analyses"]["ordination"] = ordination_results
            logger.info(f"Ordination analysis saved to {output_file}")
        except Exception as e:
            logger.error(f"Ordination analysis failed: {e}", exc_info=True)
            workflow_results["analyses"]["ordination"] = {"error": str(e)}

    # Indicator species analysis
    if args.indicators:
        if group_labels is None:
            logger.warning("Indicator species analysis requires --groups file, skipping")
        else:
            try:
                logger.info("Running indicator species analysis...")
                from metainformant.ecology.analysis.indicators import anosim, indval

                # Build distance matrix for ANOSIM
                from metainformant.ecology.analysis.ordination import distance_matrix

                dist = distance_matrix(site_abundances, method="bray_curtis")

                # IndVal
                indval_result = indval(site_abundances, group_labels)

                # ANOSIM
                anosim_result = anosim(dist, group_labels, n_permutations=999)

                indicator_results = {
                    "indval": indval_result,
                    "anosim": anosim_result,
                }

                output_file = output_dir / "indicator_analysis.json"
                dump_json(indicator_results, output_file)
                workflow_results["analyses"]["indicators"] = indicator_results
                logger.info(f"Indicator analysis saved to {output_file}")
            except Exception as e:
                logger.error(f"Indicator analysis failed: {e}", exc_info=True)
                workflow_results["analyses"]["indicators"] = {"error": str(e)}

    # Macroecology analysis
    if args.macroecology:
        try:
            logger.info("Fitting species abundance distribution models...")
            from metainformant.ecology.analysis.macroecology import compare_sad_models

            # Pool abundances across sites for SAD fitting
            pooled = [0.0] * max(len(s) for s in site_abundances)
            for site_ab in site_abundances:
                for i, a in enumerate(site_ab):
                    pooled[i] += float(a)

            pooled_positive = [a for a in pooled if a > 0]
            sad_results = compare_sad_models(pooled_positive)

            output_file = output_dir / "macroecology_analysis.json"
            dump_json(sad_results, output_file)
            workflow_results["analyses"]["macroecology"] = sad_results
            logger.info(f"Macroecology analysis saved to {output_file}")
        except Exception as e:
            logger.error(f"Macroecology analysis failed: {e}", exc_info=True)
            workflow_results["analyses"]["macroecology"] = {"error": str(e)}

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
