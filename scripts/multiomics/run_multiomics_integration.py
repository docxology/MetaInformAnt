#!/usr/bin/env python3
"""Multi-omics integration workflow orchestrator.

This script provides comprehensive orchestration for multi-omics integration workflows,
including data loading, integration, joint analysis, and correlation analysis.

Usage:
    python3 scripts/multiomics/run_multiomics_integration.py --genomics genomics.csv --transcriptomics expression.tsv --output output/multiomics/results
    python3 scripts/multiomics/run_multiomics_integration.py --genomics g.csv --transcriptomics t.tsv --joint-pca --n-components 50
    python3 scripts/multiomics/run_multiomics_integration.py --help
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Any

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core.io import dump_json, ensure_directory
from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Multi-omics integration workflow orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic integration with two omics layers
  %(prog)s --genomics genomics.csv --transcriptomics expression.tsv --output output/multiomics/basic

  # Full analysis with joint PCA
  %(prog)s --genomics g.csv --transcriptomics t.tsv --proteomics p.csv --joint-pca --n-components 50

  # Canonical correlation analysis
  %(prog)s --genomics g.csv --transcriptomics t.tsv --canonical-correlation --n-components 10
        """,
    )
    parser.add_argument(
        "--genomics",
        type=Path,
        help="Genomics data file (CSV/TSV)",
    )
    parser.add_argument(
        "--transcriptomics",
        type=Path,
        help="Transcriptomics data file (CSV/TSV)",
    )
    parser.add_argument(
        "--proteomics",
        type=Path,
        help="Proteomics data file (CSV/TSV)",
    )
    parser.add_argument(
        "--metabolomics",
        type=Path,
        help="Metabolomics data file (CSV/TSV)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/multiomics"),
        help="Output directory (default: output/multiomics)",
    )
    parser.add_argument(
        "--joint-pca",
        action="store_true",
        help="Perform joint PCA analysis",
    )
    parser.add_argument(
        "--joint-nmf",
        action="store_true",
        help="Perform joint NMF analysis",
    )
    parser.add_argument(
        "--canonical-correlation",
        action="store_true",
        help="Perform canonical correlation analysis",
    )
    parser.add_argument(
        "--n-components",
        type=int,
        default=20,
        help="Number of components for dimensionality reduction (default: 20)",
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
    """Execute multi-omics integration workflow.

    Args:
        args: Parsed command-line arguments
    """
    logger.info("Starting multi-omics integration workflow")
    logger.info(f"Output: {args.output}")

    # Collect input files
    data_dict = {}
    if args.genomics:
        data_dict["genomics"] = args.genomics
    if args.transcriptomics:
        data_dict["transcriptomics"] = args.transcriptomics
    if args.proteomics:
        data_dict["proteomics"] = args.proteomics
    if args.metabolomics:
        data_dict["metabolomics"] = args.metabolomics

    if not data_dict:
        raise ValueError("At least one omics data file must be provided")

    if args.dry_run:
        logger.info("DRY RUN - no changes will be made")
        logger.info(f"Would integrate: {list(data_dict.keys())}")
        logger.info(f"Would write to: {args.output}")
        if args.joint_pca:
            logger.info(f"Would perform joint PCA with {args.n_components} components")
        if args.joint_nmf:
            logger.info(f"Would perform joint NMF with {args.n_components} components")
        if args.canonical_correlation:
            logger.info(f"Would perform CCA with {args.n_components} components")
        return

    # Ensure output directory exists
    output_dir = ensure_directory(args.output)
    logger.info(f"Output directory: {output_dir}")

    # Integrate omics data
    logger.info("Integrating omics data...")
    from metainformant.multiomics import integrate_omics_data

    try:
        omics_data = integrate_omics_data(data_dict)
        logger.info(
            f"Integrated {omics_data.n_samples} samples across {len(omics_data.layer_names)} layers: {omics_data.layer_names}"
        )
    except Exception as e:
        logger.error(f"Failed to integrate omics data: {e}")
        raise

    # Run analyses
    workflow_results = {
        "output_dir": str(output_dir),
        "layers": omics_data.layer_names,
        "n_samples": omics_data.n_samples,
        "analyses": {},
    }

    # Joint PCA analysis
    if args.joint_pca:
        try:
            logger.info(f"Performing joint PCA with {args.n_components} components...")
            from metainformant.multiomics import joint_pca

            embeddings, loadings, variance = joint_pca(omics_data, n_components=args.n_components, standardize=True)

            results = {
                "embeddings_shape": list(embeddings.shape),
                "explained_variance": variance[:10].tolist() if len(variance) >= 10 else variance.tolist(),
                "total_variance_explained": float(sum(variance[: args.n_components])),
            }

            # Save embeddings
            import pandas as pd

            embeddings_df = pd.DataFrame(embeddings, index=omics_data.samples)
            embeddings_file = output_dir / "joint_pca_embeddings.csv"
            embeddings_df.to_csv(embeddings_file)
            logger.info(f"Joint PCA embeddings saved to {embeddings_file}")

            output_file = output_dir / "joint_pca_results.json"
            dump_json(results, output_file)
            workflow_results["analyses"]["joint_pca"] = results
        except Exception as e:
            logger.error(f"Joint PCA failed: {e}", exc_info=True)
            workflow_results["analyses"]["joint_pca"] = {"error": str(e)}

    # Joint NMF analysis
    if args.joint_nmf:
        try:
            logger.info(f"Performing joint NMF with {args.n_components} components...")
            from metainformant.multiomics import joint_nmf

            W, H = joint_nmf(omics_data, n_components=args.n_components, max_iter=200, random_state=42)

            results = {
                "sample_factors_shape": list(W.shape),
                "feature_factors": {layer: list(H[layer].shape) for layer in omics_data.layer_names if layer in H},
            }

            # Save factors
            import pandas as pd

            W_df = pd.DataFrame(W, index=omics_data.samples)
            W_file = output_dir / "joint_nmf_sample_factors.csv"
            W_df.to_csv(W_file)
            logger.info(f"Joint NMF sample factors saved to {W_file}")

            output_file = output_dir / "joint_nmf_results.json"
            dump_json(results, output_file)
            workflow_results["analyses"]["joint_nmf"] = results
        except Exception as e:
            logger.error(f"Joint NMF failed: {e}", exc_info=True)
            workflow_results["analyses"]["joint_nmf"] = {"error": str(e)}

    # Canonical correlation analysis (requires at least 2 layers)
    if args.canonical_correlation:
        if len(omics_data.layer_names) >= 2:
            try:
                logger.info(f"Performing canonical correlation analysis with {args.n_components} components...")
                from metainformant.multiomics import canonical_correlation

                layer1, layer2 = omics_data.layer_names[0], omics_data.layer_names[1]
                X_c, Y_c, X_w, Y_w, correlations = canonical_correlation(
                    omics_data, layer_pair=(layer1, layer2), n_components=args.n_components, regularization=0.01
                )

                results = {
                    "layer_pair": [layer1, layer2],
                    "canonical_correlations": (
                        correlations[:10].tolist() if len(correlations) >= 10 else correlations.tolist()
                    ),
                    "max_correlation": float(max(correlations)) if len(correlations) > 0 else 0.0,
                }

                output_file = output_dir / "canonical_correlation_results.json"
                dump_json(results, output_file)
                workflow_results["analyses"]["canonical_correlation"] = results
                logger.info(
                    f"Canonical correlation analysis complete (max correlation: {results['max_correlation']:.3f})"
                )
            except Exception as e:
                logger.error(f"Canonical correlation failed: {e}", exc_info=True)
                workflow_results["analyses"]["canonical_correlation"] = {"error": str(e)}
        else:
            logger.warning("Canonical correlation requires at least 2 omics layers")
            workflow_results["analyses"]["canonical_correlation"] = {"error": "Need at least 2 layers"}

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
