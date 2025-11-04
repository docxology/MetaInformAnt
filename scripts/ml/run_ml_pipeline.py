#!/usr/bin/env python3
"""Machine learning workflow orchestrator.

This script provides comprehensive orchestration for machine learning workflows,
including feature selection, model training, cross-validation, and evaluation.

Usage:
    python3 scripts/ml/run_ml_pipeline.py --features features.csv --labels labels.csv --output output/ml/results
    python3 scripts/ml/run_ml_pipeline.py --features X.csv --labels y.csv --classify --feature-selection
    python3 scripts/ml/run_ml_pipeline.py --help
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
        description="Machine learning workflow orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Classification with feature selection
  %(prog)s --features X.csv --labels y.csv --classify --feature-selection --n-features 50

  # Regression with cross-validation
  %(prog)s --features X.csv --labels y.csv --regress --cross-validate --cv-folds 5

  # Dimensionality reduction
  %(prog)s --features X.csv --reduce-dimensions --method pca --n-components 20
        """
    )
    parser.add_argument(
        "--features",
        type=Path,
        required=True,
        help="Feature matrix file (CSV/TSV)",
    )
    parser.add_argument(
        "--labels",
        type=Path,
        help="Labels/targets file (CSV/TSV)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/ml"),
        help="Output directory (default: output/ml)",
    )
    parser.add_argument(
        "--classify",
        action="store_true",
        help="Perform classification",
    )
    parser.add_argument(
        "--regress",
        action="store_true",
        help="Perform regression",
    )
    parser.add_argument(
        "--feature-selection",
        action="store_true",
        help="Perform feature selection",
    )
    parser.add_argument(
        "--n-features",
        type=int,
        default=100,
        help="Number of features to select (default: 100)",
    )
    parser.add_argument(
        "--reduce-dimensions",
        action="store_true",
        help="Perform dimensionality reduction",
    )
    parser.add_argument(
        "--method",
        type=str,
        default="pca",
        choices=["pca", "umap", "tsne"],
        help="Dimensionality reduction method (default: pca)",
    )
    parser.add_argument(
        "--n-components",
        type=int,
        default=20,
        help="Number of components (default: 20)",
    )
    parser.add_argument(
        "--cross-validate",
        action="store_true",
        help="Perform cross-validation",
    )
    parser.add_argument(
        "--cv-folds",
        type=int,
        default=5,
        help="Number of CV folds (default: 5)",
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
    """Execute machine learning workflow."""
    logger.info("Starting machine learning workflow")
    logger.info(f"Features: {args.features}")
    logger.info(f"Output: {args.output}")
    
    if not args.features.exists():
        raise FileNotFoundError(f"Features file not found: {args.features}")
    
    if args.dry_run:
        logger.info("DRY RUN - no changes will be made")
        return
    
    output_dir = ensure_directory(args.output)
    logger.info(f"Output directory: {output_dir}")
    
    # Load features
    logger.info(f"Loading features from {args.features}")
    try:
        import pandas as pd
        import numpy as np
        X_df = read_csv(args.features)
        X = X_df.values
        logger.info(f"Loaded features: {X.shape[0]} samples, {X.shape[1]} features")
    except Exception as e:
        logger.error(f"Failed to load features: {e}")
        raise
    
    workflow_results = {
        "input_file": str(args.features),
        "output_dir": str(output_dir),
        "feature_shape": list(X.shape),
        "analyses": {},
    }
    
    # Load labels if provided
    y = None
    if args.labels:
        if not args.labels.exists():
            raise FileNotFoundError(f"Labels file not found: {args.labels}")
        logger.info(f"Loading labels from {args.labels}")
        try:
            y_df = read_csv(args.labels)
            y = y_df.values.flatten() if len(y_df.shape) > 1 else y_df.values
            logger.info(f"Loaded labels: {len(y)} samples")
            workflow_results["labels_file"] = str(args.labels)
        except Exception as e:
            logger.error(f"Failed to load labels: {e}")
            raise
    
    # Feature selection
    if args.feature_selection and y is not None:
        try:
            logger.info(f"Selecting {args.n_features} features...")
            from metainformant.ml import select_features_univariate
            
            X_selected, selected_indices = select_features_univariate(
                X, y, k=args.n_features, method="f_score"
            )
            logger.info(f"Selected {len(selected_indices)} features")
            X = X_selected
            workflow_results["analyses"]["feature_selection"] = {
                "n_selected": len(selected_indices),
                "selected_indices": selected_indices[:100],  # Limit output size
            }
        except Exception as e:
            logger.error(f"Feature selection failed: {e}", exc_info=True)
            workflow_results["analyses"]["feature_selection"] = {"error": str(e)}
    
    # Classification
    if args.classify and y is not None:
        try:
            logger.info("Training classifier...")
            from metainformant.ml import BiologicalClassifier, evaluate_classifier
            
            classifier = BiologicalClassifier(algorithm="random_forest", random_state=42)
            classifier.fit(X, y)
            
            predictions = classifier.predict(X)
            metrics = evaluate_classifier(classifier, X, y)
            
            workflow_results["analyses"]["classification"] = {
                "algorithm": "random_forest",
                "metrics": metrics,
            }
            
            if args.cross_validate:
                from metainformant.ml import cross_validate_biological
                cv_results = cross_validate_biological(classifier, X, y, cv=args.cv_folds)
                workflow_results["analyses"]["classification"]["cv_results"] = cv_results
                logger.info(f"Cross-validation accuracy: {cv_results.get('mean_score', 0):.3f}")
        except Exception as e:
            logger.error(f"Classification failed: {e}", exc_info=True)
            workflow_results["analyses"]["classification"] = {"error": str(e)}
    
    # Regression
    if args.regress and y is not None:
        try:
            logger.info("Training regressor...")
            from metainformant.ml import BiologicalRegressor, evaluate_regressor
            
            regressor = BiologicalRegressor(algorithm="linear", random_state=42)
            regressor.fit(X, y)
            
            predictions = regressor.predict(X)
            metrics = evaluate_regressor(regressor, X, y)
            
            workflow_results["analyses"]["regression"] = {
                "algorithm": "linear",
                "metrics": metrics,
            }
        except Exception as e:
            logger.error(f"Regression failed: {e}", exc_info=True)
            workflow_results["analyses"]["regression"] = {"error": str(e)}
    
    # Dimensionality reduction
    if args.reduce_dimensions:
        try:
            logger.info(f"Reducing dimensions using {args.method}...")
            if args.method == "pca":
                from metainformant.ml import reduce_dimensions_pca
                X_reduced = reduce_dimensions_pca(X, n_components=args.n_components)
            elif args.method == "umap":
                from metainformant.ml import reduce_dimensions_umap
                X_reduced = reduce_dimensions_umap(X, n_components=args.n_components)
            elif args.method == "tsne":
                from metainformant.ml import reduce_dimensions_tsne
                X_reduced = reduce_dimensions_tsne(X, n_components=args.n_components)
            
            # Save reduced features
            import pandas as pd
            X_reduced_df = pd.DataFrame(X_reduced)
            output_file = output_dir / f"{args.method}_reduced_features.csv"
            X_reduced_df.to_csv(output_file, index=False)
            logger.info(f"Reduced features saved to {output_file}")
            
            workflow_results["analyses"]["dimensionality_reduction"] = {
                "method": args.method,
                "n_components": args.n_components,
                "output_shape": list(X_reduced.shape),
            }
        except Exception as e:
            logger.error(f"Dimensionality reduction failed: {e}", exc_info=True)
            workflow_results["analyses"]["dimensionality_reduction"] = {"error": str(e)}
    
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




