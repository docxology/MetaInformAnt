#!/usr/bin/env python3
"""Machine learning simulation script.

This script generates synthetic feature matrices, labels, and test data
for classification and regression tasks.

Usage:
    python3 scripts/simulation/simulate_ml.py --type classification --n-samples 100 --n-features 50
    python3 scripts/simulation/simulate_ml.py --type regression --n-samples 200 --n-features 30
    python3 scripts/simulation/simulate_ml.py --type features --n-samples 150 --n-features 100
"""

import argparse
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core import io, logging, paths, validation

logger = logging.get_logger(__name__)


def simulate_classification(
    output_dir: Path,
    n_samples: int,
    n_features: int,
    n_classes: int,
    seed: int,
) -> dict:
    """Simulate classification dataset."""
    logger.info(f"Generating classification data: {n_samples} samples, {n_features} features, {n_classes} classes")
    rng = random.Random(seed)
    np.random.seed(seed)

    # Generate features with class-specific patterns
    X = []
    y = []

    samples_per_class = n_samples // n_classes

    for class_idx in range(n_classes):
        # Class-specific mean
        class_mean = np.random.randn(n_features) * 2 + class_idx * 3

        for _ in range(samples_per_class):
            # Generate sample from class distribution
            sample = class_mean + np.random.randn(n_features)
            X.append(sample)
            y.append(class_idx)

    # Add remaining samples
    for _ in range(n_samples - len(X)):
        class_idx = rng.randint(0, n_classes - 1)
        class_mean = np.random.randn(n_features) * 2 + class_idx * 3
        sample = class_mean + np.random.randn(n_features)
        X.append(sample)
        y.append(class_idx)

    # Create DataFrames
    feature_names = [f"feature_{i:04d}" for i in range(n_features)]
    X_df = pd.DataFrame(X, columns=feature_names)
    X_df.insert(0, "sample_id", [f"sample_{i:04d}" for i in range(n_samples)])

    y_df = pd.DataFrame({"sample_id": [f"sample_{i:04d}" for i in range(n_samples)], "label": y})

    # Save
    X_file = output_dir / "features.csv"
    y_file = output_dir / "labels.csv"
    X_df.to_csv(X_file, index=False)
    y_df.to_csv(y_file, index=False)

    logger.info(f"Classification data saved to {X_file} and {y_file}")

    return {
        "type": "classification",
        "n_samples": n_samples,
        "n_features": n_features,
        "n_classes": n_classes,
        "features_file": str(X_file),
        "labels_file": str(y_file),
    }


def simulate_regression(
    output_dir: Path,
    n_samples: int,
    n_features: int,
    noise_level: float,
    seed: int,
) -> dict:
    """Simulate regression dataset."""
    logger.info(f"Generating regression data: {n_samples} samples, {n_features} features")
    rng = random.Random(seed)
    np.random.seed(seed)

    # Generate features
    X = np.random.randn(n_samples, n_features)

    # Generate true coefficients
    n_relevant = n_features // 2
    coefficients = np.zeros(n_features)
    relevant_features = rng.sample(range(n_features), n_relevant)
    for idx in relevant_features:
        coefficients[idx] = np.random.randn()

    # Generate target
    y = X @ coefficients + np.random.normal(0, noise_level, n_samples)

    # Create DataFrames
    feature_names = [f"feature_{i:04d}" for i in range(n_features)]
    X_df = pd.DataFrame(X, columns=feature_names)
    X_df.insert(0, "sample_id", [f"sample_{i:04d}" for i in range(n_samples)])

    y_df = pd.DataFrame(
        {
            "sample_id": [f"sample_{i:04d}" for i in range(n_samples)],
            "target": y,
        }
    )

    # Save
    X_file = output_dir / "features.csv"
    y_file = output_dir / "targets.csv"
    X_df.to_csv(X_file, index=False)
    y_df.to_csv(y_file, index=False)

    # Save true coefficients
    coef_df = pd.DataFrame(
        {
            "feature": feature_names,
            "coefficient": coefficients,
            "is_relevant": [i in relevant_features for i in range(n_features)],
        }
    )
    coef_file = output_dir / "true_coefficients.csv"
    coef_df.to_csv(coef_file, index=False)

    logger.info(f"Regression data saved to {X_file} and {y_file}")

    return {
        "type": "regression",
        "n_samples": n_samples,
        "n_features": n_features,
        "n_relevant": n_relevant,
        "features_file": str(X_file),
        "targets_file": str(y_file),
        "coefficients_file": str(coef_file),
    }


def simulate_features(
    output_dir: Path,
    n_samples: int,
    n_features: int,
    feature_types: list[str],
    seed: int,
) -> dict:
    """Simulate feature matrix with different feature types."""
    logger.info(f"Generating feature matrix: {n_samples} samples, {n_features} features")
    rng = random.Random(seed)
    np.random.seed(seed)

    # Generate features
    features = {}

    for i in range(n_features):
        feature_type = rng.choice(feature_types) if feature_types else "continuous"

        if feature_type == "continuous":
            values = np.random.normal(0, 1, n_samples)
        elif feature_type == "binary":
            values = np.random.binomial(1, 0.5, n_samples)
        elif feature_type == "categorical":
            n_categories = rng.randint(3, 6)
            values = np.random.randint(0, n_categories, n_samples)
        else:  # count
            values = np.random.poisson(5, n_samples)

        features[f"feature_{i:04d}"] = values

    # Create DataFrame
    df = pd.DataFrame(features)
    df.insert(0, "sample_id", [f"sample_{i:04d}" for i in range(n_samples)])

    csv_file = output_dir / "features.csv"
    df.to_csv(csv_file, index=False)

    logger.info(f"Feature matrix saved to {csv_file}")

    return {
        "type": "features",
        "n_samples": n_samples,
        "n_features": n_features,
        "output_file": str(csv_file),
    }


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Machine learning simulation for testing and validation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Simulate classification data
  %(prog)s --type classification --n-samples 100 --n-features 50 --n-classes 3

  # Simulate regression data
  %(prog)s --type regression --n-samples 200 --n-features 30 --noise 0.1

  # Simulate feature matrix
  %(prog)s --type features --n-samples 150 --n-features 100
        """,
    )

    parser.add_argument(
        "--type",
        required=True,
        choices=["classification", "regression", "features"],
        help="Type of simulation to run",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/simulation/ml"),
        help="Output directory (default: output/simulation/ml)",
    )
    parser.add_argument("--n-samples", type=int, default=100, help="Number of samples")
    parser.add_argument("--n-features", type=int, default=50, help="Number of features")
    parser.add_argument("--n-classes", type=int, default=3, help="Number of classes (classification type)")
    parser.add_argument("--noise", type=float, default=0.1, help="Noise level (regression type)")
    parser.add_argument(
        "--feature-types",
        nargs="+",
        choices=["continuous", "binary", "categorical", "count"],
        default=["continuous"],
        help="Feature types (features type)",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")

    args = parser.parse_args()

    if args.verbose:
        import logging as std_logging

        logger.setLevel(std_logging.DEBUG)

    output_dir = paths.ensure_directory(args.output)

    try:
        if args.type == "classification":
            results = simulate_classification(output_dir, args.n_samples, args.n_features, args.n_classes, args.seed)
        elif args.type == "regression":
            results = simulate_regression(output_dir, args.n_samples, args.n_features, args.noise, args.seed)
        elif args.type == "features":
            results = simulate_features(output_dir, args.n_samples, args.n_features, args.feature_types, args.seed)

        # Save summary
        summary_file = output_dir / "simulation_summary.json"
        io.dump_json(results, summary_file, indent=2)
        logger.info(f"Simulation complete. Summary saved to {summary_file}")

        return 0
    except Exception as e:
        logger.error(f"Simulation failed: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
