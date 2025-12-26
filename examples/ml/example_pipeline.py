#!/usr/bin/env python3
"""Machine learning pipeline example.

This example demonstrates machine learning workflows for biological data using METAINFORMANT's ML toolkit.

Usage:
    python examples/ml/example_pipeline.py

Output:
    output/examples/ml/pipeline_results.json
"""

from __future__ import annotations

import numpy as np
from pathlib import Path
from metainformant.core import io
from metainformant.ml.classification import train_ensemble_classifier, cross_validate_biological

def main():
    """Demonstrate ML pipeline."""
    # Setup output directory
    output_dir = Path("output/examples/ml")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== METAINFORMANT ML Pipeline Example ===")

    # Create simulated biological data
    np.random.seed(42)
    n_samples, n_features = 200, 50
    X = np.random.randn(n_samples, n_features)
    y = (X[:, 0] + X[:, 5] + np.random.randn(n_samples) * 0.5) > 0  # Binary classification

    print(f"✓ Created dataset: {n_samples} samples, {n_features} features")

    # Train classifier
    print("\nTraining classifier...")
    models = train_ensemble_classifier(X, y, algorithms=["random_forest"])
    model = models["random_forest"]

    # Cross-validation
    cv_results = cross_validate_biological(X, y, algorithm="random_forest", cv_folds=5)
    accuracy = cv_results.get('accuracy', cv_results.get('mean_accuracy', 0.0))
    print(f"Cross-validation accuracy: {accuracy:.3f}")

    # Save results
    results_file = output_dir / "pipeline_results.json"
    # Convert numpy arrays and scalars to Python types for JSON serialization
    json_cv_results = {}
    for key, value in cv_results.items():
        if isinstance(value, np.ndarray):
            json_cv_results[key] = value.tolist()
        elif isinstance(value, (np.integer, np.floating, np.float64, np.float32)):
            json_cv_results[key] = float(value)
        elif isinstance(value, list):
            # Convert any numpy arrays or scalars in lists
            json_cv_results[key] = []
            for item in value:
                if isinstance(item, np.ndarray):
                    json_cv_results[key].append(item.tolist())
                elif isinstance(item, (np.integer, np.floating)):
                    json_cv_results[key].append(float(item))
                else:
                    json_cv_results[key].append(item)
        else:
            json_cv_results[key] = value

    io.dump_json({
        "ml_pipeline": {
            "dataset_shape": [n_samples, n_features],
            "model_type": "random_forest",
            "cv_results": json_cv_results
        }
    }, results_file, indent=2)

    print(f"✓ Results saved to: {results_file}")
    print("\n=== ML Pipeline Example Complete ===")

if __name__ == "__main__":
    main()
