#!/usr/bin/env python3
"""Machine Learning {{analysis_type}} Example

This example demonstrates:
- Loading and preprocessing biological data
- {{analysis_type.lower()}} model training
- Performance evaluation and visualization

Usage:
    python examples/ml/example_{{name}}.py

Expected output:
    output/examples/ml/{{name}}_results.json
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from metainformant.core import io
from metainformant.ml import classification


def main():
    """Main function demonstrating ML {{analysis_type.lower()}}."""
    print("Machine Learning {{analysis_type}} Example")
    print("=" * 40)

    # Create output directory
    output_dir = Path("output/examples/ml")
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "{{name}}_results.json"

    try:
        # Generate sample biological data
        np.random.seed(42)
        n_samples, n_features = 100, 20

        # Create synthetic features (e.g., gene expression levels)
        X = np.random.randn(n_samples, n_features)

        # Create synthetic binary labels (e.g., disease status)
        y = (X[:, 0] + X[:, 1] + np.random.randn(n_samples) * 0.5 > 0).astype(int)

        print(f"Generated dataset: {n_samples} samples, {n_features} features")
        print(f"Class distribution: {np.bincount(y)}")

        # Split data
        train_size = int(0.8 * n_samples)
        X_train, X_test = X[:train_size], X[train_size:]
        y_train, y_test = y[:train_size], y[train_size:]

        # Train model
        print("Training {{analysis_type.lower()}} model...")
        model = classification.train_classifier(X_train, y_train, method="rf")

        # Make predictions
        predictions = classification.predict_with_confidence(model, X_test)

        # Evaluate performance
        accuracy = np.mean(predictions[0] == y_test)
        print(".3f")

        # Save results
        results_data = {
            "dataset_info": {
                "n_samples": n_samples,
                "n_features": n_features,
                "train_samples": len(X_train),
                "test_samples": len(X_test),
            },
            "model_performance": {"accuracy": float(accuracy), "method": "rf"},
            "predictions": {
                "true_labels": y_test.tolist(),
                "predicted_labels": predictions[0].tolist(),
                "confidences": predictions[1].tolist(),
            },
        }

        results = {
            "example": "{{name}}",
            "domain": "ml",
            "analysis_type": "{{analysis_type}}",
            "description": "Machine learning {{analysis_type.lower()}} example",
            "results": results_data,
        }

        io.dump_json(results, output_file)
        print(f"✅ Results saved to: {output_file}")

    except Exception as e:
        print(f"❌ Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
