#!/usr/bin/env python3
"""Predict outcomes from event sequences using trained model.

This script loads a trained model and makes predictions on new sequences.

Usage:
    python3 scripts/life_events/predict_outcomes.py --model output/life_events/model.json --sequences data/test_sequences.json --output output/life_events/predictions.json
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core import io, paths
from metainformant.core.utils.logging import setup_logger
from metainformant.life_events import (
    EventSequencePredictor,
    convert_sequences_to_tokens,
    load_sequences_from_json,
)

logger = setup_logger(__name__, level="INFO")


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Predict outcomes from event sequences",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--model",
        type=Path,
        required=True,
        help="Trained model file (JSON format)",
    )
    parser.add_argument(
        "--sequences",
        type=Path,
        required=True,
        help="Input sequences file (JSON format)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/life_events/predictions.json"),
        help="Output predictions file (default: output/life_events/predictions.json)",
    )
    parser.add_argument(
        "--probabilities",
        action="store_true",
        help="Include class probabilities (for classification)",
    )
    return parser.parse_args()


def main():
    """Main function."""
    args = parse_args()
    
    if not args.model.exists():
        logger.error(f"Model file not found: {args.model}")
        return 1
    
    if not args.sequences.exists():
        logger.error(f"Sequences file not found: {args.sequences}")
        return 1
    
    logger.info(f"Loading model from {args.model}")
    predictor = EventSequencePredictor.load_model(args.model)
    logger.info(f"✅ Loaded model (type: {predictor.model_type}, task: {predictor.task_type})")
    
    logger.info(f"Loading sequences from {args.sequences}")
    sequences = load_sequences_from_json(args.sequences)
    logger.info(f"✅ Loaded {len(sequences)} sequences")
    
    sequences_tokens = convert_sequences_to_tokens(sequences)
    
    logger.info("Making predictions...")
    predictions = predictor.predict(sequences_tokens)
    
    results = {
        "predictions": predictions.tolist() if hasattr(predictions, "tolist") else list(predictions),
        "n_sequences": len(sequences),
        "sequence_ids": [seq.person_id for seq in sequences],
    }
    
    if args.probabilities and predictor.task_type == "classification":
        try:
            probabilities = predictor.predict_proba(sequences_tokens)
            results["probabilities"] = probabilities.tolist()
        except Exception as e:
            logger.warning(f"Could not compute probabilities: {e}")
    
    # Save predictions
    paths.ensure_directory(args.output.parent)
    io.dump_json(results, args.output, indent=2)
    
    logger.info(f"✅ Saved predictions to {args.output}")
    
    # Print summary
    if predictor.task_type == "classification":
        unique, counts = np.unique(predictions, return_counts=True)
        logger.info("Prediction summary:")
        for label, count in zip(unique, counts):
            logger.info(f"  Class {label}: {count} sequences")
    else:
        logger.info(f"Prediction range: [{predictions.min():.3f}, {predictions.max():.3f}]")
        logger.info(f"Mean prediction: {predictions.mean():.3f}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

