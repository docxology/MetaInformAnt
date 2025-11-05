#!/usr/bin/env python3
"""Interpret model predictions.

This script analyzes model predictions to identify important events and patterns.

Usage:
    python3 scripts/life_events/interpret_predictions.py --model output/life_events/model.json --sequences data/sequences.json --output output/life_events/interpretation/
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core import io, paths
from metainformant.core.logging import setup_logger
from metainformant.life_events import (
    EventSequencePredictor,
    convert_sequences_to_tokens,
    event_importance,
    learn_event_embeddings,
    load_sequences_from_json,
    plot_prediction_importance,
)

logger = setup_logger(__name__, level="INFO")


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Interpret model predictions",
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
        default=Path("output/life_events/interpretation"),
        help="Output directory (default: output/life_events/interpretation)",
    )
    parser.add_argument(
        "--method",
        type=str,
        default="permutation",
        choices=["permutation", "gradient"],
        help="Importance computation method (default: permutation)",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=20,
        help="Number of top events to display (default: 20)",
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
    logger.info(f"✅ Loaded model")
    
    logger.info(f"Loading sequences from {args.sequences}")
    sequences = load_sequences_from_json(args.sequences)
    logger.info(f"✅ Loaded {len(sequences)} sequences")
    
    sequences_tokens = convert_sequences_to_tokens(sequences)
    
    # Get embeddings
    if hasattr(predictor, "event_embeddings") and predictor.event_embeddings:
        embeddings = predictor.event_embeddings
    else:
        logger.info("Learning embeddings for interpretation...")
        embeddings = learn_event_embeddings(sequences_tokens, embedding_dim=100)
    
    logger.info("Computing event importance...")
    importance = event_importance(
        predictor,
        sequences_tokens,
        embeddings,
        method=args.method,
    )
    
    paths.ensure_directory(args.output)
    
    # Save importance scores
    importance_file = args.output / "event_importance.json"
    io.dump_json(importance, importance_file, indent=2)
    logger.info(f"✅ Saved importance scores to {importance_file}")
    
    # Generate visualization
    logger.info("Generating importance visualization...")
    plot_prediction_importance(
        importance,
        top_n=args.top_n,
        output_path=args.output / "prediction_importance.png",
    )
    logger.info("✅ Importance visualization generated")
    
    # Print top events
    sorted_events = sorted(importance.items(), key=lambda x: x[1], reverse=True)
    logger.info(f"Top {args.top_n} most important events:")
    for i, (event, score) in enumerate(sorted_events[:args.top_n], 1):
        logger.info(f"  {i}. {event}: {score:.4f}")
    
    logger.info(f"✅ Interpretation complete. Results saved to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

