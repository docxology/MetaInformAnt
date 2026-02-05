#!/usr/bin/env python3
"""Train prediction models on event sequences.

This script trains models to predict outcomes from event sequences.

Usage:
    python3 scripts/life_events/train_model.py --sequences data/sequences.json --outcomes data/outcomes.json --output output/life_events/model.json
    python3 scripts/life_events/train_model.py --sequences data/sequences.json --outcomes data/outcomes.json --model-type lstm --epochs 20
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
    GRUSequenceModel,
    LSTMSequenceModel,
    convert_sequences_to_tokens,
    load_sequences_from_json,
)

logger = setup_logger(__name__, level="INFO")


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Train prediction models on event sequences",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--sequences",
        type=Path,
        required=True,
        help="Input sequences file (JSON format)",
    )
    parser.add_argument(
        "--outcomes",
        type=Path,
        required=True,
        help="Input outcomes file (JSON format)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/life_events/model.json"),
        help="Output model file (default: output/life_events/model.json)",
    )
    parser.add_argument(
        "--model-type",
        type=str,
        default="embedding",
        choices=["embedding", "simple", "lstm", "gru"],
        help="Model type (default: embedding)",
    )
    parser.add_argument(
        "--task-type",
        type=str,
        default="classification",
        choices=["classification", "regression"],
        help="Task type (default: classification)",
    )
    parser.add_argument(
        "--embedding-dim",
        type=int,
        default=100,
        help="Embedding dimension (default: 100)",
    )
    parser.add_argument(
        "--hidden-dim",
        type=int,
        default=64,
        help="Hidden dimension for LSTM/GRU (default: 64)",
    )
    parser.add_argument(
        "--epochs",
        type=int,
        default=10,
        help="Training epochs (default: 10)",
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=42,
        help="Random seed (default: 42)",
    )
    return parser.parse_args()


def main():
    """Main function."""
    args = parse_args()

    if not args.sequences.exists():
        logger.error(f"Sequences file not found: {args.sequences}")
        return 1

    if not args.outcomes.exists():
        logger.error(f"Outcomes file not found: {args.outcomes}")
        return 1

    logger.info(f"Loading sequences from {args.sequences}")
    sequences = load_sequences_from_json(args.sequences)
    logger.info(f"✅ Loaded {len(sequences)} sequences")

    logger.info(f"Loading outcomes from {args.outcomes}")
    outcomes_data = io.load_json(args.outcomes)
    outcomes = np.array(outcomes_data.get("outcomes", outcomes_data.get("outcome", [])))
    logger.info(f"✅ Loaded outcomes for {len(outcomes)} sequences")

    if len(outcomes) != len(sequences):
        logger.error(f"Mismatch: {len(sequences)} sequences but {len(outcomes)} outcomes")
        return 1

    # Convert to tokens
    sequences_tokens = convert_sequences_to_tokens(sequences)

    logger.info(f"Training {args.model_type} model...")

    if args.model_type == "lstm":
        model = LSTMSequenceModel(
            embedding_dim=args.embedding_dim,
            hidden_dim=args.hidden_dim,
            task_type=args.task_type,
            epochs=args.epochs,
            random_state=args.random_state,
        )
    elif args.model_type == "gru":
        model = GRUSequenceModel(
            embedding_dim=args.embedding_dim,
            hidden_dim=args.hidden_dim,
            task_type=args.task_type,
            epochs=args.epochs,
            random_state=args.random_state,
        )
    else:
        model = EventSequencePredictor(
            model_type=args.model_type,
            task_type=args.task_type,
            embedding_dim=args.embedding_dim,
            hidden_dim=args.hidden_dim,
            random_state=args.random_state,
        )

    model.fit(sequences_tokens, outcomes)

    # Save model
    paths.ensure_directory(args.output.parent)
    if hasattr(model, "save_model"):
        model.save_model(args.output)
    else:
        logger.warning("Model does not support saving, skipping...")

    logger.info(f"✅ Model trained and saved to {args.output}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
