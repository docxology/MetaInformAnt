#!/usr/bin/env python3
"""Learn event embeddings from sequences.

This script learns dense vector representations of events using Word2Vec-style methods.

Usage:
    python3 scripts/life_events/learn_embeddings.py --input data/sequences.json --output output/life_events/embeddings.json
    python3 scripts/life_events/learn_embeddings.py --input data/sequences.json --embedding-dim 200 --window-size 10
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core import io, paths
from metainformant.core.logging import setup_logger
from metainformant.life_events import (
    convert_sequences_to_tokens,
    learn_event_embeddings,
    load_sequences_from_json,
)

logger = setup_logger(__name__, level="INFO")


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Learn event embeddings from sequences",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Input sequences file (JSON format)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/life_events/embeddings.json"),
        help="Output embeddings file (default: output/life_events/embeddings.json)",
    )
    parser.add_argument(
        "--embedding-dim",
        type=int,
        default=100,
        help="Embedding dimension (default: 100)",
    )
    parser.add_argument(
        "--window-size",
        type=int,
        default=5,
        help="Context window size (default: 5)",
    )
    parser.add_argument(
        "--method",
        type=str,
        default="skipgram",
        choices=["skipgram", "cbow"],
        help="Embedding method (default: skipgram)",
    )
    parser.add_argument(
        "--epochs",
        type=int,
        default=10,
        help="Number of training epochs (default: 10)",
    )
    parser.add_argument(
        "--learning-rate",
        type=float,
        default=0.01,
        help="Learning rate (default: 0.01)",
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=42,
        help="Random seed (default: 42)",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose logging",
    )
    return parser.parse_args()


def main():
    """Main function."""
    args = parse_args()
    
    if args.verbose:
        logger.setLevel("DEBUG")
    
    if not args.input.exists():
        logger.error(f"Input file not found: {args.input}")
        return 1
    
    logger.info(f"Loading sequences from {args.input}")
    sequences = load_sequences_from_json(args.input)
    logger.info(f"✅ Loaded {len(sequences)} sequences")
    
    # Convert to tokens
    sequences_tokens = convert_sequences_to_tokens(sequences)
    
    logger.info("Learning event embeddings...")
    embeddings = learn_event_embeddings(
        sequences_tokens,
        embedding_dim=args.embedding_dim,
        window_size=args.window_size,
        method=args.method,
        epochs=args.epochs,
        learning_rate=args.learning_rate,
        random_state=args.random_state,
        verbose=args.verbose,
    )
    
    # Save embeddings
    paths.ensure_directory(args.output.parent)
    embeddings_dict = {k: v.tolist() for k, v in embeddings.items()}
    io.dump_json(embeddings_dict, args.output, indent=2)
    
    logger.info(f"✅ Learned embeddings for {len(embeddings)} events")
    logger.info(f"✅ Saved embeddings to {args.output}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

