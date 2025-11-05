#!/usr/bin/env python3
"""Export embeddings in various formats.

This script exports event embeddings to different file formats.

Usage:
    python3 scripts/life_events/export_embeddings.py --input output/life_events/embeddings.json --output output/life_events/embeddings.csv --format csv
    python3 scripts/life_events/export_embeddings.py --input output/life_events/embeddings.json --output output/life_events/embeddings.txt --format word2vec
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core import io, paths
from metainformant.core.logging import setup_logger

logger = setup_logger(__name__, level="INFO")


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Export embeddings in various formats",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Input embeddings file (JSON format)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output file path",
    )
    parser.add_argument(
        "--format",
        type=str,
        default="csv",
        choices=["csv", "word2vec", "numpy"],
        help="Output format (default: csv)",
    )
    return parser.parse_args()


def export_csv(embeddings: dict, output_path: Path):
    """Export embeddings to CSV format."""
    import csv
    
    tokens = sorted(embeddings.keys())
    embedding_dim = len(embeddings[tokens[0]])
    
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        # Header
        header = ['token'] + [f'dim_{i}' for i in range(embedding_dim)]
        writer.writerow(header)
        # Data
        for token in tokens:
            row = [token] + embeddings[token].tolist()
            writer.writerow(row)


def export_word2vec(embeddings: dict, output_path: Path):
    """Export embeddings to Word2Vec text format."""
    tokens = sorted(embeddings.keys())
    embedding_dim = len(embeddings[tokens[0]])
    
    with open(output_path, 'w') as f:
        f.write(f"{len(tokens)} {embedding_dim}\n")
        for token in tokens:
            vec_str = ' '.join(str(x) for x in embeddings[token])
            f.write(f"{token} {vec_str}\n")


def export_numpy(embeddings: dict, output_path: Path):
    """Export embeddings to NumPy format."""
    tokens = sorted(embeddings.keys())
    embedding_matrix = np.array([embeddings[t] for t in tokens])
    
    np.save(output_path, embedding_matrix)
    
    # Also save token list
    tokens_file = output_path.with_suffix('.tokens.txt')
    with open(tokens_file, 'w') as f:
        for token in tokens:
            f.write(f"{token}\n")


def main():
    """Main function."""
    args = parse_args()
    
    if not args.input.exists():
        logger.error(f"Input file not found: {args.input}")
        return 1
    
    logger.info(f"Loading embeddings from {args.input}")
    embeddings_data = io.load_json(args.input)
    
    # Convert to numpy arrays
    embeddings = {k: np.array(v) for k, v in embeddings_data.items()}
    logger.info(f"✅ Loaded embeddings for {len(embeddings)} events")
    
    paths.ensure_directory(args.output.parent)
    
    logger.info(f"Exporting embeddings to {args.format} format...")
    
    if args.format == "csv":
        export_csv(embeddings, args.output)
    elif args.format == "word2vec":
        export_word2vec(embeddings, args.output)
    elif args.format == "numpy":
        export_numpy(embeddings, args.output)
    
    logger.info(f"✅ Exported embeddings to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

