#!/usr/bin/env python3
"""RNA expression simulation script.

This script generates synthetic RNA expression data including count matrices,
differential expression patterns, batch effects, and transcript abundance.

Usage:
    python3 scripts/simulation/simulate_rna.py --type counts --num-genes 1000 --num-samples 20
    python3 scripts/simulation/simulate_rna.py --type differential --num-genes 500 --num-samples 10 --n-de-genes 50
    python3 scripts/simulation/simulate_rna.py --type batch --num-genes 2000 --num-samples 30 --n-batches 3
"""

import argparse
import random
import sys
from pathlib import Path

import pandas as pd

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core import io, logging, paths, validation
from metainformant.simulation.rna import simulate_counts_negative_binomial

logger = logging.get_logger(__name__)


def simulate_counts(
    output_dir: Path,
    num_genes: int,
    num_samples: int,
    mean_expression: float,
    dispersion: float,
    seed: int,
) -> dict:
    """Simulate basic RNA-seq count matrix.

    Args:
        output_dir: Output directory for results
        num_genes: Number of genes
        num_samples: Number of samples
        mean_expression: Mean expression level
        dispersion: Dispersion parameter (must be > 0)
        seed: Random seed for reproducibility

    Returns:
        Dictionary with simulation results and metadata

    Raises:
        ValidationError: If parameters are invalid
    """
    # Validate parameters
    validation.validate_range(num_genes, min_val=1, name="num_genes")
    validation.validate_range(num_samples, min_val=1, name="num_samples")
    validation.validate_range(mean_expression, min_val=0.0, name="mean_expression")
    validation.validate_range(dispersion, min_val=0.0, name="dispersion")

    logger.info(f"Simulating expression counts: {num_genes} genes x {num_samples} samples")
    rng = random.Random(seed)

    counts = simulate_counts_negative_binomial(
        num_genes=num_genes,
        num_samples=num_samples,
        mean_expression=mean_expression,
        dispersion=dispersion,
        rng=rng,
    )

    # Save as CSV
    df = pd.DataFrame(counts, columns=[f"sample_{i:03d}" for i in range(num_samples)])
    df.index = [f"gene_{i:05d}" for i in range(num_genes)]

    csv_file = output_dir / "expression_counts.csv"
    df.to_csv(csv_file)

    logger.info(f"Expression counts saved to {csv_file}")

    return {
        "type": "counts",
        "num_genes": num_genes,
        "num_samples": num_samples,
        "mean_expression": mean_expression,
        "dispersion": dispersion,
        "output_file": str(csv_file),
    }


def simulate_differential(
    output_dir: Path,
    num_genes: int,
    num_samples: int,
    n_de_genes: int,
    fold_change: float,
    mean_expression: float,
    dispersion: float,
    seed: int,
) -> dict:
    """Simulate differential expression between two groups.

    Args:
        output_dir: Output directory for results
        num_genes: Number of genes
        num_samples: Number of samples (split into two groups)
        n_de_genes: Number of differentially expressed genes
        fold_change: Fold change for DE genes
        mean_expression: Mean expression level
        dispersion: Dispersion parameter
        seed: Random seed for reproducibility

    Returns:
        Dictionary with simulation results and metadata

    Raises:
        ValidationError: If parameters are invalid
    """
    # Validate parameters
    validation.validate_range(num_genes, min_val=1, name="num_genes")
    validation.validate_range(num_samples, min_val=2, name="num_samples")
    validation.validate_range(n_de_genes, min_val=0, max_val=num_genes, name="n_de_genes")
    validation.validate_range(fold_change, min_val=1.0, name="fold_change")
    validation.validate_range(mean_expression, min_val=0.0, name="mean_expression")
    validation.validate_range(dispersion, min_val=0.0, name="dispersion")

    logger.info(f"Simulating differential expression: {n_de_genes} DE genes")
    rng = random.Random(seed)

    n_group1 = num_samples // 2
    n_group2 = num_samples - n_group1

    # Select DE genes
    n_de_genes = min(n_de_genes, num_genes)  # Ensure we don't exceed available genes
    de_genes = set(rng.sample(range(num_genes), n_de_genes))

    counts = []
    gene_info = []

    for gene_idx in range(num_genes):
        is_de = gene_idx in de_genes
        direction = rng.choice([-1, 1]) if is_de else 1

        # Group 1 expression
        mean1 = mean_expression
        row1 = [
            simulate_counts_negative_binomial(1, 1, mean_expression=mean1, dispersion=dispersion, rng=rng)[0][0]
            for _ in range(n_group1)
        ]

        # Group 2 expression (different if DE)
        mean2 = mean_expression * (fold_change**direction) if is_de else mean_expression
        row2 = [
            simulate_counts_negative_binomial(1, 1, mean_expression=mean2, dispersion=dispersion, rng=rng)[0][0]
            for _ in range(n_group2)
        ]

        counts.append(row1 + row2)
        gene_info.append(
            {
                "gene_id": f"gene_{gene_idx:05d}",
                "is_de": is_de,
                "fold_change": fold_change**direction if is_de else 1.0,
            }
        )

    # Save counts
    df = pd.DataFrame(counts, columns=[f"sample_{i:03d}" for i in range(num_samples)])
    df.index = [f"gene_{i:05d}" for i in range(num_genes)]
    csv_file = output_dir / "differential_expression.csv"
    df.to_csv(csv_file)

    # Save gene info
    gene_df = pd.DataFrame(gene_info)
    gene_file = output_dir / "gene_info.csv"
    gene_df.to_csv(gene_file, index=False)

    logger.info(f"Differential expression data saved to {csv_file}")

    return {
        "type": "differential",
        "num_genes": num_genes,
        "num_samples": num_samples,
        "n_de_genes": n_de_genes,
        "fold_change": fold_change,
        "output_file": str(csv_file),
        "gene_info_file": str(gene_file),
    }


def simulate_batch(
    output_dir: Path,
    num_genes: int,
    num_samples: int,
    n_batches: int,
    batch_effect_size: float,
    mean_expression: float,
    dispersion: float,
    seed: int,
) -> dict:
    """Simulate batch effects in expression data.

    Args:
        output_dir: Output directory for results
        num_genes: Number of genes
        num_samples: Number of samples
        n_batches: Number of batches
        batch_effect_size: Size of batch effects
        mean_expression: Mean expression level
        dispersion: Dispersion parameter
        seed: Random seed for reproducibility

    Returns:
        Dictionary with simulation results and metadata

    Raises:
        ValidationError: If parameters are invalid
    """
    # Validate parameters
    validation.validate_range(num_genes, min_val=1, name="num_genes")
    validation.validate_range(num_samples, min_val=1, name="num_samples")
    validation.validate_range(n_batches, min_val=1, max_val=num_samples, name="n_batches")
    validation.validate_range(batch_effect_size, min_val=0.0, name="batch_effect_size")
    validation.validate_range(mean_expression, min_val=0.0, name="mean_expression")
    validation.validate_range(dispersion, min_val=0.0, name="dispersion")

    logger.info(f"Simulating batch effects: {n_batches} batches")
    rng = random.Random(seed)

    samples_per_batch = num_samples // n_batches

    # Generate batch effects
    batch_effects = [rng.gauss(0, batch_effect_size) for _ in range(n_batches)]

    counts = []
    sample_info = []

    for sample_idx in range(num_samples):
        batch_idx = sample_idx // samples_per_batch
        if batch_idx >= n_batches:
            batch_idx = n_batches - 1

        batch_effect = batch_effects[batch_idx]
        sample_mean = mean_expression * (2**batch_effect)

        row = [
            simulate_counts_negative_binomial(1, 1, mean_expression=sample_mean, dispersion=dispersion, rng=rng)[0][0]
            for _ in range(num_genes)
        ]
        counts.append(row)

        sample_info.append(
            {
                "sample_id": f"sample_{sample_idx:03d}",
                "batch": batch_idx,
                "batch_effect": batch_effect,
            }
        )

    # Transpose to genes x samples
    df = pd.DataFrame(counts).T
    df.columns = [f"sample_{i:03d}" for i in range(num_samples)]
    df.index = [f"gene_{i:05d}" for i in range(num_genes)]

    csv_file = output_dir / "batch_expression.csv"
    df.to_csv(csv_file)

    # Save sample info
    sample_df = pd.DataFrame(sample_info)
    sample_file = output_dir / "sample_info.csv"
    sample_df.to_csv(sample_file, index=False)

    logger.info(f"Batch effect data saved to {csv_file}")

    return {
        "type": "batch",
        "num_genes": num_genes,
        "num_samples": num_samples,
        "n_batches": n_batches,
        "batch_effect_size": batch_effect_size,
        "output_file": str(csv_file),
        "sample_info_file": str(sample_file),
    }


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="RNA expression simulation for testing and validation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Simulate basic count matrix
  %(prog)s --type counts --num-genes 1000 --num-samples 20

  # Simulate differential expression
  %(prog)s --type differential --num-genes 500 --num-samples 10 --n-de-genes 50

  # Simulate batch effects
  %(prog)s --type batch --num-genes 2000 --num-samples 30 --n-batches 3
        """,
    )

    parser.add_argument(
        "--type",
        required=True,
        choices=["counts", "differential", "batch"],
        help="Type of simulation to run",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/simulation/rna"),
        help="Output directory (default: output/simulation/rna)",
    )
    parser.add_argument("--num-genes", type=int, default=1000, help="Number of genes")
    parser.add_argument("--num-samples", type=int, default=20, help="Number of samples")
    parser.add_argument("--mean-expression", type=float, default=100.0, help="Mean expression level")
    parser.add_argument("--dispersion", type=float, default=0.1, help="Dispersion parameter")
    parser.add_argument("--n-de-genes", type=int, default=50, help="Number of DE genes (differential type)")
    parser.add_argument("--fold-change", type=float, default=2.0, help="Fold change for DE genes (differential type)")
    parser.add_argument("--n-batches", type=int, default=3, help="Number of batches (batch type)")
    parser.add_argument("--batch-effect-size", type=float, default=0.5, help="Batch effect size (batch type)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")

    args = parser.parse_args()

    if args.verbose:
        import logging as std_logging

        logger.setLevel(std_logging.DEBUG)

    # Validate output directory
    output_dir = paths.ensure_directory(args.output)

    # Validate common parameters
    validation.validate_range(args.num_genes, min_val=1, name="num_genes")
    validation.validate_range(args.num_samples, min_val=1, name="num_samples")
    validation.validate_range(args.mean_expression, min_val=0.0, name="mean_expression")
    validation.validate_range(args.dispersion, min_val=0.0, name="dispersion")
    if hasattr(args, "n_de_genes"):
        validation.validate_range(args.n_de_genes, min_val=0, max_val=args.num_genes, name="n_de_genes")
    if hasattr(args, "fold_change"):
        validation.validate_range(args.fold_change, min_val=1.0, name="fold_change")
    if hasattr(args, "n_batches"):
        validation.validate_range(args.n_batches, min_val=1, max_val=args.num_samples, name="n_batches")

    try:
        if args.type == "counts":
            results = simulate_counts(
                output_dir,
                args.num_genes,
                args.num_samples,
                args.mean_expression,
                args.dispersion,
                args.seed,
            )
        elif args.type == "differential":
            results = simulate_differential(
                output_dir,
                args.num_genes,
                args.num_samples,
                args.n_de_genes,
                args.fold_change,
                args.mean_expression,
                args.dispersion,
                args.seed,
            )
        elif args.type == "batch":
            results = simulate_batch(
                output_dir,
                args.num_genes,
                args.num_samples,
                args.n_batches,
                args.batch_effect_size,
                args.mean_expression,
                args.dispersion,
                args.seed,
            )

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
