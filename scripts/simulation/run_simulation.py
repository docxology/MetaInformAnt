#!/usr/bin/env python3
"""Biological simulation workflow orchestrator.

This script provides comprehensive orchestration for biological simulation workflows,
including sequence generation, agent-based models, and expression simulation.

Usage:
    python3 scripts/simulation/run_simulation.py --model sequences --n 1000 --length 500 --output output/simulation/sequences
    python3 scripts/simulation/run_simulation.py --model agents --width 50 --height 50 --steps 100
    python3 scripts/simulation/run_simulation.py --model expression --num-genes 1000 --num-samples 20
    python3 scripts/simulation/run_simulation.py --help
"""

import argparse
import logging
import random
import sys
from pathlib import Path
from typing import Any

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core.io import dump_json, ensure_directory
from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Biological simulation workflow orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate random DNA sequences
  %(prog)s --model sequences --n 1000 --length 500 --output output/simulation/dna

  # Run agent-based model
  %(prog)s --model agents --width 50 --height 50 --num-agents 100 --steps 100

  # Simulate gene expression counts
  %(prog)s --model expression --num-genes 1000 --num-samples 20 --output output/simulation/expression
        """,
    )
    parser.add_argument(
        "--model",
        type=str,
        required=True,
        choices=["sequences", "agents", "expression"],
        help="Simulation model type",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/simulation"),
        help="Output directory (default: output/simulation)",
    )
    # Sequence simulation parameters
    parser.add_argument("--n", type=int, default=100, help="Number of sequences (sequences model)")
    parser.add_argument("--length", type=int, default=500, help="Sequence length (sequences model)")
    parser.add_argument("--gc-content", type=float, default=0.5, help="GC content (sequences model)")
    parser.add_argument("--mutations", type=int, default=0, help="Number of mutations per sequence (sequences model)")
    # Agent-based model parameters
    parser.add_argument("--width", type=int, default=50, help="Grid width (agents model)")
    parser.add_argument("--height", type=int, default=50, help="Grid height (agents model)")
    parser.add_argument("--num-agents", type=int, default=100, help="Number of agents (agents model)")
    parser.add_argument("--steps", type=int, default=100, help="Number of simulation steps (agents model)")
    # Expression simulation parameters
    parser.add_argument("--num-genes", type=int, default=1000, help="Number of genes (expression model)")
    parser.add_argument("--num-samples", type=int, default=20, help="Number of samples (expression model)")
    parser.add_argument("--mean-expression", type=float, default=100.0, help="Mean expression (expression model)")
    parser.add_argument("--dispersion", type=float, default=0.1, help="Dispersion (expression model)")
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility",
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


def run_sequence_simulation(args, output_dir: Path) -> dict[str, Any]:
    """Run sequence generation simulation."""
    logger.info(f"Generating {args.n} sequences of length {args.length}...")
    from metainformant.dna.sequences import write_fasta
    from metainformant.simulation import generate_random_dna, mutate_sequence

    rng = random.Random(args.seed)
    sequences = {}

    for i in range(args.n):
        seq = generate_random_dna(args.length, gc_content=args.gc_content, rng=rng)
        if args.mutations > 0:
            seq = mutate_sequence(seq, n_mut=args.mutations, rng=rng)
        sequences[f"seq_{i}"] = seq

    # Save sequences
    fasta_file = output_dir / "simulated_sequences.fasta"
    write_fasta(sequences, str(fasta_file))
    logger.info(f"Sequences saved to {fasta_file}")

    return {
        "n_sequences": args.n,
        "length": args.length,
        "gc_content": args.gc_content,
        "mutations_per_seq": args.mutations,
        "output_file": str(fasta_file),
    }


def run_agent_simulation(args, output_dir: Path) -> dict[str, Any]:
    """Run agent-based model simulation."""
    logger.info(f"Running agent-based model: {args.width}x{args.height} grid, {args.num_agents} agents...")
    from metainformant.simulation import GridWorld

    rng = random.Random(args.seed)
    world = GridWorld(width=args.width, height=args.height, num_agents=args.num_agents, rng=rng)

    # Run simulation
    position_history = []
    for step in range(args.steps):
        world.step()
        if step % 10 == 0:  # Record every 10 steps
            positions = world.positions()
            position_history.append({"step": step, "positions": positions[:100]})  # Limit output size

    results = {
        "grid_size": [args.width, args.height],
        "num_agents": args.num_agents,
        "steps": args.steps,
        "position_history": position_history,
        "final_positions": world.positions()[:100],  # Limit output size
    }

    output_file = output_dir / "agent_simulation.json"
    dump_json(results, output_file)
    logger.info(f"Agent simulation saved to {output_file}")

    return results


def run_expression_simulation(args, output_dir: Path) -> dict[str, Any]:
    """Run expression count simulation."""
    logger.info(f"Simulating expression: {args.num_genes} genes, {args.num_samples} samples...")
    import pandas as pd

    from metainformant.simulation import simulate_counts_negative_binomial

    rng = random.Random(args.seed)
    counts = simulate_counts_negative_binomial(
        num_genes=args.num_genes,
        num_samples=args.num_samples,
        mean_expression=args.mean_expression,
        dispersion=args.dispersion,
        rng=rng,
    )

    # Save as CSV
    df = pd.DataFrame(counts, columns=[f"sample_{i}" for i in range(args.num_samples)])
    df.index = [f"gene_{i}" for i in range(args.num_genes)]

    csv_file = output_dir / "simulated_expression.csv"
    df.to_csv(csv_file)
    logger.info(f"Expression counts saved to {csv_file}")

    return {
        "num_genes": args.num_genes,
        "num_samples": args.num_samples,
        "mean_expression": args.mean_expression,
        "dispersion": args.dispersion,
        "output_file": str(csv_file),
        "shape": list(df.shape),
    }


def run_workflow(args):
    """Execute biological simulation workflow."""
    logger.info("Starting biological simulation workflow")
    logger.info(f"Model: {args.model}")
    logger.info(f"Output: {args.output}")

    if args.dry_run:
        logger.info("DRY RUN - no changes will be made")
        return

    output_dir = ensure_directory(args.output)
    logger.info(f"Output directory: {output_dir}")

    workflow_results = {
        "model_type": args.model,
        "output_dir": str(output_dir),
        "seed": args.seed,
        "results": {},
    }

    try:
        if args.model == "sequences":
            results = run_sequence_simulation(args, output_dir)
            workflow_results["results"] = results
        elif args.model == "agents":
            results = run_agent_simulation(args, output_dir)
            workflow_results["results"] = results
        elif args.model == "expression":
            results = run_expression_simulation(args, output_dir)
            workflow_results["results"] = results

        # Save summary
        summary_file = output_dir / "workflow_summary.json"
        dump_json(workflow_results, summary_file, indent=2)
        logger.info(f"Workflow summary saved to {summary_file}")

    except Exception as e:
        logger.error(f"Simulation failed: {e}", exc_info=True)
        raise

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
