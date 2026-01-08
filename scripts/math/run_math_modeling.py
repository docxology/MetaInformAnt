#!/usr/bin/env python3
"""Mathematical biology workflow orchestrator.

This script provides comprehensive orchestration for mathematical biology workflows,
including population dynamics, epidemiology, selection models, and population genetics.

Usage:
    python3 scripts/math/run_math_modeling.py --model sir --output output/math/sir_simulation
    python3 scripts/math/run_math_modeling.py --model logistic --r 3.5 --steps 100
    python3 scripts/math/run_math_modeling.py --model selection --generations 50
    python3 scripts/math/run_math_modeling.py --help
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Any

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core.io import ensure_directory, dump_json
from metainformant.core.logging import get_logger

logger = get_logger(__name__)


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Mathematical biology workflow orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # SIR epidemiology model
  %(prog)s --model sir --beta 0.3 --gamma 0.1 --S0 990 --I0 10 --steps 100

  # Logistic growth model
  %(prog)s --model logistic --r 3.5 --x0 0.5 --steps 100

  # Selection experiment simulation
  %(prog)s --model selection --generations 50 --n 10000 --s_hat 0.6

  # Price equation analysis
  %(prog)s --model price --fitness 0.2,0.4,0.1 --trait-parent 1.0,1.2,0.9 --trait-offspring 1.25,1.35,0.95
        """
    )
    parser.add_argument(
        "--model",
        type=str,
        required=True,
        choices=["sir", "logistic", "selection", "price", "lotka-volterra"],
        help="Model type to run",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/math"),
        help="Output directory (default: output/math)",
    )
    # SIR model parameters
    parser.add_argument("--S0", type=float, default=990.0, help="Initial susceptible population (SIR)")
    parser.add_argument("--I0", type=float, default=10.0, help="Initial infected population (SIR)")
    parser.add_argument("--R0_init", type=float, default=0.0, help="Initial recovered population (SIR)")
    parser.add_argument("--beta", type=float, default=0.3, help="Transmission rate (SIR)")
    parser.add_argument("--gamma", type=float, default=0.1, help="Recovery rate (SIR)")
    # Logistic model parameters
    parser.add_argument("--r", type=float, default=3.5, help="Growth rate (logistic)")
    parser.add_argument("--x0", type=float, default=0.5, help="Initial value (logistic)")
    # Selection model parameters
    parser.add_argument("--generations", type=int, default=50, help="Number of generations (selection)")
    parser.add_argument("--n", type=int, default=10000, help="Population size (selection)")
    parser.add_argument("--s_hat", type=float, default=0.6, help="Signal fidelity (selection)")
    # General parameters
    parser.add_argument("--steps", type=int, default=100, help="Number of simulation steps")
    parser.add_argument("--dt", type=float, default=0.01, help="Time step size")
    # Price equation parameters (comma-separated lists)
    parser.add_argument("--fitness", type=str, help="Comma-separated fitness values (Price)")
    parser.add_argument("--trait-parent", type=str, help="Comma-separated parent trait values (Price)")
    parser.add_argument("--trait-offspring", type=str, help="Comma-separated offspring trait values (Price)")
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


def run_sir_model(args, output_dir: Path) -> dict[str, Any]:
    """Run SIR epidemiology model."""
    logger.info("Running SIR epidemiology model...")
    from metainformant.math import sir_step, basic_reproduction_number
    
    S, I, R = args.S0, args.I0, args.R0_init
    beta, gamma = args.beta, args.gamma
    
    results = {
        "time_points": [],
        "S": [],
        "I": [],
        "R": [],
    }
    
    for t in range(args.steps):
        results["time_points"].append(t * args.dt)
        results["S"].append(S)
        results["I"].append(I)
        results["R"].append(R)
        
        S, I, R = sir_step(S, I, R, beta, gamma, dt=args.dt)
    
    R0 = basic_reproduction_number(beta, gamma)
    results["R0"] = R0
    
    output_file = output_dir / "sir_simulation.json"
    dump_json(results, output_file)
    logger.info(f"SIR simulation saved to {output_file}")
    logger.info(f"R0 = {R0:.3f}")
    
    return results


def run_logistic_model(args, output_dir: Path) -> dict[str, Any]:
    """Run logistic growth model."""
    logger.info("Running logistic growth model...")
    from metainformant.math import logistic_map
    
    trajectory = logistic_map(r=args.r, x0=args.x0, n_iterations=args.steps)
    
    results = {
        "r": args.r,
        "x0": args.x0,
        "trajectory": trajectory,
        "final_value": trajectory[-1],
    }
    
    output_file = output_dir / "logistic_simulation.json"
    dump_json(results, output_file)
    logger.info(f"Logistic simulation saved to {output_file}")
    
    return results


def run_selection_model(args, output_dir: Path) -> dict[str, Any]:
    """Run selection experiment simulation."""
    logger.info("Running selection experiment simulation...")
    import numpy as np
    from metainformant.math.selection_experiments import simulate_generations
    
    results_obj = simulate_generations(
        generations=args.generations,
        n=args.n,
        s_hat=args.s_hat
    )
    
    results = {
        "generations": args.generations,
        "n": args.n,
        "s_hat": args.s_hat,
        "mean_s": results_obj.mean_s.tolist() if hasattr(results_obj.mean_s, 'tolist') else list(results_obj.mean_s),
        "mean_q": results_obj.mean_q.tolist() if hasattr(results_obj.mean_q, 'tolist') else list(results_obj.mean_q),
    }
    
    output_file = output_dir / "selection_simulation.json"
    dump_json(results, output_file)
    logger.info(f"Selection simulation saved to {output_file}")
    
    return results


def run_price_equation(args, output_dir: Path) -> dict[str, Any]:
    """Run Price equation analysis."""
    logger.info("Running Price equation analysis...")
    from metainformant.math import price_equation
    
    # Parse comma-separated values
    fitness = [float(x.strip()) for x in args.fitness.split(",")]
    trait_parent = [float(x.strip()) for x in args.trait_parent.split(",")]
    trait_offspring = [float(x.strip()) for x in args.trait_offspring.split(",")]
    
    if len(fitness) != len(trait_parent) or len(fitness) != len(trait_offspring):
        raise ValueError("Fitness, trait_parent, and trait_offspring must have same length")
    
    cov_term, trans_term, total = price_equation(fitness, trait_parent, trait_offspring)
    
    results = {
        "covariance_term": float(cov_term),
        "transmission_term": float(trans_term),
        "total_change": float(total),
    }
    
    output_file = output_dir / "price_equation.json"
    dump_json(results, output_file)
    logger.info(f"Price equation analysis saved to {output_file}")
    logger.info(f"Total change: {total:.6f} (covariance: {cov_term:.6f}, transmission: {trans_term:.6f})")
    
    return results


def run_lotka_volterra(args, output_dir: Path) -> dict[str, Any]:
    """Run Lotka-Volterra predator-prey model."""
    logger.info("Running Lotka-Volterra model...")
    from metainformant.math import lotka_volterra_step
    
    prey, predator = 100.0, 10.0
    alpha, beta, delta, gamma = 1.0, 0.1, 0.075, 1.5
    
    results = {
        "time_points": [],
        "prey": [],
        "predator": [],
    }
    
    for t in range(args.steps):
        results["time_points"].append(t * args.dt)
        results["prey"].append(prey)
        results["predator"].append(predator)
        
        prey, predator = lotka_volterra_step(
            prey, predator, alpha, beta, delta, gamma, dt=args.dt
        )
    
    output_file = output_dir / "lotka_volterra_simulation.json"
    dump_json(results, output_file)
    logger.info(f"Lotka-Volterra simulation saved to {output_file}")
    
    return results


def run_workflow(args):
    """Execute mathematical modeling workflow."""
    logger.info("Starting mathematical biology workflow")
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
        "results": {},
    }
    
    try:
        if args.model == "sir":
            results = run_sir_model(args, output_dir)
            workflow_results["results"] = results
        elif args.model == "logistic":
            results = run_logistic_model(args, output_dir)
            workflow_results["results"] = results
        elif args.model == "selection":
            results = run_selection_model(args, output_dir)
            workflow_results["results"] = results
        elif args.model == "price":
            if not all([args.fitness, args.trait_parent, args.trait_offspring]):
                raise ValueError("Price equation requires --fitness, --trait-parent, and --trait-offspring")
            results = run_price_equation(args, output_dir)
            workflow_results["results"] = results
        elif args.model == "lotka-volterra":
            results = run_lotka_volterra(args, output_dir)
            workflow_results["results"] = results
        
        # Save summary
        summary_file = output_dir / "workflow_summary.json"
        dump_json(workflow_results, summary_file, indent=2)
        logger.info(f"Workflow summary saved to {summary_file}")
        
    except Exception as e:
        logger.error(f"Model simulation failed: {e}", exc_info=True)
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




