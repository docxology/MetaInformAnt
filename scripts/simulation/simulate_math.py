#!/usr/bin/env python3
"""Math simulation script.

This script generates synthetic data for mathematical biology models including
coalescent simulations, selection experiments, and population genetics models.

Usage:
    python3 scripts/simulation/simulate_math.py --type coalescent --n-samples 10 --n-loci 1000
    python3 scripts/simulation/simulate_math.py --type selection --n-generations 50 --population-size 100
    python3 scripts/simulation/simulate_math.py --type popgen --n-sequences 20 --length 1000
"""

import argparse
import random
import sys
from pathlib import Path

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core import io, logging, paths, validation
from metainformant.simulation.popgen import (
    generate_population_sequences,
    generate_site_frequency_spectrum,
    generate_two_populations,
)

logger = logging.get_logger(__name__)


def simulate_coalescent(
    output_dir: Path,
    n_samples: int,
    n_loci: int,
    effective_size: float,
    seed: int,
) -> dict:
    """Simulate coalescent data."""
    logger.info(f"Generating coalescent simulation: {n_samples} samples, {n_loci} loci")
    rng = random.Random(seed)
    
    # Use population genetics simulation to approximate coalescent
    sequences = generate_population_sequences(
        n_samples,
        n_loci,
        nucleotide_diversity=0.01,
        rng=rng,
    )
    
    # Calculate site frequency spectrum
    sfs = generate_site_frequency_spectrum(sequences, rng=rng)
    
    # Save sequences
    from metainformant.dna.sequences import write_fasta
    
    sequences_dict = {f"sample_{i:04d}": seq for i, seq in enumerate(sequences)}
    fasta_file = output_dir / "coalescent_sequences.fasta"
    write_fasta(sequences_dict, str(fasta_file))
    
    # Save SFS
    sfs_file = output_dir / "site_frequency_spectrum.json"
    io.dump_json(sfs, sfs_file, indent=2)
    
    logger.info(f"Coalescent simulation saved to {fasta_file}")
    
    return {
        "type": "coalescent",
        "n_samples": n_samples,
        "n_loci": n_loci,
        "effective_size": effective_size,
        "output_file": str(fasta_file),
        "sfs_file": str(sfs_file),
    }


def simulate_selection(
    output_dir: Path,
    n_generations: int,
    population_size: int,
    selection_coefficient: float,
    seed: int,
) -> dict:
    """Simulate selection experiment."""
    logger.info(f"Generating selection experiment: {n_generations} generations")
    rng = random.Random(seed)
    
    # Simulate allele frequency change under selection
    initial_frequency = 0.1
    frequencies = [initial_frequency]
    
    for generation in range(1, n_generations + 1):
        # Wright-Fisher with selection
        current_freq = frequencies[-1]
        # Fitness: w_A = 1 + s, w_a = 1
        fitness_ratio = (1 + selection_coefficient) * current_freq / (
            (1 + selection_coefficient) * current_freq + (1 - current_freq)
        )
        
        # Sample next generation
        next_freq = rng.betavariate(
            population_size * fitness_ratio,
            population_size * (1 - fitness_ratio),
        )
        frequencies.append(next_freq)
    
    # Save trajectory
    trajectory = [
        {"generation": gen, "allele_frequency": freq} for gen, freq in enumerate(frequencies)
    ]
    
    trajectory_file = output_dir / "selection_trajectory.json"
    io.dump_json(trajectory, trajectory_file, indent=2)
    
    logger.info(f"Selection experiment saved to {trajectory_file}")
    
    return {
        "type": "selection",
        "n_generations": n_generations,
        "population_size": population_size,
        "selection_coefficient": selection_coefficient,
        "output_file": str(trajectory_file),
    }


def simulate_popgen(
    output_dir: Path,
    n_sequences: int,
    sequence_length: int,
    diversity: float,
    fst: float | None,
    seed: int,
) -> dict:
    """Simulate population genetics data."""
    logger.info(f"Generating popgen data: {n_sequences} sequences")
    rng = random.Random(seed)
    
    if fst is not None:
        # Two populations
        n_pop1 = n_sequences // 2
        n_pop2 = n_sequences - n_pop1
        pop1, pop2 = generate_two_populations(
            n_pop1, n_pop2, sequence_length, fst=fst, within_pop_diversity=diversity, rng=rng
        )
        
        from metainformant.dna.sequences import write_fasta
        
        sequences_dict = {}
        for i, seq in enumerate(pop1):
            sequences_dict[f"pop1_seq_{i:04d}"] = seq
        for i, seq in enumerate(pop2):
            sequences_dict[f"pop2_seq_{i:04d}"] = seq
        
        fasta_file = output_dir / "popgen_sequences.fasta"
        write_fasta(sequences_dict, str(fasta_file))
    else:
        # Single population
        sequences = generate_population_sequences(
            n_sequences, sequence_length, nucleotide_diversity=diversity, rng=rng
        )
        
        from metainformant.dna.sequences import write_fasta
        
        sequences_dict = {f"seq_{i:04d}": seq for i, seq in enumerate(sequences)}
        fasta_file = output_dir / "popgen_sequences.fasta"
        write_fasta(sequences_dict, str(fasta_file))
    
    logger.info(f"Popgen data saved to {fasta_file}")
    
    return {
        "type": "popgen",
        "n_sequences": n_sequences,
        "length": sequence_length,
        "diversity": diversity,
        "fst": fst,
        "output_file": str(fasta_file),
    }


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Mathematical biology simulation for testing and validation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Simulate coalescent
  %(prog)s --type coalescent --n-samples 10 --n-loci 1000

  # Simulate selection experiment
  %(prog)s --type selection --n-generations 50 --population-size 100 --selection 0.1

  # Simulate population genetics
  %(prog)s --type popgen --n-sequences 20 --length 1000 --diversity 0.01
        """,
    )
    
    parser.add_argument(
        "--type",
        required=True,
        choices=["coalescent", "selection", "popgen"],
        help="Type of simulation to run",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/simulation/math"),
        help="Output directory (default: output/simulation/math)",
    )
    parser.add_argument("--n-samples", type=int, default=10, help="Number of samples (coalescent type)")
    parser.add_argument("--n-loci", type=int, default=1000, help="Number of loci (coalescent type)")
    parser.add_argument("--effective-size", type=float, default=1000.0, help="Effective population size (coalescent type)")
    parser.add_argument("--n-generations", type=int, default=50, help="Number of generations (selection type)")
    parser.add_argument("--population-size", type=int, default=100, help="Population size (selection type)")
    parser.add_argument("--selection", type=float, default=0.1, help="Selection coefficient (selection type)")
    parser.add_argument("--n-sequences", type=int, default=20, help="Number of sequences (popgen type)")
    parser.add_argument("--length", type=int, default=1000, help="Sequence length (popgen type)")
    parser.add_argument("--diversity", type=float, default=0.01, help="Nucleotide diversity (popgen type)")
    parser.add_argument("--fst", type=float, help="Fst value for two populations (popgen type)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")
    
    args = parser.parse_args()
    
    if args.verbose:
        import logging as std_logging
        logger.setLevel(std_logging.DEBUG)
    
    output_dir = paths.ensure_directory(args.output)
    
    try:
        if args.type == "coalescent":
            results = simulate_coalescent(
                output_dir, args.n_samples, args.n_loci, args.effective_size, args.seed
            )
        elif args.type == "selection":
            results = simulate_selection(
                output_dir,
                args.n_generations,
                args.population_size,
                args.selection,
                args.seed,
            )
        elif args.type == "popgen":
            results = simulate_popgen(
                output_dir, args.n_sequences, args.length, args.diversity, args.fst, args.seed
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

