#!/usr/bin/env python3
"""DNA sequence simulation script.

This script generates synthetic DNA sequences for testing and validation,
including sequences, mutations, population data, alignment test data, and
phylogenetic tree test data.

Usage:
    python3 scripts/simulation/simulate_dna.py --type sequences --n 100 --length 1000
    python3 scripts/simulation/simulate_dna.py --type population --n 50 --length 2000 --diversity 0.01
    python3 scripts/simulation/simulate_dna.py --type alignment --n 10 --length 500
    python3 scripts/simulation/simulate_dna.py --type phylogeny --n 20 --length 1000
"""

import argparse
import random
import sys
from pathlib import Path

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core import io, logging, paths, validation
from metainformant.dna.sequences import write_fasta
from metainformant.simulation import (
    generate_population_sequences,
    generate_random_dna,
    generate_two_populations,
    mutate_sequence,
)

logger = logging.get_logger(__name__)


def simulate_sequences(
    output_dir: Path,
    n_sequences: int,
    sequence_length: int,
    gc_content: float,
    mutations: int,
    seed: int,
) -> dict:
    """Simulate basic DNA sequences.
    
    Args:
        output_dir: Output directory for results
        n_sequences: Number of sequences to generate
        sequence_length: Length of each sequence
        gc_content: GC content (0.0-1.0)
        mutations: Number of mutations per sequence
        seed: Random seed for reproducibility
        
    Returns:
        Dictionary with simulation results and metadata
        
    Raises:
        ValidationError: If parameters are invalid
    """
    # Validate parameters
    validation.validate_range(n_sequences, min_val=1, name="n_sequences")
    validation.validate_range(sequence_length, min_val=1, name="sequence_length")
    validation.validate_range(gc_content, min_val=0.0, max_val=1.0, name="gc_content")
    validation.validate_range(mutations, min_val=0, name="mutations")
    
    logger.info(f"Generating {n_sequences} DNA sequences of length {sequence_length}")
    rng = random.Random(seed)
    
    sequences = {}
    for i in range(n_sequences):
        seq = generate_random_dna(sequence_length, gc_content=gc_content, rng=rng)
        if mutations > 0:
            seq = mutate_sequence(seq, n_mut=mutations, rng=rng)
        sequences[f"seq_{i:04d}"] = seq
    
    fasta_file = output_dir / "sequences.fasta"
    write_fasta(sequences, str(fasta_file))
    logger.info(f"Sequences saved to {fasta_file}")
    
    return {
        "type": "sequences",
        "n_sequences": n_sequences,
        "length": sequence_length,
        "gc_content": gc_content,
        "mutations_per_seq": mutations,
        "output_file": str(fasta_file),
    }


def simulate_population(
    output_dir: Path,
    n_sequences: int,
    sequence_length: int,
    diversity: float | None,
    gc_content: float,
    seed: int,
) -> dict:
    """Simulate population genetics data.
    
    Args:
        output_dir: Output directory for results
        n_sequences: Number of sequences in population
        sequence_length: Length of each sequence
        diversity: Nucleotide diversity (None for default)
        gc_content: GC content (0.0-1.0)
        seed: Random seed for reproducibility
        
    Returns:
        Dictionary with simulation results and metadata
        
    Raises:
        ValidationError: If parameters are invalid
    """
    # Validate parameters
    validation.validate_range(n_sequences, min_val=2, name="n_sequences")
    validation.validate_range(sequence_length, min_val=1, name="sequence_length")
    if diversity is not None:
        validation.validate_range(diversity, min_val=0.0, name="diversity")
    validation.validate_range(gc_content, min_val=0.0, max_val=1.0, name="gc_content")
    
    logger.info(f"Generating population: {n_sequences} sequences with diversity {diversity}")
    rng = random.Random(seed)
    
    sequences_list = generate_population_sequences(
        n_sequences,
        sequence_length,
        nucleotide_diversity=diversity,
        gc_content=gc_content,
        rng=rng,
    )
    
    sequences = {f"pop_seq_{i:04d}": seq for i, seq in enumerate(sequences_list)}
    fasta_file = output_dir / "population.fasta"
    write_fasta(sequences, str(fasta_file))
    
    # Also save as JSON for metadata
    metadata = {
        "n_sequences": n_sequences,
        "sequence_length": sequence_length,
        "nucleotide_diversity": diversity,
        "gc_content": gc_content,
    }
    metadata_file = output_dir / "population_metadata.json"
    io.dump_json(metadata, metadata_file)
    
    logger.info(f"Population sequences saved to {fasta_file}")
    
    return {
        "type": "population",
        "output_file": str(fasta_file),
        "metadata_file": str(metadata_file),
        **metadata,
    }


def simulate_alignment(
    output_dir: Path,
    n_sequences: int,
    sequence_length: int,
    gc_content: float,
    mutation_rate: float,
    seed: int,
) -> dict:
    """Simulate sequences for alignment testing.
    
    Args:
        output_dir: Output directory for results
        n_sequences: Number of sequences to generate
        sequence_length: Length of each sequence
        gc_content: GC content (0.0-1.0)
        mutation_rate: Mutation rate per site (0.0-1.0)
        seed: Random seed for reproducibility
        
    Returns:
        Dictionary with simulation results and metadata
        
    Raises:
        ValidationError: If parameters are invalid
    """
    # Validate parameters
    validation.validate_range(n_sequences, min_val=2, name="n_sequences")
    validation.validate_range(sequence_length, min_val=1, name="sequence_length")
    validation.validate_range(gc_content, min_val=0.0, max_val=1.0, name="gc_content")
    validation.validate_range(mutation_rate, min_val=0.0, max_val=1.0, name="mutation_rate")
    
    logger.info(f"Generating alignment test data: {n_sequences} sequences")
    rng = random.Random(seed)
    
    # Generate reference sequence
    reference = generate_random_dna(sequence_length, gc_content=gc_content, rng=rng)
    
    # Generate related sequences with mutations
    sequences = {"reference": reference}
    n_mutations = int(sequence_length * mutation_rate)
    
    for i in range(n_sequences - 1):
        seq = mutate_sequence(reference, n_mut=n_mutations, rng=rng)
        sequences[f"aligned_{i:04d}"] = seq
    
    fasta_file = output_dir / "alignment_test.fasta"
    write_fasta(sequences, str(fasta_file))
    
    logger.info(f"Alignment test sequences saved to {fasta_file}")
    
    return {
        "type": "alignment",
        "n_sequences": n_sequences,
        "length": sequence_length,
        "mutation_rate": mutation_rate,
        "output_file": str(fasta_file),
    }


def simulate_phylogeny(
    output_dir: Path,
    n_sequences: int,
    sequence_length: int,
    gc_content: float,
    diversity: float,
    seed: int,
) -> dict:
    """Simulate sequences for phylogenetic tree construction.
    
    Args:
        output_dir: Output directory for results
        n_sequences: Number of sequences (minimum 3 for tree)
        sequence_length: Length of each sequence
        gc_content: GC content (0.0-1.0)
        diversity: Nucleotide diversity
        seed: Random seed for reproducibility
        
    Returns:
        Dictionary with simulation results and metadata
        
    Raises:
        ValidationError: If parameters are invalid
    """
    # Validate parameters
    validation.validate_range(n_sequences, min_val=3, name="n_sequences")
    validation.validate_range(sequence_length, min_val=1, name="sequence_length")
    validation.validate_range(gc_content, min_val=0.0, max_val=1.0, name="gc_content")
    validation.validate_range(diversity, min_val=0.0, name="diversity")
    
    logger.info(f"Generating phylogeny test data: {n_sequences} sequences")
    rng = random.Random(seed)
    
    # Generate population with diversity
    sequences_list = generate_population_sequences(
        n_sequences,
        sequence_length,
        nucleotide_diversity=diversity,
        gc_content=gc_content,
        rng=rng,
    )
    
    sequences = {f"phylo_seq_{i:04d}": seq for i, seq in enumerate(sequences_list)}
    fasta_file = output_dir / "phylogeny_test.fasta"
    write_fasta(sequences, str(fasta_file))
    
    logger.info(f"Phylogeny test sequences saved to {fasta_file}")
    
    return {
        "type": "phylogeny",
        "n_sequences": n_sequences,
        "length": sequence_length,
        "diversity": diversity,
        "output_file": str(fasta_file),
    }


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="DNA sequence simulation for testing and validation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate basic sequences
  %(prog)s --type sequences --n 100 --length 1000 --output output/simulation/dna

  # Generate population with diversity
  %(prog)s --type population --n 50 --length 2000 --diversity 0.01

  # Generate alignment test data
  %(prog)s --type alignment --n 10 --length 500 --mutation-rate 0.05

  # Generate phylogeny test data
  %(prog)s --type phylogeny --n 20 --length 1000 --diversity 0.02
        """,
    )
    
    parser.add_argument(
        "--type",
        required=True,
        choices=["sequences", "population", "alignment", "phylogeny"],
        help="Type of simulation to run",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/simulation/dna"),
        help="Output directory (default: output/simulation/dna)",
    )
    parser.add_argument("--n", type=int, default=100, help="Number of sequences")
    parser.add_argument("--length", type=int, default=1000, help="Sequence length")
    parser.add_argument("--gc-content", type=float, default=0.5, help="GC content (0.0-1.0)")
    parser.add_argument("--mutations", type=int, default=0, help="Number of mutations per sequence (sequences type)")
    parser.add_argument("--diversity", type=float, help="Nucleotide diversity (population/phylogeny types)")
    parser.add_argument("--mutation-rate", type=float, default=0.01, help="Mutation rate per site (alignment type)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")
    
    args = parser.parse_args()
    
    if args.verbose:
        import logging as std_logging
        logger.setLevel(std_logging.DEBUG)
    
    # Validate output directory
    output_dir = paths.ensure_directory(args.output)
    
    # Validate common parameters
    validation.validate_range(args.n, min_val=1, name="n")
    validation.validate_range(args.length, min_val=1, name="length")
    validation.validate_range(args.gc_content, min_val=0.0, max_val=1.0, name="gc_content")
    if args.mutations is not None:
        validation.validate_range(args.mutations, min_val=0, name="mutations")
    if args.diversity is not None:
        validation.validate_range(args.diversity, min_val=0.0, name="diversity")
    if hasattr(args, "mutation_rate"):
        validation.validate_range(args.mutation_rate, min_val=0.0, max_val=1.0, name="mutation_rate")
    
    try:
        if args.type == "sequences":
            results = simulate_sequences(
                output_dir, args.n, args.length, args.gc_content, args.mutations, args.seed
            )
        elif args.type == "population":
            results = simulate_population(
                output_dir, args.n, args.length, args.diversity, args.gc_content, args.seed
            )
        elif args.type == "alignment":
            results = simulate_alignment(
                output_dir, args.n, args.length, args.gc_content, args.mutation_rate, args.seed
            )
        elif args.type == "phylogeny":
            results = simulate_phylogeny(
                output_dir, args.n, args.length, args.gc_content, args.diversity or 0.01, args.seed
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

