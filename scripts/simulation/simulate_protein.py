#!/usr/bin/env python3
"""Protein sequence simulation script.

This script generates synthetic protein sequences, structure coordinates,
domain annotations, and PPI interaction data for testing and validation.

Usage:
    python3 scripts/simulation/simulate_protein.py --type sequences --n 100 --length 200
    python3 scripts/simulation/simulate_protein.py --type structure --n 10 --length 150
    python3 scripts/simulation/simulate_protein.py --type domains --n 50 --length 300
    python3 scripts/simulation/simulate_protein.py --type ppi --n 100 --interactions 200
"""

import argparse
import random
import sys
from pathlib import Path

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core import io, logging, paths, validation
from metainformant.protein.sequences import write_fasta
from metainformant.simulation.sequences import generate_random_protein, mutate_sequence

logger = logging.get_logger(__name__)


def simulate_sequences(
    output_dir: Path,
    n_sequences: int,
    sequence_length: int,
    mutations: int,
    seed: int,
) -> dict:
    """Simulate basic protein sequences.
    
    Args:
        output_dir: Output directory for results
        n_sequences: Number of sequences to generate
        sequence_length: Length of each sequence
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
    validation.validate_range(mutations, min_val=0, name="mutations")
    
    logger.info(f"Generating {n_sequences} protein sequences of length {sequence_length}")
    rng = random.Random(seed)
    
    sequences = {}
    for i in range(n_sequences):
        seq = generate_random_protein(sequence_length, rng=rng)
        if mutations > 0:
            seq = mutate_sequence(seq, n_mut=mutations, rng=rng)
        sequences[f"protein_{i:04d}"] = seq
    
    fasta_file = output_dir / "protein_sequences.fasta"
    write_fasta(sequences, str(fasta_file))
    logger.info(f"Protein sequences saved to {fasta_file}")
    
    return {
        "type": "sequences",
        "n_sequences": n_sequences,
        "length": sequence_length,
        "mutations_per_seq": mutations,
        "output_file": str(fasta_file),
    }


def simulate_structure(
    output_dir: Path,
    n_sequences: int,
    sequence_length: int,
    seed: int,
) -> dict:
    """Simulate protein structure coordinates (CA atoms)."""
    logger.info(f"Generating structure coordinates for {n_sequences} proteins")
    rng = random.Random(seed)
    
    structures = {}
    sequences = {}
    
    for i in range(n_sequences):
        seq = generate_random_protein(sequence_length, rng=rng)
        protein_id = f"protein_{i:04d}"
        sequences[protein_id] = seq
        
        # Generate CA coordinates (simplified - random walk)
        coords = []
        x, y, z = 0.0, 0.0, 0.0
        for _ in range(len(seq)):
            x += rng.gauss(0, 3.8)  # ~3.8A CA-CA distance
            y += rng.gauss(0, 3.8)
            z += rng.gauss(0, 3.8)
            coords.append({"x": x, "y": y, "z": z})
        
        structures[protein_id] = {
            "sequence": seq,
            "ca_coordinates": coords,
        }
    
    # Save sequences
    fasta_file = output_dir / "structure_sequences.fasta"
    write_fasta(sequences, str(fasta_file))
    
    # Save structures
    structure_file = output_dir / "structures.json"
    io.dump_json(structures, structure_file, indent=2)
    
    logger.info(f"Structure data saved to {structure_file}")
    
    return {
        "type": "structure",
        "n_sequences": n_sequences,
        "length": sequence_length,
        "output_file": str(structure_file),
        "sequences_file": str(fasta_file),
    }


def simulate_domains(
    output_dir: Path,
    n_sequences: int,
    sequence_length: int,
    n_domains_per_seq: int,
    seed: int,
) -> dict:
    """Simulate protein domain annotations."""
    logger.info(f"Generating domain annotations for {n_sequences} proteins")
    rng = random.Random(seed)
    
    domain_types = ["SH2", "SH3", "PH", "PDZ", "kinase", "helix", "sheet", "coiled_coil"]
    sequences = {}
    annotations = []
    
    for i in range(n_sequences):
        seq = generate_random_protein(sequence_length, rng=rng)
        protein_id = f"protein_{i:04d}"
        sequences[protein_id] = seq
        
        # Generate random domain annotations
        n_domains = rng.randint(1, n_domains_per_seq)
        domains = []
        used_positions = set()
        
        for _ in range(n_domains):
            domain_type = rng.choice(domain_types)
            domain_length = rng.randint(20, min(100, sequence_length // 2))
            start = rng.randint(0, sequence_length - domain_length)
            
            # Avoid overlapping domains
            if any(pos in used_positions for pos in range(start, start + domain_length)):
                continue
            
            used_positions.update(range(start, start + domain_length))
            domains.append({
                "domain_type": domain_type,
                "start": start,
                "end": start + domain_length,
            })
        
        annotations.append({
            "protein_id": protein_id,
            "n_domains": len(domains),
            "domains": domains,
        })
    
    # Save sequences
    fasta_file = output_dir / "domain_sequences.fasta"
    write_fasta(sequences, str(fasta_file))
    
    # Save annotations
    annotation_file = output_dir / "domain_annotations.json"
    io.dump_json(annotations, annotation_file, indent=2)
    
    logger.info(f"Domain annotations saved to {annotation_file}")
    
    return {
        "type": "domains",
        "n_sequences": n_sequences,
        "output_file": str(annotation_file),
        "sequences_file": str(fasta_file),
    }


def simulate_ppi(
    output_dir: Path,
    n_proteins: int,
    n_interactions: int,
    seed: int,
) -> dict:
    """Simulate protein-protein interaction network."""
    logger.info(f"Generating PPI network: {n_proteins} proteins, {n_interactions} interactions")
    rng = random.Random(seed)
    
    # Generate protein IDs
    proteins = [f"protein_{i:04d}" for i in range(n_proteins)]
    
    # Generate interactions
    interactions = []
    interaction_set = set()
    
    while len(interactions) < n_interactions:
        p1 = rng.choice(proteins)
        p2 = rng.choice(proteins)
        
        if p1 == p2:
            continue
        
        # Ensure undirected edge representation
        edge = tuple(sorted([p1, p2]))
        if edge in interaction_set:
            continue
        
        interaction_set.add(edge)
        interactions.append({
            "protein1": p1,
            "protein2": p2,
            "confidence": rng.random(),
            "interaction_type": rng.choice(["physical", "genetic", "regulatory"]),
        })
    
    # Save as edge list
    import pandas as pd
    df = pd.DataFrame(interactions)
    edge_file = output_dir / "ppi_edges.csv"
    df.to_csv(edge_file, index=False)
    
    logger.info(f"PPI network saved to {edge_file}")
    
    return {
        "type": "ppi",
        "n_proteins": n_proteins,
        "n_interactions": len(interactions),
        "output_file": str(edge_file),
    }


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Protein simulation for testing and validation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate protein sequences
  %(prog)s --type sequences --n 100 --length 200

  # Generate structure coordinates
  %(prog)s --type structure --n 10 --length 150

  # Generate domain annotations
  %(prog)s --type domains --n 50 --length 300 --n-domains 3

  # Generate PPI network
  %(prog)s --type ppi --n 100 --interactions 200
        """,
    )
    
    parser.add_argument(
        "--type",
        required=True,
        choices=["sequences", "structure", "domains", "ppi"],
        help="Type of simulation to run",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/simulation/protein"),
        help="Output directory (default: output/simulation/protein)",
    )
    parser.add_argument("--n", type=int, default=100, help="Number of proteins/sequences")
    parser.add_argument("--length", type=int, default=200, help="Sequence length (sequences/structure/domains)")
    parser.add_argument("--mutations", type=int, default=0, help="Number of mutations per sequence (sequences type)")
    parser.add_argument("--n-domains", type=int, default=2, help="Max domains per sequence (domains type)")
    parser.add_argument("--interactions", type=int, default=200, help="Number of interactions (ppi type)")
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
    if hasattr(args, "length"):
        validation.validate_range(args.length, min_val=1, name="length")
    if hasattr(args, "mutations"):
        validation.validate_range(args.mutations, min_val=0, name="mutations")
    if hasattr(args, "n_domains"):
        validation.validate_range(args.n_domains, min_val=1, name="n_domains")
    if hasattr(args, "interactions"):
        validation.validate_range(args.interactions, min_val=0, name="interactions")
    
    try:
        if args.type == "sequences":
            results = simulate_sequences(output_dir, args.n, args.length, args.mutations, args.seed)
        elif args.type == "structure":
            results = simulate_structure(output_dir, args.n, args.length, args.seed)
        elif args.type == "domains":
            results = simulate_domains(output_dir, args.n, args.length, args.n_domains, args.seed)
        elif args.type == "ppi":
            results = simulate_ppi(output_dir, args.n, args.interactions, args.seed)
        
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

