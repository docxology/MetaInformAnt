#!/usr/bin/env python3
"""Network simulation script.

This script generates synthetic biological networks including PPI networks,
regulatory networks, and pathway networks.

Usage:
    python3 scripts/simulation/simulate_networks.py --type ppi --n-nodes 100 --n-edges 200
    python3 scripts/simulation/simulate_networks.py --type regulatory --n-nodes 50 --n-edges 150
    python3 scripts/simulation/simulate_networks.py --type pathway --n-nodes 30 --n-edges 50
"""

import argparse
import random
import sys
from pathlib import Path

import pandas as pd

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core import io, logging, paths, validation

logger = logging.get_logger(__name__)


def simulate_ppi(
    output_dir: Path,
    n_nodes: int,
    n_edges: int,
    seed: int,
) -> dict:
    """Simulate protein-protein interaction network."""
    logger.info(f"Generating PPI network: {n_nodes} nodes, {n_edges} edges")
    rng = random.Random(seed)
    
    # Generate node IDs
    nodes = [f"protein_{i:04d}" for i in range(n_nodes)]
    
    # Generate edges
    edges = []
    edge_set = set()
    
    while len(edges) < n_edges:
        node1 = rng.choice(nodes)
        node2 = rng.choice(nodes)
        
        if node1 == node2:
            continue
        
        # Undirected edge
        edge = tuple(sorted([node1, node2]))
        if edge in edge_set:
            continue
        
        edge_set.add(edge)
        edges.append({
            "source": node1,
            "target": node2,
            "interaction_type": rng.choice(["physical", "genetic", "regulatory"]),
            "confidence": rng.random(),
            "evidence": rng.choice(["EXP", "HTP", "LTP", "IEA"]),
        })
    
    # Save edge list
    df = pd.DataFrame(edges)
    csv_file = output_dir / "ppi_network.csv"
    df.to_csv(csv_file, index=False)
    
    logger.info(f"PPI network saved to {csv_file}")
    
    return {
        "type": "ppi",
        "n_nodes": n_nodes,
        "n_edges": len(edges),
        "output_file": str(csv_file),
    }


def simulate_regulatory(
    output_dir: Path,
    n_nodes: int,
    n_edges: int,
    seed: int,
) -> dict:
    """Simulate gene regulatory network."""
    logger.info(f"Generating regulatory network: {n_nodes} nodes, {n_edges} edges")
    rng = random.Random(seed)
    
    # Split into TFs and targets
    n_tfs = n_nodes // 3
    tfs = [f"TF_{i:04d}" for i in range(n_tfs)]
    targets = [f"gene_{i:04d}" for i in range(n_tfs, n_nodes)]
    all_nodes = tfs + targets
    
    # Generate regulatory edges (TFs -> targets)
    edges = []
    edge_set = set()
    
    while len(edges) < n_edges:
        tf = rng.choice(tfs)
        target = rng.choice(targets)
        
        edge = (tf, target)
        if edge in edge_set:
            continue
        
        edge_set.add(edge)
        edges.append({
            "source": tf,
            "target": target,
            "regulation_type": rng.choice(["activation", "repression"]),
            "strength": rng.random(),
            "evidence": rng.choice(["ChIP-seq", "expression", "motif"]),
        })
    
    # Save edge list
    df = pd.DataFrame(edges)
    csv_file = output_dir / "regulatory_network.csv"
    df.to_csv(csv_file, index=False)
    
    logger.info(f"Regulatory network saved to {csv_file}")
    
    return {
        "type": "regulatory",
        "n_nodes": n_nodes,
        "n_edges": len(edges),
        "n_tfs": n_tfs,
        "output_file": str(csv_file),
    }


def simulate_pathway(
    output_dir: Path,
    n_nodes: int,
    n_edges: int,
    seed: int,
) -> dict:
    """Simulate pathway network."""
    logger.info(f"Generating pathway network: {n_nodes} nodes, {n_edges} edges")
    rng = random.Random(seed)
    
    # Generate node IDs
    nodes = [f"pathway_node_{i:04d}" for i in range(n_nodes)]
    
    # Pathway networks are often directed and hierarchical
    edges = []
    edge_set = set()
    
    while len(edges) < n_edges:
        source = rng.choice(nodes)
        target = rng.choice(nodes)
        
        if source == target:
            continue
        
        edge = (source, target)
        if edge in edge_set:
            continue
        
        edge_set.add(edge)
        edges.append({
            "source": source,
            "target": target,
            "interaction_type": rng.choice(["phosphorylation", "binding", "activation", "inhibition"]),
            "pathway": rng.choice(["pathway_A", "pathway_B", "pathway_C"]),
        })
    
    # Save edge list
    df = pd.DataFrame(edges)
    csv_file = output_dir / "pathway_network.csv"
    df.to_csv(csv_file, index=False)
    
    logger.info(f"Pathway network saved to {csv_file}")
    
    return {
        "type": "pathway",
        "n_nodes": n_nodes,
        "n_edges": len(edges),
        "output_file": str(csv_file),
    }


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Network simulation for testing and validation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Simulate PPI network
  %(prog)s --type ppi --n-nodes 100 --n-edges 200

  # Simulate regulatory network
  %(prog)s --type regulatory --n-nodes 50 --n-edges 150

  # Simulate pathway network
  %(prog)s --type pathway --n-nodes 30 --n-edges 50
        """,
    )
    
    parser.add_argument(
        "--type",
        required=True,
        choices=["ppi", "regulatory", "pathway"],
        help="Type of simulation to run",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/simulation/networks"),
        help="Output directory (default: output/simulation/networks)",
    )
    parser.add_argument("--n-nodes", type=int, default=100, help="Number of nodes")
    parser.add_argument("--n-edges", type=int, default=200, help="Number of edges")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")
    
    args = parser.parse_args()
    
    if args.verbose:
        import logging as std_logging
        logger.setLevel(std_logging.DEBUG)
    
    output_dir = paths.ensure_directory(args.output)
    
    try:
        if args.type == "ppi":
            results = simulate_ppi(output_dir, args.n_nodes, args.n_edges, args.seed)
        elif args.type == "regulatory":
            results = simulate_regulatory(output_dir, args.n_nodes, args.n_edges, args.seed)
        elif args.type == "pathway":
            results = simulate_pathway(output_dir, args.n_nodes, args.n_edges, args.seed)
        
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

