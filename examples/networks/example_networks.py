#!/usr/bin/env python3
"""Network analysis example.

This example demonstrates biological network construction and analysis using METAINFORMANT's networks toolkit.

Usage:
    python examples/networks/example_networks.py

Output:
    output/examples/networks/network_analysis.json
"""

from __future__ import annotations

from pathlib import Path
from metainformant.core import io
from metainformant.networks.graph import create_network, centrality_measures

def main():
    """Demonstrate network analysis."""
    # Setup output directory
    output_dir = Path("output/examples/networks")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== METAINFORMANT Networks Example ===")

    # Create sample biological network (protein-protein interactions)
    edges = [
        ("Protein_A", "Protein_B"),
        ("Protein_A", "Protein_C"),
        ("Protein_B", "Protein_D"),
        ("Protein_C", "Protein_D"),
        ("Protein_C", "Protein_E"),
        ("Protein_D", "Protein_F")
    ]

    print(f"✓ Creating network with {len(edges)} interactions")

    # Build network - extract unique nodes from edges
    nodes = set()
    for edge in edges:
        nodes.add(edge[0])
        nodes.add(edge[1])
    network = create_network(list(nodes), directed=False)

    # Add edges to network
    for node1, node2 in edges:
        network.add_edge(node1, node2)

    # Calculate centrality measures
    all_centrality = centrality_measures(network)
    centrality = all_centrality.get("degree", {})

    analysis_results = {
        "network_stats": {
            "nodes": network.num_nodes(),
            "edges": network.num_edges(),
            "average_degree": sum(len(network.get_neighbors(node)) for node in network.nodes) / network.num_nodes()
        },
        "centrality_measures": centrality
    }

    print(f"Network: {analysis_results['network_stats']['nodes']} nodes, {analysis_results['network_stats']['edges']} edges")

    # Save results
    results_file = output_dir / "network_analysis.json"
    io.dump_json({
        "network_analysis": analysis_results
    }, results_file, indent=2)

    print(f"✓ Results saved to: {results_file}")
    print("\n=== Networks Example Complete ===")

if __name__ == "__main__":
    main()
