#!/usr/bin/env python3
"""DNA phylogenetic analysis example.

This example demonstrates phylogenetic tree construction using METAINFORMANT's phylogeny toolkit, including neighbor-joining and UPGMA methods.

Usage:
    python examples/dna/example_phylogeny.py

Output:
    output/examples/dna/phylogeny_tree.newick
"""

from __future__ import annotations

from pathlib import Path
from metainformant.core import io
from metainformant.dna import phylogeny

def main():
    """Demonstrate DNA phylogenetic analysis."""
    # Setup output directory
    output_dir = Path("output/examples/dna")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== METAINFORMANT DNA Phylogeny Example ===")

    # 1. Create evolutionary related sequences
    print("\n1. Creating evolutionary related sequences...")

    # Simulate evolutionary relationships
    # Root sequence
    ancestral_seq = "ATCGATCGATCGATCGATCG"

    # Simulated mutations to create related sequences
    sequences = {
        "ancestral": ancestral_seq,
        "species_a": "ATCGATCGATCGATCGATCA",  # One mutation
        "species_b": "ATCGATCGATCGATCGATCG",  # Identical to ancestral
        "species_c": "ATCGATCGATCGATCGTTTG",  # Two mutations
        "species_d": "GTCGATCGATCGATCGATCG",  # One mutation at start
        "species_e": "ATCGATCGATCGATCGATCG",  # Identical to ancestral
        "species_f": "ATCGATCGATCGATCGAAAA",  # Multiple mutations at end
        "outgroup":  "GCTAGCTAGCTAGCTAGCTA"   # Very different (outgroup)
    }

    print(f"✓ Created {len(sequences)} sequences for phylogenetic analysis")

    # 2. Build neighbor-joining tree
    print("\n2. Building neighbor-joining tree...")

    nj_tree = phylogeny.neighbor_joining_tree(sequences)
    nj_newick = phylogeny.to_newick(nj_tree)
    nj_stats = phylogeny.basic_tree_stats(nj_tree)

    print("✓ Neighbor-joining tree constructed")
    print(f"  Tree statistics: {nj_stats}")

    # Save Newick format
    nj_file = output_dir / "neighbor_joining_tree.newick"
    with open(nj_file, "w") as f:
        f.write(nj_newick)

    print(f"✓ Newick tree saved to: {nj_file}")

    # 3. Build UPGMA tree
    print("\n3. Building UPGMA tree...")

    upgma_tree = phylogeny.upgma_tree(sequences)
    upgma_newick = phylogeny.to_newick(upgma_tree)
    upgma_stats = phylogeny.basic_tree_stats(upgma_tree)

    print("✓ UPGMA tree constructed")
    print(f"  Tree statistics: {upgma_stats}")

    # Save Newick format
    upgma_file = output_dir / "upgma_tree.newick"
    with open(upgma_file, "w") as f:
        f.write(upgma_newick)

    print(f"✓ Newick tree saved to: {upgma_file}")

    # 4. K-mer based tree construction
    print("\n4. Building k-mer based phylogenetic tree...")

    # Test different k values
    k_values = [2, 3, 4, 5]
    kmer_trees = {}

    for k in k_values:
        try:
            kmer_tree = phylogeny.nj_tree_from_kmer(sequences, k=k, metric="cosine")
            kmer_newick = phylogeny.to_newick(kmer_tree)
            kmer_stats = phylogeny.basic_tree_stats(kmer_tree)

            kmer_trees[f"k{k}"] = {
                "newick": kmer_newick,
                "stats": kmer_stats
            }

            print(f"✓ K-mer tree (k={k}) constructed: {kmer_stats}")

            # Save Newick format
            kmer_file = output_dir / f"kmer_tree_k{k}.newick"
            with open(kmer_file, "w") as f:
                f.write(kmer_newick)

        except Exception as e:
            print(f"✗ Failed to build k-mer tree (k={k}): {e}")
            kmer_trees[f"k{k}"] = {"error": str(e)}

    # 5. ASCII tree visualization
    print("\n5. ASCII tree visualization...")

    try:
        nj_ascii = phylogeny.to_ascii(nj_tree)
        print("Neighbor-joining tree:")
        print(nj_ascii)
        print()
    except Exception as e:
        print(f"ASCII visualization failed: {e}")
        nj_ascii = None

    try:
        upgma_ascii = phylogeny.to_ascii(upgma_tree)
        print("UPGMA tree:")
        print(upgma_ascii)
        print()
    except Exception as e:
        print(f"ASCII visualization failed: {e}")
        upgma_ascii = None

    # 6. Bootstrap support analysis
    print("\n6. Bootstrap support analysis...")

    # Note: Bootstrap analysis with small dataset may not be meaningful
    # but demonstrates the API
    try:
        bootstrap_tree = phylogeny.bootstrap_support(
            nj_tree, sequences, n_replicates=10, method="nj"
        )
        bootstrap_newick = phylogeny.to_newick(bootstrap_tree)

        bootstrap_file = output_dir / "bootstrap_tree.newick"
        with open(bootstrap_file, "w") as f:
            f.write(bootstrap_newick)

        print("✓ Bootstrap support analysis completed (10 replicates)")
        print(f"  Bootstrap tree saved to: {bootstrap_file}")

    except Exception as e:
        print(f"Bootstrap analysis failed: {e}")
        bootstrap_newick = None

    # 7. Compare tree topologies
    print("\n7. Comparing tree topologies...")

    tree_comparison = {
        "neighbor_joining": {
            "newick_length": len(nj_newick),
            "stats": nj_stats,
            "topology_summary": "Distance-based method assuming molecular clock"
        },
        "upgma": {
            "newick_length": len(upgma_newick),
            "stats": upgma_stats,
            "topology_summary": "Assumes constant molecular clock rate"
        }
    }

    # Compare Newick strings to see topological differences
    if nj_newick != upgma_newick:
        tree_comparison["topology_difference"] = "Trees have different topologies"
    else:
        tree_comparison["topology_difference"] = "Trees have identical topologies"

    print(f"  NJ vs UPGMA: {tree_comparison['topology_difference']}")

    # 8. Sequence distance analysis
    print("\n8. Analyzing sequence relationships...")

    # Calculate pairwise distances (simple Hamming distance for demonstration)
    sequence_distances = {}

    seq_ids = list(sequences.keys())
    for i, id1 in enumerate(seq_ids):
        for id2 in seq_ids[i+1:]:
            seq1 = sequences[id1]
            seq2 = sequences[id2]

            # Simple distance calculation
            distance = sum(a != b for a, b in zip(seq1, seq2))
            normalized_distance = distance / len(seq1)

            sequence_distances[f"{id1}_vs_{id2}"] = {
                "distance": distance,
                "normalized_distance": normalized_distance,
                "similarity": 1 - normalized_distance
            }

    # Find most similar and dissimilar pairs
    most_similar = min(sequence_distances.items(), key=lambda x: x[1]["distance"])
    most_dissimilar = max(sequence_distances.items(), key=lambda x: x[1]["distance"])

    relationship_analysis = {
        "most_similar_pair": {
            "pair": most_similar[0],
            "distance": most_similar[1]["distance"],
            "similarity": most_similar[1]["similarity"]
        },
        "most_dissimilar_pair": {
            "pair": most_dissimilar[0],
            "distance": most_dissimilar[1]["distance"],
            "similarity": most_dissimilar[1]["similarity"]
        },
        "all_distances": sequence_distances
    }

    print(f"  Most similar: {most_similar[0]} (distance: {most_similar[1]['distance']})")
    print(f"  Most dissimilar: {most_dissimilar[0]} (distance: {most_dissimilar[1]['distance']})")

    # 9. Create comprehensive phylogeny results
    print("\n9. Creating comprehensive phylogeny results...")

    phylogeny_results = {
        "dna_phylogeny_analysis": {
            "timestamp": "2024-12-26T10:00:00Z",
            "input_sequences": sequences,
            "tree_construction_methods": {
                "neighbor_joining": {
                    "description": "Distance-based method that doesn't assume molecular clock",
                    "algorithm": "Saitou and Nei (1987)",
                    "use_case": "General phylogenetic reconstruction"
                },
                "upgma": {
                    "description": "Assumes constant rate of evolution (molecular clock)",
                    "algorithm": "Sokal and Michener (1958)",
                    "use_case": "When evolutionary rates are expected to be constant"
                },
                "kmer_based": {
                    "description": "Uses k-mer frequencies for distance calculation",
                    "algorithm": "Neighbor-joining with k-mer distances",
                    "use_case": "Large datasets where full alignment is expensive"
                }
            },
            "results": {
                "neighbor_joining": {
                    "newick_string": nj_newick,
                    "statistics": nj_stats,
                    "ascii_visualization": nj_ascii
                },
                "upgma": {
                    "newick_string": upgma_newick,
                    "statistics": upgma_stats,
                    "ascii_visualization": upgma_ascii
                },
                "kmer_trees": kmer_trees,
                "bootstrap_analysis": {
                    "performed": bootstrap_newick is not None,
                    "replicates": 10 if bootstrap_newick else 0,
                    "newick_string": bootstrap_newick
                },
                "tree_comparison": tree_comparison,
                "sequence_relationships": relationship_analysis
            },
            "summary_metrics": {
                "sequences_analyzed": len(sequences),
                "tree_methods_tested": 3,  # NJ, UPGMA, k-mer
                "k_values_tested": len(k_values),
                "bootstrap_replicates": 10,
                "output_files_generated": [
                    "neighbor_joining_tree.newick",
                    "upgma_tree.newick",
                    "bootstrap_tree.newick"
                ] + [f"kmer_tree_k{k}.newick" for k in k_values if f"k{k}" in kmer_trees and "error" not in kmer_trees[f"k{k}"]],
                "most_similar_sequences": most_similar[0],
                "most_dissimilar_sequences": most_dissimilar[0]
            },
            "functions_demonstrated": [
                "neighbor_joining_tree", "upgma_tree", "to_newick",
                "to_ascii", "basic_tree_stats", "bootstrap_support",
                "nj_tree_from_kmer"
            ],
            "key_insights": [
                "Different tree construction methods can yield different topologies",
                "Neighbor-joining is more flexible than UPGMA for varying evolutionary rates",
                "K-mer based methods scale better for large sequence datasets",
                "Bootstrap support helps assess tree reliability",
                "Sequence similarity analysis complements tree reconstruction"
            ]
        }
    }

    results_file = output_dir / "phylogeny_analysis.json"
    io.dump_json(phylogeny_results, results_file, indent=2)

    print(f"✓ Comprehensive phylogeny analysis saved to: {results_file}")

    print("\n=== DNA Phylogeny Example Complete ===")
    print("This example demonstrated METAINFORMANT's phylogenetic analysis capabilities:")
    print("- Neighbor-joining and UPGMA tree construction")
    print("- K-mer based distance methods")
    print("- Newick format output and ASCII visualization")
    print("- Bootstrap support analysis")
    print("- Tree topology comparison and sequence relationship analysis")

    print(f"\nAll outputs saved to: {output_dir}")

if __name__ == "__main__":
    main()
