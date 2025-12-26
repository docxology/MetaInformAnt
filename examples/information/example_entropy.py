#!/usr/bin/env python3
"""Information theory entropy analysis example.

This example demonstrates information entropy analysis using METAINFORMANT's information toolkit.

Usage:
    python examples/information/example_entropy.py

Output:
    output/examples/information/entropy_analysis.json
"""

from __future__ import annotations

from pathlib import Path
from metainformant.core import io
from metainformant.information.syntactic import shannon_entropy, shannon_entropy_from_counts

def main():
    """Demonstrate entropy analysis."""
    # Setup output directory
    output_dir = Path("output/examples/information")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== METAINFORMANT Information Theory Example ===")

    # Test sequences with different complexities
    sequences = {
        "uniform": "ATCGATCGATCG",  # Equal base frequencies
        "repetitive": "AAAAAAAAAAAA",  # Low complexity
        "complex": "ATCGCGTATAGC",  # Higher complexity
        "mixed": "ATGCGTACGTTAG"  # Mixed complexity
    }

    entropy_results = {}

    for seq_name, sequence in sequences.items():
        # Calculate Shannon entropy for different k-mer sizes
        entropies = {}
        for k in [1, 2, 3]:
            # Convert sequence to list of symbols and count frequencies
            symbols = [sequence[i:i+k] for i in range(len(sequence)-k+1)]
            from collections import Counter
            symbol_counts = Counter(symbols)
            entropy = shannon_entropy_from_counts(symbol_counts)
            entropies[f"k{k}"] = entropy

        entropy_results[seq_name] = {
            "sequence": sequence,
            "length": len(sequence),
            "entropies": entropies,
            "complexity_rank": "high" if entropies["k1"] > 1.8 else "medium" if entropies["k1"] > 1.5 else "low"
        }

        print(f"  {seq_name}: H₁ = {entropies['k1']:.3f} ({entropy_results[seq_name]['complexity_rank']} complexity)")

    # Save results
    results_file = output_dir / "entropy_analysis.json"
    io.dump_json({
        "information_analysis": {
            "sequences_analyzed": len(sequences),
            "entropy_measures": entropy_results
        }
    }, results_file, indent=2)

    print(f"✓ Results saved to: {results_file}")
    print("\n=== Information Theory Example Complete ===")

if __name__ == "__main__":
    main()
