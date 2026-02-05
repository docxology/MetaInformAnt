#!/usr/bin/env python3
"""DNA {{analysis_type}} Analysis Example

This example demonstrates:
- Loading DNA sequences from FASTA files
- Performing {{analysis_type.lower()}} analysis
- Visualizing results

Usage:
    python examples/dna/example_{{name}}.py

Expected output:
    output/examples/dna/{{name}}_results.json
"""

from __future__ import annotations

from pathlib import Path

from metainformant.core import io
from metainformant.dna import sequences

# Sample DNA sequences for demonstration
SAMPLE_SEQUENCES = {"seq1": "ATCGATCGATCGATCG", "seq2": "GCTAGCTAGCTAGCTA", "seq3": "ATATATATATATATAT"}


def main():
    """Main function demonstrating DNA {{analysis_type.lower()}} analysis."""
    print("DNA {{analysis_type}} Analysis Example")
    print("=" * 40)

    # Create output directory
    output_dir = Path("output/examples/dna")
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "{{name}}_results.json"

    try:
        results_data = {}

        # Example {{analysis_type.lower()}} analysis
        for seq_id, sequence in SAMPLE_SEQUENCES.items():
            print(f"Analyzing {seq_id}...")

            # Basic sequence properties
            results_data[seq_id] = {
                "sequence": sequence,
                "length": len(sequence),
                "gc_content": sequences.gc_content(sequence),
            }

            # {{analysis_type}} specific analysis
            {{analysis_code}}

        # Save results
        results = {
            "example": "{{name}}",
            "domain": "dna",
            "analysis_type": "{{analysis_type}}",
            "description": "DNA {{analysis_type}} analysis example",
            "results": results_data,
        }

        io.dump_json(results, output_file)
        print(f"✅ Results saved to: {output_file}")

    except Exception as e:
        print(f"❌ Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
