#!/usr/bin/env python3
"""Protein sequence analysis example.

This example demonstrates protein sequence processing and analysis using METAINFORMANT's protein toolkit.

Usage:
    python examples/protein/example_sequences.py

Output:
    output/examples/protein/sequence_analysis.json
"""

from __future__ import annotations

from pathlib import Path
from metainformant.core import io

def main():
    """Demonstrate protein sequence analysis."""
    # Setup output directory
    output_dir = Path("output/examples/protein")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== METAINFORMANT Protein Sequences Example ===")

    # Create sample protein sequences
    sequences = {
        "protein_A": "MAKDVKFGNDPLAVDKGHLVQIYLGKGLTVLGGNDLRQFGGGK",
        "protein_B": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
        "protein_C": "MKVLWAALLVTFLAGCQAKVEQAVETEPEPELRQQTEWQSGQRWELALGRFWDYLRWVQTLSEQVQEELLSSQVTQELRALMDETMKELKAYKSELEEQLTPVAEETRARLSKELQAAQARLGADMEDVCGRLVQYRGEVQAMLGQSTEELRVRLASHLRKLRKRLLRDADDLQKRLAVYQAGAREGAERGLSAIRERLGPLVEQGRVRAATVGSLAGQPLQERAQAWGERLRARAEQDAAKGYR",
        "protein_D": "MGDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKGIIWGEDTLMEYLENPKKYIPGTKMIFAGLKKEVGEHGLGHVTRAFSDGLAHLDNLKGTFAQLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
    }

    print(f"✓ Analyzing {len(sequences)} protein sequences")

    # Basic sequence analysis
    analysis_results = {}

    for seq_id, sequence in sequences.items():
        result = {
            "sequence": sequence,
            "length": len(sequence),
            "molecular_weight": len(sequence) * 110,  # Rough estimate
            "hydrophobic_residues": sum(1 for aa in sequence if aa in "AILMFWYV"),
            "charged_residues": sum(1 for aa in sequence if aa in "DEKRH"),
            "hydrophobic_ratio": sum(1 for aa in sequence if aa in "AILMFWYV") / len(sequence),
            "isoelectric_point_estimate": 6.0  # Simplified
        }
        analysis_results[seq_id] = result

        print(f"  {seq_id}: {result['length']} aa, {result['hydrophobic_residues']} hydrophobic residues")

    # Save results
    results_file = output_dir / "sequence_analysis.json"
    io.dump_json({
        "protein_sequence_analysis": {
            "sequences_analyzed": len(sequences),
            "results": analysis_results
        }
    }, results_file, indent=2)

    print(f"✓ Results saved to: {results_file}")
    print("\n=== Protein Sequences Example Complete ===")

if __name__ == "__main__":
    main()
