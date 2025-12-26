#!/usr/bin/env python3
"""DNA sequence alignment example.

This example demonstrates sequence alignment techniques using METAINFORMANT's DNA alignment toolkit, including global and local alignment algorithms.

Usage:
    python examples/dna/example_alignment.py

Output:
    output/examples/dna/alignment_results.json
"""

from __future__ import annotations

from pathlib import Path
from metainformant.core import io
from metainformant.dna import alignment

def main():
    """Demonstrate DNA sequence alignment."""
    # Setup output directory
    output_dir = Path("output/examples/dna")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== METAINFORMANT DNA Alignment Example ===")

    # 1. Define sample sequences for alignment
    print("\n1. Preparing sample sequences for alignment...")

    # Related sequences with mutations
    reference_seq = "ATCGATCGATCGATCGATCG"  # 20bp reference
    similar_seq = "ATCGATCGATCGATCGATCA"   # One base difference at end
    different_seq = "GCTAGCTAGCTAGCTAGCTA"  # Completely different sequence
    mutated_seq = "ATCGATCGATCGATCGTTTG"   # Multiple mutations

    # Longer sequences for more interesting alignments
    seq_a = "ATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    seq_b = "ATCGATCGATCGATCGATCGATCGATCGATCGATCA"  # One SNP at end
    seq_c = "ATCGATCGATCGATCGATCGATCGATCGATCGTTTT"  # Multiple changes

    test_sequences = {
        "reference": reference_seq,
        "similar": similar_seq,
        "different": different_seq,
        "mutated": mutated_seq,
        "seq_a": seq_a,
        "seq_b": seq_b,
        "seq_c": seq_c
    }

    print(f"✓ Prepared {len(test_sequences)} test sequences for alignment")

    # 2. Global alignment examples
    print("\n2. Performing global alignments...")

    global_alignments = {}

    # Default parameters: match=1, mismatch=-1, gap=-2
    alignment_pairs = [
        ("reference", "similar"),
        ("reference", "mutated"),
        ("seq_a", "seq_b"),
        ("seq_a", "seq_c")
    ]

    for seq1_name, seq2_name in alignment_pairs:
        seq1 = test_sequences[seq1_name]
        seq2 = test_sequences[seq2_name]

        # Perform global alignment
        alignment_result = alignment.global_align(seq1, seq2)

        global_alignments[f"{seq1_name}_vs_{seq2_name}"] = {
            "sequence1": seq1,
            "sequence2": seq2,
            "aligned_seq1": alignment_result.aligned_seq1,
            "aligned_seq2": alignment_result.aligned_seq2,
            "score": alignment_result.score,
            "identity": alignment.calculate_alignment_identity(alignment_result),
            "conserved_regions": alignment.find_conserved_regions(alignment_result),
            "statistics": alignment.alignment_statistics(alignment_result)
        }

        print(f"  {seq1_name} vs {seq2_name}:")
        print(f"    Score: {alignment_result.score}")
        print(f"    Identity: {global_alignments[f'{seq1_name}_vs_{seq2_name}']['identity']:.1%}")
        print(f"    Conserved regions: {len(global_alignments[f'{seq1_name}_vs_{seq2_name}']['conserved_regions'])}")

    # 3. Local alignment examples
    print("\n3. Performing local alignments...")

    local_alignments = {}

    # Local alignment finds best matching subsequences
    local_pairs = [
        ("reference", "different"),  # Very different sequences
        ("seq_a", "seq_c"),
        ("reference", "seq_a")  # Subsequence matching
    ]

    for seq1_name, seq2_name in local_pairs:
        seq1 = test_sequences[seq1_name]
        seq2 = test_sequences[seq2_name]

        alignment_result = alignment.local_align(seq1, seq2)

        local_alignments[f"{seq1_name}_vs_{seq2_name}"] = {
            "sequence1": seq1,
            "sequence2": seq2,
            "aligned_seq1": alignment_result.aligned_seq1,
            "aligned_seq2": alignment_result.aligned_seq2,
            "score": alignment_result.score,
            "identity": alignment.calculate_alignment_identity(alignment_result),
            "conserved_regions": alignment.find_conserved_regions(alignment_result),
            "statistics": alignment.alignment_statistics(alignment_result)
        }

        print(f"  {seq1_name} vs {seq2_name} (local):")
        print(f"    Score: {alignment_result.score}")
        print(f"    Identity: {local_alignments[f'{seq1_name}_vs_{seq2_name}']['identity']:.1%}")

    # 4. Custom alignment parameters
    print("\n4. Testing custom alignment parameters...")

    custom_alignments = {}

    # Test different scoring schemes
    scoring_schemes = [
        {"name": "strict", "match": 1, "mismatch": -3, "gap": -5},
        {"name": "lenient", "match": 2, "mismatch": -1, "gap": -1},
        {"name": "nucleotide_bias", "match": 1, "mismatch": -2, "gap": -3}
    ]

    test_seq1 = test_sequences["reference"]
    test_seq2 = test_sequences["mutated"]

    for scheme in scoring_schemes:
        name = scheme["name"]
        result = alignment.global_align(
            test_seq1, test_seq2,
            match_score=scheme["match"],
            mismatch_score=scheme["mismatch"],
            gap_score=scheme["gap"]
        )

        custom_alignments[name] = {
            "scoring_scheme": scheme,
            "aligned_seq1": result.aligned_seq1,
            "aligned_seq2": result.aligned_seq2,
            "score": result.score,
            "identity": alignment.calculate_alignment_identity(result),
            "statistics": alignment.alignment_statistics(result)
        }

        print(f"  {name} scoring: score={result.score}, identity={custom_alignments[name]['identity']:.1%}")

    # 5. Conserved region analysis
    print("\n5. Analyzing conserved regions...")

    conserved_analysis = {}

    for alignment_name, align_data in global_alignments.items():
        conserved_regions = align_data["conserved_regions"]

        # Analyze each conserved region
        region_details = []
        for region_seq, start_pos, end_pos in conserved_regions:
            region_info = {
                "sequence": region_seq,
                "start_position": start_pos,
                "end_position": end_pos,
                "length": len(region_seq),
                "gc_content": region_seq.count('G') + region_seq.count('C')
            }
            region_details.append(region_info)

        conserved_analysis[alignment_name] = {
            "total_regions": len(conserved_regions),
            "regions": region_details,
            "total_conserved_bases": sum(len(r["sequence"]) for r in region_details),
            "average_region_length": sum(len(r["sequence"]) for r in region_details) / len(region_details) if region_details else 0
        }

        if conserved_regions:
            print(f"  {alignment_name}: {len(conserved_regions)} conserved regions, {conserved_analysis[alignment_name]['total_conserved_bases']} conserved bases")

    # 6. Alignment statistics comparison
    print("\n6. Comparing alignment statistics...")

    stats_comparison = {}

    all_alignments = {
        "global": global_alignments,
        "local": local_alignments,
        "custom": {f"reference_vs_mutated_{k}": v for k, v in custom_alignments.items()}
    }

    for alignment_type, alignments in all_alignments.items():
        type_stats = []

        for align_name, align_data in alignments.items():
            stats = align_data["statistics"]
            stats["alignment_name"] = align_name
            stats["identity"] = align_data["identity"]
            type_stats.append(stats)

        stats_comparison[alignment_type] = type_stats

        if type_stats:
            avg_identity = sum(s["identity"] for s in type_stats) / len(type_stats)
            print(f"  {alignment_type} alignments: average identity = {avg_identity:.1%}")

    # 7. Create comprehensive alignment results
    print("\n7. Creating comprehensive alignment results...")

    alignment_results = {
        "dna_alignment_analysis": {
            "timestamp": "2024-12-26T10:00:00Z",
            "test_sequences": test_sequences,
            "alignment_algorithms": {
                "global_alignment": {
                    "description": "Needleman-Wunsch algorithm for full sequence alignment",
                    "use_case": "Comparing complete sequences of similar length",
                    "default_parameters": {"match": 1, "mismatch": -1, "gap": -2}
                },
                "local_alignment": {
                    "description": "Smith-Waterman algorithm for subsequence alignment",
                    "use_case": "Finding conserved regions in dissimilar sequences",
                    "default_parameters": {"match": 1, "mismatch": -1, "gap": -2}
                }
            },
            "results": {
                "global_alignments": global_alignments,
                "local_alignments": local_alignments,
                "custom_parameter_alignments": custom_alignments,
                "conserved_region_analysis": conserved_analysis,
                "statistics_comparison": stats_comparison
            },
            "summary_metrics": {
                "total_sequences_tested": len(test_sequences),
                "total_alignments_performed": len(global_alignments) + len(local_alignments) + len(custom_alignments),
                "alignment_types": ["global", "local", "custom_scoring"],
                "average_global_identity": sum(a["identity"] for a in global_alignments.values()) / len(global_alignments) if global_alignments else 0,
                "total_conserved_regions_found": sum(len(a["conserved_regions"]) for a in global_alignments.values()),
                "custom_scoring_schemes_tested": len(scoring_schemes)
            },
            "functions_demonstrated": [
                "global_align", "local_align", "calculate_alignment_identity",
                "find_conserved_regions", "alignment_statistics"
            ],
            "key_insights": [
                "Global alignment works best for similar sequences of equal length",
                "Local alignment finds conserved regions in dissimilar sequences",
                "Scoring parameters significantly affect alignment results",
                "Conserved region analysis helps identify functional motifs",
                "Alignment statistics provide quantitative comparison measures"
            ]
        }
    }

    results_file = output_dir / "alignment_results.json"
    io.dump_json(alignment_results, results_file, indent=2)

    print(f"✓ Comprehensive alignment analysis saved to: {results_file}")

    print("\n=== DNA Alignment Example Complete ===")
    print("This example demonstrated METAINFORMANT's DNA sequence alignment capabilities:")
    print("- Global vs local alignment algorithms")
    print("- Custom scoring parameters and their effects")
    print("- Alignment identity and statistics calculation")
    print("- Conserved region identification and analysis")
    print("- Comparative analysis of different alignment approaches")

    print(f"\nAll outputs saved to: {output_dir}")

if __name__ == "__main__":
    main()
