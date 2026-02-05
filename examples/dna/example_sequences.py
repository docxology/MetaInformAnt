#!/usr/bin/env python3
"""DNA sequence processing basics example.

This example demonstrates fundamental DNA sequence operations using METAINFORMANT's DNA sequence analysis toolkit.

Usage:
    python examples/dna/example_sequences.py

Output:
    output/examples/dna/sequences_analysis.json
"""

from __future__ import annotations

from pathlib import Path

from metainformant.core import io
from metainformant.dna import sequences


def main():
    """Demonstrate DNA sequence processing basics."""
    # Setup output directory
    output_dir = Path("output/examples/dna")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== METAINFORMANT DNA Sequences Example ===")

    # 1. Create sample FASTA data
    print("\n1. Creating sample FASTA data...")

    sample_sequences = {
        "seq1": "ATCGATCGATCGATCGATCG",  # 20bp, GC=10/20=50%
        "seq2": "GCTAGCTAGCTAGCTAGCTA",  # 20bp, GC=10/20=50%
        "seq3": "TTTTAAAAAAAATTTTTTAA",  # 20bp, GC=0/20=0%
        "seq4": "CGCGCGCGCGCGCGCGCGCG",  # 20bp, GC=20/20=100%
        "seq5": "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",  # 40bp, GC=20/40=50%
    }

    # Write to FASTA file
    fasta_file = output_dir / "sample_sequences.fasta"
    with open(fasta_file, "w") as f:
        for seq_id, sequence in sample_sequences.items():
            f.write(f">{seq_id}\n{sequence}\n")

    print(f"✓ Created sample FASTA file: {fasta_file}")

    # 2. Read FASTA file
    print("\n2. Reading FASTA sequences...")

    loaded_sequences = sequences.read_fasta(fasta_file)
    print(f"✓ Loaded {len(loaded_sequences)} sequences from FASTA file")

    # 3. Basic sequence analysis
    print("\n3. Basic sequence analysis...")

    sequence_analysis = {}

    for seq_id, sequence in loaded_sequences.items():
        analysis = {
            "original_sequence": sequence,
            "length": sequences.sequence_length(sequence),
            "gc_content": sequences.gc_content(sequence),
            "is_valid_dna": sequences.validate_dna_sequence(sequence),
            "reverse_complement": sequences.reverse_complement(sequence),
        }

        # Verify reverse complement is valid
        analysis["reverse_complement_valid"] = sequences.validate_dna_sequence(analysis["reverse_complement"])

        sequence_analysis[seq_id] = analysis

        print(f"  {seq_id}: {analysis['length']}bp, GC={analysis['gc_content']:.1%}")

    # 4. K-mer analysis
    print("\n4. K-mer analysis...")

    kmer_results = {}

    for seq_id, sequence in loaded_sequences.items():
        # Calculate k-mer counts for k=2,3,4
        kmer_analysis = {}

        for k in [2, 3, 4]:
            counts = sequences.kmer_counts(sequence, k)
            frequencies = sequences.kmer_frequencies(sequence, k)

            kmer_analysis[f"k{k}"] = {
                "counts": counts,
                "frequencies": frequencies,
                "unique_kmers": len(counts),
                "most_common": max(counts.items(), key=lambda x: x[1]) if counts else None,
            }

        kmer_results[seq_id] = kmer_analysis

        # Show k=2 results for first sequence
        if seq_id == "seq1":
            print(f"  {seq_id} k-mer analysis:")
            for k, data in kmer_analysis.items():
                most_common = data["most_common"]
                print(f"    {k}: {data['unique_kmers']} unique, most common: {most_common[0]} ({most_common[1]}x)")

    # 5. Sequence complexity and motifs
    print("\n5. Sequence complexity and motif analysis...")

    complexity_results = {}

    for seq_id, sequence in loaded_sequences.items():
        complexity = sequences.calculate_sequence_complexity(sequence)
        entropy = sequences.calculate_sequence_entropy(sequence, k=1)

        # Find simple motifs
        motifs_to_find = ["ATCG", "AAAA", "CGCG"]
        found_motifs = sequences.find_motifs(sequence, motifs_to_find)

        # Find repeats
        repeats = sequences.find_repeats(sequence, min_length=3)

        complexity_results[seq_id] = {
            "complexity_score": complexity,
            "entropy_h1": entropy,
            "motifs_found": found_motifs,
            "repeats_found": repeats,
        }

        print(f"  {seq_id}: complexity={complexity:.3f}, entropy={entropy:.3f}")

        # Show motif results for sequences with motifs
        total_motifs = sum(len(positions) for positions in found_motifs.values())
        if total_motifs > 0:
            print(f"    Found {total_motifs} motif occurrences")

    # 6. ORF and codon analysis
    print("\n6. Open reading frame and codon analysis...")

    orf_results = {}

    for seq_id, sequence in loaded_sequences.items():
        # Find ORFs (minimum length 9 for demonstration)
        orfs = sequences.find_orfs(sequence, min_length=9)

        # Codon usage (requires sequence length divisible by 3)
        codon_usage = {}
        if len(sequence) % 3 == 0:
            codon_usage = sequences.calculate_codon_usage(sequence)

        # Find start and stop codons
        start_positions = sequences.find_start_codons(sequence)
        stop_positions = sequences.find_stop_codons(sequence)

        orf_results[seq_id] = {
            "open_reading_frames": orfs,
            "codon_usage": codon_usage,
            "start_codons_at": start_positions,
            "stop_codons_at": stop_positions,
            "sequence_length_mod3": len(sequence) % 3,
        }

        print(f"  {seq_id}: {len(orfs)} ORFs, {len(start_positions)} start codons, {len(stop_positions)} stop codons")

    # 7. Sequence bias and skew analysis
    print("\n7. Sequence bias and skew analysis...")

    bias_results = {}

    for seq_id, sequence in loaded_sequences.items():
        gc_skew = sequences.calculate_gc_skew(sequence)
        at_skew = sequences.calculate_at_skew(sequence)
        bias = sequences.detect_sequence_bias(sequence)

        # Find palindromes
        palindromes = sequences.find_palindromes(sequence, min_length=4)

        # Calculate melting temperature
        try:
            tm = sequences.calculate_melting_temperature(sequence)
        except Exception:
            tm = None

        bias_results[seq_id] = {
            "gc_skew": gc_skew,
            "at_skew": at_skew,
            "base_bias": bias,
            "palindromes": palindromes,
            "melting_temperature_celsius": tm,
        }

        print(f"  {seq_id}: GC skew={gc_skew:.3f}, AT skew={at_skew:.3f}")

        if tm is not None:
            print(f"    Melting temperature: {tm:.1f}°C")

    # 8. Sequence complementarity
    print("\n8. Sequence complementarity analysis...")

    complementarity_scores = {}

    seq_ids = list(loaded_sequences.keys())
    for i, seq_id1 in enumerate(seq_ids):
        for seq_id2 in seq_ids[i + 1 :]:
            seq1 = loaded_sequences[seq_id1]
            seq2 = loaded_sequences[seq_id2]

            # Only calculate complementarity for sequences of equal length
            if len(seq1) == len(seq2):
                score = sequences.dna_complementarity_score(seq1, seq2)
                complementarity_scores[f"{seq_id1}_vs_{seq_id2}"] = score
            else:
                complementarity_scores[f"{seq_id1}_vs_{seq_id2}"] = None

            # Only show high complementarity scores for comparable sequences
            if score is not None and score > 0.5:
                print(f"  {seq_id1} vs {seq_id2}: complementarity = {score:.3f}")

    # 9. Create comprehensive analysis summary
    print("\n9. Creating comprehensive analysis summary...")

    analysis_summary = {
        "dna_sequences_analysis": {
            "timestamp": "2024-12-26T10:00:00Z",
            "input_file": str(fasta_file.relative_to(output_dir)),
            "sequences_analyzed": len(loaded_sequences),
            "analyses_performed": [
                "sequence_validation",
                "gc_content_calculation",
                "reverse_complement",
                "kmer_analysis",
                "sequence_complexity",
                "motif_detection",
                "orf_finding",
                "codon_usage",
                "sequence_bias",
                "palindrome_detection",
                "melting_temperature",
                "complementarity_scoring",
            ],
            "results": {
                "sequence_analysis": sequence_analysis,
                "kmer_analysis": kmer_results,
                "complexity_analysis": complexity_results,
                "orf_analysis": orf_results,
                "bias_analysis": bias_results,
                "complementarity_scores": complementarity_scores,
            },
            "summary_statistics": {
                "total_sequences": len(loaded_sequences),
                "total_bases": sum(len(seq) for seq in loaded_sequences.values()),
                "average_gc_content": sum(sequence_analysis[s]["gc_content"] for s in sequence_analysis)
                / len(sequence_analysis),
                "sequences_with_orfs": sum(1 for r in orf_results.values() if r["open_reading_frames"]),
                "palindromes_found": sum(len(r["palindromes"]) for r in bias_results.values()),
            },
            "functions_demonstrated": [
                "read_fasta",
                "reverse_complement",
                "gc_content",
                "validate_dna_sequence",
                "kmer_counts",
                "kmer_frequencies",
                "calculate_sequence_complexity",
                "find_motifs",
                "find_orfs",
                "calculate_codon_usage",
                "detect_sequence_bias",
                "calculate_gc_skew",
                "calculate_at_skew",
                "find_palindromes",
                "calculate_melting_temperature",
                "dna_complementarity_score",
            ],
        }
    }

    results_file = output_dir / "sequences_analysis.json"
    io.dump_json(analysis_summary, results_file, indent=2)

    print(f"✓ Comprehensive analysis saved to: {results_file}")

    print("\n=== DNA Sequences Example Complete ===")
    print("This example demonstrated METAINFORMANT's comprehensive DNA sequence analysis capabilities:")
    print("- Sequence I/O and validation")
    print("- Basic sequence statistics (GC content, length)")
    print("- K-mer analysis and complexity measures")
    print("- Motif and repeat detection")
    print("- ORF finding and codon analysis")
    print("- Sequence bias and skew calculations")
    print("- Melting temperature and complementarity")

    print(f"\nAll outputs saved to: {output_dir}")


if __name__ == "__main__":
    main()
