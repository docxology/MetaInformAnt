"""Contamination detection and analysis for biological data."""

from __future__ import annotations

import re
from collections import Counter
from typing import Dict, List, Set, Tuple

import numpy as np
from Bio import SeqIO


def detect_cross_species_contamination(
    sequences: List[str],
    reference_genomes: Dict[str, str],
    threshold: float = 0.1,
) -> Dict[str, List[str]]:
    """
    Detect potential cross-species contamination by comparing sequences
    against multiple reference genomes.

    Args:
        sequences: List of DNA sequences to check
        reference_genomes: Dict mapping species names to reference sequences
        threshold: Minimum similarity threshold for contamination detection

    Returns:
        Dict mapping sequence indices to list of potential contaminant species
    """
    contamination_results = {}

    for i, seq in enumerate(sequences):
        potential_contaminants = []

        for species, ref_seq in reference_genomes.items():
            # Simple similarity check (in production, use more sophisticated alignment)
            similarity = _calculate_sequence_similarity(seq, ref_seq)
            if similarity > threshold:
                potential_contaminants.append(f"{species} ({similarity:.2f})")

        if potential_contaminants:
            contamination_results[str(i)] = potential_contaminants

    return contamination_results


def detect_rrna_contamination(
    sequences: List[str],
    rrna_patterns: List[str] | None = None,
) -> Dict[str, float]:
    """
    Detect ribosomal RNA contamination in sequencing data.

    Args:
        sequences: List of sequences to check
        rrna_patterns: List of rRNA sequence patterns (auto-generated if None)

    Returns:
        Dict mapping sequence indices to contamination scores
    """
    if rrna_patterns is None:
        # Common rRNA patterns (simplified)
        rrna_patterns = [
            "GGAAGGAGCAGTG",
            "GGAAAGAGCAGTG",
            "AGGAAGGAGCAGT",
            "AAGGAAGGAGCAG",
        ]

    contamination_scores = {}

    for i, seq in enumerate(sequences):
        score = 0.0

        for pattern in rrna_patterns:
            if pattern in seq:
                score += 1.0

        # Normalize by sequence length
        if len(seq) > 0:
            score = score / len(seq) * 1000  # Scale for readability

        if score > 0:
            contamination_scores[str(i)] = score

    return contamination_scores


def detect_mycoplasma_contamination(
    sequences: List[str],
    mycoplasma_genome: str | None = None,
) -> Dict[str, bool]:
    """
    Detect mycoplasma contamination in cell culture samples.

    Args:
        sequences: List of sequences to check
        mycoplasma_genome: Mycoplasma reference genome sequence

    Returns:
        Dict mapping sequence indices to contamination status
    """
    if mycoplasma_genome is None:
        # Simplified mycoplasma detection pattern
        mycoplasma_pattern = "TTAAATTTAAATTTAAATTT"

    contamination_results = {}

    for i, seq in enumerate(sequences):
        if mycoplasma_genome:
            similarity = _calculate_sequence_similarity(seq, mycoplasma_genome)
            is_contaminated = similarity > 0.8
        else:
            # Pattern-based detection
            is_contaminated = mycoplasma_pattern in seq

        if is_contaminated:
            contamination_results[str(i)] = True

    return contamination_results


def detect_adapter_contamination(
    sequences: List[str],
    adapter_sequences: List[str] | None = None,
) -> Dict[str, List[str]]:
    """
    Detect adapter sequence contamination.

    Args:
        sequences: List of sequences to check
        adapter_sequences: List of known adapter sequences

    Returns:
        Dict mapping sequence indices to list of detected adapters
    """
    if adapter_sequences is None:
        # Common Illumina adapters
        adapter_sequences = [
            "AGATCGGAAGAG",  # TruSeq Universal Adapter
            "TGGAATTCTCGG",  # Small RNA 3' Adapter
            "CTGTCTCTTATA",  # Nextera Transposase
        ]

    contamination_results = {}

    for i, seq in enumerate(sequences):
        detected_adapters = []

        for adapter in adapter_sequences:
            if adapter in seq:
                detected_adapters.append(adapter)

        if detected_adapters:
            contamination_results[str(i)] = detected_adapters

    return contamination_results


def detect_vector_contamination(
    sequences: List[str],
    vector_sequences: Dict[str, str] | None = None,
) -> Dict[str, List[str]]:
    """
    Detect plasmid or viral vector contamination.

    Args:
        sequences: List of sequences to check
        vector_sequences: Dict mapping vector names to sequences

    Returns:
        Dict mapping sequence indices to list of potential vector contaminants
    """
    if vector_sequences is None:
        # Common vector elements (simplified)
        vector_sequences = {
            "pUC19": "GGCCGCTCTAGAACTAGTGGATC",
            "pBR322": "CCTGCAGGTCGACTCTAGAGGAT",
            "pET": "CATATGCGGTGTGAAATACCGC",
        }

    contamination_results = {}

    for i, seq in enumerate(sequences):
        potential_vectors = []

        for vector_name, vector_seq in vector_sequences.items():
            similarity = _calculate_sequence_similarity(seq, vector_seq)
            if similarity > 0.7:
                potential_vectors.append(f"{vector_name} ({similarity:.2f})")

        if potential_vectors:
            contamination_results[str(i)] = potential_vectors

    return contamination_results


def _calculate_sequence_similarity(seq1: str, seq2: str) -> float:
    """Calculate simple sequence similarity (Jaccard index of k-mers)."""
    if not seq1 or not seq2:
        return 0.0

    k = min(6, len(seq1), len(seq2))  # Use 6-mers or shorter if sequences are small

    def get_kmers(seq: str, k: int) -> Set[str]:
        kmers = set()
        for i in range(len(seq) - k + 1):
            kmers.add(seq[i:i+k])
        return kmers

    kmers1 = get_kmers(seq1, k)
    kmers2 = get_kmers(seq2, k)

    if not kmers1 or not kmers2:
        return 0.0

    intersection = len(kmers1.intersection(kmers2))
    union = len(kmers1.union(kmers2))

    return intersection / union if union > 0 else 0.0


def generate_contamination_report(
    contamination_results: Dict[str, Dict[str, List[str]]],
    sample_metadata: Dict[str, str] | None = None,
) -> str:
    """
    Generate a comprehensive contamination analysis report.

    Args:
        contamination_results: Results from various contamination detection methods
        sample_metadata: Optional metadata about samples

    Returns:
        Formatted contamination report
    """
    report_lines = []
    report_lines.append("METAINFORMANT Contamination Analysis Report")
    report_lines.append("=" * 50)
    report_lines.append("")

    total_samples = 0
    contaminated_samples = 0

    for method, results in contamination_results.items():
        if results:
            report_lines.append(f"{method.upper()} Contamination:")
            for sample_id, contaminants in results.items():
                total_samples += 1
                if contaminants:
                    contaminated_samples += 1
                    report_lines.append(f"  Sample {sample_id}: {', '.join(contaminants)}")
            report_lines.append("")

    report_lines.append("Summary:")
    report_lines.append(f"  Total samples analyzed: {total_samples}")
    report_lines.append(f"  Samples with contamination: {contaminated_samples}")

    if total_samples > 0:
        contamination_rate = (contaminated_samples / total_samples) * 100
        report_lines.append(f"  Contamination rate: {contamination_rate:.1f}%")

    return "\n".join(report_lines)
