"""DNA mutation detection and analysis utilities.

This module provides tools for detecting, classifying, and analyzing mutations
in DNA sequences, including point mutations, insertions, deletions, and
evolutionary rate calculations.
"""

from __future__ import annotations

import random
from typing import Dict, List, Tuple

from metainformant.core import logging

logger = logging.get_logger(__name__)


def calculate_mutation_rate(ancestral: str, derived: str) -> float:
    """Calculate mutation rate between ancestral and derived sequences.

    Args:
        ancestral: Ancestral DNA sequence
        derived: Derived DNA sequence

    Returns:
        Mutation rate (mutations per site)

    Raises:
        ValueError: If sequences have different lengths

    Example:
        >>> ancestral = "ATCG"
        >>> derived = "ATCG"
        >>> rate = calculate_mutation_rate(ancestral, derived)
        >>> rate == 0.0
        True
    """
    if len(ancestral) != len(derived):
        raise ValueError("Sequences must have equal length")

    if not ancestral:
        return 0.0

    mutations = 0
    for a, d in zip(ancestral.upper(), derived.upper()):
        if a != d:
            mutations += 1

    return mutations / len(ancestral)


def classify_mutations(ancestral: str, derived: str) -> Dict[str, int]:
    """Classify mutations by type between ancestral and derived sequences.

    Args:
        ancestral: Ancestral DNA sequence
        derived: Derived DNA sequence

    Returns:
        Dictionary with mutation counts by type

    Raises:
        ValueError: If sequences have different lengths

    Example:
        >>> ancestral = "ATCG"
        >>> derived = "AGCT"
        >>> mutations = classify_mutations(ancestral, derived)
        >>> "transitions" in mutations
        True
    """
    if len(ancestral) != len(derived):
        raise ValueError("Sequences must have equal length")

    mutation_types = {
        'transitions': 0,
        'transversions': 0,
        'synonymous': 0,
        'nonsynonymous': 0,
        'total': 0
    }

    # Transition pairs: A<->G, C<->T
    transition_pairs = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}

    # Codon translation table (simplified)
    genetic_code = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'G', 'GCC': 'G', 'GCA': 'G', 'GCG': 'G',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

    for i, (a, d) in enumerate(zip(ancestral.upper(), derived.upper())):
        if a != d:
            mutation_types['total'] += 1

            # Classify transition vs transversion
            if (a, d) in transition_pairs:
                mutation_types['transitions'] += 1
            else:
                mutation_types['transversions'] += 1

            # Check if this affects protein coding (simplified)
            # Look at codon context
            codon_start = (i // 3) * 3
            if codon_start + 3 <= len(ancestral):
                ancestral_codon = ancestral[codon_start:codon_start+3].upper()
                derived_codon = ancestral_codon[:i % 3] + d + ancestral_codon[i % 3 + 1:]

                ancestral_aa = genetic_code.get(ancestral_codon, 'X')
                derived_aa = genetic_code.get(derived_codon, 'X')

                if ancestral_aa == derived_aa:
                    mutation_types['synonymous'] += 1
                else:
                    mutation_types['nonsynonymous'] += 1

    return mutation_types


def generate_point_mutations(sequence: str, num_mutations: int,
                           mutation_rate: float = 0.001) -> str:
    """Generate point mutations in a DNA sequence.

    Args:
        sequence: Original DNA sequence
        num_mutations: Number of mutations to introduce
        mutation_rate: Probability of mutation at each site

    Returns:
        Mutated DNA sequence

    Example:
        >>> seq = "ATCGATCG"
        >>> mutated = generate_point_mutations(seq, num_mutations=2)
        >>> len(mutated) == len(seq)
        True
    """
    if not sequence:
        return sequence

    nucleotides = ['A', 'C', 'G', 'T']
    mutated = list(sequence.upper())

    mutations_applied = 0

    for i in range(len(mutated)):
        if mutations_applied >= num_mutations:
            break

        if random.random() < mutation_rate:
            current = mutated[i]
            # Choose different nucleotide
            possible = [n for n in nucleotides if n != current]
            if possible:
                mutated[i] = random.choice(possible)
                mutations_applied += 1

    return ''.join(mutated)


def simulate_sequence_evolution(sequence: str, generations: int,
                              mutation_rate: float = 0.001) -> List[str]:
    """Simulate sequence evolution under mutation pressure.

    Args:
        sequence: Starting DNA sequence
        generations: Number of generations to simulate
        mutation_rate: Mutation rate per site per generation

    Returns:
        List of sequences for each generation

    Example:
        >>> seq = "ATCG"
        >>> evolution = simulate_sequence_evolution(seq, generations=3)
        >>> len(evolution) == 4  # Including original
        True
    """
    if not sequence or generations < 0:
        return [sequence]

    evolutionary_lineage = [sequence]

    for gen in range(generations):
        # Mutate the current sequence
        current_seq = evolutionary_lineage[-1]
        mutated_seq = generate_point_mutations(current_seq, num_mutations=1,
                                             mutation_rate=mutation_rate)
        evolutionary_lineage.append(mutated_seq)

    return evolutionary_lineage


def find_mutational_hotspots(sequence: str, window_size: int = 10) -> List[Tuple[int, float]]:
    """Identify potential mutational hotspots in a sequence.

    Args:
        sequence: DNA sequence to analyze
        window_size: Size of sliding window for analysis

    Returns:
        List of tuples (position, hotspot_score)

    Example:
        >>> seq = "ATCGATCGATCG"
        >>> hotspots = find_mutational_hotspots(seq, window_size=4)
        >>> isinstance(hotspots, list)
        True
    """
    if not sequence or window_size <= 0:
        return []

    hotspots = []

    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size].upper()

        # Calculate GC content (simple hotspot predictor)
        gc_count = window.count('G') + window.count('C')
        gc_content = gc_count / window_size

        # Calculate repeat content
        repeat_score = 0
        for j in range(window_size - 1):
            if window[j] == window[j + 1]:
                repeat_score += 1

        # Combined hotspot score
        hotspot_score = (gc_content * 0.6) + (repeat_score / window_size * 0.4)

        hotspots.append((i, hotspot_score))

    return hotspots


def calculate_substitution_matrix(seq1: str, seq2: str) -> Dict[Tuple[str, str], int]:
    """Calculate nucleotide substitution matrix between two sequences.

    Args:
        seq1: First DNA sequence
        seq2: Second DNA sequence

    Returns:
        Dictionary mapping (from_nuc, to_nuc) pairs to counts

    Example:
        >>> seq1 = "ATCG"
        >>> seq2 = "AGCT"
        >>> matrix = calculate_substitution_matrix(seq1, seq2)
        >>> ("T", "G") in matrix
        True
    """
    if len(seq1) != len(seq2):
        return {}

    substitutions = {}

    for a, b in zip(seq1.upper(), seq2.upper()):
        if a != b:
            pair = (a, b)
            substitutions[pair] = substitutions.get(pair, 0) + 1

    return substitutions


def detect_selection_signatures(sequence: str, reference_sequences: List[str]) -> Dict[str, float]:
    """Detect signatures of natural selection in a sequence.

    Args:
        sequence: Query sequence
        reference_sequences: List of reference sequences

    Returns:
        Dictionary with selection metrics

    Example:
        >>> query = "ATCGATCG"
        >>> refs = ["ATCGATCG", "ATCGATCG"]
        >>> signatures = detect_selection_signatures(query, refs)
        >>> "diversity" in signatures
        True
    """
    if not reference_sequences:
        return {'diversity': 0.0, 'conservation': 1.0}

    metrics = {
        'diversity': 0.0,
        'conservation': 0.0,
        'synonymous_sites': 0,
        'nonsynonymous_sites': 0
    }

    # Simple diversity calculation
    diversities = []
    for ref in reference_sequences:
        if len(ref) == len(sequence):
            div = calculate_mutation_rate(sequence, ref)
            diversities.append(div)

    if diversities:
        metrics['diversity'] = sum(diversities) / len(diversities)

        # Conservation is inverse of diversity
        metrics['conservation'] = 1.0 - metrics['diversity']

    return metrics


def generate_mutant_library(sequence: str, positions: List[int],
                          mutations: List[str]) -> List[str]:
    """Generate a library of mutant sequences.

    Args:
        sequence: Original DNA sequence
        positions: List of positions to mutate
        mutations: List of mutation strings (one per position)

    Returns:
        List of mutant sequences

    Example:
        >>> seq = "ATCG"
        >>> mutants = generate_mutant_library(seq, [1, 3], ["G", "T"])
        >>> len(mutants) == 2
        True
    """
    if len(positions) != len(mutations):
        raise ValueError("Positions and mutations lists must have equal length")

    mutants = []

    for pos, mutation in zip(positions, mutations):
        if 0 <= pos < len(sequence):
            mutant = list(sequence)
            mutant[pos] = mutation.upper()
            mutants.append(''.join(mutant))

    return mutants


def analyze_mutation_spectrum(sequence: str, reference: str) -> Dict[str, any]:
    """Analyze the spectrum of mutations between sequences.

    Args:
        sequence: Derived sequence
        reference: Reference sequence

    Returns:
        Dictionary with mutation spectrum analysis

    Example:
        >>> ref = "ATCGATCG"
        >>> seq = "AGCGATCG"
        >>> spectrum = analyze_mutation_spectrum(seq, ref)
        >>> "mutation_types" in spectrum
        True
    """
    if len(sequence) != len(reference):
        raise ValueError("Sequences must have equal length")

    mutation_types = classify_mutations(reference, sequence)

    # Calculate mutation density
    mutation_density = mutation_types['total'] / len(sequence) if sequence else 0

    # Calculate transition/transversion ratio
    ti_tv_ratio = (mutation_types['transitions'] / mutation_types['transversions']
                  if mutation_types['transversions'] > 0 else float('inf'))

    spectrum = {
        'mutation_types': mutation_types,
        'mutation_density': mutation_density,
        'ti_tv_ratio': ti_tv_ratio,
        'sequence_length': len(sequence),
        'substitution_matrix': calculate_substitution_matrix(reference, sequence)
    }

    return spectrum








