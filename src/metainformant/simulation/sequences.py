"""Sequence simulation utilities for generating synthetic biological sequences.

This module provides functions for generating random DNA and protein sequences,
introducing mutations, and simulating sequence evolution over generations.
All functions support reproducible results through random seed control.
"""

from __future__ import annotations

import random
from typing import Dict, List, Optional, Tuple, Any
from collections import Counter

from metainformant.core import logging, validation, errors

logger = logging.get_logger(__name__)

# Standard genetic code and amino acids
DNA_BASES = "ATCG"
RNA_BASES = "AUCG"
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"  # Standard 20 amino acids

# Genetic code mapping (simplified)
GENETIC_CODE: Dict[str, str] = {
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
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}


def generate_random_dna(length: int, *, gc_content: float = 0.5, rng: random.Random | None = None) -> str:
    """Generate a random DNA sequence with specified GC content.

    Args:
        length: Length of the sequence to generate
        gc_content: Target GC content (0.0 to 1.0)
        rng: Random number generator for reproducibility

    Returns:
        Random DNA sequence as string

    Raises:
        ValueError: If length <= 0 or gc_content not in [0, 1]
    """
    validation.validate_range(length, min_val=1, name="length")
    validation.validate_range(gc_content, min_val=0.0, max_val=1.0, name="gc_content")

    if rng is None:
        rng = random.Random()

    # Calculate base frequencies
    gc_freq = gc_content
    at_freq = 1.0 - gc_content

    # Split GC and AT evenly
    g_freq = c_freq = gc_freq / 2
    a_freq = t_freq = at_freq / 2

    bases = ['A', 'T', 'G', 'C']
    weights = [a_freq, t_freq, g_freq, c_freq]

    sequence = rng.choices(bases, weights=weights, k=length)
    return ''.join(sequence)


def mutate_sequence(seq: str, n_mut: int, *, rng: random.Random | None = None) -> str:
    """Introduce point mutations into a sequence.

    Args:
        seq: Input sequence
        n_mut: Number of mutations to introduce
        rng: Random number generator for reproducibility

    Returns:
        Mutated sequence

    Raises:
        ValueError: If n_mut < 0 or n_mut > len(seq)
        TypeError: If seq is not a string
    """
    validation.validate_type(seq, str, "seq")
    validation.validate_range(n_mut, min_val=0, name="n_mut")
    validation.validate_range(n_mut, min_val=0, max_val=len(seq), name="n_mut")

    if rng is None:
        rng = random.Random()

    if n_mut == 0:
        return seq

    sequence = list(seq)
    positions = rng.sample(range(len(sequence)), n_mut)

    for pos in positions:
        current_base = sequence[pos]

        if current_base in DNA_BASES:
            # DNA mutation
            possible_mutations = [b for b in DNA_BASES if b != current_base]
        elif current_base in AMINO_ACIDS:
            # Protein mutation - mutate to any other amino acid
            possible_mutations = [aa for aa in AMINO_ACIDS if aa != current_base]
        else:
            # Unknown character - replace with random DNA base
            possible_mutations = list(DNA_BASES)

        sequence[pos] = rng.choice(possible_mutations)

    return ''.join(sequence)


def generate_random_protein(length: int, *, rng: random.Random | None = None) -> str:
    """Generate a random protein sequence.

    Args:
        length: Length of the protein sequence to generate
        rng: Random number generator for reproducibility

    Returns:
        Random protein sequence as string

    Raises:
        ValueError: If length <= 0
    """
    validation.validate_range(length, min_val=1, name="length")

    if rng is None:
        rng = random.Random()

    sequence = rng.choices(list(AMINO_ACIDS), k=length)
    return ''.join(sequence)


def evolve_sequence(sequence: str, generations: int, *,
                   mutation_rate: float = 0.001, rng: random.Random | None = None) -> str:
    """Evolve a sequence over multiple generations with mutations.

    Args:
        sequence: Initial sequence
        generations: Number of generations to evolve
        mutation_rate: Probability of mutation per base per generation
        rng: Random number generator for reproducibility

    Returns:
        Evolved sequence

    Raises:
        ValueError: If generations < 0 or mutation_rate not in [0, 1]
        TypeError: If sequence is not a string
    """
    validation.validate_type(sequence, str, "sequence")
    validation.validate_range(generations, min_val=0, name="generations")
    validation.validate_range(mutation_rate, min_val=0.0, max_val=1.0, name="mutation_rate")

    if rng is None:
        rng = random.Random()

    current_seq = sequence

    for gen in range(generations):
        # Calculate expected number of mutations
        expected_mutations = len(current_seq) * mutation_rate

        # Use Poisson distribution for actual number of mutations
        n_mutations = rng.poisson(expected_mutations)

        if n_mutations > 0:
            current_seq = mutate_sequence(current_seq, min(n_mutations, len(current_seq)), rng=rng)

    return current_seq


def translate_dna_to_protein(dna_sequence: str, *, frame: int = 0) -> str:
    """Translate a DNA sequence to protein using the genetic code.

    Args:
        dna_sequence: DNA sequence to translate
        frame: Reading frame (0, 1, or 2)

    Returns:
        Translated protein sequence

    Raises:
        ValueError: If frame not in [0, 1, 2] or sequence invalid
    """
    validation.validate_range(frame, min_val=0, max_val=2, name="frame")
    validation.validate_type(dna_sequence, str, "dna_sequence")

    # Validate DNA sequence
    if not all(base in DNA_BASES for base in dna_sequence.upper()):
        raise errors.ValidationError(f"Invalid DNA sequence: contains non-DNA bases")

    sequence = dna_sequence.upper()[frame:]
    protein = []

    for i in range(0, len(sequence) - len(sequence) % 3, 3):
        codon = sequence[i:i+3]
        if len(codon) == 3:
            amino_acid = GENETIC_CODE.get(codon, 'X')  # X for unknown codons
            protein.append(amino_acid)
            if amino_acid == '*':  # Stop codon
                break

    return ''.join(protein)


def reverse_transcribe_protein_to_dna(protein_sequence: str, *,
                                    rng: random.Random | None = None) -> str:
    """Reverse translate a protein sequence to DNA (one possible codon usage).

    Args:
        protein_sequence: Protein sequence to reverse translate
        rng: Random number generator for reproducibility

    Returns:
        DNA sequence encoding the protein

    Raises:
        TypeError: If protein_sequence is not a string
    """
    validation.validate_type(protein_sequence, str, "protein_sequence")

    if rng is None:
        rng = random.Random()

    dna_sequence = []

    for aa in protein_sequence.upper():
        if aa == '*':
            # Stop codon - choose randomly
            stop_codons = ['TAA', 'TAG', 'TGA']
            codon = rng.choice(stop_codons)
        elif aa in AMINO_ACIDS:
            # Find all codons for this amino acid
            possible_codons = [codon for codon, amino in GENETIC_CODE.items() if amino == aa]
            codon = rng.choice(possible_codons) if possible_codons else 'NNN'
        else:
            # Unknown amino acid
            codon = 'NNN'

        dna_sequence.append(codon)

    return ''.join(dna_sequence)


def generate_coding_sequence(length: int, *, gc_content: float = 0.5,
                           rng: random.Random | None = None) -> Tuple[str, str]:
    """Generate a random coding DNA sequence and its protein translation.

    Args:
        length: Length of the DNA sequence (will be multiple of 3)
        gc_content: Target GC content for DNA
        rng: Random number generator for reproducibility

    Returns:
        Tuple of (dna_sequence, protein_sequence)

    Raises:
        ValueError: If length not divisible by 3
    """
    if length % 3 != 0:
        raise errors.ValidationError(f"DNA coding sequence length must be divisible by 3, got {length}")

    dna_seq = generate_random_dna(length, gc_content=gc_content, rng=rng)
    protein_seq = translate_dna_to_protein(dna_seq)

    return dna_seq, protein_seq


def calculate_sequence_similarity(seq1: str, seq2: str) -> float:
    """Calculate sequence similarity (fraction of identical positions).

    Args:
        seq1: First sequence
        seq2: Second sequence

    Returns:
        Similarity score (0.0 to 1.0)

    Raises:
        TypeError: If sequences are not strings
        ValueError: If sequences have different lengths
    """
    validation.validate_type(seq1, str, "seq1")
    validation.validate_type(seq2, str, "seq2")

    if len(seq1) != len(seq2):
        raise errors.ValidationError(f"Sequences must have equal length: {len(seq1)} vs {len(seq2)}")

    if len(seq1) == 0:
        return 1.0

    matches = sum(a == b for a, b in zip(seq1, seq2))
    return matches / len(seq1)


def generate_sequence_family(ancestor: str, n_descendants: int, generations: int, *,
                           mutation_rate: float = 0.001, rng: random.Random | None = None) -> List[str]:
    """Generate a family of related sequences from a common ancestor.

    Args:
        ancestor: Ancestral sequence
        n_descendants: Number of descendant sequences to generate
        generations: Number of generations of evolution
        mutation_rate: Mutation rate per base per generation
        rng: Random number generator for reproducibility

    Returns:
        List of descendant sequences (including ancestor as first element)

    Raises:
        ValueError: If n_descendants < 1 or generations < 0
    """
    validation.validate_type(ancestor, str, "ancestor")
    validation.validate_range(n_descendants, min_val=1, name="n_descendants")
    validation.validate_range(generations, min_val=0, name="generations")

    if rng is None:
        rng = random.Random()

    descendants = [ancestor]  # Include ancestor

    for _ in range(n_descendants):
        descendant = evolve_sequence(ancestor, generations,
                                   mutation_rate=mutation_rate, rng=rng)
        descendants.append(descendant)

    return descendants


def analyze_sequence_divergence(sequences: List[str]) -> Dict[str, Any]:
    """Analyze sequence divergence in a set of related sequences.

    Args:
        sequences: List of sequences to analyze

    Returns:
        Dictionary with divergence statistics

    Raises:
        ValueError: If fewer than 2 sequences provided
    """
    validation.validate_type(sequences, list, "sequences")
    if len(sequences) < 2:
        raise errors.ValidationError(f"At least 2 sequences required for divergence analysis, got {len(sequences)}")

    # Check all sequences have same length
    lengths = [len(seq) for seq in sequences]
    if len(set(lengths)) != 1:
        raise errors.ValidationError(f"All sequences must have equal length, got lengths: {lengths}")

    n_seq = len(sequences)
    similarities = []

    # Calculate pairwise similarities
    for i in range(n_seq):
        for j in range(i + 1, n_seq):
            sim = calculate_sequence_similarity(sequences[i], sequences[j])
            similarities.append(sim)

    # Calculate statistics
    mean_similarity = sum(similarities) / len(similarities)
    divergences = [1.0 - sim for sim in similarities]
    mean_divergence = sum(divergences) / len(divergences)

    # Count variable positions
    seq_length = lengths[0]
    variable_positions = 0

    for pos in range(seq_length):
        bases_at_pos = [seq[pos] for seq in sequences]
        if len(set(bases_at_pos)) > 1:
            variable_positions += 1

    return {
        "num_sequences": n_seq,
        "sequence_length": seq_length,
        "mean_similarity": mean_similarity,
        "mean_divergence": mean_divergence,
        "variable_positions": variable_positions,
        "variable_fraction": variable_positions / seq_length if seq_length > 0 else 0,
        "pairwise_similarities": similarities,
        "pairwise_divergences": divergences,
    }


def simulate_gene_duplication(original_gene: str, n_copies: int, *,
                            divergence_time: int = 1000, mutation_rate: float = 1e-8,
                            rng: random.Random | None = None) -> List[str]:
    """Simulate gene duplication and subsequent divergence.

    Args:
        original_gene: Original gene sequence
        n_copies: Number of duplicate copies to generate
        divergence_time: Generations since duplication
        mutation_rate: Mutation rate per base per generation
        rng: Random number generator for reproducibility

    Returns:
        List of diverged gene copies

    Raises:
        ValueError: If n_copies < 1 or divergence_time < 0
    """
    validation.validate_type(original_gene, str, "original_gene")
    validation.validate_range(n_copies, min_val=1, name="n_copies")
    validation.validate_range(divergence_time, min_val=0, name="divergence_time")

    if rng is None:
        rng = random.Random()

    copies = []
    for _ in range(n_copies):
        # Each copy evolves independently from the original
        copy = evolve_sequence(original_gene, divergence_time,
                             mutation_rate=mutation_rate, rng=rng)
        copies.append(copy)

    return copies





