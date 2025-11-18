from __future__ import annotations

from collections import Counter
from typing import Dict

from Bio import SeqIO


def read_fasta(path: str) -> Dict[str, str]:
    """Read a FASTA file into a dictionary of id -> sequence string."""
    records = SeqIO.parse(path, "fasta")
    seqs: Dict[str, str] = {}
    for rec in records:
        seqs[rec.id] = str(rec.seq)
    return seqs


def reverse_complement(seq: str) -> str:
    """Compute reverse complement of a DNA sequence.
    
    Args:
        seq: DNA sequence string
        
    Returns:
        Reverse complement sequence
    """
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def gc_content(seq: str) -> float:
    """Calculate GC content of a DNA sequence."""
    if not seq:
        return 0.0

    # Convert to uppercase for consistent counting
    upper = seq.upper()

    # Count GC bases efficiently
    gc_count = upper.count("G") + upper.count("C")

    # Handle non-standard characters by counting only ACGT
    total_valid = sum(upper.count(base) for base in "ACGT")

    return gc_count / total_valid if total_valid > 0 else 0.0


def kmer_counts(seq: str, k: int) -> Dict[str, int]:
    """Count k-mers (substrings of length k) in a sequence.
    
    Args:
        seq: DNA sequence string
        k: K-mer length
        
    Returns:
        Dictionary mapping k-mer to count
    """
    if k <= 0 or len(seq) < k:
        return {}
    return dict(Counter(seq[i : i + k] for i in range(0, len(seq) - k + 1)))


def kmer_frequencies(seq: str, k: int) -> Dict[str, float]:
    """Calculate k-mer frequencies (normalized counts).
    
    Args:
        seq: DNA sequence string
        k: K-mer length
        
    Returns:
        Dictionary mapping k-mer to frequency (0.0 to 1.0)
    """
    counts = kmer_counts(seq, k)
    total = sum(counts.values())
    if total == 0:
        return {}
    return {kmer: cnt / total for kmer, cnt in counts.items()}


def sequence_length(seq: str) -> int:
    """Get the length of a DNA sequence, ignoring case and non-standard characters."""
    return len(seq.replace(" ", "").replace("\n", "").replace("\t", ""))


def validate_dna_sequence(seq: str) -> bool:
    """Validate that a string contains only valid DNA characters (ACGTacgtNn-).

    Args:
        seq: DNA sequence to validate

    Returns:
        True if sequence contains only valid DNA characters
    """
    valid_chars = set("ACGTacgtNn-")
    return all(char in valid_chars for char in seq)


def dna_complementarity_score(seq1: str, seq2: str) -> float:
    """Calculate complementarity score between two DNA sequences.

    Args:
        seq1: First DNA sequence
        seq2: Second DNA sequence

    Returns:
        Complementarity score (0-1, where 1 is perfectly complementary)
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length")

    complement_map = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", "a": "t", "t": "a", "c": "g", "g": "c", "n": "n"}

    matches = 0
    total = len(seq1)

    for i in range(total):
        char1, char2 = seq1[i].upper(), seq2[i].upper()
        if char1 in complement_map and complement_map[char1] == char2:
            matches += 1

    return matches / total if total > 0 else 0.0


def find_repeats(seq: str, min_length: int = 3) -> Dict[str, list[int]]:
    """Find repeated substrings in a DNA sequence.

    Args:
        seq: DNA sequence to analyze
        min_length: Minimum length of repeats to find

    Returns:
        Dictionary mapping repeat sequences to their positions
    """
    repeats = {}
    seq_len = len(seq)

    for length in range(min_length, seq_len // 2 + 1):
        for i in range(seq_len - length + 1):
            substring = seq[i:i + length]
            if substring not in repeats:
                positions = []
                for j in range(seq_len - length + 1):
                    if seq[j:j + length] == substring:
                        positions.append(j)
                if len(positions) > 1:
                    repeats[substring] = positions

    return repeats


def find_motifs(seq: str, motif_patterns: list[str]) -> Dict[str, list[int]]:
    """Find multiple motif patterns in a DNA sequence.

    Args:
        seq: DNA sequence to search
        motif_patterns: List of motif patterns (supports IUPAC codes)

    Returns:
        Dictionary mapping motif patterns to their positions
    """
    from .motifs import find_motif_positions

    results = {}
    for pattern in motif_patterns:
        positions = find_motif_positions(seq, pattern)
        if positions:
            results[pattern] = positions

    return results


def calculate_sequence_complexity(seq: str) -> float:
    """Calculate sequence complexity based on k-mer diversity.

    Args:
        seq: DNA sequence to analyze

    Returns:
        Complexity score (0-1, where 1 is maximum diversity)
    """
    if len(seq) < 2:
        return 0.0

    # Use 2-mer diversity as a measure of complexity
    k = 2
    kmers = kmer_counts(seq, k)
    total_kmers = len(seq) - k + 1

    if total_kmers == 0:
        return 0.0

    # Calculate Shannon entropy
    import math
    entropy = 0.0
    for count in kmers.values():
        if count > 0:
            p = count / total_kmers
            entropy -= p * math.log2(p)

    # Normalize by maximum possible entropy
    max_entropy = math.log2(len(kmers)) if kmers else 0
    return entropy / max_entropy if max_entropy > 0 else 0.0


def find_orfs(seq: str, min_length: int = 30) -> list[tuple[int, int, str]]:
    """Find open reading frames in a DNA sequence.

    Args:
        seq: DNA sequence to analyze
        min_length: Minimum ORF length in nucleotides

    Returns:
        List of (start, end, frame) tuples for ORFs
    """
    from .translation import find_orfs as find_translation_orfs

    orfs = find_translation_orfs(seq, min_aa=min_length // 3)
    return [(orf.start, orf.end, orf.frame) for orf in orfs]


def calculate_sequence_entropy(seq: str, k: int = 1) -> float:
    """Calculate Shannon entropy of k-mer frequencies.

    Args:
        seq: DNA sequence to analyze
        k: K-mer size for entropy calculation

    Returns:
        Shannon entropy value
    """
    import math

    kmers = kmer_counts(seq, k)
    total = sum(kmers.values())

    if total == 0:
        return 0.0

    entropy = 0.0
    for count in kmers.values():
        if count > 0:
            p = count / total
            entropy -= p * math.log2(p)

    # Normalize by the maximum possible entropy for the observed k-mers so that
    # the value is in [0, 1] regardless of k. This makes entropy values
    # comparable across different k and ensures that more complex k-mers do not
    # trivially appear to have lower entropy due to fewer observations.
    max_entropy = math.log2(len(kmers)) if kmers else 0.0
    return entropy / max_entropy if max_entropy > 0 else 0.0


def detect_sequence_bias(seq: str) -> Dict[str, float]:
    """Detect nucleotide composition bias in a sequence.

    Args:
        seq: DNA sequence to analyze

    Returns:
        Dictionary with bias statistics
    """
    if not seq:
        return {"gc_content": 0.0, "at_content": 0.0, "purine_content": 0.0, "pyrimidine_content": 0.0}

    upper_seq = seq.upper()
    total = len(upper_seq)

    gc_count = upper_seq.count("G") + upper_seq.count("C")
    at_count = upper_seq.count("A") + upper_seq.count("T")
    purine_count = upper_seq.count("A") + upper_seq.count("G")
    pyrimidine_count = upper_seq.count("C") + upper_seq.count("T")

    return {
        "gc_content": gc_count / total,
        "at_content": at_count / total,
        "purine_content": purine_count / total,
        "pyrimidine_content": pyrimidine_count / total,
        "total_bases": total
    }


def calculate_gc_skew(seq: str) -> float:
    """Calculate GC vs AT skew: (GC - AT) / (GC + AT).
    
    Args:
        seq: DNA sequence to analyze
        
    Returns:
        GC skew value (-1 to 1)
    """
    if not seq:
        return 0.0

    upper_seq = seq.upper()
    gc_count = upper_seq.count("G") + upper_seq.count("C")
    at_count = upper_seq.count("A") + upper_seq.count("T")

    if gc_count + at_count == 0:
        return 0.0

    return (gc_count - at_count) / (gc_count + at_count)


def calculate_at_skew(seq: str) -> float:
    """Calculate AT skew: (A - T) / (A + T).
    
    Args:
        seq: DNA sequence to analyze
        
    Returns:
        AT skew value (-1 to 1)
    """
    if not seq:
        return 0.0
        
    upper_seq = seq.upper()
    a_count = upper_seq.count('A')
    t_count = upper_seq.count('T')
    
    if a_count + t_count == 0:
        return 0.0
        
    return (a_count - t_count) / (a_count + t_count)


def find_palindromes(seq: str, min_length: int = 4) -> list[tuple[str, int, int]]:
    """Find palindromic sequences in DNA.
    
    Args:
        seq: DNA sequence to search
        min_length: Minimum length of palindromes to find
        
    Returns:
        List of (palindrome, start, end) tuples
    """
    palindromes = []
    
    for i in range(len(seq) - min_length + 1):
        for j in range(i + min_length, len(seq) + 1):
            substring = seq[i:j]
            if len(substring) >= min_length:
                # Check if it's a palindrome
                rev_comp = reverse_complement(substring)
                if substring == rev_comp:
                    palindromes.append((substring, i, j))
                    break  # Don't find overlapping palindromes
    
    return palindromes


def calculate_melting_temperature(seq: str, method: str = "wallace") -> float:
    """Calculate DNA melting temperature.
    
    Args:
        seq: DNA sequence
        method: Method to use ('wallace' or 'enhanced')
        
    Returns:
        Melting temperature in Celsius
    """
    if not seq:
        return 0.0
        
    # Count nucleotides
    a_count = seq.upper().count('A')
    t_count = seq.upper().count('T')
    g_count = seq.upper().count('G')
    c_count = seq.upper().count('C')
    
    if method == "wallace":
        # Wallace rule: Tm = 2*(A+T) + 4*(G+C)
        return 2 * (a_count + t_count) + 4 * (g_count + c_count)
    elif method == "enhanced":
        # Enhanced formula for longer sequences
        total = len(seq)
        if total <= 14:
            return 2 * (a_count + t_count) + 4 * (g_count + c_count)
        else:
            return 64.9 + 41 * (g_count + c_count - 16.4) / total
    else:
        raise ValueError(f"Unknown method: {method}")


def calculate_codon_usage(seq: str) -> dict[str, float]:
    """Calculate codon usage frequencies.
    
    Args:
        seq: DNA sequence (must be divisible by 3)
        
    Returns:
        Dictionary mapping codons to their frequencies
    """
    from collections import Counter
    
    if len(seq) % 3 != 0:
        raise ValueError("Sequence length must be divisible by 3")

    codons = [seq[i:i+3].upper() for i in range(0, len(seq), 3)]
    # Exclude stop codons from usage statistics so that codon usage focuses on
    # sense codons only.
    stop_codons = {"TAA", "TAG", "TGA"}
    sense_codons = [c for c in codons if c not in stop_codons]

    if not sense_codons:
        return {}

    codon_counts = Counter(sense_codons)
    total_codons = len(sense_codons)

    return {codon: count / total_codons for codon, count in codon_counts.items()}


def find_start_codons(seq: str) -> list[int]:
    """Find positions of ATG start codons.
    
    Args:
        seq: DNA sequence to search
        
    Returns:
        List of 0-based positions of ATG codons
    """
    positions = []
    for i in range(len(seq) - 2):
        if seq[i:i+3].upper() == 'ATG':
            positions.append(i)
    return positions


def find_stop_codons(seq: str) -> list[int]:
    """Find positions of stop codons (TAA, TAG, TGA).
    
    Args:
        seq: DNA sequence to search
        
    Returns:
        List of 0-based positions of stop codons
    """
    positions = []
    stop_codons = ['TAA', 'TAG', 'TGA']
    
    for i in range(len(seq) - 2):
        codon = seq[i:i+3].upper()
        if codon in stop_codons:
            positions.append(i)
    
    return positions
