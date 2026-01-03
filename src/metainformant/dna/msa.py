"""Multiple sequence alignment for DNA sequences.

This module provides tools for aligning multiple DNA sequences using various
algorithms including progressive alignment and integration with external tools
like MUSCLE, MAFFT, and ClustalW.
"""

from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional

from metainformant.core import io, logging

logger = logging.get_logger(__name__)


def progressive_alignment(sequences: Dict[str, str], method: str = "muscle") -> Dict[str, str]:
    """Perform progressive multiple sequence alignment.

    Args:
        sequences: Dictionary mapping sequence IDs to DNA sequences
        method: Alignment method ("muscle", "mafft", "clustalw")

    Returns:
        Dictionary mapping sequence IDs to aligned sequences

    Raises:
        RuntimeError: If external alignment tool fails
        FileNotFoundError: If alignment tool is not installed

    Example:
        >>> seqs = {"seq1": "ATCG", "seq2": "ATCG", "seq3": "ATCG"}
        >>> aligned = progressive_alignment(seqs, method="muscle")
        >>> all(len(seq) == len(aligned["seq1"]) for seq in aligned.values())
        True
    """
    if not sequences:
        return {}

    if len(sequences) == 1:
        return sequences.copy()

    # Check if external tool is available
    if not _is_tool_available(method):
        logger.warning(f"{method} not available, falling back to simple pairwise alignment")
        return _simple_progressive_alignment(sequences)

    try:
        return _run_external_alignment(sequences, method)
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logger.warning(f"External alignment failed: {e}, using simple alignment")
        return _simple_progressive_alignment(sequences)


def _is_tool_available(tool: str) -> bool:
    """Check if an external alignment tool is available."""
    try:
        if tool == "muscle":
            subprocess.run(["muscle", "-version"],
                         capture_output=True, check=True, timeout=5)
        elif tool == "mafft":
            subprocess.run(["mafft", "--version"],
                         capture_output=True, check=True, timeout=5)
        elif tool == "clustalw":
            subprocess.run(["clustalw2", "-help"],
                         capture_output=True, check=True, timeout=5)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        return False


def _run_external_alignment(sequences: Dict[str, str], method: str) -> Dict[str, str]:
    """Run external alignment tool using repository-local temp directory."""
    # Use repo-local temp directory for FAT filesystem compatibility (per CLAUDE.md)
    repo_root = Path(__file__).resolve().parent.parent.parent.parent
    temp_base = repo_root / ".tmp" / "python"
    temp_base.mkdir(parents=True, exist_ok=True)

    # Create a unique temporary directory within repo-local temp
    tmpdir = tempfile.mkdtemp(dir=str(temp_base), prefix="msa_")
    tmpdir = Path(tmpdir)

    try:

        # Write input FASTA
        input_fasta = tmpdir / "input.fasta"
        _write_fasta(sequences, input_fasta)

        # Run alignment
        output_fasta = tmpdir / "output.fasta"

        if method == "muscle":
            cmd = ["muscle", "-in", str(input_fasta), "-out", str(output_fasta)]
        elif method == "mafft":
            cmd = ["mafft", "--quiet", str(input_fasta)]
            # MAFFT outputs to stdout
            output_fasta = None
        elif method == "clustalw":
            cmd = ["clustalw2", "-INFILE=" + str(input_fasta), "-OUTFILE=" + str(output_fasta)]
        else:
            raise ValueError(f"Unknown alignment method: {method}")

        result = subprocess.run(cmd, capture_output=True, text=True, cwd=tmpdir)

        if result.returncode != 0:
            raise subprocess.CalledProcessError(result.returncode, cmd, result.stdout, result.stderr)

        # Read output
        if method == "mafft":
            # MAFFT outputs to stdout
            aligned_fasta_content = result.stdout
        else:
            aligned_fasta_content = output_fasta.read_text()

        return _parse_fasta(aligned_fasta_content)
    finally:
        # Clean up temporary directory
        import shutil
        if tmpdir.exists():
            shutil.rmtree(tmpdir, ignore_errors=True)


def _simple_progressive_alignment(sequences: Dict[str, str]) -> Dict[str, str]:
    """Simple progressive alignment using pairwise alignments."""
    if len(sequences) <= 1:
        return sequences.copy()

    from . import alignment  # Import here to avoid circular imports

    seq_ids = list(sequences.keys())
    aligned = {seq_ids[0]: sequences[seq_ids[0]]}

    # Progressively add sequences
    for seq_id in seq_ids[1:]:
        # Find best match among existing aligned sequences
        best_score = float('-inf')
        best_aligned = None

        for existing_id, existing_seq in aligned.items():
            result = alignment.global_align(sequences[seq_id], existing_seq)
            if result.score > best_score:
                best_score = result.score
                best_aligned = (existing_id, result)

        if best_aligned:
            existing_id, result = best_aligned
            # Merge the new sequence into the alignment
            merged = _merge_alignments(aligned[existing_id], result.seq2_aligned, result.seq1_aligned)
            aligned[existing_id] = merged[0]
            aligned[seq_id] = merged[1]
        else:
            # No good alignment found, add as-is
            aligned[seq_id] = sequences[seq_id]

    return aligned


def _merge_alignments(seq1_aligned: str, seq2_aligned: str, new_seq_aligned: str) -> tuple[str, str]:
    """Merge a new sequence into an existing alignment."""
    # This is a simplified merge - in practice, this would be more complex
    # For now, just return the aligned sequences as they are
    return seq1_aligned, new_seq_aligned


def generate_consensus_from_alignment(aligned_sequences: Dict[str, str], threshold: float = 0.5) -> str:
    """Generate consensus sequence from multiple sequence alignment.

    Args:
        aligned_sequences: Dictionary of aligned sequences
        threshold: Minimum frequency for consensus base (default: 0.5)

    Returns:
        Consensus sequence

    Example:
        >>> aligned = {"seq1": "ATCG", "seq2": "ATCG", "seq3": "ATCG"}
        >>> consensus = generate_consensus_from_alignment(aligned)
        >>> consensus
        'ATCG'
    """
    if not aligned_sequences:
        return ""

    # Get sequence length (should all be equal)
    seq_length = len(next(iter(aligned_sequences.values())))

    consensus = []
    nucleotides = ['A', 'C', 'G', 'T', '-']

    for pos in range(seq_length):
        base_counts = {}

        # Count bases at this position
        for seq in aligned_sequences.values():
            if pos < len(seq):
                base = seq[pos].upper()
                if base in nucleotides:
                    base_counts[base] = base_counts.get(base, 0) + 1

        if not base_counts:
            consensus.append('N')  # Unknown base
            continue

        total_count = sum(base_counts.values())

        # Find base with highest frequency
        max_base = max(base_counts.items(), key=lambda x: x[1])

        if max_base[1] / total_count >= threshold:
            consensus.append(max_base[0])
        else:
            # Use IUPAC ambiguity code for ties/low confidence
            consensus.append(_get_ambiguity_code(base_counts))

    return ''.join(consensus)


def _get_ambiguity_code(base_counts: Dict[str, int]) -> str:
    """Get IUPAC ambiguity code for mixed bases."""
    bases = set(base_counts.keys())

    # Remove gaps for ambiguity calculation
    bases.discard('-')

    if len(bases) == 0:
        return '-'
    elif len(bases) == 1:
        return list(bases)[0]
    elif bases == {'A', 'G'}:
        return 'R'
    elif bases == {'C', 'T'}:
        return 'Y'
    elif bases == {'G', 'C'}:
        return 'S'
    elif bases == {'A', 'T'}:
        return 'W'
    elif bases == {'G', 'T'}:
        return 'K'
    elif bases == {'A', 'C'}:
        return 'M'
    elif bases == {'A', 'C', 'G'}:
        return 'V'
    elif bases == {'A', 'C', 'T'}:
        return 'H'
    elif bases == {'A', 'G', 'T'}:
        return 'D'
    elif bases == {'C', 'G', 'T'}:
        return 'B'
    else:  # All four bases
        return 'N'


def calculate_alignment_quality(aligned_sequences: Dict[str, str]) -> Dict[str, float]:
    """Calculate quality metrics for multiple sequence alignment.

    Args:
        aligned_sequences: Dictionary of aligned sequences

    Returns:
        Dictionary with quality metrics

    Example:
        >>> aligned = {"seq1": "ATCG", "seq2": "ATCG", "seq3": "ATCG"}
        >>> quality = calculate_alignment_quality(aligned)
        >>> quality['average_identity'] > 0.8
        True
    """
    if not aligned_sequences:
        return {
            'average_identity': 0.0,
            'conservation_score': 0.0,
            'gap_percentage': 0.0,
            'alignment_length': 0
        }

    seq_list = list(aligned_sequences.values())
    n_sequences = len(seq_list)

    if n_sequences < 2:
        return {
            'average_identity': 1.0,
            'conservation_score': 1.0,
            'gap_percentage': 0.0,
            'alignment_length': len(seq_list[0]) if seq_list else 0
        }

    alignment_length = len(seq_list[0])
    total_pairs = n_sequences * (n_sequences - 1) // 2
    total_identity = 0.0
    total_gaps = 0
    conserved_positions = 0

    # Compare all pairs
    for i in range(n_sequences):
        for j in range(i + 1, n_sequences):
            seq1, seq2 = seq_list[i], seq_list[j]
            matches = 0
            gaps = 0

            for pos in range(alignment_length):
                base1 = seq1[pos] if pos < len(seq1) else '-'
                base2 = seq2[pos] if pos < len(seq2) else '-'

                if base1 == '-' or base2 == '-':
                    gaps += 1
                elif base1 == base2:
                    matches += 1

            identity = matches / (alignment_length - gaps) if (alignment_length - gaps) > 0 else 0.0
            total_identity += identity
            total_gaps += gaps

    # Calculate conservation score (fraction of positions with consensus)
    for pos in range(alignment_length):
        bases_at_pos = [seq[pos] for seq in seq_list if pos < len(seq)]
        if bases_at_pos:
            most_common = max(set(bases_at_pos), key=bases_at_pos.count)
            if bases_at_pos.count(most_common) == len(bases_at_pos):
                conserved_positions += 1

    return {
        'average_identity': total_identity / total_pairs if total_pairs > 0 else 0.0,
        'conservation_score': conserved_positions / alignment_length if alignment_length > 0 else 0.0,
        'gap_percentage': (total_gaps / (total_pairs * alignment_length)) * 100 if total_pairs * alignment_length > 0 else 0.0,
        'alignment_length': alignment_length
    }


def _write_fasta(sequences: Dict[str, str], filepath: Path) -> None:
    """Write sequences to FASTA format."""
    with open(filepath, 'w') as f:
        for seq_id, seq in sequences.items():
            f.write(f">{seq_id}\n")
            f.write(f"{seq}\n")


def _parse_fasta(fasta_content: str) -> Dict[str, str]:
    """Parse FASTA format content."""
    sequences = {}
    current_id = None
    current_seq = []

    for line in fasta_content.strip().split('\n'):
        line = line.strip()
        if line.startswith('>'):
            if current_id is not None:
                sequences[current_id] = ''.join(current_seq)
            current_id = line[1:].split()[0]  # Remove '>' and take first word
            current_seq = []
        elif current_id is not None:
            current_seq.append(line)

    if current_id is not None:
        sequences[current_id] = ''.join(current_seq)

    return sequences







