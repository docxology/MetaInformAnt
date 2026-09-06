"""Multiple sequence alignment for DNA sequences.

This module provides tools for aligning multiple DNA sequences using various
algorithms including progressive alignment and integration with external tools
like MUSCLE, MAFFT, and ClustalW.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def align_msa(sequences: Dict[str, str], method: str = "auto") -> Dict[str, str]:
    """Align multiple DNA sequences with an external tool or local fallback.

    Args:
        sequences: Dictionary mapping sequence IDs to DNA sequences.
        method: ``"auto"`` chooses the first available supported CLI tool,
            otherwise one of ``"muscle"``, ``"mafft"``, ``"clustalo"``, or
            ``"clustalw"``.

    Returns:
        Dictionary mapping sequence IDs to equal-length aligned sequences.
    """
    if method == "auto":
        for candidate in ("muscle", "mafft", "clustalo", "clustalw"):
            if _is_tool_available(candidate):
                return _pad_alignment_to_common_length(_run_external_alignment(sequences, candidate))
        return _pad_alignment_to_common_length(_simple_progressive_alignment(sequences))

    return _pad_alignment_to_common_length(progressive_alignment(sequences, method=method))


def align_with_cli(sequences: Dict[str, str], tool: str = "muscle") -> Dict[str, str]:
    """Align sequences with a named external alignment command.

    Unlike :func:`align_msa`, this explicit CLI helper does not silently fall
    back when the requested command is unavailable.
    """
    if not _is_tool_available(tool):
        raise FileNotFoundError(f"{tool} alignment tool is not available on PATH")
    return _pad_alignment_to_common_length(_run_external_alignment(sequences, tool))


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
    supported = ("muscle", "mafft", "clustalo", "clustalw")
    if method not in supported:
        raise ValueError(f"Unknown alignment method: {method}")

    if not sequences:
        return {}

    if len(sequences) == 1:
        return sequences.copy()

    if not _is_tool_available(method):
        logger.warning(f"{method} not available, falling back to simple pairwise alignment")
        return _pad_alignment_to_common_length(_simple_progressive_alignment(sequences))

    try:
        return _pad_alignment_to_common_length(_run_external_alignment(sequences, method))
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logger.warning(f"External alignment failed: {e}, using simple alignment")
        return _pad_alignment_to_common_length(_simple_progressive_alignment(sequences))


def _is_tool_available(tool: str) -> bool:
    """Check if an external alignment tool is available."""
    try:
        if tool == "muscle":
            subprocess.run(["muscle", "-version"], capture_output=True, check=True, timeout=5)
        elif tool == "mafft":
            subprocess.run(["mafft", "--version"], capture_output=True, check=True, timeout=5)
        elif tool == "clustalo":
            subprocess.run(["clustalo", "--version"], capture_output=True, check=True, timeout=5)
        elif tool == "clustalw":
            subprocess.run(["clustalw2", "-help"], capture_output=True, check=True, timeout=5)
        else:
            return False
        return True
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        return False


def _run_external_alignment(sequences: Dict[str, str], method: str) -> Dict[str, str]:
    """Run external alignment tool using repository-local temp directory."""
    # Use repo-local temp directory for FAT filesystem compatibility (per CLAUDE.md).
    # parents: [0] alignment/, [1] dna/, [2] metainformant/, [3] src/, [4] repo root.
    repo_root = Path(__file__).resolve().parents[4]
    temp_base = repo_root / ".tmp" / "python"
    temp_base.mkdir(parents=True, exist_ok=True)

    # Create a unique temporary directory within repo-local temp
    tmpdir = Path(tempfile.mkdtemp(dir=str(temp_base), prefix="msa_"))

    try:
        # Write input FASTA
        input_fasta = tmpdir / "input.fasta"
        _write_fasta(sequences, input_fasta)

        # Run alignment
        output_fasta = tmpdir / "output.fasta"

        if method == "mafft":
            cmd = ["mafft", "--quiet", str(input_fasta)]
            # MAFFT outputs to stdout
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=tmpdir)
        elif method == "muscle":
            result = _run_muscle(input_fasta, output_fasta, tmpdir)
        elif method == "clustalo":
            cmd = ["clustalo", "-i", str(input_fasta), "-o", str(output_fasta), "--force"]
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=tmpdir)
        elif method == "clustalw":
            cmd = ["clustalw2", "-INFILE=" + str(input_fasta), "-OUTFILE=" + str(output_fasta)]
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=tmpdir)
        else:
            raise ValueError(f"Unknown alignment method: {method}")

        if result.returncode != 0:
            raise subprocess.CalledProcessError(
                result.returncode,
                result.args,
                result.stdout,
                result.stderr,
            )

        # Read output
        if method == "mafft":
            # MAFFT outputs to stdout
            aligned_fasta_content = result.stdout
        else:
            aligned_fasta_content = output_fasta.read_text()

        return _pad_alignment_to_common_length(_parse_fasta(aligned_fasta_content))
    finally:
        # Clean up temporary directory
        if tmpdir.exists():
            shutil.rmtree(tmpdir, ignore_errors=True)


def _run_muscle(input_fasta: Path, output_fasta: Path, tmpdir: Path) -> subprocess.CompletedProcess[str]:
    """Run MUSCLE with v5 syntax first, then v3 syntax for older installs."""
    version_5_cmd = ["muscle", "-align", str(input_fasta), "-output", str(output_fasta)]
    result = subprocess.run(version_5_cmd, capture_output=True, text=True, cwd=tmpdir)
    if result.returncode == 0:
        return result

    version_3_cmd = ["muscle", "-in", str(input_fasta), "-out", str(output_fasta)]
    return subprocess.run(version_3_cmd, capture_output=True, text=True, cwd=tmpdir)


def _simple_progressive_alignment(sequences: Dict[str, str]) -> Dict[str, str]:
    """Align every sequence to the first sequence, propagating gap columns.

    Each sequence is pairwise-aligned against the (ungapped) first sequence;
    gap columns opened in the reference are inserted into every previously
    aligned sequence so all columns stay homologous.
    """
    if len(sequences) <= 1:
        return sequences.copy()

    from ..alignment import pairwise as alignment  # Import here to avoid circular imports

    seq_ids = list(sequences.keys())
    ref_id = seq_ids[0]
    ref_ungapped = sequences[ref_id]

    # One column per alignment position; each column holds one char per
    # sequence added so far. ``col_is_ref[k]`` marks whether column k holds
    # a character of the (ungapped) reference sequence or an inserted gap.
    columns: List[List[str]] = [[char] for char in ref_ungapped]
    col_is_ref: List[bool] = [True] * len(ref_ungapped)
    aligned_ids = [ref_id]

    for seq_id in seq_ids[1:]:
        result = alignment.global_align(ref_ungapped, sequences[seq_id])
        ref_aln, new_aln = result.seq1_aligned, result.seq2_aligned

        new_columns: List[List[str]] = []
        new_col_is_ref: List[bool] = []
        k = 0
        for ref_char, new_char in zip(ref_aln, new_aln):
            if ref_char == "-":
                # Gap opened in the reference: every previously aligned
                # sequence gets a gap; the new sequence keeps its base.
                new_columns.append(["-"] * len(aligned_ids) + [new_char])
                new_col_is_ref.append(False)
            else:
                # Previously inserted gap columns anchored before this
                # reference character are kept; the new sequence has no
                # counterpart for them, so it receives a gap.
                while not col_is_ref[k]:
                    gap_column = columns[k]
                    gap_column.append("-")
                    new_columns.append(gap_column)
                    new_col_is_ref.append(False)
                    k += 1
                column = columns[k]
                column.append(new_char)
                new_columns.append(column)
                new_col_is_ref.append(True)
                k += 1
        while k < len(columns):
            gap_column = columns[k]
            gap_column.append("-")
            new_columns.append(gap_column)
            new_col_is_ref.append(False)
            k += 1
        columns = new_columns
        col_is_ref = new_col_is_ref
        aligned_ids.append(seq_id)

    return {seq_id: "".join(columns[pos][i] for pos in range(len(columns))) for i, seq_id in enumerate(aligned_ids)}


def _pad_alignment_to_common_length(sequences: Dict[str, str]) -> Dict[str, str]:
    """Right-pad aligned sequences with gaps so every sequence has equal length."""
    if not sequences:
        return {}

    max_length = max(len(seq) for seq in sequences.values())
    return {seq_id: seq.ljust(max_length, "-") for seq_id, seq in sequences.items()}


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
    nucleotides = ["A", "C", "G", "T", "-"]

    for pos in range(seq_length):
        base_counts: Dict[str, int] = {}

        # Count bases at this position
        for seq in aligned_sequences.values():
            if pos < len(seq):
                base = seq[pos].upper()
                if base in nucleotides:
                    base_counts[base] = base_counts.get(base, 0) + 1

        if not base_counts:
            consensus.append("N")  # Unknown base
            continue

        total_count = sum(base_counts.values())

        # Find base with highest frequency
        max_base = max(base_counts.items(), key=lambda x: x[1])

        if max_base[1] / total_count >= threshold:
            consensus.append(max_base[0])
        else:
            # Use IUPAC ambiguity code for ties/low confidence
            consensus.append(_get_ambiguity_code(base_counts))

    return "".join(consensus)


def _get_ambiguity_code(base_counts: Dict[str, int]) -> str:
    """Get IUPAC ambiguity code for mixed bases."""
    bases = set(base_counts.keys())

    # Remove gaps for ambiguity calculation
    bases.discard("-")

    if len(bases) == 0:
        return "-"
    elif len(bases) == 1:
        return list(bases)[0]
    elif bases == {"A", "G"}:
        return "R"
    elif bases == {"C", "T"}:
        return "Y"
    elif bases == {"G", "C"}:
        return "S"
    elif bases == {"A", "T"}:
        return "W"
    elif bases == {"G", "T"}:
        return "K"
    elif bases == {"A", "C"}:
        return "M"
    elif bases == {"A", "C", "G"}:
        return "V"
    elif bases == {"A", "C", "T"}:
        return "H"
    elif bases == {"A", "G", "T"}:
        return "D"
    elif bases == {"C", "G", "T"}:
        return "B"
    else:  # All four bases
        return "N"


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
        return {"average_identity": 0.0, "conservation_score": 0.0, "gap_percentage": 0.0, "alignment_length": 0}

    seq_list = list(aligned_sequences.values())
    n_sequences = len(seq_list)

    if n_sequences < 2:
        return {
            "average_identity": 1.0,
            "conservation_score": 1.0,
            "gap_percentage": 0.0,
            "alignment_length": len(seq_list[0]) if seq_list else 0,
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
                base1 = seq1[pos] if pos < len(seq1) else "-"
                base2 = seq2[pos] if pos < len(seq2) else "-"

                if base1 == "-" or base2 == "-":
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
        "average_identity": total_identity / total_pairs if total_pairs > 0 else 0.0,
        "conservation_score": conserved_positions / alignment_length if alignment_length > 0 else 0.0,
        "gap_percentage": (
            (total_gaps / (total_pairs * alignment_length)) * 100 if total_pairs * alignment_length > 0 else 0.0
        ),
        "alignment_length": alignment_length,
    }


def _write_fasta(sequences: Dict[str, str], filepath: Path) -> None:
    """Write sequences to FASTA format."""
    with open(filepath, "w") as f:
        for seq_id, seq in sequences.items():
            f.write(f">{seq_id}\n")
            f.write(f"{seq}\n")


def _parse_fasta(fasta_content: str) -> Dict[str, str]:
    """Parse FASTA format content."""
    sequences = {}
    current_id = None
    current_seq: List[str] = []

    for line in fasta_content.strip().split("\n"):
        line = line.strip()
        if line.startswith(">"):
            if current_id is not None:
                sequences[current_id] = "".join(current_seq)
            current_id = line[1:].split()[0]  # Remove '>' and take first word
            current_seq = []
        elif current_id is not None:
            current_seq.append(line)

    if current_id is not None:
        sequences[current_id] = "".join(current_seq)

    return sequences
