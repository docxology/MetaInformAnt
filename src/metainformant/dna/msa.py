from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Dict

from Bio.Align import MultipleSeqAlignment, PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def align_msa(id_to_seq: Dict[str, str], method: str = "auto") -> Dict[str, str]:
    """Very small MSA utility using progressive pairwise alignment.

    This avoids external binary dependencies (MUSCLE/Clustal) while keeping tests lightweight.
    Not intended for large datasets.
    """
    # Start with the first sequence as reference
    items = list(id_to_seq.items())
    if not items:
        return {}
    ref_id, ref_seq = items[0]
    aligned: Dict[str, str] = {ref_id: ref_seq}

    aligner = PairwiseAligner()
    aligner.mode = "global"

    def pad_to_length(s: str, n: int) -> str:
        return s + ("-" * (n - len(s)))

    for seq_id, seq in items[1:]:
        # Align to the current consensus (use reference-guided progressive strategy)
        alignment = aligner.align(ref_seq, seq)[0]
        a_ref = str(alignment.aligned[0])
        # PairwiseAligner doesn't directly give gapped strings; rebuild from coordinates
        # Fallback: compute naive edit script by aligning using scores and traceback via format
        # Use alignment.format() to get a three-line gapped representation
        lines = alignment.format().splitlines()
        gapped_ref = lines[0].strip()
        gapped_seq = lines[2].strip()

        # Update previous aligned sequences to include any gaps introduced in gapped_ref
        new_aligned: Dict[str, str] = {}
        ref_pointer = 0
        for i, ch in enumerate(gapped_ref):
            if ch != "-":
                ref_pointer += 1
        # For simplicity, pad all existing sequences with '-' to match new length
        new_len = len(gapped_ref)
        for k, v in aligned.items():
            new_aligned[k] = pad_to_length(v, new_len)
        new_aligned[seq_id] = gapped_seq
        aligned = new_aligned
        ref_seq = gapped_ref

    # Ensure equal lengths
    max_len = max(len(s) for s in aligned.values())
    for k in list(aligned.keys()):
        aligned[k] = pad_to_length(aligned[k], max_len)
    return aligned


def align_with_cli(id_to_seq: Dict[str, str], tool: str = "muscle") -> Dict[str, str]:
    """Align with external CLI (MUSCLE/Clustal Omega) if present.

    Requires the tool to be available in PATH. Writes temporary FASTA and reads back.
    """
    exe = shutil.which(tool)
    if exe is None:
        raise RuntimeError(f"{tool} not found in PATH")
    tmpdir = Path.cwd() / ".tmp_msa"
    tmpdir.mkdir(exist_ok=True)
    in_fa = tmpdir / "in.fasta"
    out_fa = tmpdir / "out.fasta"
    # write
    with in_fa.open("w") as fh:
        for k, v in id_to_seq.items():
            fh.write(f">{k}\n{v}\n")
    # call
    if tool.lower() == "muscle":
        cmd = [exe, "-align", str(in_fa), "-output", str(out_fa)]
    elif tool.lower() in {"clustalo", "clustalomega"}:
        cmd = [exe, "-i", str(in_fa), "-o", str(out_fa), "--force"]
    else:
        raise RuntimeError(f"Unsupported tool: {tool}")
    subprocess.run(cmd, check=True)
    # read
    from Bio import SeqIO

    out: Dict[str, str] = {}
    for rec in SeqIO.parse(str(out_fa), "fasta"):
        out[rec.id] = str(rec.seq)
    return out
