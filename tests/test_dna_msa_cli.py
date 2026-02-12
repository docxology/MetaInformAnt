from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from metainformant.dna.alignment import msa
from metainformant.dna.sequence import core as sequences


@pytest.mark.external_tool
@pytest.mark.skipif(shutil.which("muscle") is None, reason="MUSCLE not installed")
def test_muscle_cli_alignment_roundtrip(tmp_path: Path) -> None:
    fasta_path = Path(__file__).parent / "data" / "dna" / "toy.fasta"
    if not fasta_path.exists():
        pytest.skip("Test data file not found (tests/data/dna/toy.fasta)")
    seqs = sequences.read_fasta(str(fasta_path))
    aln = msa.align_with_cli(seqs, tool="muscle")
    assert set(aln.keys()) == set(seqs.keys())
    lengths = {len(s) for s in aln.values()}
    assert len(lengths) == 1
