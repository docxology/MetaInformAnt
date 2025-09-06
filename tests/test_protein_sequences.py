from __future__ import annotations

from pathlib import Path


def test_parse_fasta_and_composition(tmp_path: Path):
    from metainformant.protein.sequences import (
        calculate_aa_composition,
        is_valid_protein_sequence,
        kmer_frequencies,
        parse_fasta,
    )

    fasta = ">seq1\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQANN\n>seq2\nGASGDLGKK\n"
    fasta_path = tmp_path / "toy.faa"
    fasta_path.write_text(fasta)

    records = parse_fasta(fasta_path)
    assert set(records.keys()) == {"seq1", "seq2"}
    assert records["seq1"].startswith("MKTAYI")

    comp = calculate_aa_composition(records["seq2"])  # length 9
    assert abs(sum(comp.values()) - 1.0) < 1e-6
    # Ensure expected residues are present
    for aa in {"G", "A", "S", "D", "L", "K"}:
        assert aa in comp

    assert is_valid_protein_sequence("MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTI") is True
    assert is_valid_protein_sequence("MTEY*KL") is False

    km = kmer_frequencies(records["seq1"], k=2)
    assert sum(km.values()) == len(records["seq1"]) - 1
    assert all(len(k) == 2 for k in km)
