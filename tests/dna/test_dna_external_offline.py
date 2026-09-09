"""Offline tests for DNA external/codon utilities (no network, no mocks)."""

from __future__ import annotations

import pytest

from metainformant.dna.expression import codon
from metainformant.dna.external.entrez import _parse_feature_table

GENBANK_SNIPPET = """\
LOCUS       NC_000913            5386 bp    DNA     circular BCT 05-NOV-2020
FEATURES             Location/Qualifiers
     source          1..5386
                     /organism="Escherichia coli K-12"
                     /mol_type="genomic DNA"
     CDS             190..255
                     /gene="thrA"
                     /product="aspartokinase I"
                     /codon_start=1
                     /translation="MKKIAAIV"
ORIGIN
        1 atgaaattaa
//"""


def test_parse_feature_table_extracts_features_and_qualifiers() -> None:
    """Indented feature rows and deeper qualifier rows parse into the expected structure."""
    features = _parse_feature_table(GENBANK_SNIPPET)
    assert [f["type"] for f in features] == ["source", "CDS"]

    assert features[0]["location"] == "1..5386"
    assert features[0]["qualifiers"] == {
        "organism": "Escherichia coli K-12",
        "mol_type": "genomic DNA",
    }

    assert features[1]["location"] == "190..255"
    assert features[1]["qualifiers"] == {
        "gene": "thrA",
        "product": "aspartokinase I",
        "codon_start": "1",
        "translation": "MKKIAAIV",
    }


def test_parse_feature_table_joins_wrapped_qualifiers() -> None:
    """Continuation lines of a wrapped qualifier value are joined onto it."""
    text = """\
FEATURES             Location/Qualifiers
     CDS             10..20
                     /product="some very long product
                     name spanning lines"
"""
    features = _parse_feature_table(text)
    assert features[0]["qualifiers"]["product"] == "some very long product name spanning lines"


def test_parse_feature_table_without_features_section() -> None:
    """A record with no FEATURES section yields no features."""
    assert _parse_feature_table("LOCUS       ABC\n//") == []


def test_cai_rejects_sequence_length_not_multiple_of_three() -> None:
    """cai raises ValueError on bad length like its sibling functions."""
    with pytest.raises(ValueError):
        codon.cai("ATGGCCATTGTAATGGGCC")  # 19 bases


def test_cai_output_pinned() -> None:
    """CAI values are pinned so the O(1) synonymous-frequency refactor stays exact."""
    assert codon.cai("ATGGCCATTGTAATGGGCC" * 3) == pytest.approx(0.9858965833364024)
    assert codon.cai("ATGGCTATGGCTATGGCT") == pytest.approx(0.9517476831808721)


def test_calculate_enc_uses_input_sequence_not_reference() -> None:
    """ENC is computed from the input sequence's codon frequencies; reference_usage is inert."""
    seq = "ATGGCCATGGCTATGGCCATGGCT"  # balanced Ala GCC/GCT use: F = 0.5 -> ENC 2.0
    flat_reference = {c: 1.0 / 64 for c in codon.GENETIC_CODE}
    skewed_reference = {"GCC": 0.95, "GCT": 0.05, "ATG": 1.0}

    from_input = codon.calculate_enc(seq)
    assert from_input == pytest.approx(2.0)
    assert codon.calculate_enc(seq, reference_usage=flat_reference) == pytest.approx(from_input)
    assert codon.calculate_enc(seq, reference_usage=skewed_reference) == pytest.approx(from_input)

    # A reference-only computation (the old behaviour) would read the flat table
    # as unbiased (~61) for this sequence; the input-driven value differs.
    assert from_input != pytest.approx(61.0)


def test_codon_usage_matches_manual_counts() -> None:
    """codon_usage frequencies equal manually counted codons."""
    assert codon.codon_usage("ATGATGGCC") == {"ATG": 2 / 3, "GCC": 1 / 3}
    assert codon.codon_usage("ATGGCCATT") == {"ATG": 1 / 3, "GCC": 1 / 3, "ATT": 1 / 3}
    assert codon.codon_usage("") == {}
