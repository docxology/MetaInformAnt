"""Offline regression tests for DNA-RNA integration utilities."""

from __future__ import annotations

import pytest

from metainformant.dna.integration import rna


def test_predict_transcription_start_sites_single_box_reported_once() -> None:
    """A single TATA box yields exactly one TSS entry at box start + box length."""
    seq = "A" * 30 + "TATA" + "A" * 20
    tss = rna.predict_transcription_start_sites(seq)
    assert tss == [(34, 1.0)]


def test_predict_transcription_start_sites_two_boxes() -> None:
    """Two distinct TATA boxes yield exactly two entries, in position order."""
    seq = "TATA" + "A" * 40 + "TATA" + "A" * 30
    tss = rna.predict_transcription_start_sites(seq)
    assert tss == [(4, 0.5), (48, 1.0)]


def test_predict_transcription_start_sites_case_insensitive() -> None:
    """Matching is case-insensitive."""
    seq = "a" * 30 + "tata" + "a" * 20
    assert rna.predict_transcription_start_sites(seq) == [(34, 1.0)]


def test_correlate_gc_content_with_expression() -> None:
    """Per-gene GC content is correlated against expression (not a constant vector)."""
    dna_features = {
        "g1": {"gc_content": 0.3},
        "g2": {"gc_content": 0.5},
        "g3": {"gc_content": 0.7},
    }
    rna_expression = {"g1": 10.0, "g2": 50.0, "g3": 90.0}
    result = rna.correlate_dna_with_rna_expression(dna_features, rna_expression)
    assert abs(result["gc_expression"]) > 0.9
    assert result["gc_expression"] != 0.0


def test_correlate_gc_content_from_per_gene_sequences() -> None:
    """Per-gene GC content can also be derived from sequences."""
    dna_features = {"gc_free": "AAATTTAAATTT", "gc_rich": "GCGCGCGCGCGC"}
    rna_expression = {"gc_free": 5.0, "gc_rich": 100.0}
    result = rna.correlate_dna_with_rna_expression(dna_features, rna_expression)
    assert result["gc_expression"] == pytest.approx(1.0)


def test_correlate_empty_inputs_return_neutral_shape() -> None:
    """Empty inputs return the empty correlation dict without raising."""
    assert rna.correlate_dna_with_rna_expression({}, {}) == {}
    assert rna.correlate_dna_with_rna_expression({"g1": {"gc_content": 0.5}}, {}) == {}
    assert rna.correlate_dna_with_rna_expression({}, {"g1": 1.0}) == {}


def test_predict_gene_function_empty_sequence_is_neutral() -> None:
    """Empty input returns the same result shape with neutral values."""
    result = rna.predict_gene_function_from_sequence("")
    assert set(result) == {
        "gc_content",
        "orf_count",
        "protein_features",
        "regulatory_elements",
        "predicted_function",
    }
    assert result["gc_content"] == 0.0
    assert result["orf_count"] == 0
    assert result["protein_features"]["cai"] == 0.0
    assert result["protein_features"]["codon_usage"] == {}
    assert result["protein_features"]["protein_length"] == 0
    assert result["regulatory_elements"]["tata_box"] == []
    assert result["predicted_function"] == "unknown"


def test_predict_gene_function_whitespace_sequence_is_neutral() -> None:
    """Whitespace-only input does not raise and yields neutral gc_content."""
    result = rna.predict_gene_function_from_sequence("   \n\t")
    assert result["gc_content"] == 0.0


def test_predict_gene_function_counts_lowercase_gc() -> None:
    """GC counting is case-insensitive (9 of 18 bases are G/C)."""
    result = rna.predict_gene_function_from_sequence("atggccattgtaatgggc")
    assert result["gc_content"] == pytest.approx(0.5)
