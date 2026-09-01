"""Regression tests for dna.expression fixes (R4_DNA lane, 2026-09-01).

Covers three verified bugs:
1. translation.translate / translate_dna ignored the ``genetic_code``
   parameter (always used the standard table).
2. translation.optimize_codons compared DNA input codons against the
   RNA-keyed genetic code, raising ``ValueError`` for every valid input,
   and treated missing target usages as frequency 1.0 (making the target
   profile meaningless) instead of leaving unstated codons unchanged.
3. codon.back_translate returned RNA codons despite a DNA docstring.
"""

from __future__ import annotations

from metainformant.dna.expression import codon, translation


def test_translate_honors_genetic_code_table_2() -> None:
    """Vertebrate mitochondrial code: AGA is a stop, AUA is Met."""
    # AUG AGA: standard -> M R ; table 2 -> M *
    assert translation.translate("AUGAGA", genetic_code=2) == "M*"
    assert translation.translate("AUGAGA", genetic_code=1) == "MR"


def test_translate_dna_honors_genetic_code_table_2() -> None:
    assert translation.translate_dna("ATGAGA", genetic_code=2) == "M*"
    # ATA (DNA) -> AUA (RNA) -> Met under table 2; AAA -> Lys unchanged
    assert translation.translate_dna("ATAAAA", genetic_code=2) == "MK"
    # Standard table still default: ATA -> Ile
    assert translation.translate_dna("ATAAAA") == "IK"


def test_translate_honors_invertebrate_mito_table_5() -> None:
    # Table 5: AGA/AGG -> Ser, UGA -> Trp
    assert translation.translate("AGAUGG", genetic_code=5) == "SW"


def test_optimize_codons_accepts_dna_input() -> None:
    """DNA codons must be valid input (regression: ValueError on 'ATG')."""
    result = translation.optimize_codons("ATGGCCAAA", {"GCA": 1.0})
    assert result == "ATGGCAAAA"


def test_optimize_codons_uses_target_frequencies() -> None:
    # Two Ala codons present in target with distinct frequencies
    result = translation.optimize_codons("ATGGCTAAA", {"GCA": 0.9, "GCC": 0.1})
    assert result == "ATGGCAAAA"
    result = translation.optimize_codons("ATGGCAAAA", {"GCA": 0.1, "GCC": 0.9})
    assert result == "ATGGCCAAA"


def test_optimize_codons_keeps_codon_without_target_usage() -> None:
    # Empty target: no synonym has stated usage -> sequence unchanged
    assert translation.optimize_codons("ATGGCCAAA", {}) == "ATGGCCAAA"


def test_optimize_codons_preserves_stops() -> None:
    assert translation.optimize_codons("TAACCC", {"TGG": 1.0}) == "TAACCC"


def test_codon_back_translate_returns_dna() -> None:
    """Default preferences use RNA bases internally; output must be DNA."""
    assert codon.back_translate("MA") == "ATGGCT"
    assert "U" not in codon.back_translate("MAVW")


def test_codon_back_translate_accepts_dna_preferences() -> None:
    assert codon.back_translate("MA", {"M": "ATG", "A": "GCA"}) == "ATGGCA"
    # RNA-keyed preferences still normalize to DNA output
    assert codon.back_translate("MA", {"M": "AUG", "A": "GCA"}) == "ATGGCA"


def test_codon_back_translate_unknown_residue() -> None:
    assert codon.back_translate("X") == "NNN"
