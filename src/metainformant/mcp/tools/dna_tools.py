"""DNA MCP tools: translation, GC/codon composition (pure functions, read-only)."""

from __future__ import annotations

from typing import Any


def _handle_translate(sequence: str, genetic_code: int = 1, to_stop: bool = False) -> dict:
    from metainformant.dna.translation import translate_dna

    protein = translate_dna(sequence, genetic_code=genetic_code, to_stop=to_stop)
    return {"protein": protein, "length_aa": len(protein.replace("*", ""))}


def _handle_composition(sequence: str) -> dict:
    from metainformant.dna.sequence.composition import (
        codon_frequencies,
        gc_content,
        gc_skew,
        melting_temperature,
        nucleotide_frequencies,
    )

    return {
        "gc_content": gc_content(sequence),
        "gc_skew": gc_skew(sequence),
        "melting_temperature": melting_temperature(sequence),
        "nucleotide_frequencies": nucleotide_frequencies(sequence),
        "codon_frequencies": codon_frequencies(sequence),
        "length": len(sequence),
    }


TOOL_SPEC: dict[str, Any] = {
    "name": "dna_translate",
    "description": "Translate a DNA sequence to protein (optionally stopping at first stop codon).",
    "input_schema": {
        "type": "object",
        "properties": {
            "sequence": {"type": "string"},
            "genetic_code": {"type": "integer"},
            "to_stop": {"type": "boolean"},
        },
        "required": ["sequence"],
    },
    "handler": _handle_translate,
    "writes": "read-only",
}

COMP_SPEC: dict[str, Any] = {
    "name": "dna_composition",
    "description": "GC content/skew, melting temperature, nucleotide and codon frequencies for a DNA sequence.",
    "input_schema": {
        "type": "object",
        "properties": {"sequence": {"type": "string"}},
        "required": ["sequence"],
    },
    "handler": _handle_composition,
    "writes": "read-only",
}

ALL_SPECS: list[dict[str, Any]] = [TOOL_SPEC, COMP_SPEC]
