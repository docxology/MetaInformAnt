"""Protein MCP tools: pairwise alignment identity (read-only)."""

from __future__ import annotations

from typing import Any


def _handle_align(seq1: str, seq2: str, match: int = 1, mismatch: int = -1, gap: int = -2) -> dict:
    from metainformant.protein.sequence.alignment import calculate_alignment_identity, global_align

    result = global_align(seq1, seq2, match=match, mismatch=mismatch, gap=gap)
    identity = calculate_alignment_identity(result)
    return {"score": result.get("score"), "identity": identity, "aligned_length": len(result.get("aligned_seq1", ""))}


TOOL_SPEC: dict[str, Any] = {
    "name": "protein_align_identity",
    "description": "Global (Needleman-Wunsch) protein alignment with identity score.",
    "input_schema": {
        "type": "object",
        "properties": {
            "seq1": {"type": "string"},
            "seq2": {"type": "string"},
            "match": {"type": "integer"},
            "mismatch": {"type": "integer"},
            "gap": {"type": "integer"},
        },
        "required": ["seq1", "seq2"],
    },
    "handler": _handle_align,
    "writes": "read-only",
}


def _handle_composition(sequence: str) -> dict:
    from metainformant.protein.sequence.sequences import (
        amino_acid_composition,
        molecular_weight,
    )

    return {
        "amino_acid_composition": amino_acid_composition(sequence),
        "molecular_weight": molecular_weight(sequence),
        "length": len(sequence),
    }


COMP_SPEC: dict[str, Any] = {
    "name": "protein_composition",
    "description": "Amino-acid composition and molecular weight for a protein sequence.",
    "input_schema": {
        "type": "object",
        "properties": {"sequence": {"type": "string"}},
        "required": ["sequence"],
    },
    "handler": _handle_composition,
    "writes": "read-only",
}

ALL_SPECS: list[dict[str, Any]] = [TOOL_SPEC, COMP_SPEC]
