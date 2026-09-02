"""Math/popgen MCP tool: Hardy-Weinberg and Fst summaries (descriptive)."""

from __future__ import annotations

from typing import Any


def _handle_popgen_summary(allele_freq_p: float) -> dict:
    """Descriptive population-genetics summary from an allele frequency p."""
    from metainformant.math.popgen import fst_from_freqs, hardy_weinberg_allele_freqs

    aa, het, qq = hardy_weinberg_allele_freqs(allele_freq_p)
    return {
        "freq_AA": aa,
        "freq_Aa": het,
        "freq_aa": qq,
        "expected_heterozygosity_hw": het,
        "fst_example_two_loci": fst_from_freqs([allele_freq_p], [1 - allele_freq_p]),
        "note": "descriptive summary only; no inferential tests",
    }


TOOL_SPEC: dict[str, Any] = {
    "name": "math_popgen_summary",
    "description": "Hardy-Weinberg decomposition and Fst-style descriptive summary from an allele frequency.",
    "input_schema": {
        "type": "object",
        "properties": {"allele_freq_p": {"type": "number"}},
        "required": ["allele_freq_p"],
    },
    "handler": _handle_popgen_summary,
    "writes": "read-only",
}

ALL_SPECS: list[dict[str, Any]] = [TOOL_SPEC]
