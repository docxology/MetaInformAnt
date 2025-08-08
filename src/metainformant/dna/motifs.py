from __future__ import annotations

import re
from typing import Dict, List


_IUPAC: Dict[str, str] = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "U": "U",
    "R": "[AG]",
    "Y": "[CT]",
    "S": "[GC]",
    "W": "[AT]",
    "K": "[GT]",
    "M": "[AC]",
    "B": "[CGT]",
    "D": "[AGT]",
    "H": "[ACT]",
    "V": "[ACG]",
    "N": "[ACGT]",
}


def _iupac_to_regex(pattern: str) -> str:
    regex = "".join(_IUPAC.get(ch.upper(), re.escape(ch)) for ch in pattern)
    return regex


def find_motif_positions(seq: str, motif: str) -> List[int]:
    """Return 0-based start positions where motif (IUPAC) matches.

    Supports standard IUPAC ambiguity codes; case-insensitive.
    """
    if not motif:
        return []
    regex = _iupac_to_regex(motif)
    pat = re.compile(regex, re.IGNORECASE)
    return [m.start() for m in pat.finditer(seq)]


