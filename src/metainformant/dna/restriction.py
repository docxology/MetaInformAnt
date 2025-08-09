from __future__ import annotations

import re
from typing import Dict, List

from .motifs import _iupac_to_regex


def find_restriction_sites(seq: str, enzyme_to_motif: Dict[str, str]) -> Dict[str, List[int]]:
    """Find 0-based start positions of enzyme motifs in a sequence.

    Motifs may include IUPAC ambiguity codes. Only the forward strand is
    searched.
    """
    results: Dict[str, List[int]] = {}
    if not seq:
        return {name: [] for name in enzyme_to_motif}
    for name, motif in enzyme_to_motif.items():
        if not motif:
            results[name] = []
            continue
        regex = _iupac_to_regex(motif)
        pat = re.compile(regex, re.IGNORECASE)
        results[name] = [m.start() for m in pat.finditer(seq)]
    return results


