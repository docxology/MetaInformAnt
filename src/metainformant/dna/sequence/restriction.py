"""Restriction enzyme analysis and virtual digestion utilities.

This module provides tools for analyzing restriction enzyme recognition sites,
performing virtual DNA digestion, and analyzing restriction fragment patterns.
"""

from __future__ import annotations

from typing import Dict, List, Tuple

from metainformant.core import logging

logger = logging.get_logger(__name__)


# Common restriction enzymes and their recognition sites
RESTRICTION_ENZYMES = {
    # Type II restriction enzymes
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "PstI": "CTGCAG",
    "SalI": "GTCGAC",
    "XbaI": "TCTAGA",
    "SacI": "GAGCTC",
    "KpnI": "GGTACC",
    "SmaI": "CCCGGG",
    "XhoI": "CTCGAG",
    "NotI": "GCGGCCGC",
    "EcoRV": "GATATC",
    "PvuII": "CAGCTG",
    "NdeI": "CATATG",
    "NcoI": "CCATGG",
    "BglII": "AGATCT",
    "ClaI": "ATCGAT",
    "SpeI": "ACTAGT",
    "ApaI": "GGGCCC",
    "Blunt_example": "CCCGGG",  # SmaI produces blunt ends
    # Additional enzymes
    "AluI": "AGCT",
    "HaeIII": "GGCC",
    "RsaI": "GTAC",
    "TaqI": "TCGA",
    "MspI": "CCGG",
    "HpaII": "CCGG",
}


def find_restriction_sites(seq: str, enzymes: List[str]) -> Dict[str, List[int]]:
    """Find restriction sites for specified enzymes in a DNA sequence.

    Args:
        seq: DNA sequence to analyze
        enzymes: List of enzyme names

    Returns:
        Dictionary mapping enzyme names to lists of cut positions

    Example:
        >>> seq = "GAATTCGGATCC"
        >>> sites = find_restriction_sites(seq, ["EcoRI", "BamHI"])
        >>> "EcoRI" in sites
        True
    """
    if not seq:
        return {}

    results = {}

    for enzyme in enzymes:
        if enzyme not in RESTRICTION_ENZYMES:
            logger.warning(f"Unknown enzyme: {enzyme}")
            continue

        recognition_site = RESTRICTION_ENZYMES[enzyme]
        positions = []

        # Find all occurrences of recognition site
        start = 0
        while True:
            pos = seq.upper().find(recognition_site, start)
            if pos == -1:
                break
            positions.append(pos)
            start = pos + 1

        results[enzyme] = positions

    return results


def virtual_digest(seq: str, enzyme: str) -> List[str]:
    """Perform virtual restriction digest of a DNA sequence.

    Args:
        seq: DNA sequence to digest
        enzyme: Restriction enzyme name

    Returns:
        List of DNA fragments after digestion

    Raises:
        ValueError: If enzyme is not recognized

    Example:
        >>> seq = "GAATTCATCG"
        >>> fragments = virtual_digest(seq, "EcoRI")
        >>> len(fragments) >= 1
        True
    """
    if enzyme not in RESTRICTION_ENZYMES:
        raise ValueError(f"Unknown enzyme: {enzyme}")

    if not seq:
        return [seq]

    recognition_site = RESTRICTION_ENZYMES[enzyme]
    site_length = len(recognition_site)

    fragments = []
    last_cut = 0

    # Find all cut positions
    positions = find_restriction_sites(seq, [enzyme])[enzyme]

    # Add sequence end
    positions.append(len(seq))

    # Extract fragments
    for cut_pos in positions:
        fragment = seq[last_cut:cut_pos]
        if fragment:  # Skip empty fragments
            fragments.append(fragment)
        last_cut = cut_pos + site_length  # Skip the recognition site

    return fragments


def calculate_fragment_sizes(seq: str, enzyme: str) -> List[int]:
    """Calculate sizes of restriction fragments.

    Args:
        seq: DNA sequence
        enzyme: Restriction enzyme name

    Returns:
        List of fragment sizes in base pairs

    Example:
        >>> seq = "GAATTCATCG"
        >>> sizes = calculate_fragment_sizes(seq, "EcoRI")
        >>> all(size > 0 for size in sizes)
        True
    """
    fragments = virtual_digest(seq, enzyme)
    return [len(fragment) for fragment in fragments]


def double_digest(seq: str, enzyme1: str, enzyme2: str) -> List[str]:
    """Perform double restriction digest.

    Args:
        seq: DNA sequence
        enzyme1: First restriction enzyme
        enzyme2: Second restriction enzyme

    Returns:
        List of DNA fragments after double digestion

    Example:
        >>> seq = "GAATTCGGATCCATCG"
        >>> fragments = double_digest(seq, "EcoRI", "BamHI")
        >>> len(fragments) >= 1
        True
    """
    if enzyme1 not in RESTRICTION_ENZYMES or enzyme2 not in RESTRICTION_ENZYMES:
        raise ValueError("Unknown enzyme(s)")

    # Find all cut sites from both enzymes
    sites1 = find_restriction_sites(seq, [enzyme1])[enzyme1]
    sites2 = find_restriction_sites(seq, [enzyme2])[enzyme2]

    # Combine and sort all cut positions
    all_sites = []
    enzyme_lengths = {}

    for pos in sites1:
        all_sites.append((pos, enzyme1))
        enzyme_lengths[pos] = len(RESTRICTION_ENZYMES[enzyme1])

    for pos in sites2:
        all_sites.append((pos, enzyme2))
        enzyme_lengths[pos] = len(RESTRICTION_ENZYMES[enzyme2])

    # Sort by position
    all_sites.sort(key=lambda x: x[0])

    # Extract fragments
    fragments = []
    last_cut = 0

    for cut_pos, enzyme in all_sites:
        fragment = seq[last_cut:cut_pos]
        if fragment:
            fragments.append(fragment)
        last_cut = cut_pos + enzyme_lengths[cut_pos]

    # Add remaining fragment
    if last_cut < len(seq):
        fragment = seq[last_cut:]
        if fragment:
            fragments.append(fragment)

    return fragments


def find_unique_sites(seq: str, enzymes: List[str]) -> Dict[str, int]:
    """Find enzymes that cut the sequence at unique positions.

    Args:
        seq: DNA sequence
        enzymes: List of enzyme names to test

    Returns:
        Dictionary mapping enzyme names to number of cut sites

    Example:
        >>> seq = "GAATTCGGATCC"
        >>> sites = find_unique_sites(seq, ["EcoRI", "BamHI"])
        >>> all(count >= 0 for count in sites.values())
        True
    """
    all_sites = find_restriction_sites(seq, enzymes)
    return {enzyme: len(positions) for enzyme, positions in all_sites.items()}


def analyze_restriction_map(seq: str, enzymes: List[str]) -> Dict[str, any]:
    """Generate a complete restriction map analysis.

    Args:
        seq: DNA sequence
        enzymes: List of enzymes to analyze

    Returns:
        Dictionary with comprehensive restriction analysis

    Example:
        >>> seq = "GAATTCGGATCCATCG"
        >>> analysis = analyze_restriction_map(seq, ["EcoRI", "BamHI"])
        >>> "sites" in analysis
        True
    """
    sites = find_restriction_sites(seq, enzymes)

    analysis = {"sequence_length": len(seq), "sites": sites, "fragments": {}, "total_cuts": {}, "fragment_sizes": {}}

    for enzyme in enzymes:
        if enzyme in sites:
            fragments = virtual_digest(seq, enzyme)
            analysis["fragments"][enzyme] = fragments
            analysis["total_cuts"][enzyme] = len(sites[enzyme])
            analysis["fragment_sizes"][enzyme] = [len(f) for f in fragments]

    return analysis


def find_restriction_pattern(seq: str, enzyme: str) -> str:
    """Generate a visual restriction pattern string.

    Args:
        seq: DNA sequence
        enzyme: Restriction enzyme

    Returns:
        String showing cut positions (e.g., "|---|---|---|")

    Example:
        >>> seq = "GAATTCATCG"
        >>> pattern = find_restriction_pattern(seq, "EcoRI")
        >>> "|" in pattern
        True
    """
    if enzyme not in RESTRICTION_ENZYMES:
        return ""

    sites = find_restriction_sites(seq, [enzyme])[enzyme]

    if not sites:
        return "|" + "-" * len(seq) + "|"

    pattern = "|"
    last_pos = 0

    for site in sites:
        # Add dashes for uncut region
        pattern += "-" * (site - last_pos)
        pattern += "|"  # Cut marker
        last_pos = site + len(RESTRICTION_ENZYMES[enzyme])

    # Add remaining sequence
    if last_pos < len(seq):
        pattern += "-" * (len(seq) - last_pos)

    pattern += "|"

    return pattern


def is_palindromic_site(recognition_site: str) -> bool:
    """Check if a restriction site is palindromic.

    Args:
        recognition_site: DNA sequence

    Returns:
        True if the site is palindromic

    Example:
        >>> is_palindromic_site("GAATTC")  # EcoRI site
        True
        >>> is_palindromic_site("CATATG")  # Not palindromic
        False
    """
    if not recognition_site:
        return True

    # Check if sequence equals its reverse complement
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    reverse_complement = "".join(complement.get(c, c) for c in reversed(recognition_site.upper()))

    return recognition_site.upper() == reverse_complement


def get_enzyme_properties() -> Dict[str, Dict[str, any]]:
    """Get properties of all known restriction enzymes.

    Returns:
        Dictionary with enzyme properties
    """
    properties = {}

    for enzyme, site in RESTRICTION_ENZYMES.items():
        properties[enzyme] = {
            "recognition_site": site,
            "site_length": len(site),
            "is_palindromic": is_palindromic_site(site),
            "cuts_blunt": site == site[::-1],  # Simple check for blunt ends
            "cut_position": len(site) // 2,  # Approximate cut position
        }

    return properties


def find_compatible_enzymes(seq: str, min_sites: int = 1, max_sites: int = 5) -> List[str]:
    """Find restriction enzymes that cut the sequence an appropriate number of times.

    Args:
        seq: DNA sequence
        min_sites: Minimum number of cut sites
        max_sites: Maximum number of cut sites

    Returns:
        List of compatible enzyme names

    Example:
        >>> seq = "GAATTCGGATCCATCGATCG"
        >>> enzymes = find_compatible_enzymes(seq, min_sites=1, max_sites=3)
        >>> isinstance(enzymes, list)
        True
    """
    compatible = []

    for enzyme in RESTRICTION_ENZYMES.keys():
        sites = find_restriction_sites(seq, [enzyme])[enzyme]
        if min_sites <= len(sites) <= max_sites:
            compatible.append(enzyme)

    return compatible


def simulate_cloning(
    seq: str,
    insert_enzyme1: str,
    insert_enzyme2: str,
    vector_enzyme1: str,
    vector_enzyme2: str,
    vector_sequence: str | None = None,
) -> Dict[str, any]:
    """Simulate molecular cloning with restriction enzymes.

    Args:
        seq: Insert DNA sequence
        insert_enzyme1: First enzyme for insert digestion
        insert_enzyme2: Second enzyme for insert digestion
        vector_enzyme1: First enzyme for vector digestion
        vector_enzyme2: Second enzyme for vector digestion
        vector_sequence: Optional vector DNA sequence for full simulation

    Returns:
        Dictionary with cloning simulation results including:
        - insert_fragments: Fragments from insert digestion
        - vector_fragments: Fragments from vector digestion (if vector provided)
        - compatible_enzymes: Enzymes used
        - ligation_possible: Whether ligation is theoretically possible
        - enzyme_compatibility: Analysis of enzyme sticky end compatibility

    Example:
        >>> insert = "GAATTCATCGGGATCC"
        >>> results = simulate_cloning(insert, "EcoRI", "BamHI", "EcoRI", "BamHI")
        >>> "insert_fragments" in results
        True
    """
    # Digest insert
    insert_fragments = double_digest(seq, insert_enzyme1, insert_enzyme2)

    # Analyze enzyme compatibility
    insert_enzymes = {insert_enzyme1, insert_enzyme2}
    vector_enzymes = {vector_enzyme1, vector_enzyme2}

    # Check if enzymes create compatible sticky ends
    enzyme_compatibility = {
        "insert_vector_compatible": insert_enzyme1 == vector_enzyme1 and insert_enzyme2 == vector_enzyme2,
        "insert_enzymes": list(insert_enzymes),
        "vector_enzymes": list(vector_enzymes),
    }

    result = {
        "insert_fragments": insert_fragments,
        "compatible_enzymes": [insert_enzyme1, insert_enzyme2],
        "enzyme_compatibility": enzyme_compatibility,
    }

    # Digest vector if sequence provided
    if vector_sequence:
        vector_fragments = double_digest(vector_sequence, vector_enzyme1, vector_enzyme2)
        result["vector_fragments"] = vector_fragments
        result["ligation_possible"] = (
            len(insert_fragments) > 0
            and len(vector_fragments) > 0
            and enzyme_compatibility["insert_vector_compatible"]
        )
    else:
        # Without vector sequence, we can only analyze the insert
        result["vector_fragments"] = None
        result["ligation_possible"] = len(insert_fragments) > 0 and enzyme_compatibility["insert_vector_compatible"]
        result["note"] = "Vector sequence not provided - full simulation requires vector_sequence parameter"

    return result
