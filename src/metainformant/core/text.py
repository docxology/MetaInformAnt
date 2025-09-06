from __future__ import annotations

import re
from pathlib import Path

_WHITESPACE_RE = re.compile(r"\s+")
_SLUG_INVALID_RE = re.compile(r"[^a-z0-9-]+")


def normalize_whitespace(s: str) -> str:
    """Collapse all whitespace to single spaces and strip ends."""
    return _WHITESPACE_RE.sub(" ", s).strip()


def slugify(s: str) -> str:
    """Lowercase, normalize spaces to dashes, and drop invalid URL chars."""
    s = normalize_whitespace(s).lower().replace(" ", "-")
    s = _SLUG_INVALID_RE.sub("", s)
    s = re.sub(r"-+", "-", s)
    return s.strip("-")


def safe_filename(name: str) -> str:
    """Make a safe filename, preserving extension when present."""
    p = Path(name)
    stem = slugify(p.stem)
    suffix = p.suffix
    return f"{stem}{suffix}"


def clean_whitespace(text: str) -> str:
    """Clean excessive whitespace from text.

    Args:
        text: Text to clean

    Returns:
        Text with normalized whitespace
    """
    return normalize_whitespace(text)


def remove_control_chars(text: str) -> str:
    """Remove control characters from text.

    Args:
        text: Text to clean

    Returns:
        Text without control characters
    """
    import unicodedata

    return "".join(char for char in text if unicodedata.category(char)[0] != "C" or char in "\n\t ")


def standardize_gene_name(gene_name: str) -> str:
    """Standardize gene name format.

    Args:
        gene_name: Raw gene name

    Returns:
        Standardized gene name
    """
    # Remove common prefixes/suffixes and standardize case
    name = gene_name.strip().upper()
    # Remove common separators
    name = re.sub(r"[-_\.]", "", name)
    return name


def format_species_name(species_name: str) -> str:
    """Format species name in proper binomial format.

    Args:
        species_name: Raw species name

    Returns:
        Properly formatted species name
    """
    parts = species_name.lower().split()
    if len(parts) >= 2:
        return f"{parts[0].capitalize()} {parts[1].lower()}"
    elif len(parts) == 1:
        return parts[0].capitalize()
    return species_name


def clean_sequence_id(sequence_id: str) -> str:
    """Clean sequence ID to extract main identifier.

    Args:
        sequence_id: Raw sequence ID (e.g., from FASTA header)

    Returns:
        Cleaned sequence identifier
    """
    # Remove leading '>' if present
    clean_id = sequence_id.lstrip(">")

    # Extract main identifier from complex headers
    # Pattern: gi|number|ref|accession|description
    if "|" in clean_id:
        parts = clean_id.split("|")
        # Look for RefSeq pattern
        for i, part in enumerate(parts):
            if part in ["ref", "gb", "emb", "dbj"] and i + 1 < len(parts):
                return parts[i + 1].split()[0]

    # Remove everything after first space or bracket
    clean_id = re.split(r"[\s\[\]]", clean_id)[0]

    return clean_id
