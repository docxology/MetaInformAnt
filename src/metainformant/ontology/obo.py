"""OBO format parsing and processing.

This module provides functionality to parse OBO (Open Biological and Biomedical Ontologies)
format files and convert them into METAINFORMANT ontology structures.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Any, Iterator, Optional
from metainformant.core import logging, io, validation
from .types import Ontology, Term, Relationship, create_term, create_relationship, create_ontology

logger = logging.get_logger(__name__)


def parse_obo(path: str | Path) -> Ontology:
    """Parse an OBO format file into an Ontology object.

    Args:
        path: Path to the OBO file to parse.

    Returns:
        Ontology object containing all parsed terms and relationships.

    Raises:
        FileNotFoundError: If the OBO file does not exist.
        ValueError: If the file format is invalid.

    Examples:
        >>> ontology = parse_obo("data/go.obo")
        >>> len(ontology)
        50000
    """
    path = validation.validate_path_exists(Path(path))

    logger.info(f"Parsing OBO file: {path}")

    terms = {}
    relationships = []
    header_metadata = {}

    # Read file content
    content = io.read_text(path)
    lines = content.splitlines()

    # Parse header
    header_metadata = _parse_header(lines)

    # Parse stanzas
    stanza_lines = []
    i = 0

    while i < len(lines):
        line = lines[i].strip()

        if line.startswith("[") and line.endswith("]"):
            # Process previous stanza if exists
            if stanza_lines:
                stanza_type = _get_stanza_type(stanza_lines)
                if stanza_type == "Term":
                    term = _parse_term_stanza(stanza_lines)
                    if term:
                        terms[term.id] = term
                elif stanza_type == "Typedef":
                    # Handle typedefs (relationship types)
                    typedef_info = _parse_typedef_stanza(stanza_lines)
                    if typedef_info:
                        header_metadata[f"typedef:{typedef_info['id']}"] = typedef_info
                elif stanza_type == "Instance":
                    # Handle instances if present
                    instance_info = _parse_instance_stanza(stanza_lines)
                    if instance_info:
                        header_metadata[f"instance:{instance_info['id']}"] = instance_info

            # Start new stanza
            stanza_lines = [line]
        elif line and not line.startswith("!"):
            # Add non-comment lines to current stanza
            stanza_lines.append(line)

        i += 1

    # Process final stanza
    if stanza_lines:
        stanza_type = _get_stanza_type(stanza_lines)
        if stanza_type == "Term":
            term = _parse_term_stanza(stanza_lines)
            if term:
                terms[term.id] = term

    # Build relationships from term relationships
    relationships = _extract_relationships_from_terms(terms)

    logger.info(f"Parsed {len(terms)} terms and {len(relationships)} relationships")

    return create_ontology(terms, relationships, **header_metadata)


def _parse_header(lines: List[str]) -> Dict[str, Any]:
    """Parse the OBO file header."""
    header = {}
    i = 0

    while i < len(lines) and not lines[i].startswith("["):
        line = lines[i].strip()
        if line and not line.startswith("!"):
            if ":" in line:
                key, value = line.split(":", 1)
                key = key.strip()
                value = value.strip()
                header[key] = value
        i += 1

    return header


def _get_stanza_type(stanza_lines: List[str]) -> Optional[str]:
    """Extract stanza type from stanza lines."""
    if stanza_lines and stanza_lines[0].startswith("[") and stanza_lines[0].endswith("]"):
        return stanza_lines[0][1:-1]
    return None


def _parse_term_stanza(stanza_lines: List[str]) -> Optional[Term]:
    """Parse a [Term] stanza into a Term object."""
    term_data = {}

    for line in stanza_lines[1:]:  # Skip [Term] line
        line = line.strip()
        if not line or line.startswith("!"):
            continue

        if ":" not in line:
            continue

        key, value = line.split(":", 1)
        key = key.strip()
        value = value.strip()

        # Handle multi-line values and qualifiers
        if key in ["name", "namespace", "def", "comment"]:
            # Remove quotes if present
            value = value.strip('"')
            term_data[key] = value
        elif key == "id":
            term_data["id"] = value
        elif key == "is_obsolete":
            term_data["is_obsolete"] = value.lower() == "true"
        elif key.startswith("synonym"):
            if "synonyms" not in term_data:
                term_data["synonyms"] = []
            # Extract synonym text from quotes
            synonym_match = re.search(r'"([^"]*)"', value)
            if synonym_match:
                term_data["synonyms"].append(synonym_match.group(1))
        elif key.startswith("xref"):
            if "xrefs" not in term_data:
                term_data["xrefs"] = []
            term_data["xrefs"].append(value)
        elif key in [
            "is_a",
            "part_of",
            "regulates",
            "negatively_regulates",
            "positively_regulates",
            "has_part",
            "occurs_in",
            "happens_during",
            "ends_during",
        ]:
            if "relationships" not in term_data:
                term_data["relationships"] = []
            term_data["relationships"].append((key, value))

    if "id" not in term_data:
        return None

    # Create term
    term = create_term(
        id=term_data["id"],
        name=term_data.get("name"),
        definition=term_data.get("def"),
        namespace=term_data.get("namespace"),
        synonyms=term_data.get("synonyms", []),
        xrefs=term_data.get("xrefs", []),
        is_obsolete=term_data.get("is_obsolete", False),
    )

    # Store relationships for later processing
    term.metadata["relationships"] = term_data.get("relationships", [])

    return term


def _parse_typedef_stanza(stanza_lines: List[str]) -> Optional[Dict[str, Any]]:
    """Parse a [Typedef] stanza."""
    typedef_data = {}

    for line in stanza_lines[1:]:
        line = line.strip()
        if not line or line.startswith("!"):
            continue

        if ":" not in line:
            continue

        key, value = line.split(":", 1)
        key = key.strip()
        value = value.strip()

        if key == "id":
            typedef_data["id"] = value
        elif key == "name":
            typedef_data["name"] = value
        elif key == "def":
            typedef_data["definition"] = value.strip('"')

    return typedef_data if "id" in typedef_data else None


def _parse_instance_stanza(stanza_lines: List[str]) -> Optional[Dict[str, Any]]:
    """Parse an [Instance] stanza."""
    instance_data = {}

    for line in stanza_lines[1:]:
        line = line.strip()
        if not line or line.startswith("!"):
            continue

        if ":" not in line:
            continue

        key, value = line.split(":", 1)
        key = key.strip()
        value = value.strip()

        if key == "id":
            instance_data["id"] = value
        elif key == "name":
            instance_data["name"] = value

    return instance_data if "id" in instance_data else None


def _extract_relationships_from_terms(terms: Dict[str, Term]) -> List[Relationship]:
    """Extract relationships from term metadata."""
    relationships = []

    for term in terms.terms.values():
        if "relationships" in term.metadata:
            for rel_type, target_id in term.metadata["relationships"]:
                relationship = create_relationship(source=term.id, target=target_id, relation_type=rel_type)
                relationships.append(relationship)

    return relationships


def validate_obo_format(path: str | Path) -> tuple[bool, List[str]]:
    """Validate OBO file format.

    Args:
        path: Path to the OBO file to validate.

    Returns:
        Tuple of (is_valid, error_messages).
    """
    path = validation.validate_path_exists(Path(path))

    errors = []

    try:
        content = io.read_text(path)
        lines = content.splitlines()

        # Check for header
        has_header = False
        for line in lines[:10]:  # Check first 10 lines for format-version
            if line.startswith("format-version:"):
                has_header = True
                break

        if not has_header:
            errors.append("Missing format-version in header")

        # Check for at least one [Term] stanza
        has_term = False
        for line in lines:
            if line.strip() == "[Term]":
                has_term = True
                break

        if not has_term:
            errors.append("No [Term] stanzas found")

        # Try parsing to catch other issues
        try:
            parse_obo(path)
        except Exception as e:
            errors.append(f"Parsing failed: {str(e)}")

    except Exception as e:
        errors.append(f"File reading failed: {str(e)}")

    return len(errors) == 0, errors


def get_obo_statistics(path: str | Path) -> Dict[str, Any]:
    """Get statistics about an OBO file.

    Args:
        path: Path to the OBO file.

    Returns:
        Dictionary with file statistics.
    """
    path = validation.validate_path_exists(Path(path))

    content = io.read_text(path)
    lines = content.splitlines()

    stats = {
        "total_lines": len(lines),
        "term_count": 0,
        "typedef_count": 0,
        "instance_count": 0,
        "relationship_count": 0,
        "obsolete_count": 0,
    }

    current_stanza = None

    for line in lines:
        line = line.strip()

        if line.startswith("[") and line.endswith("]"):
            current_stanza = line[1:-1]
            if current_stanza == "Term":
                stats["term_count"] += 1
            elif current_stanza == "Typedef":
                stats["typedef_count"] += 1
            elif current_stanza == "Instance":
                stats["instance_count"] += 1
        elif current_stanza == "Term" and line.startswith("is_obsolete: true"):
            stats["obsolete_count"] += 1
        elif current_stanza == "Term" and ":" in line:
            key = line.split(":", 1)[0].strip()
            if key in ["is_a", "part_of", "regulates", "has_part", "occurs_in"]:
                stats["relationship_count"] += 1

    return stats
