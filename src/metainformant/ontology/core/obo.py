"""OBO format parsing and processing.

This module provides functionality to parse OBO (Open Biological and Biomedical Ontologies)
format files and convert them into METAINFORMANT ontology structures.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

from metainformant.core.data import validation
from metainformant.core.utils import logging

from .types import Ontology, Relationship, Term, create_ontology, create_relationship, create_term

logger = logging.get_logger(__name__)

# OBO flat file format 1.4, "Escape characters": a backslash followed by one of these stands for the character itself,
# and \n, \t and \W stand for a newline, a tab and a space. Any other backslash sequence is left exactly as written.
_OBO_ESCAPED_PUNCTUATION = frozenset(':,"\\()[]{}')
_OBO_ESCAPED_WHITESPACE = {"n": "\n", "t": "\t", "W": " "}


def _unescape_obo(text: str) -> str:
    """Decode OBO backslash escapes (``\\"``, ``\\\\``, ``\\,``, ``\\n``, ``\\t``, ``\\W`` ...)."""
    if "\\" not in text:
        return text
    decoded: List[str] = []
    i = 0
    while i < len(text):
        char = text[i]
        if char == "\\" and i + 1 < len(text):
            following = text[i + 1]
            if following in _OBO_ESCAPED_WHITESPACE:
                decoded.append(_OBO_ESCAPED_WHITESPACE[following])
                i += 2
                continue
            if following in _OBO_ESCAPED_PUNCTUATION:
                decoded.append(following)
                i += 2
                continue
        decoded.append(char)
        i += 1
    return "".join(decoded)


def quote_obo_string(text: str) -> str:
    """Render ``text`` as an OBO quoted string, surrounding quotes included.

    Backslashes and double quotes are escaped, and line breaks and tabs are written as ``\\n`` and ``\\t`` so the
    value stays on one line; :func:`parse_obo_definition` reads the result back to the original text.
    """
    escaped = text.replace("\r\n", "\n").replace("\r", "\n")
    escaped = escaped.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n").replace("\t", "\\t")
    return f'"{escaped}"'


def format_obo_xref_list(xrefs: Iterable[str]) -> str:
    """Render dbxrefs as the bracketed list that follows a definition, e.g. ``[GOC:mah, PMID:123]``.

    Characters that would end an item or the list (backslash, comma, closing bracket, double quote) are escaped, so
    :func:`parse_obo_definition` reads back exactly the items given here.
    """
    items = [
        xref.replace("\\", "\\\\").replace(",", "\\,").replace("]", "\\]").replace('"', '\\"').replace("\n", "\\n")
        for xref in xrefs
    ]
    return "[" + ", ".join(items) + "]"


def _closing_quote_index(value: str) -> int:
    """Index of the quote that closes the quoted string opening at ``value[0]``, or -1 if there is none.

    An unescaped quote only closes the string when what follows it is the end of the value, the ``[`` of the dbxref
    list, the ``{`` of trailing modifiers or a ``!`` comment. Any other unescaped quote is kept as literal text, so
    definitions written without escaping (as older versions of this package did) still read back whole.
    """
    start = 1
    while True:
        index = value.find('"', start)
        if index < 0:
            return -1
        # A quote is escaped when an odd number of backslashes immediately precede it
        backslashes = 0
        while index - 1 - backslashes > 0 and value[index - 1 - backslashes] == "\\":
            backslashes += 1
        if backslashes % 2 == 0:
            tail = value[index + 1 :].lstrip()
            if not tail or tail[0] in "[{!":
                return index
        start = index + 1


def _parse_obo_xref_list(rest: str) -> List[str]:
    """Parse the dbxref list at the start of ``rest`` (which begins with ``[``) into its items.

    Items are separated by commas outside quotes and brackets; an item may carry a quoted description.
    """
    end = rest.find("]")
    plain = rest[1:] if end < 0 else rest[1:end]
    if '"' not in plain and "\\" not in plain:  # common case: no quotes or escapes, so the first "]" ends the list
        return [item.strip() for item in plain.split(",") if item.strip()]

    items: List[str] = []
    current: List[str] = []
    in_quote = False
    i = 1
    while i < len(rest):
        char = rest[i]
        if char == "\\" and i + 1 < len(rest):
            current.append(rest[i : i + 2])  # keep an escape pair together; decoded below
            i += 2
            continue
        if char == '"':
            in_quote = not in_quote
        elif not in_quote and char == "]":
            break
        elif not in_quote and char == ",":
            items.append("".join(current))
            current = []
            i += 1
            continue
        current.append(char)
        i += 1
    items.append("".join(current))
    return [_unescape_obo(item.strip()) for item in items if item.strip()]


def parse_obo_definition(value: str) -> Tuple[str, List[str]]:
    """Parse the value of an OBO ``def:`` tag into ``(definition_text, dbxrefs)``.

    The value is a quoted string, then a bracketed dbxref list, then optional ``{...}`` trailing modifiers and a
    ``! comment``, for example ``"A process that \\"x\\" does." [GOC:mah, PMID:123]``. The text is returned without
    its quotes, with escape sequences decoded; the dbxrefs are returned as a list, so neither leaks into the
    definition. A ``!`` inside the quoted text is part of the text, not the start of a comment.

    Malformed values are read leniently rather than rejected: a value that is not quoted is returned as written
    (minus a trailing ``! comment``), and a quoted string that is never closed runs to the end of the value.

    Args:
        value: Everything after ``def:`` on the line.

    Returns:
        Tuple of the definition text and the list of dbxrefs (empty when none are given).

    Examples:
        >>> parse_obo_definition('"A process that \\\\"x\\\\" does." [GOC:mah, PMID:123]')
        ('A process that "x" does.', ['GOC:mah', 'PMID:123'])
    """
    value = value.strip()
    if not value.startswith('"'):
        return value.split(" ! ", 1)[0].strip(), []

    close = _closing_quote_index(value)
    if close < 0:
        return _unescape_obo(value[1:].rstrip().rstrip('"')), []

    text = _unescape_obo(value[1:close])
    rest = value[close + 1 :].lstrip()
    xrefs = _parse_obo_xref_list(rest) if rest.startswith("[") else []
    return text, xrefs


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
    header_metadata: Dict[str, Any] = {}

    # Read file content
    content = Path(path).read_text(encoding="utf-8")
    lines = content.splitlines()

    # Parse header
    header_metadata = _parse_header(lines)

    # Parse stanzas
    stanza_lines: List[str] = []
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
        elif stanza_type == "Typedef":
            typedef_info = _parse_typedef_stanza(stanza_lines)
            if typedef_info:
                header_metadata[f"typedef:{typedef_info['id']}"] = typedef_info
        elif stanza_type == "Instance":
            instance_info = _parse_instance_stanza(stanza_lines)
            if instance_info:
                header_metadata[f"instance:{instance_info['id']}"] = instance_info

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
    term_data: Dict[str, Any] = {}

    for line in stanza_lines[1:]:  # Skip [Term] line
        line = line.strip()
        if not line or line.startswith("!"):
            continue

        if ":" not in line:
            continue

        key, value = line.split(":", 1)
        key = key.strip()
        raw_value = value.strip()
        # Remove trailing comments
        if " ! " in value:
            value = value.split(" ! ", 1)[0]
        value = value.strip()

        # Handle multi-line values and qualifiers
        if key == "def":
            # The definition is a quoted string followed by a dbxref list; parse it from the raw value because a
            # " ! " inside the quotes is text, not a trailing comment.
            term_data["def"], definition_xrefs = parse_obo_definition(raw_value)
            if definition_xrefs:
                term_data["def_xrefs"] = definition_xrefs
        elif key in ["name", "namespace", "comment"]:
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
        elif key == "alt_id":
            if "alt_ids" not in term_data:
                term_data["alt_ids"] = []
            term_data["alt_ids"].append(value)
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

    if "alt_ids" in term_data:
        term.alt_ids = term_data["alt_ids"]

    # Evidence cited by the definition (the dbxref list that follows it), kept apart from the definition text
    if "def_xrefs" in term_data:
        term.metadata["definition_xrefs"] = term_data["def_xrefs"]

    # Store relationships for later processing
    term.metadata["relationships"] = term_data.get("relationships", [])

    return term


def _parse_typedef_stanza(stanza_lines: List[str]) -> Optional[Dict[str, Any]]:
    """Parse a [Typedef] stanza."""
    typedef_data: Dict[str, Any] = {}

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
            typedef_data["definition"], definition_xrefs = parse_obo_definition(value)
            if definition_xrefs:
                typedef_data["definition_xrefs"] = definition_xrefs

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

    for term in terms.values():
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
        content = Path(path).read_text(encoding="utf-8")
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

    content = Path(path).read_text(encoding="utf-8")
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
