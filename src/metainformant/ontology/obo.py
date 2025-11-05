from __future__ import annotations

from pathlib import Path
from typing import Iterable, Iterator, List

from metainformant.core.errors import ValidationError, IOError as CoreIOError
from metainformant.core.logging import get_logger
from metainformant.core.paths import expand_and_resolve

from .types import Ontology, Term

logger = get_logger(__name__)


def _iter_stanzas(lines: Iterable[str]) -> Iterator[List[str]]:
    stanza: List[str] = []
    for raw in lines:
        line = raw.strip("\n")
        if not line:
            # keep blank lines inside stanza
            stanza.append("")
            continue
        if line.startswith("[") and line.endswith("]"):
            # yield previous stanza if any
            if stanza:
                yield stanza
                stanza = []
            stanza.append(line)
        else:
            stanza.append(line)
    if stanza:
        yield stanza


def parse_obo(path: str | Path) -> Ontology:
    """Parse OBO (Open Biological and Biomedical Ontologies) format file.
    
    Parses a subset of OBO format sufficient for Gene Ontology and similar
    hierarchical ontologies. Supports essential fields for term representation
    and hierarchy traversal.
    
    Args:
        path: Path to OBO format file (typically .obo extension)
        
    Returns:
        Ontology object containing all parsed terms and their relationships.
        Each term includes identifier, name, namespace, definition, alternative
        IDs, and parent relationships (is_a).
        
    Raises:
        IOError: If file does not exist or cannot be read
        ValidationError: If file is malformed or contains invalid data
        
    Supported OBO Fields:
        - id: Term identifier (e.g., "GO:0008150")
        - name: Human-readable term name
        - namespace: Ontology namespace/category
        - def: Term definition (with optional citations)
        - alt_id: Alternative identifiers for the term
        - is_a: Parent term relationships (creates hierarchy)
        
    Examples:
        >>> onto = parse_obo("go-basic.obo")
        >>> onto.num_terms()
        45000
        >>> onto.has_term("GO:0008150")
        True
        >>> term = onto.terms["GO:0008150"]
        >>> term.name
        'biological_process'
        
    Note:
        This is a lightweight parser for standard OBO files. Complex
        relationship types beyond is_a, qualifiers, and advanced OBO
        features may require specialized parsing libraries.
        
    References:
        OBO Format Specification: https://owlcollab.github.io/oboformat/doc/GO.format.html
    """
    p = expand_and_resolve(path)
    
    if not p.exists():
        raise CoreIOError(f"OBO file not found: {path}")
    
    if not p.is_file():
        raise CoreIOError(f"Path is not a file: {path}")
    
    if p.stat().st_size == 0:
        raise ValidationError(f"OBO file is empty: {path}")
    
    logger.info(f"Parsing OBO file: {path} (size: {p.stat().st_size} bytes)")
    
    ontology = Ontology()
    line_num = 0
    term_count = 0
    errors: List[str] = []
    
    try:
        with p.open("rt", encoding="utf-8") as fh:
            for stanza in _iter_stanzas(fh):
                line_num += len(stanza)
                if not stanza:
                    continue
                if stanza[0] != "[Term]":
                    continue

                term_id: str | None = None
                name: str | None = None
                namespace: str | None = None
                definition: str | None = None
                alt_ids: List[str] = []
                parents: List[str] = []
                relationships: dict[str, List[str]] = {}
                synonyms: List[str] = []
                xrefs: List[str] = []
                subsets: List[str] = []

                for line in stanza[1:]:
                    if not line:
                        continue
                    if line.startswith("id: "):
                        term_id = line[4:].strip()
                    elif line.startswith("name: "):
                        name = line[6:].strip()
                    elif line.startswith("namespace: "):
                        namespace = line[11:].strip()
                    elif line.startswith("def: "):
                        definition = line[5:].strip()
                    elif line.startswith("alt_id: "):
                        alt_ids.append(line[8:].strip())
                    elif line.startswith("is_a: "):
                        # Extract term ID (before space or !)
                        parent_line = line[6:].strip()
                        parent_id = parent_line.split()[0].split("!")[0].strip()
                        parents.append(parent_id)
                    elif line.startswith("synonym: "):
                        # Parse synonym: "text" [SYNONYM_TYPE] [optional citations]
                        synonym_line = line[9:].strip()
                        # Extract quoted text (first quoted string)
                        if '"' in synonym_line:
                            start = synonym_line.find('"') + 1
                            end = synonym_line.find('"', start)
                            if end > start:
                                synonym_text = synonym_line[start:end]
                                synonyms.append(synonym_text)
                    elif line.startswith("xref: "):
                        # Parse xref: database:identifier or database:identifier "description"
                        xref_line = line[6:].strip()
                        # Extract database:identifier (before space or quote)
                        xref_parts = xref_line.split()[0].split('"')[0].strip()
                        if xref_parts:
                            xrefs.append(xref_parts)
                    elif line.startswith("subset: "):
                        # Parse subset: subset_name
                        subset_line = line[8:].strip()
                        subset_name = subset_line.split()[0].strip()
                        if subset_name:
                            subsets.append(subset_name)
                    elif ": " in line and not line.startswith("comment: "):
                        # Parse other relationship types (part_of, regulates, has_part, etc.)
                        # Format: relationship_type: term_id [optional description]
                        parts = line.split(": ", 1)
                        if len(parts) == 2:
                            rel_type = parts[0].strip()
                            rel_value = parts[1].strip()
                            # Extract term ID (before space or !)
                            rel_term_id = rel_value.split()[0].split("!")[0].strip()
                            if rel_type not in ("id", "name", "namespace", "def", "alt_id", "is_a", "synonym", "xref", "subset", "comment"):
                                # This is a relationship type
                                if rel_type not in relationships:
                                    relationships[rel_type] = []
                                relationships[rel_type].append(rel_term_id)

                if term_id and name:
                    try:
                        term = Term(
                            term_id=term_id,
                            name=name,
                            namespace=namespace,
                            definition=definition,
                            alt_ids=alt_ids,
                            is_a_parents=parents,
                            relationships=relationships,
                            synonyms=synonyms,
                            xrefs=xrefs,
                            subsets=subsets,
                        )
                        ontology.add_term(term)
                        term_count += 1
                    except ValidationError as e:
                        errors.append(f"Line ~{line_num}: {e}")
                elif term_id or name:
                    errors.append(f"Line ~{line_num}: Term missing required field (id or name)")
                    if term_id:
                        errors.append(f"  Term ID: {term_id}")
                    if name:
                        errors.append(f"  Name: {name}")
        
        if term_count == 0:
            raise ValidationError(f"No valid terms found in OBO file: {path}")
        
        logger.info(f"Parsed {term_count} terms from {path}")
        
        if errors:
            logger.warning(f"Encountered {len(errors)} parsing errors (non-fatal)")
            for error in errors[:10]:  # Log first 10 errors
                logger.warning(f"  {error}")
            if len(errors) > 10:
                logger.warning(f"  ... and {len(errors) - 10} more errors")
        
    except UnicodeDecodeError as e:
        raise ValidationError(f"OBO file encoding error in {path}: {e}") from e
    except Exception as e:
        raise ValidationError(f"Error parsing OBO file {path} at line ~{line_num}: {e}") from e
    
    return ontology
