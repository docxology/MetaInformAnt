"""Ontology serialization utilities.

Provides functions to save and load ontology objects to/from JSON format,
preserving all term attributes and relationships.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict

from metainformant.core import io as core_io
from metainformant.core.errors import ValidationError, IOError as CoreIOError
from metainformant.core.logging import get_logger
from metainformant.core.paths import expand_and_resolve, prepare_file_path

from .types import Ontology, Term

logger = get_logger(__name__)


def save_ontology(onto: Ontology, dest: str | Path) -> Path:
    """Save ontology to JSON file.
    
    Serializes an Ontology object to JSON format, preserving all term
    attributes including relationships, synonyms, xrefs, and subsets.
    
    Args:
        onto: Ontology object to serialize
        dest: Destination path for JSON file
        
    Returns:
        Path to the created JSON file
        
    Raises:
        IOError: If file cannot be written
        
    Examples:
        >>> onto = load_go_obo("go.obo")
        >>> path = save_ontology(onto, "go_saved.json")
        >>> path.exists()
        True
    """
    dest = expand_and_resolve(dest)
    prepare_file_path(dest)
    
    logger.info(f"Saving ontology with {onto.num_terms()} terms to {dest}")
    
    # Serialize terms
    terms_data: Dict[str, Dict[str, Any]] = {}
    for term_id, term in onto.terms.items():
        terms_data[term_id] = {
            "term_id": term.term_id,
            "name": term.name,
            "namespace": term.namespace,
            "definition": term.definition,
            "alt_ids": term.alt_ids,
            "is_a_parents": term.is_a_parents,
            "relationships": term.relationships,
            "synonyms": term.synonyms,
            "xrefs": term.xrefs,
            "subsets": term.subsets,
        }
    
    # Serialize parent/child mappings (convert sets to lists for JSON)
    parents_of_data = {
        term_id: list(parent_set) for term_id, parent_set in onto.parents_of.items()
    }
    children_of_data = {
        term_id: list(child_set) for term_id, child_set in onto.children_of.items()
    }
    
    data = {
        "terms": terms_data,
        "parents_of": parents_of_data,
        "children_of": children_of_data,
    }
    
    try:
        core_io.dump_json(data, dest, indent=2)
        logger.info(f"Successfully saved ontology to {dest}")
    except Exception as e:
        raise CoreIOError(f"Failed to save ontology to {dest}: {e}") from e
    
    return dest


def load_ontology(path: str | Path) -> Ontology:
    """Load ontology from JSON file.
    
    Deserializes an Ontology object from JSON format, reconstructing all
    term attributes and relationships.
    
    Args:
        path: Path to JSON file
        
    Returns:
        Ontology object reconstructed from JSON
        
    Raises:
        IOError: If file does not exist or cannot be read
        ValidationError: If file format is invalid
        
    Examples:
        >>> onto = save_ontology(original_onto, "go.json")
        >>> loaded = load_ontology("go.json")
        >>> loaded.num_terms() == original_onto.num_terms()
        True
    """
    path = expand_and_resolve(path)
    
    if not path.exists():
        raise CoreIOError(f"Ontology file not found: {path}")
    
    if not path.is_file():
        raise CoreIOError(f"Path is not a file: {path}")
    
    logger.info(f"Loading ontology from {path}")
    
    try:
        data = core_io.load_json(path)
    except Exception as e:
        raise ValidationError(f"Failed to read JSON from {path}: {e}") from e
    
    # Validate structure
    if not isinstance(data, dict):
        raise ValidationError(f"Invalid ontology JSON format: expected dict, got {type(data)}")
    
    if "terms" not in data:
        raise ValidationError("Invalid ontology JSON format: missing 'terms' key")
    
    # Reconstruct ontology
    onto = Ontology()
    terms_data = data.get("terms", {})
    
    for term_id, term_data in terms_data.items():
        if not isinstance(term_data, dict):
            raise ValidationError(f"Invalid term data for {term_id}: expected dict")
        
        # Extract term fields
        term = Term(
            term_id=term_data.get("term_id", term_id),
            name=term_data.get("name", ""),
            namespace=term_data.get("namespace"),
            definition=term_data.get("definition"),
            alt_ids=term_data.get("alt_ids", []),
            is_a_parents=term_data.get("is_a_parents", []),
            relationships=term_data.get("relationships", {}),
            synonyms=term_data.get("synonyms", []),
            xrefs=term_data.get("xrefs", []),
            subsets=term_data.get("subsets", []),
        )
        
        try:
            onto.add_term(term)
        except ValidationError as e:
            logger.warning(f"Skipping invalid term {term_id}: {e}")
            continue
    
    # Reconstruct parent/child mappings if present
    parents_of_data = data.get("parents_of", {})
    children_of_data = data.get("children_of", {})
    
    for term_id, parent_list in parents_of_data.items():
        if term_id in onto.terms:
            onto.parents_of[term_id] = set(parent_list)
    
    for term_id, child_list in children_of_data.items():
        if term_id in onto.terms:
            onto.children_of[term_id] = set(child_list)
    
    logger.info(f"Loaded {onto.num_terms()} terms from {path}")
    
    return onto

