"""AntWiki phenotype data loading utilities.

This module provides functions for loading and validating AntWiki phenotype data
from JSON files, including morphological measurements and behavioral traits.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ..core.errors import IOError as CoreIOError, ValidationError
from ..core.io import load_json
from ..core.logging import get_logger

logger = get_logger(__name__)


def _validate_antwiki_entry(entry: dict[str, Any], index: int) -> None:
    """Validate a single AntWiki entry has expected structure.
    
    Args:
        entry: Dictionary entry to validate
        index: Index of entry in list (for error messages)
        
    Raises:
        ValidationError: If entry structure is invalid
    """
    if not isinstance(entry, dict):
        raise ValidationError(f"Entry at index {index} is not a dictionary: {type(entry).__name__}")
    
    # Check for species identifier (either 'species' or 'taxon')
    if "species" not in entry and "taxon" not in entry:
        raise ValidationError(
            f"Entry at index {index} missing required field: expected 'species' or 'taxon'"
        )
    
    # Validate measurements if present
    if "measurements" in entry:
        if not isinstance(entry["measurements"], dict):
            raise ValidationError(
                f"Entry at index {index} has invalid 'measurements' field: expected dict, got {type(entry['measurements']).__name__}"
            )
    
    # Validate traits if present
    if "traits" in entry:
        if not isinstance(entry["traits"], list):
            raise ValidationError(
                f"Entry at index {index} has invalid 'traits' field: expected list, got {type(entry['traits']).__name__}"
            )


def load_antwiki_json(path: Path, validate: bool = True) -> list[dict[str, Any]]:
    """Load AntWiki phenotype data from JSON file.
    
    Parses JSON files containing AntWiki species phenotype information,
    including morphological measurements and behavioral traits. Uses core.io
    utilities for robust file handling with gzip support and proper error
    reporting.
    
    Args:
        path: Path to JSON file containing AntWiki data
        validate: If True, validate that entries have expected structure (species/taxon, measurements, traits)
        
    Returns:
        List of dictionaries, each representing a species or observation.
        If input is a single dictionary, returns list with one element.
        
    Raises:
        FileNotFoundError: If the file does not exist
        CoreIOError: If file read fails or JSON parsing fails
        ValidationError: If data structure is invalid (not list or dict) or entries don't have expected structure
        
    Examples:
        >>> from pathlib import Path
        >>> data = load_antwiki_json(Path("antwiki_data.json"))
        >>> len(data)
        150
        >>> data[0].keys()
        dict_keys(['species', 'measurements', 'traits'])
        
    Note:
        AntWiki provides comprehensive ant species phenotype databases.
        This loader supports both single-object and array JSON formats.
    """
    logger.debug(f"Loading AntWiki JSON from {path}")
    
    try:
        data = load_json(path)
    except FileNotFoundError:
        logger.error(f"AntWiki JSON file not found: {path}")
        raise
    except Exception as e:
        logger.error(f"Failed to load AntWiki JSON from {path}: {e}")
        raise CoreIOError(f"Failed to load JSON from {path}: {e}") from e
    
    # Validate data structure
    if isinstance(data, list):
        logger.debug(f"Loaded list format with {len(data)} entries")
        entries = data
    elif isinstance(data, dict):
        logger.debug("Loaded single dict format, converting to list")
        entries = [data]
    else:
        # Invalid data structure
        logger.error(f"Invalid data structure in {path}: expected list or dict, got {type(data).__name__}")
        raise ValidationError(
            f"Invalid AntWiki data format in {path}: expected list or dict, got {type(data).__name__}"
        )
    
    # Validate entry structures if requested
    if validate:
        logger.debug(f"Validating {len(entries)} entries")
        for i, entry in enumerate(entries):
            try:
                _validate_antwiki_entry(entry, i)
            except ValidationError as e:
                logger.error(f"Validation failed for entry {i}: {e}")
                raise
    
    logger.info(f"Successfully loaded {len(entries)} AntWiki entries from {path}")
    return entries
