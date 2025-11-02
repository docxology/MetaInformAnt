from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def load_antwiki_json(path: Path) -> list[dict[str, Any]]:
    """Load AntWiki phenotype data from JSON file.
    
    Parses JSON files containing AntWiki species phenotype information,
    including morphological measurements and behavioral traits.
    
    Args:
        path: Path to JSON file containing AntWiki data
        
    Returns:
        List of dictionaries, each representing a species or observation.
        If input is a single dictionary, returns list with one element.
        Returns empty list for invalid input.
        
    Examples:
        >>> data = load_antwiki_json(Path("antwiki_data.json"))
        >>> len(data)
        150
        >>> data[0].keys()
        dict_keys(['species', 'measurements', 'traits'])
        
    Note:
        AntWiki provides comprehensive ant species phenotype databases.
        This loader supports both single-object and array JSON formats.
    """
    data = json.loads(path.read_text())
    if isinstance(data, list):
        return data
    if isinstance(data, dict):
        return [data]
    return []
