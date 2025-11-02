from __future__ import annotations

from pathlib import Path
from typing import Optional

from metainformant.core import io as core_io
from metainformant.core.paths import expand_and_resolve

from .obo import parse_obo
from .types import Ontology


def count_go_scripts(go_dir: Path) -> int:
    """Count Gene Ontology annotation script files in a directory.
    
    Args:
        go_dir: Directory path containing GO annotation files
        
    Returns:
        Integer count of annotation script files found
    """
    return sum(1 for p in go_dir.glob("*.py") if p.is_file())


def load_go_obo(path: str | Path) -> Ontology:
    """Load Gene Ontology from an OBO format file.
    
    Parses a Gene Ontology (GO) OBO file and constructs an Ontology object
    representing the GO term hierarchy and relationships.
    
    Args:
        path: Path to OBO format file (typically go.obo or go-basic.obo)
        
    Returns:
        Ontology object containing all GO terms and their relationships
        
    Examples:
        >>> onto = load_go_obo("data/go.obo")
        >>> onto.num_terms()
        45000
        >>> onto.has_term("GO:0008150")
        True
        
    Note:
        This is a lightweight reader suitable for standard GO OBO files.
        For advanced features (relationship types beyond is_a, complex qualifiers),
        consider specialized OBO parsing libraries.
    """
    return parse_obo(path)


def write_go_summary(onto: Ontology, dest: str | Path | None = None) -> Path:
    """Write a JSON summary of Gene Ontology statistics to file.
    
    Creates a summary document containing ontology metrics such as number
    of terms. Writes to output/ontology/go_summary.json by default.
    
    Args:
        onto: Ontology object to summarize
        dest: Optional destination path for summary file. If None, writes
            to output/ontology/go_summary.json
            
    Returns:
        Path to the created summary file
        
    Examples:
        >>> onto = load_go_obo("go-basic.obo")
        >>> summary_path = write_go_summary(onto)
        >>> summary_path.exists()
        True
        >>> import json
        >>> data = json.loads(summary_path.read_text())
        >>> "num_terms" in data
        True
    """
    if dest is None:
        dest = Path("output/ontology/go_summary.json")
    dest = expand_and_resolve(dest)
    core_io.dump_json({"num_terms": onto.num_terms()}, dest, indent=2)
    return dest
