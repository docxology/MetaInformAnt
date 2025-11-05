from __future__ import annotations

from pathlib import Path
from typing import Any, List, Optional, Set

from metainformant.core import io as core_io
from metainformant.core.errors import ValidationError
from metainformant.core.logging import get_logger
from metainformant.core.paths import expand_and_resolve, prepare_file_path

from .obo import parse_obo
from .query import get_roots, get_leaves
from .types import Ontology

logger = get_logger(__name__)


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
        
    Raises:
        IOError: If file does not exist or cannot be read
        ValidationError: If file is malformed or contains invalid data
        
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
    logger.info(f"Loading Gene Ontology from {path}")
    onto = parse_obo(path)
    logger.info(f"Loaded {onto.num_terms()} GO terms")
    return onto


def write_go_summary(onto: Ontology, dest: str | Path | None = None) -> Path:
    """Write a JSON summary of Gene Ontology statistics to file.
    
    Creates a summary document containing ontology metrics such as number
    of terms, namespaces, root terms, and depth information. Writes to
    output/ontology/go_summary.json by default.
    
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
        >>> "namespaces" in data
        True
    """
    if dest is None:
        dest = Path("output/ontology/go_summary.json")
    dest = expand_and_resolve(dest)
    prepare_file_path(dest)
    
    # Collect namespace statistics
    namespaces: dict[str, int] = {}
    for term in onto.terms.values():
        ns = term.namespace or "unknown"
        namespaces[ns] = namespaces.get(ns, 0) + 1
    
    # Get root terms (terms with no parents)
    root_terms = get_roots(onto)
    
    # Get leaf terms (terms with no children)
    leaf_terms = get_leaves(onto)
    
    # Calculate approximate depth (max path length from any root to any leaf)
    max_depth = 0
    from .query import path_to_root
    for leaf in list(leaf_terms)[:100]:  # Sample up to 100 leaves for performance
        try:
            path = path_to_root(onto, leaf)
            max_depth = max(max_depth, len(path) - 1)
        except Exception:
            pass
    
    summary = {
        "num_terms": onto.num_terms(),
        "namespaces": namespaces,
        "num_roots": len(root_terms),
        "num_leaves": len(leaf_terms),
        "max_depth": max_depth,
    }
    
    logger.info(f"Writing GO summary to {dest}")
    core_io.dump_json(summary, dest, indent=2)
    return dest


def validate_go_ontology(onto: Ontology) -> tuple[bool, List[str]]:
    """Validate Gene Ontology structure and integrity.
    
    Checks for common issues in GO ontologies:
    - Cycles in relationships
    - Orphaned terms
    - Missing required fields
    - Namespace consistency
    
    Args:
        onto: Ontology object to validate
        
    Returns:
        Tuple of (is_valid, list_of_errors)
        
    Examples:
        >>> onto = load_go_obo("go.obo")
        >>> is_valid, errors = validate_go_ontology(onto)
        >>> is_valid
        True
    """
    errors: List[str] = []
    
    # Use Ontology.validate() for basic integrity checks
    is_valid, validation_errors = onto.validate()
    errors.extend(validation_errors)
    
    # Check for expected GO namespaces
    expected_namespaces = {"biological_process", "molecular_function", "cellular_component"}
    found_namespaces = {term.namespace for term in onto.terms.values() if term.namespace}
    unexpected = found_namespaces - expected_namespaces
    if unexpected and len(expected_namespaces & found_namespaces) == 0:
        errors.append(f"Warning: No standard GO namespaces found. Found: {unexpected}")
    
    # Check for terms with empty names (should not happen in GO)
    empty_names = [tid for tid, term in onto.terms.items() if not term.name]
    if empty_names:
        errors.append(f"Found {len(empty_names)} terms with empty names")
    
    return len(errors) == 0, errors


def enrich_genes(
    genes: List[str],
    background: List[str] | None,
    onto: Ontology,
    gene_to_terms: dict[str, Set[str]],
) -> dict[str, Any]:
    """Gene enrichment analysis using Gene Ontology.
    
    **Note**: This function requires the `scipy` package for statistical tests.
    Install with: `pip install scipy`
    
    Performs Fisher's exact test or hypergeometric test for GO term enrichment.
    Requires gene-to-term mappings (typically from GAF/GPAD annotation files).
    
    Args:
        genes: List of gene identifiers to test for enrichment
        background: Optional background gene set. If None, uses all genes in
            gene_to_terms as background
        onto: Ontology object containing GO terms
        gene_to_terms: Dictionary mapping gene_id -> set of GO term IDs
        
    Returns:
        Dictionary with enrichment results:
        - 'terms': List of enriched GO terms with statistics
        - 'method': Statistical test used
        - 'correction': Multiple testing correction method
        
    Raises:
        ImportError: If scipy is not installed
        ValueError: If gene_to_terms is empty or invalid
        
    Examples:
        >>> # Requires scipy and gene-to-term mappings
        >>> gene_to_terms = {"GENE1": {"GO:0008150"}, "GENE2": {"GO:0008150"}}
        >>> genes = ["GENE1", "GENE2"]
        >>> onto = load_go_obo("go.obo")
        >>> # result = enrich_genes(genes, None, onto, gene_to_terms)
    
    Note:
        This function is a placeholder. Full implementation requires:
        - scipy for statistical tests (fisher_exact, hypergeom)
        - GAF/GPAD file parsers for gene-to-term mappings
        - Multiple testing correction (FDR, Bonferroni)
    """
    try:
        import scipy.stats  # noqa: F401
    except ImportError:
        raise ImportError(
            "scipy is required for enrichment analysis. "
            "Install with: pip install scipy"
        )
    
    if not gene_to_terms:
        raise ValueError("gene_to_terms cannot be empty")
    
    if not genes:
        raise ValueError("genes list cannot be empty")
    
    # Placeholder implementation
    logger.warning(
        "enrich_genes is a placeholder. Full implementation requires "
        "scipy.stats for statistical tests and GAF/GPAD parsers for "
        "gene-to-term mappings."
    )
    
    return {
        "terms": [],
        "method": "fisher_exact",
        "correction": "fdr",
        "note": "Placeholder implementation - requires scipy and annotation parsers",
    }


def semantic_similarity(
    onto: Ontology,
    term1: str,
    term2: str,
    method: str = "resnik",
    gene_to_terms: dict[str, Set[str]] | None = None,
) -> float:
    """Calculate semantic similarity between two GO terms.
    
    **Note**: This function requires the `scipy` package and gene annotation data.
    Install with: `pip install scipy`
    
    Implements semantic similarity measures such as Resnik, Lin, or Jiang-Conrath.
    Requires information content calculations from gene annotation frequencies.
    
    Args:
        onto: Ontology object containing GO terms
        term1: First GO term identifier
        term2: Second GO term identifier
        method: Similarity method: "resnik", "lin", or "jiang_conrath"
        gene_to_terms: Optional dictionary mapping gene_id -> set of GO term IDs.
            If None, information content cannot be calculated and returns 0.0.
            
    Returns:
        Similarity score between 0.0 and 1.0 (or -inf to inf for Resnik)
        
    Raises:
        ImportError: If scipy is not installed
        ValueError: If term1 or term2 not found, or invalid method
        
    Examples:
        >>> onto = load_go_obo("go.obo")
        >>> # Requires scipy and gene annotations
        >>> # similarity = semantic_similarity(onto, "GO:0008150", "GO:0009987")
    
    Note:
        This function is a placeholder. Full implementation requires:
        - scipy for mathematical operations
        - Information content calculation from gene annotation frequencies
        - Lowest common ancestor (LCA) finding
        - Gene annotation data (GAF/GPAD format)
    """
    try:
        import scipy  # noqa: F401
    except ImportError:
        raise ImportError(
            "scipy is required for semantic similarity. "
            "Install with: pip install scipy"
        )
    
    if not onto.has_term(term1):
        raise ValueError(f"Term {term1} not found in ontology")
    if not onto.has_term(term2):
        raise ValueError(f"Term {term2} not found in ontology")
    
    if method not in ("resnik", "lin", "jiang_conrath"):
        raise ValueError(f"Invalid method: {method}. Must be 'resnik', 'lin', or 'jiang_conrath'")
    
    if gene_to_terms is None:
        logger.warning(
            "semantic_similarity requires gene_to_terms for information content. "
            "Returning 0.0 as placeholder."
        )
        return 0.0
    
    # Placeholder implementation
    logger.warning(
        "semantic_similarity is a placeholder. Full implementation requires "
        "information content calculations and LCA finding."
    )
    
    return 0.0
