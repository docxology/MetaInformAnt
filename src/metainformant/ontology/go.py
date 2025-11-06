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
    method: str = "fisher_exact",
    correction: str = "fdr",
    alpha: float = 0.05,
    propagate_annotations: bool = True,
) -> dict[str, Any]:
    """Gene enrichment analysis using Gene Ontology.
    
    **Note**: This function requires the `scipy` package for statistical tests.
    Install with: `pip install scipy`
    
    Performs Fisher's exact test for GO term enrichment. Annotations are
    automatically propagated to ancestor terms (if a gene is annotated with a term,
    it's also considered annotated with all parent terms).
    
    Args:
        genes: List of gene identifiers to test for enrichment
        background: Optional background gene set. If None, uses all genes in
            gene_to_terms as background
        onto: Ontology object containing GO terms
        gene_to_terms: Dictionary mapping gene_id -> set of GO term IDs
        method: Statistical test method ("fisher_exact" only currently supported)
        correction: Multiple testing correction ("fdr", "bonferroni", or "none")
        alpha: Significance level (default: 0.05)
        propagate_annotations: If True, propagate annotations to ancestor terms
        
    Returns:
        Dictionary with enrichment results:
        - 'terms': List of enriched GO terms with statistics, sorted by p-value
        - 'method': Statistical test used
        - 'correction': Multiple testing correction method
        - 'n_genes': Number of genes in query set
        - 'n_background': Number of genes in background
        - 'n_tests': Number of GO terms tested
        - 'n_significant': Number of significantly enriched terms
        
    Raises:
        ImportError: If scipy is not installed
        ValueError: If gene_to_terms is empty or invalid
        
    Examples:
        >>> # Requires scipy and gene-to-term mappings
        >>> gene_to_terms = {"GENE1": {"GO:0008150"}, "GENE2": {"GO:0008150"}}
        >>> genes = ["GENE1", "GENE2"]
        >>> onto = load_go_obo("go.obo")
        >>> result = enrich_genes(genes, None, onto, gene_to_terms)
        >>> len(result["terms"]) > 0
        True
    """
    try:
        import scipy.stats
    except ImportError:
        raise ImportError(
            "scipy is required for enrichment analysis. "
            "Install with: pip install scipy"
        )
    
    if not gene_to_terms:
        raise ValueError("gene_to_terms cannot be empty")
    
    if not genes:
        raise ValueError("genes list cannot be empty")
    
    if method != "fisher_exact":
        raise ValueError(f"Method '{method}' not supported. Use 'fisher_exact'")
    
    if correction not in ("fdr", "bonferroni", "none"):
        raise ValueError(f"Correction '{correction}' not supported. Use 'fdr', 'bonferroni', or 'none'")
    
    # Determine background
    if background is None:
        background = list(gene_to_terms.keys())
    
    # Convert to sets for efficient operations
    genes_set = set(genes)
    background_set = set(background)
    
    # Validate that all query genes are in background
    missing_genes = genes_set - background_set
    if missing_genes:
        logger.warning(f"Query genes not in background: {len(missing_genes)} genes")
        # Remove missing genes from query
        genes_set = genes_set & background_set
    
    if not genes_set:
        raise ValueError("No valid query genes found in background")
    
    # Propagate annotations to ancestors if requested
    from .query import ancestors
    
    gene_to_all_terms: dict[str, Set[str]] = {}
    for gene_id, terms in gene_to_terms.items():
        all_terms = set(terms)
        if propagate_annotations:
            # Add all ancestor terms for each annotation
            for term_id in terms:
                if onto.has_term(term_id):
                    ancestor_set = ancestors(onto, term_id)
                    all_terms.update(ancestor_set)
        gene_to_all_terms[gene_id] = all_terms
    
    # Count genes per term in query and background
    n_query = len(genes_set)
    n_background = len(background_set)
    
    # Get all terms that appear in background
    all_terms_in_background: Set[str] = set()
    for gene_id in background_set:
        if gene_id in gene_to_all_terms:
            all_terms_in_background.update(gene_to_all_terms[gene_id])
    
    # Filter to terms that exist in ontology
    terms_to_test = [term_id for term_id in all_terms_in_background if onto.has_term(term_id)]
    
    if not terms_to_test:
        logger.warning("No valid GO terms found in background annotations")
        return {
            "terms": [],
            "method": method,
            "correction": correction,
            "n_genes": n_query,
            "n_background": n_background,
            "n_tests": 0,
            "n_significant": 0,
        }
    
    # Perform enrichment test for each term
    enrichment_results = []
    pvalues = []
    
    for term_id in terms_to_test:
        # Count genes annotated with this term in query and background
        query_annotated = sum(1 for g in genes_set if g in gene_to_all_terms and term_id in gene_to_all_terms[g])
        query_not_annotated = n_query - query_annotated
        
        background_annotated = sum(1 for g in background_set if g in gene_to_all_terms and term_id in gene_to_all_terms[g])
        background_not_annotated = n_background - background_annotated
        
        # Skip if no genes annotated in background
        if background_annotated == 0:
            continue
        
        # Fisher's exact test: 2x2 contingency table
        #                Annotated  Not Annotated
        # Query         a          b
        # Background    c          d
        a = query_annotated
        b = query_not_annotated
        c = background_annotated - query_annotated
        d = background_not_annotated - query_not_annotated
        
        # Ensure all counts are non-negative
        if a < 0 or b < 0 or c < 0 or d < 0:
            continue
        
        # Perform Fisher's exact test
        try:
            table = [[a, b], [c, d]]
            odds_ratio, p_value = scipy.stats.fisher_exact(table, alternative="greater")
        except Exception as e:
            logger.debug(f"Fisher's exact test failed for {term_id}: {e}")
            continue
        
        # Get term information
        term = onto.terms.get(term_id)
        term_name = term.name if term else "Unknown"
        term_namespace = term.namespace if term else None
        
        enrichment_results.append({
            "term_id": term_id,
            "term_name": term_name,
            "namespace": term_namespace,
            "query_annotated": a,
            "query_total": n_query,
            "background_annotated": background_annotated,
            "background_total": n_background,
            "odds_ratio": float(odds_ratio),
            "p_value": float(p_value),
        })
        pvalues.append(float(p_value))
    
    if not enrichment_results:
        return {
            "terms": [],
            "method": method,
            "correction": correction,
            "n_genes": n_query,
            "n_background": n_background,
            "n_tests": 0,
            "n_significant": 0,
        }
    
    # Apply multiple testing correction
    if correction == "fdr":
        from metainformant.gwas.correction import fdr_correction
        fdr_result = fdr_correction(pvalues, alpha=alpha)
        corrected_pvalues = fdr_result.get("corrected_pvalues", pvalues)
        significant_indices = set(fdr_result.get("significant_indices", []))
    elif correction == "bonferroni":
        from metainformant.gwas.correction import bonferroni_correction
        bonf_result = bonferroni_correction(pvalues, alpha=alpha)
        corrected_alpha = bonf_result.get("corrected_alpha", alpha)
        significant_indices = {i for i, p in enumerate(pvalues) if p < corrected_alpha}
        corrected_pvalues = [min(p * len(pvalues), 1.0) for p in pvalues]
    else:  # none
        corrected_pvalues = pvalues
        significant_indices = {i for i, p in enumerate(pvalues) if p < alpha}
    
    # Add corrected p-values to results
    for i, result in enumerate(enrichment_results):
        result["corrected_p_value"] = corrected_pvalues[i]
        result["significant"] = i in significant_indices
    
    # Sort by p-value (ascending)
    enrichment_results.sort(key=lambda x: x["p_value"])
    
    return {
        "terms": enrichment_results,
        "method": method,
        "correction": correction,
        "n_genes": n_query,
        "n_background": n_background,
        "n_tests": len(enrichment_results),
        "n_significant": len(significant_indices),
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
    
    # Calculate information content for all terms
    from metainformant.information.semantic import information_content_from_annotations
    
    # Convert gene_to_terms format for information module
    term_annotations = {gene_id: terms for gene_id, terms in gene_to_terms.items()}
    term_ic = information_content_from_annotations(term_annotations)
    
    # Get term ICs
    ic1 = term_ic.get(term1, 0.0)
    ic2 = term_ic.get(term2, 0.0)
    
    if ic1 == 0.0 or ic2 == 0.0:
        logger.debug(f"Term {term1} or {term2} has zero information content")
        return 0.0
    
    # Find most informative common ancestor (MICA)
    from .query import common_ancestors
    
    common_ancs = common_ancestors(onto, term1, term2)
    
    if not common_ancs:
        return 0.0
    
    # Find MICA (term with highest IC among common ancestors)
    mica_ic = max(term_ic.get(anc, 0.0) for anc in common_ancs)
    
    if method == "resnik":
        # Resnik similarity: IC of MICA
        return float(mica_ic)
    elif method == "lin":
        # Lin similarity: 2 * IC(MICA) / (IC(term1) + IC(term2))
        if (ic1 + ic2) == 0:
            return 0.0
        return float(2.0 * mica_ic / (ic1 + ic2))
    elif method == "jiang_conrath":
        # Jiang-Conrath distance: IC(term1) + IC(term2) - 2 * IC(MICA)
        # Convert distance to similarity: 1 / (1 + distance)
        distance = ic1 + ic2 - 2.0 * mica_ic
        if distance < 0:
            distance = 0.0
        similarity = 1.0 / (1.0 + distance)
        return float(similarity)
    else:
        raise ValueError(f"Unknown method: {method}")
