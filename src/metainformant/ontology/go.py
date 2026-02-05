"""Gene Ontology analysis and enrichment.

This module provides functionality for working with Gene Ontology (GO) data,
including enrichment analysis, semantic similarity calculations, and GO-specific utilities.
"""

from __future__ import annotations

import math
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import pandas as pd

from metainformant.core import io, logging, paths, validation

from .obo import parse_obo
from .query import ancestors, descendants
from .types import Ontology, Term, create_ontology

logger = logging.get_logger(__name__)


def count_go_scripts(go_dir: Path) -> int:
    """Count the number of GO analysis scripts in a directory.

    Args:
        go_dir: Directory containing GO analysis scripts.

    Returns:
        Number of GO-related script files.

    Examples:
        >>> count_go_scripts(Path("scripts/go_analysis"))
        15
    """
    go_dir = validation.validate_path_is_dir(go_dir)

    script_extensions = [".py", ".r", ".sh", ".pl"]
    go_files = []

    for ext in script_extensions:
        go_files.extend(paths.find_files_by_extension(go_dir, ext.lstrip(".")))

    # Filter for GO-related files
    go_related = [f for f in go_files if "go" in f.name.lower() or "ontology" in f.name.lower()]

    logger.info(f"Found {len(go_related)} GO-related scripts in {go_dir}")
    return len(go_related)


def load_go_obo(path: str | Path) -> Ontology:
    """Load Gene Ontology from OBO file.

    Args:
        path: Path to the GO OBO file.

    Returns:
        Ontology object containing GO terms and relationships.

    Raises:
        FileNotFoundError: If the OBO file does not exist.
        ValueError: If the ontology is not a valid GO ontology.

    Examples:
        >>> go_ontology = load_go_obo("data/go.obo")
        >>> len(go_ontology)
        50000
    """
    ontology = parse_obo(path)

    # Validate it's a GO ontology
    if not any("Gene Ontology" in str(value) for value in ontology.metadata.values()):
        logger.warning("Loaded ontology may not be Gene Ontology - missing GO metadata")

    logger.info(f"Loaded GO ontology with {len(ontology)} terms")
    return ontology


def write_go_summary(onto: Ontology, dest: str | Path | None = None) -> Path:
    """Write a summary of Gene Ontology content.

    Args:
        onto: GO ontology to summarize.
        dest: Destination path for the summary file. If None, writes to output/go_summary.txt.

    Returns:
        Path to the written summary file.

    Examples:
        >>> summary_path = write_go_summary(go_ontology)
        >>> print(f"Summary written to: {summary_path}")
    """
    if dest is None:
        output_dir = paths.ensure_directory(Path("output"))
        dest = output_dir / "go_summary.txt"
    else:
        dest = Path(dest)
        dest.parent.mkdir(parents=True, exist_ok=True)

    # Count terms by namespace
    namespace_counts = defaultdict(int)
    obsolete_count = 0

    for term in onto.terms.values():
        if term.namespace:
            namespace_counts[term.namespace] += 1
        if term.is_obsolete:
            obsolete_count += 1

    # Count relationships by type
    relationship_counts = defaultdict(int)
    for rel in onto.relationships:
        relationship_counts[rel.relation_type] += 1

    # Get roots and leaves
    roots = onto.get_roots()
    leaves = onto.get_leaves()

    # Write summary
    content = f"""Gene Ontology Summary
{'='*50}

Basic Statistics:
- Total terms: {len(onto)}
- Obsolete terms: {obsolete_count}
- Active terms: {len(onto) - obsolete_count}
- Root terms: {len(roots)}
- Leaf terms: {len(leaves)}
- Relationships: {len(onto.relationships)}

Terms by Namespace:
"""

    for namespace, count in sorted(namespace_counts.items()):
        content += f"- {namespace}: {count}\n"

    content += f"\nRelationships by Type:\n"
    for rel_type, count in sorted(relationship_counts.items()):
        content += f"- {rel_type}: {count}\n"

    content += f"\nRoot Terms:\n"
    for root in sorted(list(roots)[:10]):  # Show first 10 roots
        term = onto.get_term(root)
        if term:
            content += f"- {root}: {term.name}\n"
    if len(roots) > 10:
        content += f"... and {len(roots) - 10} more\n"

    # Write to file
    io.write_text(content, dest)
    logger.info(f"GO summary written to {dest}")

    return dest


def validate_go_ontology(onto: Ontology) -> tuple[bool, List[str]]:
    """Validate Gene Ontology structure and content.

    Args:
        onto: Ontology to validate.

    Returns:
        Tuple of (is_valid, list_of_issues).

    Examples:
        >>> is_valid, issues = validate_go_ontology(go_ontology)
        >>> if not is_valid:
        ...     print("Validation issues:", issues)
    """
    issues = []

    # Check for required namespaces
    required_namespaces = {"biological_process", "molecular_function", "cellular_component"}
    found_namespaces = set()

    for term in onto.terms.values():
        if term.namespace:
            found_namespaces.add(term.namespace)

    missing_namespaces = required_namespaces - found_namespaces
    if missing_namespaces:
        issues.append(f"Missing required GO namespaces: {missing_namespaces}")

    # Check for root terms
    roots = onto.get_roots()
    if len(roots) < 3:
        issues.append(f"Expected at least 3 root terms, found {len(roots)}")

    # Check for relationships
    if len(onto.relationships) == 0:
        issues.append("No relationships found in ontology")

    # Check for GO ID format
    go_id_pattern = re.compile(r"GO:\d{7}")
    invalid_ids = []

    for term_id in onto.terms.keys():
        if not go_id_pattern.match(term_id):
            invalid_ids.append(term_id)

    if invalid_ids:
        issues.append(f"Found {len(invalid_ids)} terms with invalid GO IDs (first 5: {invalid_ids[:5]})")

    # Check for orphan terms (terms not connected to roots)
    if roots:
        connected_terms = set()
        for root in roots:
            connected_terms.update(onto.get_descendants(root))
            connected_terms.update(onto.get_ancestors(root))
        connected_terms.update(roots)

        orphan_terms = set(onto.terms.keys()) - connected_terms
        if orphan_terms:
            issues.append(f"Found {len(orphan_terms)} orphan terms not connected to root terms")

    is_valid = len(issues) == 0

    if is_valid:
        logger.info("GO ontology validation passed")
    else:
        logger.warning(f"GO ontology validation found {len(issues)} issues")

    return is_valid, issues


def enrich_genes(
    genes: List[str],
    background: List[str] | None,
    annotations: Dict[str, Set[str]],
    alpha: float = 0.05,
    method: str = "fisher",
) -> pd.DataFrame:
    """Perform Gene Ontology enrichment analysis.

    Args:
        genes: List of gene identifiers to test for enrichment.
        background: Background gene list. If None, uses all genes in annotations.
        annotations: Dictionary mapping GO term IDs to sets of gene IDs.
        alpha: Significance threshold for enrichment.
        method: Statistical method ("fisher" or "hypergeometric").

    Returns:
        DataFrame with enrichment results sorted by p-value.

    Raises:
        ValueError: If method is not supported or inputs are invalid.

    Examples:
        >>> results = enrich_genes(
        ...     genes=['GENE1', 'GENE2', 'GENE3'],
        ...     background=None,
        ...     annotations=go_annotations,
        ...     alpha=0.05,
        ...     method='fisher'
        ... )
        >>> results.head()
    """
    validation.validate_not_empty(genes, "genes")
    validation.validate_not_none(annotations, "annotations")

    if method not in ["fisher", "hypergeometric"]:
        raise ValueError(f"Unsupported method: {method}. Use 'fisher' or 'hypergeometric'.")

    # Prepare background
    if background is None:
        background = set()
        for gene_set in annotations.values():
            background.update(gene_set)
        background = list(background)

    background_set = set(background)
    gene_set = set(genes)

    # Validate gene sets
    invalid_genes = gene_set - background_set
    if invalid_genes:
        logger.warning(f"{len(invalid_genes)} genes not found in background: {list(invalid_genes)[:5]}...")

    # Remove genes not in background
    gene_set = gene_set & background_set

    if len(gene_set) == 0:
        logger.warning("No valid genes remaining after filtering against background")
        return pd.DataFrame()

    logger.info(f"Testing enrichment for {len(gene_set)} genes against {len(background_set)} background genes")

    results = []

    for go_term, annotated_genes in annotations.items():
        annotated_set = set(annotated_genes)

        # Calculate contingency table
        # [[in_query_and_term, in_query_not_term],
        #  [not_in_query_and_term, not_in_query_not_term]]

        in_query_and_term = len(gene_set & annotated_set)
        in_query_not_term = len(gene_set - annotated_set)
        not_in_query_and_term = len((background_set - gene_set) & annotated_set)
        not_in_query_not_term = len((background_set - gene_set) - annotated_set)

        # Skip if no genes in this term
        if in_query_and_term == 0:
            continue

        # Calculate enrichment
        if method == "fisher":
            p_value = _fisher_exact_test(
                in_query_and_term, in_query_not_term, not_in_query_and_term, not_in_query_not_term
            )
        else:  # hypergeometric
            p_value = _hypergeometric_test(in_query_and_term, len(annotated_set), len(gene_set), len(background_set))

        # Calculate enrichment ratio
        expected = (len(annotated_set) * len(gene_set)) / len(background_set)
        enrichment = in_query_and_term / expected if expected > 0 else float("inf")

        results.append(
            {
                "go_term": go_term,
                "p_value": p_value,
                "enrichment": enrichment,
                "observed": in_query_and_term,
                "expected": expected,
                "term_size": len(annotated_set),
                "query_size": len(gene_set),
            }
        )

    # Convert to DataFrame and sort
    df = pd.DataFrame(results)
    if len(df) > 0:
        df = df.sort_values("p_value")

        # Apply multiple testing correction (Bonferroni)
        df["bonferroni_p"] = df["p_value"] * len(df)
        df["bonferroni_p"] = df["bonferroni_p"].clip(upper=1.0)

        # Filter by significance
        significant = df[df["bonferroni_p"] <= alpha].copy()
        logger.info(f"Found {len(significant)} significantly enriched GO terms (alpha={alpha})")
    else:
        logger.warning("No enrichment results generated")

    return df


def semantic_similarity(
    term1: str, term2: str, term_ic: Dict[str, float], hierarchy: Dict[str, Set[str]], method: str = "resnik"
) -> float:
    """Calculate semantic similarity between two GO terms.

    Args:
        term1: First GO term ID.
        term2: Second GO term ID.
        term_ic: Dictionary mapping term IDs to information content values.
        hierarchy: Dictionary mapping term IDs to sets of ancestor terms.
        method: Similarity method ("resnik", "lin", "jiang-conrath").

    Returns:
        Similarity score between 0 and 1.

    Raises:
        ValueError: If method is not supported or terms are invalid.

    Examples:
        >>> similarity = semantic_similarity(
        ...     "GO:0008150", "GO:0009987",
        ...     term_ic=ic_dict,
        ...     hierarchy=hierarchy_dict,
        ...     method="resnik"
        ... )
    """
    if method not in ["resnik", "lin", "jiang-conrath"]:
        raise ValueError(f"Unsupported method: {method}")

    if term1 not in term_ic or term2 not in term_ic:
        logger.warning(f"One or both terms not found in IC dictionary: {term1}, {term2}")
        return 0.0

    # Find most informative common ancestor (MICA)
    ancestors1 = hierarchy.get(term1, set()) | {term1}
    ancestors2 = hierarchy.get(term2, set()) | {term2}

    common_ancestors = ancestors1 & ancestors2

    if not common_ancestors:
        return 0.0

    # Find MICA
    mica = max(common_ancestors, key=lambda x: term_ic.get(x, 0.0))
    mica_ic = term_ic.get(mica, 0.0)

    if method == "resnik":
        # Resnik similarity: IC of MICA
        similarity = mica_ic

    elif method == "lin":
        # Lin similarity: 2 * IC(MICA) / (IC(term1) + IC(term2))
        ic1 = term_ic[term1]
        ic2 = term_ic[term2]
        if ic1 + ic2 > 0:
            similarity = 2 * mica_ic / (ic1 + ic2)
        else:
            similarity = 0.0

    else:  # jiang-conrath
        # Jiang-Conrath distance: IC(term1) + IC(term2) - 2 * IC(MICA)
        ic1 = term_ic[term1]
        ic2 = term_ic[term2]
        distance = ic1 + ic2 - 2 * mica_ic

        # Convert distance to similarity (assuming max distance is 2 * max IC)
        max_ic = max(term_ic.values()) if term_ic else 1.0
        max_distance = 2 * max_ic
        similarity = 1.0 - (distance / max_distance) if max_distance > 0 else 0.0
        similarity = max(0.0, min(1.0, similarity))  # Clamp to [0, 1]

    return similarity


def _fisher_exact_test(a: int, b: int, c: int, d: int) -> float:
    """Calculate Fisher's exact test p-value."""
    # Simplified implementation - in practice, would use scipy.stats.fisher_exact
    try:
        from scipy.stats import fisher_exact

        _, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")
        return p_value
    except ImportError:
        # Fallback calculation using hypergeometric distribution
        return _hypergeometric_test(a, a + c, a + b, a + b + c + d)


def _hypergeometric_test(k: int, K: int, n: int, N: int) -> float:
    """Calculate hypergeometric test p-value."""
    try:
        from scipy.stats import hypergeom

        # P(X >= k) where X ~ Hypergeometric(N, K, n)
        p_value = hypergeom.sf(k - 1, N, K, n)
        return p_value
    except ImportError:
        # Very simplified approximation
        from math import exp, lgamma

        try:
            # Log probability calculation
            log_p = (
                lgamma(K + 1)
                + lgamma(N - K + 1)
                + lgamma(n + 1)
                + lgamma(N - n + 1)
                - lgamma(N + 1)
                - lgamma(k + 1)
                - lgamma(K - k + 1)
                - lgamma(n - k + 1)
                - lgamma(N - K - n + k + 1)
            )
            return min(1.0, exp(log_p))
        except (ValueError, OverflowError):
            return 1.0  # Conservative fallback


def calculate_term_ic(onto: Ontology, term_counts: Dict[str, int], total_annotations: int) -> Dict[str, float]:
    """Calculate information content for ontology terms.

    Args:
        onto: Ontology object.
        term_counts: Dictionary mapping term IDs to annotation counts.
        total_annotations: Total number of annotations across all terms.

    Returns:
        Dictionary mapping term IDs to information content values.

    Examples:
        >>> term_ic = calculate_term_ic(ontology, annotation_counts, total_genes)
        >>> term_ic['GO:0008150']  # Information content of biological_process
    """
    term_ic = {}

    for term_id in onto.terms.keys():
        count = term_counts.get(term_id, 0)
        if count > 0 and total_annotations > 0:
            # IC = -log2(p) where p is probability of term annotation
            p = count / total_annotations
            ic = -math.log2(p)
            term_ic[term_id] = ic
        else:
            term_ic[term_id] = 0.0

    logger.info(f"Calculated IC for {len(term_ic)} terms")
    return term_ic


def build_hierarchy_dict(onto: Ontology) -> Dict[str, Set[str]]:
    """Build hierarchy dictionary for semantic similarity calculations.

    Args:
        onto: Ontology object.

    Returns:
        Dictionary mapping term IDs to sets of ancestor terms.

    Examples:
        >>> hierarchy = build_hierarchy_dict(ontology)
        >>> ancestors = hierarchy['GO:0008150']  # Ancestors of biological_process
    """
    hierarchy = {}

    for term_id in onto.terms.keys():
        ancestors_set = onto.get_ancestors(term_id)
        hierarchy[term_id] = ancestors_set

    logger.info(f"Built hierarchy dictionary for {len(hierarchy)} terms")
    return hierarchy


# Import re at module level for GO ID validation
import re
