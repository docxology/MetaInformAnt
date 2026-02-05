"""Taxonomic classification for amplicon sequences.

Implements naive Bayes and alignment-based taxonomic classification methods
for assigning taxonomy to OTU/ASV representative sequences. Includes
hierarchical taxonomy tree construction and bootstrap confidence estimation.

Taxonomic ranks follow the standard hierarchy:
    Domain > Phylum > Class > Order > Family > Genus > Species
"""

from __future__ import annotations

import math
import os
import random
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

_ENV_PREFIX = "META_"

TAXONOMIC_RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]


@dataclass
class TaxonomyAssignment:
    """Taxonomic assignment for a single sequence.

    Attributes:
        sequence_id: Identifier of the classified sequence.
        lineage: Full taxonomic lineage as list of (rank, name) tuples.
        confidence: Confidence scores per rank (0.0 to 1.0).
        method: Classification method used.
        best_hit_id: ID of the best matching reference (if alignment-based).
        best_hit_identity: Sequence identity to best hit.
    """

    sequence_id: str
    lineage: list[tuple[str, str]]
    confidence: dict[str, float] = field(default_factory=dict)
    method: str = "naive_bayes"
    best_hit_id: str | None = None
    best_hit_identity: float = 0.0


@dataclass
class TaxonomyNode:
    """Node in a hierarchical taxonomy tree.

    Attributes:
        name: Taxon name at this node.
        rank: Taxonomic rank (domain, phylum, etc.).
        count: Number of sequences assigned to this taxon.
        children: Child nodes in the hierarchy.
        parent: Parent node reference (or None for root).
    """

    name: str
    rank: str
    count: int = 0
    children: dict[str, TaxonomyNode] = field(default_factory=dict)
    parent: TaxonomyNode | None = field(default=None, repr=False)

    def add_child(self, name: str, rank: str) -> TaxonomyNode:
        """Add or get a child node."""
        if name not in self.children:
            child = TaxonomyNode(name=name, rank=rank, parent=self)
            self.children[name] = child
        return self.children[name]

    def get_lineage(self) -> list[tuple[str, str]]:
        """Get full lineage from root to this node."""
        lineage: list[tuple[str, str]] = []
        node: TaxonomyNode | None = self
        while node is not None and node.rank != "root":
            lineage.append((node.rank, node.name))
            node = node.parent
        return list(reversed(lineage))

    def total_descendants(self) -> int:
        """Count total sequences in this subtree."""
        total = self.count
        for child in self.children.values():
            total += child.total_descendants()
        return total

    def to_dict(self) -> dict[str, Any]:
        """Convert tree to nested dictionary representation."""
        result: dict[str, Any] = {
            "name": self.name,
            "rank": self.rank,
            "count": self.count,
            "total": self.total_descendants(),
        }
        if self.children:
            result["children"] = {name: child.to_dict() for name, child in self.children.items()}
        return result


def _build_kmer_profiles(
    reference_db: dict[str, str],
    k: int = 8,
) -> dict[str, set[str]]:
    """Build k-mer profiles for reference sequences.

    Args:
        reference_db: Mapping of reference IDs to sequences.
        k: K-mer size.

    Returns:
        Mapping of reference IDs to sets of k-mers.
    """
    profiles: dict[str, set[str]] = {}
    for ref_id, seq in reference_db.items():
        seq_upper = seq.upper()
        kmers = set()
        for i in range(len(seq_upper) - k + 1):
            kmer = seq_upper[i : i + k]
            if all(c in "ACGT" for c in kmer):
                kmers.add(kmer)
        profiles[ref_id] = kmers
    return profiles


def _naive_bayes_classify(
    query_seq: str,
    reference_db: dict[str, str],
    reference_taxonomy: dict[str, list[tuple[str, str]]],
    k: int = 8,
    bootstrap_n: int = 100,
    bootstrap_threshold: float = 0.8,
) -> TaxonomyAssignment:
    """Classify a query sequence using naive Bayes k-mer approach.

    The RDP-style naive Bayes classifier:
    1. Extract all k-mers from the query sequence.
    2. For each reference taxon at each rank, compute the posterior probability
       based on k-mer frequencies.
    3. Assign the query to the taxon with highest posterior at each rank.
    4. Use bootstrap resampling of k-mers for confidence estimation.

    Args:
        query_seq: Query nucleotide sequence.
        reference_db: Reference sequences keyed by ID.
        reference_taxonomy: Taxonomy for each reference sequence.
        k: K-mer size for classification.
        bootstrap_n: Number of bootstrap iterations for confidence.
        bootstrap_threshold: Minimum bootstrap ratio for confident assignment.

    Returns:
        TaxonomyAssignment with lineage and per-rank confidence scores.
    """
    query_upper = query_seq.upper()
    query_kmers: list[str] = []
    for i in range(len(query_upper) - k + 1):
        kmer = query_upper[i : i + k]
        if all(c in "ACGT" for c in kmer):
            query_kmers.append(kmer)

    if not query_kmers:
        return TaxonomyAssignment(
            sequence_id="",
            lineage=[(rank, "unclassified") for rank in TAXONOMIC_RANKS],
            confidence={rank: 0.0 for rank in TAXONOMIC_RANKS},
            method="naive_bayes",
        )

    # Build per-taxon k-mer probability models at each rank
    # For each rank, group references by their taxon at that rank
    rank_taxon_kmers: dict[str, dict[str, Counter[str]]] = {}
    for rank in TAXONOMIC_RANKS:
        rank_taxon_kmers[rank] = defaultdict(Counter)

    for ref_id, lineage in reference_taxonomy.items():
        if ref_id not in reference_db:
            continue
        ref_seq = reference_db[ref_id].upper()
        ref_kmer_set: list[str] = []
        for i in range(len(ref_seq) - k + 1):
            kmer = ref_seq[i : i + k]
            if all(c in "ACGT" for c in kmer):
                ref_kmer_set.append(kmer)

        lineage_dict = {rank: name for rank, name in lineage}
        for rank in TAXONOMIC_RANKS:
            taxon_name = lineage_dict.get(rank, "unclassified")
            for kmer in ref_kmer_set:
                rank_taxon_kmers[rank][taxon_name][kmer] += 1

    # Classify at each rank using log-likelihood
    assigned_lineage: list[tuple[str, str]] = []
    confidence_scores: dict[str, float] = {}

    for rank in TAXONOMIC_RANKS:
        taxon_models = rank_taxon_kmers[rank]
        if not taxon_models:
            assigned_lineage.append((rank, "unclassified"))
            confidence_scores[rank] = 0.0
            continue

        # Compute total k-mers per taxon for probability estimation
        taxon_totals: dict[str, int] = {}
        all_kmers_vocab: set[str] = set()
        for taxon, kmer_counts in taxon_models.items():
            taxon_totals[taxon] = sum(kmer_counts.values())
            all_kmers_vocab.update(kmer_counts.keys())

        vocab_size = len(all_kmers_vocab) if all_kmers_vocab else 1

        # Log-likelihood for each taxon
        taxon_scores: dict[str, float] = {}
        total_refs = sum(taxon_totals.values())

        for taxon, kmer_counts in taxon_models.items():
            # Prior: fraction of references in this taxon
            prior = math.log(taxon_totals[taxon] / total_refs) if total_refs > 0 else 0.0
            log_likelihood = prior
            total_t = taxon_totals[taxon]

            for qk in query_kmers:
                # Laplace-smoothed probability
                count = kmer_counts.get(qk, 0)
                prob = (count + 0.5) / (total_t + 0.5 * vocab_size)
                log_likelihood += math.log(prob)

            taxon_scores[taxon] = log_likelihood

        # Best taxon
        best_taxon = max(taxon_scores, key=taxon_scores.get)  # type: ignore[arg-type]

        # Bootstrap confidence: resample k-mers and reclassify
        bootstrap_hits = 0
        if bootstrap_n > 0 and len(query_kmers) > 0:
            for _ in range(bootstrap_n):
                # Resample k-mers with replacement
                resampled = random.choices(query_kmers, k=len(query_kmers))
                resample_scores: dict[str, float] = {}
                for taxon, kmer_counts in taxon_models.items():
                    prior = math.log(taxon_totals[taxon] / total_refs) if total_refs > 0 else 0.0
                    ll = prior
                    total_t = taxon_totals[taxon]
                    for qk in resampled:
                        count = kmer_counts.get(qk, 0)
                        prob = (count + 0.5) / (total_t + 0.5 * vocab_size)
                        ll += math.log(prob)
                    resample_scores[taxon] = ll
                resample_best = max(resample_scores, key=resample_scores.get)  # type: ignore[arg-type]
                if resample_best == best_taxon:
                    bootstrap_hits += 1
            confidence = bootstrap_hits / bootstrap_n
        else:
            confidence = 0.0

        assigned_lineage.append((rank, best_taxon))
        confidence_scores[rank] = confidence

    return TaxonomyAssignment(
        sequence_id="",
        lineage=assigned_lineage,
        confidence=confidence_scores,
        method="naive_bayes",
    )


def _blast_classify(
    query_seq: str,
    reference_db: dict[str, str],
    reference_taxonomy: dict[str, list[tuple[str, str]]],
    top_n: int = 5,
) -> TaxonomyAssignment:
    """Classify using alignment-based (BLAST-like) approach.

    Finds the top N most similar reference sequences using k-mer overlap,
    then assigns taxonomy based on the consensus of top hits.

    Args:
        query_seq: Query sequence.
        reference_db: Reference sequences.
        reference_taxonomy: Taxonomy per reference.
        top_n: Number of top hits to consider.

    Returns:
        TaxonomyAssignment with consensus lineage.
    """
    query_upper = query_seq.upper()
    k = 11  # Longer k-mers for more specific matching

    query_kmers: set[str] = set()
    for i in range(len(query_upper) - k + 1):
        kmer = query_upper[i : i + k]
        if all(c in "ACGT" for c in kmer):
            query_kmers.add(kmer)

    if not query_kmers:
        return TaxonomyAssignment(
            sequence_id="",
            lineage=[(rank, "unclassified") for rank in TAXONOMIC_RANKS],
            confidence={rank: 0.0 for rank in TAXONOMIC_RANKS},
            method="blast",
        )

    # Score each reference by k-mer overlap
    ref_scores: list[tuple[str, float]] = []
    for ref_id, ref_seq in reference_db.items():
        ref_upper = ref_seq.upper()
        ref_kmers: set[str] = set()
        for i in range(len(ref_upper) - k + 1):
            kmer = ref_upper[i : i + k]
            if all(c in "ACGT" for c in kmer):
                ref_kmers.add(kmer)
        if not ref_kmers:
            continue
        shared = len(query_kmers & ref_kmers)
        total = len(query_kmers | ref_kmers)
        jaccard = shared / total if total > 0 else 0.0
        ref_scores.append((ref_id, jaccard))

    ref_scores.sort(key=lambda x: -x[1])
    top_hits = ref_scores[:top_n]

    if not top_hits:
        return TaxonomyAssignment(
            sequence_id="",
            lineage=[(rank, "unclassified") for rank in TAXONOMIC_RANKS],
            confidence={rank: 0.0 for rank in TAXONOMIC_RANKS},
            method="blast",
        )

    best_hit_id = top_hits[0][0]
    best_hit_identity = top_hits[0][1]

    # Consensus taxonomy from top hits
    assigned_lineage: list[tuple[str, str]] = []
    confidence_scores: dict[str, float] = {}

    for rank_idx, rank in enumerate(TAXONOMIC_RANKS):
        taxon_votes: Counter[str] = Counter()
        total_weight = 0.0
        for ref_id, score in top_hits:
            if ref_id in reference_taxonomy:
                lineage = reference_taxonomy[ref_id]
                lineage_dict = {r: n for r, n in lineage}
                taxon_name = lineage_dict.get(rank, "unclassified")
                taxon_votes[taxon_name] += score  # Weight by similarity
                total_weight += score

        if taxon_votes and total_weight > 0:
            best_taxon = taxon_votes.most_common(1)[0][0]
            confidence = taxon_votes[best_taxon] / total_weight
        else:
            best_taxon = "unclassified"
            confidence = 0.0

        assigned_lineage.append((rank, best_taxon))
        confidence_scores[rank] = confidence

    return TaxonomyAssignment(
        sequence_id="",
        lineage=assigned_lineage,
        confidence=confidence_scores,
        method="blast",
        best_hit_id=best_hit_id,
        best_hit_identity=best_hit_identity,
    )


def classify_taxonomy(
    sequences: dict[str, str],
    reference_db: dict[str, str],
    reference_taxonomy: dict[str, list[tuple[str, str]]] | None = None,
    method: str = "naive_bayes",
    confidence_threshold: float = 0.8,
    k: int = 8,
    bootstrap_n: int = 100,
) -> list[TaxonomyAssignment]:
    """Classify sequences taxonomically against a reference database.

    Supports two classification methods:
    - naive_bayes: RDP-style k-mer based Bayesian classifier with bootstrap confidence.
    - blast: Alignment-based classification using k-mer overlap and consensus of top hits.

    Args:
        sequences: Dictionary mapping sequence IDs to nucleotide sequences.
        reference_db: Reference database mapping reference IDs to sequences.
        reference_taxonomy: Taxonomy for each reference sequence, as a dict
            mapping reference IDs to lists of (rank, name) tuples.
            If None, a minimal taxonomy is inferred from reference IDs.
        method: Classification method - "naive_bayes" (default) or "blast".
        confidence_threshold: Minimum confidence for reporting classification.
            Below this threshold, the taxon is reported as "unclassified".
        k: K-mer size for naive Bayes classifier (default 8).
        bootstrap_n: Number of bootstrap replicates for confidence (default 100).

    Returns:
        List of TaxonomyAssignment objects, one per input sequence.

    Raises:
        ValueError: If sequences or reference_db is empty, or method is unknown.

    Examples:
        >>> query = {"q1": "ATCGATCGATCG"}
        >>> ref_db = {"r1": "ATCGATCGATCG", "r2": "TTTTTTTTTTTT"}
        >>> ref_tax = {
        ...     "r1": [("domain", "Bacteria"), ("phylum", "Firmicutes")],
        ...     "r2": [("domain", "Bacteria"), ("phylum", "Proteobacteria")],
        ... }
        >>> results = classify_taxonomy(query, ref_db, ref_tax)
        >>> len(results)
        1
    """
    if not sequences:
        raise ValueError("Input sequences dictionary must not be empty")
    if not reference_db:
        raise ValueError("Reference database must not be empty")
    if method not in ("naive_bayes", "blast"):
        raise ValueError(f"Unknown classification method: {method}. Use 'naive_bayes' or 'blast'.")

    logger.info(f"Classifying {len(sequences)} sequences using {method} against {len(reference_db)} references")

    # Build default taxonomy if not provided
    if reference_taxonomy is None:
        reference_taxonomy = {}
        for ref_id in reference_db:
            # Attempt to parse taxonomy from ID (e.g., "Bacteria;Firmicutes;...")
            parts = ref_id.split(";")
            lineage: list[tuple[str, str]] = []
            for i, part in enumerate(parts):
                if i < len(TAXONOMIC_RANKS):
                    lineage.append((TAXONOMIC_RANKS[i], part.strip()))
            if not lineage:
                lineage = [("domain", "unclassified")]
            reference_taxonomy[ref_id] = lineage

    assignments: list[TaxonomyAssignment] = []

    for seq_id, seq in sequences.items():
        if method == "naive_bayes":
            assignment = _naive_bayes_classify(
                query_seq=seq,
                reference_db=reference_db,
                reference_taxonomy=reference_taxonomy,
                k=k,
                bootstrap_n=bootstrap_n,
                bootstrap_threshold=confidence_threshold,
            )
        else:
            assignment = _blast_classify(
                query_seq=seq,
                reference_db=reference_db,
                reference_taxonomy=reference_taxonomy,
            )

        assignment.sequence_id = seq_id

        # Apply confidence threshold: mask low-confidence assignments
        filtered_lineage: list[tuple[str, str]] = []
        for rank, name in assignment.lineage:
            if assignment.confidence.get(rank, 0.0) >= confidence_threshold:
                filtered_lineage.append((rank, name))
            else:
                filtered_lineage.append((rank, "unclassified"))
        assignment.lineage = filtered_lineage

        assignments.append(assignment)

    logger.info(f"Classification complete: {len(assignments)} sequences classified")
    return assignments


def build_taxonomy_tree(classifications: list[TaxonomyAssignment]) -> TaxonomyNode:
    """Build a hierarchical taxonomy tree from classification results.

    Constructs a tree where each node represents a taxon at a specific rank,
    with count tracking for the number of sequences assigned.

    Args:
        classifications: List of TaxonomyAssignment objects.

    Returns:
        Root TaxonomyNode of the hierarchical tree.

    Examples:
        >>> from metainformant.metagenomics.amplicon.taxonomy import TaxonomyAssignment
        >>> assignments = [
        ...     TaxonomyAssignment("s1", [("domain", "Bacteria"), ("phylum", "Firmicutes")]),
        ...     TaxonomyAssignment("s2", [("domain", "Bacteria"), ("phylum", "Proteobacteria")]),
        ... ]
        >>> tree = build_taxonomy_tree(assignments)
        >>> tree.total_descendants()
        2
    """
    root = TaxonomyNode(name="root", rank="root")

    for assignment in classifications:
        current_node = root
        for rank, name in assignment.lineage:
            if name == "unclassified":
                # Create unclassified branch but don't traverse further
                child = current_node.add_child(f"unclassified_{rank}", rank)
                child.count += 1
                break
            child = current_node.add_child(name, rank)
            current_node = child
        # Increment the leaf count
        current_node.count += 1

    logger.info(
        f"Taxonomy tree built: {root.total_descendants()} sequences, "
        f"{len(root.children)} top-level taxa"
    )
    return root


def calculate_confidence(
    assignments: list[TaxonomyAssignment],
    min_confidence: float = 0.0,
) -> dict[str, dict[str, float]]:
    """Calculate and summarize confidence metrics for taxonomy assignments.

    Computes per-rank confidence statistics across all assignments:
    mean confidence, median confidence, fraction above threshold.

    Args:
        assignments: List of TaxonomyAssignment objects with confidence scores.
        min_confidence: Minimum confidence threshold for "confident" assignments.

    Returns:
        Dictionary mapping rank names to confidence statistics:
        {rank: {"mean": float, "median": float, "fraction_confident": float, "n": int}}

    Examples:
        >>> from metainformant.metagenomics.amplicon.taxonomy import TaxonomyAssignment
        >>> assignments = [
        ...     TaxonomyAssignment("s1", [], confidence={"domain": 0.95, "phylum": 0.80}),
        ...     TaxonomyAssignment("s2", [], confidence={"domain": 0.90, "phylum": 0.60}),
        ... ]
        >>> stats = calculate_confidence(assignments, min_confidence=0.8)
        >>> stats["domain"]["mean"]
        0.925
    """
    if not assignments:
        return {}

    # Collect per-rank confidence values
    rank_confidences: dict[str, list[float]] = defaultdict(list)
    for assignment in assignments:
        for rank, conf in assignment.confidence.items():
            rank_confidences[rank].append(conf)

    # Compute statistics
    result: dict[str, dict[str, float]] = {}
    for rank, values in rank_confidences.items():
        n = len(values)
        mean_conf = sum(values) / n if n > 0 else 0.0
        sorted_vals = sorted(values)
        if n % 2 == 0 and n > 0:
            median_conf = (sorted_vals[n // 2 - 1] + sorted_vals[n // 2]) / 2.0
        elif n > 0:
            median_conf = sorted_vals[n // 2]
        else:
            median_conf = 0.0
        fraction_confident = sum(1 for v in values if v >= min_confidence) / n if n > 0 else 0.0

        result[rank] = {
            "mean": mean_conf,
            "median": median_conf,
            "fraction_confident": fraction_confident,
            "n": float(n),
        }

    return result


__all__ = [
    "TAXONOMIC_RANKS",
    "TaxonomyAssignment",
    "TaxonomyNode",
    "build_taxonomy_tree",
    "calculate_confidence",
    "classify_taxonomy",
]
