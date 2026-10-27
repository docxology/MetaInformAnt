"""Taxonomic classification for amplicon sequences.

Implements naive Bayes and alignment-based taxonomic classification methods
for assigning taxonomy to OTU/ASV representative sequences. Includes
hierarchical taxonomy tree construction and bootstrap confidence estimation.

Taxonomic ranks follow the standard hierarchy:
    Domain > Phylum > Class > Order > Family > Genus > Species
"""

from __future__ import annotations

import math
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


def _nb_winner(
    kmers: list[str],
    log_prior: dict[str, float],
    kmer_logp: dict[str, dict[str, float]],
) -> str | None:
    """Return the best taxon for a k-mer sample, or None when it cannot be decided.

    Args:
        kmers: K-mers of the query (or of a bootstrap sample of it).
        log_prior: Log prior of each taxon.
        kmer_logp: Log-probability of each informative k-mer under each taxon. K-mers
            missing from this mapping are absent from every reference and are ignored.

    Returns:
        The taxon with the highest log posterior, or None if the sample has no
        informative k-mer or the two best taxa tie.
    """
    known = [kmer for kmer in kmers if kmer in kmer_logp]
    if not known:
        return None
    scores = dict(log_prior)
    for kmer in known:
        for taxon, logp in kmer_logp[kmer].items():
            scores[taxon] += logp
    ranked = sorted(scores.items(), key=lambda item: item[1], reverse=True)
    if len(ranked) > 1 and ranked[0][1] - ranked[1][1] <= 1e-12 * max(1.0, abs(ranked[0][1])):
        return None
    return ranked[0][0]


def _naive_bayes_classify(
    query_seq: str,
    reference_db: dict[str, str],
    reference_taxonomy: dict[str, list[tuple[str, str]]],
    k: int = 8,
    bootstrap_n: int = 100,
    confidence_threshold: float = 0.8,
    seed: int | None = 42,
) -> TaxonomyAssignment:
    """Classify a query sequence using naive Bayes k-mer approach.

    The RDP-style naive Bayes classifier:
    1. Extract all k-mers from the query sequence.
    2. For each reference taxon at each rank, compute the posterior probability
       based on k-mer frequencies. K-mers absent from every reference carry no
       taxonomic information and are ignored; a query with no known k-mer, or one
       whose best two taxa tie, cannot be classified.
    3. Assign the query to the taxon with highest posterior at each rank.
    4. Bootstrap: draw ``ceil(n_kmers / k)`` k-mers with replacement (about the
       number of non-overlapping words, as in the RDP classifier) ``bootstrap_n``
       times and reclassify; the confidence is the fraction of trials that agree
       with the full-query assignment. Trials without a known k-mer, or with a tie,
       count as disagreeing.
    5. Abstain: a rank is reported as "unclassified" when its bootstrap confidence is
       below ``confidence_threshold``, and so is every rank below it (as in RDP, the
       lineage is cut at the first rank that is not confidently assigned).

    Args:
        query_seq: Query nucleotide sequence.
        reference_db: Reference sequences keyed by ID.
        reference_taxonomy: Taxonomy for each reference sequence.
        k: K-mer size for classification.
        bootstrap_n: Number of bootstrap iterations for confidence.
        confidence_threshold: Minimum bootstrap confidence to report a rank (default 0.8,
            the RDP default for full-length sequences).
        seed: Seed of the bootstrap random generator (None draws fresh entropy).

    Returns:
        TaxonomyAssignment with lineage and per-rank bootstrap confidence scores
        (the confidence of ranks below the first abstention is still reported).
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

    rng = random.Random(seed)
    subsample_size = max(1, math.ceil(len(query_kmers) / k))
    distinct_query_kmers = set(query_kmers)

    # Classify at each rank using log-likelihood
    assigned_lineage: list[tuple[str, str]] = []
    confidence_scores: dict[str, float] = {}
    abstained = False

    for rank in TAXONOMIC_RANKS:
        taxon_models = rank_taxon_kmers[rank]
        taxon_totals = {taxon: sum(counts.values()) for taxon, counts in taxon_models.items()}
        total_refs = sum(taxon_totals.values())
        if not taxon_models or total_refs == 0:
            assigned_lineage.append((rank, "unclassified"))
            confidence_scores[rank] = 0.0
            abstained = True
            continue

        all_kmers_vocab: set[str] = set()
        for kmer_counts in taxon_models.values():
            all_kmers_vocab.update(kmer_counts.keys())
        vocab_size = len(all_kmers_vocab)

        # Prior: fraction of reference k-mers in each taxon
        log_prior = {taxon: math.log(total / total_refs) for taxon, total in taxon_totals.items()}

        # Laplace-smoothed log-probability of each informative query k-mer under each taxon.
        # K-mers absent from every taxon are uninformative and are left out.
        kmer_logp: dict[str, dict[str, float]] = {}
        for qk in distinct_query_kmers & all_kmers_vocab:
            kmer_logp[qk] = {
                taxon: math.log((kmer_counts.get(qk, 0) + 0.5) / (taxon_totals[taxon] + 0.5 * vocab_size))
                for taxon, kmer_counts in taxon_models.items()
            }

        best_taxon = _nb_winner(query_kmers, log_prior, kmer_logp)

        # Bootstrap confidence: resample k-mers and reclassify
        if best_taxon is None or bootstrap_n <= 0:
            confidence = 0.0
        else:
            bootstrap_hits = 0
            for _ in range(bootstrap_n):
                resampled = rng.choices(query_kmers, k=subsample_size)
                if _nb_winner(resampled, log_prior, kmer_logp) == best_taxon:
                    bootstrap_hits += 1
            confidence = bootstrap_hits / bootstrap_n

        confidence_scores[rank] = confidence
        if abstained or best_taxon is None or best_taxon == "unclassified" or confidence < confidence_threshold:
            abstained = True
            assigned_lineage.append((rank, "unclassified"))
        else:
            assigned_lineage.append((rank, best_taxon))

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
        taxon_votes: dict[str, float] = {}
        total_weight = 0.0
        for ref_id, score in top_hits:
            if ref_id in reference_taxonomy:
                lineage = reference_taxonomy[ref_id]
                lineage_dict = {r: n for r, n in lineage}
                taxon_name = lineage_dict.get(rank, "unclassified")
                taxon_votes[taxon_name] = taxon_votes.get(taxon_name, 0.0) + score
                total_weight += score

        if taxon_votes and total_weight > 0:
            best_taxon = max(taxon_votes, key=lambda name: taxon_votes[name])
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
    seed: int | None = 42,
) -> list[TaxonomyAssignment]:
    """Classify sequences taxonomically against a reference database.

    Supports two classification methods:
    - naive_bayes: RDP-style k-mer based Bayesian classifier. It abstains (reports
      "unclassified") at every rank whose bootstrap confidence is below
      ``confidence_threshold`` and at all ranks below it.
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
        seed: Seed of the naive Bayes bootstrap (default 42, so results are
            reproducible; None draws fresh entropy). Every query uses its own
            generator started from this seed, so a result does not depend on the
            other queries in the batch.

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
                confidence_threshold=confidence_threshold,
                seed=seed,
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
                # Create unclassified branch but don't traverse further; the sequence is
                # counted there and nowhere else
                current_node = current_node.add_child(f"unclassified_{rank}", rank)
                break
            current_node = current_node.add_child(name, rank)
        # Count the sequence once, at the deepest node it reached
        current_node.count += 1

    logger.info(f"Taxonomy tree built: {root.total_descendants()} sequences, " f"{len(root.children)} top-level taxa")
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
            "n": n,
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
