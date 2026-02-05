"""Community profiling from shotgun metagenomic reads.

Implements taxonomic profiling using k-mer-based classification (Kraken-style),
relative abundance calculation, and k-mer index construction. The approach
classifies reads by matching k-mers against a reference database and using
a lowest common ancestor (LCA) algorithm for taxonomic assignment.

K-mer profiling approach:
1. Build a k-mer index from reference genomes, associating each k-mer
   with its taxonomic lineage.
2. For each read, extract k-mers and look up their taxonomy.
3. Use LCA to assign the read to the most specific taxon supported by
   the majority of its k-mers.
4. Aggregate assignments into community profiles with relative abundances.
"""

from __future__ import annotations

import os
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

_ENV_PREFIX = "META_"


@dataclass
class TaxonProfile:
    """Profile for a single taxon in the community.

    Attributes:
        taxon_id: Unique taxon identifier.
        name: Taxon name.
        rank: Taxonomic rank.
        lineage: Full taxonomic lineage as list of (rank, name) tuples.
        read_count: Number of reads assigned.
        relative_abundance: Fraction of total classified reads.
        kmer_count: Number of unique k-mers matched.
    """

    taxon_id: str
    name: str
    rank: str
    lineage: list[tuple[str, str]] = field(default_factory=list)
    read_count: int = 0
    relative_abundance: float = 0.0
    kmer_count: int = 0


@dataclass
class CommunityProfile:
    """Full community taxonomic profile.

    Attributes:
        taxa: List of taxon profiles.
        total_reads: Total number of reads profiled.
        classified_reads: Number of reads successfully classified.
        unclassified_reads: Number of reads that could not be classified.
        classification_rate: Fraction of reads classified.
        rank_summaries: Per-rank abundance summaries.
    """

    taxa: list[TaxonProfile]
    total_reads: int = 0
    classified_reads: int = 0
    unclassified_reads: int = 0
    classification_rate: float = 0.0
    rank_summaries: dict[str, dict[str, float]] = field(default_factory=dict)


@dataclass
class KmerIndex:
    """K-mer index mapping k-mers to taxonomic lineages.

    Attributes:
        kmer_to_taxon: Mapping of k-mer strings to taxon IDs.
        taxon_lineages: Mapping of taxon IDs to full lineages.
        k: K-mer size used in the index.
        total_kmers: Total number of unique k-mers indexed.
        num_taxa: Number of taxa represented.
    """

    kmer_to_taxon: dict[str, str]
    taxon_lineages: dict[str, list[tuple[str, str]]]
    k: int
    total_kmers: int = 0
    num_taxa: int = 0


def build_kmer_index(
    reference_sequences: dict[str, str],
    taxonomy: dict[str, list[tuple[str, str]]] | None = None,
    k: int = 31,
) -> KmerIndex:
    """Build a k-mer index from reference sequences for taxonomic profiling.

    For each reference sequence, extracts all k-mers and associates them
    with the sequence's taxonomic lineage. K-mers shared between multiple
    taxa are assigned to their lowest common ancestor (LCA).

    Args:
        reference_sequences: Dictionary mapping reference IDs to DNA sequences.
        taxonomy: Taxonomy for each reference. Maps reference IDs to lists
            of (rank, name) tuples. If None, each reference ID is treated
            as its own taxon.
        k: K-mer size (default 31). Larger values give higher specificity
            but require more memory and may miss divergent sequences.

    Returns:
        KmerIndex ready for use in profiling.

    Raises:
        ValueError: If reference_sequences is empty or k < 11.

    Examples:
        >>> refs = {"ref1": "ATCGATCGATCGATCGATCGATCGATCGATCG"}
        >>> tax = {"ref1": [("domain", "Bacteria"), ("phylum", "Firmicutes")]}
        >>> index = build_kmer_index(refs, tax, k=11)
        >>> index.total_kmers > 0
        True
    """
    if not reference_sequences:
        raise ValueError("Reference sequences must not be empty")
    if k < 11:
        raise ValueError(f"K-mer size must be >= 11, got {k}")

    logger.info(f"Building k-mer index (k={k}) from {len(reference_sequences)} references")

    # Default taxonomy if not provided
    if taxonomy is None:
        taxonomy = {ref_id: [("species", ref_id)] for ref_id in reference_sequences}

    # Extract k-mers and associate with taxa
    kmer_taxa: dict[str, set[str]] = defaultdict(set)  # kmer -> set of taxon_ids

    for ref_id, seq in reference_sequences.items():
        seq_upper = seq.upper()
        for i in range(len(seq_upper) - k + 1):
            kmer = seq_upper[i : i + k]
            if all(c in "ACGT" for c in kmer):
                kmer_taxa[kmer].add(ref_id)

    # Resolve shared k-mers using LCA
    kmer_to_taxon: dict[str, str] = {}
    taxon_lineages: dict[str, list[tuple[str, str]]] = {}

    for kmer, taxa_ids in kmer_taxa.items():
        if len(taxa_ids) == 1:
            taxon_id = next(iter(taxa_ids))
            kmer_to_taxon[kmer] = taxon_id
            if taxon_id not in taxon_lineages:
                taxon_lineages[taxon_id] = taxonomy.get(taxon_id, [("species", taxon_id)])
        else:
            # Compute LCA
            lca_lineage = _compute_lca(
                [taxonomy.get(tid, [("species", tid)]) for tid in taxa_ids]
            )
            lca_id = _lineage_to_id(lca_lineage)
            kmer_to_taxon[kmer] = lca_id
            if lca_id not in taxon_lineages:
                taxon_lineages[lca_id] = lca_lineage

    # Ensure all reference taxa are in lineages
    for ref_id in reference_sequences:
        if ref_id not in taxon_lineages:
            taxon_lineages[ref_id] = taxonomy.get(ref_id, [("species", ref_id)])

    index = KmerIndex(
        kmer_to_taxon=kmer_to_taxon,
        taxon_lineages=taxon_lineages,
        k=k,
        total_kmers=len(kmer_to_taxon),
        num_taxa=len(taxon_lineages),
    )

    logger.info(f"K-mer index built: {index.total_kmers} k-mers, {index.num_taxa} taxa")
    return index


def _compute_lca(
    lineages: list[list[tuple[str, str]]],
) -> list[tuple[str, str]]:
    """Compute the lowest common ancestor of multiple lineages.

    Walks down the lineage from domain to species, stopping at the
    first rank where the lineages disagree.

    Args:
        lineages: List of taxonomic lineages.

    Returns:
        LCA lineage (the shared prefix of all input lineages).
    """
    if not lineages:
        return []
    if len(lineages) == 1:
        return lineages[0]

    # Convert lineages to dicts keyed by rank
    rank_order = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    lineage_dicts = []
    for lin in lineages:
        ld: dict[str, str] = {}
        for rank, name in lin:
            ld[rank] = name
        lineage_dicts.append(ld)

    lca: list[tuple[str, str]] = []
    for rank in rank_order:
        names = set()
        for ld in lineage_dicts:
            if rank in ld:
                names.add(ld[rank])
            else:
                names.add(None)  # type: ignore[arg-type]
        if len(names) == 1 and None not in names:
            lca.append((rank, names.pop()))
        else:
            break

    return lca


def _lineage_to_id(lineage: list[tuple[str, str]]) -> str:
    """Convert a lineage to a string identifier.

    Args:
        lineage: List of (rank, name) tuples.

    Returns:
        String identifier like "Bacteria;Firmicutes;Bacilli".
    """
    return ";".join(name for _, name in lineage) if lineage else "unclassified"


def profile_community(
    reads: dict[str, str] | list[str],
    database: KmerIndex | None = None,
    reference_sequences: dict[str, str] | None = None,
    reference_taxonomy: dict[str, list[tuple[str, str]]] | None = None,
    k: int = 31,
    min_kmer_hits: int = 2,
    confidence_threshold: float = 0.5,
) -> CommunityProfile:
    """Profile the taxonomic composition of a metagenomic sample.

    Classifies reads by matching k-mers against a reference index and
    assigning taxonomy via LCA-based consensus. Each read is classified
    to the most specific taxon supported by >= confidence_threshold of
    its matched k-mers.

    Args:
        reads: Sequencing reads as dict (id->sequence) or list of sequences.
        database: Pre-built KmerIndex. If None, one is built from
            reference_sequences and reference_taxonomy.
        reference_sequences: Reference sequences for building index.
            Required if database is None.
        reference_taxonomy: Taxonomy for references. Required if database is None.
        k: K-mer size (default 31). Ignored if database is provided.
        min_kmer_hits: Minimum k-mer matches required to classify a read.
        confidence_threshold: Fraction of k-mers that must agree on
            taxonomy for classification (default 0.5).

    Returns:
        CommunityProfile with per-taxon abundance data.

    Raises:
        ValueError: If no reads provided, or no database/references available.

    Examples:
        >>> reads = ["ATCGATCGATCGATCGATCGATCGATCGATCG"]
        >>> refs = {"ref1": "ATCGATCGATCGATCGATCGATCGATCGATCG"}
        >>> tax = {"ref1": [("domain", "Bacteria"), ("phylum", "Firmicutes")]}
        >>> profile = profile_community(reads, reference_sequences=refs,
        ...                             reference_taxonomy=tax, k=11)
        >>> profile.total_reads
        1
    """
    if not reads:
        raise ValueError("Input reads must not be empty")

    # Build index if not provided
    if database is None:
        if reference_sequences is None:
            raise ValueError("Either database or reference_sequences must be provided")
        database = build_kmer_index(reference_sequences, reference_taxonomy, k=k)

    k = database.k

    # Normalize reads
    if isinstance(reads, dict):
        read_list = list(reads.values())
    else:
        read_list = list(reads)

    logger.info(f"Profiling {len(read_list)} reads against index with {database.total_kmers} k-mers")

    # Classify each read
    taxon_counts: Counter[str] = Counter()
    classified = 0
    unclassified = 0

    for read in read_list:
        read_upper = read.upper()
        # Extract k-mers from read
        kmer_hits: Counter[str] = Counter()
        total_kmers_in_read = 0

        for i in range(len(read_upper) - k + 1):
            kmer = read_upper[i : i + k]
            if all(c in "ACGT" for c in kmer):
                total_kmers_in_read += 1
                if kmer in database.kmer_to_taxon:
                    taxon_id = database.kmer_to_taxon[kmer]
                    kmer_hits[taxon_id] += 1

        total_hits = sum(kmer_hits.values())

        if total_hits < min_kmer_hits:
            unclassified += 1
            continue

        # Find the most specific classification with sufficient confidence
        # First, aggregate hits at each taxonomic level
        if kmer_hits:
            best_taxon = kmer_hits.most_common(1)[0][0]
            best_count = kmer_hits[best_taxon]

            # Check confidence: fraction of hits supporting this taxon
            confidence = best_count / total_hits if total_hits > 0 else 0.0

            if confidence >= confidence_threshold:
                taxon_counts[best_taxon] += 1
                classified += 1
            else:
                # Try LCA of all hit taxa
                hit_lineages = []
                for tid in kmer_hits:
                    if tid in database.taxon_lineages:
                        hit_lineages.append(database.taxon_lineages[tid])
                if hit_lineages:
                    lca = _compute_lca(hit_lineages)
                    lca_id = _lineage_to_id(lca)
                    taxon_counts[lca_id] += 1
                    classified += 1
                else:
                    unclassified += 1
        else:
            unclassified += 1

    # Build community profile
    total = len(read_list)
    taxa_profiles: list[TaxonProfile] = []

    for taxon_id, count in taxon_counts.most_common():
        lineage = database.taxon_lineages.get(taxon_id, [])
        if not lineage:
            # Parse from taxon_id
            parts = taxon_id.split(";")
            rank_order = ["domain", "phylum", "class", "order", "family", "genus", "species"]
            lineage = [(rank_order[i] if i < len(rank_order) else "unknown", p) for i, p in enumerate(parts)]

        name = lineage[-1][1] if lineage else taxon_id
        rank = lineage[-1][0] if lineage else "unknown"

        profile = TaxonProfile(
            taxon_id=taxon_id,
            name=name,
            rank=rank,
            lineage=lineage,
            read_count=count,
            relative_abundance=count / classified if classified > 0 else 0.0,
        )
        taxa_profiles.append(profile)

    # Build rank summaries
    rank_summaries: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    for tp in taxa_profiles:
        for rank, name in tp.lineage:
            rank_summaries[rank][name] += tp.relative_abundance

    community = CommunityProfile(
        taxa=taxa_profiles,
        total_reads=total,
        classified_reads=classified,
        unclassified_reads=unclassified,
        classification_rate=classified / total if total > 0 else 0.0,
        rank_summaries=dict(rank_summaries),
    )

    logger.info(
        f"Profiling complete: {classified}/{total} reads classified "
        f"({community.classification_rate:.1%}), {len(taxa_profiles)} taxa detected"
    )
    return community


def calculate_relative_abundance(
    profile: CommunityProfile,
    rank: str | None = None,
    min_abundance: float = 0.0,
) -> dict[str, float]:
    """Calculate relative abundance from a community profile.

    Normalizes taxon counts to relative abundances summing to 1.0,
    optionally filtered to a specific taxonomic rank.

    Args:
        profile: CommunityProfile to process.
        rank: Taxonomic rank to aggregate at (e.g., "phylum", "genus").
            If None, uses the most specific rank for each taxon.
        min_abundance: Minimum relative abundance to include (default 0.0).
            Taxa below this threshold are grouped as "Other".

    Returns:
        Dictionary mapping taxon names to relative abundances (sum = 1.0).

    Examples:
        >>> # Given a CommunityProfile with taxa
        >>> profile = CommunityProfile(
        ...     taxa=[TaxonProfile("t1", "Firmicutes", "phylum", read_count=70),
        ...           TaxonProfile("t2", "Proteobacteria", "phylum", read_count=30)],
        ...     classified_reads=100
        ... )
        >>> abundances = calculate_relative_abundance(profile)
        >>> abs(sum(abundances.values()) - 1.0) < 0.01
        True
    """
    if not profile.taxa:
        return {}

    # Aggregate by rank if specified
    if rank:
        rank_counts: Counter[str] = Counter()
        for taxon in profile.taxa:
            taxon_name = "unclassified"
            for r, name in taxon.lineage:
                if r == rank:
                    taxon_name = name
                    break
            rank_counts[taxon_name] += taxon.read_count

        total = sum(rank_counts.values())
        if total == 0:
            return {}
        abundances = {name: count / total for name, count in rank_counts.items()}
    else:
        total = sum(t.read_count for t in profile.taxa)
        if total == 0:
            return {}
        abundances = {t.name: t.read_count / total for t in profile.taxa}

    # Apply minimum abundance filter
    if min_abundance > 0:
        other_abundance = 0.0
        filtered: dict[str, float] = {}
        for name, abund in abundances.items():
            if abund >= min_abundance:
                filtered[name] = abund
            else:
                other_abundance += abund
        if other_abundance > 0:
            filtered["Other"] = other_abundance
        abundances = filtered

    # Renormalize to ensure sum = 1.0
    total_abund = sum(abundances.values())
    if total_abund > 0:
        abundances = {name: a / total_abund for name, a in abundances.items()}

    return abundances


__all__ = [
    "CommunityProfile",
    "KmerIndex",
    "TaxonProfile",
    "build_kmer_index",
    "calculate_relative_abundance",
    "profile_community",
]
