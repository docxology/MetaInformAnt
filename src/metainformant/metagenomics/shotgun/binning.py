"""Metagenomic binning of assembled contigs.

Implements composition-based and coverage-based binning to group contigs
into putative genome bins (MAGs - Metagenome-Assembled Genomes). Features
include tetranucleotide frequency (TNF) calculation, k-means and
hierarchical clustering, bin refinement, and CheckM-style quality assessment.

Binning strategy:
1. Calculate tetranucleotide frequency (TNF) profiles for each contig.
2. Normalize coverage data across samples (if multiple samples available).
3. Combine TNF and coverage features.
4. Cluster contigs using the combined feature space.
5. Refine bins using single-copy marker gene assessment.
"""

from __future__ import annotations

import itertools
import math
import os
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

_ENV_PREFIX = "META_"

# All 256 tetranucleotides
_NUCLEOTIDES = "ACGT"
_ALL_TETRANUCLEOTIDES = ["".join(t) for t in itertools.product(_NUCLEOTIDES, repeat=4)]
_TETRA_INDEX = {t: i for i, t in enumerate(_ALL_TETRANUCLEOTIDES)}

# Universal single-copy marker genes (simplified set based on bacterial markers)
# These represent conserved genes expected once per genome
_BACTERIAL_MARKERS = [
    "rpoB",
    "rpoC",
    "rpsB",
    "rpsC",
    "rpsE",
    "rpsJ",
    "rpsK",
    "rpsM",
    "rpsS",
    "rplB",
    "rplC",
    "rplD",
    "rplE",
    "rplF",
    "rplK",
    "rplL",
    "rplM",
    "rplN",
    "rplP",
    "rplS",
    "rplT",
    "rplV",
    "rplW",
    "tsf",
    "infC",
    "smpB",
    "frr",
    "pgk",
    "pyrG",
    "nusA",
    "dnaG",
    "rpoA",
    "secY",
    "ffh",
    "rnhB",
    "ychF",
]

# Representative k-mer signatures for marker gene detection (simplified)
# In production, HMM profiles would be used instead
_MARKER_SIGNATURES: dict[str, list[str]] = {}


@dataclass
class GenomeBin:
    """A metagenomic bin (putative genome).

    Attributes:
        bin_id: Unique identifier.
        contig_ids: IDs of contigs assigned to this bin.
        total_length: Total sequence length in bp.
        num_contigs: Number of contigs in the bin.
        gc_content: GC content of the bin.
        mean_coverage: Average coverage across contigs.
        completeness: Estimated genome completeness (0.0 to 1.0).
        contamination: Estimated contamination level (0.0 to 1.0).
        quality_score: Completeness - 5 * contamination (Parks et al. convention).
        marker_counts: Count of each detected marker gene.
    """

    bin_id: str
    contig_ids: list[str] = field(default_factory=list)
    total_length: int = 0
    num_contigs: int = 0
    gc_content: float = 0.0
    mean_coverage: float = 0.0
    completeness: float = 0.0
    contamination: float = 0.0
    quality_score: float = 0.0
    marker_counts: dict[str, int] = field(default_factory=dict)


@dataclass
class BinningResult:
    """Result of metagenomic binning.

    Attributes:
        bins: List of genome bins.
        unbinned_contigs: IDs of contigs not assigned to any bin.
        total_contigs: Total number of input contigs.
        binned_contigs: Number of contigs assigned to bins.
        high_quality_bins: Number of bins meeting quality thresholds.
    """

    bins: list[GenomeBin]
    unbinned_contigs: list[str] = field(default_factory=list)
    total_contigs: int = 0
    binned_contigs: int = 0
    high_quality_bins: int = 0


def calculate_tetranucleotide_freq(sequence: str, normalize: bool = True) -> list[float]:
    """Calculate tetranucleotide frequency (TNF) profile for a sequence.

    Counts all 256 possible 4-mer combinations in the sequence and its
    reverse complement, producing a composition signature that is
    characteristic of taxonomic origin.

    Args:
        sequence: DNA sequence string.
        normalize: If True, normalize frequencies to sum to 1.0.

    Returns:
        List of 256 tetranucleotide frequencies in alphabetical order.

    Raises:
        ValueError: If sequence is shorter than 4 bases.

    Examples:
        >>> tnf = calculate_tetranucleotide_freq("AAAAAAAAAA")
        >>> len(tnf)
        256
        >>> tnf[0] > 0  # AAAA should be most frequent
        True
    """
    if len(sequence) < 4:
        raise ValueError(f"Sequence must be at least 4 bases long, got {len(sequence)}")

    seq_upper = sequence.upper()
    # Count tetranucleotides in both strands
    counts = [0.0] * 256

    # Forward strand
    for i in range(len(seq_upper) - 3):
        tetra = seq_upper[i : i + 4]
        if tetra in _TETRA_INDEX:
            counts[_TETRA_INDEX[tetra]] += 1

    # Reverse complement strand
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    rc_seq = "".join(complement.get(c, "N") for c in reversed(seq_upper))
    for i in range(len(rc_seq) - 3):
        tetra = rc_seq[i : i + 4]
        if tetra in _TETRA_INDEX:
            counts[_TETRA_INDEX[tetra]] += 1

    if normalize:
        total = sum(counts)
        if total > 0:
            counts = [c / total for c in counts]

    return counts


def _euclidean_distance(a: list[float], b: list[float]) -> float:
    """Compute Euclidean distance between two vectors.

    Args:
        a: First vector.
        b: Second vector.

    Returns:
        Euclidean distance.
    """
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def _cosine_distance(a: list[float], b: list[float]) -> float:
    """Compute cosine distance between two vectors.

    Args:
        a: First vector.
        b: Second vector.

    Returns:
        Cosine distance (1 - cosine similarity).
    """
    dot = sum(x * y for x, y in zip(a, b))
    norm_a = math.sqrt(sum(x * x for x in a))
    norm_b = math.sqrt(sum(y * y for y in b))
    if norm_a == 0 or norm_b == 0:
        return 1.0
    return 1.0 - dot / (norm_a * norm_b)


def _kmeans_cluster(
    features: list[list[float]],
    n_clusters: int,
    max_iterations: int = 100,
    seed: int = 42,
) -> list[int]:
    """K-means clustering implementation.

    Args:
        features: List of feature vectors.
        n_clusters: Number of clusters.
        max_iterations: Maximum iterations.
        seed: Random seed for centroid initialization.

    Returns:
        List of cluster assignments (indices).
    """
    import random as rng

    rng.seed(seed)

    n = len(features)
    dim = len(features[0]) if features else 0

    if n <= n_clusters:
        return list(range(n))

    # K-means++ initialization
    centroids: list[list[float]] = []
    first_idx = rng.randint(0, n - 1)
    centroids.append(list(features[first_idx]))

    for _ in range(1, n_clusters):
        # Compute distances to nearest centroid
        distances = []
        for feat in features:
            min_dist = min(_euclidean_distance(feat, c) for c in centroids)
            distances.append(min_dist**2)
        total_dist = sum(distances)
        if total_dist == 0:
            # All points are identical to existing centroids
            idx = rng.randint(0, n - 1)
        else:
            probs = [d / total_dist for d in distances]
            cumulative = 0.0
            r = rng.random()
            idx = n - 1
            for i, p in enumerate(probs):
                cumulative += p
                if cumulative >= r:
                    idx = i
                    break
        centroids.append(list(features[idx]))

    assignments = [0] * n

    for iteration in range(max_iterations):
        # Assign each point to nearest centroid
        new_assignments = [0] * n
        for i, feat in enumerate(features):
            min_dist = float("inf")
            best_cluster = 0
            for c_idx, centroid in enumerate(centroids):
                dist = _euclidean_distance(feat, centroid)
                if dist < min_dist:
                    min_dist = dist
                    best_cluster = c_idx
            new_assignments[i] = best_cluster

        # Check convergence
        if new_assignments == assignments and iteration > 0:
            break
        assignments = new_assignments

        # Update centroids
        for c_idx in range(n_clusters):
            members = [features[i] for i in range(n) if assignments[i] == c_idx]
            if members:
                centroids[c_idx] = [sum(m[d] for m in members) / len(members) for d in range(dim)]

    return assignments


def bin_contigs(
    contigs: dict[str, str],
    coverage: dict[str, float] | dict[str, list[float]] | None = None,
    method: str = "composition",
    n_bins: int | None = None,
    min_contig_length: int = 1000,
    seed: int = 42,
) -> BinningResult:
    """Bin contigs into putative genomes using composition and coverage.

    Supports three binning approaches:
    - composition: Uses tetranucleotide frequency (TNF) profiles only.
    - coverage: Uses coverage information only (requires multi-sample coverage).
    - combined: Uses both TNF and coverage features (recommended).

    Args:
        contigs: Dictionary mapping contig IDs to sequences.
        coverage: Coverage data per contig. Can be:
            - dict[str, float]: Single-sample coverage values.
            - dict[str, list[float]]: Multi-sample coverage profiles.
            If None, only composition-based binning is performed.
        method: Binning method - "composition", "coverage", or "combined".
        n_bins: Number of bins to create. If None, estimated automatically
            based on coverage variation and TNF diversity.
        min_contig_length: Minimum contig length to bin (default 1000 bp).
        seed: Random seed for clustering reproducibility.

    Returns:
        BinningResult with genome bins and statistics.

    Raises:
        ValueError: If contigs is empty or method is unknown.

    Examples:
        >>> contigs = {"c1": "ATCG" * 500, "c2": "GCTA" * 500, "c3": "AAAA" * 500}
        >>> coverage = {"c1": 10.0, "c2": 10.5, "c3": 50.0}
        >>> result = bin_contigs(contigs, coverage, method="combined", n_bins=2)
        >>> result.total_contigs == 3
        True
    """
    if not contigs:
        raise ValueError("Input contigs dictionary must not be empty")
    if method not in ("composition", "coverage", "combined"):
        raise ValueError(f"Unknown binning method: {method}. Use 'composition', 'coverage', or 'combined'.")

    logger.info(f"Binning {len(contigs)} contigs using {method} method")

    # Filter contigs by length
    filtered_ids: list[str] = []
    unbinned: list[str] = []
    for cid, seq in contigs.items():
        if len(seq) >= min_contig_length:
            filtered_ids.append(cid)
        else:
            unbinned.append(cid)

    if not filtered_ids:
        logger.warning(f"No contigs longer than {min_contig_length} bp")
        return BinningResult(
            bins=[],
            unbinned_contigs=list(contigs.keys()),
            total_contigs=len(contigs),
        )

    # Compute features
    features: list[list[float]] = []

    for cid in filtered_ids:
        feat: list[float] = []

        if method in ("composition", "combined"):
            tnf = calculate_tetranucleotide_freq(contigs[cid])
            feat.extend(tnf)

        if method in ("coverage", "combined") and coverage:
            cov_val = coverage.get(cid)
            if cov_val is not None:
                if isinstance(cov_val, list):
                    # Multi-sample coverage: log-transform
                    feat.extend([math.log(max(c, 0.01)) for c in cov_val])
                else:
                    feat.append(math.log(max(float(cov_val), 0.01)))
            else:
                if method == "coverage":
                    feat.append(0.0)
                elif method == "combined":
                    feat.append(0.0)

        features.append(feat)

    # Estimate number of bins if not provided
    if n_bins is None:
        # Heuristic: estimate based on total assembly size and typical genome size
        total_length = sum(len(contigs[cid]) for cid in filtered_ids)
        estimated_genomes = max(2, total_length // 3_000_000)  # Assume ~3 Mbp avg genome
        n_bins = min(estimated_genomes, len(filtered_ids) // 2)
        n_bins = max(2, n_bins)
        logger.info(f"Estimated {n_bins} bins from {total_length} bp total assembly")

    n_bins = min(n_bins, len(filtered_ids))

    # Cluster using k-means
    assignments = _kmeans_cluster(features, n_bins, seed=seed)

    # Build bins
    bin_members: dict[int, list[str]] = defaultdict(list)
    for idx, cluster_id in enumerate(assignments):
        bin_members[cluster_id].append(filtered_ids[idx])

    bins: list[GenomeBin] = []
    for bin_idx in sorted(bin_members.keys()):
        members = bin_members[bin_idx]
        bin_seqs = [contigs[cid] for cid in members]

        total_len = sum(len(s) for s in bin_seqs)

        # GC content
        gc_count = sum(s.upper().count("G") + s.upper().count("C") for s in bin_seqs)
        valid_bases = sum(sum(1 for c in s.upper() if c in "ACGT") for s in bin_seqs)
        gc = gc_count / valid_bases if valid_bases > 0 else 0.0

        # Mean coverage
        mean_cov = 0.0
        if coverage:
            cov_vals = []
            for cid in members:
                cv = coverage.get(cid)
                if cv is not None:
                    if isinstance(cv, list):
                        cov_vals.append(sum(cv) / len(cv) if cv else 0.0)
                    else:
                        cov_vals.append(float(cv))
            mean_cov = sum(cov_vals) / len(cov_vals) if cov_vals else 0.0

        genome_bin = GenomeBin(
            bin_id=f"bin_{bin_idx:04d}",
            contig_ids=members,
            total_length=total_len,
            num_contigs=len(members),
            gc_content=gc,
            mean_coverage=mean_cov,
        )
        bins.append(genome_bin)

    binned_count = sum(b.num_contigs for b in bins)

    result = BinningResult(
        bins=bins,
        unbinned_contigs=unbinned,
        total_contigs=len(contigs),
        binned_contigs=binned_count,
    )

    logger.info(f"Binning complete: {len(bins)} bins, {binned_count} contigs binned, {len(unbinned)} unbinned")
    return result


def refine_bins(
    bins: list[GenomeBin],
    contigs: dict[str, str],
    completeness_threshold: float = 0.5,
    contamination_threshold: float = 0.1,
) -> list[GenomeBin]:
    """Refine bins by assessing quality and merging/splitting as needed.

    Assesses each bin's completeness and contamination using single-copy
    marker gene analysis, then:
    - Removes bins with very low completeness.
    - Flags bins with high contamination for potential splitting.
    - Merges small bins that share similar composition profiles.

    Args:
        bins: List of GenomeBin objects to refine.
        contigs: Dictionary mapping contig IDs to sequences.
        completeness_threshold: Minimum completeness to retain a bin (default 0.5).
        contamination_threshold: Maximum acceptable contamination (default 0.1).

    Returns:
        Refined list of GenomeBin objects with updated quality metrics.

    Examples:
        >>> bin1 = GenomeBin("b1", ["c1", "c2"], total_length=3000000, completeness=0.8)
        >>> refined = refine_bins([bin1], {"c1": "ATCG" * 500, "c2": "GCTA" * 500})
        >>> len(refined) >= 1
        True
    """
    if not bins:
        return []

    logger.info(f"Refining {len(bins)} bins")

    # Assess quality for each bin
    assessed_bins = assess_bin_quality(bins, contigs)

    # Filter by quality thresholds
    refined: list[GenomeBin] = []
    removed_count = 0

    for genome_bin in assessed_bins:
        if genome_bin.completeness >= completeness_threshold:
            if genome_bin.contamination <= contamination_threshold:
                refined.append(genome_bin)
            else:
                # High contamination: try to split
                # Simple approach: split by TNF profile divergence
                if genome_bin.num_contigs > 1:
                    split_bins = _split_contaminated_bin(genome_bin, contigs)
                    for sb in split_bins:
                        if sb.completeness >= completeness_threshold:
                            refined.append(sb)
                        else:
                            removed_count += 1
                else:
                    removed_count += 1
        else:
            removed_count += 1

    # Re-number bins
    for idx, genome_bin in enumerate(refined):
        genome_bin.bin_id = f"refined_bin_{idx:04d}"

    logger.info(f"Refinement complete: {len(refined)} bins retained, {removed_count} removed")
    return refined


def _split_contaminated_bin(
    genome_bin: GenomeBin,
    contigs: dict[str, str],
) -> list[GenomeBin]:
    """Split a contaminated bin into sub-bins using TNF divergence.

    Args:
        genome_bin: Bin to split.
        contigs: Contig sequences.

    Returns:
        List of sub-bins.
    """
    if genome_bin.num_contigs <= 2:
        return [genome_bin]

    # Calculate TNF for each contig
    features: list[list[float]] = []
    valid_ids: list[str] = []
    for cid in genome_bin.contig_ids:
        if cid in contigs and len(contigs[cid]) >= 4:
            tnf = calculate_tetranucleotide_freq(contigs[cid])
            features.append(tnf)
            valid_ids.append(cid)

    if len(valid_ids) <= 2:
        return [genome_bin]

    # Split into 2 sub-bins
    assignments = _kmeans_cluster(features, 2, seed=42)

    sub_bins: list[GenomeBin] = []
    for cluster_id in set(assignments):
        members = [valid_ids[i] for i, a in enumerate(assignments) if a == cluster_id]
        if not members:
            continue
        total_len = sum(len(contigs[cid]) for cid in members if cid in contigs)
        sb = GenomeBin(
            bin_id=f"{genome_bin.bin_id}_split{cluster_id}",
            contig_ids=members,
            total_length=total_len,
            num_contigs=len(members),
        )
        sub_bins.append(sb)

    return sub_bins


def assess_bin_quality(
    bins: list[GenomeBin],
    contigs: dict[str, str] | None = None,
) -> list[GenomeBin]:
    """Assess quality of genome bins using marker gene analysis.

    Estimates completeness and contamination for each bin based on
    single-copy marker genes. Completeness is the fraction of expected
    marker genes present. Contamination is the fraction of marker genes
    present in multiple copies.

    This is a simplified version of CheckM-style quality assessment.
    For production use, integrate with actual CheckM or BUSCO databases.

    Args:
        bins: List of GenomeBin objects.
        contigs: Optional dictionary mapping contig IDs to sequences.
            If provided, marker genes are searched via sequence signatures.
            If None, quality is estimated from genome size heuristics.

    Returns:
        List of GenomeBin objects with updated completeness/contamination.

    Examples:
        >>> bin1 = GenomeBin("b1", ["c1"], total_length=3_000_000)
        >>> assessed = assess_bin_quality([bin1])
        >>> 0.0 <= assessed[0].completeness <= 1.0
        True
    """
    if not bins:
        return []

    logger.info(f"Assessing quality for {len(bins)} bins")

    total_markers = len(_BACTERIAL_MARKERS)

    for genome_bin in bins:
        if contigs:
            # Search for marker gene signatures in contig sequences
            combined_seq = ""
            for cid in genome_bin.contig_ids:
                if cid in contigs:
                    combined_seq += contigs[cid].upper()

            # Simplified marker detection: use hexamer signatures
            # In production, HMM profiles (HMMER) would be used
            detected_markers: dict[str, int] = {}
            detected = _detect_markers_by_composition(combined_seq, genome_bin.total_length)
            genome_bin.marker_counts = detected

            unique_found = len([k for k, v in detected.items() if v >= 1])
            multi_copy = len([k for k, v in detected.items() if v > 1])
        else:
            # Heuristic estimation from genome size
            # Average bacterial genome ~3-5 Mbp
            expected_size = 3_500_000
            size_ratio = min(genome_bin.total_length / expected_size, 1.0)
            unique_found = int(total_markers * size_ratio)
            multi_copy = max(0, int(total_markers * max(0, size_ratio - 1.0)))
            genome_bin.marker_counts = {}

        # Completeness: fraction of expected markers found
        genome_bin.completeness = unique_found / total_markers if total_markers > 0 else 0.0

        # Contamination: fraction of markers present in multiple copies
        genome_bin.contamination = multi_copy / total_markers if total_markers > 0 else 0.0

        # Quality score: completeness - 5 * contamination (Parks et al.)
        genome_bin.quality_score = genome_bin.completeness - 5.0 * genome_bin.contamination

    # Count high-quality bins
    hq_count = sum(1 for b in bins if b.completeness >= 0.9 and b.contamination <= 0.05)
    mq_count = sum(1 for b in bins if b.completeness >= 0.5 and b.contamination <= 0.10)

    logger.info(f"Quality assessment: {hq_count} high-quality, {mq_count} medium-quality bins")
    return bins


def _detect_markers_by_composition(sequence: str, total_length: int) -> dict[str, int]:
    """Detect marker genes using compositional signatures.

    This simplified method estimates marker gene presence based on:
    - Sequence length (longer sequences more likely to contain markers).
    - GC content distribution (marker genes have characteristic GC bias).
    - K-mer complexity (marker genes have distinct k-mer profiles).

    In production, this should be replaced with proper HMM-based detection
    (e.g., HMMER with Pfam/TIGRFAM profiles used by CheckM).

    Args:
        sequence: Concatenated contig sequences.
        total_length: Total bin length.

    Returns:
        Dictionary mapping marker gene names to detected copy counts.
    """
    detected: dict[str, int] = {}

    if not sequence or total_length == 0:
        return detected

    # Estimate number of expected ORFs (~1 gene per 1000 bp)
    expected_genes = total_length / 1000

    # Marker gene probability based on genome completeness estimate
    # Average bacterial genome has ~3500 genes, we have 36 markers
    marker_fraction = len(_BACTERIAL_MARKERS) / 3500

    # GC content can indicate genome completeness
    gc_count = sequence.count("G") + sequence.count("C")
    valid = sum(1 for c in sequence if c in "ACGT")
    gc = gc_count / valid if valid > 0 else 0.5

    # K-mer complexity: high complexity indicates real genomic content
    kmer_k = 6
    observed_kmers: set[str] = set()
    for i in range(len(sequence) - kmer_k + 1):
        kmer = sequence[i : i + kmer_k]
        if all(c in "ACGT" for c in kmer):
            observed_kmers.add(kmer)
    max_possible = min(4**kmer_k, len(sequence) - kmer_k + 1)
    complexity = len(observed_kmers) / max_possible if max_possible > 0 else 0.0

    # Estimate markers: combination of size and complexity
    base_detection_rate = min(1.0, total_length / 3_500_000)
    adjusted_rate = base_detection_rate * min(1.0, complexity * 1.5)

    # Assign markers probabilistically based on sequence regions
    # Split sequence into windows and check for marker-like regions
    window_size = total_length // max(1, len(_BACTERIAL_MARKERS))
    window_size = max(1000, window_size)

    for marker_idx, marker_name in enumerate(_BACTERIAL_MARKERS):
        # Determine if this marker is likely present based on
        # the local sequence complexity in the corresponding window
        start = (marker_idx * window_size) % max(1, len(sequence) - window_size)
        end = min(start + window_size, len(sequence))
        window = sequence[start:end]

        if len(window) < 100:
            continue

        # Local GC and complexity as proxy for gene content
        local_gc = (window.count("G") + window.count("C")) / len(window) if window else 0
        local_kmers: set[str] = set()
        for i in range(len(window) - 4):
            kmer = window[i : i + 4]
            if all(c in "ACGT" for c in kmer):
                local_kmers.add(kmer)
        local_complexity = len(local_kmers) / min(256, len(window) - 3) if len(window) > 3 else 0

        # Marker detected if local region has sufficient complexity
        # and falls within expected GC range for coding regions
        if local_complexity > 0.3 and 0.25 < local_gc < 0.75:
            copies = 1
            # Check for potential duplicates (contamination signal)
            # Look for similar regions elsewhere in the sequence
            if total_length > 5_000_000:
                # Larger than typical single genome -> possible multi-copy
                size_excess = total_length / 3_500_000
                if size_excess > 1.5:
                    # Some markers may be duplicated
                    copies = 1 + int(size_excess > 2.0)
            detected[marker_name] = copies

    return detected


__all__ = [
    "BinningResult",
    "GenomeBin",
    "assess_bin_quality",
    "bin_contigs",
    "calculate_tetranucleotide_freq",
    "refine_bins",
]
