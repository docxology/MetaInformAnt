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

Bin quality (completeness and contamination) is only ever computed from
single-copy marker gene annotations supplied by the caller (``marker_genes``: the
marker genes found on each contig, e.g. from an HMMER search). Without them the
values are reported as unavailable (``None``); nothing is estimated from sequence
composition or genome size.
"""

from __future__ import annotations

import itertools
import math
import random
from collections import Counter, defaultdict
from collections.abc import Collection, Mapping, Sequence
from dataclasses import dataclass, field

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

_ENV_PREFIX = "META_"

# All 256 tetranucleotides
_NUCLEOTIDES = "ACGT"
_ALL_TETRANUCLEOTIDES = ["".join(t) for t in itertools.product(_NUCLEOTIDES, repeat=4)]
_TETRA_INDEX = {t: i for i, t in enumerate(_ALL_TETRANUCLEOTIDES)}

# Universal single-copy marker genes (simplified set based on bacterial markers)
# These represent conserved genes expected once per genome. Gene symbols only: annotations
# must use the same names. Lineage-specific sets (CheckM, BUSCO, Pfam/TIGRFAM IDs) are
# passed through ``marker_set``.
BACTERIAL_SINGLE_COPY_MARKERS: tuple[str, ...] = (
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
)
_BACTERIAL_MARKERS = list(BACTERIAL_SINGLE_COPY_MARKERS)


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
        completeness: Fraction of the expected single-copy markers present (0.0 to 1.0),
            or None while unavailable (no marker annotations were supplied).
        contamination: Extra copies of expected markers divided by the number of expected
            markers (CheckM definition; 0.0 or more, and above 1.0 for a bin with
            several copies of most markers), or None while unavailable.
        quality_score: Completeness - 5 * contamination (Parks et al. convention), or
            None while unavailable.
        marker_counts: Copy number of each expected marker found in the bin.
    """

    bin_id: str
    contig_ids: list[str] = field(default_factory=list)
    total_length: int = 0
    num_contigs: int = 0
    gc_content: float = 0.0
    mean_coverage: float = 0.0
    completeness: float | None = None
    contamination: float | None = None
    quality_score: float | None = None
    marker_counts: dict[str, int] = field(default_factory=dict)


@dataclass
class BinningResult:
    """Result of metagenomic binning.

    Attributes:
        bins: List of genome bins.
        unbinned_contigs: IDs of contigs not assigned to any bin.
        total_contigs: Total number of input contigs.
        binned_contigs: Number of contigs assigned to bins.
        high_quality_bins: Number of bins with completeness >= 0.9 and contamination
            <= 0.05, or None when bin quality is unavailable (``bin_contigs`` was not
            given ``marker_genes``).
    """

    bins: list[GenomeBin]
    unbinned_contigs: list[str] = field(default_factory=list)
    total_contigs: int = 0
    binned_contigs: int = 0
    high_quality_bins: int | None = None


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
    rng = random.Random(seed)

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
    marker_genes: Mapping[str, Sequence[str]] | None = None,
    marker_set: Collection[str] | None = None,
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
        marker_genes: Optional single-copy marker gene annotations per contig
            (contig ID -> marker names found on it, repeated for multiple copies).
            When given, every bin's completeness and contamination are computed
            (see ``assess_bin_quality``) and ``high_quality_bins`` is filled in.
        marker_set: Expected marker genes for ``marker_genes`` (default
            ``BACTERIAL_SINGLE_COPY_MARKERS``).

    Returns:
        BinningResult with genome bins and statistics. Bin completeness, contamination
        and ``high_quality_bins`` are None unless ``marker_genes`` is given.

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
            high_quality_bins=0,
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

    if marker_genes is not None:
        assess_bin_quality(bins, marker_genes=marker_genes, marker_set=marker_set)

    result = BinningResult(
        bins=bins,
        unbinned_contigs=unbinned,
        total_contigs=len(contigs),
        binned_contigs=binned_count,
        high_quality_bins=count_high_quality_bins(bins),
    )

    logger.info(f"Binning complete: {len(bins)} bins, {binned_count} contigs binned, {len(unbinned)} unbinned")
    return result


def refine_bins(
    bins: list[GenomeBin],
    contigs: dict[str, str],
    completeness_threshold: float = 0.5,
    contamination_threshold: float = 0.1,
    marker_genes: Mapping[str, Sequence[str]] | None = None,
    marker_set: Collection[str] | None = None,
) -> list[GenomeBin]:
    """Refine bins by assessing quality and merging/splitting as needed.

    Assesses each bin's completeness and contamination from single-copy marker
    gene annotations (``assess_bin_quality``), then:
    - Removes bins with very low completeness.
    - Splits bins with high contamination by TNF profile; the sub-bins are assessed
      with the same annotations and kept when their completeness reaches the threshold.
    - Merging of small bins is not implemented.

    Args:
        bins: List of GenomeBin objects to refine.
        contigs: Dictionary mapping contig IDs to sequences.
        completeness_threshold: Minimum completeness to retain a bin (default 0.5).
        contamination_threshold: Maximum acceptable contamination (default 0.1).
        marker_genes: Single-copy marker gene annotations per contig (contig ID -> marker
            names found on it). Required whenever ``bins`` is not empty: without
            them neither completeness nor contamination exists to refine on.
        marker_set: Expected marker genes (default ``BACTERIAL_SINGLE_COPY_MARKERS``).

    Returns:
        Refined list of GenomeBin objects with updated quality metrics.

    Raises:
        ValueError: If ``bins`` is not empty and ``marker_genes`` is None.

    Examples:
        >>> bin1 = GenomeBin("b1", ["c1", "c2"], total_length=3000000)
        >>> markers = {"c1": list(BACTERIAL_SINGLE_COPY_MARKERS)}
        >>> refined = refine_bins([bin1], {"c1": "ATCG" * 500, "c2": "GCTA" * 500}, marker_genes=markers)
        >>> len(refined) >= 1
        True
    """
    if not bins:
        return []
    if marker_genes is None:
        raise ValueError(
            "refine_bins requires marker_genes (single-copy marker gene annotations per contig): "
            "completeness and contamination cannot be estimated without them"
        )

    logger.info(f"Refining {len(bins)} bins")

    # Assess quality for each bin
    assessed_bins = assess_bin_quality(bins, contigs, marker_genes=marker_genes, marker_set=marker_set)

    # Filter by quality thresholds
    refined: list[GenomeBin] = []
    removed_count = 0

    for genome_bin in assessed_bins:
        completeness, contamination = _quality_values(genome_bin)
        if completeness >= completeness_threshold:
            if contamination <= contamination_threshold:
                refined.append(genome_bin)
            else:
                # High contamination: try to split
                # Simple approach: split by TNF profile divergence
                if genome_bin.num_contigs > 1:
                    split_bins = _split_contaminated_bin(genome_bin, contigs)
                    assess_bin_quality(split_bins, marker_genes=marker_genes, marker_set=marker_set)
                    for sb in split_bins:
                        if _quality_values(sb)[0] >= completeness_threshold:
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


def _quality_values(genome_bin: GenomeBin) -> tuple[float, float]:
    """Return ``(completeness, contamination)`` of an assessed bin.

    Raises:
        ValueError: If the bin has not been assessed.
    """
    if genome_bin.completeness is None or genome_bin.contamination is None:
        raise ValueError(f"Bin {genome_bin.bin_id!r} has no completeness/contamination; assess it with marker_genes")
    return genome_bin.completeness, genome_bin.contamination


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
    marker_genes: Mapping[str, Sequence[str]] | None = None,
    marker_set: Collection[str] | None = None,
) -> list[GenomeBin]:
    """Assess quality of genome bins from single-copy marker gene annotations.

    CheckM-style estimates for a flat marker list of ``M`` expected single-copy genes:

    - completeness = (number of expected markers present at least once in the bin) / M
    - contamination = (extra copies, i.e. sum over markers of ``copies - 1`` where
      ``copies >= 1``) / M
    - quality_score = completeness - 5 * contamination (Parks et al.)

    Marker genes must come from a real search of the contigs (for example HMMER with
    TIGRFAM/Pfam marker profiles); the search is not done here. Without
    ``marker_genes`` nothing can be said about completeness or contamination, so both
    (and ``quality_score``) are set to None for every bin: they are never estimated from
    genome size or sequence composition.

    Args:
        bins: List of GenomeBin objects (updated in place).
        contigs: Accepted for backward compatibility and unused: no marker search is
            performed on sequences.
        marker_genes: Marker gene annotations per contig (contig ID -> marker names found
            on that contig; a marker listed twice counts as two copies). Contigs that are
            not in the mapping are treated as carrying no markers.
        marker_set: Expected marker genes defining the denominator ``M``. Defaults to
            ``BACTERIAL_SINGLE_COPY_MARKERS``; pass the set matching the annotation names
            (for example a lineage-specific CheckM or BUSCO set).

    Returns:
        List of GenomeBin objects with updated completeness/contamination
        (None when ``marker_genes`` is None).

    Raises:
        ValueError: If ``marker_set`` is given but empty.

    Examples:
        >>> bin1 = GenomeBin("b1", ["c1"], total_length=3_000_000)
        >>> assessed = assess_bin_quality([bin1])
        >>> assessed[0].completeness is None
        True
        >>> markers = {"c1": list(BACTERIAL_SINGLE_COPY_MARKERS)}
        >>> assess_bin_quality([bin1], marker_genes=markers)[0].completeness
        1.0
    """
    if not bins:
        return []

    if marker_genes is None:
        logger.warning(
            f"No marker_genes supplied: completeness and contamination of {len(bins)} bins are unavailable "
            "(no marker search is performed on sequences)"
        )
        for genome_bin in bins:
            genome_bin.completeness = None
            genome_bin.contamination = None
            genome_bin.quality_score = None
            genome_bin.marker_counts = {}
        return bins

    expected = set(BACTERIAL_SINGLE_COPY_MARKERS) if marker_set is None else set(marker_set)
    if not expected:
        raise ValueError("marker_set must contain at least one marker gene")

    logger.info(f"Assessing quality for {len(bins)} bins against {len(expected)} marker genes")

    if marker_genes and not any(cid in marker_genes for b in bins for cid in b.contig_ids):
        logger.warning("None of the contig IDs in marker_genes occur in the bins; check the contig naming")

    total_markers = len(expected)

    for genome_bin in bins:
        counts: Counter[str] = Counter()
        for cid in genome_bin.contig_ids:
            annotated = marker_genes.get(cid, ())
            for marker in [annotated] if isinstance(annotated, str) else annotated:
                if marker in expected:
                    counts[marker] += 1

        extra_copies = sum(copies - 1 for copies in counts.values())

        genome_bin.marker_counts = dict(counts)
        genome_bin.completeness = len(counts) / total_markers
        genome_bin.contamination = extra_copies / total_markers
        genome_bin.quality_score = genome_bin.completeness - 5.0 * genome_bin.contamination

    hq_count = count_high_quality_bins(bins)
    mq_count = sum(
        1
        for b in bins
        if b.completeness is not None
        and b.contamination is not None
        and b.completeness >= 0.5
        and b.contamination <= 0.10
    )

    logger.info(f"Quality assessment: {hq_count} high-quality, {mq_count} medium-quality bins")
    return bins


def count_high_quality_bins(
    bins: list[GenomeBin],
    min_completeness: float = 0.9,
    max_contamination: float = 0.05,
) -> int | None:
    """Count bins with completeness >= ``min_completeness`` and contamination <= ``max_contamination``.

    Args:
        bins: Bins whose quality was assessed with ``assess_bin_quality``.
        min_completeness: Completeness threshold (default 0.9).
        max_contamination: Contamination threshold (default 0.05).

    Returns:
        The number of high-quality bins (0 for no bins), or None if any bin has no
        completeness/contamination (quality unavailable, so the count is unknown).
    """
    count = 0
    for genome_bin in bins:
        if genome_bin.completeness is None or genome_bin.contamination is None:
            return None
        if genome_bin.completeness >= min_completeness and genome_bin.contamination <= max_contamination:
            count += 1
    return count


__all__ = [
    "BACTERIAL_SINGLE_COPY_MARKERS",
    "BinningResult",
    "GenomeBin",
    "assess_bin_quality",
    "bin_contigs",
    "calculate_tetranucleotide_freq",
    "count_high_quality_bins",
    "refine_bins",
]
