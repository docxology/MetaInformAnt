"""Haplotype phasing from long-read sequencing data.

Implements read-based phasing using heterozygous variants to assign reads
to haplotypes and construct phase blocks. Uses a graph-based approach where
variants are nodes and reads connecting variants are edges.

The phasing algorithm:
1. Build a read-variant matrix from heterozygous SNPs
2. Construct a fragment graph connecting variants observed on the same reads
3. Partition the graph into two haplotypes using a max-cut approximation
4. Extend phase blocks across connected components
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Sequence

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


@dataclass
class Variant:
    """A heterozygous variant used for phasing.

    Attributes:
        chromosome: Reference chromosome.
        position: 0-based genomic position.
        ref_allele: Reference allele.
        alt_allele: Alternate allele.
        genotype: Genotype string (e.g., '0/1').
    """

    chromosome: str
    position: int
    ref_allele: str = ""
    alt_allele: str = ""
    genotype: str = "0/1"


@dataclass
class ReadVariantObservation:
    """Observation of a variant on a specific read.

    Attributes:
        read_name: Read identifier.
        variant_index: Index into the variants list.
        allele: Observed allele (0 for ref, 1 for alt).
        quality: Base quality at the variant position.
    """

    read_name: str
    variant_index: int
    allele: int
    quality: int = 30


@dataclass
class PhaseBlock:
    """A contiguous block of phased variants.

    Attributes:
        block_id: Unique block identifier.
        chromosome: Reference chromosome.
        start: Start position of the phase block.
        end: End position of the phase block.
        variants: List of variant positions in the block.
        haplotype1: Alleles on haplotype 1 (list of 0/1 per variant).
        haplotype2: Alleles on haplotype 2 (list of 0/1 per variant).
        num_variants: Number of variants in the block.
        quality: Phase quality score.
    """

    block_id: int = 0
    chromosome: str = ""
    start: int = 0
    end: int = 0
    variants: list[int] = field(default_factory=list)
    haplotype1: list[int] = field(default_factory=list)
    haplotype2: list[int] = field(default_factory=list)
    num_variants: int = 0
    quality: float = 0.0


@dataclass
class PhaseBlockStats:
    """Summary statistics for a set of phase blocks.

    Attributes:
        num_blocks: Total number of phase blocks.
        num_phased_variants: Total phased variants across all blocks.
        n50: N50 of phase block sizes (in base pairs).
        max_block_size: Largest phase block (bp).
        mean_block_size: Mean phase block size (bp).
        median_block_size: Median phase block size (bp).
        switch_error_rate: Estimated switch error rate.
        total_phased_span: Total genomic span covered by phase blocks.
    """

    num_blocks: int = 0
    num_phased_variants: int = 0
    n50: int = 0
    max_block_size: int = 0
    mean_block_size: float = 0.0
    median_block_size: float = 0.0
    switch_error_rate: float = 0.0
    total_phased_span: int = 0


def phase_reads(
    variants: Sequence[Variant | dict[str, Any]],
    reads: Sequence[dict[str, Any]],
    min_variant_quality: int = 20,
    min_base_quality: int = 10,
) -> tuple[list[PhaseBlock], dict[str, int]]:
    """Phase reads using heterozygous variants.

    Assigns reads to haplotype 1 or haplotype 2 based on their alleles at
    heterozygous variant sites. Uses a greedy graph-based approach.

    Algorithm:
    1. Extract variant observations from each read
    2. Build a weighted graph where edges connect co-observed variants
    3. Partition into two haplotypes using spectral-like max-cut
    4. Assign reads to haplotypes based on majority vote

    Args:
        variants: Sequence of Variant objects or dictionaries with keys:
            chromosome, position, ref_allele, alt_allele.
        reads: Sequence of read dictionaries, each containing:
            - read_name (str): Unique read identifier
            - variants (list): List of dicts with 'position', 'allele', 'quality'
            OR
            - aligned_pairs (list): List of (query_pos, ref_pos) tuples
            - query_sequence (str): Read sequence
        min_variant_quality: Minimum variant quality to use for phasing.
        min_base_quality: Minimum base quality at variant positions.

    Returns:
        Tuple of (phase_blocks, read_assignments) where read_assignments maps
        read_name -> haplotype (1 or 2).
    """
    # Normalize variants
    norm_variants: list[Variant] = []
    for v in variants:
        if isinstance(v, dict):
            norm_variants.append(Variant(
                chromosome=v.get("chromosome", ""),
                position=v.get("position", 0),
                ref_allele=v.get("ref_allele", ""),
                alt_allele=v.get("alt_allele", ""),
                genotype=v.get("genotype", "0/1"),
            ))
        else:
            norm_variants.append(v)

    if not norm_variants:
        return [], {}

    # Build variant position index
    var_positions: dict[tuple[str, int], int] = {}
    for i, v in enumerate(norm_variants):
        var_positions[(v.chromosome, v.position)] = i

    # Extract read-variant observations
    observations: list[ReadVariantObservation] = []
    reads_per_variant: dict[int, list[ReadVariantObservation]] = defaultdict(list)
    variants_per_read: dict[str, list[ReadVariantObservation]] = defaultdict(list)

    for read in reads:
        read_name = read.get("read_name", read.get("read_id", ""))
        if not read_name:
            continue

        # Get variant observations from the read
        read_vars = read.get("variants", [])
        if read_vars:
            for rv in read_vars:
                pos = rv.get("position", 0)
                chrom = rv.get("chromosome", read.get("reference_name", ""))
                key = (chrom, pos)
                if key in var_positions:
                    allele = int(rv.get("allele", 0))
                    quality = int(rv.get("quality", 30))
                    if quality >= min_base_quality:
                        obs = ReadVariantObservation(
                            read_name=read_name,
                            variant_index=var_positions[key],
                            allele=allele,
                            quality=quality,
                        )
                        observations.append(obs)
                        reads_per_variant[var_positions[key]].append(obs)
                        variants_per_read[read_name].append(obs)

    if not observations:
        logger.warning("No read-variant observations found for phasing")
        return [], {}

    # Build fragment connection graph
    # Edge weight = number of reads supporting the connection, with sign
    # indicating whether alleles are in phase (+) or anti-phase (-)
    edge_weights: dict[tuple[int, int], float] = defaultdict(float)

    for read_name, read_obs_list in variants_per_read.items():
        if len(read_obs_list) < 2:
            continue

        # For each pair of variants on this read, add an edge
        for i in range(len(read_obs_list)):
            for j in range(i + 1, len(read_obs_list)):
                vi = read_obs_list[i].variant_index
                vj = read_obs_list[j].variant_index
                if vi > vj:
                    vi, vj = vj, vi

                # Same allele = in-phase (+1), different allele = anti-phase (-1)
                ai = read_obs_list[i].allele
                aj = read_obs_list[j].allele
                weight = 1.0 if ai == aj else -1.0

                edge_weights[(vi, vj)] += weight

    # Graph-based phasing: partition variants into two haplotypes
    # Using a greedy max-cut approximation
    num_vars = len(norm_variants)
    haplotype_assignment = _greedy_max_cut(num_vars, edge_weights)

    # Build phase blocks from connected components
    # Find connected components of the variant graph
    components = _find_connected_components(num_vars, set(edge_weights.keys()))

    phase_blocks: list[PhaseBlock] = []
    for block_id, component in enumerate(components):
        if len(component) < 2:
            continue

        sorted_indices = sorted(component)
        block_variants = [norm_variants[i] for i in sorted_indices]

        chrom = block_variants[0].chromosome
        start = min(v.position for v in block_variants)
        end = max(v.position for v in block_variants)

        h1 = [haplotype_assignment[i] for i in sorted_indices]
        h2 = [1 - h for h in h1]

        # Quality based on consistency of edge weights
        block_quality = _compute_block_quality(sorted_indices, edge_weights, haplotype_assignment)

        phase_blocks.append(PhaseBlock(
            block_id=block_id,
            chromosome=chrom,
            start=start,
            end=end,
            variants=[v.position for v in block_variants],
            haplotype1=h1,
            haplotype2=h2,
            num_variants=len(sorted_indices),
            quality=block_quality,
        ))

    # Assign reads to haplotypes based on majority vote
    read_assignments: dict[str, int] = {}
    for read_name, read_obs_list in variants_per_read.items():
        h1_support = 0
        h2_support = 0
        for obs in read_obs_list:
            if haplotype_assignment[obs.variant_index] == obs.allele:
                h1_support += 1
            else:
                h2_support += 1

        if h1_support > h2_support:
            read_assignments[read_name] = 1
        elif h2_support > h1_support:
            read_assignments[read_name] = 2
        # If tied, leave unassigned

    logger.info(
        "Phased %d variants into %d blocks, assigned %d reads to haplotypes",
        sum(b.num_variants for b in phase_blocks),
        len(phase_blocks),
        len(read_assignments),
    )

    return phase_blocks, read_assignments


def build_haplotype_blocks(
    phased_variants: Sequence[dict[str, Any]],
    max_gap: int = 300000,
) -> list[PhaseBlock]:
    """Construct haplotype phase blocks from pre-phased variant data.

    Groups consecutive phased variants into blocks, breaking at gaps
    larger than max_gap or at phase set boundaries.

    Args:
        phased_variants: Sequence of phased variant dictionaries, each containing:
            - chromosome (str)
            - position (int)
            - haplotype (int): 1 or 2
            - allele (int): 0 (ref) or 1 (alt)
            - phase_set (int, optional): Phase set ID for grouping
        max_gap: Maximum gap between variants before starting a new block.

    Returns:
        List of PhaseBlock objects.
    """
    if not phased_variants:
        return []

    # Group by chromosome and phase_set
    groups: dict[tuple[str, int], list[dict[str, Any]]] = defaultdict(list)
    for v in phased_variants:
        chrom = v.get("chromosome", "")
        ps = v.get("phase_set", 0)
        groups[(chrom, ps)].append(v)

    blocks: list[PhaseBlock] = []
    block_id = 0

    for (chrom, ps), group_vars in groups.items():
        sorted_vars = sorted(group_vars, key=lambda v: v["position"])

        current_block_vars: list[dict[str, Any]] = [sorted_vars[0]]

        for i in range(1, len(sorted_vars)):
            gap = sorted_vars[i]["position"] - sorted_vars[i - 1]["position"]
            if gap > max_gap:
                # Start a new block
                if len(current_block_vars) >= 2:
                    blocks.append(_make_phase_block(block_id, chrom, current_block_vars))
                    block_id += 1
                current_block_vars = [sorted_vars[i]]
            else:
                current_block_vars.append(sorted_vars[i])

        if len(current_block_vars) >= 2:
            blocks.append(_make_phase_block(block_id, chrom, current_block_vars))
            block_id += 1

    return blocks


def _make_phase_block(block_id: int, chrom: str, variants: list[dict[str, Any]]) -> PhaseBlock:
    """Create a PhaseBlock from a list of variant dictionaries."""
    positions = [v["position"] for v in variants]
    h1 = []
    h2 = []
    for v in variants:
        hp = v.get("haplotype", 1)
        allele = v.get("allele", 0)
        if hp == 1:
            h1.append(allele)
            h2.append(1 - allele)
        else:
            h1.append(1 - allele)
            h2.append(allele)

    return PhaseBlock(
        block_id=block_id,
        chromosome=chrom,
        start=min(positions),
        end=max(positions),
        variants=positions,
        haplotype1=h1,
        haplotype2=h2,
        num_variants=len(positions),
        quality=1.0,
    )


def tag_reads_by_haplotype(
    reads: Sequence[dict[str, Any]],
    phase_blocks: Sequence[PhaseBlock],
) -> list[dict[str, Any]]:
    """Assign reads to haplotypes based on overlap with phase blocks.

    For each read, determines which phase block it overlaps, then assigns
    it to haplotype 1 or 2 based on allele concordance with the phase block.

    Args:
        reads: Sequence of read dictionaries with keys:
            - read_name (str)
            - reference_name (str): Chromosome
            - reference_start (int): Alignment start
            - reference_end (int): Alignment end
            - variants (list, optional): Observed variants
        phase_blocks: Sequence of PhaseBlock objects.

    Returns:
        List of read dictionaries with added 'HP' (haplotype) and 'PS'
        (phase set) tags.
    """
    # Build interval index for phase blocks
    block_index: dict[str, list[PhaseBlock]] = defaultdict(list)
    for block in phase_blocks:
        block_index[block.chromosome].append(block)

    # Sort blocks by start position for each chromosome
    for chrom in block_index:
        block_index[chrom].sort(key=lambda b: b.start)

    tagged_reads: list[dict[str, Any]] = []

    for read in reads:
        read_dict = dict(read)
        chrom = read_dict.get("reference_name", "")
        read_start = read_dict.get("reference_start", 0)
        read_end = read_dict.get("reference_end", read_start)

        assigned_hp = 0
        assigned_ps = 0

        if chrom in block_index:
            for block in block_index[chrom]:
                # Check overlap
                if read_start < block.end and read_end > block.start:
                    # Read overlaps this phase block
                    # Determine haplotype from variant observations
                    read_vars = read_dict.get("variants", [])
                    if read_vars:
                        h1_support = 0
                        h2_support = 0

                        for rv in read_vars:
                            pos = rv.get("position", 0)
                            allele = rv.get("allele", 0)

                            if pos in block.variants:
                                idx = block.variants.index(pos)
                                if idx < len(block.haplotype1):
                                    if block.haplotype1[idx] == allele:
                                        h1_support += 1
                                    elif block.haplotype2[idx] == allele:
                                        h2_support += 1

                        if h1_support > h2_support:
                            assigned_hp = 1
                            assigned_ps = block.block_id
                        elif h2_support > h1_support:
                            assigned_hp = 2
                            assigned_ps = block.block_id
                    break  # Only use the first overlapping block

        read_dict["HP"] = assigned_hp
        read_dict["PS"] = assigned_ps
        tagged_reads.append(read_dict)

    tagged_count = sum(1 for r in tagged_reads if r["HP"] > 0)
    logger.info("Tagged %d/%d reads with haplotype assignments", tagged_count, len(tagged_reads))

    return tagged_reads


def calculate_phase_block_stats(
    blocks: Sequence[PhaseBlock],
) -> PhaseBlockStats:
    """Calculate summary statistics for a set of phase blocks.

    Computes N50, size distribution, and switch error rate estimates
    for a collection of phase blocks.

    Args:
        blocks: Sequence of PhaseBlock objects.

    Returns:
        PhaseBlockStats with comprehensive statistics.
    """
    if not blocks:
        return PhaseBlockStats()

    sizes = [b.end - b.start for b in blocks if b.end > b.start]
    if not sizes:
        return PhaseBlockStats(num_blocks=len(blocks))

    num_variants = sum(b.num_variants for b in blocks)
    total_span = sum(sizes)

    sorted_sizes = sorted(sizes)
    n = len(sorted_sizes)

    # Mean and median
    mean_size = sum(sorted_sizes) / n
    if n % 2 == 0:
        median_size = (sorted_sizes[n // 2 - 1] + sorted_sizes[n // 2]) / 2.0
    else:
        median_size = float(sorted_sizes[n // 2])

    # N50 calculation
    n50 = 0
    cumulative = 0
    threshold = total_span * 0.5
    for size in sorted(sizes, reverse=True):
        cumulative += size
        if cumulative >= threshold:
            n50 = size
            break

    # Estimate switch error rate from quality scores
    # Higher quality blocks have lower switch error rates
    total_quality = sum(b.quality * b.num_variants for b in blocks if b.num_variants > 0)
    total_vars_quality = sum(b.num_variants for b in blocks if b.num_variants > 0)
    if total_vars_quality > 0:
        mean_quality = total_quality / total_vars_quality
        # Convert quality to approximate switch error rate
        # Using a simple model: SER ~ 1 / (1 + quality * 10)
        switch_error_rate = 1.0 / (1.0 + mean_quality * 10.0)
    else:
        switch_error_rate = 0.0

    return PhaseBlockStats(
        num_blocks=len(blocks),
        num_phased_variants=num_variants,
        n50=n50,
        max_block_size=max(sizes),
        mean_block_size=mean_size,
        median_block_size=median_size,
        switch_error_rate=switch_error_rate,
        total_phased_span=total_span,
    )


# --- Internal helper functions ---


def _greedy_max_cut(
    num_nodes: int,
    edge_weights: dict[tuple[int, int], float],
) -> list[int]:
    """Greedy max-cut approximation for graph bipartitioning.

    Assigns each node to partition 0 or 1 to maximize the total weight of
    edges crossing the cut. Uses a greedy approach that considers each node
    and assigns it to the partition that maximizes the cut.

    Args:
        num_nodes: Number of nodes.
        edge_weights: Dictionary mapping (i, j) -> weight.

    Returns:
        List of partition assignments (0 or 1) for each node.
    """
    assignment = [0] * num_nodes

    # Build adjacency list
    adj: dict[int, list[tuple[int, float]]] = defaultdict(list)
    for (i, j), w in edge_weights.items():
        adj[i].append((j, w))
        adj[j].append((i, w))

    # Greedy assignment
    # Start with the first node in partition 0
    for node in range(1, num_nodes):
        # Calculate benefit of placing in partition 0 vs 1
        benefit_0 = 0.0  # Total weight of cut edges if in partition 0
        benefit_1 = 0.0  # Total weight of cut edges if in partition 1

        for neighbor, weight in adj[node]:
            if neighbor >= node:
                continue  # Only consider already-assigned nodes

            if assignment[neighbor] == 0:
                benefit_1 += weight  # In partition 1, this crosses the cut
            else:
                benefit_0 += weight  # In partition 0, this crosses the cut

        # For phasing, negative weights mean anti-phase (should be in different partitions)
        # Positive weights mean in-phase (should be in same partition)
        # We want to maximize: sum of |w| for edges crossing the cut if w < 0,
        # plus sum of |w| for edges not crossing the cut if w > 0
        score_0 = 0.0
        score_1 = 0.0

        for neighbor, weight in adj[node]:
            if neighbor >= node:
                continue

            if weight > 0:
                # In-phase: want same partition
                if assignment[neighbor] == 0:
                    score_0 += weight
                else:
                    score_1 += weight
            else:
                # Anti-phase: want different partition
                if assignment[neighbor] == 0:
                    score_1 += abs(weight)
                else:
                    score_0 += abs(weight)

        assignment[node] = 0 if score_0 >= score_1 else 1

    return assignment


def _find_connected_components(
    num_nodes: int,
    edges: set[tuple[int, int]],
) -> list[list[int]]:
    """Find connected components in an undirected graph.

    Args:
        num_nodes: Number of nodes.
        edges: Set of (i, j) edges.

    Returns:
        List of connected components, each a list of node indices.
    """
    adj: dict[int, set[int]] = defaultdict(set)
    nodes_with_edges: set[int] = set()

    for i, j in edges:
        adj[i].add(j)
        adj[j].add(i)
        nodes_with_edges.add(i)
        nodes_with_edges.add(j)

    visited: set[int] = set()
    components: list[list[int]] = []

    for node in sorted(nodes_with_edges):
        if node in visited:
            continue

        # BFS
        component: list[int] = []
        queue = [node]
        visited.add(node)

        while queue:
            current = queue.pop(0)
            component.append(current)

            for neighbor in adj[current]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)

        components.append(sorted(component))

    return components


def _compute_block_quality(
    variant_indices: list[int],
    edge_weights: dict[tuple[int, int], float],
    assignment: list[int],
) -> float:
    """Compute phase block quality based on consistency of edges.

    Quality is the fraction of edges that are consistent with the phasing
    (in-phase variants in same partition, anti-phase in different partitions).
    """
    consistent = 0
    total = 0

    for i in range(len(variant_indices)):
        for j in range(i + 1, len(variant_indices)):
            vi = variant_indices[i]
            vj = variant_indices[j]
            key = (min(vi, vj), max(vi, vj))

            if key in edge_weights:
                weight = edge_weights[key]
                same_partition = assignment[vi] == assignment[vj]
                total += 1

                if (weight > 0 and same_partition) or (weight < 0 and not same_partition):
                    consistent += 1

    if total == 0:
        return 0.0

    return consistent / total
