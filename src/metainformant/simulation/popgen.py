"""Population genetics synthetic data generation.

This module provides methods to generate synthetic population genetics data
with specified properties, including diversity, Fst, allele frequencies,
demographic scenarios, and evolutionary models.
"""

from __future__ import annotations

import math
import random
from collections.abc import Sequence
from typing import Any

from ..simulation.sequences import _DNA, generate_random_dna, mutate_sequence


def generate_population_sequences(
    n_sequences: int,
    sequence_length: int,
    *,
    nucleotide_diversity: float | None = None,
    wattersons_theta: float | None = None,
    reference_sequence: str | None = None,
    mutation_rate: float = 0.001,
    gc_content: float = 0.5,
    rng: random.Random | None = None,
) -> list[str]:
    """Generate a population of sequences with specified diversity.
    
    Generates sequences that approximate target nucleotide diversity (π) or
    Watterson's theta (θ_W). Uses a coalescent-inspired approach to introduce
    mutations that create realistic polymorphism patterns.
    
    Args:
        n_sequences: Number of sequences to generate
        sequence_length: Length of each sequence
        nucleotide_diversity: Target π (average pairwise differences per site).
            If provided, mutations are introduced to approximate this value.
        wattersons_theta: Target θ_W (Watterson's estimator). Alternative to
            nucleotide_diversity. If both provided, nucleotide_diversity takes precedence.
        reference_sequence: Starting sequence. If None, generates random sequence.
        mutation_rate: Per-site mutation rate (used if diversity targets not specified)
        gc_content: GC content for reference sequence (if generating new one)
        rng: Random number generator
    
    Returns:
        List of sequences with specified diversity properties
    
    Examples:
        >>> seqs = generate_population_sequences(
        ...     n_sequences=10,
        ...     sequence_length=1000,
        ...     nucleotide_diversity=0.01
        ... )
        >>> len(seqs)
        10
        >>> len(seqs[0])
        1000
    """
    r = rng or random
    
    # Generate or use reference sequence
    if reference_sequence is None:
        reference_sequence = generate_random_dna(
            sequence_length, gc_content=gc_content, rng=r
        )
    else:
        if len(reference_sequence) != sequence_length:
            raise ValueError(
                f"Reference sequence length {len(reference_sequence)} != "
                f"sequence_length {sequence_length}"
            )
    
    # Calculate target number of mutations
    if nucleotide_diversity is not None:
        # π ≈ average number of differences per site
        # For n sequences, total pairwise differences ≈ n(n-1)/2 * π * L
        # Each mutation creates ~n differences (one per sequence)
        # So target mutations ≈ n(n-1)/2 * π * L / n = (n-1)/2 * π * L
        target_mutations = int((n_sequences - 1) / 2 * nucleotide_diversity * sequence_length)
    elif wattersons_theta is not None:
        # θ_W ≈ S / (a₁ * L) where S is segregating sites
        # a₁ = Σ(1/i) for i=1 to n-1
        a1 = sum(1.0 / i for i in range(1, n_sequences))
        target_segregating = int(wattersons_theta * a1 * sequence_length)
        # Approximate: each mutation creates one segregating site
        target_mutations = target_segregating
    else:
        # Use mutation_rate to determine mutations
        # Expected mutations per sequence = mutation_rate * L
        target_mutations = int(mutation_rate * sequence_length * n_sequences)
    
    # Generate sequences by introducing mutations
    sequences = [reference_sequence]
    
    # Introduce mutations progressively
    for _ in range(n_sequences - 1):
        # Start from a randomly chosen existing sequence
        parent_seq = r.choice(sequences)
        
        # Introduce mutations to create new sequence
        # Number of mutations per new sequence
        mutations_per_seq = max(1, target_mutations // (n_sequences - 1))
        
        new_seq = mutate_sequence(parent_seq, n_mut=mutations_per_seq, rng=r)
        sequences.append(new_seq)
    
    # If we need more diversity, add additional mutations
    current_pi = _estimate_pi(sequences)
    if nucleotide_diversity is not None and current_pi < nucleotide_diversity * 0.9:
        # Add more mutations to reach target
        additional_mutations = int(
            (nucleotide_diversity - current_pi) * sequence_length * n_sequences / 2
        )
        for _ in range(additional_mutations):
            seq_idx = r.randint(0, len(sequences) - 1)
            pos = r.randint(0, sequence_length - 1)
            seq_list = list(sequences[seq_idx])
            bases = [b for b in _DNA if b != seq_list[pos]]
            if bases:
                seq_list[pos] = r.choice(bases)
                sequences[seq_idx] = "".join(seq_list)
    
    return sequences


def generate_two_populations(
    n_pop1: int,
    n_pop2: int,
    sequence_length: int,
    *,
    fst: float = 0.1,
    within_pop_diversity: float = 0.01,
    reference_sequence: str | None = None,
    gc_content: float = 0.5,
    rng: random.Random | None = None,
) -> tuple[list[str], list[str]]:
    """Generate two populations with specified Fst.
    
    Creates two populations that are differentiated by the specified Fst value.
    Higher Fst means more differentiation between populations.
    
    Args:
        n_pop1: Number of sequences in population 1
        n_pop2: Number of sequences in population 2
        sequence_length: Length of each sequence
        fst: Target Fst value (0-1). Higher values = more differentiation
        within_pop_diversity: Target nucleotide diversity within each population
        reference_sequence: Starting sequence (ancestral)
        gc_content: GC content for reference sequence
        rng: Random number generator
    
    Returns:
        Tuple of (pop1_sequences, pop2_sequences)
    
    Examples:
        >>> pop1, pop2 = generate_two_populations(
        ...     n_pop1=10,
        ...     n_pop2=10,
        ...     sequence_length=1000,
        ...     fst=0.2
        ... )
        >>> len(pop1), len(pop2)
        (10, 10)
    """
    r = rng or random
    
    # Generate ancestral sequence
    if reference_sequence is None:
        reference_sequence = generate_random_dna(
            sequence_length, gc_content=gc_content, rng=r
        )
    
    # Generate population 1
    pop1 = generate_population_sequences(
        n_pop1,
        sequence_length,
        nucleotide_diversity=within_pop_diversity,
        reference_sequence=reference_sequence,
        rng=r,
    )
    
    # Create population 2 with differentiation
    # Fst = variance_between / (variance_between + variance_within)
    # To achieve target Fst, we need to introduce population-specific mutations
    
    # Calculate number of fixed differences needed
    # Fst ≈ fixed_differences / (fixed_differences + shared_polymorphisms)
    # For simplicity, introduce mutations that are fixed in one population
    
    num_fixed_differences = int(fst * sequence_length * 0.1)  # Rough estimate
    
    # Generate population 2 starting from ancestral
    pop2_seed = reference_sequence
    for _ in range(num_fixed_differences):
        pos = r.randint(0, sequence_length - 1)
        seq_list = list(pop2_seed)
        bases = [b for b in _DNA if b != seq_list[pos]]
        if bases:
            seq_list[pos] = r.choice(bases)
            pop2_seed = "".join(seq_list)
    
    # Generate population 2 from differentiated seed
    pop2 = generate_population_sequences(
        n_pop2,
        sequence_length,
        nucleotide_diversity=within_pop_diversity,
        reference_sequence=pop2_seed,
        rng=r,
    )
    
    return pop1, pop2


def generate_genotype_matrix(
    n_individuals: int,
    n_sites: int,
    *,
    allele_frequencies: Sequence[float] | None = None,
    min_maf: float = 0.05,
    max_maf: float = 0.5,
    hwe: bool = True,
    ploidy: int = 2,
    rng: random.Random | None = None,
) -> list[list[int]]:
    """Generate a genotype matrix with specified allele frequencies.
    
    Creates a genotype matrix where each individual is represented by a row
    and each site by a column. Genotypes can be generated under Hardy-Weinberg
    equilibrium or with specified allele frequencies.
    
    Args:
        n_individuals: Number of individuals (rows)
        n_sites: Number of sites (columns)
        allele_frequencies: Optional list of allele frequencies per site.
            If provided, length must equal n_sites. If None, frequencies are
            randomly sampled between min_maf and max_maf.
        min_maf: Minimum minor allele frequency (if generating frequencies)
        max_maf: Maximum minor allele frequency (if generating frequencies)
        hwe: If True, genotypes follow Hardy-Weinberg equilibrium
        ploidy: Ploidy level (2 for diploid, 1 for haploid)
        rng: Random number generator
    
    Returns:
        Genotype matrix as list of lists. For diploid: 0=AA, 1=AB, 2=BB.
        For haploid: 0=reference, 1=alternate.
    
    Examples:
        >>> genotypes = generate_genotype_matrix(
        ...     n_individuals=10,
        ...     n_sites=5,
        ...     allele_frequencies=[0.2, 0.3, 0.4, 0.1, 0.5]
        ... )
        >>> len(genotypes)
        10
        >>> len(genotypes[0])
        5
    """
    r = rng or random
    
    # Generate or use allele frequencies
    if allele_frequencies is None:
        freqs = [r.uniform(min_maf, max_maf) for _ in range(n_sites)]
    else:
        if len(allele_frequencies) != n_sites:
            raise ValueError(
                f"allele_frequencies length {len(allele_frequencies)} != n_sites {n_sites}"
            )
        freqs = list(allele_frequencies)
    
    # Generate genotypes
    genotypes = []
    for _ in range(n_individuals):
        individual = []
        for freq in freqs:
            if ploidy == 2:
                # Diploid: Hardy-Weinberg equilibrium
                if hwe:
                    # Genotype frequencies: p², 2pq, q²
                    p = freq
                    q = 1 - p
                    rand = r.random()
                    if rand < p * p:
                        genotype = 0  # AA (homozygous reference)
                    elif rand < p * p + 2 * p * q:
                        genotype = 1  # AB (heterozygous)
                    else:
                        genotype = 2  # BB (homozygous alternate)
                else:
                    # Random assignment based on allele frequency
                    # Sample two alleles independently
                    allele1 = 1 if r.random() < freq else 0
                    allele2 = 1 if r.random() < freq else 0
                    genotype = allele1 + allele2
            else:
                # Haploid
                genotype = 1 if r.random() < freq else 0
            
            individual.append(genotype)
        genotypes.append(individual)
    
    return genotypes


def simulate_bottleneck_population(
    n_sequences: int,
    sequence_length: int,
    *,
    pre_bottleneck_diversity: float = 0.01,
    bottleneck_size: int = 5,
    bottleneck_duration: int = 10,
    recovery_generations: int = 20,
    mutation_rate: float = 0.001,
    rng: random.Random | None = None,
) -> list[str]:
    """Simulate a population that went through a bottleneck.
    
    Models a population that experienced a reduction in size (bottleneck),
    then recovery. This creates a characteristic signature: reduced diversity,
    excess of rare alleles, and negative Tajima's D.
    
    Args:
        n_sequences: Number of sequences to generate
        sequence_length: Length of each sequence
        pre_bottleneck_diversity: Nucleotide diversity before bottleneck
        bottleneck_size: Effective population size during bottleneck
        bottleneck_duration: Number of generations at bottleneck size
        recovery_generations: Number of generations since bottleneck
        mutation_rate: Per-site mutation rate
        rng: Random number generator
    
    Returns:
        List of sequences reflecting bottleneck signature
    
    Examples:
        >>> seqs = simulate_bottleneck_population(
        ...     n_sequences=20,
        ...     sequence_length=1000,
        ...     bottleneck_size=5,
        ...     bottleneck_duration=10
        ... )
        >>> len(seqs)
        20
    """
    r = rng or random
    
    # Generate ancestral population with high diversity
    ancestral = generate_population_sequences(
        bottleneck_size * 2,  # Larger ancestral population
        sequence_length,
        nucleotide_diversity=pre_bottleneck_diversity,
        mutation_rate=mutation_rate,
        rng=r,
    )
    
    # Bottleneck: sample only bottleneck_size sequences
    bottleneck_samples = r.sample(ancestral, min(bottleneck_size, len(ancestral)))
    
    # During bottleneck: introduce mutations in small population
    bottleneck_sequences = []
    for seq in bottleneck_samples:
        # High mutation rate during bottleneck (more generations)
        mutations = int(mutation_rate * sequence_length * bottleneck_duration)
        mutated = mutate_sequence(seq, n_mut=mutations, rng=r)
        bottleneck_sequences.append(mutated)
    
    # Recovery: expand from bottleneck
    # Create new sequences by mutating from bottleneck survivors
    final_sequences = list(bottleneck_sequences)
    
    # Expand to n_sequences
    while len(final_sequences) < n_sequences:
        parent = r.choice(final_sequences)
        mutations = int(mutation_rate * sequence_length * recovery_generations / n_sequences)
        new_seq = mutate_sequence(parent, n_mut=mutations, rng=r)
        final_sequences.append(new_seq)
    
    # Ensure we have exactly n_sequences
    return final_sequences[:n_sequences]


def simulate_population_expansion(
    n_sequences: int,
    sequence_length: int,
    *,
    initial_diversity: float = 0.005,
    expansion_factor: float = 10.0,
    growth_rate: float = 0.1,
    mutation_rate: float = 0.001,
    rng: random.Random | None = None,
) -> list[str]:
    """Simulate a population that underwent recent expansion.
    
    Models exponential population growth, which creates excess of rare alleles
    and negative Tajima's D (similar to bottleneck but different mechanism).
    
    Args:
        n_sequences: Number of sequences to generate
        sequence_length: Length of each sequence
        initial_diversity: Nucleotide diversity in small ancestral population
        expansion_factor: How much population grew (e.g., 10x)
        growth_rate: Per-generation growth rate
        mutation_rate: Per-site mutation rate
        rng: Random number generator
    
    Returns:
        List of sequences reflecting expansion signature
    
    Examples:
        >>> seqs = simulate_population_expansion(
        ...     n_sequences=20,
        ...     sequence_length=1000,
        ...     expansion_factor=10.0
        ... )
        >>> len(seqs)
        20
    """
    r = rng or random
    
    # Small ancestral population
    initial_size = max(2, int(n_sequences / expansion_factor))
    
    ancestral = generate_population_sequences(
        initial_size,
        sequence_length,
        nucleotide_diversity=initial_diversity,
        mutation_rate=mutation_rate,
        rng=r,
    )
    
    # Simulate expansion: more individuals, more recent mutations
    # Recent expansion means many rare (young) alleles
    
    # Calculate generations since expansion started
    generations = int(math.log(expansion_factor) / growth_rate) if growth_rate > 0 else 10
    
    # Expand population by adding new sequences derived from ancestral
    expanded = list(ancestral)
    
    while len(expanded) < n_sequences:
        # Choose parent from existing sequences
        parent = r.choice(expanded)
        
        # Introduce mutations (recent expansion = many new mutations)
        mutations = int(mutation_rate * sequence_length * generations / n_sequences)
        new_seq = mutate_sequence(parent, n_mut=mutations, rng=r)
        expanded.append(new_seq)
    
    return expanded[:n_sequences]


def generate_site_frequency_spectrum(
    sample_size: int,
    n_sites: int,
    *,
    theta: float = 0.01,
    folded: bool = True,
    rng: random.Random | None = None,
) -> list[int]:
    """Generate a site frequency spectrum with specified properties.
    
    Creates a site frequency spectrum (SFS) under the standard neutral model
    or with specified theta. The SFS describes the distribution of allele
    frequencies across polymorphic sites.
    
    Args:
        sample_size: Number of sampled sequences (n)
        n_sites: Number of polymorphic sites to generate
        theta: Population mutation parameter θ = 4Neμ (for diploid)
        folded: If True, return folded SFS (minor allele frequencies only)
        rng: Random number generator
    
    Returns:
        List of counts per frequency bin. For folded SFS, length is n//2.
        For unfolded SFS, length is n-1.
    
    Examples:
        >>> sfs = generate_site_frequency_spectrum(
        ...     sample_size=10,
        ...     n_sites=100,
        ...     theta=0.01
        ... )
        >>> len(sfs)
        5  # Folded SFS for n=10
        >>> sum(sfs)
        100  # Total polymorphic sites
    """
    r = rng or random
    
    if folded:
        bins = sample_size // 2
    else:
        bins = sample_size - 1
    
    sfs = [0] * bins
    
    # Under standard neutral model, expected frequency of i derived alleles is
    # proportional to 1/i for unfolded SFS
    # For folded: proportional to 1/i + 1/(n-i) for i < n/2
    
    # Generate sites with frequencies following neutral model
    for _ in range(n_sites):
        if folded:
            # Folded: sample minor allele frequency
            # Frequency distribution: higher weight for rare alleles
            freq = r.randint(1, bins)
            # Weight by inverse frequency (more rare alleles)
            if r.random() < 1.0 / freq:
                sfs[freq - 1] += 1
            else:
                # Still add to a bin, but with probability
                sfs[freq - 1] += 1
        else:
            # Unfolded: sample derived allele count
            # Under neutral model: E[SFS[i]] ∝ 1/i
            freq = r.randint(1, bins + 1)
            # Weight by inverse frequency
            if r.random() < 1.0 / freq:
                sfs[freq - 1] += 1
            else:
                sfs[freq - 1] += 1
    
    return sfs


def generate_linkage_disequilibrium_data(
    n_individuals: int,
    n_sites: int,
    *,
    r_squared_target: float = 0.5,
    recombination_rate: float = 0.01,
    allele_frequencies: Sequence[float] | None = None,
    rng: random.Random | None = None,
) -> list[list[int]]:
    """Generate genotype data with specified linkage disequilibrium.
    
    Creates genotypes with linkage disequilibrium (LD) between nearby sites.
    LD is controlled by recombination rate: lower recombination = higher LD.
    
    Args:
        n_individuals: Number of individuals
        n_sites: Number of sites
        r_squared_target: Target r² value (LD measure)
        recombination_rate: Recombination rate between sites (c)
        allele_frequencies: Optional allele frequencies per site
        rng: Random number generator
    
    Returns:
        Genotype matrix with specified LD patterns
    
    Examples:
        >>> genotypes = generate_linkage_disequilibrium_data(
        ...     n_individuals=100,
        ...     n_sites=10,
        ...     r_squared_target=0.5,
        ...     recombination_rate=0.01
        ... )
        >>> len(genotypes)
        100
    """
    r = rng or random
    
    # Generate allele frequencies if not provided
    if allele_frequencies is None:
        freqs = [r.uniform(0.1, 0.5) for _ in range(n_sites)]
    else:
        freqs = list(allele_frequencies)
    
    # Generate haplotypes with LD
    # LD decay: r² ≈ e^(-2 * c * distance)
    # For nearby sites, create correlation
    
    haplotypes = []
    for _ in range(n_individuals * 2):  # Two haplotypes per individual
        haplotype = []
        prev_allele = 1 if r.random() < freqs[0] else 0
        
        for site_idx, freq in enumerate(freqs):
            if site_idx == 0:
                allele = prev_allele
            else:
                # Calculate LD based on distance and recombination
                distance = 1.0  # Assume sites are 1 unit apart
                ld_decay = math.exp(-2 * recombination_rate * distance)
                
                # With probability based on LD, keep same allele as previous
                if r.random() < ld_decay * r_squared_target:
                    allele = prev_allele
                else:
                    # Independent sampling
                    allele = 1 if r.random() < freq else 0
            
            haplotype.append(allele)
            prev_allele = allele
        
        haplotypes.append(haplotype)
    
    # Convert haplotypes to diploid genotypes
    genotypes = []
    for i in range(0, len(haplotypes), 2):
        haplotype1 = haplotypes[i]
        haplotype2 = haplotypes[i + 1]
        genotype = [h1 + h2 for h1, h2 in zip(haplotype1, haplotype2)]
        genotypes.append(genotype)
    
    return genotypes


def _estimate_pi(sequences: Sequence[str]) -> float:
    """Estimate nucleotide diversity from sequences."""
    if len(sequences) < 2:
        return 0.0
    
    L = min(len(s) for s in sequences)
    if L == 0:
        return 0.0
    
    n = len(sequences)
    total_diff = 0.0
    num_pairs = 0
    
    for i in range(n):
        for j in range(i + 1, n):
            diff = sum(1 for a, b in zip(sequences[i][:L], sequences[j][:L]) if a != b)
            total_diff += diff / L
            num_pairs += 1
    
    return total_diff / num_pairs if num_pairs > 0 else 0.0

