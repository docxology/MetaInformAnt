"""Population genetics simulation utilities for generating synthetic genomic data.

This module provides functions for simulating various population genetic scenarios
including population bottlenecks, expansions, migration, and natural selection.
All functions support reproducible results through random seed control.
"""

from __future__ import annotations

import random
from typing import Dict, List, Optional, Tuple, Any
import numpy as np

from metainformant.core import logging, validation, errors

logger = logging.get_logger(__name__)


def generate_population_sequences(
    n_sequences: int,
    length: int | None = None,
    *,
    theta: float | None = None,
    n_sites: int = 1000,
    rng: random.Random | None = None,
    # Parameter aliases for backward compatibility
    sequence_length: int | None = None,
    random_seed: int | None = None,
    mutation_rate: float | None = None,
    nucleotide_diversity: float | None = None,
    wattersons_theta: float | None = None,
    reference_sequence: str | None = None,
) -> List[str]:
    """Generate a population of DNA sequences with neutral mutations.

    Args:
        n_sequences: Number of sequences to generate
        length: Length of each sequence
        theta: Population mutation parameter (4*N*mu)
        n_sites: Number of polymorphic sites to simulate
        rng: Random number generator for reproducibility
        sequence_length: Alias for length (backward compatibility)
        random_seed: Random seed for reproducibility (creates rng if provided)
        mutation_rate: Per-site mutation rate (alternative to theta)
        nucleotide_diversity: Target nucleotide diversity (pi)
        wattersons_theta: Watterson's theta estimate
        reference_sequence: Reference sequence to use as ancestor

    Returns:
        List of DNA sequences

    Raises:
        ValueError: If parameters are invalid
    """
    # Handle parameter aliases
    if sequence_length is not None and length is None:
        length = sequence_length
    if reference_sequence is not None:
        length = len(reference_sequence)
    if length is None:
        raise ValueError("length (or sequence_length) is required")
    if random_seed is not None and rng is None:
        rng = random.Random(random_seed)

    # Determine theta from various parameters
    if theta is None:
        if mutation_rate is not None:
            theta = 4 * n_sequences * mutation_rate
        elif nucleotide_diversity is not None:
            theta = nucleotide_diversity
        elif wattersons_theta is not None:
            theta = wattersons_theta
        else:
            theta = 0.01  # Default

    validation.validate_range(n_sequences, min_val=2, name="n_sequences")
    validation.validate_range(length, min_val=1, name="length")
    validation.validate_range(theta, min_val=0.0, name="theta")
    validation.validate_range(n_sites, min_val=1, name="n_sites")

    if rng is None:
        rng = random.Random()

    np.random.seed(rng.randint(0, 2**32))

    # Generate or use ancestral sequence
    if reference_sequence is not None:
        ancestral_seq = reference_sequence
    else:
        ancestral_seq = "".join(rng.choices("ATCG", k=length))

    # Generate mutation positions
    mutation_positions = rng.sample(range(length), min(n_sites, length))

    sequences = []

    for _ in range(n_sequences):
        seq = list(ancestral_seq)

        # Apply mutations at polymorphic sites
        mutation_prob = theta / (4 * n_sequences) if theta > 0 else 0.01
        for pos in mutation_positions:
            if rng.random() < mutation_prob:
                # Mutate to different base
                current_base = seq[pos]
                possible_mutations = [b for b in "ATCG" if b != current_base]
                seq[pos] = rng.choice(possible_mutations)

        sequences.append("".join(seq))

    return sequences


def generate_two_populations(
    n_pop1: int,
    n_pop2: int,
    length: int | None = None,
    *,
    f_st: float | None = None,
    fst: float | None = None,
    theta: float = 0.01,
    n_sites: int = 1000,
    rng: random.Random | None = None,
    sequence_length: int | None = None,
) -> Tuple[List[str], List[str]]:
    """Generate two populations with specified differentiation (F_ST).

    Args:
        n_pop1: Size of population 1
        n_pop2: Size of population 2
        length: Sequence length
        f_st: Desired F_ST value (0-1)
        fst: Alias for f_st (lowercase)
        theta: Population mutation parameter
        n_sites: Number of polymorphic sites
        rng: Random number generator for reproducibility
        sequence_length: Alias for length

    Returns:
        Tuple of (population1_sequences, population2_sequences)

    Raises:
        ValueError: If parameters are invalid
    """
    # Handle parameter aliases
    if sequence_length is not None and length is None:
        length = sequence_length
    if fst is not None and f_st is None:
        f_st = fst
    if f_st is None:
        f_st = 0.1
    if length is None:
        raise ValueError("length (or sequence_length) is required")

    validation.validate_range(n_pop1, min_val=2, name="n_pop1")
    validation.validate_range(n_pop2, min_val=2, name="n_pop2")
    validation.validate_range(length, min_val=1, name="length")
    validation.validate_range(f_st, min_val=0.0, max_val=1.0, name="f_st")
    validation.validate_range(theta, min_val=0.0, name="theta")
    validation.validate_range(n_sites, min_val=1, name="n_sites")

    if rng is None:
        rng = random.Random()

    np.random.seed(rng.randint(0, 2**32))

    # Generate ancestral population
    ancestral_pop = generate_population_sequences(max(n_pop1, n_pop2), length, theta=theta, n_sites=n_sites, rng=rng)

    # Split into two populations
    pop1_indices = rng.sample(range(len(ancestral_pop)), n_pop1)
    pop2_indices = rng.sample(range(len(ancestral_pop)), n_pop2)

    pop1 = [ancestral_pop[i] for i in pop1_indices]
    pop2 = [ancestral_pop[i] for i in pop2_indices]

    # Apply additional differentiation to achieve target F_ST
    if f_st > 0:
        # Add population-specific mutations
        n_diff_sites = int(f_st * n_sites)
        diff_positions = rng.sample(range(length), n_diff_sites)

        for pos in diff_positions:
            # Mutate some individuals in pop2
            for i in range(len(pop2)):
                if rng.random() < 0.5:  # 50% chance to mutate
                    seq_list = list(pop2[i])
                    current_base = seq_list[pos]
                    possible_mutations = [b for b in "ATCG" if b != current_base]
                    seq_list[pos] = rng.choice(possible_mutations)
                    pop2[i] = "".join(seq_list)

    return pop1, pop2


def generate_genotype_matrix(
    n_individuals: int,
    n_snps: int | None = None,
    *,
    n_sites: int | None = None,
    maf_min: float = 0.05,
    maf_max: float = 0.5,
    hwe_deviation: float = 0.0,
    hwe: bool = True,
    allele_frequencies: List[float] | None = None,
    ploidy: int = 2,
    rng: random.Random | None = None,
) -> List[List[int]]:
    """Generate a genotype matrix for population genetics analysis.

    Args:
        n_individuals: Number of individuals
        n_snps: Number of SNPs (alias: n_sites)
        n_sites: Number of sites (alias for n_snps)
        maf_min: Minimum minor allele frequency
        maf_max: Maximum minor allele frequency
        hwe_deviation: Deviation from Hardy-Weinberg equilibrium (0 = perfect HWE)
        hwe: Whether to generate genotypes under HWE (True) or not (False)
        allele_frequencies: List of allele frequencies for each site
        ploidy: Ploidy level (1=haploid, 2=diploid)
        rng: Random number generator for reproducibility

    Returns:
        Genotype matrix as list of lists (n_individuals x n_sites) with values 0-ploidy

    Raises:
        ValueError: If parameters are invalid
    """
    # Handle parameter aliases
    if n_sites is not None and n_snps is None:
        n_snps = n_sites
    if n_snps is None:
        raise ValueError("n_snps (or n_sites) is required")

    # Use allele_frequencies to infer n_snps if provided
    if allele_frequencies is not None:
        n_snps = len(allele_frequencies)

    validation.validate_range(n_individuals, min_val=1, name="n_individuals")
    validation.validate_range(n_snps, min_val=1, name="n_snps")
    validation.validate_range(maf_min, min_val=0.0, max_val=0.5, name="maf_min")
    validation.validate_range(maf_max, min_val=maf_min, max_val=0.5, name="maf_max")
    validation.validate_range(hwe_deviation, min_val=0.0, max_val=1.0, name="hwe_deviation")
    validation.validate_range(ploidy, min_val=1, max_val=2, name="ploidy")

    # Convert hwe boolean to hwe_deviation
    if not hwe:
        hwe_deviation = 0.3  # Some deviation from HWE

    if rng is None:
        rng = random.Random()

    np.random.seed(rng.randint(0, 2**32))

    genotype_matrix = []

    for ind_idx in range(n_individuals):
        ind_genotypes = []
        for snp_idx in range(n_snps):
            # Get allele frequency
            if allele_frequencies is not None:
                maf = allele_frequencies[snp_idx]
            else:
                maf = rng.uniform(maf_min, maf_max)

            # Generate genotype based on ploidy
            if ploidy == 1:
                # Haploid: 0 or 1
                genotype = 1 if rng.random() < maf else 0
            else:
                # Diploid: 0, 1, or 2
                if hwe_deviation == 0:
                    # Perfect HWE
                    p_aa = (1 - maf) ** 2
                    p_ab = 2 * maf * (1 - maf)
                    p_bb = maf**2

                    genotype = rng.choices([0, 1, 2], weights=[p_aa, p_ab, p_bb])[0]
                else:
                    # Deviate from HWE
                    if rng.random() < hwe_deviation:
                        genotype = rng.choice([0, 2])  # Favor homozygotes
                    else:
                        p_aa = (1 - maf) ** 2
                        p_ab = 2 * maf * (1 - maf)
                        p_bb = maf**2
                        genotype = rng.choices([0, 1, 2], weights=[p_aa, p_ab, p_bb])[0]

            ind_genotypes.append(genotype)
        genotype_matrix.append(ind_genotypes)

    return genotype_matrix


def simulate_bottleneck_population(
    initial_size: int | None = None,
    bottleneck_size: int | None = None,
    final_size: int | None = None,
    generations: int | None = None,
    *,
    mutation_rate: float = 1e-8,
    rng: random.Random | None = None,
    # Alternative interface for sequence generation
    n_sequences: int | None = None,
    sequence_length: int | None = None,
    bottleneck_duration: int | None = None,
    pre_bottleneck_diversity: float | None = None,
) -> List[str] | Dict[str, Any]:
    """Simulate a population bottleneck followed by recovery.

    Can return either sequences (if n_sequences provided) or statistics dict.

    Args:
        initial_size: Initial population size
        bottleneck_size: Population size during bottleneck
        final_size: Final population size after recovery
        generations: Total generations to simulate
        mutation_rate: Per-base mutation rate
        rng: Random number generator for reproducibility
        n_sequences: Number of sequences to generate (alternative interface)
        sequence_length: Length of each sequence
        bottleneck_duration: Duration of bottleneck in generations
        pre_bottleneck_diversity: Initial nucleotide diversity

    Returns:
        List of sequences if n_sequences provided, otherwise dictionary with simulation results

    Raises:
        ValueError: If parameters are invalid
    """
    if rng is None:
        rng = random.Random()

    # Alternative interface: return sequences
    if n_sequences is not None:
        if sequence_length is None:
            raise ValueError("sequence_length required when using n_sequences interface")
        if bottleneck_size is None:
            bottleneck_size = max(2, n_sequences // 4)
        if bottleneck_duration is None:
            bottleneck_duration = 10

        # Generate pre-bottleneck population
        diversity = pre_bottleneck_diversity if pre_bottleneck_diversity else 0.01
        pre_seqs = generate_population_sequences(
            n_sequences=bottleneck_size,
            sequence_length=sequence_length,
            mutation_rate=diversity,
            rng=rng,
        )

        # Expand from bottleneck survivors
        sequences = []
        for _ in range(n_sequences):
            # Sample from bottleneck survivors
            base_seq = rng.choice(pre_seqs)
            # Add some mutations
            seq_list = list(base_seq)
            for i in range(len(seq_list)):
                if rng.random() < diversity / 10:  # Reduced diversity post-bottleneck
                    current = seq_list[i]
                    seq_list[i] = rng.choice([b for b in "ATCG" if b != current])
            sequences.append("".join(seq_list))

        return sequences

    # Original interface: return statistics
    if initial_size is None:
        initial_size = 1000
    if bottleneck_size is None:
        bottleneck_size = 100
    if final_size is None:
        final_size = initial_size
    if generations is None:
        generations = 100

    validation.validate_range(initial_size, min_val=2, name="initial_size")
    validation.validate_range(bottleneck_size, min_val=2, name="bottleneck_size")
    validation.validate_range(final_size, min_val=2, name="final_size")
    validation.validate_range(generations, min_val=1, name="generations")
    validation.validate_range(mutation_rate, min_val=0.0, name="mutation_rate")

    np.random.seed(rng.randint(0, 2**32))

    # Simulate population size changes
    bottleneck_start = generations // 3
    bottleneck_end = 2 * generations // 3

    population_sizes = []

    for gen in range(generations):
        if gen < bottleneck_start:
            size = initial_size
        elif gen < bottleneck_end:
            size = bottleneck_size
        else:
            size = final_size
        population_sizes.append(size)

    # Simulate genetic diversity changes (simplified)
    initial_diversity = 1.0
    diversity_trajectory = []

    for size in population_sizes:
        diversity = initial_diversity * (size / initial_size)
        diversity = max(diversity, 0.01)
        diversity_trajectory.append(diversity)

    # Simulate mutations
    total_mutations = 0
    mutation_trajectory = []

    for gen, size in enumerate(population_sizes):
        expected_mutations = size * 1000 * mutation_rate
        mutations_this_gen = rng.poisson(expected_mutations)
        total_mutations += mutations_this_gen
        mutation_trajectory.append(total_mutations)

    return {
        "population_sizes": population_sizes,
        "diversity_trajectory": diversity_trajectory,
        "mutation_trajectory": mutation_trajectory,
        "bottleneck_start": bottleneck_start,
        "bottleneck_end": bottleneck_end,
        "total_mutations": total_mutations,
        "final_diversity": diversity_trajectory[-1],
    }


def simulate_population_expansion(
    initial_size: int | None = None,
    final_size: int | None = None,
    expansion_time: int | None = None,
    *,
    mutation_rate: float = 1e-8,
    rng: random.Random | None = None,
    # Alternative interface for sequence generation
    n_sequences: int | None = None,
    sequence_length: int | None = None,
    expansion_factor: float | None = None,
) -> List[str] | Dict[str, Any]:
    """Simulate exponential population expansion.

    Can return either sequences (if n_sequences provided) or statistics dict.

    Args:
        initial_size: Initial population size
        final_size: Final population size
        expansion_time: Generations over which expansion occurs
        mutation_rate: Per-base mutation rate
        rng: Random number generator for reproducibility
        n_sequences: Number of sequences to generate (alternative interface)
        sequence_length: Length of each sequence
        expansion_factor: Factor by which population expanded

    Returns:
        List of sequences if n_sequences provided, otherwise dictionary with simulation results

    Raises:
        ValueError: If parameters are invalid
    """
    if rng is None:
        rng = random.Random()

    # Alternative interface: return sequences
    if n_sequences is not None:
        if sequence_length is None:
            raise ValueError("sequence_length required when using n_sequences interface")
        if expansion_factor is None:
            expansion_factor = 10.0

        # Start with small founding population
        founding_size = max(2, int(n_sequences / expansion_factor))

        # Generate founding sequences
        founding_seqs = generate_population_sequences(
            n_sequences=founding_size,
            sequence_length=sequence_length,
            mutation_rate=0.01,
            rng=rng,
        )

        # Expand population
        sequences = []
        for _ in range(n_sequences):
            base_seq = rng.choice(founding_seqs)
            # Add some mutations during expansion
            seq_list = list(base_seq)
            for i in range(len(seq_list)):
                if rng.random() < 0.001:  # Some new mutations
                    current = seq_list[i]
                    seq_list[i] = rng.choice([b for b in "ATCG" if b != current])
            sequences.append("".join(seq_list))

        return sequences

    # Original interface: return statistics
    if initial_size is None:
        initial_size = 100
    if final_size is None:
        final_size = 10000
    if expansion_time is None:
        expansion_time = 100

    validation.validate_range(initial_size, min_val=2, name="initial_size")
    validation.validate_range(final_size, min_val=initial_size, name="final_size")
    validation.validate_range(expansion_time, min_val=1, name="expansion_time")
    validation.validate_range(mutation_rate, min_val=0.0, name="mutation_rate")

    np.random.seed(rng.randint(0, 2**32))

    # Exponential growth model
    growth_rate = (np.log(final_size) - np.log(initial_size)) / expansion_time

    population_sizes = []
    diversity_trajectory = []
    mutation_trajectory = []

    total_mutations = 0
    current_diversity = 1.0

    for gen in range(expansion_time + 1):
        size = int(initial_size * np.exp(growth_rate * gen))
        size = min(size, final_size)
        population_sizes.append(size)

        diversity_increase = (size - initial_size) / (final_size - initial_size) * 0.1
        current_diversity = min(1.0, current_diversity + diversity_increase)
        diversity_trajectory.append(current_diversity)

        expected_mutations = size * 1000 * mutation_rate
        mutations_this_gen = rng.poisson(expected_mutations)
        total_mutations += mutations_this_gen
        mutation_trajectory.append(total_mutations)

    return {
        "population_sizes": population_sizes,
        "diversity_trajectory": diversity_trajectory,
        "mutation_trajectory": mutation_trajectory,
        "growth_rate": growth_rate,
        "expansion_time": expansion_time,
        "total_mutations": total_mutations,
        "final_diversity": diversity_trajectory[-1],
    }


def generate_site_frequency_spectrum(
    n_samples: int | None = None,
    n_sites: int | None = None,
    *,
    sample_size: int | None = None,
    folded: bool = False,
    demographic_model: str = "constant",
    parameters: Dict[str, float] | None = None,
    rng: random.Random | None = None,
) -> List[int]:
    """Generate a site frequency spectrum under different demographic models.

    Args:
        n_samples: Number of sampled chromosomes (alias: sample_size)
        n_sites: Number of polymorphic sites
        sample_size: Alias for n_samples
        folded: Whether to return folded SFS (default False)
        demographic_model: Demographic model ("constant", "expansion", "bottleneck")
        parameters: Model-specific parameters
        rng: Random number generator for reproducibility

    Returns:
        Site frequency spectrum as list of counts

    Raises:
        ValueError: If parameters are invalid
    """
    # Handle parameter aliases
    if sample_size is not None and n_samples is None:
        n_samples = sample_size
    if n_samples is None:
        raise ValueError("n_samples (or sample_size) is required")
    if n_sites is None:
        raise ValueError("n_sites is required")

    validation.validate_range(n_samples, min_val=2, name="n_samples")
    validation.validate_range(n_sites, min_val=1, name="n_sites")

    if demographic_model not in ["constant", "expansion", "bottleneck"]:
        raise errors.ValidationError(f"Unknown demographic model: {demographic_model}")

    if rng is None:
        rng = random.Random()

    np.random.seed(rng.randint(0, 2**32))

    if parameters is None:
        parameters = {}

    # Unfolded SFS has n-1 bins
    sfs = [0] * (n_samples - 1)

    for site in range(n_sites):
        if demographic_model == "constant":
            # Neutral coalescent - 1/i distribution
            weights = [1.0 / i for i in range(1, n_samples)]
            total = sum(weights)
            weights = [w / total for w in weights]
            freq = rng.choices(range(1, n_samples), weights=weights)[0]
        elif demographic_model == "expansion":
            # Favor low frequency variants
            alpha = parameters.get("alpha", 1.0)
            freq = min(rng.zipf(alpha), n_samples - 1)
        elif demographic_model == "bottleneck":
            # Favor high frequency variants
            bottleneck_strength = parameters.get("bottleneck_strength", 0.1)
            if rng.random() < bottleneck_strength:
                freq = rng.randint(n_samples // 2, n_samples - 1)
            else:
                freq = rng.randint(1, max(2, n_samples // 2))

        # Ensure valid frequency
        freq = max(1, min(freq, n_samples - 1))

        sfs[freq - 1] += 1

    # Fold SFS if requested
    if folded:
        folded_sfs = []
        mid = n_samples // 2
        for i in range(mid):
            # Combine frequency class i and n-1-i (complement frequencies)
            j = n_samples - 2 - i
            if i < len(sfs) and j < len(sfs) and i != j:
                folded_sfs.append(sfs[i] + sfs[j])
            elif i < len(sfs):
                # Middle bin (only when i == j) or edge case
                folded_sfs.append(sfs[i])
        return folded_sfs

    return sfs


def generate_linkage_disequilibrium_data(
    n_individuals: int,
    n_snps: int | None = None,
    *,
    n_sites: int | None = None,
    recombination_rate: float = 1e-8,
    selection_coefficient: float = 0.0,
    r_squared_target: float | None = None,
    allele_frequencies: List[float] | None = None,
    rng: random.Random | None = None,
) -> List[List[int]]:
    """Generate genotype data with linkage disequilibrium patterns.

    Args:
        n_individuals: Number of individuals
        n_snps: Number of SNPs (alias: n_sites)
        n_sites: Alias for n_snps
        recombination_rate: Recombination rate between adjacent SNPs
        selection_coefficient: Selection coefficient (0 = neutral)
        r_squared_target: Target r-squared for LD between adjacent sites
        allele_frequencies: List of allele frequencies per site
        rng: Random number generator for reproducibility

    Returns:
        Genotype matrix as list of lists

    Raises:
        ValueError: If parameters are invalid
    """
    # Handle parameter aliases
    if n_sites is not None and n_snps is None:
        n_snps = n_sites
    if allele_frequencies is not None:
        n_snps = len(allele_frequencies)
    if n_snps is None:
        raise ValueError("n_snps (or n_sites) is required")

    validation.validate_range(n_individuals, min_val=2, name="n_individuals")
    validation.validate_range(n_snps, min_val=2, name="n_snps")
    validation.validate_range(recombination_rate, min_val=0.0, name="recombination_rate")

    if rng is None:
        rng = random.Random()

    np.random.seed(rng.randint(0, 2**32))

    # Determine linkage strength from r_squared_target if provided
    linkage_strength = 0.5  # default
    if r_squared_target is not None:
        # r_squared_target controls how correlated adjacent SNPs are
        linkage_strength = r_squared_target

    genotype_matrix = np.zeros((n_individuals, n_snps), dtype=int)

    for ind in range(n_individuals):
        for pos in range(n_snps):
            if allele_frequencies is not None:
                # Use provided allele frequency
                af = allele_frequencies[pos]
            else:
                # Use default 0.5
                af = 0.5

            if pos == 0:
                # First SNP based on allele frequency
                allele = 1 if rng.random() < af else 0
            else:
                # Subsequent SNPs linked to previous with strength based on r_squared_target
                prev_allele = genotype_matrix[ind, pos - 1]

                if selection_coefficient == 0:
                    # Neutral LD - use linkage_strength
                    if rng.random() < linkage_strength:
                        allele = prev_allele  # Maintain LD
                    else:
                        # Generate based on allele frequency
                        allele = 1 if rng.random() < af else 0
                else:
                    # Selection maintains LD
                    if rng.random() < abs(selection_coefficient):
                        allele = prev_allele
                    else:
                        allele = 1 if rng.random() < af else 0

            genotype_matrix[ind, pos] = allele

    # Return as List[List[int]]
    return [[int(genotype_matrix[i, j]) for j in range(n_snps)] for i in range(n_individuals)]


def simulate_admixture(
    n_populations: int,
    population_sizes: List[int],
    admixture_proportions: np.ndarray,
    generations: int,
    *,
    migration_rate: float = 0.01,
    rng: random.Random | None = None,
) -> Dict[str, Any]:
    """Simulate population admixture.

    Args:
        n_populations: Number of ancestral populations
        population_sizes: Size of each ancestral population
        admixture_proportions: Mixing proportions (n_populations x n_populations)
        generations: Generations since admixture
        migration_rate: Migration rate between populations
        rng: Random number generator for reproducibility

    Returns:
        Dictionary with admixture simulation results

    Raises:
        ValueError: If parameters are invalid
    """
    validation.validate_range(n_populations, min_val=2, name="n_populations")
    validation.validate_type(population_sizes, list, "population_sizes")
    validation.validate_range(generations, min_val=0, name="generations")
    validation.validate_range(migration_rate, min_val=0.0, max_val=1.0, name="migration_rate")

    if len(population_sizes) != n_populations:
        raise errors.ValidationError(f"population_sizes must have length {n_populations}")

    if admixture_proportions.shape != (n_populations, n_populations):
        raise errors.ValidationError(f"admixture_proportions must have shape {(n_populations, n_populations)}")

    if rng is None:
        rng = random.Random()

    np.random.seed(rng.randint(0, 2**32))

    # Initialize ancestral populations
    ancestral_frequencies = []
    for i in range(n_populations):
        freq = rng.uniform(0.1, 0.9)  # Allele frequency in ancestral population
        ancestral_frequencies.append(freq)

    # Simulate admixture over generations
    current_frequencies = np.array(ancestral_frequencies)
    frequency_trajectory = [current_frequencies.copy()]

    for gen in range(generations):
        # Apply migration/admixture
        new_frequencies = admixture_proportions @ current_frequencies

        # Add some drift
        for i in range(n_populations):
            drift_effect = rng.normal(0, 0.01)
            new_frequencies[i] += drift_effect
            new_frequencies[i] = np.clip(new_frequencies[i], 0.0, 1.0)

        current_frequencies = new_frequencies
        frequency_trajectory.append(current_frequencies.copy())

    # Calculate admixture proportions over time
    admixture_trajectory = []
    for freqs in frequency_trajectory:
        admixture_trajectory.append(freqs / np.sum(freqs))

    return {
        "ancestral_frequencies": ancestral_frequencies,
        "frequency_trajectory": frequency_trajectory,
        "admixture_trajectory": admixture_trajectory,
        "final_frequencies": current_frequencies,
        "generations": generations,
        "migration_rate": migration_rate,
    }


def simulate_selection(
    genotype_matrix: np.ndarray,
    fitness_effects: np.ndarray,
    generations: int,
    *,
    selection_strength: float = 0.1,
    rng: random.Random | None = None,
) -> Dict[str, Any]:
    """Simulate natural selection on genotypes.

    Args:
        genotype_matrix: Initial genotype matrix (individuals x SNPs)
        fitness_effects: Fitness effect for each genotype (3 values per SNP: AA, AB, BB)
        generations: Generations to simulate
        selection_strength: Overall strength of selection
        rng: Random number generator for reproducibility

    Returns:
        Dictionary with selection simulation results

    Raises:
        ValueError: If parameters are invalid
    """
    validation.validate_type(genotype_matrix, np.ndarray, "genotype_matrix")
    validation.validate_type(fitness_effects, np.ndarray, "fitness_effects")
    validation.validate_range(generations, min_val=1, name="generations")
    validation.validate_range(selection_strength, min_val=0.0, name="selection_strength")

    n_individuals, n_snps = genotype_matrix.shape

    if fitness_effects.shape != (n_snps, 3):
        raise errors.ValidationError(f"fitness_effects must have shape ({n_snps}, 3)")

    if rng is None:
        rng = random.Random()

    np.random.seed(rng.randint(0, 2**32))

    # Track allele frequencies over time
    allele_frequencies = []
    current_genotypes = genotype_matrix.copy()

    for gen in range(generations + 1):
        # Calculate current allele frequencies
        freqs = []
        for snp in range(n_snps):
            genotypes = current_genotypes[:, snp]
            alt_allele_count = np.sum(genotypes)  # 0=AA, 1=AB, 2=BB
            freq = alt_allele_count / (2 * n_individuals)
            freqs.append(freq)
        allele_frequencies.append(freqs)

        if gen == generations:
            break

        # Calculate fitness for each individual
        fitness_values = np.ones(n_individuals)

        for snp in range(n_snps):
            genotypes = current_genotypes[:, snp]
            fitness_snp = np.zeros(n_individuals)

            for i, genotype in enumerate(genotypes):
                fitness_snp[i] = fitness_effects[snp, genotype]

            fitness_values *= 1 + selection_strength * fitness_snp

        # Normalize fitness (Wright-Fisher model)
        fitness_values /= np.mean(fitness_values)

        # Generate next generation
        next_genotypes = np.zeros_like(current_genotypes)

        for snp in range(n_snps):
            # Sample genotypes proportional to fitness
            genotypes = current_genotypes[:, snp]

            # Create weighted sampling
            weights = fitness_values
            selected_indices = rng.choices(range(n_individuals), weights=weights, k=n_individuals)

            next_genotypes[:, snp] = genotypes[selected_indices]

        current_genotypes = next_genotypes

    return {
        "initial_genotypes": genotype_matrix,
        "final_genotypes": current_genotypes,
        "allele_frequencies": allele_frequencies,
        "fitness_effects": fitness_effects,
        "generations": generations,
        "selection_strength": selection_strength,
    }
