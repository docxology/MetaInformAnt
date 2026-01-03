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


def generate_population_sequences(n_sequences: int, length: int, *,
                                theta: float = 0.01, n_sites: int = 1000,
                                rng: random.Random | None = None) -> List[str]:
    """Generate a population of DNA sequences with neutral mutations.

    Args:
        n_sequences: Number of sequences to generate
        length: Length of each sequence
        theta: Population mutation parameter (4*N*mu)
        n_sites: Number of polymorphic sites to simulate
        rng: Random number generator for reproducibility

    Returns:
        List of DNA sequences

    Raises:
        ValueError: If parameters are invalid
    """
    validation.validate_range(n_sequences, min_val=2, name="n_sequences")
    validation.validate_range(length, min_val=1, name="length")
    validation.validate_range(theta, min_val=0.0, name="theta")
    validation.validate_range(n_sites, min_val=1, name="n_sites")

    if rng is None:
        rng = random.Random()

    np.random.seed(rng.randint(0, 2**32))

    # Generate ancestral sequence
    ancestral_seq = ''.join(rng.choices("ATCG", k=length))

    # Generate mutation positions
    mutation_positions = rng.sample(range(length), min(n_sites, length))

    sequences = []

    for _ in range(n_sequences):
        seq = list(ancestral_seq)

        # Apply mutations at polymorphic sites
        for pos in mutation_positions:
            if rng.random() < theta / (4 * n_sequences):  # Simplified mutation probability
                # Mutate to different base
                current_base = seq[pos]
                possible_mutations = [b for b in "ATCG" if b != current_base]
                seq[pos] = rng.choice(possible_mutations)

        sequences.append(''.join(seq))

    return sequences


def generate_two_populations(n_pop1: int, n_pop2: int, length: int, *,
                           f_st: float = 0.1, theta: float = 0.01,
                           n_sites: int = 1000, rng: random.Random | None = None) -> Tuple[List[str], List[str]]:
    """Generate two populations with specified differentiation (F_ST).

    Args:
        n_pop1: Size of population 1
        n_pop2: Size of population 2
        length: Sequence length
        f_st: Desired F_ST value (0-1)
        theta: Population mutation parameter
        n_sites: Number of polymorphic sites
        rng: Random number generator for reproducibility

    Returns:
        Tuple of (population1_sequences, population2_sequences)

    Raises:
        ValueError: If parameters are invalid
    """
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
    ancestral_pop = generate_population_sequences(
        max(n_pop1, n_pop2), length, theta=theta, n_sites=n_sites, rng=rng
    )

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
                    pop2[i] = ''.join(seq_list)

    return pop1, pop2


def generate_genotype_matrix(n_individuals: int, n_snps: int, *,
                           maf_min: float = 0.05, maf_max: float = 0.5,
                           hwe_deviation: float = 0.0,
                           rng: random.Random | None = None) -> np.ndarray:
    """Generate a genotype matrix for population genetics analysis.

    Args:
        n_individuals: Number of individuals
        n_snps: Number of SNPs
        maf_min: Minimum minor allele frequency
        maf_max: Maximum minor allele frequency
        hwe_deviation: Deviation from Hardy-Weinberg equilibrium (0 = perfect HWE)
        rng: Random number generator for reproducibility

    Returns:
        Genotype matrix (n_individuals x n_snps) with values 0, 1, 2

    Raises:
        ValueError: If parameters are invalid
    """
    validation.validate_range(n_individuals, min_val=2, name="n_individuals")
    validation.validate_range(n_snps, min_val=1, name="n_snps")
    validation.validate_range(maf_min, min_val=0.0, max_val=0.5, name="maf_min")
    validation.validate_range(maf_max, min_val=maf_min, max_val=0.5, name="maf_max")
    validation.validate_range(hwe_deviation, min_val=0.0, max_val=1.0, name="hwe_deviation")

    if rng is None:
        rng = random.Random()

    np.random.seed(rng.randint(0, 2**32))

    genotype_matrix = np.zeros((n_individuals, n_snps), dtype=int)

    for snp_idx in range(n_snps):
        # Generate allele frequency
        maf = rng.uniform(maf_min, maf_max)

        # Generate genotypes
        for ind_idx in range(n_individuals):
            if hwe_deviation == 0:
                # Perfect HWE
                p_aa = (1 - maf) ** 2
                p_ab = 2 * maf * (1 - maf)
                p_bb = maf ** 2

                genotype = rng.choices([0, 1, 2], weights=[p_aa, p_ab, p_bb])[0]
            else:
                # Deviate from HWE (e.g., due to selection, assortative mating)
                # Simplified deviation: increase homozygosity
                homozygosity_bias = hwe_deviation

                if rng.random() < homozygosity_bias:
                    genotype = rng.choice([0, 2])  # Favor homozygotes
                else:
                    genotype = rng.choices([0, 1, 2], weights=[(1-maf)**2, 2*maf*(1-maf), maf**2])[0]

            genotype_matrix[ind_idx, snp_idx] = genotype

    return genotype_matrix


def simulate_bottleneck_population(initial_size: int, bottleneck_size: int,
                                 final_size: int, generations: int, *,
                                 mutation_rate: float = 1e-8,
                                 rng: random.Random | None = None) -> Dict[str, Any]:
    """Simulate a population bottleneck followed by recovery.

    Args:
        initial_size: Initial population size
        bottleneck_size: Population size during bottleneck
        final_size: Final population size after recovery
        generations: Total generations to simulate
        mutation_rate: Per-base mutation rate
        rng: Random number generator for reproducibility

    Returns:
        Dictionary with simulation results including population sizes and genetic diversity

    Raises:
        ValueError: If parameters are invalid
    """
    validation.validate_range(initial_size, min_val=2, name="initial_size")
    validation.validate_range(bottleneck_size, min_val=2, name="bottleneck_size")
    validation.validate_range(final_size, min_val=2, name="final_size")
    validation.validate_range(generations, min_val=1, name="generations")
    validation.validate_range(mutation_rate, min_val=0.0, name="mutation_rate")

    if rng is None:
        rng = random.Random()

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
    # In a real simulation, this would track actual alleles
    initial_diversity = 1.0
    diversity_trajectory = []

    for size in population_sizes:
        # Diversity decreases during bottleneck due to drift
        diversity = initial_diversity * (size / initial_size)
        diversity = max(diversity, 0.01)  # Minimum diversity
        diversity_trajectory.append(diversity)

    # Simulate mutations
    total_mutations = 0
    mutation_trajectory = []

    for gen, size in enumerate(population_sizes):
        # Mutations proportional to population size and mutation rate
        expected_mutations = size * 1000 * mutation_rate  # Assuming 1000bp genome
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


def simulate_population_expansion(initial_size: int, final_size: int,
                                expansion_time: int, *,
                                mutation_rate: float = 1e-8,
                                rng: random.Random | None = None) -> Dict[str, Any]:
    """Simulate exponential population expansion.

    Args:
        initial_size: Initial population size
        final_size: Final population size
        expansion_time: Generations over which expansion occurs
        mutation_rate: Per-base mutation rate
        rng: Random number generator for reproducibility

    Returns:
        Dictionary with simulation results

    Raises:
        ValueError: If parameters are invalid
    """
    validation.validate_range(initial_size, min_val=2, name="initial_size")
    validation.validate_range(final_size, min_val=initial_size, name="final_size")
    validation.validate_range(expansion_time, min_val=1, name="expansion_time")
    validation.validate_range(mutation_rate, min_val=0.0, name="mutation_rate")

    if rng is None:
        rng = random.Random()

    np.random.seed(rng.randint(0, 2**32))

    # Exponential growth model
    growth_rate = (np.log(final_size) - np.log(initial_size)) / expansion_time

    population_sizes = []
    diversity_trajectory = []
    mutation_trajectory = []

    total_mutations = 0
    current_diversity = 1.0

    for gen in range(expansion_time + 1):
        # Population size at this generation
        size = int(initial_size * np.exp(growth_rate * gen))
        size = min(size, final_size)  # Cap at final size
        population_sizes.append(size)

        # Diversity increases slightly with population size
        diversity_increase = (size - initial_size) / (final_size - initial_size) * 0.1
        current_diversity = min(1.0, current_diversity + diversity_increase)
        diversity_trajectory.append(current_diversity)

        # Mutations
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


def generate_site_frequency_spectrum(n_samples: int, n_sites: int, *,
                                   demographic_model: str = "constant",
                                   parameters: Dict[str, float] | None = None,
                                   rng: random.Random | None = None) -> np.ndarray:
    """Generate a site frequency spectrum under different demographic models.

    Args:
        n_samples: Number of sampled chromosomes
        n_sites: Number of polymorphic sites
        demographic_model: Demographic model ("constant", "expansion", "bottleneck")
        parameters: Model-specific parameters
        rng: Random number generator for reproducibility

    Returns:
        Site frequency spectrum (array of length n_samples-1)

    Raises:
        ValueError: If parameters are invalid
    """
    validation.validate_range(n_samples, min_val=2, name="n_samples")
    validation.validate_range(n_sites, min_val=1, name="n_sites")

    if demographic_model not in ["constant", "expansion", "bottleneck"]:
        raise errors.ValidationError(f"Unknown demographic model: {demographic_model}")

    if rng is None:
        rng = random.Random()

    np.random.seed(rng.randint(0, 2**32))

    if parameters is None:
        parameters = {}

    sfs = np.zeros(n_samples - 1)

    for site in range(n_sites):
        if demographic_model == "constant":
            # Neutral coalescent
            freq = rng.randint(1, n_samples)  # Uniform frequency
        elif demographic_model == "expansion":
            # Favor low frequency variants (recent expansion)
            alpha = parameters.get("alpha", 1.0)
            freq = rng.zipf(alpha)  # Zipf distribution favors low frequencies
            freq = min(freq, n_samples - 1)
        elif demographic_model == "bottleneck":
            # Favor high frequency variants (old bottleneck)
            bottleneck_strength = parameters.get("bottleneck_strength", 0.1)
            if rng.random() < bottleneck_strength:
                freq = rng.randint(n_samples//2, n_samples-1)  # High frequency
            else:
                freq = rng.randint(1, n_samples//2)  # Low frequency

        # Ensure valid frequency
        freq = max(1, min(freq, n_samples - 1))

        sfs[freq - 1] += 1  # SFS is 1-indexed

    return sfs


def generate_linkage_disequilibrium_data(n_individuals: int, n_snps: int, *,
                                       recombination_rate: float = 1e-8,
                                       selection_coefficient: float = 0.0,
                                       rng: random.Random | None = None) -> Tuple[np.ndarray, np.ndarray]:
    """Generate genotype data with linkage disequilibrium patterns.

    Args:
        n_individuals: Number of individuals
        n_snps: Number of SNPs
        recombination_rate: Recombination rate between adjacent SNPs
        selection_coefficient: Selection coefficient (0 = neutral)
        rng: Random number generator for reproducibility

    Returns:
        Tuple of (genotype_matrix, ld_matrix)

    Raises:
        ValueError: If parameters are invalid
    """
    validation.validate_range(n_individuals, min_val=2, name="n_individuals")
    validation.validate_range(n_snps, min_val=2, name="n_snps")
    validation.validate_range(recombination_rate, min_val=0.0, name="recombination_rate")

    if rng is None:
        rng = random.Random()

    np.random.seed(rng.randint(0, 2**32))

    # Generate haplotype blocks
    block_size = max(2, int(1.0 / recombination_rate) if recombination_rate > 0 else n_snps)
    n_blocks = (n_snps + block_size - 1) // block_size

    genotype_matrix = np.zeros((n_individuals, n_snps), dtype=int)

    for block_start in range(0, n_snps, block_size):
        block_end = min(block_start + block_size, n_snps)
        block_size_actual = block_end - block_start

        # Generate haplotypes for this block
        haplotypes = []

        for _ in range(n_individuals):
            # Generate haplotype with LD structure
            haplotype = []

            for pos in range(block_size_actual):
                if pos == 0:
                    # First SNP independent
                    allele = rng.choice([0, 1])
                else:
                    # Subsequent SNPs linked to previous
                    prev_allele = haplotype[-1]

                    if selection_coefficient == 0:
                        # Neutral LD decay
                        linkage_strength = np.exp(-recombination_rate * pos)
                        if rng.random() < linkage_strength:
                            allele = prev_allele  # Maintain LD
                        else:
                            allele = 1 - prev_allele  # Break LD
                    else:
                        # Selection maintains LD
                        if rng.random() < abs(selection_coefficient):
                            allele = prev_allele
                        else:
                            allele = 1 - prev_allele

                haplotype.append(allele)

            haplotypes.append(haplotype)

        # Convert haplotypes to genotypes (simplified)
        for ind in range(n_individuals):
            for pos in range(block_size_actual):
                genotype_matrix[ind, block_start + pos] = haplotypes[ind][pos]

    # Calculate LD matrix (simplified D')
    ld_matrix = np.zeros((n_snps, n_snps))

    for i in range(n_snps):
        for j in range(i + 1, n_snps):
            # Calculate correlation between SNPs
            snp_i = genotype_matrix[:, i]
            snp_j = genotype_matrix[:, j]

            # Simple correlation coefficient as LD measure
            if np.std(snp_i) > 0 and np.std(snp_j) > 0:
                corr = np.corrcoef(snp_i, snp_j)[0, 1]
                ld_matrix[i, j] = abs(corr)
                ld_matrix[j, i] = abs(corr)

    return genotype_matrix, ld_matrix


def simulate_admixture(n_populations: int, population_sizes: List[int],
                     admixture_proportions: np.ndarray, generations: int, *,
                     migration_rate: float = 0.01,
                     rng: random.Random | None = None) -> Dict[str, Any]:
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


def simulate_selection(genotype_matrix: np.ndarray, fitness_effects: np.ndarray,
                      generations: int, *, selection_strength: float = 0.1,
                      rng: random.Random | None = None) -> Dict[str, Any]:
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

            fitness_values *= (1 + selection_strength * fitness_snp)

        # Normalize fitness (Wright-Fisher model)
        fitness_values /= np.mean(fitness_values)

        # Generate next generation
        next_genotypes = np.zeros_like(current_genotypes)

        for snp in range(n_snps):
            # Sample genotypes proportional to fitness
            genotypes = current_genotypes[:, snp]

            # Create weighted sampling
            weights = fitness_values
            selected_indices = rng.choices(
                range(n_individuals),
                weights=weights,
                k=n_individuals
            )

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





