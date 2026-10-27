# Simulation: Population Genetics

Synthetic data generation for population genetics analysis.

**Functions**: `generate_population_sequences`, `generate_two_populations`, `generate_genotype_matrix`, `simulate_bottleneck_population`, `simulate_population_expansion`, `generate_site_frequency_spectrum`, `generate_linkage_disequilibrium_data`.

## Overview

This module provides comprehensive synthetic data generation for population genetics analysis. You can generate:

- **Sequences with a chosen mutation level** (a per-site probability that sets π and θ_W)
- **Two populations with a target Fst** (Balding-Nichols allele frequencies)
- **Genotype matrices** with Hardy-Weinberg equilibrium
- **Demographic scenarios** (bottlenecks, expansions)
- **Site frequency spectra** (constant, expansion and bottleneck shapes)
- **Linkage disequilibrium** data

## Functions

### `generate_population_sequences(n_sequences, length=None, *, ...)`

Generate a population of sequences by mutating one ancestral sequence independently in each sequence.

**Parameters**:
- `n_sequences`: Number of sequences to generate (at least 2)
- `length`: Length of each sequence (alias `sequence_length`); required unless `reference_sequence` is given
- `n_sites`: Number of candidate polymorphic positions, chosen at random (default 1000, capped at the sequence length); only these positions can vary
- `nucleotide_diversity`: Per-sequence probability that a candidate position differs from the ancestral base. Despite the name it is used as this mutation probability, not as a target π (see the note below)
- `wattersons_theta`: Used the same way as `nucleotide_diversity` when that is not given
- `mutation_rate`: Per-site mutation rate; it is converted to θ = 4 · `n_sequences` · `mutation_rate`, which gives the same probability
- `theta`: Population mutation parameter; the probability is θ / (4 · `n_sequences`). The default θ is 0.01
- `reference_sequence`: Ancestral sequence. Without one, a random A/C/G/T sequence (about 50% GC) is drawn; a reference replaces `length` with its own length
- `rng`: `random.Random` for reproducibility (or `random_seed`, an integer)

**Returns**: List of `n_sequences` strings of equal length

**Note**: Two sequences differ at a position when exactly one of them changed, so for probability *p* the expected π is about 2*p*(1 − *p*) × `min(n_sites, length) / length`. With the default `n_sites` and `length=1000`, `nucleotide_diversity=0.01` gives π ≈ 0.02.

**Example**:
```python-snippet
from metainformant.simulation.models.popgen import generate_population_sequences

# Each sequence differs from the ancestor at about 1% of the 1000 positions
sequences = generate_population_sequences(
    n_sequences=10,
    sequence_length=1000,
    nucleotide_diversity=0.01
)

# Analyze the generated data
from metainformant.dna.population import nucleotide_diversity
pi = nucleotide_diversity(sequences)
print(f"Observed π: {pi:.4f}")  # about 0.02, not 0.01 (see the note above)
```

### `generate_two_populations(n_pop1, n_pop2, length=None, *, f_st=None, fst=None, ...)`

Generate two populations differentiated by a target Fst under the Balding-Nichols model.

Each polymorphic site gets an ancestral derived-allele frequency *p* drawn uniformly from [0.05, 0.95]. Each population's frequency is then drawn from Beta(*p*(1 − F)/F, (1 − *p*)(1 − F)/F), so that Hudson's Fst equals the target F in expectation. Every sequence carries the derived allele at each site independently with its population's frequency, and every other position keeps the ancestral base. `fst=0` gives identical frequencies and `fst=1` fixes each population for one allele.

**Parameters**:
- `n_pop1`: Number of sequences in population 1 (at least 2)
- `n_pop2`: Number of sequences in population 2 (at least 2)
- `length`: Length of each sequence (alias `sequence_length`); required
- `fst`: Target Fst between 0 and 1 (alias `f_st`; default 0.1)
- `theta`: Per-site population mutation rate; it sets the default number of polymorphic sites (default 0.01)
- `n_sites`: Number of polymorphic sites, at most `length`. The default is Watterson's expectation `theta * length * a_n`, where `a_n` is the harmonic sum over the `n_pop1 + n_pop2` sequences (at least 1)
- `rng`: `random.Random` for reproducibility

**Returns**: Tuple of (pop1_sequences, pop2_sequences)

Hudson's estimate is noisy when few sites are polymorphic. The default gives 42 polymorphic sites for two populations of 20 sequences at length 1000, and the estimate then varies between runs with a standard deviation of about 0.04 at Fst = 0.2; raise `n_sites` to tighten it.

**Example**:
```python-snippet
from metainformant.simulation.models.popgen import generate_two_populations
from metainformant.dna.population import hudson_fst

# Generate two populations with Fst = 0.2
pop1, pop2 = generate_two_populations(
    n_pop1=10,
    n_pop2=10,
    sequence_length=1000,
    fst=0.2,
    n_sites=500
)

# Verify Fst
fst_observed = hudson_fst(pop1, pop2)
print(f"Observed Fst: {fst_observed:.4f}")  # about 0.2, varying by roughly ±0.02 between runs
```

### `generate_genotype_matrix(n_individuals, n_snps=None, *, ...)`

Generate a genotype matrix with specified allele frequencies.

**Parameters**:
- `n_individuals`: Number of individuals (rows)
- `n_snps`: Number of sites (columns); alias `n_sites`. Taken from the length of `allele_frequencies` when that is given
- `allele_frequencies`: Optional list of alternate-allele frequencies per site
- `maf_min`: Minimum minor allele frequency (0 to 0.5, default 0.05); sites without a given frequency draw theirs uniformly between `maf_min` and `maf_max`
- `maf_max`: Maximum minor allele frequency (default 0.5)
- `hwe`: If True (default), genotypes follow Hardy-Weinberg equilibrium. If False, 30% of genotypes are replaced by a random homozygote (this overrides `hwe_deviation`)
- `hwe_deviation`: Probability between 0 and 1 that a genotype is replaced by a random homozygote (default 0)
- `ploidy`: Ploidy level (2 for diploid, 1 for haploid)
- `rng`: `random.Random` for reproducibility

**Returns**: Genotype matrix as a list of lists with one row per individual and one column per site. For diploid: 0=AA, 1=AB, 2=BB (copies of the alternate allele); for haploid: 0 or 1.

**Example**:
```python-snippet
from metainformant.simulation.models.popgen import generate_genotype_matrix
from metainformant.dna.population import allele_frequencies

# Generate genotypes with specified frequencies
genotypes = generate_genotype_matrix(
    n_individuals=1000,
    n_sites=10,
    allele_frequencies=[0.2, 0.3, 0.4, 0.1, 0.5, 0.3, 0.2, 0.4, 0.1, 0.6],
    hwe=True
)

# Verify frequencies. The matrix has one row per individual, but
# allele_frequencies expects one row per locus, so transpose it first.
freqs = allele_frequencies(list(zip(*genotypes)))
print(f"Allele frequencies: {[round(f, 2) for f in freqs]}")
```

### `simulate_bottleneck_population(n_sequences=..., sequence_length=..., ...)`

Simulate a population that went through a bottleneck. Pass `n_sequences` and `sequence_length` as keyword arguments to get sequences; the function's first four positional parameters belong to the trajectory form described below.

**Parameters**:
- `n_sequences`: Number of sequences to generate
- `sequence_length`: Length of each sequence (required with `n_sequences`)
- `pre_bottleneck_diversity`: Per-site mutation probability of the founder sequences (default 0.01); descendants gain new mutations at one tenth of it
- `bottleneck_size`: Number of founder sequences that survive the bottleneck (default `max(2, n_sequences // 4)`)
- `bottleneck_duration`: Accepted for compatibility; the sequence form ignores it
- `rng`: `random.Random` for reproducibility

**Returns**: List of sequences: each one is a randomly chosen founder with a few new mutations

**Trajectory form**: called without `n_sequences`, the function takes `initial_size`, `bottleneck_size`, `final_size`, `generations` and `mutation_rate` instead and returns a dictionary with `population_sizes`, `diversity_trajectory`, `mutation_trajectory`, `bottleneck_start`, `bottleneck_end`, `total_mutations` and `final_diversity`. The bottleneck occupies the middle third of `generations`.

**Example**:
```python-snippet
from metainformant.simulation.models.popgen import simulate_bottleneck_population
from metainformant.dna.population import nucleotide_diversity, tajimas_d

# 20 sequences resampled from 5 bottleneck survivors
sequences = simulate_bottleneck_population(
    n_sequences=20,
    sequence_length=1000,
    pre_bottleneck_diversity=0.01,
    bottleneck_size=5
)

# Survivors are resampled with few new mutations; this heuristic does not
# simulate a coalescent, so the sign of D is not guaranteed.
pi = nucleotide_diversity(sequences)
d = tajimas_d(sequences)
print(f"π: {pi:.4f}, Tajima's D: {d:.4f}")
```

### `simulate_population_expansion(n_sequences=..., sequence_length=..., ...)`

Simulate a population that underwent recent expansion. Pass `n_sequences` and `sequence_length` as keyword arguments to get sequences; the function's first three positional parameters belong to the trajectory form described below.

**Parameters**:
- `n_sequences`: Number of sequences to generate
- `sequence_length`: Length of each sequence (required with `n_sequences`)
- `expansion_factor`: How much the population grew (default 10.0); the founding population has `max(2, int(n_sequences / expansion_factor))` sequences, drawn with a per-site mutation probability of 0.01
- `rng`: `random.Random` for reproducibility

**Returns**: List of sequences: each one is a randomly chosen founder with new mutations at probability 0.001 per site

**Trajectory form**: called without `n_sequences`, the function takes `initial_size`, `final_size`, `expansion_time` and `mutation_rate` instead and returns a dictionary with `population_sizes`, `diversity_trajectory`, `mutation_trajectory`, `growth_rate`, `expansion_time`, `total_mutations` and `final_diversity`. The `growth_rate` is the per-generation exponential rate, (ln `final_size` − ln `initial_size`) / `expansion_time`, so it is an output, not an input.

**Example**:
```python-snippet
from metainformant.simulation.models.popgen import simulate_population_expansion

# Simulate 10x expansion
sequences = simulate_population_expansion(
    n_sequences=20,
    sequence_length=1000,
    expansion_factor=10.0
)

# Descendants share founder haplotypes, which tends to give positive D here,
# unlike the excess of rare variants a real expansion leaves.
from metainformant.dna.population import tajimas_d
d = tajimas_d(sequences)
print(f"Tajima's D: {d:.4f}")
```

### `generate_site_frequency_spectrum(n_samples=None, n_sites=None, *, ...)`

Generate a site frequency spectrum with specified properties.

**Parameters**:
- `n_samples`: Number of sampled sequences (alias `sample_size`; at least 2)
- `n_sites`: Number of polymorphic sites to generate (required)
- `folded`: If True, return the folded SFS (minor allele counts only)
- `demographic_model`: `"constant"` (default; counts follow the neutral 1/i expectation), `"expansion"` or `"bottleneck"`; anything else raises `ValidationError`
- `parameters`: Optional dictionary of model parameters: `alpha` (default 1.0) for `"expansion"`, a Zipf exponent that only takes effect above 1, and `bottleneck_strength` (default 0.1) for `"bottleneck"`
- `rng`: `random.Random` for reproducibility

**Returns**: List of site counts per frequency class. The unfolded spectrum has `n_samples - 1` bins (bin *i* holds sites where *i* + 1 samples carry the derived allele); the folded spectrum has `n_samples // 2` bins (bin *i* holds sites with minor allele count *i* + 1).

**Example**:
```python-snippet
from metainformant.simulation.models.popgen import generate_site_frequency_spectrum

# Generate folded SFS for 10 samples, 100 sites
sfs = generate_site_frequency_spectrum(
    sample_size=10,
    n_sites=100,
    folded=True
)

# The constant model has more rare alleles (first bins)
print(f"SFS: {sfs}")  # sums to 100, with the largest counts in the first bins
```

### `generate_linkage_disequilibrium_data(n_individuals, n_snps=None, *, ...)`

Generate genotype data in which neighbouring sites are correlated.

Sites are generated left to right: each site copies the allele of the previous site with a fixed probability and is otherwise drawn from its own allele frequency.

**Parameters**:
- `n_individuals`: Number of individuals (at least 2)
- `n_snps`: Number of sites (at least 2); alias `n_sites`. Taken from the length of `allele_frequencies` when that is given
- `r_squared_target`: The copy probability (default 0.5). Despite the name it is not r² itself: with equal allele frequencies the correlation r between adjacent sites equals it, so the observed r² is about its square
- `selection_coefficient`: Default 0; when non-zero, its absolute value replaces the copy probability
- `recombination_rate`: Accepted and validated, but it does not change the output
- `allele_frequencies`: Optional frequency of allele 1 per site (default 0.5 at every site)
- `rng`: `random.Random` for reproducibility

**Returns**: List of lists with one row per individual and one column per site, holding 0/1 alleles (not 0/1/2 genotypes)

**Example**:
```python-snippet
from metainformant.simulation.models.popgen import generate_linkage_disequilibrium_data

# Generate genotypes with LD
genotypes = generate_linkage_disequilibrium_data(
    n_individuals=1000,
    n_sites=10,
    r_squared_target=0.5
)

# Analyze LD
from metainformant.math.core.utilities import r_squared
# Calculate LD between adjacent sites
ld = r_squared([row[0] for row in genotypes], [row[1] for row in genotypes])
print(f"Observed r²: {ld:.4f}")  # about 0.25: r is about 0.5, so r² is about 0.25
```

## Use Cases

### 1. Testing Population Genetics Methods

Generate synthetic data to test analysis methods:

```python-snippet
from metainformant.simulation.models.popgen import generate_population_sequences
from metainformant.dna.population.analysis import calculate_summary_statistics

# Generate data with known properties
sequences = generate_population_sequences(
    n_sequences=20,
    sequence_length=1000,
    nucleotide_diversity=0.01
)

# Test analysis methods
stats = calculate_summary_statistics(sequences=sequences)
print(f"π: {stats['nucleotide_diversity']:.4f}")
print(f"S: {stats['segregating_sites']}")
print(f"θ_W: {stats['watterson_theta']:.2f}")  # S / a_n for the whole alignment, not per site
```

### 2. Comparing Demographic Models

Compare different demographic scenarios:

```python-snippet
from metainformant.simulation.models.popgen import (
    simulate_bottleneck_population,
    simulate_population_expansion,
)
from metainformant.dna.population.analysis import neutrality_test_suite

# Bottleneck scenario
bottleneck_seqs = simulate_bottleneck_population(
    n_sequences=20,
    sequence_length=1000,
    bottleneck_size=5
)

# Expansion scenario
expansion_seqs = simulate_population_expansion(
    n_sequences=20,
    sequence_length=1000,
    expansion_factor=10.0
)

# Compare neutrality tests
bottleneck_results = neutrality_test_suite(bottleneck_seqs)
expansion_results = neutrality_test_suite(expansion_seqs)

print(f"Bottleneck D: {bottleneck_results['tajima_d']:.4f}")
print(f"Expansion D: {expansion_results['tajima_d']:.4f}")
# These generators perturb founder sequences rather than simulating a
# coalescent, so they do not reproduce textbook neutrality signatures. Use a
# coalescent simulator when a test's null distribution matters.
```

### 3. Generating Training Data

Create datasets for machine learning:

```python-snippet
from metainformant.simulation.models.popgen import generate_genotype_matrix

# Generate large genotype matrix
genotypes = generate_genotype_matrix(
    n_individuals=1000,
    n_sites=10000,
    maf_min=0.05,
    maf_max=0.5,
    hwe=True
)

# Use for GWAS simulation, PCA, etc.
from metainformant.gwas.analysis.structure import compute_pca
pca_result = compute_pca(genotypes, n_components=10)
print(pca_result["status"], len(pca_result["pcs"]))  # success 1000 (one score row per individual)
```

### 4. Simulating Fst Scenarios

Generate populations with different levels of differentiation:

```python-snippet
from metainformant.simulation.models.popgen import generate_two_populations
from metainformant.dna.population import hudson_fst

# Low differentiation
pop1_low, pop2_low = generate_two_populations(
    n_pop1=10, n_pop2=10, sequence_length=1000, fst=0.05, n_sites=500
)
fst_low = hudson_fst(pop1_low, pop2_low)

# High differentiation
pop1_high, pop2_high = generate_two_populations(
    n_pop1=10, n_pop2=10, sequence_length=1000, fst=0.5, n_sites=500
)
fst_high = hudson_fst(pop1_high, pop2_high)

print(f"Low Fst: {fst_low:.4f}")    # about 0.05
print(f"High Fst: {fst_high:.4f}")  # about 0.5
```

## Integration with Analysis

### Full Workflow Example

```python-snippet
from metainformant.simulation.models.popgen import generate_two_populations
from metainformant.dna.population import hudson_fst
from metainformant.dna.population.analysis import calculate_summary_statistics, compare_populations
from metainformant.math.population_genetics.demography import bottleneck_effective_size

# Generate populations
pop1, pop2 = generate_two_populations(
    n_pop1=20,
    n_pop2=20,
    sequence_length=1000,
    fst=0.2,
    n_sites=500
)

# Analyze: summarize each population, then compare the summaries
comparison = compare_populations(
    calculate_summary_statistics(sequences=pop1),
    calculate_summary_statistics(sequences=pop2)
)

print(f"Fst: {hudson_fst(pop1, pop2):.4f}")
print(f"Diversity difference (pop2 - pop1): {comparison['nucleotide_diversity_difference']:.4f}")

# Compare to demographic model
ne_bottleneck = bottleneck_effective_size(
    pre_bottleneck_size=10000,
    bottleneck_size=100,
    post_bottleneck_size=10000
)
print(f"Effective size under bottleneck: {ne_bottleneck:.2f}")
```

## Best Practices

### 1. Use Seeds for Reproducibility

```python
import random
from metainformant.simulation.models.popgen import generate_population_sequences

rng = random.Random(42)

sequences = generate_population_sequences(
    n_sequences=10,
    sequence_length=1000,
    mutation_rate=0.01,
    rng=rng
)
```

### 2. Validate Generated Data

Always check that generated data matches expectations. For `generate_population_sequences` the expectation is the one in the note above: with probability *p* = 0.01 over 1000 positions, π is about 2*p*(1 − *p*) ≈ 0.0198.

```python-snippet
from metainformant.simulation.models.popgen import generate_population_sequences
from metainformant.dna.population import nucleotide_diversity

p = 0.01
sequences = generate_population_sequences(
    n_sequences=40,
    sequence_length=1000,
    nucleotide_diversity=p
)

# Verify
pi_expected = 2 * p * (1 - p)
pi_observed = nucleotide_diversity(sequences)
print(f"Expected π: {pi_expected:.4f}, Observed π: {pi_observed:.4f}")
assert abs(pi_observed - pi_expected) < 0.005  # Allow small variance
```

### 3. Share an Ancestor Between Populations

Each call draws its own random ancestral sequence, so two separate calls give unrelated, unaligned populations. To simulate several populations from one ancestor, pass the same `reference_sequence`:

```python
from metainformant.simulation.models.popgen import generate_population_sequences
from metainformant.dna.population import nucleotide_diversity

# One shared ancestral sequence
ancestor = generate_population_sequences(n_sequences=2, sequence_length=1000)[0]

# Two derived populations with different amounts of change
low = generate_population_sequences(
    n_sequences=20,
    reference_sequence=ancestor,
    nucleotide_diversity=0.005
)
high = generate_population_sequences(
    n_sequences=20,
    reference_sequence=ancestor,
    nucleotide_diversity=0.02
)
print(f"π low: {nucleotide_diversity(low):.4f}, π high: {nucleotide_diversity(high):.4f}")
```

`simulate_bottleneck_population` and `simulate_population_expansion` take no `reference_sequence`; they draw their own founders, so they cannot be started from an existing ancestor.

## Limitations

1. **Simplified models**: These simulations use simplified models. For publication-quality analysis, consider more sophisticated simulators (e.g., msprime, simuPOP).

2. **No coalescent genealogy**: Sequences are drawn independently from one ancestral sequence, and the bottleneck and expansion generators resample founders, so none of them reproduce coalescent haplotype structure. Only `generate_two_populations` (Balding-Nichols allele frequencies) and `generate_site_frequency_spectrum` (the neutral 1/i expectation) follow a population genetic model, and only at the level of allele frequencies.

3. **Linkage**: `generate_linkage_disequilibrium_data` copies each site from the previous one with a fixed probability. It does not model recombination (`recombination_rate` does not change the output) or haplotype blocks. For realistic LD patterns, use specialized tools.

4. **Selection**: The generators above are neutral. `simulate_selection(genotype_matrix, fitness_effects, generations, ...)` in the same module applies fitness-proportional resampling to an existing genotype matrix (a NumPy array, with `fitness_effects` of shape `(n_sites, 3)`), and `simulate_admixture(n_populations, population_sizes, admixture_proportions, generations, ...)` mixes ancestral allele frequencies; neither generates sequences.

## References

- **Hudson, R. R. (2002)**. Generating samples under a Wright-Fisher neutral model of genetic variation. *Bioinformatics*, 18(2), 337-338.

- **Hernandez, R. D. (2008)**. A flexible forward simulator for populations subject to selection and demography. *Bioinformatics*, 24(23), 2786-2787.

- **Kelleher, J., et al. (2016)**. Efficient coalescent simulation and genealogical analysis for large sample sizes. *PLoS Computational Biology*, 12(5), e1004842.

## See Also

- [`dna.population`](../dna/population.md) - Population genetics analysis
- [`dna.population_analysis`](../dna/population_workflows.md) - Analysis workflows
- [`simulation.sequences`](sequences.md) - Basic sequence generation
- [`math.demography`](../math/demography.md) - Demographic models
