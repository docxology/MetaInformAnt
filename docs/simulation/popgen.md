# Simulation: Population Genetics

Synthetic data generation for population genetics analysis.

**Functions**: `generate_population_sequences`, `generate_two_populations`, `generate_genotype_matrix`, `simulate_bottleneck_population`, `simulate_population_expansion`, `generate_site_frequency_spectrum`, `generate_linkage_disequilibrium_data`. The coalescent engine behind the sequence generators: `PopulationEpoch`, `simulate_coalescent_spectra`, `simulate_coalescent_haplotypes`, `haplotypes_to_sequences`.

## Overview

This module provides synthetic data generation for population genetics analysis. You can generate:

- **Neutral coalescent samples**: aligned sequences at a chosen per-site θ, with the neutral site frequency spectrum, π and Tajima's D
- **Bottleneck and expansion samples**: the same coalescent under a demographic history, with the loss of rare variants (positive Tajima's D) after a bottleneck and the excess of rare variants (negative D) after an expansion
- **Two populations with a target Fst** (Balding-Nichols allele frequencies)
- **Genotype matrices** with Hardy-Weinberg equilibrium
- **Site frequency spectra** (the neutral 1/i expectation, plus heuristic expansion and bottleneck shapes)
- **Linkage disequilibrium** data with a chosen adjacent-site r²
- **Coalescent building blocks**: demographic histories, site frequency spectra of many replicates, and haplotype matrices

### Models at a glance

| Function | Model |
|----------|-------|
| `generate_population_sequences` | Neutral Kingman coalescent, constant size, infinite sites, no recombination |
| `simulate_bottleneck_population` (sequences) | The same coalescent with a bottleneck history |
| `simulate_population_expansion` (sequences) | The same coalescent with exponential growth |
| `simulate_bottleneck_population`, `simulate_population_expansion` (statistics dictionaries) | Deterministic size trajectory with heuristic diversity and Poisson mutation counts; not a coalescent |
| `generate_two_populations` | Balding-Nichols allele frequencies |
| `generate_genotype_matrix` | Independent genotypes at per-site allele frequencies |
| `generate_site_frequency_spectrum` | `"constant"`: independent draws from the neutral 1/i expectation; `"expansion"` and `"bottleneck"`: heuristic shapes, not coalescent-derived |
| `generate_linkage_disequilibrium_data` | First-order Markov chain along the sites with exact conditional probabilities |
| `PopulationEpoch`, `simulate_coalescent_spectra`, `simulate_coalescent_haplotypes`, `haplotypes_to_sequences` | The coalescent engine under the sequence generators |

## The Coalescent Model

`generate_population_sequences` and the sequence forms of `simulate_bottleneck_population` and `simulate_population_expansion` draw their samples from the Kingman coalescent with infinite-sites mutation and no recombination:

- **One genealogy per sample.** All sites share a single tree. Every variable site is one mutation on one branch of it and is carried by exactly the sequences below that branch, so any two variable sites pass the four-gamete test (at most three of the four allele combinations occur).
- **Time and scale.** Time runs backwards from the sample in units of `2 * N_ref` generations, where `N_ref` is the reference population size that defines θ = 4 · `N_ref` · μ, the population mutation parameter of the whole sequence. A pair of lineages coalesces at rate 1 / size(*t*), with size(*t*) relative to `N_ref`, and mutations fall on every branch at rate θ / 2 per unit time. Size changes are described by `PopulationEpoch` histories, which rescale the coalescent rates through time exactly.
- **Neutral expectations at constant size.** For *n* sequences, E[S] = θ · a_n segregating sites, E[π] = θ pairwise differences and E[ξ_i] = θ / *i* mutations carried by *i* sequences, where a_n is the sum of 1 / *i* for *i* < *n*. Tajima's D is slightly below 0 on average (about −0.1 for 20 sequences; the standardized statistic is slightly biased).
- **One locus is one draw.** The statistics of a single sample scatter widely around those means. For 20 sequences with a per-site θ of 0.01 over 1000 bp, π has a mean of 0.010 and a standard deviation of about 0.005 across seeds. Average over seeded replicates when an expectation matters.
- **Mutations.** Variable sites are spread uniformly over the A/C/G/T positions of the ancestral sequence, no site mutates twice, and the derived base at a site is drawn uniformly from the three other bases (no transition/transversion bias).

## Functions

### `generate_population_sequences(n_sequences, length=None, *, ...)`

Generate aligned DNA sequences as one neutral coalescent sample: constant population size, infinite-sites mutation and no recombination (see [The Coalescent Model](#the-coalescent-model)).

**Parameters**:
- `n_sequences`: Number of sequences to generate (at least 2)
- `length`: Length of each sequence (alias `sequence_length`); required unless `reference_sequence` is given
- `nucleotide_diversity`, `wattersons_theta`, `theta`, `mutation_rate`: Four ways to set the per-site population mutation parameter θ_site (4 · *N* · μ per site). When several are given, the first in this order decides and the others are ignored:
  1. `nucleotide_diversity`: θ_site equals the per-site π asked for, because E[π] = θ_site
  2. `wattersons_theta`: θ_site equals the per-site Watterson's θ asked for, because E[θ_W] = θ_site
  3. `theta`: used as θ_site (per site, as in `generate_two_populations`)
  4. `mutation_rate`: the per-site mutation rate μ per generation, with the population size taken to be `n_sequences`, so θ_site = 4 · `n_sequences` · `mutation_rate`
  5. none of the four: θ_site = 0.01

  The sample is drawn with θ = θ_site · *L*, for the *L* A/C/G/T positions of the ancestral sequence. The targets are expectations over genealogies: the π and θ_W of a realised sample scatter around them.
- `reference_sequence`: Ancestral sequence; it replaces `length` with its own length. Only its A/C/G/T positions (either case) can mutate, and a derived base keeps the case of the base it replaces. Without one, a random A/C/G/T sequence (about 50% GC) is drawn
- `n_sites`: Accepted and validated (at least 1) but not used. The number of variable sites is an outcome of the genealogy and θ (mean θ · a_n), not an input, and the variable sites are spread uniformly over the sequence
- `rng`: `random.Random` for reproducibility (or `random_seed`, an integer). NumPy's global random state is not touched

**Returns**: List of `n_sequences` strings of equal length. At small θ a sample can be monomorphic (all sequences identical): with θ = θ_site · *L* for the whole sequence, the probability is the product of *j* / (*j* + θ) over *j* = 1 to `n_sequences` − 1.

**Raises**: `ValueError` when `length` is missing, the mutation parameter is not finite, or the infinite-sites model cannot hold the variation; `ValidationError` for a numeric parameter out of range.

**Choosing θ**: The expected number of variable sites is θ_site · a_n · *L*, so a target of *S* variable sites needs θ_site = *S* / (a_n · *L*) (a_n is 3.55 for 20 sequences). Keep θ_site · a_n well below 1: each variable site sits at its own position, so a request that expects more variable sites than the sequence has A/C/G/T positions raises `ValueError`. For example `mutation_rate=0.01` with 10 sequences means θ_site = 0.4, about 113 expected variable sites on 100 positions. The drawn number of variable sites is checked as well, so a request near the limit can fail for some seeds and succeed for others.

**Example**:
```python-snippet
from metainformant.dna.population import nucleotide_diversity, segregating_sites
from metainformant.simulation.models.popgen import generate_population_sequences

# One neutral sample: 10 sequences of 1000 bp with an expected π of 0.01 per site
sequences = generate_population_sequences(
    n_sequences=10,
    sequence_length=1000,
    nucleotide_diversity=0.01,
    random_seed=1,
)

# Analyze the generated data
pi = nucleotide_diversity(sequences)
print(f"Observed π: {pi:.4f}")  # scatters around 0.01: one genealogy, standard deviation about 0.005
print(f"Segregating sites: {segregating_sites(sequences)}")  # mean θ · a_n = 10 · 2.83 = 28
```

Average over seeds to see the expectations:

```python-snippet
import random
from statistics import mean

from metainformant.dna.population import nucleotide_diversity, tajimas_d
from metainformant.simulation.models.popgen import generate_population_sequences

samples = [
    generate_population_sequences(20, 1000, nucleotide_diversity=0.01, rng=random.Random(seed))
    for seed in range(150)
]
print(f"Mean π: {mean(nucleotide_diversity(s) for s in samples):.4f}")  # 0.0100
print(f"Mean D: {mean(tajimas_d(s) for s in samples):.2f}")  # -0.04
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

Simulate a population that went through a bottleneck. Pass `n_sequences` and `sequence_length` as keyword arguments to get sequences; the function's first four positional parameters belong to the statistics form described below.

The sequences are a neutral coalescent sample (see [The Coalescent Model](#the-coalescent-model)) from a population of `initial_size` individuals that shrank to `bottleneck_size` for the last `bottleneck_duration` generations. The sample is drawn at the end of the bottleneck, before any recovery. In units of `2 * initial_size` generations the bottleneck has relative size `bottleneck_size / initial_size` and lasts `bottleneck_duration / (2 * initial_size)`; its strength is set mainly by `bottleneck_duration / (2 * bottleneck_size)`, the coalescent time it spans. Lineages coalesce quickly inside the bottleneck and leave a few long ancestral branches, so diversity falls, rare variants are lost and Tajima's D is positive on average, the textbook signature of a bottleneck observed shortly afterwards. A severe bottleneck can leave the sample monomorphic.

**Parameters**:
- `n_sequences`: Number of sequences to generate (at least 2); passing it selects the sequence form
- `sequence_length`: Length of each sequence (required with `n_sequences`)
- `bottleneck_size`: Population size during the bottleneck, at least 2 and at most `initial_size` (default `max(2, n_sequences // 4)`)
- `bottleneck_duration`: Generations the bottleneck lasts, at least 0 (default 10). A duration of 0 means no bottleneck: the sample is the constant-size sample of `generate_population_sequences` with `theta=pre_bottleneck_diversity`
- `initial_size`: Population size before the bottleneck; it sets the time scale (default `max(1000, 10 * bottleneck_size)`, at least `bottleneck_size`)
- `pre_bottleneck_diversity`: Per-site population mutation parameter of the ancestral population (4 · *N* · μ per site, its expected nucleotide diversity; default 0.01). A value of 0 is honoured and gives a monomorphic sample
- `rng`: `random.Random` for reproducibility

`final_size`, `generations` and `mutation_rate` belong to the statistics form and have no effect on sequences: the sample is drawn before any recovery, `bottleneck_duration` sets the length of the bottleneck, and `pre_bottleneck_diversity` already fixes 4 · *N* · μ.

**Returns**: List of sequences

**Raises**: `ValueError` if `sequence_length` is missing or the infinite-sites model cannot hold the variation (see `generate_population_sequences`); `ValidationError` if a numeric parameter is out of range, including `bottleneck_size` outside 2 to `initial_size`, a negative `bottleneck_duration` and a negative `pre_bottleneck_diversity`.

**Statistics form**: called without `n_sequences`, the function takes `initial_size` (default 1000), `bottleneck_size` (default 100), `final_size` (default `initial_size`), `generations` (default 100) and `mutation_rate` (default 1e-8) instead and returns a dictionary with `population_sizes`, `diversity_trajectory`, `mutation_trajectory`, `bottleneck_start`, `bottleneck_end`, `total_mutations` and `final_diversity`. This is a deterministic size trajectory with heuristic diversity and Poisson mutation counts; it is not a coalescent simulation. The bottleneck occupies the middle third of `generations`.

**Example**:
```python-snippet
import random
from statistics import mean

from metainformant.dna.population import nucleotide_diversity, tajimas_d
from metainformant.simulation.models.popgen import simulate_bottleneck_population

# Samples of 20 sequences drawn at the end of a 10-generation bottleneck to
# 5 individuals in a population of 1000 (the defaults of the sequence form)
samples = [
    simulate_bottleneck_population(
        n_sequences=20,
        sequence_length=1000,
        initial_size=1000,
        bottleneck_size=5,
        bottleneck_duration=10,
        pre_bottleneck_diversity=0.01,
        rng=random.Random(seed),
    )
    for seed in range(150)
]

# One sample is one genealogy, so average over seeds
print(f"Mean π: {mean(nucleotide_diversity(s) for s in samples):.4f}")  # 0.0038, against 0.0100 without the bottleneck
print(f"Mean D: {mean(tajimas_d(s) for s in samples):+.2f}")  # +0.92
```

### `simulate_population_expansion(n_sequences=..., sequence_length=..., ...)`

Simulate a population that underwent recent expansion. Pass `n_sequences` and `sequence_length` as keyword arguments to get sequences; the function's first three positional parameters belong to the statistics form described below.

The sequences are a neutral coalescent sample (see [The Coalescent Model](#the-coalescent-model)) from a population of `initial_size` individuals that grew exponentially by a factor of `expansion_factor` over the last `expansion_time` generations; the sample is drawn when the growth ends. In units of `2 * initial_size` generations the expansion lasts `expansion_time / (2 * initial_size)` and its growth rate is `log(expansion_factor)` divided by that. The large recent population barely coalesces, so external branches dominate the genealogy: rare variants are in excess and Tajima's D is negative on average (about −0.9 for 20 sequences with the defaults), the textbook signature of a recent expansion. The ancestral population has the fixed per-site population mutation parameter 0.01 (4 · *N* · μ per site), so π rises with the expansion (0.013 with the defaults).

**Parameters**:
- `n_sequences`: Number of sequences to generate (at least 2); passing it selects the sequence form
- `sequence_length`: Length of each sequence (required with `n_sequences`)
- `initial_size`: Population size before the expansion; it sets the time scale (default 100, at least 2)
- `expansion_factor`: Factor by which the population grew (default 10.0, at least 1; 1 means no expansion and gives the constant-size sample of `generate_population_sequences` with θ_site = 0.01)
- `final_size`: Another way to set the factor, `final_size / initial_size`; it takes precedence over `expansion_factor` and must be at least `initial_size`
- `expansion_time`: Generations over which the expansion occurs (default 100, at least 1)
- `rng`: `random.Random` for reproducibility

`mutation_rate` belongs to the statistics form and has no effect on sequences, whose ancestral 4 · *N* · μ per site is fixed at 0.01.

**Returns**: List of sequences

**Raises**: `ValueError` if `sequence_length` is missing or the infinite-sites model cannot hold the variation (see `generate_population_sequences`); `ValidationError` if a numeric parameter is out of range, including an `initial_size` below 2, a `final_size` below `initial_size`, an `expansion_factor` below 1 and an `expansion_time` below 1.

**Statistics form**: called without `n_sequences`, the function takes `initial_size` (default 100), `final_size` (default 10000), `expansion_time` (default 100) and `mutation_rate` (default 1e-8) instead and returns a dictionary with `population_sizes`, `diversity_trajectory`, `mutation_trajectory`, `growth_rate`, `expansion_time`, `total_mutations` and `final_diversity`. This is a deterministic size trajectory with heuristic diversity and Poisson mutation counts; it is not a coalescent simulation. The `growth_rate` is the per-generation exponential rate, (ln `final_size` − ln `initial_size`) / `expansion_time`, so it is an output, not an input.

**Example**:
```python-snippet
import random
from statistics import mean

from metainformant.dna.population import nucleotide_diversity, tajimas_d
from metainformant.simulation.models.popgen import simulate_population_expansion

# Samples of 20 sequences after a 10x expansion of a population of 100 over 100 generations
samples = [
    simulate_population_expansion(
        n_sequences=20,
        sequence_length=1000,
        expansion_factor=10.0,
        rng=random.Random(seed),
    )
    for seed in range(150)
]

print(f"Mean π: {mean(nucleotide_diversity(s) for s in samples):.4f}")  # 0.0131, against 0.0100 without the expansion
print(f"Mean D: {mean(tajimas_d(s) for s in samples):+.2f}")  # -0.89
```

### `generate_site_frequency_spectrum(n_samples=None, n_sites=None, *, ...)`

Generate a site frequency spectrum with specified properties.

**Parameters**:
- `n_samples`: Number of sampled sequences (alias `sample_size`; at least 2)
- `n_sites`: Number of polymorphic sites to generate (required)
- `folded`: If True, return the folded SFS (minor allele counts only)
- `demographic_model`: `"constant"` (default; each site's class is an independent draw from the neutral 1/i expectation), `"expansion"` or `"bottleneck"`; anything else raises `ValidationError`
- `parameters`: Optional dictionary of model parameters: `alpha` (default 1.0) for `"expansion"`, a Zipf exponent that only takes effect above 1, and `bottleneck_strength` (default 0.1) for `"bottleneck"`
- `rng`: `random.Random` for reproducibility. It also seeds NumPy's global generator, which the `"expansion"` model's Zipf draw uses

The `"expansion"` and `"bottleneck"` shapes are heuristic: a Zipf draw that favors rare classes, and a mixture of two uniform draws, from the upper half of the classes with probability `bottleneck_strength` and from the lower half otherwise. They are not derived from a coalescent history. For a spectrum under a demographic history, use `simulate_coalescent_spectra` with `epochs`.

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

Generate 0/1 allele data in which neighbouring sites are correlated.

Sites are generated left to right as a first-order Markov chain, one row per individual. Each site keeps its own allele frequency exactly, and each adjacent pair gets the disequilibrium D = r · √(p1 (1 − p1) p2 (1 − p2)) with r = √`r_squared_target`, so the expected r² between adjacent sites equals `r_squared_target`. The probability of allele 1 given the previous site's allele is an exact function of the two frequencies and r. For sites with equal frequencies the correlation between sites *k* apart falls to `r ** k`.

**Parameters**:
- `n_individuals`: Number of individuals (at least 2)
- `n_snps`: Number of sites (at least 2); alias `n_sites`. Taken from the length of `allele_frequencies` when that is given
- `r_squared_target`: Target r² between adjacent sites, within [0, 1] (default 0.25). A sample r² scatters around the target by about 2r(1 − r²)/√`n_individuals`. A pair whose allele frequencies cannot reach the target (a rare allele next to a common one limits r²) is generated at the largest attainable r², and a warning is logged
- `selection_coefficient`: Default 0. When non-zero, its absolute value replaces √`r_squared_target` as the correlation r between adjacent sites (capped at 1); a heuristic coupling, not a model of selection
- `recombination_rate`: Accepted and validated (it must not be negative) but it does not change the output. A first-order chain cannot separate the decay of LD with distance from the LD of adjacent sites, so the model has no recombination step: adjacent LD is set by `r_squared_target` alone and decays geometrically with distance
- `allele_frequencies`: Frequency of allele 1 at each site, each within [0, 1] (default 0.5 at every site); a frequency of 0 or 1 gives a monomorphic site
- `rng`: `random.Random` for reproducibility

**Returns**: List of lists with one row per individual and one column per site, holding 0/1 alleles (not 0/1/2 genotypes)

**Raises**: `ValueError` if neither `n_snps` nor `n_sites` is given; `ValidationError` if a parameter is out of range or `r_squared_target`, an allele frequency or `selection_coefficient` is not a valid number.

**Example**:
```python-snippet
from metainformant.math.core.utilities import r_squared
from metainformant.simulation.models.popgen import generate_linkage_disequilibrium_data

# Generate alleles with LD between adjacent sites
genotypes = generate_linkage_disequilibrium_data(
    n_individuals=1000,
    n_sites=10,
    r_squared_target=0.5
)

# Calculate LD between adjacent sites
ld = r_squared([row[0] for row in genotypes], [row[1] for row in genotypes])
print(f"Observed r²: {ld:.4f}")  # about 0.5, scattering by about 0.02 between runs
```

## Coalescent Engine

The sequence generators are built on four public pieces that can be used directly: a demographic history (`PopulationEpoch`), the site frequency spectra of many replicates (`simulate_coalescent_spectra`), a haplotype matrix (`simulate_coalescent_haplotypes`) and a writer that turns a haplotype matrix into DNA (`haplotypes_to_sequences`). The three functions take an `rng` that accepts a `random.Random`, a `numpy.random.Generator` or an integer seed; `None` draws fresh entropy. None of them touches NumPy's global random state.

### `PopulationEpoch(start, size, growth_rate=0.0)`

One epoch of a piecewise-exponential demographic history, read backwards from the sample. Pass a list of epochs as `epochs=` to the simulators; omit it for a constant-size population. Time `t` is in units of `2 * N_ref` generations backwards from the sample (`t = 0`), and sizes are relative to the reference size `N_ref` that defines θ = 4 · `N_ref` · μ.

**Fields**:
- `start`: Backward time at which the epoch begins; the first epoch starts at 0. An epoch lasts until the next one begins (for ever, in the last)
- `size`: Population size at `start`, relative to `N_ref` (positive)
- `growth_rate`: Forward-in-time exponential growth rate per unit time (default 0). Inside the epoch the relative size at time `t` is `size * exp(-growth_rate * (t - start))`, so 0 keeps the size constant, a positive rate describes a population that was smaller further back (an expansion) and a negative rate one that was larger (a decline)

A history must have finite values, a first epoch starting at 0, strictly increasing starts, positive sizes and a last epoch with `growth_rate >= 0` (otherwise lineages may never coalesce); anything else raises `ValueError` when it is used.

**Example**:
```python-snippet
import math

from metainformant.simulation.models.popgen import PopulationEpoch, simulate_coalescent_spectra

# A bottleneck of relative size 0.05 that began 0.06 and ended 0.01 before the sample
bottleneck = [PopulationEpoch(0.0, 1.0), PopulationEpoch(0.01, 0.05), PopulationEpoch(0.06, 1.0)]

# An instantaneous 100-fold expansion 0.5 before the sample
expansion = [PopulationEpoch(0.0, 100.0), PopulationEpoch(0.5, 1.0)]

# An exponential 10-fold expansion that took 0.1 and ended at the sample
growth = [PopulationEpoch(0.0, 10.0, math.log(10.0) / 0.1), PopulationEpoch(0.1, 1.0)]

spectra = simulate_coalescent_spectra(n_samples=20, n_replicates=1000, theta=10.0, epochs=growth, rng=1)
```

The sign of Tajima's D after a bottleneck depends on when the population is sampled: positive soon after the bottleneck (a sample drawn at its end, as `simulate_bottleneck_population` does) and negative once the population has recovered and grown.

### `simulate_coalescent_spectra(n_samples, n_replicates, *, theta=None, segregating_sites=None, epochs=None, rng=None)`

Simulate the site frequency spectra of many neutral coalescent samples. Each replicate draws a genealogy of `n_samples` sequences (constant size, or the history `epochs`) and places infinite-sites mutations on its branches. Give exactly one of `theta` and `segregating_sites`:

- `theta`: mutations are Poisson with mean θ / 2 per unit of branch length, so under constant size E[S] = θ · a_n and E[ξ_i] = θ / *i*
- `segregating_sites`: every replicate gets exactly that many mutations, each on a branch chosen with probability proportional to its length, on a genealogy drawn without weighting (the approach of `ms -s`; Hudson 2002). This is the distribution of the spectrum given the number of segregating sites when θ is integrated out with the scale-invariant prior 1 / θ, so no mutation rate has to be estimated or assumed. For a fixed θ the genealogy given *S* would be weighted by `T ** S * exp(-theta * T / 2)`, with *T* its total branch length; the difference is small, and p-values built on this null are calibrated from a few tens of mutations upward. This is the null behind `p_value_method="simulation"` in `dna.population.analysis` (see [Neutrality tests with a simulated null](#neutrality-tests-with-a-simulated-null))

**Parameters**:
- `n_samples`: Number of sampled sequences (at least 2)
- `n_replicates`: Number of independent samples to simulate (at least 1). Large requests are simulated in vectorised blocks, so memory stays bounded
- `theta`: Population mutation parameter of the whole sequence, 4 · `N_ref` · μ
- `segregating_sites`: Number of mutations in every replicate (a non-negative integer)
- `epochs`: Demographic history (see `PopulationEpoch`); constant size when omitted
- `rng`: `random.Random`, `numpy.random.Generator` or integer seed

**Returns**: Integer array of shape `(n_replicates, n_samples - 1)` whose entry `[r, i - 1]` counts the mutations of replicate `r` carried by exactly `i` of the sequences (ancestral state known)

**Raises**: `ValueError` if neither or both of `theta` and `segregating_sites` are given, or the epochs or the number of segregating sites are invalid; `ValidationError` for `n_samples` or `n_replicates` out of range.

**Example**:
```python-snippet
from metainformant.simulation.models.popgen import simulate_coalescent_spectra

# Poisson mutation load: E[xi_i] = theta / i
spectra = simulate_coalescent_spectra(n_samples=20, n_replicates=2000, theta=10.0, rng=1)
print(spectra.shape)  # (2000, 19)
print(spectra.mean(axis=0)[:3])  # about [10.0, 5.0, 3.3]

# Conditioned on the number of mutations: every replicate has exactly 30
conditioned = simulate_coalescent_spectra(20, 5, segregating_sites=30, rng=1)
print(conditioned.sum(axis=1))  # [30 30 30 30 30]
```

### `simulate_coalescent_haplotypes(n_samples, theta, *, epochs=None, rng=None)`

Simulate one neutral coalescent sample with infinite-sites mutation as a haplotype matrix. A genealogy of `n_samples` sequences is drawn (constant size, or the history `epochs`) and every branch of length *l* receives a Poisson number of mutations with mean θ · *l* / 2. A mutation is carried by exactly the sequences below its branch, and no site mutates twice. There is no recombination, so any two columns pass the four-gamete test.

**Parameters**:
- `n_samples`: Number of sampled sequences (at least 2)
- `theta`: Population mutation parameter of the whole sequence, 4 · `N_ref` · μ, with `N_ref` the reference size that sets the time unit of `epochs`; it must be finite and non-negative
- `epochs`: Demographic history (see `PopulationEpoch`); constant size when omitted
- `rng`: `random.Random`, `numpy.random.Generator` or integer seed

**Returns**: `uint8` array of shape `(n_samples, S)`: entry `[i, j]` is 1 when sequence `i` carries the derived allele of segregating site `j` and 0 for the ancestral allele. Columns are in random order and none is invariant; *S* is random (mean θ · a_n under constant size).

**Raises**: `ValueError` if `theta` is not finite or the epochs are invalid; `ValidationError` if `n_samples` is below 2 or `theta` is negative.

### `haplotypes_to_sequences(haplotypes, length=None, *, reference_sequence=None, rng=None)`

Write an infinite-sites haplotype matrix as aligned DNA sequences. Every column becomes one segregating site at its own position, drawn uniformly without replacement from the A/C/G/T positions of the ancestral sequence; sequences with a 0 there keep the ancestral base and sequences with a 1 carry a derived base drawn uniformly from the three other bases. Positions and bases are seeded through `rng`.

**Parameters**:
- `haplotypes`: Matrix of shape `(n_sequences, n_sites)` of 0 (ancestral) and 1 (derived) alleles, as returned by `simulate_coalescent_haplotypes`
- `length`: Sequence length; required without `reference_sequence`
- `reference_sequence`: Ancestral sequence, which sets the length. Positions holding anything but A/C/G/T (either case) never mutate, and a derived base keeps the case of the base it replaces. Without one, `length` bases are drawn uniformly from A/C/G/T
- `rng`: `random.Random`, `numpy.random.Generator` or integer seed

**Returns**: One string per row of `haplotypes`

**Raises**: `ValueError` if the matrix is malformed, neither `length` nor `reference_sequence` is given (or they disagree), or there are more segregating sites than A/C/G/T positions.

**Example**: a coalescent sample under a chosen demographic history, written on a chosen ancestor:
```python-snippet
import math
import random

from metainformant.simulation.models.popgen import (
    PopulationEpoch,
    haplotypes_to_sequences,
    simulate_coalescent_haplotypes,
)
from metainformant.simulation.models.sequences import generate_random_dna

ancestor = generate_random_dna(500, rng=random.Random(1))
growth = [PopulationEpoch(0.0, 10.0, math.log(10.0) / 0.1), PopulationEpoch(0.1, 1.0)]

# theta is for the whole sequence: 0.01 per site over 500 sites
haplotypes = simulate_coalescent_haplotypes(n_samples=20, theta=5.0, epochs=growth, rng=2)
sequences = haplotypes_to_sequences(haplotypes, reference_sequence=ancestor, rng=2)
print(haplotypes.shape)  # (20, S): one column per segregating site
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

Compare different demographic scenarios. One sample is one genealogy, so average over seeds:

```python-snippet
import random
from statistics import mean

from metainformant.dna.population import tajimas_d
from metainformant.simulation.models.popgen import (
    simulate_bottleneck_population,
    simulate_population_expansion,
)

seeds = range(150)

# Bottleneck scenario: sampled at the end of the bottleneck
bottleneck_d = [
    tajimas_d(simulate_bottleneck_population(n_sequences=20, sequence_length=1000, rng=random.Random(seed)))
    for seed in seeds
]

# Expansion scenario: sampled when the growth ends
expansion_d = [
    tajimas_d(simulate_population_expansion(n_sequences=20, sequence_length=1000, rng=random.Random(seed)))
    for seed in seeds
]

print(f"Bottleneck D: {mean(bottleneck_d):+.2f}")  # +0.92
print(f"Expansion D: {mean(expansion_d):+.2f}")  # -0.89
# The neutral constant-size sample of generate_population_sequences averages about -0.04
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

### Neutrality tests with a simulated null

The neutrality tests in `dna.population.analysis` (`calculate_tajima_d`, `calculate_fu_li_d`, `calculate_fu_li_f`, `calculate_fay_wu_h` and `neutrality_test_suite`) take the keyword-only options `p_value_method` (`"normal"`, the default, or `"simulation"`), `n_simulations` (default 2000) and `random_state`. With `"simulation"` the p-value is a Monte Carlo p-value against `simulate_coalescent_spectra(..., segregating_sites=S)`, where *S* is the number of mutations in the sample. The null is a constant-size neutral population, so a small p-value on a bottleneck or expansion sample reflects its demography. See [Neutrality p-values](../dna/population.md#neutrality-p-values-analysis) for the exact definition.

```python-snippet
import random

from metainformant.dna.population.analysis import calculate_tajima_d, neutrality_test_suite
from metainformant.simulation.models.popgen import simulate_population_expansion

sequences = simulate_population_expansion(
    n_sequences=40, sequence_length=2000, expansion_factor=100.0, rng=random.Random(1)
)

# Normal approximation (the default)
d, p_normal = calculate_tajima_d(sequences)

# Monte Carlo p-value against the neutral constant-size coalescent; the seed makes it reproducible
d, p_simulated = calculate_tajima_d(sequences, p_value_method="simulation", n_simulations=2000, random_state=1)

suite = neutrality_test_suite(sequences, p_value_method="simulation", random_state=1)
print(suite["p_value_method"], suite["n_simulations"])  # simulation 2000
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
    mutation_rate=1e-4,  # θ per site = 4 * 10 * 1e-4 = 0.004
    rng=rng
)
```

The same seed gives the same sequences. The sequence generators draw only from `rng`, so seeding does not depend on, or change, NumPy's global random state.

### 2. Validate Generated Data

Always check that generated data matches expectations. One sample is one genealogy, so compare the mean over seeded replicates with the target, within its standard error:

```python-snippet
import random
from statistics import mean, stdev

from metainformant.dna.population import nucleotide_diversity
from metainformant.simulation.models.popgen import generate_population_sequences

target = 0.01
values = [
    nucleotide_diversity(
        generate_population_sequences(
            n_sequences=40, sequence_length=1000, nucleotide_diversity=target, rng=random.Random(seed)
        )
    )
    for seed in range(100)
]

# Verify
pi_observed = mean(values)
standard_error = stdev(values) / len(values) ** 0.5
print(f"Expected π: {target:.4f}, Observed π: {pi_observed:.4f} ± {standard_error:.4f}")
assert abs(pi_observed - target) < 4 * standard_error
```

### 3. Share an Ancestor Between Populations

Each call draws its own random ancestral sequence, so two separate calls give unrelated, unaligned populations. To simulate several populations from one ancestral sequence, pass the same `reference_sequence`:

```python
from metainformant.simulation.models.popgen import generate_population_sequences
from metainformant.simulation.models.sequences import generate_random_dna
from metainformant.dna.population import nucleotide_diversity

# One shared ancestral sequence
ancestor = generate_random_dna(1000)

# Two populations with different amounts of variation
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
print(f"π low: {nucleotide_diversity(low):.4f}, π high: {nucleotide_diversity(high):.4f}")  # each is one draw around 0.005 and 0.02
```

Each call draws its own genealogy, so the populations share the ancestral bases but no population history. Use `generate_two_populations` for a target Fst.

`simulate_bottleneck_population` and `simulate_population_expansion` take no `reference_sequence`; they draw their own ancestral sequence. To start a demographic history from an existing ancestor, write `simulate_coalescent_haplotypes(..., epochs=...)` onto it with `haplotypes_to_sequences(..., reference_sequence=ancestor)`.

## Limitations

1. **Simplified models**: These simulations use compact models for tests, benchmarks and teaching. For publication-quality analysis, consider more sophisticated simulators (e.g., msprime, simuPOP).

2. **One genealogy, no recombination**: The coalescent samples carry a single genealogy for the whole sequence (infinite sites, no recombination, one panmictic population). Variable sites are not independent loci: windows along a sequence share the tree, so they are not replicates, and any two variable sites pass the four-gamete test. Demographic histories are piecewise constant or exponential changes of one population size, with no migration or structure. Draw one sample per seed and average over seeds when a distribution matters.

3. **Linkage**: `generate_linkage_disequilibrium_data` is a first-order Markov chain along the sites. Adjacent-site r² is set by `r_squared_target`, and for equal allele frequencies the correlation decays geometrically with distance. It does not model recombination (`recombination_rate` does not change the output) or haplotype blocks. For realistic LD patterns, use specialized tools.

4. **Selection**: The generators above are neutral. `simulate_selection(genotype_matrix, fitness_effects, generations, ...)` in the same module applies fitness-proportional resampling to an existing genotype matrix (a NumPy array, with `fitness_effects` of shape `(n_sites, 3)`), and `simulate_admixture(n_populations, population_sizes, admixture_proportions, generations, ...)` mixes ancestral allele frequencies; neither generates sequences.

5. **Heuristic generators**: The `"expansion"` and `"bottleneck"` models of `generate_site_frequency_spectrum` and the statistics-dictionary forms of `simulate_bottleneck_population` and `simulate_population_expansion` are heuristic, not coalescent simulations. `generate_two_populations` and `generate_genotype_matrix` model allele frequencies, not genealogies.

## References

- **Kingman, J. F. C. (1982)**. The coalescent. *Stochastic Processes and their Applications*, 13(3), 235-248.

- **Hudson, R. R. (1990)**. Gene genealogies and the coalescent process. *Oxford Surveys in Evolutionary Biology*, 7, 1-44.

- **Griffiths, R. C., & Tavaré, S. (1994)**. Sampling theory for neutral alleles in a varying environment. *Philosophical Transactions of the Royal Society B*, 344(1310), 403-410.

- **Watterson, G. A. (1975)**. On the number of segregating sites in genetical models without recombination. *Theoretical Population Biology*, 7(2), 256-276.

- **Fu, Y.-X. (1995)**. Statistical properties of segregating sites. *Theoretical Population Biology*, 48(2), 172-197.

- **Hudson, R. R. (2002)**. Generating samples under a Wright-Fisher neutral model of genetic variation. *Bioinformatics*, 18(2), 337-338.

- **Hernandez, R. D. (2008)**. A flexible forward simulator for populations subject to selection and demography. *Bioinformatics*, 24(23), 2786-2787.

- **Kelleher, J., et al. (2016)**. Efficient coalescent simulation and genealogical analysis for large sample sizes. *PLoS Computational Biology*, 12(5), e1004842.

## See Also

- [`dna.population`](../dna/population.md) - Population genetics analysis
- [`dna.population_analysis`](../dna/population_workflows.md) - Analysis workflows
- [`simulation.sequences`](sequences.md) - Basic sequence generation
- [`math.coalescent`](../math/coalescent.md) - Analytical coalescent expectations
- [`math.demography`](../math/demography.md) - Demographic models
