# Specification: population_genetics

## 🎯 Scope
Population genetics submodule.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 9 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `coalescent.py`
- `core.py`
- `demography.py`
- `effective_size.py`
- `fst.py`
- `ld.py`
- `selection.py`
- `statistics.py`

### Behavior Contracts

#### `fst.py`
- `fst_from_allele_freqs(pop1_freqs, pop2_freqs=None)` — moment F_ST pooled over
  loci: `sum_j var_p_j / sum_j (var_p_j + mean_i p_ij (1 - p_ij))`. For two
  populations this is algebraically identical to `(Ht - Hs) / Ht`, and the
  single-list `[p1, p2]` call equals the two-list call for one locus.
- `fst_from_allele_freq_matrix(pop_freqs)` — same estimator for k >= 2
  populations; `var_p` uses the unbiased `k/(k-1)` sample variance for `k > 2`
  and the population variance for `k == 2`.
- `weirs_fst(population_counts)` — Weir & Cockerham (1984) single-locus
  variance components (`a`, `b`, `c`) from per-population haplotype (gene-copy)
  counts; `c = 0` because haplotype counts carry no diploid heterozygosity, and
  `theta = sum_u a_u / sum_u (a_u + b_u + c_u)` clamped to [0, 1]. The legacy
  `(counts, labels)` CRC32-jitter signature was removed.

#### `ld.py`
- `ld_coefficients(pA, pa, pB, pb, pAB)` — frequency mode returns signed `(D, D')`.
- `ld_coefficients(genotypes)` — each row is one PHASED haplotype `[g1, g2]`;
  returns `{"D", "D_prime", "r_squared"}`.
- `ld_coefficients(genotypes, phased=False)` — each row is an unphased diploid
  `[[a1, a2], [b1, b2]]`; haplotype frequencies are recovered by a deterministic
  EM routine before computing D, D' and r^2.
- `ld_decay_r2(distances, r_squared_values, max_distance=None)` — bounded least-squares
  fit of `r^2 = a * exp(-k * d)` (`a >= 0`, `k >= 0`) over every pair with equal weight;
  nothing is binned and the rate comes back per unit of distance. Returns
  `{"decay_rate": k, "half_decay_distance": ln 2 / k, "initial_r2": a}`. Flat or rising
  data report no decay (`decay_rate` 0.0, `half_decay_distance` `inf`, `initial_r2` the
  mean r^2); fewer than three distinct distances after `max_distance` return only
  `decay_rate` and `half_decay_distance`. `k` is a descriptive constant, not an estimate
  of `4 * Ne * c`. `ValueError` for unequal lengths, non-finite values or a negative
  distance; `RuntimeError` if the optimiser does not converge. A scalar first argument
  is the closed form `ld_decay_r2(r2_start, recombination_rate=c, generations=t)` =
  `r2_start * (1 - c) ** (2 * t)`.

#### `effective_size.py`
- `effective_size_sex_ratio(male_size, female_size, sex_ratio=0.5)` — always
  `4 * Nm * Nf / (Nm + Nf)`. `sex_ratio` must lie strictly between 0 and 1 but does not
  change the result, because the breeding sex ratio is already given by `Nm` and `Nf`.

#### `statistics.py`
- `effective_size_from_family_size_variance(family_sizes, Vk=None)` — with `Vk`,
  `family_sizes` is the census size N (an int or a float; a bool is rejected) and the
  result is `(4N - 2) / (Vk + 2)`. Without `Vk`, a list, tuple or array of offspring
  counts gives `k / (V - 1)` with `V` the population variance (`inf` when `V <= 1`, 0.0 for
  an empty list). Wrong types raise `TypeError`; `N <= 0`, a non-finite value or a
  negative `Vk` raise `ValueError`.
- `calculate_confidence_intervals(values, confidence_level=0.95, method="bootstrap")` —
  the bootstrap draws from a private `numpy.random.RandomState(42)`: the interval for a
  list is always the same and numpy's global random state is never reseeded or advanced.
  A list and the equivalent 1-D array give the same interval. A dictionary maps names to
  observation lists and returns one interval per name; a single point value raises
  `ValueError` (no sampling information).
- `compare_statistics(stats1, stats2, test_type="ttest")` — runs exactly the requested
  test from `COMPARISON_TESTS` (`"ttest"` Welch, `"mannwhitney"` two-sided, `"kruskal"`);
  unknown types and untestable input raise `ValueError` rather than returning p = 1.0.
  Dictionaries of summary statistics support only `"ttest"`.

#### `coalescent.py`
- `simulate_coalescent(n_samples, effective_size, mutation_rate, *, random_state=None)` —
  `random_state` is an int, a `numpy.random.Generator` or a `numpy.random.RandomState`.
  A seed or generator reproduces the genealogy whatever numpy's global state is, and a
  seeded run leaves that state unchanged. `None` draws from numpy's global generator.
