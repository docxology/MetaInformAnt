# Specification: models

## 🎯 Scope
Simulation models subpackage.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 5 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `agents.py`
- `popgen.py`
- `rna.py`
- `sequences.py`

### Function contracts (`popgen.py`, behavior pinned by tests)

Pinned by `tests/simulation/test_simulation_popgen.py` and
`tests/simulation/test_simulation_coalescent.py`; the linkage-disequilibrium generator by
`tests/simulation/test_simulation_r3_fixes.py`.

- **Model.** The sequence generators sample the neutral Kingman coalescent with
  infinite-sites mutation and no recombination: one genealogy per call, every variable
  site one mutation on one branch. Time runs backwards from the sample in units of
  `2 * N_ref` generations, a pair of lineages coalesces at rate `1 / size(t)` and
  mutations fall on every branch at rate `theta / 2`, with `theta = 4 * N_ref * mu` for
  the whole sequence. Under constant size `E[S] = theta * a_n`, `E[pi] = theta`,
  `E[xi_i] = theta / i` and `P(S = 0)` is the product of `j / (j + theta)` over `j < n`.
- `PopulationEpoch(start, size, growth_rate=0.0)`: frozen dataclass; the relative size at
  backward time `t` inside an epoch is `size * exp(-growth_rate * (t - start))`, and the
  coalescent rates are rescaled through time exactly (integrated hazard of `1 / size`). A
  history needs finite values, a first epoch at 0, strictly increasing starts, positive
  sizes and a last epoch with `growth_rate >= 0`; otherwise `ValueError`.
- `simulate_coalescent_spectra(n_samples, n_replicates, *, theta=None,
  segregating_sites=None, epochs=None, rng=None)`: integer array of shape
  `(n_replicates, n_samples - 1)`, entry `[r, i - 1]` the mutations of replicate `r`
  carried by `i` sequences. Exactly one of `theta` (Poisson load) and `segregating_sites`
  (exactly that many mutations per replicate, each on an unweighted genealogy in
  proportion to branch length: the `ms -s` approach) is required, else `ValueError`.
  `rng` is a `random.Random`, a `numpy.random.Generator`, an integer or `None`; numpy's
  global state is never used. Large requests run in blocks without changing the model.
- `simulate_coalescent_haplotypes(n_samples, theta, *, epochs=None, rng=None)`: `uint8`
  matrix `(n_samples, S)` of derived-allele indicators; any two columns pass the
  four-gamete test; it samples the same model as the spectrum simulator.
- `haplotypes_to_sequences(haplotypes, length=None, *, reference_sequence=None, rng=None)`:
  one segregating site per column at distinct positions drawn without replacement from
  the A/C/G/T positions of the ancestral sequence; the derived base is uniform among the
  other three and keeps the case of the base it replaces; other positions never mutate.
  `ValueError` for a malformed matrix, a missing or disagreeing `length` and
  `reference_sequence`, or more sites than A/C/G/T positions.
- `generate_population_sequences`: one neutral coalescent sample. The per-site theta is
  the first given of `nucleotide_diversity`, `wattersons_theta`, `theta` and
  `mutation_rate` (`4 * n_sequences * mutation_rate`), default 0.01, applied to the A/C/G/T
  sites of the ancestral sequence. A sample can be monomorphic at small theta. `n_sites`
  is validated (at least 1) and unused. A constant-size request that expects more
  variable sites than A/C/G/T positions, or a draw that needs more, raises `ValueError`.
  Numpy's global state is untouched.
- `simulate_bottleneck_population` with `n_sequences`: a sample drawn at the end of the
  bottleneck, epochs `[PopulationEpoch(0, bottleneck_size / initial_size),
  PopulationEpoch(bottleneck_duration / (2 * initial_size), 1)]`, per-site theta
  `pre_bottleneck_diversity` (default 0.01; 0 is honoured). `bottleneck_duration=0` gives
  the constant-size sample. `ValidationError` unless `2 <= bottleneck_size <=
  initial_size` and `bottleneck_duration >= 0`. `final_size`, `generations` and
  `mutation_rate` do not affect sequences. Mean Tajima's D is positive (about +0.9 for 20
  sequences of 1000 bp).
- `simulate_population_expansion` with `n_sequences`: a sample drawn when exponential
  growth ends, epochs `[PopulationEpoch(0, f, log(f) / T), PopulationEpoch(T, 1)]` with
  `T = expansion_time / (2 * initial_size)` and `f = expansion_factor` (`final_size /
  initial_size` when `final_size` is given, which wins), ancestral per-site theta fixed at
  0.01. `expansion_factor=1.0` gives the constant-size sample; `ValidationError` for an
  `expansion_factor` below 1, a `final_size` below `initial_size` or an `expansion_time`
  below 1. `mutation_rate` does not affect sequences. Mean Tajima's D is negative (about
  -0.9).
- The statistics forms of both functions (no `n_sequences`) are deterministic size
  trajectories with heuristic diversity and Poisson mutation counts, not coalescent
  simulations; they seed numpy's global generator from `rng`.
- `generate_linkage_disequilibrium_data`: a first-order Markov chain along the sites with
  exact conditional probabilities. Marginal allele frequencies are exact and adjacent-site
  r-squared equals `r_squared_target` (default 0.25), capped at the attainable value with
  a logged warning when the frequencies cannot reach it. A non-zero
  `selection_coefficient` replaces the correlation `r` by its absolute value.
  `recombination_rate` is validated and does not change the output.
- `generate_site_frequency_spectrum`: `"constant"` draws each site's class from the
  neutral `1 / i` expectation; `"expansion"` and `"bottleneck"` are heuristic shapes, not
  coalescent-derived. `generate_two_populations` (Balding-Nichols) is a model of allele
  frequencies.
