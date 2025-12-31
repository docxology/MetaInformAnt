# AI Agents in Mathematical Biology Development

This document outlines AI assistance in developing METAINFORMANT's mathematical and theoretical biology capabilities.

## AI Contributions

### Mathematical Biology Architecture
**Code Assistant Agent** designed:
- Comprehensive mathematical biology framework
- Theoretical model implementations
- Population genetics algorithms
- Statistical analysis utilities

### Analysis Components
**Code Assistant Agent** contributed to:
- Population genetics models
- Coalescent theory implementations
- Epidemiological modeling
- Evolutionary game theory
- Quantitative genetics
- Statistical utility functions (correlation, regression, Fisher's exact test)
- Caching optimizations for expensive computations (Tajima constants, entropy calculations)
- NumPy vectorization for performance improvements

### Recent Functional Improvements (December 2024)

**Code Assistant Agent** implemented key enhancements:

#### Performance Optimizations
- **Vectorized Correlation**: Replaced manual correlation calculation with NumPy's optimized `corrcoef()`
- **Vectorized Regression**: Enhanced linear regression using NumPy's `polyfit()` for better performance
- **Caching Enhancements**: Improved caching for expensive mathematical computations

#### Code Quality Improvements
- **Modern Type Hints**: Updated all type annotations to use `|` union syntax
- **Input Validation**: Added validation using `metainformant.core.validation`
- **Error Handling**: Standardized error messages and exception handling
- **Documentation**: Enhanced docstrings with examples, references, and performance notes

#### Testing Enhancements
- **Edge Case Testing**: Added tests for error conditions and edge cases
- **Validation Testing**: Verified input validation and error handling
- **Performance Testing**: Benchmarking of vectorized operations

### Quality Assurance
**Documentation Agent** assisted with:
- Mathematical biology documentation
- Algorithm explanation and derivation
- Usage examples and mathematical context
- Integration guides for quantitative workflows

## Development Approach

- **Mathematical Rigor**: AI helped implement accurate mathematical models
- **Computational Efficiency**: Optimized algorithms for large datasets
- **Theoretical Accuracy**: Verified implementations against literature
- **Extensibility**: Framework for adding new mathematical methods

## Quality Assurance

- Human oversight ensures mathematical accuracy and biological relevance
- AI assistance accelerates development while maintaining standards
- Comprehensive testing validates mathematical functionality

This mathematical biology infrastructure provides a solid foundation for quantitative biological analysis.

## Complete Function Signatures

### Population Genetics (`popgen.py`)
- `hardy_weinberg_genotype_freqs(allele_a_frequency: float) -> tuple[float, float, float]`
- `selection_update(allele_a_frequency: float, selection_coefficient: float, dominance_coefficient: float = 0.5) -> float`
- `mutation_update(allele_a_frequency: float, mu: float, nu: float) -> float`
- `fixation_probability(allele_a_frequency: float, effective_population_size: float, selection_coefficient: float = 0.0) -> float`
- `heterozygosity_decay(initial_heterozygosity: float, effective_population_size: float, generations: int) -> float`
- `inbreeding_coefficient(effective_population_size: float, generations: int) -> float`
- `equilibrium_heterozygosity_infinite_alleles(effective_population_size: float, mutation_rate: float) -> float`
- `island_model_update(local_frequency: float, migration_rate: float, migrant_pool_frequency: float) -> float`
- `mutation_selection_balance_recessive(mutation_rate: float, selection_coefficient: float) -> float`
- `mutation_selection_balance_dominant(mutation_rate: float, selection_coefficient: float) -> float`
- `effective_population_size_from_heterozygosity(observed_heterozygosity: float, mutation_rate: float) -> float`
- `inbreeding_coefficient_from_fst(fst: float, subpopulations: int = 2) -> float`
- `linkage_disequilibrium_decay_distance(r_squared: float, recombination_rate: float) -> float`
- `coalescent_time_to_mrca(sample_size: int, effective_size: float) -> float`

### Price Equation & Selection (`price.py`)
- `expectation(values: Iterable[float], weights: Iterable[float] | None = None) -> float`
- `covariance(x: Sequence[float], y: Sequence[float]) -> float`
- `variance(values: Sequence[float]) -> float`
- `correlation(x: Sequence[float], y: Sequence[float]) -> float`
- `standard_deviation(values: Sequence[float]) -> float`
- `weighted_variance(values: Sequence[float], weights: Sequence[float]) -> float`
- `weighted_covariance(x: Sequence[float], y: Sequence[float], weights: Sequence[float]) -> float`
- `weighted_correlation(x: Sequence[float], y: Sequence[float], weights: Sequence[float]) -> float`
- `relative_fitness(fitness: Sequence[float]) -> list[float]`
- `selection_differential(phenotype: Sequence[float], fitness: Sequence[float]) -> float`
- `selection_gradient(phenotype: Sequence[float], fitness: Sequence[float]) -> float`
- `selection_intensity(phenotype: Sequence[float], fitness: Sequence[float]) -> float`
- `price_equation(parent_mean: float, offspring_mean: float, parent_variance: float, offspring_variance: float, covariance: float) -> dict[str, float]`
- `delta_mean_trait(phenotype: Sequence[float], fitness: Sequence[float]) -> float`

### Drift-Diffusion Models (`ddm.py`)
- `ddm_decision_time(drift_rate: float, boundary_separation: float, starting_point: float = 0.5, non_decision_time: float = 0.0) -> float`
- `ddm_accuracy(drift_rate: float, boundary_separation: float) -> float`
- `ddm_log_likelihood(rt: float, choice: int, drift_rate: float, boundary_separation: float, starting_point: float = 0.5, non_decision_time: float = 0.0) -> float`
- `fit_ddm_parameters(rt_data: Sequence[float], choice_data: Sequence[int], bounds: dict[str, tuple[float, float]] | None = None) -> dict[str, float]`

### Linkage Disequilibrium (`ld.py`)
- `d_prime(haplotype_freqs: dict[str, float]) -> float`
- `r_squared(haplotype_freqs: dict[str, float]) -> float`
- `lod_score(haplotype_freqs: dict[str, float]) -> float`
- `lewontin_d_prime(haplotype_freqs: dict[str, float]) -> float`
- `hill_bulmer_r_squared(haplotype_freqs: dict[str, float]) -> float`
- `decay_ld_with_distance(initial_ld: float, recombination_rate: float, distance: float) -> float`

### Coalescent Theory (`coalescent.py`)
- `simulate_coalescent(n_samples: int, effective_size: float = 10000, mutation_rate: float = 1e-8) -> dict[str, Any]`
- `coalescent_tree(n_samples: int, effective_size: float) -> dict[str, Any]`
- `coalescent_time_to_mrca(n_samples: int, effective_size: float) -> float`
- `coalescent_branch_lengths(n_samples: int, effective_size: float) -> list[float]`

### Quantitative Genetics (`quantgen.py`)
- `heritability(broader_sense: bool, variance_components: dict[str, float]) -> float`
- `genetic_correlation(trait1: Sequence[float], trait2: Sequence[float], relatedness: Sequence[float]) -> float`
- `response_to_selection(selection_differential: float, heritability: float) -> float`
- `breeder_equation(selection_intensity: float, heritability: float, phenotypic_sd: float) -> float`

### Epidemiological Models (`epidemiology.py`)
- `sir_model(beta: float, gamma: float, s0: float = 0.99, i0: float = 0.01, r0: float = 0.0, days: int = 100) -> dict[str, list[float]]`
- `seir_model(beta: float, gamma: float, sigma: float, s0: float = 0.99, e0: float = 0.0, i0: float = 0.01, r0: float = 0.0, days: int = 100) -> dict[str, list[float]]`
- `basic_reproduction_number(beta: float, gamma: float) -> float`
- `herd_immunity_threshold(r0: float) -> float`

### Statistical Utilities (`__init__.py` and various modules)
- `correlation_coefficient(x: list[float], y: list[float]) -> float`
- `linear_regression(x: list[float], y: list[float]) -> tuple[float, float, float]`
- `fisher_exact_test(a: int, b: int, c: int, d: int) -> tuple[float, float]`
- `shannon_entropy(values: list[float]) -> float`
- `jensen_shannon_divergence(p: list[float], q: list[float]) -> float`

### Evolutionary Game Theory (`egt.py`)
- `replicator_dynamics(fitnesses: Sequence[Sequence[float]], initial_frequencies: Sequence[float], generations: int = 100) -> dict[str, Any]`
- `evolutionarily_stable_strategy(fitness_matrix: Sequence[Sequence[float]], strategy: int) -> bool`
- `prisoners_dilemma_payoff(cooperation_rate: float, temptation: float = 5.0, reward: float = 3.0, punishment: float = 1.0, sucker: float = 0.0) -> dict[str, float]`

### Effective Population Size (`effective_size.py`)
- `census_to_effective_size(census_size: int, variance_in_reproductive_success: float) -> float`
- `temporal_method_heterozygosity(heterozygosities: Sequence[float], generations: int) -> float`
- `linkage_disequilibrium_method(r_squared_values: Sequence[float], recombination_rates: Sequence[float]) -> float`

### Demography (`demography.py`)
- `exponential_growth_model(initial_size: float, growth_rate: float, generations: int) -> list[float]`
- `logistic_growth_model(initial_size: float, carrying_capacity: float, growth_rate: float, generations: int) -> list[float]`
- `age_structure_model(fertility_rates: Sequence[float], survival_rates: Sequence[float], initial_age_structure: Sequence[float], generations: int) -> dict[str, Any]`

### Selection Experiments (`selection_experiments/`)
- `run_selection_experiment(config: dict[str, Any]) -> dict[str, Any]`
- `analyze_selection_response(phenotypes: Sequence[Sequence[float]], generations: int) -> dict[str, Any]`
- `calculate_selection_coefficient(before_freq: float, after_freq: float, generations: int) -> float`
- `plot_selection_trajectory(results: dict[str, Any], output_path: str | Path | None = None) -> matplotlib.figure.Figure`

## Planned Function Signatures - **NOT IMPLEMENTED**

*All math functions are currently planned for future implementation*

### Population Genetics (`popgen.py`) - **NOT IMPLEMENTED**
*Planned: Population genetics models and calculations*
- `hardy_weinberg_genotype_freqs(allele_a_frequency: float) -> tuple[float, float, float]`
- `selection_update(allele_a_frequency: float, selection_coefficient: float, dominance_coefficient: float = 0.5) -> float`
- `mutation_update(allele_a_frequency: float, mu: float, nu: float) -> float`
- `fixation_probability(allele_a_frequency: float, effective_population_size: float, selection_coefficient: float = 0.0) -> float`
- `heterozygosity_decay(initial_heterozygosity: float, effective_population_size: float, generations: int) -> float`
- `inbreeding_coefficient(effective_population_size: float, generations: int) -> float`
- `equilibrium_heterozygosity_infinite_alleles(effective_population_size: float, mutation_rate: float) -> float`
- `island_model_update(local_frequency: float, migration_rate: float, migrant_pool_frequency: float) -> float`
- `mutation_selection_balance_recessive(mutation_rate: float, selection_coefficient: float) -> float`
- `mutation_selection_balance_dominant(mutation_rate: float, selection_coefficient: float) -> float`
- `effective_population_size_from_heterozygosity(observed_heterozygosity: float, mutation_rate: float) -> float`
- `inbreeding_coefficient_from_fst(fst: float, subpopulations: int = 2) -> float`
- `linkage_disequilibrium_decay_distance(r_squared: float, recombination_rate: float) -> float`
- `coalescent_time_to_mrca(sample_size: int, effective_size: float) -> float`

### Price Equation & Selection (`price.py`) - **NOT IMPLEMENTED**
*Planned: Price equation and selection analysis*
- `expectation(values: Iterable[float], weights: Iterable[float] | None = None) -> float`
- `covariance(x: Sequence[float], y: Sequence[float]) -> float`
- `variance(values: Sequence[float]) -> float`
- `correlation(x: Sequence[float], y: Sequence[float]) -> float`
- `standard_deviation(values: Sequence[float]) -> float`
- `weighted_variance(values: Sequence[float], weights: Sequence[float]) -> float`
- `weighted_covariance(x: Sequence[float], y: Sequence[float], weights: Sequence[float]) -> float`
- `weighted_correlation(x: Sequence[float], y: Sequence[float], weights: Sequence[float]) -> float`
- `relative_fitness(fitness: Sequence[float]) -> list[float]`
- `selection_differential(phenotype: Sequence[float], fitness: Sequence[float]) -> float`
- `selection_gradient(phenotype: Sequence[float], fitness: Sequence[float]) -> float`
- `selection_intensity(phenotype: Sequence[float], fitness: Sequence[float]) -> float`
- `price_equation(parent_mean: float, offspring_mean: float, parent_variance: float, offspring_variance: float, covariance: float) -> dict[str, float]`
- `delta_mean_trait(phenotype: Sequence[float], fitness: Sequence[float]) -> float`

### Drift-Diffusion Models (`ddm.py`) - **NOT IMPLEMENTED**
*Planned: Decision-making models*
- `ddm_decision_time(drift_rate: float, boundary_separation: float, starting_point: float = 0.5, non_decision_time: float = 0.0) -> float`
- `ddm_accuracy(drift_rate: float, boundary_separation: float) -> float`
- `ddm_log_likelihood(rt: float, choice: int, drift_rate: float, boundary_separation: float, starting_point: float = 0.5, non_decision_time: float = 0.0) -> float`
- `fit_ddm_parameters(rt_data: Sequence[float], choice_data: Sequence[int], bounds: dict[str, tuple[float, float]] | None = None) -> dict[str, float]`

### Linkage Disequilibrium (`ld.py`) - **NOT IMPLEMENTED**
*Planned: LD analysis and calculations*
- `d_prime(haplotype_freqs: dict[str, float]) -> float`
- `r_squared(haplotype_freqs: dict[str, float]) -> float`
- `lod_score(haplotype_freqs: dict[str, float]) -> float`
- `lewontin_d_prime(haplotype_freqs: dict[str, float]) -> float`
- `hill_bulmer_r_squared(haplotype_freqs: dict[str, float]) -> float`
- `decay_ld_with_distance(initial_ld: float, recombination_rate: float, distance: float) -> float`

### Coalescent Theory (`coalescent.py`) - **NOT IMPLEMENTED**
*Planned: Coalescent simulations and analysis*
- `simulate_coalescent(n_samples: int, effective_size: float = 10000, mutation_rate: float = 1e-8) -> dict[str, Any]`
- `coalescent_tree(n_samples: int, effective_size: float) -> dict[str, Any]`
- `coalescent_time_to_mrca(n_samples: int, effective_size: float) -> float`
- `coalescent_branch_lengths(n_samples: int, effective_size: float) -> list[float]`

### Quantitative Genetics (`quantgen.py`) - **NOT IMPLEMENTED**
*Planned: Quantitative genetics calculations*
- `heritability(broader_sense: bool, variance_components: dict[str, float]) -> float`
- `genetic_correlation(trait1: Sequence[float], trait2: Sequence[float], relatedness: Sequence[float]) -> float`
- `response_to_selection(selection_differential: float, heritability: float) -> float`
- `breeder_equation(selection_intensity: float, heritability: float, phenotypic_sd: float) -> float`

### Epidemiological Models (`epidemiology.py`) - **NOT IMPLEMENTED**
*Planned: Disease modeling and epidemiology*
- `sir_model(beta: float, gamma: float, s0: float = 0.99, i0: float = 0.01, r0: float = 0.0, days: int = 100) -> dict[str, list[float]]`
- `seir_model(beta: float, gamma: float, sigma: float, s0: float = 0.99, e0: float = 0.0, i0: float = 0.01, r0: float = 0.0, days: int = 100) -> dict[str, list[float]]`
- `basic_reproduction_number(beta: float, gamma: float) -> float`
- `herd_immunity_threshold(r0: float) -> float`

### Evolutionary Game Theory (`egt.py`) - **NOT IMPLEMENTED**
*Planned: Game theory models for evolution*
- `replicator_dynamics(fitnesses: Sequence[Sequence[float]], initial_frequencies: Sequence[float], generations: int = 100) -> dict[str, Any]`
- `evolutionarily_stable_strategy(fitness_matrix: Sequence[Sequence[float]], strategy: int) -> bool`
- `prisoners_dilemma_payoff(cooperation_rate: float, temptation: float = 5.0, reward: float = 3.0, punishment: float = 1.0, sucker: float = 0.0) -> dict[str, float]`

### Effective Population Size (`effective_size.py`) - **NOT IMPLEMENTED**
*Planned: Effective population size calculations*
- `census_to_effective_size(census_size: int, variance_in_reproductive_success: float) -> float`
- `temporal_method_heterozygosity(heterozygosities: Sequence[float], generations: int) -> float`
- `linkage_disequilibrium_method(r_squared_values: Sequence[float], recombination_rates: Sequence[float]) -> float`

### Demography (`demography.py`) - **NOT IMPLEMENTED**
*Planned: Population demography models*
- `exponential_growth_model(initial_size: float, growth_rate: float, generations: int) -> list[float]`
- `logistic_growth_model(initial_size: float, carrying_capacity: float, growth_rate: float, generations: int) -> list[float]`
- `age_structure_model(fertility_rates: Sequence[float], survival_rates: Sequence[float], initial_age_structure: Sequence[float], generations: int) -> dict[str, Any]`

### Population Dynamics (`dynamics.py`) - **NOT IMPLEMENTED**
*Planned: Population dynamics modeling*

### F-statistics (`fst.py`) - **NOT IMPLEMENTED**
*Planned: F-statistics and population differentiation*

### Population Genetics Statistics (`popgen_stats.py`) - **NOT IMPLEMENTED**
*Planned: Population genetics statistical methods*

### Selection Models (`selection.py`) - **NOT IMPLEMENTED**
*Planned: Natural selection modeling and analysis*
