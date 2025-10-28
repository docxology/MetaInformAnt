# Mathematical Biology Module

The `math` module provides theoretical and quantitative biology tools, implementing mathematical models from evolutionary theory, population genetics, and decision-making processes.

## Overview

This module contains implementations of key mathematical frameworks in theoretical biology:
- **Population Genetics**: Price equation, kin selection, multilevel selection
- **Decision Theory**: Drift-diffusion models for behavioral analysis
- **Evolutionary Dynamics**: Selection experiments and quantitative trait analysis
- **Statistical Methods**: Epidemiological models and linkage disequilibrium

## Submodules

### Population Genetics (`price.py`, `kin_selection.py`, `multilevel_selection.py`)
Mathematical models of evolutionary processes.

**Key Features:**
- Price equation decomposition of evolutionary change
- Hamilton's rule for kin selection
- Multilevel selection theory
- Heritability and response to selection calculations

**Usage:**
```python
from metainformant.math import price_equation, kin_selection_response

# Price equation analysis
cov, trans, total = price_equation.decompose_change(
    trait_values=[1.0, 1.2, 0.9],
    fitness_values=[0.2, 0.4, 0.1],
    trait_changes=[0.25, 0.35, 0.15]
)

# Kin selection analysis
response = kin_selection_response(
    relatedness=0.5,
    benefit=0.4,
    cost=0.1
)
```

### Decision Theory (`ddm.py`)
Drift-diffusion models for decision-making processes.

**Key Features:**
- Analytic accuracy calculations
- Mean decision time computation
- Parameter estimation and fitting
- Model comparison and validation

**Usage:**
```python
from metainformant.math import ddm_analytic_accuracy, ddm_mean_decision_time

# DDM analysis
accuracy = ddm_analytic_accuracy(drift_rate=0.5, threshold=1.0, noise=0.1)
decision_time = ddm_mean_decision_time(drift_rate=0.5, threshold=1.0, noise=0.1)
```

### Population Dynamics (`dynamics.py`)
Mathematical models of population dynamics and ecological processes.

**Key Features:**
- Logistic growth models
- Lotka-Volterra competition/predation
- Population stability analysis
- Phase space analysis and attractors

**Usage:**
```python
from metainformant.math import dynamics

# Logistic growth simulation
growth_curve = dynamics.logistic_growth(
    initial_population=100,
    carrying_capacity=1000,
    growth_rate=0.1,
    time_steps=100
)

# Lotka-Volterra competition
competition_result = dynamics.lotka_volterra_competition(
    species1_initial=100,
    species2_initial=80,
    alpha=1.2,  # Competition coefficient
    beta=0.8,
    time_steps=200
)
```

### Evolutionary Experiments (`selection_experiments/`)
Natural selection simulation framework with signal processing constraints.

**Key Features:**
- Structural vs. quality trait evolution
- Signal fidelity modeling
- Selection dynamics simulation
- Multi-generation evolutionary experiments

**Usage:**
```python
from metainformant.math.selection_experiments import simulate_generation, simulate_generations

# Single generation simulation
result = simulate_generation(
    s=np.random.normal(0, 1, 10000),
    s_hat=0.8
)

# Multi-generation evolution
results = simulate_generations(
    generations=50,
    n=10000,
    s_hat=0.6
)
```

### Epidemiology (`epidemiology.py`)
Mathematical models of disease dynamics.

**Key Features:**
- SIR (Susceptible-Infected-Recovered) models
- Parameter estimation from data
- Model fitting and validation
- Scenario simulation and forecasting

**Usage:**
```python
from metainformant.math import sir_model

# Run SIR simulation
results = sir_model.simulate(
    population=1000,
    initial_infected=10,
    transmission_rate=0.3,
    recovery_rate=0.1,
    days=100
)
```

### Population Genetics (`popgen.py`, `fst.py`, `ld.py`)
Population genetic statistics and analysis.

**Key Features:**
- Fst calculation and interpretation
- Linkage disequilibrium analysis
- Population structure inference
- Genetic diversity metrics

**Usage:**
```python
from metainformant.math import fst, ld

# Population differentiation
fst_value = fst.calculate_fst(pop1_alleles, pop2_alleles)

# Linkage disequilibrium
ld_coefficient = ld.calculate_ld(allele1_freq, allele2_freq, haplotype_freq)
```

### Coalescent Theory (`coalescent.py`)
Coalescent models for genealogical processes.

**Key Features:**
- Coalescent tree simulation
- Time to most recent common ancestor
- Population size estimation
- Demographic inference

**Usage:**
```python
from metainformant.math import coalescent

# Simulate coalescent tree
tree = coalescent.simulate_tree(n_samples=10, theta=2.0)
t_mrca = coalescent.time_to_mrca(tree)
```

## Theoretical Background

### Price Equation
The Price equation provides a general framework for understanding evolutionary change:
```
Δ̄z = Cov(w, z) + E(wΔz)
```
Where Δ̄z is the change in average trait value, Cov(w, z) represents natural selection, and E(wΔz) represents transmission bias.

### Kin Selection
Hamilton's rule states that a trait will be favored by selection when:
```
rB > C
```
Where r is the genetic relatedness, B is the benefit to the recipient, and C is the cost to the actor.

### Drift-Diffusion Models
DDMs model decision-making as a stochastic process where evidence accumulates over time until a threshold is reached.

## Integration with Other Modules

### With DNA Module
```python
from metainformant.dna import population
from metainformant.math import price_equation

# Theoretical analysis of DNA diversity
pi = population.nucleotide_diversity(sequences)
theoretical_change = price_equation.decompose_change(...)
```

### With Simulation Module
```python
from metainformant.math import ddm_analytic_accuracy
from metainformant.simulation import agents

# Use DDM in agent-based models
accuracy = ddm_analytic_accuracy(drift_rate, threshold, noise)
agent_decision = agents.make_decision(accuracy)
```

## Mathematical Rigor

All implementations are based on peer-reviewed mathematical literature and include:
- Proper handling of edge cases and numerical stability
- Comprehensive parameter validation
- Clear documentation of assumptions and limitations
- References to original theoretical work

## Performance Features

- **Numerical Stability**: Robust algorithms for floating-point computations
- **Vectorization**: NumPy-based implementations for efficiency
- **Caching**: Expensive computations are cached when appropriate
- **Parallel Processing**: CPU-intensive simulations support parallel execution

## Testing

Comprehensive tests ensure mathematical correctness:
- Verification against known analytical solutions
- Numerical stability testing
- Performance benchmarking
- Edge case handling validation

## Dependencies

- **Core**: NumPy for numerical computations
- **Optional**: SciPy for advanced statistical functions
- **Visualization**: Matplotlib for plotting mathematical functions

## Usage Examples

### Price Equation Analysis
```python
from metainformant.math import price_equation

# Analyze evolutionary change
trait_before = [1.0, 1.2, 0.9, 1.1, 0.8]
fitness = [0.2, 0.4, 0.1, 0.3, 0.2]
trait_after = [1.25, 1.35, 0.95, 1.15, 0.85]

selection, transmission, total_change = price_equation.decompose_change(
    trait_before, fitness, trait_after
)

print(f"Selection component: {selection:.4f}")
print(f"Transmission component: {transmission:.4f}")
print(f"Total change: {total_change:.4f}")
```

### Selection Experiments
```python
from metainformant.math.selection_experiments import simulate_generations

# Run evolutionary simulation
results = simulate_generations(
    generations=100,
    n=10000,
    s_hat=0.7,  # Signal fidelity
    fitness_fn=lambda s: s**2  # Quadratic fitness
)

# Analyze results
correlation = results.correlation_s_q()
selection_gradient = results.selection_gradient()
```

This module provides rigorous mathematical tools for theoretical biology, enabling quantitative analysis of evolutionary and behavioral processes.
