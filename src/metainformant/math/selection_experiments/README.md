# Natural Selection Experiments

This package contains a comprehensive simulation of natural selection dynamics with signal processing constraints, implementing the mathematical framework from quantitative genetics and evolutionary theory.

## Overview

This example is a computational implementation, modifying the mathematical and computational framework presented in:

### Source Paper: Qualia & Natural Selection
**Ryan Williams** (2025)  
*[Qualia & Natural Selection: Formal Constraints on the Evolution of Consciousness](https://arxiv.org/abs/2505.05480)*  
*[DOI: 10.48550/arXiv.2505.05480](https://doi.org/10.48550/arXiv.2505.05480)*

> This paper explores foundational questions about the relationship of qualia to natural selection. The primary result is a derivation of specific formal conditions under which structural systems subject to natural selection can convey consistent effects in an associated qualitative domain, placing theoretical and empirical constraints on theories of consciousness. In order to achieve this result, information-theoretic measures are developed to quantify the mutual determinability between structure and quality, quantifying fidelity between the two domains. The fidelities represented by that space are then incorporated into the Price Equation to yield key bounds on the transmission of selective effects between domains. Finally, transmission of higher-order structures between domains is explored. Placement within a broader philosophical context can be found in the companion paper Structure & Quality.

The experiments explore how natural selection operates on two types of traits:

- **Structural traits (s)**: The underlying physical characteristics of organisms
- **Quality traits (q)**: Observable signals that correlate with structural traits but may be noisy

The key insight is that selection acts on quality traits, but evolution occurs in structural traits. The relationship between these traits is mediated by signal fidelity, represented by the parameter \(\hat{s}\).

## Mathematical Framework

### Core Functions

#### Signal Mapping: \(\phi(s)\)

The relationship between structural and quality traits is defined by:

\[\phi(s) = 2s^2 + 4 + \epsilon\]

Where \(\epsilon \sim \mathcal{N}(0, \sigma^2)\) and \(\sigma^2\) is determined by the signal fidelity \(\hat{s}\):

\[\sigma^2 = \frac{\hat{s}^{-2} - 1}{\text{Var}(\bar{\phi}(s))}\]

#### Selection Dynamics

Natural selection operates through fitness-weighted reproduction:

\[w_i = f(s_i) + \eta_i\]

Where \(f(s_i)\) is the fitness function and \(\eta_i\) is environmental noise.

#### Price Equation Components

The evolutionary change is decomposed using the Price equation:

\[\Delta \bar{z} = \text{Cov}(w, z) + E(w \Delta z)\]

For our traits:

- \(\Delta \bar{s} = \text{Cov}(w, s) + E(w \Delta s)\)
- \(\Delta \bar{q} = \text{Cov}(w, q) + E(w \Delta q)\)

## Experimental Scenarios

### 1. Trait-Environment Correlation

**File**: `outputs/plot-s-vs-q.png`

- Examines the mapping between structural traits and quality signals
- Shows how signal fidelity affects the correlation \(\rho_{s,q}\)
- High \(\hat{s}\) = strong correlation, low noise in signal
- Low \(\hat{s}\) = weak correlation, high noise in signal

### 2. Selection on Correlated Traits

**File**: `outputs/plot-sq-vs-w.png`

- Demonstrates how selection acts on quality traits but affects structural traits
- Shows different selection gradients for \(s\) vs \(q\)
- Illustrates the concept of "selection on signals"

### 3. Rebound Selection

**File**: `outputs/plot-ns-rebound.png`

- Selection with strong directional bias
- Shows how structural traits can "rebound" when selection is strong
- Parameters: \(\hat{s}=0.6\), directional mutation

### 4. Inverse Selection

**File**: `outputs/plot-ns-inverse.png`

- Selection with opposing forces on signals and structures
- Demonstrates potential for evolutionary conflict
- Parameters: \(\hat{s}=0.6\), opposing selection gradients

### 5. Neutral Selection

**File**: `outputs/plot-ns.png`

- Weak selection scenario
- Shows drift-dominated dynamics
- Parameters: \(\hat{s}=0.1\), weak selection

### 6. Quality Signal Limit (QSL)

**File**: `outputs/plot-ns-qsl.png`

- Extreme scenario with very low signal fidelity
- Demonstrates breakdown of signal-based selection
- Parameters: \(\hat{s}=0.006\), very weak signal

## Key Parameters

| Parameter | Description | Range |
|-----------|-------------|-------|
| \(\hat{s}\) | Signal fidelity | 0-1 |
| \(n\) | Population size | 1,000-1,000,000 |
| \(generations\) | Number of generations | 20-100 |
| \(\Delta s\) | Structural mutation | Various distributions |
| \(w\) | Fitness function | Various forms |

## Output Metrics

Each experiment reports:

- **\(\rho_{s,q}\)**: Correlation between structural and quality traits
- **\(\rho_{w,s}\)**: Selection gradient on structural traits
- **\(\rho_{w,q}\)**: Selection gradient on quality traits
- **QSC**: Quality Selection Coefficient
- **CV**: Coefficient of variation for traits

## Usage

### Command Line

```bash
# Run all experiments with default settings (generates individual PNG plots)
python -m metainformant math selection replay

# Run experiments with custom output directory
python -m metainformant math selection replay --dest /path/to/output

# Create composite graphical abstract (all experiments combined)
python -m metainformant math selection abstract

# Create composite abstract with custom output directory
python -m metainformant math selection abstract --dest /path/to/output
```

### Python API

```python
from metainformant.math.selection_experiments import (
    simulate_generation, simulate_generations,
    display_s_vs_q, display_generation
)

# Single generation analysis
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

### Signal Processing in Evolution

The model implements the concept that biological signals (like mate choice cues, warning coloration, etc.) must balance:

1. **Fidelity**: How well the signal represents the underlying trait
2. **Cost**: Energetic cost of producing/maintaining signals
3. **Selection**: Differential survival based on signal perception

### Evolutionary Game Theory

The framework connects to evolutionary game theory through:

- Signal costs as strategy payoffs
- Receiver responses as evolutionary stable strategies
- Coevolution of signalers and receivers

### Quantitative Genetics

Links to quantitative genetics via:

- Heritability of signals vs. structures
- Selection response equations
- Genetic constraints on signal evolution

## References and Further Reading

1. **Primary Source**: Williams, R. (2025). Qualia & Natural Selection: Formal Constraints on the Evolution of Consciousness. arXiv preprint arXiv:2505.05480. [https://arxiv.org/abs/2505.05480](https://arxiv.org/abs/2505.05480)
2. **Price Equation**: Price, G. R. (1970). Selection and covariance. Nature.
3. **Signal Evolution**: Zahavi, A. (1975). Mate selection-a selection for a handicap. Journal of Theoretical Biology.
4. **Honest Signaling**: Grafen, A. (1990). Biological signals as handicaps. Journal of Theoretical Biology.
5. **Quantitative Genetics**: Falconer, D. S., & Mackay, T. F. C. (1996). Introduction to Quantitative Genetics.

## Output Files

### Individual Experiment Plots
All plots are saved as high-resolution PNG files in the `outputs/` subfolder:

- **`plot-s-vs-q.png`**: Structural trait vs quality signal relationship
- **`plot-sq-vs-w.png`**: Selection on correlated traits
- **`plot-ns-rebound.png`**: Rebound selection dynamics over generations
- **`plot-ns-inverse.png`**: Inverse selection dynamics over generations
- **`plot-ns.png`**: Neutral selection scenario
- **`plot-ns-qsl.png`**: Quality signal limit scenario

### Composite Graphical Abstract
- **`graphical_abstract.png`**: Comprehensive 8-panel figure combining all experiments and key insights

## Implementation Notes

- **Educational Purpose**: This is an educational implementation demonstrating the mathematical framework from Williams (2025)
- **Random Seed**: Set to 428 for reproducible results
- **Plot Format**: High-resolution PNG (200-300 DPI) for publication quality
- **Accessibility**: Large fonts (16-20pt), high contrast colors, clear legends
- **Memory Management**: Plots are closed after saving to prevent memory leaks
- **Vectorization**: All simulations use NumPy vectorization for efficiency
- **Grid Layouts**: Subtle grids for better data readability
- **Theoretical Fidelity**: Maintains mathematical accuracy while providing visual demonstrations

## Contributing

To extend these experiments:

1. Add new fitness functions in `model.py`
2. Create new plotting functions in `plotting.py`
3. Add experiment scenarios in `cli.py`
4. Update this README with new theoretical background

## License

This code is part of the METAINFORMANT project and follows the same licensing terms.
