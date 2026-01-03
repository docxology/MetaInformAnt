# Simulation Visualization

This document provides comprehensive documentation for simulation data visualization capabilities in METAINFORMANT, including evolutionary simulations, population dynamics, agent-based models, and synthetic data analysis.

## Overview

Simulation visualization includes specialized plots for analyzing results from various simulation types, including sequence evolution, population genetics models, agent-based ecosystems, and RNA-seq data generation. These tools create publication-quality figures for simulation-based research.

## Module Functions

### Sequence Evolution

#### Evolution Trajectories
```python
from metainformant.simulation import visualization as sim_viz

# Plot sequence evolution over generations
sequence_history = [
    "ATCGATCGATCG",  # Generation 0
    "ATCGATCGATCA",  # Generation 1
    "ATCGATCGATCT",  # Generation 2
    # ... more generations
]

ax = sim_viz.plot_sequence_evolution(sequence_history, figsize=(12, 8))
```

#### Evolution Animation
```python
# Create animated sequence evolution
fig, anim = sim_viz.animate_sequence_evolution(sequence_history, interval=500)
```

### Population Dynamics

#### Population Trajectories
```python
# Plot population size changes over time
population_history = [
    [100, 50, 30],  # Generation 0: [pop1, pop2, pop3]
    [95, 55, 32],   # Generation 1
    [102, 48, 28],  # Generation 2
    # ... more generations
]

ax = sim_viz.plot_population_dynamics_simulation(population_history)
```

#### Population Animation
```python
# Animate population dynamics
fig, anim = sim_viz.animate_population_dynamics(population_history)
```

### RNA-seq Simulation

#### Expression Analysis
```python
# Plot RNA-seq simulation results
rnaseq_data = {
    'expression_values': np.random.lognormal(0, 1, 1000),
    'library_sizes': np.random.normal(1000000, 100000, 10),
    'mean_variance': {
        'means': np.random.lognormal(0, 1, 100),
        'variances': np.random.lognormal(0, 1, 100)
    }
}

ax = sim_viz.plot_rnaseq_simulation_results(rnaseq_data)
```

### Agent-Based Models

#### Agent Trajectories
```python
# Plot agent-based simulation results
agent_data = {
    'agent_states': [
        [0, 1, 0, 1, 2],  # Time 0
        [1, 1, 1, 0, 2],  # Time 1
        [1, 0, 1, 1, 1],  # Time 2
    ],
    'agent_positions': [
        [(0, 0), (1, 1), (2, 2)],  # Final positions
    ]
}

ax = sim_viz.plot_agent_based_model_results(agent_data)
```

### Evolutionary Simulations

#### Fitness Evolution
```python
# Plot comprehensive evolutionary simulation
evolution_data = {
    'fitness_history': {
        'generations': list(range(100)),
        'mean_fitness': np.random.normal(0.5, 0.1, 100) + np.linspace(0, 1, 100),
        'max_fitness': np.random.normal(0.7, 0.05, 100) + np.linspace(0, 1.5, 100)
    },
    'diversity_history': {
        'generations': list(range(100)),
        'diversity': np.random.normal(0.8, 0.1, 100) - np.linspace(0, 0.3, 100)
    },
    'allele_frequencies': {
        'generations': list(range(100)),
        'allele_A': np.random.beta(2, 5, 100),
        'allele_B': 1 - np.random.beta(2, 5, 100)
    }
}

ax = sim_viz.plot_evolutionary_simulation_summary(evolution_data)
```

### Parameter Sensitivity

#### Sensitivity Analysis
```python
# Plot parameter sensitivity results
sensitivity_data = {
    'parameter_values': np.linspace(0.1, 1.0, 20),
    'output_values': np.linspace(0.1, 1.0, 20) + np.random.normal(0, 0.1, 20)
}

ax = sim_viz.plot_simulation_parameter_sensitivity(sensitivity_data)
```

### Validation Analysis

#### Simulation vs Observed
```python
# Compare simulated vs observed data
observed_data = np.random.normal(0.5, 0.2, 100)
simulated_data = [
    np.random.normal(0.5, 0.2, 100) for _ in range(5)  # 5 replicates
]

ax = sim_viz.plot_simulation_validation_comparison(observed_data, simulated_data)
```

### Advanced Features

#### Interactive Dashboard
```python
# Create interactive simulation results dashboard
simulation_results = {
    'time_series': {
        'time': list(range(100)),
        'values': np.random.normal(0, 1, 100)
    },
    'distributions': np.random.normal(0, 1, 1000),
    'correlations': np.random.rand(10, 10),
    'statistics': {'mean': 0.0, 'std': 1.0, 'min': -3.0, 'max': 3.0}
}

fig = sim_viz.create_interactive_simulation_dashboard(simulation_results)
```

## Integration with Simulation Module

### With Sequence Simulations
```python
from metainformant.simulation import sequences, visualization as sim_viz

# Generate and evolve sequences
ancestor = sequences.generate_random_dna(100)
evolved = sequences.evolve_sequence(ancestor, generations=50)

# Visualize evolution
sequence_history = [ancestor] + [sequences.mutate_sequence(ancestor, i) for i in range(50)]
ax = sim_viz.plot_sequence_evolution(sequence_history)
```

### With Population Genetics
```python
from metainformant.simulation import popgen, visualization as sim_viz

# Simulate population evolution
simulation_results = popgen.simulate_bottleneck_population(1000, 100, 1000, 50)

# Extract population history for visualization
population_sizes = simulation_results.get('population_sizes', [])
ax = sim_viz.plot_population_dynamics_simulation(population_sizes)
```

### With Agent-Based Models
```python
from metainformant.simulation import agents, visualization as sim_viz

# Run agent-based simulation
ecosystem = agents.create_ecosystem(50, ['predator', 'prey'], (20, 20))
simulation_data = agents.run_simulation(ecosystem, 100)

# Visualize results
agent_states = [data['agent_states'] for data in simulation_data]
agent_data = {'agent_states': agent_states}
ax = sim_viz.plot_agent_based_model_results(agent_data)
```

## Output Options

All visualization functions support:
```python
# Save static plots
ax = sim_viz.plot_population_dynamics_simulation(data, output_path="dynamics.png")

# Save animations
fig, anim = sim_viz.animate_sequence_evolution(history, output_path="evolution.gif")

# Interactive dashboards
fig = sim_viz.create_interactive_simulation_dashboard(results, output_path="dashboard.html")
```

## Data Format Requirements

- **Sequence History**: List of strings, all same length
- **Population History**: List of arrays/lists with population sizes per generation
- **Agent Data**: Dictionary with 'agent_states', 'agent_positions', etc.
- **Evolution Data**: Dictionary with 'fitness_history', 'diversity_history', etc.
- **Time Series**: Dictionary with 'time' and 'values' arrays

## Performance Considerations

- **Large Simulations**: For >1000 generations, consider subsampling for plotting
- **Complex Animations**: Reduce frame rate for smoother playback
- **Memory Usage**: Large population histories can be memory-intensive
- **Interactive Plots**: Require Plotly for web-based dashboards

## Dependencies

- **Required**: matplotlib, numpy
- **Optional**: seaborn (enhanced styling), plotly (interactive plots)
- **Animation**: matplotlib.animation (built-in) or imageio for GIF export

## Examples

### Complete Simulation Analysis Workflow
```python
from metainformant.simulation import sequences, popgen, visualization as sim_viz
import numpy as np

# Run multiple types of simulations
fig, axes = plt.subplots(2, 2, figsize=(15, 12))

# Sequence evolution
ancestor = sequences.generate_random_dna(50)
sequence_history = [ancestor]
for gen in range(20):
    mutated = sequences.mutate_sequence(sequence_history[-1], 1)
    sequence_history.append(mutated)

sim_viz.plot_sequence_evolution(sequence_history, ax=axes[0,0])

# Population dynamics
initial_pop = [100, 50, 30]
population_history = [initial_pop]
for gen in range(30):
    # Simple population change model
    new_pop = [max(0, p + np.random.normal(0, 5)) for p in population_history[-1]]
    population_history.append(new_pop)

sim_viz.plot_population_dynamics_simulation(population_history, ax=axes[0,1])

# Evolutionary fitness
generations = list(range(50))
mean_fitness = 0.5 + np.cumsum(np.random.normal(0, 0.01, 50))
max_fitness = mean_fitness + np.random.normal(0.2, 0.05, 50)

evolution_data = {
    'fitness_history': {
        'generations': generations,
        'mean_fitness': mean_fitness,
        'max_fitness': max_fitness
    }
}
sim_viz.plot_evolutionary_simulation_summary(evolution_data, ax=axes[1,0])

# Parameter sensitivity
param_values = np.linspace(0.1, 2.0, 20)
output_values = param_values * 2 + np.random.normal(0, 0.2, 20)

sensitivity_data = {
    'parameter_values': param_values,
    'output_values': output_values
}
sim_viz.plot_simulation_parameter_sensitivity(sensitivity_data, ax=axes[1,1])

plt.tight_layout()
plt.savefig("simulation_analysis.png", dpi=300, bbox_inches='tight')
```

## Color Schemes

Recommended color schemes for simulation data:
```python
# Sequence evolution
mutation_colors = ['lightgray', 'red']  # No change, mutation

# Population dynamics
population_colors = plt.cm.tab10(np.linspace(0, 1, 10))  # Up to 10 populations

# Agent states
agent_colors = {
    0: 'blue',      # State 0
    1: 'red',       # State 1
    2: 'green',     # State 2
    3: 'orange'     # State 3
}

# Fitness evolution
fitness_colors = ['blue', 'red']  # Mean, maximum fitness
```

## Troubleshooting

### Common Issues

1. **Empty Plots**: Check data dimensions and types
2. **Animation Errors**: Ensure matplotlib.animation is available
3. **Memory Errors**: Reduce simulation size or subsample data
4. **Color Cycling**: Use explicit color maps for consistent visualization

### Data Validation

- **Sequence History**: All strings must be same length
- **Population Data**: Consistent number of populations across generations
- **Time Series**: Matching lengths for time and value arrays
- **Agent Data**: Proper dictionary structure with expected keys

## Related Documentation

- **[Simulation Analysis](../simulation/)**: Core simulation functions
- **[Sequence Simulation](../simulation/sequences.md)**: DNA/protein sequence generation
- **[Population Genetics](../simulation/popgen.md)**: Population simulation models
- **[Agent-Based Models](../simulation/agents.md)**: Ecosystem simulations
- **[Visualization Integration](integration.md)**: Cross-module visualization patterns

This module provides comprehensive simulation visualization capabilities integrated with METAINFORMANT's synthetic data generation and modeling workflows.
