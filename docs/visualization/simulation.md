# Simulation Visualization

This document provides comprehensive documentation for simulation data visualization capabilities in METAINFORMANT, including evolutionary simulations, population dynamics, agent-based models, and synthetic data analysis.

## Overview

The simulation package (`metainformant.simulation`) generates data and has no plotting module of its own. Results from the various simulation types, including sequence evolution, population genetics models, agent-based ecosystems, and RNA-seq data generation, are visualized with the general plotting functions of `metainformant.visualization` (basic and analysis plots, animations, composite and interactive figures). These tools create publication-quality figures for simulation-based research.

## Module Functions

### Sequence Evolution

#### Evolution Trajectories
```python
import numpy as np
from metainformant.simulation.models import sequences
from metainformant.visualization.plots.basic import lineplot

# Evolve a sequence one generation at a time
ancestor = sequences.generate_random_dna(100)
sequence_history = [ancestor]
for _ in range(20):
    sequence_history.append(sequences.evolve_sequence(sequence_history[-1], 1, mutation_rate=0.01))

# Plot similarity to the ancestor over generations
similarity = np.array([sequences.calculate_sequence_similarity(ancestor, s) for s in sequence_history])
ax = lineplot(similarity, label="Similarity to ancestor")
ax.set_xlabel("Generation")
ax.set_ylabel("Similarity")
```

#### Evolution Animation
```python
from metainformant.visualization.plots.animations import animate_evolution

# Create animated sequence evolution (sequence_history from above)
fig, anim = animate_evolution(sequence_history, interval=500)
```

### Population Dynamics

#### Population Trajectories
```python
import pandas as pd
from metainformant.visualization.analysis.timeseries import plot_time_series

# Plot population size changes over time
population_history = [
    [100, 50, 30],  # Generation 0: [pop1, pop2, pop3]
    [95, 55, 32],   # Generation 1
    [102, 48, 28],  # Generation 2
    # ... more generations
]

data = pd.DataFrame(population_history, columns=["pop1", "pop2", "pop3"])
ax = plot_time_series(data)
```

#### Population Animation
```python
import numpy as np
from metainformant.visualization.plots.animations import animate_time_series

# Animate population dynamics (one row per population, one column per generation)
fig, anim = animate_time_series(np.array(population_history).T)
```

### RNA-seq Simulation

#### Expression Analysis
```python
from metainformant.simulation.models import rna
from metainformant.visualization.plots.basic import scatter_plot

# Simulate RNA-seq counts (samples x genes) and plot the mean-variance relationship
counts = rna.simulate_rnaseq_counts(n_genes=1000, n_samples=10)
ax = scatter_plot(counts.mean(axis=0), counts.var(axis=0), s=5, alpha=0.5)
ax.set_xlabel("Mean count")
ax.set_ylabel("Variance")
```

### Agent-Based Models

#### Agent Trajectories
```python
import pandas as pd
from metainformant.visualization.analysis.timeseries import plot_time_series

# Plot agent-based simulation results: agents of each type at every time step
agent_counts = pd.DataFrame({
    'predator': [10, 12, 15, 11, 9],
    'prey': [40, 35, 28, 33, 38],
})

ax = plot_time_series(agent_counts)
```

### Evolutionary Simulations

#### Fitness Evolution
```python
import numpy as np
import pandas as pd
from metainformant.visualization.analysis.timeseries import plot_time_series
from metainformant.visualization.dashboards.composite import multi_panel

# Plot an evolutionary simulation: fitness and diversity over generations
generations = list(range(100))
fitness_history = pd.DataFrame({
    'mean_fitness': np.random.normal(0.5, 0.1, 100) + np.linspace(0, 1, 100),
    'max_fitness': np.random.normal(0.7, 0.05, 100) + np.linspace(0, 1.5, 100)
}, index=generations)
diversity_history = pd.DataFrame({
    'diversity': np.random.normal(0.8, 0.1, 100) - np.linspace(0, 0.3, 100)
}, index=generations)

fig = multi_panel([
    {"func": plot_time_series, "kwargs": {"data": fitness_history}, "title": "Fitness"},
    {"func": plot_time_series, "kwargs": {"data": diversity_history}, "title": "Diversity"},
])
```

### Parameter Sensitivity

#### Sensitivity Analysis
```python
import numpy as np
from metainformant.visualization.plots.basic import scatter_plot

# Plot parameter sensitivity results
parameter_values = np.linspace(0.1, 1.0, 20)
output_values = np.linspace(0.1, 1.0, 20) + np.random.normal(0, 0.1, 20)

ax = scatter_plot(parameter_values, output_values)
ax.set_xlabel("Parameter value")
ax.set_ylabel("Simulation output")
```

### Validation Analysis

#### Simulation vs Observed
```python
import numpy as np
from metainformant.visualization.analysis.statistical import histogram

# Compare simulated vs observed data
observed_data = np.random.normal(0.5, 0.2, 100)
simulated_data = [
    np.random.normal(0.5, 0.2, 100) for _ in range(5)  # 5 replicates
]

ax = histogram(observed_data, bins=20, alpha=0.6, label="Observed")
for replicate in simulated_data:
    histogram(replicate, bins=20, ax=ax, histtype="step", color="gray")
ax.legend()
```

### Advanced Features

#### Interactive Dashboard
```python
import numpy as np
from metainformant.visualization.interactive_dashboards.dashboards import (
    create_interactive_scatter,
    export_to_html,
)

# Create an interactive plot of simulation results (JSON-serializable data, standalone HTML export)
time = list(range(100))
values = np.random.normal(0, 1, 100).tolist()

plot_data = create_interactive_scatter(time, values, title="Simulated time series")
export_to_html(plot_data, "dashboard.html", title="Simulation results")
```

## Integration with Simulation Module

### With Sequence Simulations
```python-snippet
from metainformant.simulation.models import sequences
from metainformant.visualization.plots.animations import animate_evolution

# Generate and evolve sequences
ancestor = sequences.generate_random_dna(100)
evolved = sequences.evolve_sequence(ancestor, generations=50)

# Visualize evolution
sequence_history = [ancestor] + [sequences.mutate_sequence(ancestor, i) for i in range(50)]
fig, anim = animate_evolution(sequence_history)
```

### With Population Genetics
```python-snippet
import pandas as pd
from metainformant.simulation.models import popgen
from metainformant.visualization.analysis.timeseries import plot_time_series

# Simulate population evolution
simulation_results = popgen.simulate_bottleneck_population(1000, 100, 1000, 50)

# Extract population history for visualization
population_sizes = simulation_results.get('population_sizes', [])
ax = plot_time_series(pd.DataFrame({'population_size': population_sizes}))
```

### With Agent-Based Models
```python-snippet
import pandas as pd
from metainformant.simulation.models import agents
from metainformant.visualization.analysis.timeseries import plot_time_series

# Run agent-based simulation
ecosystem = agents.create_ecosystem(50, ['predator', 'prey'], (20, 20))
simulation_data = agents.run_simulation(ecosystem, 100)

# Visualize results: agents of each type at every recorded step
agent_counts = pd.DataFrame(
    [data['agent_counts'] for data in simulation_data],
    index=[data['step'] for data in simulation_data],
).fillna(0)
ax = plot_time_series(agent_counts)
```

## Output Options

All visualization functions support:
```python
# Save static plots
ax = plot_time_series(data, output_path="dynamics.png")

# Save animations
fig, anim = animate_evolution(history, output_path="evolution.gif")

# Interactive plots (standalone HTML)
export_to_html(plot_data, "dashboard.html")
```

## Data Format Requirements

- **Sequence History**: List of strings, all same length (`animate_evolution`)
- **Population History**: DataFrame with one column per population, indexed by generation (`plot_time_series`); 2D array with one row per population (`animate_time_series`)
- **Agent Data**: Snapshots from `agents.run_simulation` (`step`, `n_agents`, `agent_counts`, `total_energy`, `agent_positions`)
- **Evolution Data**: DataFrame of per-generation fitness, diversity, etc., indexed by generation
- **Time Series**: DataFrame indexed by time, one column per series

## Performance Considerations

- **Large Simulations**: For >1000 generations, consider subsampling for plotting
- **Complex Animations**: Reduce frame rate for smoother playback
- **Memory Usage**: Large population histories can be memory-intensive
- **Interactive Plots**: Require Plotly for web-based dashboards

## Dependencies

- **Required**: matplotlib, numpy
- **Optional**: seaborn (enhanced styling), plotly (interactive plots)
- **Animation**: matplotlib.animation (built-in); GIF export uses Pillow

## Examples

### Complete Simulation Analysis Workflow
```python-snippet
from metainformant.simulation.models import sequences
from metainformant.visualization.analysis.timeseries import plot_time_series
from metainformant.visualization.plots.basic import lineplot, scatter_plot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Run multiple types of simulations
fig, axes = plt.subplots(2, 2, figsize=(15, 12))

# Sequence evolution
ancestor = sequences.generate_random_dna(50)
sequence_history = [ancestor]
for gen in range(20):
    mutated = sequences.mutate_sequence(sequence_history[-1], 1)
    sequence_history.append(mutated)

similarity = np.array([sequences.calculate_sequence_similarity(ancestor, s) for s in sequence_history])
lineplot(similarity, ax=axes[0,0])

# Population dynamics
initial_pop = [100, 50, 30]
population_history = [initial_pop]
for gen in range(30):
    # Simple population change model
    new_pop = [max(0, p + np.random.normal(0, 5)) for p in population_history[-1]]
    population_history.append(new_pop)

plot_time_series(pd.DataFrame(population_history, columns=["pop1", "pop2", "pop3"]), ax=axes[0,1])

# Evolutionary fitness
generations = list(range(50))
mean_fitness = 0.5 + np.cumsum(np.random.normal(0, 0.01, 50))
max_fitness = mean_fitness + np.random.normal(0.2, 0.05, 50)

fitness_history = pd.DataFrame({
    'mean_fitness': mean_fitness,
    'max_fitness': max_fitness
}, index=generations)
plot_time_series(fitness_history, ax=axes[1,0])

# Parameter sensitivity
param_values = np.linspace(0.1, 2.0, 20)
output_values = param_values * 2 + np.random.normal(0, 0.2, 20)

scatter_plot(param_values, output_values, ax=axes[1,1])

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
- **[Agent-Based Models](../simulation/index.md#agent-based-ecosystem-modeling)**: Ecosystem simulations
- **[Visualization Integration](integration.md)**: Cross-module visualization patterns

These examples combine METAINFORMANT's synthetic data generation and modeling workflows with the general visualization functions.

