# Simulation Visualization

Plotting and animation functions for simulation results including sequence evolution, RNA-seq, population dynamics, and agent-based models.

## Contents

| File | Purpose |
|------|---------|
| `visualization.py` | Plots and animations for all simulation model types |

## Key Functions

| Function | Description |
|----------|-------------|
| `plot_sequence_evolution()` | Visualize sequence divergence over evolutionary time |
| `animate_sequence_evolution()` | Animated mutation accumulation along a phylogeny |
| `plot_rnaseq_simulation_results()` | Heatmaps and distributions for simulated RNA-seq |
| `plot_population_dynamics_simulation()` | Population size trajectories over generations |
| `plot_agent_based_model_results()` | Agent counts, spatial maps, and interaction plots |
| `plot_evolutionary_simulation_summary()` | Multi-panel summary of evolutionary simulation |
| `plot_simulation_parameter_sensitivity()` | Parameter sweep sensitivity analysis |
| `animate_population_dynamics()` | Animated population size changes |
| `plot_simulation_validation_comparison()` | Simulated vs observed data comparison |
| `create_interactive_simulation_dashboard()` | Interactive Plotly dashboard for exploration |

## Usage

```python
from metainformant.simulation.visualization.visualization import (
    plot_rnaseq_simulation_results,
    plot_population_dynamics_simulation,
)

plot_rnaseq_simulation_results(counts, output_path="output/sim_rna.png")
plot_population_dynamics_simulation(trajectories, output_path="output/sim_pop.png")
```
