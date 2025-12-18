# Mathematical Biology Scripts

Mathematical modeling and theoretical biology workflow orchestrators.

## Directory Structure

```
scripts/math/
├── run_math_modeling.py           # Mathematical biology workflow orchestrator
└── README.md                      # This file
```

## Mathematical Biology Modeling (`run_math_modeling.py`)

Comprehensive mathematical biology workflow orchestrator for population dynamics, epidemiology, evolutionary models, and theoretical analysis.

**Features:**
- SIR/SEIR epidemiology models
- Logistic and exponential growth models
- Natural and artificial selection simulations
- Price equation analysis
- Population genetics models
- Stability and bifurcation analysis

**Usage:**
```bash
# SIR epidemiology model
python3 scripts/math/run_math_modeling.py --model sir --beta 0.3 --gamma 0.1 --S0 990 --I0 10 --steps 100

# Logistic growth model
python3 scripts/math/run_math_modeling.py --model logistic --r 3.5 --x0 0.5 --steps 100

# Selection experiment simulation
python3 scripts/math/run_math_modeling.py --model selection --generations 50 --n 10000 --s_hat 0.6

# Price equation analysis
python3 scripts/math/run_math_modeling.py --model price --fitness 0.2,0.4,0.1 --trait-parent 1.0,1.2,0.9 --trait-offspring 1.25,1.35,0.95
```

**Models Available:**
- `sir`: Susceptible-Infected-Recovered epidemiology
- `seir`: Susceptible-Exposed-Infected-Recovered epidemiology
- `logistic`: Logistic population growth
- `exponential`: Exponential population growth
- `selection`: Natural/artificial selection simulation
- `price`: Price equation evolutionary analysis
- `hardy_weinberg`: Hardy-Weinberg equilibrium analysis

**Options:**
- `--model`: Mathematical model to run
- `--output`: Output directory (defaults to output/math/)
- `--steps`/`--generations`: Number of simulation steps
- Model-specific parameters (beta, gamma, r, etc.)
- `--threads`: Number of threads for parallel simulations
- `--verbose`: Enable verbose logging

**Output Structure:**
```
output/math/
├── [model_name]_simulation/
│   ├── parameters.json            # Simulation parameters
│   ├── results.json               # Simulation results data
│   ├── equilibrium_analysis.json  # Stability analysis
│   ├── bifurcation_analysis.json  # Bifurcation analysis (if applicable)
│   ├── phase_portraits/           # Generated visualizations
│   │   ├── time_series.png
│   │   ├── phase_portrait.png
│   │   ├── bifurcation_diagram.png
│   │   └── stability_analysis.png
│   └── analysis_report.json       # Comprehensive analysis report
```

## Integration

Integrates with:
- **metainformant.math**: Core mathematical biology functionality
- **Scientific computing**: NumPy, SciPy, SymPy
- **Visualization**: matplotlib, seaborn for plotting
- **Core utilities**: I/O, logging, path management

## Dependencies

- **metainformant.math**: Mathematical biology module
- **NumPy/SciPy**: Numerical computing and ODE solving
- **SymPy**: Symbolic mathematics
- **matplotlib/seaborn**: Visualization support

## Related Documentation

- [Mathematical Biology Documentation](../../docs/math/README.md)
- [Population Dynamics](../../docs/math/population_dynamics.md)
- [Evolutionary Models](../../docs/math/evolutionary_models.md)
- [METAINFORMANT CLI](../../docs/cli.md)

