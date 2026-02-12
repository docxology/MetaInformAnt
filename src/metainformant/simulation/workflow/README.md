# Simulation Workflow

End-to-end simulation workflow orchestration with configuration, execution, validation, and calibration for all simulation types.

## Contents

| File | Purpose |
|------|---------|
| `workflow.py` | Simulation config, runner, validator, and calibration pipeline |

## Key Classes and Functions

| Symbol | Description |
|--------|-------------|
| `SimulationConfig` | Dataclass holding simulation type, parameters, seed, and output path |
| `create_simulation_config()` | Factory for typed simulation configurations |
| `run_simulation_workflow()` | Execute complete simulation pipeline from config |
| `run_benchmark_simulation()` | Run simulation with timing and performance metrics |
| `validate_simulation_output()` | Check simulation output against validation criteria |
| `calibrate_simulation_parameters()` | Optimize parameters to match target statistics |

## Usage

```python
from metainformant.simulation.workflow.workflow import (
    create_simulation_config,
    run_simulation_workflow,
    validate_simulation_output,
)

config = create_simulation_config("population_genetics", {"n_populations": 3})
results = run_simulation_workflow(config)
valid, issues = validate_simulation_output(results, {"min_samples": 100})
```
