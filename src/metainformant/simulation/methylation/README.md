# methylation

Sub-package of `metainformant.simulation`. See the module's
`AGENTS.md`/`README.md` for the domain overview.

`simulator.py` generates synthetic DNA methylation datasets with ground-truth
differentially methylated regions (DMRs) and computes summary statistics:

| Function | Description |
|----------|-------------|
| `simulate_methylation()` | Simulate beta values for CpG islands, gene bodies, and promoters with injected DMRs |
| `calculate_dmr_statistics()` | Per-dataset summary statistics (site/sample counts, group mean betas) |

```python
from metainformant.simulation.methylation.simulator import (
    MethylationSimulationConfig,
    calculate_dmr_statistics,
    simulate_methylation,
)

dataset = simulate_methylation(MethylationSimulationConfig(random_seed=42))
stats = calculate_dmr_statistics(dataset)
```
