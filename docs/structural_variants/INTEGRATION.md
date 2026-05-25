# Integration: structural_variants

## Related Modules
- [dna](../dna/README.md)
- [longread](../longread/README.md)
- [dna/variants](../dna/variants.md)
- [visualization](../visualization/README.md)
- [gwas](../gwas/README.md)

## Import Patterns
### Minimal (standalone)
```python
from metainformant import structural_variants
results = structural_variants.analyze(data)
```

### Full integration
```python
from metainformant.structural_variants import detect_cnvs, call_svs, annotate_svs
from metainformant.core import io, config, cache, logging
from metainformant.visualization import quickplot
```

## Data Flow Architecture
```mermaid
flowchart LR
    A[Raw Input
BAM/CRAM (alignments)] --> B[I/O
structural_variants.io]
    B --> C[Core
structural_variants.core]
    C --> D[Results
dict/DataFrame]
    D --> E[Visualization]
    D --> F[Export
JSON/CSV]
    E --> G[Plot]
    F --> H[File]
    
    style C fill:#e1f5e1
    style D fill:#fff3e0
```

## Cross-Module Workflows

### structural_variants + gwas + phenotype
```python
from metainformant import structural_variants, gwas, phenotype
prep = structural_variants.detect_cnvs(raw_data)
gwas_res = gwas.associate(prep)
pheno = phenotype.correlate(gwas_res)
visualization.manhattan(gwas_res)
```

### structural_variants + cloud
```python
from metainformant.cloud import submit_batch
job = submit_batch(
    module="structural_variants",
    parameters=dict(algorithm="accurate", workers=16)
)
results = job.wait().download()
```

### structural_variants + visualization (detailed)
```python
from metainformant.structural_variants import analyze
from metainformant.visualization import plot_heatmap, plot_timeseries, plot_network
res = analyze(data)
plot_heatmap(res.matrix)
plot_timeseries(res.series)
```

## Shared Infrastructure
| Service | Provided By | Used By |
|---------|-------------|---------|
| Config | `core.utils.config` | All modules (including structural_variants) |
| Logging | `core.utils.logging` | Structured logs across pipeline |
| I/O | `core.io` | Format-agnostic read/write |
| Caching | `core.io.cache` | Expensive computations |
| DB | `core.db` | Persistent metadata (optional) |

## Data Contract
Structural_variants produces standardized result containers compatible with downstream modules:
```python
class Result:
    values: dict            # Primary output data
    stats: dict             # Summary statistics
    metadata: dict          # Provenance, timestamps, config
    to_table() -> pd.DataFrame
    to_json() -> str
```
All modules that accept `Result` objects can operate on structural_variants's output directly.

## Multi-Module Orchestration
See [docs/agents/](../agents/MULTI_AGENT_WORKFLOWS.md) for agent-driven multi-module pipelines.
