# Phenotype Workflow

Pipeline orchestration for multi-step phenotype analysis across behavioral, chemical, morphological, electronic, and sonic domains. Supports YAML-driven configuration and structured result collection.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Re-exports `pipeline` module |
| `pipeline.py` | Pipeline configuration, execution engine, and result aggregation |

## Key Classes and Functions

| Symbol | Description |
|--------|-------------|
| `PipelineConfig` | Dataclass for pipeline settings (steps, types, paths, parameters) |
| `PipelineConfig.from_yaml()` | Load pipeline configuration from a YAML file |
| `PipelineResult` | Dataclass capturing step-by-step execution results |
| `PhenotypePipeline` | Main pipeline executor with load/validate/analyze/summarize steps |

## Usage

```python
from metainformant.phenotype.workflow.pipeline import PipelineConfig, PhenotypePipeline

config = PipelineConfig.from_yaml("config/phenotype_pipeline.yaml")
pipeline = PhenotypePipeline(config)
result = pipeline.run()
```
