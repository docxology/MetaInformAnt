### RNA: Step Runners

Each `amalgkit` subcommand has a corresponding runner in `metainformant.rna.steps` for internal orchestration. All share the same signature:

```python
def run(params: Mapping[str, Any] | None = None, *, work_dir: str | Path | None = None, log_dir: str | Path | None = None, check: bool = False) -> CompletedProcess[str]
```

Available runners: `metadata`, `integrate`, `config`, `select`, `getfastq`, `quant`, `merge`, `cstmm`, `curate`, `csca`, `sanity`.

```python
from pathlib import Path
from metainformant.rna.steps import STEP_RUNNERS

runner = STEP_RUNNERS["quant"]
res = runner({"threads": 4}, work_dir=Path("./work"), log_dir=Path("./work/logs"), check=False)
```

Related: [Workflow](./workflow.md), [amalgkit](./amalgkit.md)


