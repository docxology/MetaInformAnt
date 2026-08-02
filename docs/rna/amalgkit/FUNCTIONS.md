# Amalgkit function index

The public wrappers live in `metainformant.rna.amalgkit.amalgkit` and return
`subprocess.CompletedProcess` records with captured output and return codes.

## Current command wrappers

| Wrapper | CLI command | Scope |
|---|---|---|
| `metadata` | `metadata` | per species |
| `select` | `select` | per species |
| `getfastq` | `getfastq` | per sample/cohort |
| `integrate` | `integrate` | per species |
| `quant` | `quant` | per sample/cohort |
| `merge` | `merge` | per species |
| `wsfilter` | `wsfilter` | per species |
| `cstmm` | `cstmm` | optional cross species |
| `csfilter` | `csfilter` | optional cross species |
| `finalize` | `finalize` | per species |
| `sanity` | `sanity` | per species |

## Workflow helpers

- `metainformant.rna.engine.workflow.plan_workflow`
- `metainformant.rna.engine.workflow.execute_workflow`
- `metainformant.rna.engine.workflow.validate_workflow_config`
- `metainformant.rna.engine.workflow.validate_workflow_outputs`
- `metainformant.rna.engine.workflow_steps.validate_step_prerequisites`
- `metainformant.rna.engine.pipeline.summarize_finalize_tables`

## Command verification

```bash
uv run python - <<'PY'
from metainformant.rna.amalgkit import validate_amalgkit_version
print(validate_amalgkit_version())
PY
amalgkit --help
```

The installed command help and `src/metainformant/rna/amalgkit/commands.py`
must agree before a new option is documented.
