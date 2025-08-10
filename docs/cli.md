# CLI

Entry: `uv run python -m metainformant` or `uv run metainformant`.

```text
uv run metainformant dna fetch --assembly GCF_000001405.40
uv run metainformant rna plan --work-dir output/amalgkit/work --threads 8 --species Apis_mellifera
uv run metainformant rna run  --work-dir output/amalgkit/work --threads 8 --species Apis_mellifera --check
uv run metainformant rna run-config --config config/amalgkit_pbarbatus.yaml --check
uv run metainformant tests -q
```

Subcommands

- **dna fetch**: validates assembly accessions (see [DNA Accessions](./dna/accessions.md))
- **rna plan**: prints an ordered plan of subcommands and parameters (see [RNA Workflow](./rna/workflow.md))
- **rna run**: executes the workflow; use `--check` to stop on first failure; logs written in `work-dir/logs` (default examples place this under `output/`)
- **rna run-config**: executes the workflow from a config file under `config/`; logs and manifest written under paths specified by the config
- **tests**: runs the repo tests (see [Testing](./testing.md))

```mermaid
sequenceDiagram
  participant U as User
  participant CLI as __main__.py
  participant DNA as dna/*
  participant RNA as rna/*
  U->>CLI: uv run metainformant rna plan --work-dir W
  CLI->>RNA: plan_workflow(config)
  RNA-->>CLI: steps: [(name, params)...]
  U->>CLI: uv run metainformant rna run --work-dir W --check
  CLI->>RNA: execute_workflow(config, check=True)
  RNA->>RNA: run_amalgkit(step, params)
  RNA-->>CLI: return codes
```

See: [RNA Workflow](./rna/workflow.md), [DNA](./dna/index.md).
