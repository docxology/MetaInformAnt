### CLI

Entry: `python -m metainformant` or the `metainformant` console script.

```text
metainformant dna fetch --assembly GCF_000001405.40
metainformant rna plan --work-dir ./work --threads 8 --species Apis_mellifera
metainformant rna run  --work-dir ./work --threads 8 --species Apis_mellifera --check
metainformant tests -q
```

```mermaid
sequenceDiagram
  participant U as User
  participant CLI as __main__.py
  participant DNA as dna/*
  participant RNA as rna/*
  U->>CLI: metainformant rna plan --work-dir W
  CLI->>RNA: plan_workflow(config)
  RNA-->>CLI: steps: [(name, params)...]
  U->>CLI: metainformant rna run --work-dir W --check
  CLI->>RNA: execute_workflow(config, check=True)
  RNA->>RNA: run_amalgkit(step, params)
  RNA-->>CLI: return codes
```

See: [RNA Workflow](./rna/workflow.md), [DNA](./dna/index.md).


