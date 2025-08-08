### RNA: Overview

Thin, modular wrapper over the external `amalgkit` CLI with a plan/execute workflow.

```mermaid
flowchart LR
  A[Workflow Config] --> B[plan_workflow]
  B --> C[Ordered Steps]
  C --> D[execute_workflow]
  D --> E[amalgkit subcommands]
```

See: [amalgkit](./amalgkit.md), [Workflow](./workflow.md), [Configs](./configs.md).


