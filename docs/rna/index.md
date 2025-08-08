### RNA: Overview

Thin, modular wrapper over the external `amalgkit` CLI with a plan/execute workflow.

```mermaid
flowchart LR
  A[Workflow Config] --> B[plan_workflow]
  B --> C[Ordered Steps]
  C --> D[execute_workflow]
  D --> E[amalgkit subcommands]
  E --> F[Logs + Manifest]
```

See: [amalgkit](./amalgkit.md), [Workflow](./workflow.md), [Configs](./configs.md).
Advanced: step runners in `metainformant.rna.steps` provide a stable call surface per subcommand.
See: [RNA Steps](./steps.md)

Notes
- Ensure `amalgkit` is installed and on PATH; verify using `check_cli_available`
- `execute_workflow` writes per-step logs to `work_dir/logs` and a JSONL manifest to `work_dir/amalgkit.manifest.jsonl`
- To parameterize steps by species/tissues and output layout, see [Configs](./configs.md)


