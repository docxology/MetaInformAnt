# Agent Directives: rna

**Context**: RNA transcriptomic analysis and workflow orchestration module for METAINFORMANT.

## Capabilities

End-to-end RNA-seq pipelines: ENA/SRA metadata retrieval, FASTQ streaming, Kallisto quantification, cross-species TMM normalization, and industrial-scale orchestration (8,300+ samples across 28 Hymenoptera species).

## Subpackages

| Subpackage | Purpose |
|------------|---------|
| `amalgkit/` | Amalgkit CLI integration — command construction, workflow config, step execution |
| `engine/` | Workflow orchestration — `AmalgkitWorkflowConfig`, `execute_workflow()`, `plan_workflow()` |
| `retrieval/` | ENA/SRA data retrieval and metadata download |
| `analysis/` | Expression matrix analysis and quality assessment |
| `core/` | Shared RNA utilities and data structures |
| `splicing/` | Alternative splicing analysis |
| `deconvolution/` | Cell-type deconvolution from bulk RNA-seq |

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for domain data file I/O. Direct stdlib parsing is allowed in core, protocol adapters, subprocess/CLI glue, and narrow parser internals when covered by tests.
- Prefer absolute imports from `metainformant`
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management
- Use `build_cli_args()` / `build_amalgkit_command()` to test command construction without subprocesses
- Use `tmp_path` and real YAML configs for orchestrator tests
- Use real `AmalgkitWorkflowConfig`, `WorkflowExecutionResult`, and `WorkflowStepResult` classes

See `tests/NO_MOCKING_POLICY.md` for the full policy.

## Related Documentation

- **Module guide**: [../../../docs/rna/](../../../docs/rna/) — In-depth usage, architecture, and examples
- **API reference**: [SPEC.md](SPEC.md) — Type signatures, data structures, error codes
- **Core infrastructure**: [../core/AGENTS.md](../core/AGENTS.md) — Shared utilities (logging, config, I/O)
- **Full module index**: [../../../docs/index.md](../../../docs/index.md) — Overview of all METAINFORMANT modules
- **Dna module**: [../dna/AGENTS.md](../dna/AGENTS.md)
- **Quality module**: [../quality/AGENTS.md](../quality/AGENTS.md)
- **Visualization module**: [../visualization/AGENTS.md](../visualization/AGENTS.md)
- **Core module**: [../core/AGENTS.md](../core/AGENTS.md)
