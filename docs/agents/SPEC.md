# Specification: agent-coordination-hub

**Path**: `docs/agents/`  
**Component Type**: Documentation Hub (Meta)  
**Dependency Level**: Cross-cutting (Core documentation layer)

## Scope

This specification describes the agent coordination documentation hub — the authoritative source for all multi-agent orchestration patterns, protocols, and best practices used throughout METAINFORMANT's 28-module ecosystem.

## Architecture

### Documentation Modules

| File | Role | Audience |
|------|------|----------|
| `README.md` | Hub overview with architecture diagram | Newcomers, architects |
| `index.md` | Toctree navigation | All users |
| `AGENTS.md` | Universal agent directives (rules) | All developers |
| `ARCHITECTURE.md` | System coordination architecture | System architects |
| `ORCHESTRATION.md` | Workflow manager API reference | Pipeline engineers |
| `MULTI_AGENT_WORKFLOWS.md` | Real-world workflow examples | Workflow designers |
| `COMMUNICATION_PROTOCOLS.md` | Inter-agent messaging patterns | Module developers |
| `SAFETY.md` | Error handling, validation, recovery | All implementers |
| `BEST_PRACTICES.md` | Operational excellence guide | Operators, maintainers |
| `rules/` | 28 module-specific rule files | Domain developers |

### Reference Implementation(s)

- **Primary**: `metainformant.core.engine.workflow_manager.BasePipelineManager`
- **Secondary**: `metainformant.core.execution.workflow` (config-driven)
- **Utilities**: `metainformant.core.execution.parallel`

## Data Structures

### Core Types (referenced by docs)

```python
from metainformant.core.engine.workflow_manager import (
    BasePipelineManager,
    PipelinePhase,
    PipelineItem,
    Stage,
)

# BasePipelineManager
manager: BasePipelineManager
manager.phases: List[PipelinePhase]
manager.items: Dict[str, PipelineItem]
manager.config: Dict[str, Any]
manager.executor: ThreadPoolExecutor
manager.ui: TerminalInterface

# PipelinePhase
phase = PipelinePhase(
    name: str,
    handler: Callable[[BasePipelineManager, List[PipelineItem]], None],
    filter_fn: Callable[[PipelineItem], bool] = lambda i: i.stage == Stage.PENDING,
    color: str = CYAN,
)

# PipelineItem
item = PipelineItem(
    item_id: str,
    metadata: Dict[str, Any] = field(default_factory=dict),
    stage: Stage = Stage.PENDING,
    error: str = "",
)

# Stage enum
Stage.PENDING → Stage.RUNNING → Stage.DONE | Stage.FAILED
```

## Integration

### Source Code

| Module | Purpose |
|--------|---------|
| `src/metainformant/core/engine/workflow_manager.py` | Pipeline engine |
| `src/metainformant/core/execution/workflow.py` | Config-driven execution |
| `src/metainformant/core/execution/parallel.py` | Parallel utilities |
| `src/metainformant/core/ui/tui.py` | Terminal UI |
| `src/metainformant/core/utils/logging.py` | Structured logging |
| `src/metainformant/core/io/*.py` | Atomic file I/O |

### Module Rules

- `docs/agents/rules/{module}.md` — Module-specific coding patterns
- Reside alongside: `src/metainformant/{module}/AGENTS.md`, `docs/{module}/AGENTS.md`

### Testing

- `tests/test_workflow_manager.py` — BasePipelineManager unit tests
- `tests/test_parallel.py` — Parallel execution tests
- `tests/NO_MOCKING_POLICY.md` — Zero mock policy enforced

## Interface (for Writers)

When extending documentation:

1. **New agent coordination pattern**: Add section to `MULTI_AGENT_WORKFLOWS.md`
2. **New safety protocol**: Append to `SAFETY.md` with example
3. **Orchestration API change**: Update `ORCHESTRATION.md` and `AGENTS.md` directive set
4. **New module**: Create `docs/agents/rules/{module}.md` and cross-link

### Cross-Linking Conventions

Internal links within `docs/agents/`:

```markdown
[AGENTS.md](AGENTS.md)                # sibling
[Orchestration](ORCHESTRATION.md)      # implicit .md
[rules/core.md](rules/core.md)        # subdirectory
[../core/](../core/)                  # sibling of docs/ (go up one)
```

External cross-module links:

```markdown
[Core Infrastructure](../core/README.md)
[RNA Workflow](../rna/README.md)
[Workflow Manager API](../src/metainformant/core/engine/workflow_manager.py)
```

## Testing Policy

- **Zero Mock**: All examples and doctests must use real implementations (see `tests/NO_MOCKING_POLICY.md`)
- **Runnable code**: Example snippets must be executable as-is
- **Type-checked**: All code snippets follow project's type hint conventions

## Versioning

Documentation version aligns with project version:

```
docs/agents/ — matches metainformant.__version__
```

## Change Log

| Change | Date | File(s) |
|--------|------|---------|
| Initial comprehensive hub creation | 2025-04-27 | All new files |
| Migration from minimal rules → full hub | 2025-04-27 | README.md, index.md, AGENTS.md + 6 new |

## See Also

- [Root AGENTS.md](../../AGENTS.md) — Project-wide AI agent philosophy
- [Core Module SPEC](../core/SPEC.md) — Core utilities specification
- [RNA Module Specification](../rna/SPEC.md) — Example of complex orchestration in action
- [Documentation Improvement Plan](../../docs/DOCUMENTATION_IMPROVEMENT_PLAN.md) — Audit history
