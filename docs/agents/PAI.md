# Personal AI Infrastructure (PAI) - Agent Coordination Hub

**Path**: `docs/agents/`  
**Purpose**: Comprehensive documentation for METAINFORMANT agent coordination  
**Domain**: Documentation (Meta)  
**Parent**: `docs/`

## Overview

This directory serves as the authoritative hub for agent coordination in METAINFORMANT. It encompasses AI-assisted development guidelines, orchestration patterns, multi-agent workflows, and safety protocols that govern the 28-module bioinformatics platform.

## Virtual Hierarchy

- **Type**: Documentation Hub
- **Scope**: System-wide coordination across all modules
- **Layer**: Core documentation layer (cross-cutting)

## Maintenance Notes

| Aspect | Detail |
|--------|--------|
| **System** | Part of METAINFORMANT Core infrastructure (affects all 28 modules) |
| **Style** | Strict type hints enforced, zero mocks in tests, `uv` package management |
| **Stability** | Coordinator APIs (BasePipelineManager) stable; workflows evolve |
| **Scope** | Every module participates in multi-agent workflows to some degree |

## Documentation Structure

```text
docs/agents/
  README.md                    → Hub overview & architecture diagram
  index.md                     → Toctree navigation
  AGENTS.md                    → Universal directives for all agents
  ARCHITECTURE.md              → System coordination architecture
  ORCHESTRATION.md             → Workflow manager API reference
  MULTI_AGENT_WORKFLOWS.md     → Real-world example pipelines
  COMMUNICATION_PROTOCOLS.md   → Inter-agent messaging
  SAFETY.md                    → Error handling, validation, recovery
  BEST_PRACTICES.md            → Configuration & operations
  SPEC.md                      → Specification: cursorrules meta-documentation
  rules/                       → Module-specific coding rules
 index.md
 core.md
 rna.md
 dna.md
 {28 modules}
```

## AI Workflows (Developer)

### Modification Workflow

1. **Edit documentation** in `docs/agents/` or corresponding module
2. **Run functional tests** to ensure cross-links work:

   ```bash
   pytest tests/test_documentation.py
   # or
   uv run python scripts/verify_docs.py
   ```

3. **Update toctrees** if new `.md` files added
4. **Commit** with descriptive message linking to relevant issue

### When to Update Which File

| Change Type | File(s) to Update |
|-------------|-------------------|
| New coordination pattern | `MULTI_AGENT_WORKFLOWS.md` |
| New safety finding | `SAFETY.md` (and cross-link) |
| Modified orchestration API | `ORCHESTRATION.md` + `AGENTS.md` |
| New module added | `rules/{module}.md` + `index.md` toctree |
| Cross-module protocol change | `COMMUNICATION_PROTOCOLS.md` + affected module AGENTS.md |

### Documentation Updates

- **Architectural patterns change**: Update `ARCHITECTURE.md`
- **New workflow example**: Add to `MULTI_AGENT_WORKFLOWS.md`
- **Protocol change**: Update `COMMUNICATION_PROTOCOLS.md` and relevant module docs
- **Config schema change**: Update `BEST_PRACTICES.md` and module SPEC

## Cross-References

- [Agent Directives](AGENTS.md) — Rules all agents must follow
- [Architecture](ARCHITECTURE.md) — Coordination system design
- [Orchestration API](ORCHESTRATION.md) — Workflow manager reference
- [Workflow Examples](MULTI_AGENT_WORKFLOWS.md) — Practical templates
- [Communication](COMMUNICATION_PROTOCOLS.md) — Metadata and config patterns
- [Safety](SAFETY.md) — Error handling, validation, recovery
- [Best Practices](BEST_PRACTICES.md) — Operational excellence
- [Module Rules](rules/index.md) — Domain-specific coding conventions

## Integration Points

- **Source coordination code**: `src/metainformant/core/engine/workflow_manager.py`
- **Parallel execution**: `src/metainformant/core/execution/parallel.py`
- **Config workflows**: `src/metainformant/core/execution/workflow.py`
- **Module AGENTS.md**: Per-module rules in `src/*/AGENTS.md` and `docs/*/AGENTS.md`
- **No-mocking policy**: `tests/NO_MOCKING_POLICY.md` (enforced in tests)

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0.0 | 2025-04-27 | Initial comprehensive coordination hub (Hermes Agent) |

---

*This PAI document is part of the METAINFORMANT Core layer. Updates here affect all modules.*
