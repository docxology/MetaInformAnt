# Agent Directives: rna

## 🤖 Role

Specialized agent context for the `rna` component.

## 🛠️ Tools & Capabilities

- **Context**: RNA transcriptomic analysis and workflow orchestration module for METAINFORMANT.
- **Pattern**: Source Code Pattern

## ⚠️ Rules & Constraints

- **Imports**: Prefer absolute imports from `metainformant`.
- **I/O**: Use `metainformant.core.io` for all file operations.
- **Logging**: Use `metainformant.core.utils.logging`.

## 📋 Code Quality Policy

All code in this module **MUST** be:

1. **Functional** — Every function performs real work with real data. No placeholder implementations.
2. **Modular** — Clear separation of concerns across `engine/`, `amalgkit/`, `retrieval/`, `analysis/`, `splicing/`.
3. **Tested** — Comprehensive tests using real implementations, real configs, and real filesystem.
4. **Documented** — README.md, AGENTS.md, SPEC.md at every directory level. Docstrings on all public functions.

## 🚫 NO_MOCKING_POLICY (Mandatory)

> **NEVER use `unittest.mock`, `pytest-mock`, `MagicMock`, `patch`, or any mocking library.**

Instead:

- Use `build_cli_args()` / `build_amalgkit_command()` to test command construction without running subprocesses.
- Use `tmp_path` and real YAML configs for orchestrator tests.
- Use real `AmalgkitWorkflowConfig`, `WorkflowExecutionResult`, and `WorkflowStepResult` classes.
- Use real `pd.DataFrame` for data pipeline tests.
- Use `try/except` for optional network tests (ENA API).

See `tests/NO_MOCKING_POLICY.md` for the full policy.
