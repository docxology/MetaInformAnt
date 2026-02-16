# Personal AI Infrastructure (PAI) - Root

## 🧠 Context & Intent

This is the root of the **METAINFORMANT** executable.

- **Intent**: Provide a comprehensive, modular, production-ready bioinformatics toolkit.
- **State**: As of Feb 2026, the repo is in a "Stabilized" state with high modularity and production-ready RNA-seq pipeline.

## 🏗️ Virtual Hierarchy

- `src/`: Core business logic.
- `scripts/`: Thin orchestration layer.
- `docs/`: User-facing documentation.
- `tests/`: No-mock verification suite.

## 📝 Maintenance Notes

- **Do not use mocks**: Always prefer real implementations.
- **Use `uv`**: Strict dependency management.
- **Output Isolation**: All artifacts go to `output/`.

## 🔄 AI Workflows

- **Refactoring**: When modifying `src/`, always run `scripts/run_all_scripts.py` to verify regressions.
- **Documentation**: Keep `AGENTS.md` and `README.md` synchronized with code changes.
