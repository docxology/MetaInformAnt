# AI Agent Guidance

This file is operational guidance for AI-assisted work in METAINFORMANT. Treat the checked-out code, tests, and domain documentation as the source of truth; do not preserve stale generated inventories when they disagree with the repo.

## Project Context

METAINFORMANT is a Python 3.11+ bioinformatics toolkit with modules under `src/metainformant/` for genomics, transcriptomics, GWAS, visualization, machine learning, ecology, ontology, and related domains.

Use the nearest nested `AGENTS.md` for folder-specific rules. Important entry points:

- `src/AGENTS.md` and `src/metainformant/AGENTS.md` for package work
- `docs/AGENTS.md` for documentation changes
- `scripts/AGENTS.md` for workflow scripts
- `config/AGENTS.md` for YAML/TOML/JSON configuration
- `tests/AGENTS.md` for test conventions

## Development Rules

- Use `uv` for dependency management and command execution.
- Keep runtime outputs under `output/` or a documented temporary/cache directory.
- Follow the real-implementation testing policy: exercise real implementations with small deterministic data.
- Prefer absolute imports from `metainformant` inside package code.
- Keep public APIs synchronized with `README.md`, domain docs, and `SPEC.md` files when behavior changes.

## Current Core Utility Paths

Use the current package layout in examples and new code:

- I/O helpers: `metainformant.core.io`
- Path helpers: `metainformant.core.io.paths` or re-exports from `metainformant.core.io`
- Configuration helpers: `metainformant.core.utils.config`
- Logging helpers: `metainformant.core.utils.logging`
- Validation helpers: `metainformant.core.data.validation`

Avoid obsolete imports such as the old `core.config` and `core.paths` module paths or the old top-level DNA `sequences` alias in new code. Compatibility shims for those paths are deliberately tested so existing imports keep working.

## Documentation Rules

- Documentation should describe the current checkout, not historical audit status or aspirational completion claims.
- Keep examples runnable or label them clearly as pseudocode.
- Cross-link to existing local files only; verify relative paths after moving docs.
- Keep README files concise and directory-specific. Use `docs/index.md` for broad navigation.
- Regenerate Cursor skills after adding, moving, or deleting `AGENTS.md` files:

  ```bash
  uv run python scripts/package/generate_cursor_skills.py
  uv run python scripts/package/generate_cursor_skills.py --check
  ```

## AI Assistance Disclosure

AI tools may be used for code, documentation, review, and test generation. Human maintainers remain responsible for validating correctness, scientific claims, dependency behavior, and reproducibility before results are used in research workflows.
