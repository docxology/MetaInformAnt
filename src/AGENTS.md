# AI Agents in Source Code

This document captures AI agent guidance for the METAINFORMANT source code directory.

## Role

The `src/` directory contains the main Python package source code. AI agents working here should understand the module organization and follow project conventions.

## Structure

```
src/
├── metainformant/          # Main package (see src/metainformant/AGENTS.md)
│   ├── core/               # Shared infrastructure
│   ├── dna/                # DNA sequence analysis
│   ├── rna/                # RNA-seq workflows
│   ├── protein/            # Protein analysis
│   ├── gwas/               # GWAS pipeline
│   └── ...                 # 15+ domain modules
└── metainformant.egg-info/ # Package metadata (auto-generated)
```

## Key Rules

### Package Management
- **ALWAYS use `uv`** for package operations
- Never use `pip` directly
- Editable install: `uv pip install -e .`

### NO MOCKING Policy
- All functions must use real implementations
- Never return dummy/placeholder data
- When external dependencies unavailable, raise errors or skip gracefully

### Core Utilities
Always use `metainformant.core` for:
- **I/O**: `from metainformant.core import io` (never `import json` directly)
- **Paths**: `from metainformant.core import paths`
- **Config**: `from metainformant.core import config`
- **Logging**: `from metainformant.core.logging import get_logger`

### Output Directory
- All execution outputs go to `output/`
- Never write to `src/` or create files in source directories

### Type Hints
- Python 3.11+ required
- All functions must have type hints
- Use `from __future__ import annotations`

## Module Navigation

Each module in `src/metainformant/` has its own:
- `README.md` - User documentation
- `AGENTS.md` - AI agent guidance
- `SPEC.md` - Technical specification
- `__init__.py` - Public API exports

## Related

- [Main Package](metainformant/AGENTS.md)
- [Tests](../tests/AGENTS.md)
- [Documentation](../docs/AGENTS.md)
