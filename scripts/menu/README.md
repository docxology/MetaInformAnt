# MetaInformAnt Interactive Menu

An interactive CLI menu system for navigating and launching the MetaInformAnt toolkit's workflows and scripts.

## Purpose

This module provides a user-friendly, text-based interface to:
- Browse available scripts across all biological domains (DNA, RNA, GWAS, etc.)
- Launch workflows with dependency checking
- Install optional dependencies via `uv` as needed

## Key Components

| File | Description |
|------|-------------|
| [run_menu.py](run_menu.py) | Main entry point with `main()` and interactive loop |

### Core Functions

- `check_uv_available()`: Verifies the `uv` package manager is in `PATH`.
- `check_optional_dependencies()`: Returns availability of each optional group.
- `install_optional_dependencies()`: Uses `uv` to install missing groups.
- `ensure_dependencies()`: Prompts the user and orchestrates installs.
- `run_interactive_menu()`: The main event loop.

## Usage

```bash
python3 scripts/menu/run_menu.py
```

## Related Documentation

- **Parent**: [scripts/README.md](../README.md) - Overview of all MetaInformAnt scripts.
- **Spec**: [SPEC.md](SPEC.md) - Technical specification for the menu system.
- **Agents**: [AGENTS.md](AGENTS.md) - AI development contributions.
- **Core Module**: [src/metainformant/core/README.md](../../src/metainformant/core/README.md)
