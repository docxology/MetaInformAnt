# AGENTS: Menu System

AI-assisted implementation of the interactive CLI menu for MetaInformAnt.

## AI Contributions

The **Code Assistant Agent** developed:
- The `MenuSystem` class architecture for hierarchical navigation.
- Dynamic dependency detection using the optional groups from `pyproject.toml`.
- `uv`-based interactive installation prompts.
- Error recovery and graceful exit handling.

## Function Index

| Function | Signature |
|----------|-----------|
| `check_uv_available` | `() -> bool` |
| `check_optional_dependencies` | `() -> dict[str, bool]` |
| `install_optional_dependencies` | `(missing_groups: list[str]) -> bool` |
| `ensure_dependencies` | `() -> bool` |
| `run_interactive_menu` | `(menu_system: MenuSystem) -> None` |
| `handle_menu_action` | `(item, menu_system: MenuSystem) -> str | None` |
| `main` | `() -> int` |

## Related Documentation

- **README**: [README.md](README.md) - Usage guide.
- **SPEC**: [SPEC.md](SPEC.md) - Technical specification.
- **Parent**: [scripts/AGENTS.md](../AGENTS.md)
