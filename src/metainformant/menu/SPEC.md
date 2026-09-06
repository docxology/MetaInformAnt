# Specification: menu

## Scope
Interactive menu and CLI interface for METAINFORMANT.

## Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## Data Structures
- **Sub-packages**: core, ui
- **Key Concepts**: Interactive CLI menus, configuration UI

## API Definition
### Exports
- `core.discovery` — `ScriptInfo`, `discover_scripts`, `generate_menu_from_scripts`, script metadata extraction
- `core.executor` — `validate_script_executable`, `execute_script`, `prompt_for_args`
- `ui.navigation` — `Menu`, `MenuItem`, `MenuHistory`, `MenuSystem`, navigation helpers
- `ui.display` — `format_menu`, `show_menu`, `get_choice`, `format_breadcrumb`, `clear_screen`
