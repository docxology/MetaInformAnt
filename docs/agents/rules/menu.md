# Menu Module Rules

## Purpose
Interactive command-line menu system for discovering and executing METAINFORMANT analysis workflows via script discovery and categorization.

## Source Structure
```text
src/metainformant/menu/
  app.py              # MenuApp, run_menu_app: interactive loop and action dispatch
  core/
    discovery.py      # ScriptInfo, script discovery, metadata extraction, categorization, menu generation
    executor.py       # Script execution (Python/Bash), argument prompting, validation
  ui/
    display.py        # Menu formatting, choice input, breadcrumb formatting, screen clearing
    navigation.py     # Menu data structures (MenuItem, Menu, MenuSystem, MenuHistory) and navigation
```

## Dependencies
- **Required**: Standard library only (`ast`, `dataclasses`, `pathlib`, `subprocess`, `sys`, `typing`)
- The menu imports no other `metainformant` package, including `metainformant.core`

## Import Patterns
```python
from metainformant.menu.core import discovery, executor
from metainformant.menu.ui import display, navigation
from metainformant.menu.app import MenuApp, run_menu_app
```

## Configuration
- No environment prefix needed (internal tool); the menu reads no environment variables
- `discovery.discover_scripts(repo_root)` scans `<repo_root>/scripts/` recursively for `*.py` and `*.sh` files and groups them by category
- The category is the first path component after `scripts/` (normally the domain directory); scripts under `archive/` or `test_examples/` get the category `other`
- It skips `__pycache__` paths, Python files whose name contains `test_`, and files whose name starts with `_` (an `__init__.py` is not skipped)

## Execution
- `executor.execute_script(path, args=None)` runs `.py` files with `sys.executable` and `.sh` files with `bash`, in the script's own directory, and returns the exit code (`1` for a missing file or an unsupported suffix)
- `execute_bash_script` sets the script's mode to `0o755` before running it
- `validate_script_executable` accepts any existing `.py` or `.sh` file; it does not check the executable bit

## Integration
- **Menu → scripts**: discovers and launches the scripts in any `scripts/<category>/` directory, whichever module the category belongs to
- **Menu → Core**: none; execution uses `subprocess` directly, not `metainformant.core`

## Testing
- Test script discovery against real script files written under `tmp_path` (see `tests/menu/test_menu_discovery.py`)
- Test menu navigation logic with real MenuItem structures
- All test outputs to `tmp_path`
