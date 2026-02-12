# Menu Core

Script discovery, metadata extraction, and execution engine for the interactive menu system.

## Contents

| File | Purpose |
|------|---------|
| `discovery.py` | Script scanning, categorization, metadata extraction, menu generation |
| `executor.py` | Script validation and execution for Python and Bash scripts |

## Key Classes and Functions

| Symbol | Description |
|--------|-------------|
| `ScriptInfo` | Dataclass holding script path, description, args, and dependencies |
| `discover_scripts()` | Scan repository for runnable scripts by category |
| `extract_script_metadata()` | Extract description, arguments, and deps from a script |
| `categorize_script()` | Classify script by domain (rna, gwas, simulation, etc.) |
| `generate_menu_from_scripts()` | Build Menu objects from discovered scripts |
| `validate_script_executable()` | Check script has execution permissions |
| `execute_script()` | Run Python or Bash script with arguments |
| `execute_python_script()` | Execute a Python script via subprocess |
| `execute_bash_script()` | Execute a Bash script via subprocess |

## Usage

```python
from metainformant.menu.core.discovery import discover_scripts, generate_menu_from_scripts
from metainformant.menu.core.executor import execute_script

scripts = discover_scripts(repo_root)
menus = generate_menu_from_scripts(scripts)
exit_code = execute_script(script_path, args=["--input", "data.csv"])
```
