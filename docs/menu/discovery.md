# Module Discovery System

The discovery system scans the repository `scripts/` directory to find Python and
bash scripts, extract their metadata, and generate navigable menu structures. It
powers the interactive menu by automatically cataloging every runnable script
along with its description, category, and argument requirements.

## Concepts

### Script Discovery

Discovery starts at the repository root and recursively walks `scripts/`. Each
`.py` and `.sh` file is inspected for metadata. Files inside `__pycache__`,
`archive`, or `test_examples` directories are skipped, as are files whose names
start with `_` (except `__init__.py`). The result is a dictionary mapping
category strings (derived from the first subdirectory under `scripts/`) to lists
of `ScriptInfo` objects sorted by name.

### Metadata Extraction

For Python scripts, the module docstring is parsed with `ast` and the first line
becomes the description. Argparse `add_argument` calls are walked to collect
required and optional argument names. Bash scripts are scanned for the first
meaningful comment line (length > 10 characters, skipping the shebang) whose
first sentence becomes the description.

### Menu Generation

`generate_menu_from_scripts` converts the categorized `ScriptInfo` dictionary
into a hierarchy of `Menu` objects. A root menu lists every category as a
submenu item. Each category menu lists its scripts as individual `MenuItem`
entries whose `action` field encodes the script path.

## Function Reference

### `discover_scripts(repo_root: Path) -> dict[str, list[ScriptInfo]]`

Walk `repo_root/scripts/` and return all discovered scripts grouped by category.

### `extract_script_metadata(script_path: Path) -> ScriptInfo`

Extract a `ScriptInfo` dataclass from a single script file. Delegates to
`_extract_python_metadata` or `_extract_bash_metadata` depending on the suffix.

### `categorize_script(script_path: Path) -> str`

Determine the category string for a script by inspecting its directory ancestry.
Returns the first path component after `scripts/`, or `"other"` when the script
sits directly in `scripts/` or inside an excluded directory.

### `generate_menu_from_scripts(scripts: dict[str, list[ScriptInfo]]) -> dict[str, Menu]`

Build a complete menu hierarchy from a discovery result. Returns a dictionary
keyed by menu ID (`"root"` for the top-level menu, `"menu_<category>"` for each
category submenu).

## Data Classes

### `ScriptInfo`

| Field           | Type         | Description                                 |
|-----------------|--------------|---------------------------------------------|
| `path`          | `Path`       | Absolute path to the script file            |
| `name`          | `str`        | Stem of the file name                       |
| `description`   | `str`        | First line of the docstring or comment       |
| `category`      | `str`        | Category derived from directory path         |
| `script_type`   | `str`        | `"python"` or `"bash"`                      |
| `required_args` | `list[str]`  | Argument names marked as required            |
| `optional_args` | `list[str]`  | Argument names that are optional             |

## Code Examples

```python
from pathlib import Path
from metainformant.menu.core.discovery import (
    discover_scripts,
    extract_script_metadata,
    generate_menu_from_scripts,
)

# Discover all scripts in the repository
repo = Path("/path/to/metainformant")
scripts = discover_scripts(repo)

for category, script_list in scripts.items():
    print(f"{category}: {len(script_list)} scripts")
    for info in script_list:
        print(f"  {info.name} ({info.script_type}) - {info.description}")
        if info.required_args:
            print(f"    required: {info.required_args}")

# Extract metadata from a single script
single = extract_script_metadata(repo / "scripts" / "rna" / "quantify.py")
print(single.category)  # "rna"

# Build navigable menus from discovered scripts
menus = generate_menu_from_scripts(scripts)
root_menu = menus["root"]
print(root_menu.title)  # "METAINFORMANT Scripts"
for item in root_menu.items:
    print(f"  {item.label}: {item.action}")
```

## Configuration

The discovery system has no external configuration files. Behavior is controlled
through the arguments passed to each function:

| Parameter  | Default | Effect                                            |
|------------|---------|---------------------------------------------------|
| `repo_root`| --      | Root directory containing the `scripts/` folder   |

Discovery filters are hard-coded: `archive` and `test_examples` directories are
excluded, `__pycache__` is skipped, and files prefixed with `_` (other than
`__init__.py`) are ignored.

## Import Path

```python
from metainformant.menu.core.discovery import discover_scripts
from metainformant.menu.core.discovery import ScriptInfo
```
