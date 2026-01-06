# Menu Module

The `menu` module provides an interactive text-based menu system for accessing and executing scripts throughout the METAINFORMANT repository. It includes menu display, navigation, script discovery, and execution capabilities.

## Overview

The menu module enables users to:
- Browse available scripts through interactive menus
- Navigate hierarchical menu structures with breadcrumbs
- Automatically discover scripts in the `scripts/` directory
- Execute scripts with argument prompting
- Access scripts via the `metainformant.sh` orchestrator

## Architecture

```mermaid
graph TB
    subgraph "Menu Module"
        Display[display.py<br/>Menu Rendering]
        Nav[navigation.py<br/>Navigation State]
        Disc[discovery.py<br/>Script Discovery]
        Exec[executor.py<br/>Script Execution]
    end
    
    subgraph "User Interface"
        Orchestrator[metainformant.sh<br/>Orchestrator]
        CLI[Command Line]
    end
    
    subgraph "Scripts"
        Scripts[scripts/*/<br/>Domain Scripts]
    end
    
    CLI --> Orchestrator
    Orchestrator --> Display
    Orchestrator --> Disc
    Display --> Nav
    Nav --> Exec
    Disc --> Scripts
    Exec --> Scripts
```

## Components

### Display Module (`display.py`)

Functions for rendering menus and handling user input:

- `format_menu(items, title="", width=80) -> str`: Format menu items into displayable string
- `show_menu(items, title="") -> None`: Display menu to user
- `get_choice(prompt, valid_choices=None) -> str`: Get validated user input
- `format_breadcrumb(path) -> str`: Format navigation breadcrumb
- `clear_screen() -> None`: Clear terminal screen

### Navigation Module (`navigation.py`)

Classes and functions for menu navigation:

- `MenuItem`: Single menu item with action and metadata
- `Menu`: Menu container with items
- `MenuSystem`: Navigation system with state management
- `MenuHistory`: Breadcrumb tracking
- `navigate_to_submenu(menu_system, submenu_id) -> None`: Navigate to submenu
- `go_back(menu_system) -> bool`: Navigate back in history
- `get_current_menu(menu_system) -> Menu`: Get current menu

### Discovery Module (`discovery.py`)

Script discovery and menu generation:

- `discover_scripts(repo_root) -> dict[str, list[ScriptInfo]]`: Discover all scripts
- `extract_script_metadata(script_path) -> ScriptInfo`: Extract metadata from script
- `categorize_script(script_path) -> str`: Categorize script by directory
- `generate_menu_from_scripts(scripts) -> dict[str, Menu]`: Generate menu structure
- `ScriptInfo`: Script metadata dataclass

### Executor Module (`executor.py`)

Script execution utilities:

- `execute_script(script_path, args=None) -> int`: Execute script (Python or bash)
- `execute_python_script(script_path, args=None) -> int`: Execute Python script
- `execute_bash_script(script_path, args=None) -> int`: Execute bash script
- `prompt_for_args(script_info) -> list[str]`: Prompt user for script arguments
- `validate_script_executable(script_path) -> bool`: Validate script exists

## Usage

### Interactive Menu

Launch the interactive menu system:

```bash
./metainformant.sh
# or
./metainformant.sh --menu
```

### Programmatic Usage

```python
from pathlib import Path
from metainformant.menu import discover_scripts, generate_menu_from_scripts
from metainformant.menu.display import show_menu, get_choice
from metainformant.menu.navigation import MenuSystem

# Discover scripts
repo_root = Path(".")
scripts = discover_scripts(repo_root)

# Generate menu structure
menus = generate_menu_from_scripts(scripts)

# Create menu system
menu_system = MenuSystem(menus=menus, current_menu_id="root")

# Display menu
current_menu = menu_system.get_current_menu()
show_menu(current_menu.items, current_menu.title)

# Get user choice
choice = get_choice("Select option: ")
```

### Direct Script Execution

Execute scripts directly via orchestrator:

```bash
# Execute specific script
./metainformant.sh rna run_workflow.py --config config.yaml

# Show category menu
./metainformant.sh rna

# List all scripts
./metainformant.sh --list
```

## Menu Structure

Menus are organized hierarchically:

1. **Root Menu**: Lists all script categories
2. **Category Menus**: Lists scripts within a category
3. **Script Execution**: Executes selected script

### Menu Item Actions

Menu items can have three types of actions:

1. **Submenu Navigation**: `"submenu:<menu_id>"` - Navigate to submenu
2. **Script Execution**: `"script:<path>"` - Execute script file
3. **Callable**: Python callable function - Execute function

## Script Discovery

The discovery system automatically:

- Scans `scripts/` directory recursively
- Extracts metadata from Python docstrings
- Extracts metadata from bash script comments
- Categorizes scripts by directory structure
- Skips test files and cache directories

### Metadata Extraction

**Python Scripts:**
- Extracts module docstring as description
- Parses argparse arguments (if present)
- Identifies required vs optional arguments

**Bash Scripts:**
- Extracts first meaningful comment as description
- Identifies script type from extension

## Error Handling

The menu system handles errors gracefully:

- Invalid menu selections show error messages
- Script execution errors are reported with exit codes
- Navigation errors return to previous menu
- Keyboard interrupts (Ctrl+C) exit cleanly

## Examples

### Creating Custom Menu

```python
from metainformant.menu.navigation import Menu, MenuItem, MenuSystem

# Create custom menu
items = [
    MenuItem(
        id="option1",
        label="Run Analysis",
        description="Execute analysis workflow",
        action="script:scripts/analysis/run.py"
    ),
    MenuItem(
        id="option2",
        label="View Results",
        description="Display analysis results",
        action=lambda: print("Results")
    )
]

menu = Menu(id="custom", title="Custom Menu", items=items)
menus = {"custom": menu}
system = MenuSystem(menus=menus, current_menu_id="custom")
```

### Discovering Scripts

```python
from pathlib import Path
from metainformant.menu import discover_scripts

repo_root = Path(".")
scripts = discover_scripts(repo_root)

for category, script_list in scripts.items():
    print(f"{category}:")
    for script in script_list:
        print(f"  {script.name}: {script.description}")
```

## Integration

The menu module integrates with:

- **Orchestrator Script** (`metainformant.sh`): Provides CLI interface
- **Script Discovery**: Automatically finds scripts in repository
- **Core Utilities**: Uses `metainformant.core` for path handling

## Testing

Comprehensive test suite covers:

- Menu display and formatting
- Navigation state management
- Script discovery and metadata extraction
- Script execution and error handling

Run tests with:

```bash
pytest tests/test_menu_*.py
```

## See Also

- [`metainformant.sh`](../../../metainformant.sh): Top-level orchestrator script
- [`scripts/README.md`](../../../scripts/README.md): Scripts directory documentation
- [`docs/ORCHESTRATION.md`](../../../docs/ORCHESTRATION.md): Orchestration guide




