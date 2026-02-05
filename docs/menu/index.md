# Menu Module

Interactive command-line menu system for discovering and executing METAINFORMANT analysis workflows.

## Submodules

### Core (`core/`)
- **discovery.py** - Script discovery, metadata extraction, categorization, menu generation
- **executor.py** - Script execution (Python/Bash), argument prompting, validation

### UI (`ui/`)
- **display.py** - Menu formatting, screen clearing, breadcrumb navigation, user choice handling
- **navigation.py** - Menu data structures (MenuItem, Menu, MenuSystem, MenuHistory)

## Usage

```python
from metainformant.menu.core import discovery, executor
from metainformant.menu.ui import display, navigation
```

The menu system scans the `scripts/` directory for available analysis workflows and presents them in an interactive terminal interface with categorized navigation.
