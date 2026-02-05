# Menu Module Rules

## Purpose
Interactive command-line menu system for discovering and executing METAINFORMANT analysis workflows via script discovery and categorization.

## Source Structure
```
src/metainformant/menu/
├── core/
│   ├── discovery.py   # Script discovery, metadata extraction, categorization
│   └── executor.py    # Script execution (Python/Bash), argument prompting, validation
└── ui/
    ├── display.py     # Menu formatting, screen clearing, breadcrumb navigation
    └── navigation.py  # Menu data structures (MenuItem, Menu, MenuSystem, MenuHistory)
```

## Dependencies
- **Required**: Standard library only (pathlib, subprocess, shutil)

## Import Patterns
```python
from metainformant.menu.core import discovery, executor
from metainformant.menu.ui import display, navigation
```

## Configuration
- No environment prefix needed (internal tool)
- Scans `scripts/` directory for available workflows

## Integration
- **Menu → All Modules**: Discovers and launches scripts from any module
- **Menu → Core**: Uses core utilities for script execution

## Testing
- Test script discovery against actual `scripts/` directory
- Test menu navigation logic with real MenuItem structures
- All test outputs to `tmp_path`
