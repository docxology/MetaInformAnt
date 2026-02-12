# Menu UI

Terminal-based menu display, navigation, and user interaction for the interactive script launcher.

## Contents

| File | Purpose |
|------|---------|
| `display.py` | Menu formatting, rendering, user input, and screen management |
| `navigation.py` | Menu hierarchy, history tracking, and navigation state machine |

## Key Classes and Functions

| Symbol | Description |
|--------|-------------|
| `MenuItem` | Dataclass for a single menu entry with label and action |
| `Menu` | Collection of MenuItems with title and metadata |
| `MenuSystem` | Stateful menu controller with history and navigation |
| `MenuHistory` | Stack-based navigation history for back/forward |
| `format_menu()` | Render menu items as formatted terminal text |
| `show_menu()` | Display menu and wait for user selection |
| `get_choice()` | Prompt user for validated input |
| `navigate_to_submenu()` | Push current menu and enter submenu |
| `go_back()` | Pop navigation history and return to parent |
| `format_breadcrumb()` | Render navigation path as breadcrumb string |

## Usage

```python
from metainformant.menu.ui.navigation import MenuSystem, navigate_to_submenu
from metainformant.menu.ui.display import show_menu, format_breadcrumb

menu_system = MenuSystem(root_menu)
navigate_to_submenu(menu_system, "rna_workflows")
breadcrumb = format_breadcrumb(["Home", "RNA", "Workflows"])
```
