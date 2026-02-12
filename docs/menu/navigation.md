# Interactive Menu Navigation

The navigation module manages the state of the interactive menu system. It
provides data classes for menu items, menus, and navigation history, plus a
`MenuSystem` class that tracks the current position and supports forward/back
traversal with breadcrumb display.

## Concepts

### Menu Hierarchy

A menu is a container holding an ordered list of `MenuItem` entries. Each item
carries an `action` string that describes what happens when selected. Actions
follow a simple protocol: `"submenu:<id>"` navigates to another menu,
`"script:<path>"` executes a script, or a Python callable is invoked directly.
Menus link to their parent via `parent_id`, forming a tree rooted at the menu
whose `id` is `"root"`.

### Navigation History

`MenuHistory` records every menu the user enters as a stack of `(menu_id, label)`
pairs. The breadcrumb path (e.g., "Home > RNA > Workflows") is derived from the
label stack and rendered by the display module. Calling `go_back` pops the current
menu and returns to the previous one.

### Display and Input

The companion `display` module provides `format_menu` for rendering items as a
numbered list, `show_menu` for printing to stdout, `get_choice` for validated
user input, `format_breadcrumb` for path rendering, and `clear_screen` using
ANSI escape codes.

## Function Reference

### `navigate_to_submenu(menu_system: MenuSystem, submenu_id: str) -> None`

Navigate the menu system to the specified submenu. Raises `KeyError` if the
submenu ID does not exist in the registered menus.

### `go_back(menu_system: MenuSystem) -> bool`

Move back one level in the navigation history. Returns `True` on success or
`False` when already at the root.

### `get_current_menu(menu_system: MenuSystem) -> Menu`

Return the `Menu` object for the currently active position.

### `format_menu(items: list[MenuItem], title: str, width: int) -> str`

Render menu items into a formatted string with numbered entries, descriptions
truncated to fit `width`, and a "0. Back/Exit" footer.

### `show_menu(items: list[MenuItem], title: str) -> None`

Print the formatted menu to stdout.

### `get_choice(prompt: str, valid_choices: list[str] | None) -> str`

Prompt the user for input. When `valid_choices` is provided, loops until the
user enters a value from the list. Raises `KeyboardInterrupt` or `EOFError` on
Ctrl+C or end-of-input.

### `format_breadcrumb(path: list[str]) -> str`

Join a list of menu labels with " > " separators. Returns `"Home"` for an empty
path.

### `clear_screen() -> None`

Clear the terminal using ANSI escape codes when stdout is a TTY, or print
newlines as a fallback.

## Data Classes

### `MenuItem`

| Field         | Type                  | Description                              |
|---------------|-----------------------|------------------------------------------|
| `id`          | `str`                 | Unique identifier                        |
| `label`       | `str`                 | Display text                             |
| `description` | `str`                 | Optional longer description              |
| `action`      | `str \| Callable`     | Navigation target or callable            |
| `enabled`     | `bool`                | Whether the item is selectable           |

### `Menu`

| Field       | Type              | Description                                |
|-------------|-------------------|--------------------------------------------|
| `id`        | `str`             | Unique menu identifier                     |
| `title`     | `str`             | Display title                              |
| `items`     | `list[MenuItem]`  | Ordered list of items                      |
| `parent_id` | `str \| None`     | ID of parent menu, `None` for root         |

### `MenuSystem`

| Field            | Type                | Description                           |
|------------------|---------------------|---------------------------------------|
| `menus`          | `dict[str, Menu]`   | All registered menus by ID            |
| `current_menu_id`| `str`               | ID of the active menu                 |
| `history`        | `MenuHistory`       | Navigation breadcrumb stack           |

## Code Examples

```python
from metainformant.menu.ui.navigation import (
    Menu, MenuItem, MenuSystem, navigate_to_submenu, go_back, get_current_menu,
)
from metainformant.menu.ui.display import format_menu, format_breadcrumb, show_menu

# Build a two-level menu
root_items = [
    MenuItem(id="rna", label="RNA Analysis", description="RNA-seq tools", action="submenu:menu_rna"),
]
rna_items = [
    MenuItem(id="quant", label="Quantify", description="Run quantification", action="script:quant.py"),
]
menus = {
    "root": Menu(id="root", title="Main Menu", items=root_items),
    "menu_rna": Menu(id="menu_rna", title="RNA Scripts", items=rna_items, parent_id="root"),
}
system = MenuSystem(menus=menus, current_menu_id="root")

# Display, navigate forward, then back
show_menu(get_current_menu(system).items, title="Main Menu")
navigate_to_submenu(system, "menu_rna")
print(format_breadcrumb(system.history.get_path()))  # "Main Menu > RNA Scripts"
go_back(system)
```

## Configuration

No external config files. Behavior is controlled by constructor arguments:
`menus` (registered menu set), `current_menu_id` (start position), `width` on
`format_menu` (default 80), and `enabled` on `MenuItem` (set `False` to hide).

## Import Paths

```python
from metainformant.menu.ui.navigation import MenuSystem, Menu, MenuItem, MenuHistory
from metainformant.menu.ui.display import format_menu, show_menu, get_choice
```
