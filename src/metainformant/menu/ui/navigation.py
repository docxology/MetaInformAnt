"""Interactive menu navigation system.

This module provides classes and functions for managing menu navigation state,
including submenu navigation, breadcrumbs, and history tracking.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable


@dataclass
class MenuItem:
    """Single menu item with action and metadata."""

    id: str
    label: str
    description: str = ""
    action: str | Callable = ""  # "submenu:<id>", "script:<path>", or callable
    enabled: bool = True


@dataclass
class Menu:
    """Menu container with items and metadata."""

    id: str
    title: str
    items: list[MenuItem]
    parent_id: str | None = None


@dataclass
class MenuHistory:
    """Navigation history for breadcrumb tracking."""

    path: list[str] = field(default_factory=list)
    menu_ids: list[str] = field(default_factory=list)

    def push(self, menu_id: str, label: str) -> None:
        """Add a menu to the history.

        Args:
            menu_id: ID of the menu
            label: Display label for breadcrumb
        """
        self.menu_ids.append(menu_id)
        self.path.append(label)

    def pop(self) -> tuple[str, str] | None:
        """Remove and return the last menu from history.

        Returns:
            Tuple of (menu_id, label) or None if history is empty
        """
        if not self.menu_ids:
            return None
        menu_id = self.menu_ids.pop()
        label = self.path.pop()
        return (menu_id, label)

    def get_path(self) -> list[str]:
        """Get current breadcrumb path.

        Returns:
            List of menu labels
        """
        return self.path.copy()

    def clear(self) -> None:
        """Clear navigation history."""
        self.path.clear()
        self.menu_ids.clear()


@dataclass
class MenuSystem:
    """Menu navigation system with state management."""

    menus: dict[str, Menu]
    current_menu_id: str
    history: MenuHistory = field(default_factory=MenuHistory)

    def __post_init__(self) -> None:
        """Initialize history with root menu."""
        if self.current_menu_id in self.menus:
            root_menu = self.menus[self.current_menu_id]
            self.history.push(self.current_menu_id, root_menu.title)

    def get_current_menu(self) -> Menu:
        """Get the currently active menu.

        Returns:
            Current menu object

        Raises:
            KeyError: If current menu ID is invalid
        """
        return self.menus[self.current_menu_id]

    def navigate_to(self, menu_id: str) -> bool:
        """Navigate to a specific menu.

        Args:
            menu_id: ID of menu to navigate to

        Returns:
            True if navigation successful, False otherwise
        """
        if menu_id not in self.menus:
            return False

        menu = self.menus[menu_id]
        self.current_menu_id = menu_id
        self.history.push(menu_id, menu.title)
        return True

    def go_back(self) -> bool:
        """Navigate back to previous menu.

        Returns:
            True if navigation successful, False if at root
        """
        result = self.history.pop()
        if result is None:
            return False

        # Pop current menu, get previous
        if self.history.menu_ids:
            prev_menu_id = self.history.menu_ids[-1]
            self.current_menu_id = prev_menu_id
            return True

        return False


def navigate_to_submenu(menu_system: MenuSystem, submenu_id: str) -> None:
    """Navigate to a submenu.

    Args:
        menu_system: Menu system instance
        submenu_id: ID of submenu to navigate to

    Raises:
        KeyError: If submenu doesn't exist
    """
    if not menu_system.navigate_to(submenu_id):
        raise KeyError(f"Submenu '{submenu_id}' not found")


def go_back(menu_system: MenuSystem) -> bool:
    """Navigate back in menu history.

    Args:
        menu_system: Menu system instance

    Returns:
        True if navigation successful, False if at root
    """
    return menu_system.go_back()


def get_current_menu(menu_system: MenuSystem) -> Menu:
    """Get the currently active menu.

    Args:
        menu_system: Menu system instance

    Returns:
        Current menu object
    """
    return menu_system.get_current_menu()
