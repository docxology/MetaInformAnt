"""Menu display and formatting utilities.

This module provides functions for rendering menus, handling user input,
and formatting menu output with consistent styling.
"""

from __future__ import annotations

import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .navigation import MenuItem


def format_menu(items: list[MenuItem], title: str = "", width: int = 80) -> str:
    """Format a list of menu items into a displayable menu string.

    Args:
        items: List of menu items to display
        title: Optional title for the menu
        width: Maximum width for menu display

    Returns:
        Formatted menu string ready for display
    """
    if not items:
        return "No menu items available.\n"

    lines: list[str] = []
    if title:
        title_line = f" {title} ".center(width, "=")
        lines.append(title_line)
        lines.append("")

    for i, item in enumerate(items, 1):
        if not item.enabled:
            continue

        label = item.label
        if item.description:
            desc_width = width - len(str(i)) - len(label) - 6
            if desc_width > 20:
                description = item.description[:desc_width]
                if len(item.description) > desc_width:
                    description += "..."
                lines.append(f"  {i}. {label:<30} {description}")
            else:
                lines.append(f"  {i}. {label}")
        else:
            lines.append(f"  {i}. {label}")

    lines.append("")
    lines.append("  0. Back/Exit")
    lines.append("")

    return "\n".join(lines)


def show_menu(items: list[MenuItem], title: str = "") -> None:
    """Display a menu to the user.

    Args:
        items: List of menu items to display
        title: Optional title for the menu
    """
    menu_text = format_menu(items, title)
    print(menu_text)


def get_choice(prompt: str, valid_choices: list[str] | None = None) -> str:
    """Get user choice with validation.

    Args:
        prompt: Prompt text to display
        valid_choices: Optional list of valid choice strings

    Returns:
        User's choice string

    Raises:
        KeyboardInterrupt: If user interrupts (Ctrl+C)
        EOFError: If input stream ends
    """
    while True:
        try:
            choice = input(f"{prompt} ").strip()
            if not choice:
                continue

            if valid_choices is None:
                return choice

            if choice in valid_choices:
                return choice

            print(f"Invalid choice. Valid options: {', '.join(valid_choices)}")
        except (KeyboardInterrupt, EOFError):
            print("\n")
            raise


def format_breadcrumb(path: list[str]) -> str:
    """Format a breadcrumb path for display.

    Args:
        path: List of menu IDs representing the navigation path

    Returns:
        Formatted breadcrumb string (e.g., "Home > RNA > Workflows")
    """
    if not path:
        return "Home"

    return " > ".join(path)


def clear_screen() -> None:
    """Clear the terminal screen.

    Uses ANSI escape codes for cross-platform compatibility.
    Falls back to printing newlines if ANSI codes aren't supported.
    """
    if sys.stdout.isatty():
        print("\033[2J\033[H", end="")
    else:
        print("\n" * 50)




