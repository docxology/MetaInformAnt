"""Interactive menu and CLI interface for METAINFORMANT.

This module provides the interactive command-line interface and menu system,
allowing users to discover and execute analysis workflows.
"""

from __future__ import annotations

from .ui.display import clear_screen, format_breadcrumb, format_menu, get_choice, show_menu
from .core.executor import execute_bash_script, execute_python_script, execute_script, prompt_for_args, validate_script_executable
from .ui.navigation import Menu, MenuHistory, MenuItem, MenuSystem, get_current_menu, go_back, navigate_to_submenu
from .core.discovery import discover_scripts, extract_script_metadata, categorize_script, generate_menu_from_scripts, ScriptInfo

__all__ = [
    # Data structures
    "MenuItem",
    "Menu",
    "MenuSystem",
    "MenuHistory",
    "ScriptInfo",
    # Display functions
    "show_menu",
    "format_menu",
    "get_choice",
    "format_breadcrumb",
    "clear_screen",
    # Navigation functions
    "navigate_to_submenu",
    "go_back",
    "get_current_menu",
    # Discovery functions
    "discover_scripts",
    "extract_script_metadata",
    "categorize_script",
    "generate_menu_from_scripts",
    # Execution functions
    "execute_script",
    "execute_python_script",
    "execute_bash_script",
    "prompt_for_args",
    "validate_script_executable",
]




