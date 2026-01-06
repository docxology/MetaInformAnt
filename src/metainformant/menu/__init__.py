"""Menu system for METAINFORMANT script orchestration.

This module provides interactive text menus for accessing scripts and workflows
through a unified interface. It includes menu display, navigation, script discovery,
and execution capabilities.

Key Components:
- Menu display and formatting
- Interactive navigation with submenus and breadcrumbs
- Automatic script discovery from scripts directory
- Script execution with argument prompting
"""

from __future__ import annotations

from .discovery import ScriptInfo, categorize_script, discover_scripts, extract_script_metadata, generate_menu_from_scripts
from .display import clear_screen, format_breadcrumb, format_menu, get_choice, show_menu
from .executor import execute_bash_script, execute_python_script, execute_script, prompt_for_args, validate_script_executable
from .navigation import Menu, MenuHistory, MenuItem, MenuSystem, get_current_menu, go_back, navigate_to_submenu

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




