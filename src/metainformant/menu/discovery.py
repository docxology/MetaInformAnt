"""Script discovery and menu generation.

This module provides functions for discovering scripts in the scripts directory,
extracting metadata, and generating menu structures from discovered scripts.
"""

from __future__ import annotations

import ast
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .navigation import Menu, MenuItem


@dataclass
class ScriptInfo:
    """Metadata about a discovered script."""

    path: Path
    name: str
    description: str
    category: str
    script_type: str  # "python", "bash"
    required_args: list[str] = field(default_factory=list)
    optional_args: list[str] = field(default_factory=list)


def categorize_script(script_path: Path) -> str:
    """Categorize a script based on its directory path.

    Args:
        script_path: Path to the script file

    Returns:
        Category string (e.g., "rna", "gwas", "core")
    """
    parts = script_path.parts
    if "scripts" in parts:
        scripts_idx = parts.index("scripts")
        if scripts_idx + 1 < len(parts):
            category = parts[scripts_idx + 1]
            # Skip archive and test directories
            if category not in ("archive", "test_examples"):
                return category
    return "other"


def extract_script_metadata(script_path: Path) -> ScriptInfo:
    """Extract metadata from a script file.

    Args:
        script_path: Path to the script file

    Returns:
        ScriptInfo with extracted metadata
    """
    script_type = "python" if script_path.suffix == ".py" else "bash"
    name = script_path.stem
    category = categorize_script(script_path)

    description = ""
    required_args: list[str] = []
    optional_args: list[str] = []

    try:
        if script_type == "python":
            description, required_args, optional_args = _extract_python_metadata(script_path)
        else:
            description = _extract_bash_metadata(script_path)
    except Exception:
        # If extraction fails, use defaults
        pass

    return ScriptInfo(
        path=script_path,
        name=name,
        description=description or f"{script_type.title()} script",
        category=category,
        script_type=script_type,
        required_args=required_args,
        optional_args=optional_args,
    )


def _extract_python_metadata(script_path: Path) -> tuple[str, list[str], list[str]]:
    """Extract metadata from a Python script.

    Args:
        script_path: Path to Python script

    Returns:
        Tuple of (description, required_args, optional_args)
    """
    description = ""
    required_args: list[str] = []
    optional_args: list[str] = []

    try:
        content = script_path.read_text(encoding="utf-8")
        tree = ast.parse(content, filename=str(script_path))

        # Extract module docstring
        if ast.get_docstring(tree):
            description = ast.get_docstring(tree).split("\n")[0].strip()

        # Find argparse parser to extract arguments
        for node in ast.walk(tree):
            if isinstance(node, ast.Call):
                if isinstance(node.func, ast.Attribute):
                    if node.func.attr == "add_argument":
                        # Try to extract argument info
                        for keyword in node.keywords:
                            if keyword.arg in ("dest", "name"):
                                if isinstance(keyword.value, ast.Constant):
                                    arg_name = keyword.value.value
                                    # Check if required
                                    is_required = True
                                    for kw in node.keywords:
                                        if kw.arg == "required" and isinstance(kw.value, ast.Constant):
                                            is_required = kw.value.value
                                            break

                                    if is_required:
                                        required_args.append(arg_name)
                                    else:
                                        optional_args.append(arg_name)
    except Exception:
        pass

    return (description, required_args, optional_args)


def _extract_bash_metadata(script_path: Path) -> str:
    """Extract description from a bash script.

    Args:
        script_path: Path to bash script

    Returns:
        Description string
    """
    try:
        content = script_path.read_text(encoding="utf-8")
        # Look for comment-based description
        lines = content.split("\n")
        for i, line in enumerate(lines[:20]):  # Check first 20 lines
            line = line.strip()
            if line.startswith("#") and len(line) > 2:
                # Skip shebang and empty comments
                if i == 0 and line.startswith("#!"):
                    continue
                desc = line.lstrip("#").strip()
                if desc and len(desc) > 10:  # Meaningful description
                    return desc.split(".")[0]  # First sentence
    except Exception:
        pass

    return ""


def discover_scripts(repo_root: Path) -> dict[str, list[ScriptInfo]]:
    """Discover all scripts in the scripts directory.

    Args:
        repo_root: Root directory of the repository

    Returns:
        Dictionary mapping category to list of ScriptInfo objects
    """
    scripts_dir = repo_root / "scripts"
    if not scripts_dir.exists():
        return {}

    scripts_by_category: dict[str, list[ScriptInfo]] = {}

    # Discover Python scripts
    for py_script in scripts_dir.rglob("*.py"):
        # Skip __pycache__ and test files in test directories
        if "__pycache__" in py_script.parts or "test_" in py_script.name:
            continue
        # Skip template files
        if py_script.name.startswith("_") and py_script.name != "__init__.py":
            continue

        try:
            script_info = extract_script_metadata(py_script)
            category = script_info.category
            if category not in scripts_by_category:
                scripts_by_category[category] = []
            scripts_by_category[category].append(script_info)
        except Exception:
            continue

    # Discover bash scripts
    for bash_script in scripts_dir.rglob("*.sh"):
        # Skip common utilities
        if bash_script.name.startswith("_"):
            continue

        try:
            script_info = extract_script_metadata(bash_script)
            category = script_info.category
            if category not in scripts_by_category:
                scripts_by_category[category] = []
            scripts_by_category[category].append(script_info)
        except Exception:
            continue

    # Sort scripts within each category
    for category in scripts_by_category:
        scripts_by_category[category].sort(key=lambda x: x.name)

    return scripts_by_category


def generate_menu_from_scripts(scripts: dict[str, list[ScriptInfo]]) -> dict[str, Menu]:
    """Generate menu structure from discovered scripts.

    Args:
        scripts: Dictionary mapping category to list of ScriptInfo

    Returns:
        Dictionary mapping menu ID to Menu object
    """
    from .navigation import Menu, MenuItem

    menus: dict[str, Menu] = {}

    # Create root menu with categories
    root_items: list[MenuItem] = []
    for category in sorted(scripts.keys()):
        category_label = category.replace("_", " ").title()
        root_items.append(
            MenuItem(
                id=f"category_{category}",
                label=category_label,
                description=f"Scripts in {category_label} category",
                action=f"submenu:menu_{category}",
            )
        )

    menus["root"] = Menu(id="root", title="METAINFORMANT Scripts", items=root_items)

    # Create category menus
    for category, script_list in scripts.items():
        category_label = category.replace("_", " ").title()
        menu_items: list[MenuItem] = []

        for script_info in script_list:
            script_name = script_info.name.replace("_", " ").title()
            menu_items.append(
                MenuItem(
                    id=f"script_{script_info.path.stem}",
                    label=script_name,
                    description=script_info.description,
                    action=f"script:{script_info.path}",
                )
            )

        menus[f"menu_{category}"] = Menu(
            id=f"menu_{category}",
            title=f"{category_label} Scripts",
            items=menu_items,
            parent_id="root",
        )

    return menus



