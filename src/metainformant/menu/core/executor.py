"""Script execution utilities.

This module provides functions for executing scripts with proper environment
setup, argument handling, and error reporting.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .discovery import ScriptInfo


def validate_script_executable(script_path: Path) -> bool:
    """Validate that a script exists and is executable.

    Args:
        script_path: Path to script file

    Returns:
        True if script is valid and executable, False otherwise
    """
    if not script_path.exists():
        return False

    if script_path.suffix == ".py":
        return True  # Python scripts are executed via interpreter

    if script_path.suffix == ".sh":
        # Check if bash is available
        return os.access(script_path, os.X_OK) or True  # Can always chmod

    return False


def prompt_for_args(script_info: ScriptInfo) -> list[str]:
    """Prompt user for script arguments.

    Args:
        script_info: ScriptInfo with argument metadata

    Returns:
        List of argument strings
    """
    args: list[str] = []

    # Prompt for required arguments
    for arg_name in script_info.required_args:
        while True:
            value = input(f"Enter {arg_name}: ").strip()
            if value:
                args.extend([f"--{arg_name}", value])
                break
            print(f"{arg_name} is required.")

    # Prompt for optional arguments
    for arg_name in script_info.optional_args:
        value = input(f"Enter {arg_name} (optional, press Enter to skip): ").strip()
        if value:
            args.extend([f"--{arg_name}", value])

    return args


def execute_python_script(script_path: Path, args: list[str] | None = None) -> int:
    """Execute a Python script.

    Args:
        script_path: Path to Python script
        args: Optional list of command-line arguments

    Returns:
        Exit code from script execution
    """
    if not script_path.exists():
        print(f"Error: Script not found: {script_path}")
        return 1

    cmd = [sys.executable, str(script_path)]
    if args:
        cmd.extend(args)

    try:
        result = subprocess.run(cmd, check=False, cwd=script_path.parent)
        return result.returncode
    except Exception as e:
        print(f"Error executing script: {e}")
        return 1


def execute_bash_script(script_path: Path, args: list[str] | None = None) -> int:
    """Execute a bash script.

    Args:
        script_path: Path to bash script
        args: Optional list of command-line arguments

    Returns:
        Exit code from script execution
    """
    if not script_path.exists():
        print(f"Error: Script not found: {script_path}")
        return 1

    # Ensure script is executable
    script_path.chmod(0o755)

    cmd = ["bash", str(script_path)]
    if args:
        cmd.extend(args)

    try:
        result = subprocess.run(cmd, check=False, cwd=script_path.parent)
        return result.returncode
    except Exception as e:
        print(f"Error executing script: {e}")
        return 1


def execute_script(script_path: Path, args: list[str] | None = None) -> int:
    """Execute a script (Python or bash).

    Args:
        script_path: Path to script file
        args: Optional list of command-line arguments

    Returns:
        Exit code from script execution
    """
    if script_path.suffix == ".py":
        return execute_python_script(script_path, args)
    elif script_path.suffix == ".sh":
        return execute_bash_script(script_path, args)
    else:
        print(f"Error: Unsupported script type: {script_path.suffix}")
        return 1




