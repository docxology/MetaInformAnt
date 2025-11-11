"""Environment and dependency checking functions.

This module provides functions for checking tool availability, dependencies,
and environment validation for RNA-seq workflows.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from ..core.logging import get_logger

logger = get_logger(__name__)


def check_amalgkit() -> tuple[bool, str]:
    """Check if amalgkit is available and get version.
    
    Returns:
        Tuple of (is_available: bool, message: str)
    """
    amalgkit_path = shutil.which("amalgkit")
    if not amalgkit_path:
        return False, "Not found on PATH"
    
    try:
        result = subprocess.run(
            ["amalgkit", "-h"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0:
            output = result.stdout or result.stderr
            # Extract version or help text
            for line in output.split("\n")[:5]:
                if "version" in line.lower():
                    return True, line.strip()
            return True, f"Found at {amalgkit_path}"
        else:
            return False, f"Error: {result.stderr}"
    except Exception as e:
        return False, f"Error checking version: {e}"


def check_sra_toolkit() -> tuple[bool, str]:
    """Check if SRA Toolkit is installed.
    
    Returns:
        Tuple of (is_available: bool, message: str)
    """
    fasterq_path = shutil.which("fasterq-dump")
    if not fasterq_path:
        return False, "fasterq-dump not found"
    
    try:
        result = subprocess.run(
            ["fasterq-dump", "--version"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0:
            for line in result.stdout.split("\n"):
                if "fasterq-dump" in line.lower():
                    return True, line.strip()
            return True, f"Found at {fasterq_path}"
        else:
            return False, f"Error: {result.stderr}"
    except Exception as e:
        return False, f"Error checking version: {e}"


def check_kallisto() -> tuple[bool, str]:
    """Check if kallisto is installed.
    
    Returns:
        Tuple of (is_available: bool, message: str)
    """
    kallisto_path = shutil.which("kallisto")
    if not kallisto_path:
        return False, "kallisto not found"
    
    try:
        result = subprocess.run(
            ["kallisto", "version"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0:
            return True, result.stdout.strip()
        else:
            return False, f"Error: {result.stderr}"
    except Exception as e:
        return False, f"Error checking version: {e}"


def check_metainformant() -> tuple[bool, str]:
    """Check if metainformant package is installed.
    
    Returns:
        Tuple of (is_available: bool, message: str)
    """
    try:
        import metainformant
        version = getattr(metainformant, "__version__", "unknown")
        return True, f"v{version}"
    except ImportError as e:
        return False, f"Import error: {e}"


def check_virtual_env() -> tuple[bool, str]:
    """Check if running inside a virtual environment.
    
    Returns:
        Tuple of (is_in_venv: bool, message: str)
    """
    import os
    venv_path = os.environ.get("VIRTUAL_ENV")
    
    if venv_path:
        return True, f"Active: {venv_path}"
    
    # Check if sys.prefix differs from sys.base_prefix
    if sys.prefix != sys.base_prefix:
        return True, f"Active: {sys.prefix}"
    
    return False, "Not in virtual environment"


def check_rscript() -> tuple[bool, str]:
    """Check if Rscript is available.
    
    Returns:
        Tuple of (is_available: bool, message: str)
    """
    rscript_path = shutil.which("Rscript")
    if not rscript_path:
        return False, "Rscript not found"
    
    try:
        result = subprocess.run(
            ["Rscript", "--version"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0:
            version_info = result.stderr.strip() if result.stderr else result.stdout.strip()
            return True, version_info
        else:
            return False, f"Error: {result.stderr}"
    except Exception as e:
        return False, f"Error checking version: {e}"


def check_dependencies() -> dict[str, tuple[bool, str]]:
    """Check all required dependencies for RNA-seq workflows.
    
    Returns:
        Dictionary mapping dependency name -> (is_available: bool, message: str)
    """
    return {
        "virtual_env": check_virtual_env(),
        "metainformant": check_metainformant(),
        "amalgkit": check_amalgkit(),
        "sra_toolkit": check_sra_toolkit(),
        "kallisto": check_kallisto(),
        "rscript": check_rscript(),
    }


def validate_environment() -> dict[str, Any]:
    """Comprehensive environment validation.
    
    Returns:
        Dictionary with validation results:
        - all_passed: bool
        - dependencies: dict mapping name -> (is_available, message)
        - recommendations: list of strings with recommendations
    """
    deps = check_dependencies()
    all_passed = all(available for available, _ in deps.values())
    
    recommendations = []
    if not deps["virtual_env"][0]:
        recommendations.append("Activate virtual environment: source .venv/bin/activate")
    if not deps["metainformant"][0]:
        recommendations.append("Install metainformant: uv pip install -e .")
    if not deps["amalgkit"][0]:
        recommendations.append("Install amalgkit: uv pip install git+https://github.com/kfuku52/amalgkit")
    if not deps["sra_toolkit"][0]:
        recommendations.append("Install SRA Toolkit: sudo apt-get install sra-toolkit")
    if not deps["kallisto"][0]:
        recommendations.append("Install kallisto: sudo apt-get install kallisto")
    if not deps["rscript"][0]:
        recommendations.append("Install R: sudo apt-get install r-base")
    
    return {
        "all_passed": all_passed,
        "dependencies": deps,
        "recommendations": recommendations,
    }


