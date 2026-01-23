"""RNA environment validation and dependency checking.

This module provides comprehensive tools for validating the runtime environment
and checking for required/optional external dependencies and CLI tools needed
for RNA-seq analysis workflows.

Main Functions:
    - check_amalgkit: Check if amalgkit CLI is available
    - check_sra_toolkit: Check if SRA Toolkit is available
    - check_kallisto: Check if Kallisto is available
    - check_metainformant: Check if MetaInformant is installed
    - check_virtual_env: Check if running in a virtual environment
    - check_rscript: Check if R/Rscript is available
    - check_dependencies: Check all dependencies at once
    - validate_environment: Comprehensive environment validation

Example:
    >>> from metainformant.rna import environment
    >>> ok, msg = environment.check_amalgkit()
    >>> if ok:
    ...     print(f"Amalgkit available: {msg}")
    ...
    >>> deps = environment.check_dependencies()
    >>> validation = environment.validate_environment()
"""

from __future__ import annotations

import os
import subprocess
import sys
from typing import Any, Dict, Tuple

from metainformant.core import logging

logger = logging.get_logger(__name__)


def check_amalgkit() -> Tuple[bool, str]:
    """Check if amalgkit CLI tool is available.

    Attempts to run 'amalgkit --version' to verify installation and availability.

    Returns:
        Tuple of (available, message) where available is True if command succeeds,
        and message contains version info on success or error description on failure.
    """
    try:
        result = subprocess.run(
            ["amalgkit", "--version"],
            capture_output=True,
            text=True,
            timeout=10,
        )

        if result.returncode == 0:
            version = result.stdout.strip() or result.stderr.strip()
            return True, version
        else:
            return False, f"Command failed with return code {result.returncode}"

    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError) as e:
        return False, str(e)


def check_sra_toolkit() -> Tuple[bool, str]:
    """Check if SRA Toolkit is available.

    Attempts to run 'fastq-dump --version' to verify SRA Toolkit installation.
    Falls back to checking 'prefetch' if fastq-dump not available.

    Returns:
        Tuple of (available, message) where available is True if SRA tools found,
        and message contains version info on success or error description on failure.
    """
    # Try fastq-dump first (standard SRA tool)
    for cmd in ["fastq-dump", "prefetch"]:
        try:
            result = subprocess.run(
                [cmd, "--version"],
                capture_output=True,
                text=True,
                timeout=10,
            )

            if result.returncode == 0:
                version = result.stdout.strip() or result.stderr.strip()
                return True, f"{cmd}: {version}"

        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError):
            continue

    return False, "SRA Toolkit not found (fastq-dump/prefetch not available)"


def check_kallisto() -> Tuple[bool, str]:
    """Check if Kallisto is available.

    Attempts to run 'kallisto version' to verify Kallisto installation.

    Returns:
        Tuple of (available, message) where available is True if Kallisto found,
        and message contains version info on success or error description on failure.
    """
    try:
        result = subprocess.run(
            ["kallisto", "version"],
            capture_output=True,
            text=True,
            timeout=10,
        )

        if result.returncode == 0:
            version = result.stdout.strip() or result.stderr.strip()
            return True, version
        else:
            return False, f"Command failed with return code {result.returncode}"

    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError) as e:
        return False, str(e)


def check_metainformant() -> Tuple[bool, str]:
    """Check if MetaInformant is installed and importable.

    Attempts to import metainformant and get version information.
    This check should always succeed in normal test/runtime environments.

    Returns:
        Tuple of (True, version_info) since MetaInformant should be available
        when this code is running.
    """
    try:
        import metainformant

        version = getattr(metainformant, "__version__", "unknown")
        return True, f"v{version}"

    except ImportError:
        return False, "MetaInformant not importable"


def check_virtual_env() -> Tuple[bool, str]:
    """Check if running in a Python virtual environment.

    Uses multiple heuristics to detect virtual environment:
    - VIRTUAL_ENV environment variable
    - sys.prefix != sys.base_prefix (PEP 405)
    - .venv or venv directory detection

    Returns:
        Tuple of (in_venv, message) where in_venv is True if running in venv.
    """
    in_venv = False
    venv_path = None

    # Check environment variable
    if os.environ.get("VIRTUAL_ENV"):
        in_venv = True
        venv_path = os.environ.get("VIRTUAL_ENV")

    # Check sys.prefix vs base_prefix (PEP 405)
    elif sys.prefix != sys.base_prefix:
        in_venv = True
        venv_path = sys.prefix

    if in_venv and venv_path:
        return True, f"Virtual environment at {venv_path}"
    elif in_venv:
        return True, "Virtual environment detected"
    else:
        return False, "Not running in virtual environment"


def check_rscript() -> Tuple[bool, str]:
    """Check if R/Rscript is available.

    Attempts to run 'Rscript --version' to verify R installation.

    Returns:
        Tuple of (available, message) where available is True if Rscript found,
        and message contains version info on success or error description on failure.
    """
    try:
        result = subprocess.run(
            ["Rscript", "--version"],
            capture_output=True,
            text=True,
            timeout=10,
        )

        # Rscript --version outputs to stderr
        output = result.stderr.strip() or result.stdout.strip()

        if result.returncode == 0 and output:
            return True, output
        elif result.returncode == 0:
            return True, "Rscript available"
        else:
            return False, f"Command failed with return code {result.returncode}"

    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError) as e:
        return False, str(e)


def check_dependencies() -> Dict[str, Tuple[bool, str]]:
    """Check all key dependencies at once.

    Runs all individual dependency checks and returns results in a dictionary.

    Returns:
        Dictionary with keys (str) for each dependency and values (bool, str) tuples:
        - 'virtual_env': Virtual environment status
        - 'metainformant': MetaInformant installation status
        - 'amalgkit': Amalgkit CLI availability
        - 'sra_toolkit': SRA Toolkit availability
        - 'kallisto': Kallisto availability
        - 'rscript': R/Rscript availability

    Example:
        >>> deps = check_dependencies()
        >>> if not deps['amalgkit'][0]:
        ...     print(f"Install amalgkit: {deps['amalgkit'][1]}")
    """
    return {
        "virtual_env": check_virtual_env(),
        "metainformant": check_metainformant(),
        "amalgkit": check_amalgkit(),
        "sra_toolkit": check_sra_toolkit(),
        "kallisto": check_kallisto(),
        "rscript": check_rscript(),
    }


def validate_environment() -> Dict[str, Any]:
    """Comprehensive environment validation with recommendations.

    Checks all dependencies and provides actionable recommendations for any
    missing tools or configurations.

    Returns:
        Dictionary with keys:
        - 'all_passed' (bool): True if all checks passed
        - 'dependencies' (dict): Results from check_dependencies()
        - 'recommendations' (list[str]): Installation/configuration recommendations

    Example:
        >>> validation = validate_environment()
        >>> if not validation['all_passed']:
        ...     print("Environment issues detected:")
        ...     for rec in validation['recommendations']:
        ...         print(f"  - {rec}")
    """
    deps = check_dependencies()

    # Determine overall status
    all_passed = all(available for available, _ in deps.values())

    # Build recommendations for missing dependencies
    recommendations = []

    if not deps["virtual_env"][0]:
        recommendations.append("Virtual environment recommended: python -m venv .venv && source .venv/bin/activate")

    if not deps["amalgkit"][0]:
        recommendations.append("Install amalgkit: uv pip install amalgkit")

    if not deps["sra_toolkit"][0]:
        recommendations.append("Install SRA Toolkit: apt-get install sra-toolkit or conda install sra-tools")

    if not deps["kallisto"][0]:
        recommendations.append("Install Kallisto: apt-get install kallisto or conda install kallisto")

    if not deps["rscript"][0]:
        recommendations.append("Install R: apt-get install r-base or follow https://cloud.r-project.org")

    logger.debug(f"Environment validation: all_passed={all_passed}, issues={len(recommendations)}")

    return {
        "all_passed": all_passed,
        "dependencies": deps,
        "recommendations": recommendations,
    }
