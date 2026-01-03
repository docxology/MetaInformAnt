"""Dependency checking utilities for RNA-seq workflow.

This module provides utilities for checking external tool availability
and system requirements for RNA-seq analysis.
"""

from __future__ import annotations

import shutil
import subprocess
from typing import Dict, List, Tuple

from metainformant.core import logging

logger = logging.get_logger(__name__)


def check_amalgkit_availability() -> Tuple[bool, str]:
    """Check if amalgkit is available on the system.

    Returns:
        Tuple of (available, version_or_error_message)
    """
    try:
        result = subprocess.run(['amalgkit', '--version'],
                              capture_output=True, text=True, timeout=10)

        if result.returncode == 0:
            version = result.stdout.strip() or result.stderr.strip()
            return True, version
        else:
            return False, f"Command failed with return code {result.returncode}"

    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError) as e:
        return False, str(e)


def check_quantification_tools() -> Dict[str, Tuple[bool, str]]:
    """Check availability of quantification tools.

    Returns:
        Dictionary mapping tool names to (available, version) tuples
    """
    tools = {}

    # Check kallisto
    try:
        result = subprocess.run(['kallisto', 'version'],
                              capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            version = result.stdout.strip() or result.stderr.strip()
            tools['kallisto'] = (True, version)
        else:
            tools['kallisto'] = (False, f"Return code {result.returncode}")
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError) as e:
        tools['kallisto'] = (False, str(e))

    # Check salmon
    try:
        result = subprocess.run(['salmon', '--version'],
                              capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            version = result.stdout.strip() or result.stderr.strip()
            tools['salmon'] = (True, version)
        else:
            tools['salmon'] = (False, f"Return code {result.returncode}")
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError) as e:
        tools['salmon'] = (False, str(e))

    return tools


def check_system_requirements() -> Dict[str, any]:
    """Check system requirements for RNA-seq analysis.

    Returns:
        Dictionary with system requirement check results
    """
    requirements = {
        'memory_gb': check_memory_gb(),
        'cpu_cores': check_cpu_cores(),
        'disk_space_gb': check_disk_space_gb(),
        'python_version': check_python_version(),
        'required_packages': check_required_packages()
    }

    return requirements


def check_memory_gb() -> float:
    """Check available system memory in GB."""
    try:
        import psutil
        memory = psutil.virtual_memory()
        return memory.available / (1024 ** 3)  # Convert to GB
    except ImportError:
        # Fallback to reading /proc/meminfo on Linux
        try:
            with open('/proc/meminfo', 'r') as f:
                for line in f:
                    if line.startswith('MemAvailable'):
                        kb = int(line.split()[1])
                        return kb / (1024 ** 2)  # Convert to GB
        except (FileNotFoundError, ValueError):
            pass

    return 0.0  # Unknown


def check_cpu_cores() -> int:
    """Check number of available CPU cores."""
    try:
        import multiprocessing
        return multiprocessing.cpu_count()
    except ImportError:
        return 1  # Conservative estimate


def check_disk_space_gb() -> float:
    """Check available disk space in GB."""
    try:
        import psutil
        disk = psutil.disk_usage('/')
        return disk.free / (1024 ** 3)  # Convert to GB
    except ImportError:
        return 0.0  # Unknown


def check_python_version() -> str:
    """Check Python version."""
    import sys
    return f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"


def check_required_packages() -> Dict[str, Tuple[bool, str]]:
    """Check if required Python packages are available.

    Returns:
        Dictionary mapping package names to (available, version) tuples
    """
    required_packages = [
        'numpy', 'pandas', 'scipy', 'matplotlib', 'seaborn',
        'requests', 'biopython', 'metainformant'
    ]

    packages = {}

    for package in required_packages:
        try:
            module = __import__(package)
            version = getattr(module, '__version__', 'unknown')
            packages[package] = (True, version)
        except ImportError:
            packages[package] = (False, 'not installed')

    return packages


def get_dependency_report() -> Dict[str, any]:
    """Generate a comprehensive dependency report.

    Returns:
        Dictionary with complete dependency analysis
    """
    report = {
        'amalgkit': check_amalgkit_availability(),
        'quantification_tools': check_quantification_tools(),
        'system_requirements': check_system_requirements(),
        'timestamp': __import__('time').time()
    }

    # Overall assessment
    amalgkit_ok = report['amalgkit'][0]
    quant_tools = [tool for tool, (avail, _) in report['quantification_tools'].items() if avail]
    system_ok = (report['system_requirements']['memory_gb'] >= 8 and
                 report['system_requirements']['cpu_cores'] >= 4)

    report['overall_status'] = 'READY' if (amalgkit_ok and quant_tools and system_ok) else 'ISSUES'

    return report


def print_dependency_report() -> None:
    """Print a formatted dependency report to console."""
    report = get_dependency_report()

    print("RNA-seq Workflow Dependency Report")
    print("=" * 40)

    # Amalgkit status
    amalgkit_ok, amalgkit_info = report['amalgkit']
    status = "✓ Available" if amalgkit_ok else "✗ Missing"
    print(f"Amalgkit: {status} ({amalgkit_info})")

    # Quantification tools
    print("\nQuantification Tools:")
    for tool, (avail, version) in report['quantification_tools'].items():
        status = "✓" if avail else "✗"
        print(f"  {tool}: {status} {version}")

    # System requirements
    sys_req = report['system_requirements']
    print(f"\nSystem Requirements:")
    print(f"  Memory: {sys_req['memory_gb']:.1f} GB")
    print(f"  CPU Cores: {sys_req['cpu_cores']}")
    print(f"  Disk Space: {sys_req['disk_space_gb']:.1f} GB")
    print(f"  Python: {sys_req['python_version']}")

    print(f"\nOverall Status: {report['overall_status']}")


def ensure_dependencies() -> bool:
    """Ensure all critical dependencies are available.

    Returns:
        True if all dependencies are satisfied
    """
    report = get_dependency_report()

    if report['overall_status'] == 'ISSUES':
        print("Dependency issues detected:")
        print_dependency_report()
        return False

    return True







