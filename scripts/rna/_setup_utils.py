"""
Common setup utilities for RNA scripts.

Provides automatic virtual environment setup, dependency installation, and environment validation.
"""

import os
import sys
import shutil
import subprocess
from pathlib import Path

# Add src to path before imports
repo_root = Path(__file__).parent.parent.parent.resolve()
sys.path.insert(0, str(repo_root / "src"))

from metainformant.core.logging import get_logger

logger = get_logger("rna_setup")


def setup_venv_and_dependencies(auto_setup: bool = True) -> bool:
    """
    Automatically setup virtual environment and install dependencies if needed.
    
    Args:
        auto_setup: If True, automatically create venv and install dependencies
        
    Returns:
        True if venv exists and is ready, False otherwise
    """
    venv_python = repo_root / ".venv" / "bin" / "python3"
    venv_dir = repo_root / ".venv"
    
    # Check if venv exists
    if venv_python.exists():
        logger.info(f"✓ Virtual environment found at {venv_dir}")
        return True
    
    if not auto_setup:
        logger.error("Virtual environment not found!")
        logger.error("Run: python3 -m venv .venv && source .venv/bin/activate && pip install -e . && pip install git+https://github.com/kfuku52/amalgkit")
        return False
    
    # Auto-setup venv
    logger.info("Setting up virtual environment...")
    logger.info(f"Creating venv at {venv_dir}")
    
    try:
        # Create venv (may fail on some systems, but that's OK - user can create manually)
        result = subprocess.run([sys.executable, "-m", "venv", str(venv_dir)], 
                              capture_output=True, text=True)
        if result.returncode != 0:
            logger.warning(f"Could not auto-create venv: {result.stderr}")
            logger.warning("Please create manually: python3 -m venv .venv")
            return False
        
        logger.info("✓ Virtual environment created")
        
        # Install dependencies
        logger.info("Installing dependencies...")
        
        # Install metainformant
        result = subprocess.run([str(venv_python), "-m", "pip", "install", "-e", str(repo_root)], 
                              capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"Failed to install metainformant: {result.stderr}")
            return False
        logger.info("✓ metainformant installed")
        
        # Install amalgkit
        result = subprocess.run([str(venv_python), "-m", "pip", "install", "git+https://github.com/kfuku52/amalgkit"], 
                              capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"Failed to install amalgkit: {result.stderr}")
            return False
        logger.info("✓ amalgkit installed")
        
        logger.info("✓ Setup complete!")
        return True
        
    except Exception as e:
        logger.error(f"Setup failed: {e}")
        logger.info("Please setup manually:")
        logger.info("  1. python3 -m venv .venv")
        logger.info("  2. source .venv/bin/activate")
        logger.info("  3. pip install -e .")
        logger.info("  4. pip install git+https://github.com/kfuku52/amalgkit")
        return False


def ensure_venv_activated(auto_setup: bool = True) -> bool:
    """
    Automatically activate virtual environment if needed.
    
    Args:
        auto_setup: If True, automatically setup venv if missing
        
    Returns:
        True if venv is activated, False if setup needed
    """
    venv_python = repo_root / ".venv" / "bin" / "python3"
    venv_dir = repo_root / ".venv"
    
    current_python = Path(sys.executable)
    
    # Check if already using venv Python
    try:
        current_python.relative_to(repo_root / ".venv")
        # Already in venv - ensure environment variables are set
        if "VIRTUAL_ENV" not in os.environ:
            os.environ["VIRTUAL_ENV"] = str(venv_dir)
            venv_bin = str(venv_dir / "bin")
            if venv_bin not in os.environ.get("PATH", ""):
                os.environ["PATH"] = f"{venv_bin}:{os.environ.get('PATH', '')}"
        return True
    except ValueError:
        pass
    
    # Check if venv exists
    if not venv_python.exists():
        if auto_setup:
            logger.info("Virtual environment not found - setting up automatically...")
            if not setup_venv_and_dependencies(auto_setup=True):
                logger.error("Failed to setup virtual environment")
                return False
        else:
            logger.error("Virtual environment not found!")
            logger.error("Run: python3 -m venv .venv && source .venv/bin/activate && pip install -e . && pip install git+https://github.com/kfuku52/amalgkit")
            return False
    
    # Activate venv by re-executing with venv Python
    new_env = os.environ.copy()
    new_env["VIRTUAL_ENV"] = str(venv_dir)
    venv_bin = str(venv_dir / "bin")
    new_env["PATH"] = f"{venv_bin}:{new_env.get('PATH', '')}"
    new_env.pop("PYTHONHOME", None)
    
    logger.info("🔄 Auto-activating virtual environment...")
    logger.info(f"Current Python: {current_python}")
    logger.info(f"Venv Python: {venv_python}")
    
    # Re-execute with venv Python
    os.execve(str(venv_python), [str(venv_python)] + sys.argv, new_env)
    # Never returns if venv activation succeeds


def check_environment() -> tuple[bool, list[str], list[str]]:
    """
    Check that all required tools are available.
    
    Returns:
        Tuple of (success: bool, missing: list[str], warnings: list[str])
    """
    missing = []
    warnings = []
    
    # Check virtual environment
    if os.environ.get("VIRTUAL_ENV") is None:
        venv_python = repo_root / ".venv" / "bin" / "python3"
        if not venv_python.exists():
            missing.append("Virtual environment (.venv)")
    
    # Check amalgkit
    if not shutil.which("amalgkit"):
        missing.append("amalgkit")
    
    # Check other tools (warnings, not critical)
    optional_tools = {
        "fasterq-dump": "SRA Toolkit",
        "kallisto": "kallisto",
        "fastp": "fastp",
        "seqkit": "seqkit",
        "wget": "wget (for ENA downloads)",
    }
    
    for tool, name in optional_tools.items():
        if not shutil.which(tool):
            warnings.append(f"{name} ({tool})")
    
    return len(missing) == 0, missing, warnings


def check_environment_or_exit(auto_setup: bool = True) -> None:
    """
    Check environment and exit if critical dependencies are missing.
    
    Args:
        auto_setup: If True, automatically setup missing components
    """
    # Ensure venv is activated
    if not ensure_venv_activated(auto_setup=auto_setup):
        if not auto_setup:
            sys.exit(1)
    
    # Check dependencies
    success, missing, warnings = check_environment()
    
    if not success:
        logger.error("=" * 80)
        logger.error("ENVIRONMENT CHECK FAILED")
        logger.error("=" * 80)
        logger.error("Missing required components:")
        for item in missing:
            logger.error(f"  ❌ {item}")
        
        if auto_setup:
            logger.info("Attempting automatic setup...")
            if "Virtual environment (.venv)" in missing:
                if not setup_venv_and_dependencies(auto_setup=True):
                    sys.exit(1)
                # Re-check after setup
                success, missing, warnings = check_environment()
        
        if not success:
            logger.error("\nSetup instructions:")
            logger.error("  1. python3 -m venv .venv")
            logger.error("  2. source .venv/bin/activate")
            logger.error("  3. pip install -e .")
            logger.error("  4. pip install git+https://github.com/kfuku52/amalgkit")
            sys.exit(1)
    
    if warnings:
        logger.warning("Optional tools not found (workflow may fail at certain steps):")
        for warning in warnings:
            logger.warning(f"  ⚠️  {warning}")
    
    logger.info("✓ Environment check passed")

