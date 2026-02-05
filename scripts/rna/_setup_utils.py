"""
Common setup utilities for RNA scripts.

Provides automatic virtual environment setup, dependency installation, and environment validation.
Uses uv for all environment management to handle filesystem limitations (e.g., ext6 without symlinks).
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional

# Add src to path before imports
repo_root = Path(__file__).parent.parent.parent.resolve()
src_path = repo_root / "src"
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

try:
    from metainformant.core.utils.logging import get_logger
except ImportError:
    # If import fails, we might need to debug sys.path or structure
    # Fallback logger or print
    import logging as std_logging

    def get_logger(name):
        return std_logging.getLogger(name)


try:
    from metainformant.core.io.disk import get_recommended_temp_dir
except ImportError:
    # Fallback if disk module moved or missing
    def get_recommended_temp_dir(repo_root=None):
        import tempfile

        return Path(tempfile.gettempdir())


logger = get_logger("rna_setup")


# Check if uv is available
def _check_uv() -> bool:
    """Check if uv is available."""
    return shutil.which("uv") is not None


def _find_venv_location() -> tuple[Path, Path]:
    """
    Find suitable location for virtual environment.

    Tries .venv in repo root first. If that fails due to symlink issues,
    tries /tmp/metainformant_venv and creates symlink reference.

    Returns:
        Tuple of (venv_dir, venv_python)
    """
    venv_dir = repo_root / ".venv"
    venv_python = venv_dir / "bin" / "python3"

    # If venv exists and is valid, use it
    if venv_python.exists():
        return venv_dir, venv_python

    # Try alternative location if .venv creation will fail
    alt_venv_dir = Path("/tmp/metainformant_venv")
    alt_venv_python = alt_venv_dir / "bin" / "python3"

    # If alternative exists, use it
    if alt_venv_python.exists():
        logger.info(f"Using alternative venv location: {alt_venv_dir}")
        return alt_venv_dir, alt_venv_python

    # Default to repo root
    return venv_dir, venv_python


def setup_venv_and_dependencies(auto_setup: bool = True) -> bool:
    """
    Automatically setup virtual environment and install dependencies using uv.

    Uses uv for all operations to handle filesystem limitations (e.g., ext6 without symlinks).

    Args:
        auto_setup: If True, automatically create venv and install dependencies

    Returns:
        True if venv exists and is ready, False otherwise
    """
    venv_dir, venv_python = _find_venv_location()

    # Check if venv exists and is valid
    if venv_python.exists():
        logger.info(f"✓ Virtual environment found at {venv_dir}")
        # Verify amalgkit and pyyaml are installed
        try:
            amalgkit_result = subprocess.run(
                [str(venv_python), "-c", "import amalgkit; print('ok')"], capture_output=True, text=True, timeout=5
            )
            yaml_result = subprocess.run(
                [str(venv_python), "-c", "import yaml; print('ok')"], capture_output=True, text=True, timeout=5
            )
            if amalgkit_result.returncode == 0 and yaml_result.returncode == 0:
                return True
            else:
                logger.info("Venv exists but missing dependencies, installing...")
        except Exception:
            logger.info("Venv exists but verifying dependencies...")

    if not auto_setup:
        logger.error("Virtual environment not found!")
        if _check_uv():
            logger.error(
                "Run: uv venv .venv && uv pip install -e . --python .venv/bin/python3 && uv pip install git+https://github.com/kfuku52/amalgkit --python .venv/bin/python3"
            )
        else:
            logger.error(
                "Run: python3 -m venv .venv && source .venv/bin/activate && pip install -e . && pip install git+https://github.com/kfuku52/amalgkit"
            )
        return False

    if not _check_uv():
        logger.error("uv not found! Please install uv first:")
        logger.error("  curl -LsSf https://astral.sh/uv/install.sh | sh")
        logger.error("Or use: pip install uv")
        return False

    # Auto-setup venv with uv
    logger.info("Setting up virtual environment with uv...")

    # Try creating venv in repo root first
    try:
        if venv_dir.exists():
            logger.info(f"Clearing existing venv at {venv_dir}...")
            result = subprocess.run(
                ["uv", "venv", "--clear", str(venv_dir)], capture_output=True, text=True, timeout=30
            )
        else:
            logger.info(f"Creating venv at {venv_dir}...")
            result = subprocess.run(["uv", "venv", str(venv_dir)], capture_output=True, text=True, timeout=30)

        if result.returncode != 0:
            # Try alternative location if repo root fails (symlink issues)
            logger.warning(f"Could not create venv at {venv_dir}: {result.stderr}")
            logger.info("Trying alternative location: /tmp/metainformant_venv")
            alt_venv_dir = Path("/tmp/metainformant_venv")
            result = subprocess.run(["uv", "venv", str(alt_venv_dir)], capture_output=True, text=True, timeout=30)
            if result.returncode == 0:
                venv_dir = alt_venv_dir
                venv_python = alt_venv_dir / "bin" / "python3"
                logger.info(f"✓ Virtual environment created at {venv_dir}")
            else:
                logger.error(f"Failed to create venv: {result.stderr}")
                return False
        else:
            venv_python = venv_dir / "bin" / "python3"
            logger.info("✓ Virtual environment created")
    except subprocess.TimeoutExpired:
        logger.error("Venv creation timed out")
        return False
    except Exception as e:
        logger.error(f"Venv creation failed: {e}")
        return False

    # Install dependencies with uv pip
    logger.info("Installing dependencies with uv pip...")

    try:
        # Get temp directory (prefers external drive if available)
        temp_dir = get_recommended_temp_dir(repo_root)
        logger.info(f"Using temporary directory: {temp_dir}")

        # Install metainformant
        # Use /tmp/uv-cache for cache to avoid symlink issues on ext6 filesystem
        # The output/.uv-cache location doesn't support symlinks needed by uv
        uv_cache_dir = Path("/tmp/uv-cache")
        uv_cache_dir.mkdir(parents=True, exist_ok=True)
        logger.info("Installing metainformant...")
        env = os.environ.copy()
        env["UV_CACHE_DIR"] = str(uv_cache_dir)
        env["TMPDIR"] = str(temp_dir)
        result = subprocess.run(
            ["uv", "pip", "install", "-e", str(repo_root), "--python", str(venv_python)],
            capture_output=True,
            text=True,
            timeout=300,
            env=env,
        )
        if result.returncode != 0:
            logger.error(f"Failed to install metainformant: {result.stderr}")
            return False
        logger.info("✓ metainformant installed")

        # Install amalgkit
        # Use /tmp/uv-cache for cache to avoid symlink issues on ext6 filesystem
        temp_dir = get_recommended_temp_dir(repo_root)
        uv_cache_dir = Path("/tmp/uv-cache")
        uv_cache_dir.mkdir(parents=True, exist_ok=True)
        logger.info("Installing amalgkit...")
        env = os.environ.copy()
        env["UV_CACHE_DIR"] = str(uv_cache_dir)
        env["TMPDIR"] = str(temp_dir)
        result = subprocess.run(
            ["uv", "pip", "install", "git+https://github.com/kfuku52/amalgkit", "--python", str(venv_python)],
            capture_output=True,
            text=True,
            timeout=300,
            env=env,
        )
        if result.returncode != 0:
            logger.error(f"Failed to install amalgkit: {result.stderr}")
            return False
        logger.info("✓ amalgkit installed")

        logger.info("✓ Setup complete!")

        # Update venv location tracking
        if venv_dir != repo_root / ".venv":
            logger.info(f"Note: Using venv at {venv_dir} (repo .venv has symlink limitations)")
            logger.info(f"Scripts will use this venv automatically")

        return True

    except subprocess.TimeoutExpired:
        logger.error("Package installation timed out")
        return False
    except Exception as e:
        logger.error(f"Setup failed: {e}")
        if _check_uv():
            logger.info("Please setup manually:")
            logger.info("  1. uv venv .venv")
            logger.info("  2. uv pip install -e . --python .venv/bin/python3")
            logger.info("  3. uv pip install git+https://github.com/kfuku52/amalgkit --python .venv/bin/python3")
        else:
            logger.info("Please setup manually:")
            logger.info("  1. python3 -m venv .venv")
            logger.info("  2. source .venv/bin/activate")
            logger.info("  3. pip install -e .")
            logger.info("  4. pip install git+https://github.com/kfuku52/amalgkit")
        return False


def ensure_venv_activated(auto_setup: bool = True) -> bool:
    """
    Automatically activate virtual environment if needed.

    Uses uv-based venv discovery to handle filesystem limitations.

    Args:
        auto_setup: If True, automatically setup venv if missing

    Returns:
        True if venv is activated, False if setup needed
    """
    venv_dir, venv_python = _find_venv_location()
    current_python = Path(sys.executable)

    # Check if already using venv Python (check both locations)
    try:
        current_python.relative_to(venv_dir)
        # Already in venv - ensure environment variables are set
        if "VIRTUAL_ENV" not in os.environ:
            os.environ["VIRTUAL_ENV"] = str(venv_dir)
            venv_bin = str(venv_dir / "bin")
            if venv_bin not in os.environ.get("PATH", ""):
                os.environ["PATH"] = f"{venv_bin}:{os.environ.get('PATH', '')}"
        return True
    except ValueError:
        # Also check alternative location
        alt_venv_dir = Path("/tmp/metainformant_venv")
        try:
            current_python.relative_to(alt_venv_dir)
            if "VIRTUAL_ENV" not in os.environ:
                os.environ["VIRTUAL_ENV"] = str(alt_venv_dir)
                venv_bin = str(alt_venv_dir / "bin")
                if venv_bin not in os.environ.get("PATH", ""):
                    os.environ["PATH"] = f"{venv_bin}:{os.environ.get('PATH', '')}"
            return True
        except ValueError:
            pass

    # Check if venv exists
    if not venv_python.exists():
        if auto_setup:
            logger.info("Virtual environment not found - setting up automatically with uv...")
            if not setup_venv_and_dependencies(auto_setup=True):
                logger.error("Failed to setup virtual environment")
                return False
            # Re-find venv location after setup
            venv_dir, venv_python = _find_venv_location()
        else:
            logger.error("Virtual environment not found!")
            if _check_uv():
                logger.error(
                    "Run: uv venv .venv && uv pip install -e . --python .venv/bin/python3 && uv pip install git+https://github.com/kfuku52/amalgkit --python .venv/bin/python3"
                )
            else:
                logger.error(
                    "Run: python3 -m venv .venv && source .venv/bin/activate && pip install -e . && pip install git+https://github.com/kfuku52/amalgkit"
                )
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
    logger.info(f"Venv Location: {venv_dir}")

    # Re-execute with venv Python
    os.execve(str(venv_python), [str(venv_python)] + sys.argv, new_env)
    # Never returns if venv activation succeeds


def check_environment() -> tuple[bool, list[str], list[str]]:
    """
    Check that all required tools are available.

    Uses uv-based venv discovery to handle filesystem limitations.

    Returns:
        Tuple of (success: bool, missing: list[str], warnings: list[str])
    """
    missing = []
    warnings = []

    # Check virtual environment (check both locations)
    venv_dir, venv_python = _find_venv_location()
    if os.environ.get("VIRTUAL_ENV") is None:
        if not venv_python.exists():
            missing.append("Virtual environment (.venv)")

    # Check amalgkit (check in venv first, then PATH)
    amalgkit_found = False
    if venv_python.exists():
        try:
            result = subprocess.run(
                [str(venv_python), "-c", "import amalgkit; print('ok')"], capture_output=True, text=True, timeout=5
            )
            if result.returncode == 0:
                amalgkit_found = True
        except Exception:
            pass

    if not amalgkit_found and not shutil.which("amalgkit"):
        missing.append("amalgkit")

    # Check pyyaml (required for config parsing)
    if venv_python.exists():
        try:
            result = subprocess.run(
                [str(venv_python), "-c", "import yaml; print('ok')"], capture_output=True, text=True, timeout=5
            )
            if result.returncode != 0:
                missing.append("pyyaml")
        except Exception:
            missing.append("pyyaml")

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
            logger.error("\nSetup instructions (using uv):")
            if _check_uv():
                logger.error("  1. uv venv .venv")
                logger.error("  2. uv pip install -e . --python .venv/bin/python3")
                logger.error("  3. uv pip install git+https://github.com/kfuku52/amalgkit --python .venv/bin/python3")
            else:
                logger.error("  First install uv:")
                logger.error("    curl -LsSf https://astral.sh/uv/install.sh | sh")
                logger.error("  Then:")
                logger.error("  1. uv venv .venv")
                logger.error("  2. uv pip install -e . --python .venv/bin/python3")
                logger.error("  3. uv pip install git+https://github.com/kfuku52/amalgkit --python .venv/bin/python3")
            sys.exit(1)

    if warnings:
        logger.warning("Optional tools not found (workflow may fail at certain steps):")
        for warning in warnings:
            logger.warning(f"  ⚠️  {warning}")

    logger.info("✓ Environment check passed")
