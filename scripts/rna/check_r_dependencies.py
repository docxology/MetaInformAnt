#!/usr/bin/env python3
"""Check and optionally install R dependencies for amalgkit workflows.

This script checks for R/Rscript availability and attempts installation if needed.
"""

import logging
import os
import shutil
import subprocess
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)s | %(message)s')
logger = logging.getLogger(__name__)


def check_rscript_available() -> bool:
    """Check if Rscript is available in PATH."""
    return shutil.which("Rscript") is not None


def install_r_debian() -> bool:
    """Attempt to install R on Debian/Ubuntu systems."""
    logger.info("Attempting to install R via apt...")
    
    try:
        # Check if we have sudo privileges
        result = subprocess.run(
            ["sudo", "-n", "true"],
            capture_output=True,
            timeout=5
        )
        
        if result.returncode != 0:
            logger.warning("No sudo privileges available - cannot auto-install R")
            return False
        
        # Update package lists
        logger.info("Updating package lists...")
        subprocess.run(
            ["sudo", "apt-get", "update", "-qq"],
            check=True,
            capture_output=True,
            timeout=120
        )
        
        # Install R
        logger.info("Installing r-base...")
        subprocess.run(
            ["sudo", "apt-get", "install", "-y", "-qq", "r-base"],
            check=True,
            timeout=300
        )
        
        logger.info("✅ R installation completed successfully")
        return True
        
    except subprocess.TimeoutExpired:
        logger.error("Installation timed out")
        return False
    except subprocess.CalledProcessError as e:
        logger.error(f"Installation failed: {e}")
        return False
    except Exception as e:
        logger.error(f"Unexpected error during installation: {e}")
        return False


def check_and_install_r(auto_install: bool = False) -> bool:
    """Check for R and optionally attempt installation.
    
    Args:
        auto_install: If True, attempt automatic installation on supported systems
        
    Returns:
        True if R is available (either pre-installed or successfully installed)
    """
    logger.info("=" * 80)
    logger.info("CHECKING R DEPENDENCIES")
    logger.info("=" * 80)
    
    # Check if Rscript is already available
    if check_rscript_available():
        rscript_path = shutil.which("Rscript")
        logger.info(f"✅ Rscript found: {rscript_path}")
        
        # Get R version
        try:
            result = subprocess.run(
                ["Rscript", "--version"],
                capture_output=True,
                text=True,
                timeout=5
            )
            version_info = result.stderr.strip() if result.stderr else result.stdout.strip()
            logger.info(f"   Version: {version_info}")
        except Exception as e:
            logger.warning(f"Could not determine R version: {e}")
        
        return True
    
    logger.warning("⚠️  Rscript not found in PATH")
    
    if not auto_install:
        logger.info("")
        logger.info("To install R manually:")
        logger.info("  Debian/Ubuntu: sudo apt-get install r-base")
        logger.info("  Fedora/RHEL:   sudo dnf install R")
        logger.info("  macOS:         brew install r")
        logger.info("")
        logger.info("Or run with --install flag to attempt automatic installation")
        return False
    
    # Attempt automatic installation
    logger.info("Attempting automatic installation...")
    
    # Detect OS
    if Path("/etc/debian_version").exists():
        logger.info("Detected Debian/Ubuntu system")
        return install_r_debian()
    else:
        logger.warning("Unsupported OS for automatic installation")
        logger.info("Please install R manually for your system")
        return False


def main():
    """Main entry point."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Check and optionally install R dependencies for amalgkit"
    )
    parser.add_argument(
        "--install",
        action="store_true",
        help="Attempt automatic installation if R is not found"
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress informational output"
    )
    
    args = parser.parse_args()
    
    if args.quiet:
        logging.getLogger().setLevel(logging.WARNING)
    
    # Check and optionally install
    r_available = check_and_install_r(auto_install=args.install)
    
    logger.info("=" * 80)
    if r_available:
        logger.info("✅ R DEPENDENCIES: SATISFIED")
        logger.info("=" * 80)
        return 0
    else:
        logger.warning("⚠️  R DEPENDENCIES: NOT SATISFIED")
        logger.info("=" * 80)
        return 1


if __name__ == "__main__":
    sys.exit(main())



