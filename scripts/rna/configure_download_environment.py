#!/usr/bin/env python3
"""Configure download environment for external drive usage.

Sets up SRA toolkit, temp directories, and environment variables to use
external drive instead of /tmp for all download operations.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.core.logging import get_logger
from metainformant.core.paths import expand_and_resolve

logger = get_logger("configure_download_env")


def find_repo_root(start_path: Path) -> Path:
    """Find repository root by looking for markers."""
    current = start_path.resolve()
    while current != current.parent:
        if any((current / marker).exists() for marker in [".git", "pyproject.toml", ".cursorrules"]):
            return current
        current = current.parent
    return start_path.resolve()


def configure_vdb_config(sra_temp_dir: Path) -> bool:
    """Configure vdb-config to use external drive for SRA repository."""
    logger.info("Configuring vdb-config for external drive...")
    
    try:
        # Set repository root
        cmd1 = [
            "vdb-config",
            "-s",
            f"/repository/user/main/public/root={sra_temp_dir}",
        ]
        result1 = subprocess.run(cmd1, capture_output=True, text=True, check=False)
        
        # Set flat volume
        cmd2 = [
            "vdb-config",
            "-s",
            f"/repository/user/main/public/apps/sra/volumes/sraFlat={sra_temp_dir}",
        ]
        result2 = subprocess.run(cmd2, capture_output=True, text=True, check=False)
        
        if result1.returncode == 0 and result2.returncode == 0:
            logger.info(f"✓ vdb-config configured to use: {sra_temp_dir}")
            return True
        else:
            logger.warning(f"vdb-config configuration had issues (codes: {result1.returncode}, {result2.returncode})")
            if result1.stderr:
                logger.debug(f"vdb-config stderr: {result1.stderr}")
            return False
    except FileNotFoundError:
        logger.warning("vdb-config not found in PATH - skipping configuration")
        return False
    except Exception as e:
        logger.error(f"Error configuring vdb-config: {e}")
        return False


def setup_environment_variables(repo_root: Path, work_dir: Path) -> dict[str, str]:
    """Set up environment variables for external drive temp usage."""
    logger.info("Setting up environment variables...")
    
    # Create temp directories on external drive
    sra_temp_dir = work_dir / "fastq" / "temp" / "sra"
    general_temp_dir = repo_root / "output" / ".tmp"
    
    sra_temp_dir.mkdir(parents=True, exist_ok=True)
    general_temp_dir.mkdir(parents=True, exist_ok=True)
    
    env_vars = {
        "TMPDIR": str(general_temp_dir),
        "TEMP": str(general_temp_dir),
        "TMP": str(general_temp_dir),
        "NCBI_SRA_REPOSITORY": str(sra_temp_dir),
    }
    
    logger.info(f"✓ Temp directories configured:")
    logger.info(f"  TMPDIR: {general_temp_dir}")
    logger.info(f"  NCBI_SRA_REPOSITORY: {sra_temp_dir}")
    
    return env_vars


def verify_disk_space(path: Path) -> tuple[bool, str]:
    """Verify sufficient disk space available."""
    try:
        import shutil
        total, used, free = shutil.disk_usage(path)
        free_gb = free / (1024**3)
        
        if free_gb < 50:
            return False, f"Only {free_gb:.1f}GB free - need at least 50GB for downloads"
        return True, f"{free_gb:.1f}GB free (sufficient)"
    except Exception as e:
        return False, f"Could not check disk space: {e}"


def main() -> int:
    """Main entry point."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Configure download environment for external drive",
    )
    parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="Path to workflow config file",
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Apply configuration (set environment variables)",
    )
    
    args = parser.parse_args()
    
    config_path = expand_and_resolve(args.config)
    if not config_path.exists():
        logger.error(f"Config file not found: {config_path}")
        return 1
    
    repo_root = find_repo_root(config_path)
    logger.info(f"Repository root: {repo_root}")
    
    # Load config to get work_dir
    try:
        from metainformant.rna.workflow import load_workflow_config
        
        config = load_workflow_config(config_path)
        work_dir = config.work_dir
        logger.info(f"Work directory: {work_dir}")
    except Exception as e:
        logger.error(f"Failed to load config: {e}")
        return 1
    
    # Verify disk space
    ok, msg = verify_disk_space(work_dir)
    if not ok:
        logger.warning(f"⚠ Disk space warning: {msg}")
    else:
        logger.info(f"✓ Disk space: {msg}")
    
    # Configure vdb-config
    sra_temp_dir = work_dir / "fastq" / "temp" / "sra"
    sra_temp_dir.mkdir(parents=True, exist_ok=True)
    configure_vdb_config(sra_temp_dir)
    
    # Set up environment variables
    env_vars = setup_environment_variables(repo_root, work_dir)
    
    if args.apply:
        logger.info("Applying environment variables to current shell...")
        for key, value in env_vars.items():
            os.environ[key] = value
            logger.info(f"  export {key}={value}")
        logger.info("\n✓ Environment variables set for current process")
        logger.info("To make permanent, add these to your shell profile:")
        for key, value in env_vars.items():
            print(f"export {key}={value}")
    else:
        logger.info("\nEnvironment variables to set:")
        for key, value in env_vars.items():
            print(f"export {key}={value}")
        logger.info("\nRun with --apply to set for current process")
    
    logger.info("\n✓ Configuration complete!")
    return 0


if __name__ == "__main__":
    sys.exit(main())

