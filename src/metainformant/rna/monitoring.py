"""Workflow monitoring and status checking functions.

This module provides comprehensive functions for checking workflow status,
sample completion, and progress tracking across all RNA-seq workflows.
"""

from __future__ import annotations

import subprocess
import time
from collections.abc import Mapping
from pathlib import Path
from typing import Any

# Handle both relative and absolute imports
try:
    from ...core.io import read_delimited
    from ...core.logging import get_logger
except ImportError:
    from metainformant.core.io import read_delimited
    from metainformant.core.logging import get_logger

logger = get_logger(__name__)


def count_quantified_samples(config_path: Path) -> tuple[int, int]:
    """Count quantified and total samples for a species.
    
    Args:
        config_path: Path to species workflow config file
        
    Returns:
        Tuple of (quantified_count, total_count)
    """
    try:
        from metainformant.rna.workflow import load_workflow_config
        cfg = load_workflow_config(config_path)
        quant_dir = Path(cfg.per_step.get("quant", {}).get("out_dir", cfg.work_dir / "quant"))
        
        # Count quantified
        quantified = 0
        if quant_dir.exists():
            quantified = len([
                d for d in quant_dir.iterdir()
                if d.is_dir() and (d / "abundance.tsv").exists()
            ])
        
        # Count total from metadata
        metadata_file = cfg.work_dir / "metadata" / "metadata.tsv"
        if not metadata_file.exists():
            metadata_file = cfg.work_dir / "metadata" / "metadata.filtered.tissue.tsv"
        
        total = 0
        if metadata_file.exists():
            try:
                rows = list(read_delimited(metadata_file, delimiter="\t"))
                total = len([r for r in rows if r.get("run")])
            except Exception:
                pass
        
        return quantified, total
    except Exception:
        return 0, 0


def get_sample_status(config_path: Path, sample_id: str) -> dict[str, Any]:
    """Get detailed status for a single sample.
    
    Args:
        config_path: Path to species workflow config file
        sample_id: SRA accession ID (e.g., "SRR1234567")
        
    Returns:
        Dictionary with status information:
        - quantified: bool
        - has_fastq: bool
        - has_sra: bool
        - is_downloading: bool
        - status: str ("quantified", "downloading", "has_fastq", "has_sra", "undownloaded")
    """
    try:
        from metainformant.rna.workflow import load_workflow_config
        cfg = load_workflow_config(config_path)
        fastq_dir = Path(cfg.per_step.get("getfastq", {}).get("out_dir", cfg.work_dir / "fastq"))
        quant_dir = Path(cfg.per_step.get("quant", {}).get("out_dir", cfg.work_dir / "quant"))
        
        # Check if quantified
        abundance_file = quant_dir / sample_id / "abundance.tsv"
        is_quantified = abundance_file.exists()
        
        # Check for FASTQ/SRA files
        has_fastq = False
        has_sra = False
        
        # Check getfastq subdirectory
        sample_dir_getfastq = fastq_dir / "getfastq" / sample_id
        if sample_dir_getfastq.exists():
            has_fastq = any(sample_dir_getfastq.glob("*.fastq*"))
            has_sra = any(sample_dir_getfastq.glob("*.sra"))
        
        # Check direct structure
        sample_dir_direct = fastq_dir / sample_id
        if sample_dir_direct.exists():
            has_fastq = has_fastq or any(sample_dir_direct.glob("*.fastq*"))
            has_sra = has_sra or any(sample_dir_direct.glob("*.sra"))
        
        # Check if actively downloading
        is_downloading = sample_id in check_active_downloads()
        
        # Determine status
        if is_quantified:
            status = "quantified"
        elif is_downloading:
            status = "downloading"
        elif has_fastq:
            status = "has_fastq"
        elif has_sra:
            status = "has_sra"
        else:
            status = "undownloaded"
        
        return {
            "quantified": is_quantified,
            "has_fastq": has_fastq,
            "has_sra": has_sra,
            "is_downloading": is_downloading,
            "status": status,
        }
    except Exception as e:
        logger.debug(f"Error getting sample status: {e}")
        return {
            "quantified": False,
            "has_fastq": False,
            "has_sra": False,
            "is_downloading": False,
            "status": "error",
        }


def analyze_species_status(config_path: Path) -> dict[str, Any]:
    """Comprehensive analysis of species workflow status.
    
    Args:
        config_path: Path to species workflow config file
        
    Returns:
        Dictionary with comprehensive status information:
        - total_in_metadata: int
        - quantified: int
        - quantified_and_deleted: int
        - quantified_not_deleted: int
        - downloading: int
        - failed_download: int
        - undownloaded: int
        - categories: dict mapping category -> list of sample_ids
    """
    try:
        from metainformant.rna.workflow import load_workflow_config
        cfg = load_workflow_config(config_path)
        fastq_dir = Path(cfg.per_step.get("getfastq", {}).get("out_dir", cfg.work_dir / "fastq"))
        quant_dir = Path(cfg.per_step.get("quant", {}).get("out_dir", cfg.work_dir / "quant"))
        
        # Get metadata
        metadata_file = cfg.work_dir / "metadata" / "metadata.tsv"
        if not metadata_file.exists():
            metadata_file = cfg.work_dir / "metadata" / "metadata.filtered.tissue.tsv"
        
        if not metadata_file.exists():
            return {
                "total_in_metadata": 0,
                "quantified": 0,
                "quantified_and_deleted": 0,
                "quantified_not_deleted": 0,
                "downloading": 0,
                "failed_download": 0,
                "undownloaded": 0,
                "categories": {},
            }
        
        rows = list(read_delimited(metadata_file, delimiter="\t"))
        sample_ids = [row.get("run") for row in rows if row.get("run")]
        
        # Check active downloads
        active_downloads = check_active_downloads()
        
        # Categorize samples
        categories = {
            "quantified_and_deleted": [],
            "quantified_not_deleted": [],
            "downloading": [],
            "failed_download": [],
            "undownloaded": [],
        }
        
        for sample_id in sample_ids:
            if not sample_id:
                continue
            
            # Check if quantified
            abundance_file = quant_dir / sample_id / "abundance.tsv"
            is_quantified = abundance_file.exists()
            
            # Check for FASTQ/SRA files
            has_fastq = False
            has_sra = False
            
            sample_dir_getfastq = fastq_dir / "getfastq" / sample_id
            if sample_dir_getfastq.exists():
                has_fastq = any(sample_dir_getfastq.glob("*.fastq*"))
                has_sra = any(sample_dir_getfastq.glob("*.sra"))
            
            sample_dir_direct = fastq_dir / sample_id
            if sample_dir_direct.exists():
                has_fastq = has_fastq or any(sample_dir_direct.glob("*.fastq*"))
                has_sra = has_sra or any(sample_dir_direct.glob("*.sra"))
            
            has_files = has_fastq or has_sra
            
            # Categorize
            if sample_id in active_downloads:
                categories["downloading"].append(sample_id)
            elif is_quantified and not has_files:
                categories["quantified_and_deleted"].append(sample_id)
            elif is_quantified and has_files:
                categories["quantified_not_deleted"].append(sample_id)
            elif has_files:
                categories["failed_download"].append(sample_id)
            else:
                categories["undownloaded"].append(sample_id)
        
        return {
            "total_in_metadata": len(sample_ids),
            "quantified": len(categories["quantified_and_deleted"]) + len(categories["quantified_not_deleted"]),
            "quantified_and_deleted": len(categories["quantified_and_deleted"]),
            "quantified_not_deleted": len(categories["quantified_not_deleted"]),
            "downloading": len(categories["downloading"]),
            "failed_download": len(categories["failed_download"]),
            "undownloaded": len(categories["undownloaded"]),
            "categories": categories,
        }
    except Exception as e:
        logger.warning(f"Error analyzing species status: {e}")
        return {
            "total_in_metadata": 0,
            "quantified": 0,
            "quantified_and_deleted": 0,
            "quantified_not_deleted": 0,
            "downloading": 0,
            "failed_download": 0,
            "undownloaded": 0,
            "categories": {},
        }


def find_unquantified_samples(config_path: Path) -> list[str]:
    """Find all unquantified samples.
    
    Args:
        config_path: Path to species workflow config file
        
    Returns:
        List of sample IDs that are not quantified
    """
    try:
        from metainformant.rna.workflow import load_workflow_config
        cfg = load_workflow_config(config_path)
        quant_dir = Path(cfg.per_step.get("quant", {}).get("out_dir", cfg.work_dir / "quant"))
        
        # Get metadata
        metadata_file = cfg.work_dir / "metadata" / "metadata.tsv"
        if not metadata_file.exists():
            metadata_file = cfg.work_dir / "metadata" / "metadata.filtered.tissue.tsv"
        
        if not metadata_file.exists():
            return []
        
        rows = list(read_delimited(metadata_file, delimiter="\t"))
        sample_ids = [row.get("run") for row in rows if row.get("run")]
        
        # Filter to unquantified
        unquantified = []
        for sample_id in sample_ids:
            if not sample_id:
                continue
            abundance_file = quant_dir / sample_id / "abundance.tsv"
            if not abundance_file.exists():
                unquantified.append(sample_id)
        
        return unquantified
    except Exception as e:
        logger.warning(f"Error finding unquantified samples: {e}")
        return []


def check_active_downloads() -> set[str]:
    """Check for samples currently being downloaded.
    
    Returns:
        Set of sample IDs that are actively downloading
    """
    active_samples = set()
    
    # Check running amalgkit getfastq processes
    try:
        result = subprocess.run(
            ["ps", "aux"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        import re
        for line in result.stdout.split("\n"):
            if "amalgkit getfastq" in line:
                # Extract all SRR IDs from the line
                matches = re.findall(r"SRR\d+", line)
                active_samples.update(matches)
                
                # Also check for --id parameter
                if "--id" in line:
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if part == "--id" and i + 1 < len(parts):
                            active_samples.add(parts[i + 1])
    except Exception:
        pass
    
    # Also check for recent fastq directory activity
    try:
        from ...core.paths import get_repo_root
        repo_root = get_repo_root()
        fastq_base = repo_root / "output" / "amalgkit"
        
        # Find directories modified in last 10 minutes
        current_time = time.time()
        for species_dir in fastq_base.glob("*/fastq"):
            if species_dir.exists():
                # Check for recently modified sample directories
                for sample_dir in species_dir.glob("getfastq/SRR*"):
                    if sample_dir.exists():
                        # Check if directory was modified recently
                        mtime = sample_dir.stat().st_mtime
                        if current_time - mtime < 600:  # 10 minutes
                            sample_id = sample_dir.name
                            active_samples.add(sample_id)
    except Exception:
        pass
    
    return active_samples


def check_workflow_progress(config_path: Path) -> dict[str, Any]:
    """Get workflow progress summary for a species.
    
    Args:
        config_path: Path to species workflow config file
        
    Returns:
        Dictionary with progress information:
        - quantified: int
        - total: int
        - percentage: float
        - remaining: int
        - downloading: int (number of samples currently downloading)
        - has_files: int (number of samples with downloaded files but not quantified)
    """
    quantified, total = count_quantified_samples(config_path)
    percentage = (quantified / total * 100) if total > 0 else 0.0
    remaining = total - quantified
    
    # Get downloading status
    active_downloads = check_active_downloads()
    downloading = len(active_downloads)
    
    # Count samples with files but not quantified
    try:
        from metainformant.rna.workflow import load_workflow_config
        cfg = load_workflow_config(config_path)
        fastq_dir = Path(cfg.per_step.get("getfastq", {}).get("out_dir", cfg.work_dir / "fastq"))
        quant_dir = Path(cfg.per_step.get("quant", {}).get("out_dir", cfg.work_dir / "quant"))
        
        has_files = 0
        getfastq_dir = fastq_dir / "getfastq"
        if not getfastq_dir.exists():
            getfastq_dir = fastq_dir
        
        if getfastq_dir.exists():
            for sample_dir in getfastq_dir.iterdir():
                if sample_dir.is_dir() and sample_dir.name.startswith("SRR"):
                    # Check if has files but not quantified
                    has_fastq = any(sample_dir.glob("*.fastq*"))
                    has_sra = any(sample_dir.glob("*.sra"))
                    if (has_fastq or has_sra):
                        abundance_file = quant_dir / sample_dir.name / "abundance.tsv"
                        if not abundance_file.exists() and sample_dir.name not in active_downloads:
                            has_files += 1
    except Exception:
        has_files = 0
    
    return {
        "quantified": quantified,
        "total": total,
        "percentage": percentage,
        "remaining": remaining,
        "downloading": downloading,
        "has_files": has_files,
    }


def assess_all_species_progress(
    config_dir: Path,
    *,
    repo_root: Path | None = None,
) -> dict[str, dict[str, Any]]:
    """Assess progress for all species in config directory.

    Args:
        config_dir: Directory containing amalgkit config files
        repo_root: Repository root directory (optional)

    Returns:
        Dictionary mapping species_id -> progress information
    """
    if repo_root is None:
        repo_root = Path.cwd()

    results = {}

    for config_file in sorted(config_dir.glob("amalgkit_*.yaml")):
        if "template" in config_file.stem.lower() or "test" in config_file.stem.lower():
            continue

        species_id = config_file.stem.replace("amalgkit_", "")

        try:
            progress = check_workflow_progress(config_file)
            results[species_id] = progress
        except Exception as e:
            logger.warning(f"Error assessing {species_id}: {e}")
            results[species_id] = {
                "quantified": 0,
                "total": 0,
                "percentage": 0.0,
                "remaining": 0,
                "error": str(e),
            }

    return results


def initialize_progress_tracking(
    config_path: Path,
    *,
    tracker=None,
) -> dict[str, Any]:
    """Initialize progress tracking for a species.

    Args:
        config_path: Path to species workflow config file
        tracker: Optional progress tracker instance

    Returns:
        Dictionary with initialization results
    """
    try:
        from metainformant.rna.workflow import load_workflow_config

        cfg = load_workflow_config(config_path)
        metadata_file = cfg.work_dir / "metadata" / "metadata.tsv"
        if not metadata_file.exists():
            metadata_file = cfg.work_dir / "metadata" / "metadata.filtered.tissue.tsv"

        if not metadata_file.exists():
            return {"success": False, "error": "No metadata file found"}

        rows = list(read_delimited(metadata_file, delimiter="\t"))
        sample_ids = [row.get("run") for row in rows if row.get("run")]

        if not sample_ids:
            return {"success": False, "error": "No samples in metadata"}

        # If tracker is provided, initialize it
        if tracker:
            species_id = config_path.stem.replace("amalgkit_", "")
            if species_id not in tracker.state:
                tracker.initialize_species(species_id, len(sample_ids), sample_ids)
                tracker._save_state()

        return {
            "success": True,
            "total_samples": len(sample_ids),
            "sample_ids": sample_ids,
        }
    except Exception as e:
        logger.error(f"Error initializing progress tracking: {e}")
        return {"success": False, "error": str(e)}

