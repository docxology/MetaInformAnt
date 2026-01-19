"""RNA-seq workflow monitoring and status tracking.

This module provides functions for monitoring workflow progress,
checking completion status, and analyzing workflow performance.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from metainformant.core import io, logging

logger = logging.get_logger(__name__)


def check_workflow_progress(config_path: Path) -> Dict[str, Any]:
    """Check the progress of RNA-seq workflow.

    Args:
        config_path: Path to workflow configuration file

    Returns:
        Dictionary with progress information
    """
    from metainformant.rna.engine.workflow import load_workflow_config
    try:
        config = load_workflow_config(config_path)
        status = analyze_species_status(config.work_dir)
        
        # Add required top-level keys for tests
        progress_info = status.copy()
        progress_info["total"] = status.get("total_in_metadata", 0)
        progress_info["quantified"] = status.get("quantified", 0)
        
        total = progress_info["total"]
        if total > 0:
            progress_info["percentage"] = (progress_info["quantified"] / total) * 100.0
        else:
            progress_info["percentage"] = 0.0
            
        progress_info["remaining"] = progress_info["total"] - progress_info["quantified"]
        progress_info["has_files"] = progress_info["quantified"] > 0 or progress_info["percentage"] > 0
        
        # downloading is expected by tests
        active = check_active_downloads()
        progress_info["downloading"] = len(active)
        
        return progress_info
    except Exception as e:
        logger.error(f"Failed to check workflow progress: {e}")
        return {
            "quantified": 0, "total": 0, "percentage": 0.0, 
            "remaining": 0, "downloading": 0, "has_files": False
        }


def assess_species_progress(species_dir: Path | str) -> Dict[str, Any]:
    """Alias for analyze_species_status."""
    return analyze_species_status(species_dir)


def count_quantified_samples(config_path: Path | str) -> Tuple[int, int]:
    """Count quantified samples in a workflow.

    Args:
        config_path: Path to workflow configuration file

    Returns:
        Tuple of (quantified_count, total_count)
    """
    from metainformant.rna.engine.workflow import load_workflow_config
    try:
        path = Path(config_path)
        if path.is_dir():
            # If passed a directory, assume it's the work_dir
            work_dir = path
            total = 0 # Cannot know total without config/metadata
        else:
            config = load_workflow_config(path)
            work_dir = config.work_dir
            total = 0
            metadata_file = work_dir / "metadata" / "metadata.tsv"
            if not metadata_file.exists():
                 metadata_file = work_dir / "metadata" / "metadata_selected.tsv"
                 
            if metadata_file.exists():
                import csv
                with open(metadata_file, "r") as f:
                    total = sum(1 for _ in csv.DictReader(f, delimiter="\t"))
        
        quant_dir = work_dir / "quant"
        quantified = 0
        if quant_dir.exists():
            for sample_dir in quant_dir.iterdir():
                if sample_dir.is_dir() and (sample_dir / "abundance.tsv").exists():
                    quantified += 1
        
        return quantified, total
    except Exception as e:
        logger.error(f"Failed to count quantified samples: {e}")
        return 0, 0


def analyze_species_status(species_dir: Path | str) -> Dict[str, Any]:
    """Analyze the status of a species in RNA-seq workflow.

    Args:
        species_dir: Species directory or config file

    Returns:
        Dictionary with status information
    """
    path = Path(species_dir)
    if path.is_file():
         from metainformant.rna.engine.workflow import load_workflow_config
         try:
             config = load_workflow_config(path)
             species_dir = config.work_dir
         except Exception:
             species_dir = path.parent
    else:
        species_dir = path

    status = {
        'species': species_dir.name,
        'completed': False,
        'failed': False,
        'steps_completed': [],
        'steps_failed': [],
        'last_updated': None,
        'progress_percentage': 0.0,
        'total_in_metadata': 0,
        'quantified': 0,
        'categories': {
            'quantified_and_deleted': 0,
            'undownloaded': 0,
        }
    }
    
    # Try to find metadata to get total count
    metadata_file = species_dir / "metadata" / "metadata.tsv"
    if not metadata_file.exists():
         metadata_file = species_dir / "metadata" / "metadata_selected.tsv"
         
    if metadata_file.exists():
        import csv
        try:
            with open(metadata_file, "r") as f:
                reader = csv.DictReader(f, delimiter="\t")
                rows = list(reader)
                status['total_in_metadata'] = len(rows)
                status['categories']['undownloaded'] = len(rows)
        except Exception:
            pass

    # Check for completion markers
    sanity_file = species_dir / "sanity_check.txt"
    if sanity_file.exists():
        status['completed'] = True
        status['progress_percentage'] = 100.0

        # Try to get completion time
        try:
            mtime = sanity_file.stat().st_mtime
            status['last_updated'] = mtime
        except Exception:
            pass
    else:
        # Check individual step completion
        workflow_steps = [
            'metadata', 'integrate', 'config', 'select',
            'getfastq', 'quant', 'merge', 'cstmm', 'curate', 'csca'
        ]

        completed_steps = 0
        total_steps = len(workflow_steps)

        for step in workflow_steps:
            # Check for step-specific output files
            step_outputs = _get_step_output_files(species_dir, step)
            if any(output.exists() for output in step_outputs):
                status['steps_completed'].append(step)
                completed_steps += 1
            else:
                # Check for failure indicators
                step_logs = species_dir / "logs" / f"{step}.log"
                if step_logs.exists():
                    try:
                        with open(step_logs, 'r') as f:
                            log_content = f.read()
                            if 'error' in log_content.lower() or 'failed' in log_content.lower():
                                status['steps_failed'].append(step)
                    except Exception:
                        pass

        status['progress_percentage'] = (completed_steps / total_steps) * 100.0
        
        # Count quantified samples
        quant_dir = species_dir / "quant"
        q_count = 0
        if quant_dir.exists():
            for sample_dir in quant_dir.glob("*"):
                if sample_dir.is_dir() and (sample_dir / "abundance.tsv").exists():
                    q_count += 1
        status['quantified'] = q_count

        # Mark as failed if any critical steps failed
        critical_steps = ['metadata', 'quant', 'merge']
        if any(step in status['steps_failed'] for step in critical_steps):
            status['failed'] = True

    return status


def assess_all_species_progress(config_dir: Path) -> Dict[str, Any]:
    """Assess progress across all species in the workflow.

    Args:
        config_dir: Directory containing species configuration files

    Returns:
        Comprehensive progress assessment mapping species names to status
    """
    results = {}
    if not config_dir.exists():
        return results

    # Find all amalgkit yaml configs
    config_files = list(config_dir.glob("amalgkit_*.yaml"))
    for cfg_file in config_files:
        # Extract species name from filename: amalgkit_{species}.yaml
        name = cfg_file.stem.replace("amalgkit_", "")
        if name == "test": # Skip generic test configs
            continue
            
        try:
            from metainformant.rna.engine.workflow import load_workflow_config
            config = load_workflow_config(cfg_file)
            status = analyze_species_status(config.work_dir)
            results[name] = status
        except Exception as e:
            logger.error(f"Failed to assess progress for {name}: {e}")
            results[name] = {"error": str(e), "overall_status": "error"}

    return results


def get_sample_status(config_path: Path, sample_id: str) -> Dict[str, Any]:
    """Get status of a specific sample.

    Args:
        config_path: Path to workflow configuration file
        sample_id: Sample identifier (e.g., SRR ID)

    Returns:
        Dictionary with sample status
    """
    from metainformant.rna.engine.workflow import load_workflow_config
    try:
        config = load_workflow_config(config_path)
        status = {
            "sample_id": sample_id,
            "quantified": False,
            "has_fastq": False,
            "has_sra": False,
            "status": "undownloaded"
        }
        
        # Check quant
        quant_file = config.work_dir / "quant" / sample_id / "abundance.tsv"
        if quant_file.exists():
            status["quantified"] = True
            status["status"] = "quantified"
            return status
            
        # Check fastq
        fastq_dir = config.work_dir / "fastq" / "getfastq" / sample_id
        if fastq_dir.exists() and any(fastq_dir.glob("*.fastq*")):
            status["has_fastq"] = True
            status["status"] = "has_fastq"
            return status
            
        # Check SRA
        sra_dir = config.work_dir / "fastq" / "getfastq"
        if (sra_dir / f"{sample_id}.sra").exists() or (sra_dir / f"{sample_id}.sra.part").exists():
            status["has_sra"] = True
            status["status"] = "downloaded"
            
        return status
    except Exception as e:
        logger.error(f"Failed to get sample status for {sample_id}: {e}")
        return {"error": str(e)}


def find_unquantified_samples(config_path: Path) -> List[str]:
    """Find samples that have not been quantified.

    Args:
        config_path: Path to workflow configuration file

    Returns:
        List of unquantified sample IDs
    """
    from metainformant.rna.engine.workflow import load_workflow_config
    try:
        config = load_workflow_config(config_path)
        metadata_file = config.work_dir / "metadata" / "metadata.tsv"
        if not metadata_file.exists():
             metadata_file = config.work_dir / "metadata" / "metadata_selected.tsv"
             
        sample_ids = []
        if metadata_file.exists():
            import csv
            with open(metadata_file, "r") as f:
                reader = csv.DictReader(f, delimiter="\t")
                sample_ids = [row["run"] for row in reader if "run" in row]
        
        unquantified = []
        quant_dir = config.work_dir / "quant"
        for sid in sample_ids:
            if not (quant_dir / sid / "abundance.tsv").exists():
                unquantified.append(sid)
                
        return unquantified
    except Exception as e:
        logger.error(f"Failed to find unquantified samples: {e}")
        return []


def check_active_downloads() -> set[str]:
    """Check for active amalgkit download processes.

    Returns:
        Set of sample IDs currently being downloaded
    """
    import subprocess
    active = set()
    try:
        # Check for fasterq-dump processes
        output = subprocess.check_output(["ps", "aux"], text=True)
        for line in output.splitlines():
             if "fasterq-dump" in line or "amalgkit getfastq" in line:
                  # Try to extract sample ID
                  import re
                  m = re.search(r"SRR\d+", line)
                  if m:
                       active.add(m.group(0))
    except Exception as e:
        logger.error(f"Error checking active downloads: {e}")
    return active


def initialize_progress_tracking(config_path: Path) -> Dict[str, Any]:
    """Initialize progress tracking for a workflow.

    Args:
        config_path: Path to workflow configuration file

    Returns:
        Dictionary with initialization result
    """
    from metainformant.rna.engine.workflow import load_workflow_config
    try:
        config = load_workflow_config(config_path)
        metadata_file = config.work_dir / "metadata" / "metadata.tsv"
        if not metadata_file.exists():
             metadata_file = config.work_dir / "metadata" / "metadata_selected.tsv"
             
        sample_ids = []
        if metadata_file.exists():
            import csv
            with open(metadata_file, "r") as f:
                reader = csv.DictReader(f, delimiter="\t")
                sample_ids = [row["run"] for row in reader if "run" in row]
        else:
            return {"success": False, "error": "Metadata not found"}
            
        return {
            "success": True,
            "total_samples": len(sample_ids),
            "sample_ids": sample_ids
        }
    except Exception as e:
        logger.error(f"Failed to initialize progress tracking: {e}")
        return {"success": False, "error": str(e)}


def _get_step_output_files(species_dir: Path, step: str) -> List[Path]:
    """Get expected output files for a workflow step.

    Args:
        species_dir: Species directory
        step: Workflow step name

    Returns:
        List of expected output file paths
    """
    step_outputs = {
        'metadata': [species_dir / "metadata" / "metadata.tsv", species_dir / "metadata.tsv"],
        'integrate': [species_dir / "metadata" / "metadata.tsv"],
        'config': [species_dir / "config" / "amalgkit_config.tsv"],
        'select': [species_dir / "metadata" / "metadata_selected.tsv"],
        'getfastq': [species_dir / "fastq" / "getfastq"],
        'quant': [species_dir / "quant"],
        'merge': [species_dir / "merge"],
        'cstmm': [species_dir / "cstmm"],
        'curate': [species_dir / "curate"],
        'csca': [species_dir / "csca"],
        'sanity': [species_dir / "sanity_check.txt"]
    }

    return step_outputs.get(step, [])


def get_sample_progress_report(config_path: Path) -> Dict[str, Any]:
    """Get detailed per-sample progress report for a workflow.
    
    Args:
        config_path: Path to workflow configuration file
        
    Returns:
        Dictionary with per-sample pipeline status
    """
    from metainformant.rna.engine.workflow import load_workflow_config
    from metainformant.rna.analysis.validation import validate_all_samples
    
    try:
        config = load_workflow_config(config_path)
        validation_result = validate_all_samples(config)
        
        return {
            'config_path': str(config_path),
            'work_dir': str(config.work_dir),
            'total_samples': validation_result.get('total_samples', 0),
            'summary': validation_result.get('summary', {}),
            'per_sample': validation_result.get('per_sample', {}),
            'missing_stages': validation_result.get('missing_stages', {})
        }
    except Exception as e:
        logger.error(f"Failed to generate progress report: {e}")
        return {
            'error': str(e),
            'config_path': str(config_path)
        }







