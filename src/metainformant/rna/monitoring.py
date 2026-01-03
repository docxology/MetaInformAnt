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


def check_workflow_progress(work_dir: Path, species: Optional[str] = None) -> Dict[str, Any]:
    """Check the progress of RNA-seq workflow.

    Args:
        work_dir: Workflow working directory
        species: Specific species to check (None for all)

    Returns:
        Dictionary with progress information
    """
    progress_info = {
        'species': {},
        'overall': {
            'total_species': 0,
            'completed_species': 0,
            'failed_species': 0
        }
    }

    # Find species directories or config files
    species_dirs = []
    if species:
        species_dir = work_dir / species
        if species_dir.exists():
            species_dirs = [species_dir]
    else:
        # Look for species subdirectories
        species_dirs = [d for d in work_dir.iterdir() if d.is_dir()]

    for species_dir in species_dirs:
        species_name = species_dir.name
        progress_info['species'][species_name] = analyze_species_status(species_dir)

    # Calculate overall statistics
    progress_info['overall']['total_species'] = len(progress_info['species'])
    progress_info['overall']['completed_species'] = sum(
        1 for status in progress_info['species'].values()
        if status.get('completed', False)
    )
    progress_info['overall']['failed_species'] = sum(
        1 for status in progress_info['species'].values()
        if status.get('failed', False)
    )

    return progress_info


def analyze_species_status(species_dir: Path) -> Dict[str, Any]:
    """Analyze the status of a specific species workflow.

    Args:
        species_dir: Species-specific working directory

    Returns:
        Dictionary with species status information
    """
    status = {
        'species': species_dir.name,
        'completed': False,
        'failed': False,
        'steps_completed': [],
        'steps_failed': [],
        'last_updated': None,
        'progress_percentage': 0.0
    }

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

        # Mark as failed if any critical steps failed
        critical_steps = ['metadata', 'quant', 'merge']
        if any(step in status['steps_failed'] for step in critical_steps):
            status['failed'] = True

    return status


def assess_all_species_progress(work_dir: Path) -> Dict[str, Any]:
    """Assess progress across all species in the workflow.

    Args:
        work_dir: Main workflow working directory

    Returns:
        Comprehensive progress assessment
    """
    assessment = {
        'timestamp': None,  # Would be set to current time
        'overall_status': 'unknown',
        'species_summary': {},
        'resource_usage': {},
        'recommendations': []
    }

    progress = check_workflow_progress(work_dir)

    # Determine overall status
    if progress['overall']['total_species'] == 0:
        assessment['overall_status'] = 'no_species'
    elif progress['overall']['failed_species'] > 0:
        assessment['overall_status'] = 'has_failures'
    elif progress['overall']['completed_species'] == progress['overall']['total_species']:
        assessment['overall_status'] = 'completed'
    else:
        assessment['overall_status'] = 'in_progress'

    assessment['species_summary'] = progress

    # Generate recommendations
    if progress['overall']['failed_species'] > 0:
        assessment['recommendations'].append(
            f"Check logs for {progress['overall']['failed_species']} failed species"
        )

    incomplete_species = progress['overall']['total_species'] - progress['overall']['completed_species']
    if incomplete_species > 0:
        assessment['recommendations'].append(
            f"Resume workflow for {incomplete_species} incomplete species"
        )

    return assessment


def initialize_progress_tracking(work_dir: Path) -> bool:
    """Initialize progress tracking for a workflow.

    Args:
        work_dir: Workflow working directory

    Returns:
        True if initialization successful, False otherwise
    """
    try:
        # Create progress tracking directory
        progress_dir = work_dir / ".progress"
        progress_dir.mkdir(parents=True, exist_ok=True)

        # Initialize progress file
        progress_file = progress_dir / "workflow_progress.json"
        initial_progress = {
            'initialized': True,
            'start_time': None,  # Would be set to current timestamp
            'species': {},
            'overall': {
                'total_steps': 11,  # Standard amalgkit steps
                'completed_steps': 0
            }
        }

        io.dump_json(initial_progress, progress_file)
        logger.info(f"Initialized progress tracking in {progress_dir}")

        return True

    except Exception as e:
        logger.error(f"Failed to initialize progress tracking: {e}")
        return False


def _get_step_output_files(species_dir: Path, step: str) -> List[Path]:
    """Get expected output files for a workflow step.

    Args:
        species_dir: Species directory
        step: Workflow step name

    Returns:
        List of expected output file paths
    """
    step_outputs = {
        'metadata': [species_dir / "metadata.tsv"],
        'integrate': [species_dir / "integrated_metadata.tsv"],
        'config': [species_dir / "amalgkit_config.tsv"],
        'select': [species_dir / "selected_samples.tsv"],
        'getfastq': [species_dir / "fastq" / ".complete"],  # Completion marker
        'quant': [species_dir / "quant" / ".complete"],
        'merge': [species_dir / "expression_matrix.tsv"],
        'cstmm': [species_dir / "expression_matrix_cstmm.tsv"],
        'curate': [species_dir / "expression_matrix_curated.tsv"],
        'csca': [species_dir / "csca_results.tsv"],
        'sanity': [species_dir / "sanity_check.txt"]
    }

    return step_outputs.get(step, [])







