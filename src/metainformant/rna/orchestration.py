"""RNA-seq workflow orchestration utilities.

This module provides high-level orchestration functions for managing
complex RNA-seq workflows across multiple species and samples.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

from metainformant.core import logging

logger = logging.get_logger(__name__)


def run_workflow_for_species(config_path: str | Path,
                           species: Optional[str] = None,
                           steps: Optional[List[str]] = None,
                           **kwargs: Any) -> Dict[str, Any]:
    """Run RNA-seq workflow for a specific species.

    Args:
        config_path: Path to workflow configuration file
        species: Species name (if None, inferred from config)
        steps: Specific steps to run (if None, run all)
        **kwargs: Additional workflow options

    Returns:
        Workflow execution results
    """
    from metainformant.rna.workflow import load_workflow_config, execute_workflow

    logger.info(f"Running workflow for species: {species}")

    # Load configuration
    config = load_workflow_config(config_path)

    # Override species if specified
    if species:
        config.species_list = [species]

    # Execute workflow
    return_codes = execute_workflow(config, **kwargs)

    results = {
        'species': species or config.species_list[0],
        'config_path': str(config_path),
        'return_codes': return_codes,
        'success': all(rc == 0 for rc in return_codes),
        'steps_executed': len(return_codes)
    }

    logger.info(f"Workflow completed for {results['species']}: {'SUCCESS' if results['success'] else 'FAILED'}")

    return results


def cleanup_unquantified_samples(config_path: str | Path) -> Dict[str, Any]:
    """Clean up samples that failed quantification.

    Args:
        config_path: Path to workflow configuration file

    Returns:
        Cleanup operation results
    """
    from metainformant.rna.workflow import load_workflow_config
    from metainformant.rna.cleanup import cleanup_unquantified_samples as cleanup_func

    config = load_workflow_config(config_path)

    logger.info(f"Cleaning up unquantified samples in {config.work_dir}")

    cleaned_samples = cleanup_func(config.work_dir)

    results = {
        'work_dir': str(config.work_dir),
        'cleaned_samples': cleaned_samples,
        'count': len(cleaned_samples)
    }

    logger.info(f"Cleaned up {len(cleaned_samples)} unquantified samples")

    return results


def monitor_workflows(work_dir: Path) -> Dict[str, Any]:
    """Monitor the status of all workflows in a directory.

    Args:
        work_dir: Main workflow directory

    Returns:
        Monitoring results for all workflows
    """
    from metainformant.rna.monitoring import assess_all_species_progress

    logger.info(f"Monitoring workflows in {work_dir}")

    assessment = assess_all_species_progress(work_dir)

    # Add summary logging
    overall = assessment.get('overall_status', 'unknown')
    total = assessment.get('species_summary', {}).get('overall', {}).get('total_species', 0)
    completed = assessment.get('species_summary', {}).get('overall', {}).get('completed_species', 0)

    logger.info(f"Workflow status: {overall} ({completed}/{total} species completed)")

    return assessment


def discover_species_configs(config_dir: Path) -> List[Dict[str, Any]]:
    """Discover available species configuration files.

    Args:
        config_dir: Directory containing configuration files

    Returns:
        List of species configuration information
    """
    configs = []

    if not config_dir.exists():
        return configs

    # Look for YAML files with amalgkit prefix
    for config_file in config_dir.glob("amalgkit_*.yaml"):
        try:
            species_name = config_file.stem.replace("amalgkit_", "")

            # Skip template and test configs
            if "template" in species_name.lower() or "test" in species_name.lower():
                continue

            config_info = {
                'species': species_name,
                'config_file': str(config_file),
                'path': config_file
            }

            configs.append(config_info)

        except Exception as e:
            logger.warning(f"Error processing config {config_file}: {e}")

    logger.info(f"Discovered {len(configs)} species configuration files")

    return configs



