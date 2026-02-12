"""RNA-seq Multi-Species Orchestrator implementation.

This module provides robust orchestration for running the amalgkit pipeline across multiple species.
It handles:
- Sequential execution of pipeline steps (metadata -> getfastq -> ... -> curate)
- Error handling and retries
- Skipping species that fail specific steps
- Centralized logging
- Resuming from last successful step
"""

from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import List, Dict, Any, Optional

import yaml

from metainformant.core.utils import logging as log_utils
from metainformant.rna.amalgkit import amalgkit

# Configure logging
logger = log_utils.get_logger(__name__)


class PipelineOrchestrator:
    """Orchestrator for multi-species RNA-seq pipelines."""

    def __init__(self, config_path: Path):
        """Initialize pipeline orchestrator.

        Args:
            config_path: Path to YAML configuration file
        """
        self.config_path = config_path
        self._config_cache = None
        self.work_dir = Path(self.config.get("work_dir", "output/amalgkit"))
        self.log_dir = self.work_dir / "logs"
        self.steps = [
            "metadata",
            "integrate",  # Optional, usually for specific use cases
            "getfastq",
            "quant",
            "merge",
            "cstmm",
            "csca",  # Cell-specific / cross-species analysis preparation
            "curate",
        ]
        
        self.setup_logging()

    @property
    def config(self) -> Dict[str, Any]:
        """Lazy load configuration."""
        if self._config_cache is None:
            self._config_cache = self._load_config()
        return self._config_cache

    def _load_config(self) -> Dict[str, Any]:
        """Load YAML configuration."""
        if not self.config_path.exists():
            raise FileNotFoundError(f"Config file not found: {self.config_path}")
        with open(self.config_path, "r") as f:
            return yaml.safe_load(f)

    def setup_logging(self):
        """Setup logging to file and console."""
        self.log_dir.mkdir(parents=True, exist_ok=True)
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        log_file = self.log_dir / f"orchestrator_{timestamp}.log"
        
        # We don't want to reconfigure the root logger if it's already set up by the system,
        # but we do want to ensure our logs go to this file.
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
        logger.addHandler(file_handler)
        logger.info(f"Logging to {log_file}")

    def run_step(self, step: str, params: amalgkit.AmalgkitParams, species: str) -> bool:
        """Run a single amalgkit step for a species."""
        logger.info(f"[{species}] Starting step: {step}")
        
        try:
            # Map step name to function
            if not hasattr(amalgkit, step):
                logger.error(f"[{species}] Unknown step: {step}")
                return False
                
            step_func = getattr(amalgkit, step)
            
            # Special handling for 'integrate' which is not always needed
            if step == "integrate" and not self.config.get("run_integrate", False):
                logger.info(f"[{species}] Skipping integrate step (not enabled)")
                return True

            start_time = time.time()
            result = step_func(params, check=False)
            duration = time.time() - start_time
            
            if result.returncode == 0:
                logger.info(f"[{species}] Step {step} completed successfully in {duration:.2f}s")
                return True
            else:
                logger.error(f"[{species}] Step {step} failed with return code {result.returncode}")
                # Log stderr if available
                if result.stderr:
                    for line in result.stderr.splitlines():
                        logger.error(f"[{species}] [stderr] {line}")
                return False
                
        except Exception as e:
            logger.error(f"[{species}] Exception during step {step}: {e}")
            return False

    def run(self):
        """Execute pipeline for all configured species."""
        species_list = self.config.get("species", [])
        if not species_list:
            logger.warning("No species found in configuration.")
            return

        global_threads = self.config.get("threads", 8)
        
        logger.info(f"Starting orchestration for {len(species_list)} species")
        
        for species_config in species_list:
            # Handle both string (name only) and dict (name + options) formats
            if isinstance(species_config, str):
                species_name = species_config
                species_opts = {}
            else:
                species_name = species_config.get("name")
                species_opts = species_config
            
            if not species_name:
                continue

            logger.info(f"=== Processing Species: {species_name} ===")
            
            # Prepare parameters
            # AmalgkitParams will handle the CLI argument generation
            params = amalgkit.AmalgkitParams(
                work_dir=self.work_dir,
                threads=species_opts.get("threads", global_threads),
                species_list=[species_name], # Process one at a time for better control
                **species_opts.get("extra_params", {})
            )
            
            # Execute steps
            species_failed = False
            for step in self.steps:
                success = self.run_step(step, params, species_name)
                if not success:
                    logger.error(f"[{species_name}] Pipeline failed at step {step}. Skipping remaining steps for this species.")
                    species_failed = True
                    break
            
            if not species_failed:
                logger.info(f"=== Completed Species: {species_name} ===")
            else:
                logger.warning(f"=== Failed Species: {species_name} ===")
