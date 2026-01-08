"""Config-driven workflow execution and orchestration for METAINFORMANT.

This module provides utilities for running bioinformatics workflows based on
configuration files, with step validation, error recovery, and progress tracking.
"""

from __future__ import annotations

import importlib
import inspect
import time
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Union

from metainformant.core import config, io, paths
from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


class WorkflowStep:
    """Represents a single step in a workflow."""

    def __init__(
        self,
        name: str,
        function: Callable,
        config: Dict[str, Any],
        depends_on: Optional[List[str]] = None,
    ):
        """Initialize workflow step.

        Args:
            name: Step name
            function: Function to execute
            config: Step configuration
            depends_on: Names of steps this step depends on
        """
        self.name = name
        self.function = function
        self.config = config
        self.depends_on = depends_on or []
        self.status = "pending"  # pending, running, completed, failed
        self.result = None
        self.error = None
        self.start_time = None
        self.end_time = None

    def execute(self, **kwargs: Any) -> Any:
        """Execute the workflow step.

        Args:
            **kwargs: Additional arguments to pass to function

        Returns:
            Function result

        Raises:
            Exception: If step execution fails
        """
        self.status = "running"
        self.start_time = time.time()

        try:
            logger.info(f"Executing workflow step: {self.name}")

            # Merge step config with kwargs
            step_kwargs = {**self.config, **kwargs}

            # Execute function
            result = self.function(**step_kwargs)

            self.status = "completed"
            self.result = result
            self.end_time = time.time()

            logger.info(f"Completed workflow step: {self.name}")
            return result

        except Exception as e:
            self.status = "failed"
            self.error = str(e)
            self.end_time = time.time()

            logger.error(f"Workflow step failed: {self.name} - {e}")
            raise

    def duration(self) -> Optional[float]:
        """Get step execution duration in seconds."""
        if self.start_time and self.end_time:
            return self.end_time - self.start_time
        return None


class BaseWorkflowOrchestrator:
    """Base class for workflow orchestration.

    This class provides the foundation for running bioinformatics workflows
    with dependency management, error handling, and progress tracking.
    """

    def __init__(
        self,
        config: Dict[str, Any],
        working_dir: Optional[Union[str, Path]] = None,
    ):
        """Initialize workflow orchestrator.

        Args:
            config: Workflow configuration
            working_dir: Working directory for outputs
        """
        self.config = config
        self.working_dir = Path(working_dir) if working_dir else Path("output/workflow")
        self.working_dir.mkdir(parents=True, exist_ok=True)

        self.steps: Dict[str, WorkflowStep] = {}
        self.execution_order: List[str] = []
        self.results: Dict[str, Any] = {}

        # Initialize from config
        self._load_steps_from_config()

    def _load_steps_from_config(self) -> None:
        """Load workflow steps from configuration."""
        workflow_steps = self.config.get("steps", {})

        for step_name, step_config in workflow_steps.items():
            # Get function reference
            function_ref = step_config.get("function")
            if isinstance(function_ref, str):
                function = self._resolve_function(function_ref)
            else:
                function = function_ref

            if not callable(function):
                raise ValueError(f"Invalid function for step {step_name}: {function_ref}")

            # Create step
            depends_on = step_config.get("depends_on", [])
            step = WorkflowStep(
                name=step_name,
                function=function,
                config=step_config.get("config", {}),
                depends_on=depends_on,
            )

            self.steps[step_name] = step

        # Determine execution order based on dependencies
        self._resolve_execution_order()

    def _resolve_function(self, function_ref: str) -> Callable:
        """Resolve function reference from string.

        Args:
            function_ref: Function reference in format "module:function"

        Returns:
            Callable function

        Raises:
            ValueError: If function cannot be resolved
        """
        try:
            module_name, func_name = function_ref.split(":", 1)
            module = importlib.import_module(module_name)
            function = getattr(module, func_name)

            if not callable(function):
                raise ValueError(f"{function_ref} is not callable")

            return function

        except (ImportError, AttributeError, ValueError) as e:
            raise ValueError(f"Cannot resolve function {function_ref}: {e}") from e

    def _resolve_execution_order(self) -> None:
        """Resolve execution order based on dependencies."""
        # Simple topological sort
        visited = set()
        temp_visited = set()
        order = []

        def visit(step_name: str) -> None:
            if step_name in temp_visited:
                raise ValueError(f"Circular dependency detected involving {step_name}")
            if step_name in visited:
                return

            temp_visited.add(step_name)

            # Visit dependencies first
            step = self.steps[step_name]
            for dep in step.depends_on:
                if dep not in self.steps:
                    raise ValueError(f"Step {step_name} depends on unknown step {dep}")
                visit(dep)

            temp_visited.remove(step_name)
            visited.add(step_name)
            order.append(step_name)

        # Visit all steps
        for step_name in self.steps:
            if step_name not in visited:
                visit(step_name)

        self.execution_order = order

    def execute_step(self, step_name: str, **kwargs: Any) -> Any:
        """Execute a single workflow step.

        Args:
            step_name: Name of step to execute
            **kwargs: Additional arguments

        Returns:
            Step result

        Raises:
            ValueError: If step dependencies not satisfied
        """
        if step_name not in self.steps:
            raise ValueError(f"Unknown step: {step_name}")

        step = self.steps[step_name]

        # Check dependencies
        for dep in step.depends_on:
            if dep not in self.results:
                raise ValueError(f"Step {step_name} dependency {dep} not satisfied")

        # Execute step
        result = step.execute(**kwargs)
        self.results[step_name] = result

        return result

    def run_workflow(self) -> Dict[str, Any]:
        """Run complete workflow.

        Returns:
            Dictionary with workflow results and metadata

        Raises:
            RuntimeError: If workflow execution fails
        """
        logger.info(f"Starting workflow execution with {len(self.steps)} steps")

        start_time = time.time()
        failed_steps = []

        try:
            # Execute steps in order
            for step_name in self.execution_order:
                try:
                    self.execute_step(step_name)
                except Exception as e:
                    logger.error(f"Step {step_name} failed: {e}")
                    failed_steps.append(step_name)

                    # Check if we should continue or stop
                    if self.config.get("stop_on_failure", True):
                        raise RuntimeError(f"Workflow failed at step {step_name}: {e}") from e

            execution_time = time.time() - start_time

            # Collect results
            workflow_results = {
                "success": len(failed_steps) == 0,
                "total_steps": len(self.steps),
                "completed_steps": len(self.results),
                "failed_steps": failed_steps,
                "execution_time": execution_time,
                "results": self.results.copy(),
                "step_details": {
                    name: {
                        "status": step.status,
                        "duration": step.duration(),
                        "error": step.error,
                    }
                    for name, step in self.steps.items()
                },
            }

            if workflow_results["success"]:
                logger.info(f"Workflow completed successfully in {execution_time:.2f}s")
            else:
                logger.warning(f"Workflow completed with {len(failed_steps)} failed steps")

            return workflow_results

        except Exception as e:
            execution_time = time.time() - start_time
            logger.error(f"Workflow execution failed after {execution_time:.2f}s: {e}")

            return {
                "success": False,
                "error": str(e),
                "execution_time": execution_time,
                "completed_steps": len(self.results),
                "failed_steps": failed_steps + ["workflow_orchestration"],
                "results": self.results.copy(),
            }

    def validate_workflow_config(self) -> Tuple[bool, List[str]]:
        """Validate workflow configuration.

        Returns:
            Tuple of (is_valid, error_messages)
        """
        errors = []

        # Check required fields
        if "steps" not in self.config:
            errors.append("Workflow config missing 'steps' section")
            return False, errors

        # Validate each step
        for step_name, step_config in self.config["steps"].items():
            # Check function
            if "function" not in step_config:
                errors.append(f"Step {step_name} missing 'function' field")
                continue

            # Try to resolve function
            try:
                function_ref = step_config["function"]
                if isinstance(function_ref, str):
                    self._resolve_function(function_ref)
            except ValueError as e:
                errors.append(f"Step {step_name} has invalid function: {e}")

            # Check dependencies
            depends_on = step_config.get("depends_on", [])
            for dep in depends_on:
                if dep not in self.config["steps"]:
                    errors.append(f"Step {step_name} depends on unknown step {dep}")

        # Check for circular dependencies
        try:
            self._resolve_execution_order()
        except ValueError as e:
            errors.append(f"Dependency error: {e}")

        return len(errors) == 0, errors


def download_and_process_data(
    url: str,
    processor: Callable,
    output_dir: Union[str, Path],
    **kwargs: Any,
) -> Any:
    """Download data and process it with a given function.

    Args:
        url: URL to download data from
        processor: Function to process downloaded data
        output_dir: Directory to save outputs
        **kwargs: Additional arguments for processor

    Returns:
        Result of processing function
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Downloading data from {url}")

    # Download data
    data = io.download_json(url)

    # Process data
    logger.info("Processing downloaded data")
    result = processor(data, output_dir=output_dir, **kwargs)

    logger.info("Data processing completed")
    return result


def validate_config_file(
    config_path: Union[str, Path],
    schema_path: Optional[Union[str, Path]] = None,
) -> Tuple[bool, List[str]]:
    """Validate configuration file.

    Args:
        config_path: Path to configuration file
        schema_path: Optional path to JSON schema file

    Returns:
        Tuple of (is_valid, error_messages)
    """
    errors = []

    try:
        # Load config
        raw_config = config.load_mapping_from_file(config_path)

        # Basic structure validation
        if not isinstance(raw_config, dict):
            errors.append("Configuration must be a dictionary")
            return False, errors

        # Check for required workflow fields
        if "steps" in raw_config:
            # This is a workflow config
            if not isinstance(raw_config["steps"], dict):
                errors.append("'steps' must be a dictionary")

            for step_name, step_config in raw_config["steps"].items():
                if not isinstance(step_config, dict):
                    errors.append(f"Step '{step_name}' must be a dictionary")
                    continue

                if "function" not in step_config:
                    errors.append(f"Step '{step_name}' missing 'function' field")

        # Schema validation if provided
        if schema_path:
            try:
                from metainformant.core.utils.errors import validate_json_schema
                validate_json_schema(raw_config, schema_path)
            except Exception as e:
                errors.append(f"Schema validation failed: {e}")

    except Exception as e:
        errors.append(f"Configuration loading failed: {e}")
        return False, errors

    return len(errors) == 0, errors


def create_sample_config(
    output_path: Union[str, Path],
    sample_type: str = "basic",
) -> None:
    """Create a sample workflow configuration file.

    Args:
        output_path: Path to save sample config
        sample_type: Type of sample config to create
    """
    if sample_type == "basic":
        sample_config = {
            "name": "Sample Bioinformatics Workflow",
            "description": "Basic workflow template",
            "steps": {
                "data_download": {
                    "function": "metainformant.core.workflow:download_and_process_data",
                    "config": {
                        "url": "https://example.com/data.json",
                        "processor": "my_module.process_data",
                    },
                },
                "analysis": {
                    "function": "my_module.analyze_data",
                    "depends_on": ["data_download"],
                    "config": {
                        "method": "default",
                        "output_format": "json",
                    },
                },
            },
        }
    else:
        sample_config = {
            "name": "Advanced Bioinformatics Workflow",
            "description": "Advanced workflow with multiple analysis steps",
            "steps": {},
        }

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    io.dump_json(sample_config, output_path)
    logger.info(f"Created sample workflow config at {output_path}")


def run_config_based_workflow(
    config_path: Union[str, Path],
    **kwargs: Any,
) -> Dict[str, Any]:
    """Run a workflow from configuration file.

    Args:
        config_path: Path to workflow configuration
        **kwargs: Additional workflow arguments

    Returns:
        Workflow execution results

    Raises:
        RuntimeError: If workflow fails
    """
    # Load and validate config
    is_valid, errors = validate_config_file(config_path)
    if not is_valid:
        raise RuntimeError(f"Invalid workflow config: {errors}")

    # Load config
    workflow_config = config.load_mapping_from_file(config_path)

    # Create orchestrator
    orchestrator = BaseWorkflowOrchestrator(workflow_config, **kwargs)

    # Run workflow
    results = orchestrator.run_workflow()

    if not results["success"]:
        raise RuntimeError(f"Workflow execution failed: {results.get('error', 'Unknown error')}")

    return results






