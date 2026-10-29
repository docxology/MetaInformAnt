#!/usr/bin/env python3
"""Run a real configuration-driven workflow over project metadata.

Usage:
    uv run python examples/core/example_workflow.py

Output:
    output/examples/core/workflow_config.json
    output/examples/core/workflow_results.json
"""

from __future__ import annotations

from pathlib import Path

from metainformant.core.execution import workflow
from metainformant.core.io import io


def main():
    """Read project metadata in two steps and join their actual results."""
    repository = Path(__file__).resolve().parents[2]
    output_dir = Path("output/examples/core")
    output_dir.mkdir(parents=True, exist_ok=True)

    workflow_config = {
        "steps": {
            "project": {
                "function": "metainformant.core.utils.config:load_mapping_from_file",
                "params": {"config_path": str(repository / "pyproject.toml")},
            },
            "examples": {
                "function": "metainformant.core.utils.config:load_mapping_from_file",
                "params": {"config_path": str(repository / "examples" / "dependencies.json")},
            },
            "summary": {
                "function": "builtins:dict",
                "depends_on": ["project", "examples"],
            },
        },
    }
    config_file = output_dir / "workflow_config.json"
    io.dump_json(workflow_config, config_file)

    is_valid, errors = workflow.validate_config_file(config_file)
    if not is_valid:
        raise ValueError(f"Invalid workflow configuration: {errors}")

    execution = workflow.run_config_based_workflow(config_file)
    if not execution["success"]:
        raise RuntimeError(f"Workflow execution failed: {execution}")

    joined = execution["results"]["summary"]
    project = joined["project"]["project"]
    observations = {
        "project_name": project["name"],
        "project_version": project["version"],
        "documented_examples": len(joined["examples"]["example_dependencies"]),
    }
    results_file = output_dir / "workflow_results.json"
    io.dump_json({"workflow_execution": execution, "observations": observations}, results_file)
    print(f"Completed {len(execution['execution_order'])} workflow steps")
    print(f"Read project {observations['project_name']} {observations['project_version']}")
    print(f"Found {observations['documented_examples']} documented examples")
    print(f"Results: {results_file}")


if __name__ == "__main__":
    main()
