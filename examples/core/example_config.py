#!/usr/bin/env python3
"""Configuration management example.

This example demonstrates how to load configuration files with environment variable overrides using METAINFORMANT's core configuration utilities.

Usage:
    python examples/core/example_config.py

Output:
    output/examples/core/config_example.json
"""

from __future__ import annotations

from pathlib import Path

from metainformant.core import io, paths
from metainformant.core.config import apply_env_overrides, load_mapping_from_file


def main():
    """Demonstrate configuration loading with environment overrides."""
    # Setup output directory
    output_dir = Path("output/examples/core")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create a sample config file for demonstration
    sample_config = {
        "threads": 4,
        "work_dir": "output/example_workflow",
        "log_dir": "output/example_workflow/logs",
        "max_memory_gb": 8,
        "debug": False,
        "analysis": {"method": "fast", "quality_threshold": 0.95, "save_intermediates": True},
    }

    # Save sample config
    config_file = output_dir / "sample_config.yaml"
    import yaml

    with open(config_file, "w") as f:
        yaml.dump(sample_config, f, default_flow_style=False)

    print("=== METAINFORMANT Configuration Example ===")
    print(f"Created sample config file: {config_file}")

    # Load configuration from file
    print("\n1. Loading configuration from file...")
    config = load_mapping_from_file(config_file)
    print(f"Loaded config: {config}")

    # Apply environment variable overrides
    print("\n2. Applying environment overrides...")
    print("You can override these settings using environment variables:")
    print("  export CORE_THREADS=8")
    print("  export CORE_WORK_DIR=/tmp/custom_work")
    print("  export CORE_LOG_DIR=/tmp/custom_logs")

    # Apply overrides (will use defaults if env vars not set)
    config_with_overrides = apply_env_overrides(config, prefix="CORE")
    print(f"Config with overrides: {config_with_overrides}")

    # Show what changed
    changes = {}
    for key, value in config_with_overrides.items():
        if key not in config or config[key] != value:
            changes[key] = {"from": config.get(key), "to": value}

    if changes:
        print("\nChanges from environment overrides:")
        for key, change in changes.items():
            print(f"  {key}: {change['from']} -> {change['to']}")
    else:
        print("\nNo environment overrides applied (set CORE_* variables to see changes)")

    # Demonstrate configuration validation
    print("\n3. Configuration validation...")

    # Check required fields
    required_fields = ["threads", "work_dir"]
    missing_fields = [field for field in required_fields if field not in config_with_overrides]

    if missing_fields:
        print(f"Missing required fields: {missing_fields}")
    else:
        print("All required fields present ✓")

    # Validate types
    validation_results = []

    if isinstance(config_with_overrides.get("threads"), int):
        validation_results.append("threads: valid integer ✓")
    else:
        validation_results.append("threads: invalid (should be integer) ✗")

    if isinstance(config_with_overrides.get("work_dir"), str):
        validation_results.append("work_dir: valid string ✓")
    else:
        validation_results.append("work_dir: invalid (should be string) ✗")

    print("Type validation:")
    for result in validation_results:
        print(f"  {result}")

    # Save final configuration
    result_file = output_dir / "config_example.json"
    result_data = {
        "original_config": config,
        "config_with_overrides": config_with_overrides,
        "changes_applied": changes,
        "validation_results": validation_results,
        "environment_variables_used": {
            "CORE_THREADS": "int (number of threads)",
            "CORE_WORK_DIR": "str (working directory path)",
            "CORE_LOG_DIR": "str (logging directory path)",
        },
    }

    io.dump_json(result_data, result_file)

    print("4. Results saved to:")
    print(f"   {result_file}")

    print("\n=== Configuration Example Complete ===")
    print("Try setting environment variables and re-running to see overrides in action!")


if __name__ == "__main__":
    main()
