#!/usr/bin/env python3
"""Basic workflow orchestration example.

This example demonstrates METAINFORMANT's workflow utilities for config-based data processing and step-by-step execution.

Usage:
    python examples/core/example_workflow.py

Output:
    output/examples/core/workflow_results.json
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Dict
from metainformant.core import io, logging, workflow

def main():
    """Demonstrate basic workflow orchestration."""
    # Setup output directory
    output_dir = Path("output/examples/core")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== METAINFORMANT Workflow Orchestration Example ===")

    # 1. Create sample configuration
    print("\n1. Creating sample workflow configuration...")

    workflow_config = {
        "workflow": {
            "name": "dna_sequence_analysis_demo",
            "description": "Demonstration of basic workflow orchestration",
            "version": "1.0"
        },
        "steps": {
            "load_data": {
                "description": "Load sequence data from file",
                "input_file": "sample_sequences.fasta",
                "expected_count": 5
            },
            "validate_data": {
                "description": "Validate sequence integrity",
                "min_length": 10,
                "allowed_bases": ["A", "T", "C", "G"]
            },
            "analyze_sequences": {
                "description": "Calculate sequence statistics",
                "calculate_gc": True,
                "find_motifs": ["ATCG", "GATC"]
            },
            "generate_report": {
                "description": "Create analysis summary report",
                "output_format": "json",
                "include_plots": False
            }
        },
        "output": {
            "directory": str(output_dir / "workflow_output"),
            "intermediate_files": True,
            "final_report": "analysis_report.json"
        }
    }

    config_file = output_dir / "workflow_config.yaml"
    import yaml
    with open(config_file, "w") as f:
        yaml.dump(workflow_config, f)

    print(f"✓ Created workflow config: {config_file}")

    # 2. Validate configuration file
    print("\n2. Validating configuration file...")

    is_valid, validation_errors = workflow.validate_config_file(config_file)

    if is_valid:
        print("✓ Configuration file is valid")
    else:
        print(f"✗ Configuration validation failed: {validation_errors}")
        return

    # 3. Create sample config programmatically
    print("\n3. Creating sample config programmatically...")

    sample_config_file = output_dir / "sample_workflow_config.yaml"
    workflow.create_sample_config(sample_config_file, sample_type="basic")

    print(f"✓ Created sample config: {sample_config_file}")

    # 4. Simulate workflow execution
    print("\n4. Simulating workflow execution...")

    logger = logging.get_logger("workflow_demo")

    def simulate_workflow_step(step_name: str, step_config: Dict[str, Any]) -> Dict[str, Any]:
        """Simulate execution of a workflow step."""
        logger.info(f"Starting step: {step_name}")
        start_time = time.time()

        # Simulate processing time based on step
        processing_times = {
            "load_data": 0.2,
            "validate_data": 0.1,
            "analyze_sequences": 0.5,
            "generate_report": 0.3
        }

        time.sleep(processing_times.get(step_name, 0.1))

        # Generate mock results based on step
        if step_name == "load_data":
            result = {
                "sequences_loaded": step_config.get("expected_count", 5),
                "total_bases": 1250,
                "file_format": "FASTA"
            }
        elif step_name == "validate_data":
            result = {
                "sequences_validated": 5,
                "invalid_sequences": 0,
                "validation_checks": ["length", "bases", "format"]
            }
        elif step_name == "analyze_sequences":
            result = {
                "average_length": 250,
                "gc_content": 0.42,
                "motifs_found": {"ATCG": 3, "GATC": 2},
                "complexity_score": 0.85
            }
        elif step_name == "generate_report":
            result = {
                "report_sections": ["summary", "statistics", "motifs"],
                "output_format": step_config.get("output_format", "json"),
                "file_size_kb": 15.7
            }
        else:
            result = {"status": "unknown_step"}

        elapsed = time.time() - start_time
        result["execution_time"] = elapsed
        result["status"] = "success"

        logger.info(f"Step {step_name} completed in {elapsed:.2f} seconds")
        return result

    # Execute workflow steps
    step_results = {}
    total_workflow_time = 0

    for step_name, step_config in workflow_config["steps"].items():
        step_result = simulate_workflow_step(step_name, step_config)
        step_results[step_name] = step_result
        total_workflow_time += step_result["execution_time"]

    logger.info(".2f")

    # 5. Download and process data example
    print("\n5. Demonstrating download and process pattern...")

    def mock_data_processor(data: Any) -> Dict[str, Any]:
        """Mock data processor for demonstration."""
        # Simulate processing downloaded data
        return {
            "processed_items": len(data) if isinstance(data, list) else 1,
            "data_type": type(data).__name__,
            "processing_timestamp": "2024-12-26T10:00:00Z"
        }

    # Example URLs (these won't actually be downloaded in this demo)
    example_urls = [
        "https://example.com/sample_data.json",
        "https://example.com/sequences.fasta"
    ]

    download_results = {}
    for url in example_urls:
        # In real usage, this would download actual data
        # Here we simulate with mock data
        mock_data = {"sample": "data", "url": url}
        try:
            # This would normally download and process real data
            result = workflow.download_and_process_data(url, mock_data_processor, output_dir)
            download_results[url] = result
        except Exception as e:
            download_results[url] = {"error": str(e), "status": "simulated"}

    print("✓ Demonstrated download and process pattern")

    # 6. Run config-based workflow
    print("\n6. Running config-based workflow...")

    try:
        # This demonstrates the full workflow execution pattern
        workflow_result = workflow.run_config_based_workflow(
            config_file,
            custom_param="demo_value"
        )
        print("✓ Config-based workflow completed")
        config_workflow_result = workflow_result
    except Exception as e:
        print(f"Config-based workflow demo: {e}")
        config_workflow_result = {"error": str(e), "status": "demo_only"}

    # 7. Create comprehensive results summary
    print("\n7. Creating comprehensive results summary...")

    summary = {
        "workflow_orchestration_demo": {
            "timestamp": "2024-12-26T10:00:00Z",
            "configuration": {
                "config_file_created": str(config_file.relative_to(output_dir)),
                "sample_config_created": str(sample_config_file.relative_to(output_dir)),
                "config_validation_passed": is_valid
            },
            "workflow_execution": {
                "steps_executed": len(step_results),
                "total_execution_time": total_workflow_time,
                "step_results": step_results,
                "workflow_status": "completed"
            },
            "data_processing": {
                "download_operations": len(download_results),
                "download_results": download_results,
                "processing_pattern": "download_and_process_data"
            },
            "config_based_execution": {
                "attempted": True,
                "result": config_workflow_result
            },
            "features_demonstrated": [
                "Configuration file creation and validation",
                "Sample configuration generation",
                "Step-by-step workflow execution",
                "Progress tracking and timing",
                "Data download and processing patterns",
                "Config-based workflow orchestration",
                "Error handling and result aggregation"
            ],
            "workflow_patterns": {
                "step_execution": "Simulated realistic bioinformatics workflow steps",
                "data_flow": "Config -> Validation -> Processing -> Results",
                "error_handling": "Graceful handling of validation and execution errors",
                "result_aggregation": "Comprehensive summary of all workflow outputs"
            }
        }
    }

    results_file = output_dir / "workflow_results.json"
    io.dump_json(summary, results_file, indent=2)

    print(f"✓ Comprehensive results saved to: {results_file}")

    print("\n=== Workflow Orchestration Example Complete ===")
    print("This example demonstrated METAINFORMANT's workflow orchestration capabilities:")
    print("- Configuration management and validation")
    print("- Step-by-step workflow execution")
    print("- Progress tracking and result aggregation")
    print("- Data download and processing patterns")
    print("- Config-based workflow orchestration")

    print(f"\nAll outputs saved to: {output_dir}")

if __name__ == "__main__":
    main()
