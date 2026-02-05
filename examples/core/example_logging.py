#!/usr/bin/env python3
"""Logging configuration and usage example.

This example demonstrates METAINFORMANT's structured logging system with console and file output, different log levels, and metadata logging.

Usage:
    python examples/core/example_logging.py

Output:
    output/examples/core/logging_example.log
"""

from __future__ import annotations

import time
from pathlib import Path

from metainformant.core import logging


def main():
    """Demonstrate logging functionality."""
    # Setup output directory
    output_dir = Path("output/examples/core")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== METAINFORMANT Logging Example ===")

    # 1. Basic logger setup
    print("\n1. Basic logger setup...")
    logger = logging.get_logger(__name__)
    logger.info("Basic logger created and ready")

    # 2. Setup logger with file output
    print("\n2. Setting up logger with file output...")
    log_file = output_dir / "logging_example.log"
    file_logger = logging.setup_logger("example_logger", log_file=str(log_file), level="DEBUG")

    print(f"✓ Logger configured with file output: {log_file}")

    # 3. Demonstrate different log levels
    print("\n3. Demonstrating different log levels...")

    file_logger.debug("This is a DEBUG message (detailed diagnostic info)")
    file_logger.info("This is an INFO message (general information)")
    file_logger.warning("This is a WARNING message (potential issue)")
    file_logger.error("This is an ERROR message (serious problem)")
    file_logger.critical("This is a CRITICAL message (system failure)")

    # 4. Simulate a workflow with logging
    print("\n4. Simulating workflow with logging...")

    def simulate_analysis_step(step_name: str, duration: float):
        """Simulate an analysis step with logging."""
        file_logger.info(f"Starting {step_name}")
        start_time = time.time()

        # Simulate work
        time.sleep(duration)

        elapsed = time.time() - start_time
        file_logger.info(f"Step completed in {elapsed:.2f} seconds")
        return elapsed

    # Run simulated workflow
    workflow_steps = [
        ("data_loading", 0.1),
        ("data_validation", 0.05),
        ("sequence_analysis", 0.2),
        ("statistics_calculation", 0.15),
        ("report_generation", 0.1),
    ]

    total_time = 0
    for step_name, duration in workflow_steps:
        step_time = simulate_analysis_step(step_name, duration)
        total_time += step_time

    file_logger.info(".2f")

    # 5. Metadata logging
    print("\n5. Demonstrating metadata logging...")

    # Log with structured metadata
    analysis_metadata = {
        "analysis_type": "dna_sequence_analysis",
        "sequences_processed": 150,
        "average_length": 450,
        "gc_content_range": [0.35, 0.65],
        "quality_threshold": 0.95,
        "timestamp": "2024-12-26T10:00:00Z",
    }

    logging.log_with_metadata(
        file_logger, "Analysis completed successfully", analysis_metadata, level="INFO", structured=True
    )

    # 6. Environment-based configuration
    print("\n6. Environment-based configuration...")
    print("You can control logging level with environment variables:")
    print("  export CORE_LOG_LEVEL=DEBUG  # For detailed output")
    print("  export CORE_LOG_LEVEL=INFO   # For normal output")
    print("  export CORE_LOG_LEVEL=WARNING # For quiet output")

    # Configure logging from environment (will use defaults if not set)
    logging.configure_logging_from_env(default_level="INFO")

    # 7. Multiple loggers for different components
    print("\n7. Multiple loggers for different components...")

    data_logger = logging.get_logger("data_processing")
    analysis_logger = logging.get_logger("sequence_analysis")
    report_logger = logging.get_logger("report_generation")

    data_logger.info("Data preprocessing completed")
    analysis_logger.info("Sequence alignment finished")
    report_logger.info("HTML report generated")

    # 8. Performance logging
    print("\n8. Performance logging example...")

    def time_operation(operation_name: str, operation_func):
        """Time an operation and log the result."""
        start_time = time.time()
        result = operation_func()
        elapsed = time.time() - start_time

        file_logger.info(f"Operation completed in {elapsed:.3f} seconds")
        return result, elapsed

    # Example operations
    def load_sequences():
        return ["ATCG" * 25] * 100  # Simulate loading 100 sequences

    def calculate_gc_content(sequences):
        return sum(seq.count("G") + seq.count("C") for seq in sequences) / sum(len(seq) for seq in sequences)

    def find_motifs(sequences, motif="ATCG"):
        return sum(seq.count(motif) for seq in sequences)

    # Time operations
    sequences, load_time = time_operation("sequence loading", load_sequences)
    gc_content, gc_time = time_operation("GC content calculation", lambda: calculate_gc_content(sequences))
    motif_count, motif_time = time_operation("motif finding", lambda: find_motifs(sequences))

    file_logger.info(".2f")
    file_logger.info(f"Found {motif_count} motif occurrences")

    # 9. Create summary
    print("\n9. Creating logging summary...")

    # Read the log file to show what was captured
    with open(log_file, "r") as f:
        log_lines = f.readlines()

    summary = {
        "logging_demo": {
            "timestamp": "2024-12-26T10:00:00Z",
            "log_file": str(log_file.relative_to(output_dir)),
            "total_log_lines": len(log_lines),
            "log_levels_used": ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
            "features_demonstrated": [
                "File and console logging",
                "Multiple log levels",
                "Structured metadata logging",
                "Environment-based configuration",
                "Multiple component loggers",
                "Performance timing",
                "Workflow step tracking",
            ],
            "workflow_simulation": {
                "steps_completed": len(workflow_steps),
                "total_workflow_time": total_time,
                "performance_operations": {
                    "sequence_loading": load_time,
                    "gc_content_calculation": gc_time,
                    "motif_finding": motif_time,
                },
            },
            "sample_log_entries": log_lines[-5:] if len(log_lines) >= 5 else log_lines,
        }
    }

    summary_file = output_dir / "logging_demo_summary.json"
    from metainformant.core import io

    io.dump_json(summary, summary_file, indent=2)

    print(f"✓ Logging summary saved to: {summary_file}")
    print(f"✓ Complete log file: {log_file}")
    print(f"✓ Total log entries: {len(log_lines)}")

    print("\n=== Logging Example Complete ===")
    print("Check the log file to see structured logging output!")
    print("\nKey takeaways:")
    print("- Use setup_logger() for file output with proper formatting")
    print("- Use get_logger() for component-specific logging")
    print("- Use log_with_metadata() for structured data logging")
    print("- Configure logging level via CORE_LOG_LEVEL environment variable")


if __name__ == "__main__":
    main()
