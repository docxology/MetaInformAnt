#!/usr/bin/env python3
"""Demonstration script showing complete workflow with visualization.

This script demonstrates:
1. Configuration loading and saving
2. Data processing
3. Visualization generation with informative names
4. Output organization in output/ directory

Usage:
    python3 scripts/core/run_demo.py
"""

import logging
import sys
from pathlib import Path

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from metainformant.core import io
from metainformant.core.paths import ensure_directory
from metainformant.core.utils.logging import setup_logger


def main():
    """Run complete demonstration workflow."""
    # Setup logging
    logger = setup_logger(__name__, level="INFO")
    logger.info("Complete workflow demonstration starting")
    
    # 1. Define configuration
    config = {
        "workflow": "demo",
        "version": "1.0.0",
        "parameters": {
            "threshold": 0.05,
            "iterations": 100,
            "method": "example"
        },
        "input_source": "generated",
        "output_format": "json",
        "timestamp": "2025-11-03"
    }
    
    # 2. Setup output directory
    output_base = Path("output/demo")
    ensure_directory(output_base)
    logger.info(f"Output directory: {output_base}")
    
    # 3. Save configuration with informative name
    config_path = output_base / "workflow_configuration.json"
    io.dump_json(config, config_path)
    logger.info(f"✅ Configuration saved: {config_path}")
    
    # 4. Generate and save input data with informative name
    input_data = {
        "samples": ["sample_A", "sample_B", "sample_C"],
        "values": [1.5, 2.3, 1.8],
        "metadata": {"source": "simulation", "n": 3}
    }
    input_path = output_base / "input_samples.json"
    io.dump_json(input_data, input_path)
    logger.info(f"✅ Input data saved: {input_path}")
    
    # 5. Process data (example transformation)
    processed_data = {
        "samples": input_data["samples"],
        "normalized_values": [v / sum(input_data["values"]) for v in input_data["values"]],
        "statistics": {
            "mean": sum(input_data["values"]) / len(input_data["values"]),
            "min": min(input_data["values"]),
            "max": max(input_data["values"])
        }
    }
    processed_path = output_base / "processed_normalized_data.json"
    io.dump_json(processed_data, processed_path)
    logger.info(f"✅ Processed data saved: {processed_path}")
    
    # 6. Generate visualizations with informative names
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from metainformant.visualization import lineplot, heatmap
        import numpy as np
        
        # Create visualizations subdirectory
        viz_dir = output_base / "visualizations"
        ensure_directory(viz_dir)
        
        # Visualization 1: Line plot of input values
        fig1, ax1 = plt.subplots(figsize=(8, 6))
        ax1 = lineplot(None, input_data["values"], ax=ax1, style='-o', color='steelblue')
        ax1.set_xlabel("Sample Index")
        ax1.set_ylabel("Value")
        ax1.set_title("Input Sample Values")
        ax1.grid(True, alpha=0.3)
        fig1_path = viz_dir / "input_sample_values_lineplot.png"
        fig1.savefig(fig1_path, dpi=300, bbox_inches='tight')
        plt.close(fig1)
        logger.info(f"✅ Visualization saved: {fig1_path}")
        
        # Visualization 2: Bar plot of normalized values
        fig2, ax2 = plt.subplots(figsize=(8, 6))
        ax2.bar(processed_data["samples"], processed_data["normalized_values"], color='coral')
        ax2.set_xlabel("Sample")
        ax2.set_ylabel("Normalized Value")
        ax2.set_title("Normalized Sample Values (Proportions)")
        ax2.grid(True, alpha=0.3, axis='y')
        fig2_path = viz_dir / "normalized_values_barplot.png"
        fig2.savefig(fig2_path, dpi=300, bbox_inches='tight')
        plt.close(fig2)
        logger.info(f"✅ Visualization saved: {fig2_path}")
        
        # Visualization 3: Heatmap example
        matrix_data = np.array([[1.5, 2.3, 1.8], [2.1, 1.9, 2.5], [1.7, 2.2, 1.6]])
        fig3, ax3 = plt.subplots(figsize=(8, 6))
        ax3 = heatmap(matrix_data, cmap='viridis', cbar=True, ax=ax3, annot=True)
        ax3.set_title("Sample Correlation Heatmap")
        ax3.set_xlabel("Sample")
        ax3.set_ylabel("Sample")
        fig3_path = viz_dir / "sample_correlation_heatmap.png"
        fig3.savefig(fig3_path, dpi=300, bbox_inches='tight')
        plt.close(fig3)
        logger.info(f"✅ Visualization saved: {fig3_path}")
        
        # Save visualization metadata
        viz_metadata = {
            "visualizations_generated": [
                str(fig1_path.name),
                str(fig2_path.name),
                str(fig3_path.name)
            ],
            "format": "PNG",
            "dpi": 300,
            "timestamp": config["timestamp"]
        }
        viz_meta_path = viz_dir / "visualization_metadata.json"
        io.dump_json(viz_metadata, viz_meta_path)
        logger.info(f"✅ Visualization metadata saved: {viz_meta_path}")
        
    except Exception as e:
        logger.error(f"Visualization generation failed: {e}", exc_info=True)
        return 1
    
    # 7. Generate summary report
    summary = {
        "workflow_status": "completed",
        "files_generated": {
            "configuration": str(config_path),
            "input_data": str(input_path),
            "processed_data": str(processed_path),
            "visualizations": [
                str(fig1_path),
                str(fig2_path),
                str(fig3_path)
            ]
        },
        "statistics": processed_data["statistics"],
        "output_directory": str(output_base)
    }
    summary_path = output_base / "workflow_summary_report.json"
    io.dump_json(summary, summary_path)
    logger.info(f"✅ Summary report saved: {summary_path}")
    
    logger.info("=" * 60)
    logger.info("Workflow completed successfully!")
    logger.info(f"All outputs saved to: {output_base}")
    logger.info(f"  - Configuration: {config_path.name}")
    logger.info(f"  - Input data: {input_path.name}")
    logger.info(f"  - Processed data: {processed_path.name}")
    logger.info(f"  - Visualizations: {viz_dir.name}/ ({len(viz_metadata['visualizations_generated'])} files)")
    logger.info(f"  - Summary: {summary_path.name}")
    logger.info("=" * 60)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

