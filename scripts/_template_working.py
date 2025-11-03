#!/usr/bin/env python3
"""Comprehensive script template with working examples.

This demonstrates proper imports, configuration, visualization, and output handling.
"""

import logging
import sys
from pathlib import Path

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core import io, paths
from metainformant.core.logging import setup_logger

# Example of how scripts should be structured
def main():
    # Setup logging
    logger = setup_logger(__name__, level="INFO")
    logger.info("Script starting")
    
    # Ensure output directory
    output_dir = paths.ensure_output_dir(Path("output/example"))
    logger.info(f"Output directory: {output_dir}")
    
    # Save configuration
    config = {
        "input": "data/example/input.txt",
        "parameters": {"threshold": 0.05},
        "timestamp": "2025-11-03"
    }
    config_path = output_dir / "config.json"
    io.write_json(config, config_path)
    logger.info(f"Configuration saved: {config_path}")
    
    # Generate visualization (if visualization module used)
    try:
        import matplotlib.pyplot as plt
        from metainformant.visualization import lineplot
        
        # Example data
        data = [1, 4, 2, 8, 5, 7]
        ax = lineplot(None, data)
        ax.set_xlabel("Time")
        ax.set_ylabel("Value")
        ax.set_title("Example Time Series")
        
        # Save figure with informative name
        fig_path = output_dir / "time_series_visualization.png"
        ax.figure.savefig(fig_path, dpi=300, bbox_inches='tight')
        plt.close(ax.figure)
        logger.info(f"Visualization saved: {fig_path}")
    except Exception as e:
        logger.warning(f"Visualization skipped: {e}")
    
    logger.info("Script complete")
    return 0

if __name__ == "__main__":
    sys.exit(main())



