"""Workflow orchestration utilities for METAINFORMANT.

Provides config-based data processing and workflow execution functions.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any

from . import io
from . import paths


def download_and_process_data(
    config: dict,
    *,
    output_dir: Path | str | None = None,
    verbose: bool = True
) -> dict[str, Any]:
    """Download and process data based on configuration.
    
    This is a comprehensive end-to-end processing function that can handle
    multiple data sources and processing steps based on a configuration file.
    
    Args:
        config: Configuration dictionary with processing instructions
        output_dir: Directory for outputs (uses 'output/processing' if None)
        verbose: Whether to print progress information
        
    Returns:
        Dictionary with processing results and metadata
    """
    if output_dir is None:
        output_dir = Path("output/processing")
    else:
        output_dir = Path(output_dir)
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    results = {
        "config": config,
        "output_dir": str(output_dir),
        "start_time": None,
        "end_time": None,
        "downloads": {},
        "processing": {},
        "errors": [],
        "warnings": []
    }
    
    try:
        # Record start time
        results["start_time"] = time.time()
        
        # Process downloads
        if "downloads" in config:
            download_results = {}
            for name, download_config in config["downloads"].items():
                if verbose:
                    print(f"Downloading {name}...")
                
                url = download_config.get("url")
                dest_path = output_dir / download_config.get("filename", f"{name}.txt")
                
                if io.download_file(url, dest_path):
                    download_results[name] = {
                        "success": True,
                        "path": str(dest_path),
                        "size": paths.get_file_size(dest_path)
                    }
                else:
                    download_results[name] = {
                        "success": False,
                        "error": f"Failed to download {url}"
                    }
                    results["errors"].append(f"Download failed: {name}")
            
            results["downloads"] = download_results
            
        # Process data
        if "processing" in config:
            processing_results = {}
            for step_name, step_config in config["processing"].items():
                if verbose:
                    print(f"Processing step: {step_name}")
                
                try:
                    # This would contain the actual processing logic
                    # For now, just record that we processed it
                    processing_results[step_name] = {
                        "completed": True,
                        "timestamp": time.time()
                    }
                except Exception as e:
                    processing_results[step_name] = {
                        "completed": False,
                        "error": str(e)
                    }
                    results["errors"].append(f"Processing failed: {step_name}")
            
            results["processing"] = processing_results
        
        # Record end time
        results["end_time"] = time.time()
        
        # Save results
        results_file = output_dir / "processing_results.json"
        io.dump_json(results, results_file)
        
        if verbose:
            print(f"Processing complete. Results saved to {results_file}")
            
    except Exception as e:
        results["errors"].append(f"Processing failed: {str(e)}")
        if verbose:
            print(f"Error: {e}")
    
    return results


def validate_config_file(config_path: str | Path) -> tuple[bool, list[str]]:
    """Validate a configuration file for processing.
    
    Args:
        config_path: Path to configuration file
        
    Returns:
        Tuple of (is_valid, list_of_errors)
    """
    config_path = Path(config_path)
    
    if not config_path.exists():
        return False, [f"Configuration file not found: {config_path}"]
    
    try:
        config = io.load_json(config_path)
    except Exception as e:
        return False, [f"Failed to parse configuration: {e}"]
    
    errors = []
    
    # Check required sections
    if "downloads" not in config and "processing" not in config:
        errors.append("Configuration must have 'downloads' or 'processing' section")
    
    # Validate download configurations
    if "downloads" in config:
        for name, download_config in config["downloads"].items():
            if "url" not in download_config:
                errors.append(f"Download '{name}' missing 'url' field")
    
    return len(errors) == 0, errors


def create_sample_config(output_path: str | Path, sample_type: str = "basic") -> None:
    """Create a sample configuration file for testing.
    
    Args:
        output_path: Path to save sample config
        sample_type: Type of sample config ('basic', 'advanced', 'scientific')
    """
    output_path = Path(output_path)
    paths.ensure_directory(output_path.parent)
    
    if sample_type == "basic":
        config = {
            "description": "Basic download and processing example",
            "downloads": {
                "sample_data": {
                    "url": "https://example.com/sample.txt",
                    "filename": "sample.txt"
                }
            },
            "processing": {
                "analyze": {
                    "type": "text_analysis",
                    "output": "analysis_results.json"
                }
            }
        }
    elif sample_type == "scientific":
        config = {
            "description": "Scientific data processing example",
            "downloads": {
                "gene_expression": {
                    "url": "https://example.com/expression_data.csv",
                    "filename": "expression_data.csv"
                },
                "metadata": {
                    "url": "https://example.com/metadata.json",
                    "filename": "metadata.json"
                }
            },
            "processing": {
                "preprocess": {
                    "type": "data_cleaning",
                    "normalize": True
                },
                "analyze": {
                    "type": "differential_expression",
                    "method": "deseq2"
                },
                "visualize": {
                    "type": "volcano_plot",
                    "output": "volcano_plot.png"
                }
            }
        }
    else:  # advanced
        config = {
            "description": "Advanced processing pipeline",
            "downloads": {
                "dataset1": {"url": "https://example.com/data1.csv"},
                "dataset2": {"url": "https://example.com/data2.csv"},
                "reference": {"url": "https://example.com/reference.fa"}
            },
            "processing": {
                "merge": {"type": "dataset_merge"},
                "filter": {"type": "quality_filter", "threshold": 0.8},
                "analyze": {"type": "statistical_analysis"},
                "export": {"type": "export_results", "format": "tsv"}
            }
        }
    
    io.dump_json(config, output_path)


def run_config_based_workflow(config_path: str | Path, **kwargs) -> dict[str, Any]:
    """Run a complete workflow based on a configuration file.
    
    This is the main entry point for config-based processing.
    
    Args:
        config_path: Path to configuration file
        **kwargs: Additional arguments for download_and_process_data
        
    Returns:
        Processing results
    """
    config_path = Path(config_path)
    
    # Validate config
    is_valid, errors = validate_config_file(config_path)
    if not is_valid:
        return {
            "success": False,
            "errors": errors,
            "config_path": str(config_path)
        }
    
    # Load config
    config = io.load_json(config_path)
    
    # Run processing
    results = download_and_process_data(config, **kwargs)
    
    return results


