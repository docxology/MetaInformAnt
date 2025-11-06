#!/usr/bin/env python3
"""Core utilities simulation script.

This script generates test data for core utilities including config files,
workflow test data, and I/O pattern examples.

Usage:
    python3 scripts/simulation/simulate_core.py --type config --n-configs 5
    python3 scripts/simulation/simulate_core.py --type workflow --n-steps 10
    python3 scripts/simulation/simulate_core.py --type io --n-files 20
"""

import argparse
import random
import sys
from pathlib import Path

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core import io, logging, paths, validation

logger = logging.get_logger(__name__)


def simulate_config(
    output_dir: Path,
    n_configs: int,
    config_type: str,
    seed: int,
) -> dict:
    """Simulate configuration files."""
    logger.info(f"Generating config files: {n_configs} configs, type {config_type}")
    rng = random.Random(seed)
    
    configs = []
    
    for config_idx in range(n_configs):
        if config_type == "yaml":
            config_content = f"""# Simulated config {config_idx}
work_dir: output/test_{config_idx}/
log_dir: output/test_{config_idx}/logs
threads: {rng.randint(1, 16)}
seed: {rng.randint(1, 1000)}

parameters:
  param1: {rng.uniform(0, 1):.3f}
  param2: {rng.randint(10, 100)}
  param3: "{rng.choice(['option_a', 'option_b', 'option_c'])}"

data:
  input_file: data/input_{config_idx}.csv
  output_file: output/test_{config_idx}/results.json
"""
            ext = "yaml"
        elif config_type == "json":
            config_content = {
                "work_dir": f"output/test_{config_idx}/",
                "log_dir": f"output/test_{config_idx}/logs",
                "threads": rng.randint(1, 16),
                "seed": rng.randint(1, 1000),
                "parameters": {
                    "param1": round(rng.uniform(0, 1), 3),
                    "param2": rng.randint(10, 100),
                    "param3": rng.choice(["option_a", "option_b", "option_c"]),
                },
                "data": {
                    "input_file": f"data/input_{config_idx}.csv",
                    "output_file": f"output/test_{config_idx}/results.json",
                },
            }
            ext = "json"
        else:  # toml
            config_content = f"""[general]
work_dir = "output/test_{config_idx}/"
log_dir = "output/test_{config_idx}/logs"
threads = {rng.randint(1, 16)}
seed = {rng.randint(1, 1000)}

[parameters]
param1 = {rng.uniform(0, 1):.3f}
param2 = {rng.randint(10, 100)}
param3 = "{rng.choice(['option_a', 'option_b', 'option_c'])}"

[data]
input_file = "data/input_{config_idx}.csv"
output_file = "output/test_{config_idx}/results.json"
"""
            ext = "toml"
        
        config_file = output_dir / f"config_{config_idx:03d}.{ext}"
        
        if config_type == "json":
            io.dump_json(config_content, config_file, indent=2)
        else:
            config_file.write_text(config_content)
        
        configs.append({
            "config_id": f"config_{config_idx:03d}",
            "file": str(config_file),
            "type": config_type,
        })
    
    # Save index
    index_file = output_dir / "config_index.json"
    io.dump_json(configs, index_file, indent=2)
    
    logger.info(f"Config files saved to {output_dir}")
    
    return {
        "type": "config",
        "n_configs": n_configs,
        "config_type": config_type,
        "index_file": str(index_file),
    }


def simulate_workflow(
    output_dir: Path,
    n_steps: int,
    seed: int,
) -> dict:
    """Simulate workflow test data."""
    logger.info(f"Generating workflow test data: {n_steps} steps")
    rng = random.Random(seed)
    
    step_types = ["input", "process", "output", "validation"]
    
    workflow = {
        "workflow_id": "test_workflow",
        "steps": [],
    }
    
    for step_idx in range(n_steps):
        step_type = rng.choice(step_types)
        step = {
            "step_id": f"step_{step_idx:03d}",
            "step_type": step_type,
            "name": f"Step {step_idx}",
            "dependencies": [],
        }
        
        # Add dependencies (earlier steps)
        if step_idx > 0:
            n_deps = rng.randint(0, min(3, step_idx))
            deps = rng.sample(range(step_idx), n_deps)
            step["dependencies"] = [f"step_{d:03d}" for d in deps]
        
        # Add step-specific data
        if step_type == "input":
            step["input_file"] = f"data/input_{step_idx}.csv"
        elif step_type == "process":
            step["function"] = f"process_{rng.choice(['transform', 'filter', 'aggregate'])}"
            step["parameters"] = {"param1": rng.uniform(0, 1)}
        elif step_type == "output":
            step["output_file"] = f"output/step_{step_idx}_result.json"
        else:  # validation
            step["validation_type"] = rng.choice(["check", "assert", "verify"])
        
        workflow["steps"].append(step)
    
    # Save workflow
    workflow_file = output_dir / "workflow.json"
    io.dump_json(workflow, workflow_file, indent=2)
    
    logger.info(f"Workflow data saved to {workflow_file}")
    
    return {
        "type": "workflow",
        "n_steps": n_steps,
        "output_file": str(workflow_file),
    }


def simulate_io(
    output_dir: Path,
    n_files: int,
    file_types: list[str],
    seed: int,
) -> dict:
    """Simulate I/O test data."""
    logger.info(f"Generating I/O test data: {n_files} files")
    rng = random.Random(seed)
    
    files_created = []
    
    for file_idx in range(n_files):
        file_type = rng.choice(file_types) if file_types else "json"
        file_name = f"test_data_{file_idx:03d}.{file_type}"
        file_path = output_dir / file_name
        
        if file_type == "json":
            data = {
                "id": file_idx,
                "data": [rng.randint(1, 100) for _ in range(10)],
                "metadata": {"created": "2024-01-01", "version": "1.0"},
            }
            io.dump_json(data, file_path, indent=2)
        elif file_type == "csv":
            import pandas as pd
            
            df = pd.DataFrame({
                "col1": [rng.randint(1, 100) for _ in range(10)],
                "col2": [rng.uniform(0, 1) for _ in range(10)],
                "col3": [rng.choice(["A", "B", "C"]) for _ in range(10)],
            })
            df.to_csv(file_path, index=False)
        elif file_type == "txt":
            content = "\n".join([f"Line {i}: {rng.randint(1, 100)}" for i in range(10)])
            file_path.write_text(content)
        else:  # jsonl
            with open(file_path, "w") as f:
                for i in range(10):
                    record = {"id": i, "value": rng.randint(1, 100)}
                    f.write(io.dumps_json(record) + "\n")
        
        files_created.append({
            "file_id": f"file_{file_idx:03d}",
            "file_path": str(file_path),
            "file_type": file_type,
        })
    
    # Save index
    index_file = output_dir / "io_index.json"
    io.dump_json(files_created, index_file, indent=2)
    
    logger.info(f"I/O test files saved to {output_dir}")
    
    return {
        "type": "io",
        "n_files": n_files,
        "file_types": file_types,
        "index_file": str(index_file),
    }


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Core utilities simulation for testing and validation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Simulate config files
  %(prog)s --type config --n-configs 5 --config-type yaml

  # Simulate workflow data
  %(prog)s --type workflow --n-steps 10

  # Simulate I/O test files
  %(prog)s --type io --n-files 20 --file-types json csv
        """,
    )
    
    parser.add_argument(
        "--type",
        required=True,
        choices=["config", "workflow", "io"],
        help="Type of simulation to run",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/simulation/core"),
        help="Output directory (default: output/simulation/core)",
    )
    parser.add_argument("--n-configs", type=int, default=5, help="Number of config files (config type)")
    parser.add_argument(
        "--config-type",
        type=str,
        choices=["yaml", "json", "toml"],
        default="yaml",
        help="Config file format (config type)",
    )
    parser.add_argument("--n-steps", type=int, default=10, help="Number of workflow steps (workflow type)")
    parser.add_argument(
        "--n-files",
        type=int,
        default=20,
        help="Number of I/O test files (io type)",
    )
    parser.add_argument(
        "--file-types",
        nargs="+",
        choices=["json", "csv", "txt", "jsonl"],
        default=["json", "csv"],
        help="File types to generate (io type)",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")
    
    args = parser.parse_args()
    
    if args.verbose:
        import logging as std_logging
        logger.setLevel(std_logging.DEBUG)
    
    output_dir = paths.ensure_directory(args.output)
    
    try:
        if args.type == "config":
            results = simulate_config(output_dir, args.n_configs, args.config_type, args.seed)
        elif args.type == "workflow":
            results = simulate_workflow(output_dir, args.n_steps, args.seed)
        elif args.type == "io":
            results = simulate_io(output_dir, args.n_files, args.file_types, args.seed)
        
        # Save summary
        summary_file = output_dir / "simulation_summary.json"
        io.dump_json(results, summary_file, indent=2)
        logger.info(f"Simulation complete. Summary saved to {summary_file}")
        
        return 0
    except Exception as e:
        logger.error(f"Simulation failed: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())

