#!/usr/bin/env python3
"""
Real-time monitoring dashboard for multi-species RNA-seq workflow.

Shows progress for all species, batch status, disk usage, and running processes.
"""

import sys
import time
from datetime import datetime
from pathlib import Path
from collections import defaultdict

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.rna import check_workflow_progress
from metainformant.rna.workflow import load_workflow_config
from glob import glob


def get_config_for_species_dir(species_dir: Path) -> Path | None:
    """Find config file for a species directory."""
    species_name = species_dir.name
    repo_root = Path(__file__).parent.parent.parent
    config_dir = repo_root / "config" / "amalgkit"
    
    config_pattern = str(config_dir / f"amalgkit_{species_name}.yaml")
    config_files = glob(config_pattern)
    if config_files:
        return Path(config_files[0])
    
    # Try alternative naming
    config_pattern = str(config_dir / f"amalgkit_{species_name.replace('_', '*')}.yaml")
    config_files = glob(config_pattern)
    if config_files:
        return Path(config_files[0])
    
    return None


def get_disk_usage(path: Path) -> tuple[float, str]:
    """Get disk usage in GB and formatted string."""
    if not path.exists():
        return 0.0, "0 B"
    
    import subprocess
    result = subprocess.run(
        ["du", "-sh", str(path)],
        capture_output=True,
        text=True
    )
    size_str = result.stdout.split()[0] if result.returncode == 0 else "?"
    
    # Try to get numeric value
    try:
        result = subprocess.run(
            ["du", "-sk", str(path)],
            capture_output=True,
            text=True
        )
        kb = float(result.stdout.split()[0])
        gb = kb / 1024 / 1024
        return gb, size_str
    except Exception:
        return 0.0, size_str


def check_running_processes() -> list[dict]:
    """Check for running workflow processes."""
    import subprocess
    result = subprocess.run(
        ["ps", "aux"],
        capture_output=True,
        text=True
    )
    
    processes = []
    for line in result.stdout.splitlines():
        if "run_multi_species" in line or "amalgkit" in line:
            if "grep" not in line:
                parts = line.split()
                if len(parts) >= 11:
                    processes.append({
                        "pid": parts[1],
                        "cpu": parts[2],
                        "mem": parts[3],
                        "time": parts[9],
                        "command": " ".join(parts[10:])[:60]
                    })
    return processes


def get_latest_batch_logs(species_dir: Path) -> dict:
    """Get info about latest batch processing."""
    logs_dir = species_dir / "logs"
    if not logs_dir.exists():
        return {}
    
    # Find latest getfastq batch log
    getfastq_logs = sorted(logs_dir.glob("*getfastq_batch*.stdout.log"))
    quant_logs = sorted(logs_dir.glob("*quant_batch*.stdout.log"))
    
    info = {}
    
    if getfastq_logs:
        latest = getfastq_logs[-1]
        info["latest_download"] = latest.stem
        info["download_time"] = datetime.fromtimestamp(latest.stat().st_mtime)
    
    if quant_logs:
        latest = quant_logs[-1]
        info["latest_quant"] = latest.stem
        info["quant_time"] = datetime.fromtimestamp(latest.stat().st_mtime)
    
    return info


def display_dashboard():
    """Display the monitoring dashboard."""
    base_dir = Path("output/amalgkit")
    
    if not base_dir.exists():
        print("❌ No amalgkit output directory found")
        return
    
    # Clear screen
    print("\033[2J\033[H", end="")
    
    print("=" * 100)
    print(f"{'METAINFORMANT MULTI-SPECIES RNA-SEQ WORKFLOW MONITOR':^100}")
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S'):^100}")
    print("=" * 100)
    print()
    
    # Species progress
    print("📊 SPECIES PROGRESS")
    print("-" * 100)
    print(f"{'Species':<25} {'Progress':<20} {'Quantified':<15} {'Disk (FASTQ)':<15} {'Disk (Quant)':<15}")
    print("-" * 100)
    
    total_quantified = 0
    total_samples = 0
    
    species_dirs = sorted(base_dir.glob("*/"))
    for species_dir in species_dirs:
        if not species_dir.is_dir():
            continue
        
        species_name = species_dir.name.capitalize()
        
        # Try to get config and use metainformant function
        config_path = get_config_for_species_dir(species_dir)
        if config_path and config_path.exists():
            progress_info = check_workflow_progress(config_path)
            quantified = progress_info["quantified"]
            total_count = progress_info["total"]
            pct = progress_info["percentage"]
        else:
            # Fallback: count directly
            quant_dir = species_dir / "quant"
            quantified = len(list(quant_dir.glob("*/abundance.tsv"))) if quant_dir.exists() else 0
            total_count = 0
            pct = 0.0
        
        total_quantified += quantified
        total_samples += total_count
        
        if total_count > 0:
            progress = f"{'█' * int(pct / 5):<20}"
            status = f"{quantified}/{total_count} ({pct:.1f}%)"
        else:
            progress = "?" * 20
            status = "No metadata"
        
        fastq_gb, fastq_str = get_disk_usage(species_dir / "fastq")
        quant_gb, quant_str = get_disk_usage(species_dir / "quant")
        
        # Color code based on progress
        if quantified == total_count and total_count > 0:
            color = "\033[92m"  # Green
        elif quantified > 0:
            color = "\033[93m"  # Yellow
        else:
            color = "\033[90m"  # Gray
        reset = "\033[0m"
        
        print(f"{color}{species_name:<25} {progress:<20} {status:<15} {fastq_str:<15} {quant_str:<15}{reset}")
    
    print("-" * 100)
    print(f"{'TOTAL':<25} {'':<20} {total_quantified}/{total_samples:<15} {'':<15} {'':<15}")
    print()
    
    # Running processes
    processes = check_running_processes()
    print("⚙️  RUNNING PROCESSES")
    print("-" * 100)
    if processes:
        for proc in processes:
            print(f"PID {proc['pid']} | CPU: {proc['cpu']}% | MEM: {proc['mem']}% | TIME: {proc['time']} | {proc['command']}")
    else:
        print("No workflow processes running")
    print()
    
    # Latest batch info
    print("📦 LATEST BATCH ACTIVITY")
    print("-" * 100)
    for species_dir in species_dirs[:4]:  # Show first 4
        if not species_dir.is_dir():
            continue
        species_name = species_dir.name.capitalize()
        batch_info = get_latest_batch_logs(species_dir)
        
        if batch_info:
            download_time = batch_info.get("download_time", "N/A")
            quant_time = batch_info.get("quant_time", "N/A")
            print(f"{species_name:<20} Download: {download_time} | Quant: {quant_time}")
        else:
            print(f"{species_name:<20} No batch logs yet")
    
    print()
    print("=" * 100)
    print("Press Ctrl+C to exit | Refreshes every 30 seconds")
    print("=" * 100)


def main():
    """Main monitoring loop."""
    try:
        while True:
            display_dashboard()
            time.sleep(30)
    except KeyboardInterrupt:
        print("\n\n👋 Monitoring stopped")
        sys.exit(0)


if __name__ == "__main__":
    main()

