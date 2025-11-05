#!/usr/bin/env python3
"""Check workflow status and restart if needed."""
import subprocess
import sys
from pathlib import Path
from datetime import datetime
from glob import glob

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.rna import count_quantified_samples, check_workflow_progress
from metainformant.core.logging import get_logger

logger = get_logger("check_restart")

repo_root = Path(__file__).parent.parent.parent
config_dir = repo_root / 'config' / 'amalgkit'
log_dir = repo_root / 'output'

# Discover config files dynamically
def get_species_configs():
    """Discover all species config files."""
    if not config_dir.exists():
        config_dir_alt = repo_root / 'config'
        if config_dir_alt.exists():
            config_pattern = str(config_dir_alt / 'amalgkit_*.yaml')
        else:
            return {}
    else:
        config_pattern = str(config_dir / 'amalgkit_*.yaml')
    
    config_files = sorted(glob(config_pattern))
    species_configs = {}
    
    for config_file in config_files:
        path = Path(config_file)
        if "template" in path.stem.lower():
            continue
        
        species_code = path.stem.replace("amalgkit_", "")
        species_configs[species_code] = path.name
    
    return species_configs

species_configs = get_species_configs()

def check_running_processes():
    """Check if workflow processes are running."""
    try:
        result = subprocess.run(
            ['pgrep', '-f', 'workflow_ena_integrated'],
            capture_output=True,
            text=True
        )
        pids = result.stdout.strip().split('\n') if result.stdout.strip() else []
        return len([p for p in pids if p])
    except:
        return 0

def check_log_activity(species):
    """Check if log file has been modified recently (within last 6 hours)."""
    log_pattern = f'workflow_{species}_*.log'
    logs = sorted(log_dir.glob(log_pattern), key=lambda p: p.stat().st_mtime, reverse=True)
    if not logs:
        return False, None
    
    latest_log = logs[0]
    mtime = latest_log.stat().st_mtime
    now = datetime.now().timestamp()
    hours_ago = (now - mtime) / 3600
    
    return hours_ago < 6, latest_log

def restart_workflow(species, config_file):
    """Restart workflow for a species."""
    config_path = config_dir / config_file
    if not config_path.exists():
        print(f"  ⚠️  Config not found: {config_path}")
        return False
    
    workflow_script = repo_root / 'scripts' / 'rna' / 'workflow_ena_integrated.py'
    if not workflow_script.exists():
        print(f"  ⚠️  Workflow script not found: {workflow_script}")
        return False
    
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = log_dir / f'workflow_{species}_restarted_{timestamp}.log'
    
    cmd = [
        sys.executable,
        str(workflow_script),
        '--config', str(config_path),
        '--batch-size', '12',
        '--threads', '12'
    ]
    
    print(f"  🚀 Restarting workflow for {species}...")
    print(f"     Command: {' '.join(cmd)}")
    print(f"     Log: {log_file}")
    
    try:
        with open(log_file, 'w') as log_f:
            process = subprocess.Popen(
                cmd,
                stdout=log_f,
                stderr=subprocess.STDOUT,
                cwd=str(repo_root)
            )
        print(f"  ✅ Started (PID: {process.pid})")
        return True
    except Exception as e:
        print(f"  ❌ Failed to start: {e}")
        return False

def main():
    print("=" * 80)
    print("WORKFLOW STATUS CHECK AND RESTART")
    print("=" * 80)
    print(f"\n⏰ {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    # Check running processes
    running_count = check_running_processes()
    print(f"Running workflow processes: {running_count}")
    
    # Check each species
    print("\n" + "=" * 80)
    print("SPECIES STATUS")
    print("=" * 80)
    
    needs_restart = []
    
    for species, config_file in species_configs.items():
        print(f"\n{species.upper()}:")
        
        # Use metainformant function to get progress
        config_path = config_dir / config_file
        if config_path.exists():
            progress = check_workflow_progress(config_path)
            quantified = progress["quantified"]
            total = progress["total"]
            pct = progress["percentage"]
            print(f"  Quantified: {quantified}/{total} ({pct:.1f}%)")
        else:
            print(f"  ⚠️  Config not found: {config_path}")
            needs_restart.append((species, config_file))
            continue
        
        # Check log activity
        is_active, latest_log = check_log_activity(species)
        if is_active:
            print(f"  ✅ Recent activity (log: {latest_log.name})")
        else:
            if latest_log:
                hours = (datetime.now().timestamp() - latest_log.stat().st_mtime) / 3600
                print(f"  ⚠️  No recent activity ({hours:.1f} hours ago)")
            else:
                print(f"  ⚠️  No log file found")
            needs_restart.append((species, config_file))
    
    # Restart if needed
    if needs_restart:
        print("\n" + "=" * 80)
        print("RESTARTING WORKFLOWS")
        print("=" * 80)
        for species, config_file in needs_restart:
            restart_workflow(species, config_file)
    else:
        print("\n✅ All workflows appear to be running")
    
    print("\n" + "=" * 80)
    print("COMPLETE")
    print("=" * 80)

if __name__ == '__main__':
    main()

