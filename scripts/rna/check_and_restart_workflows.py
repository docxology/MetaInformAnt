#!/usr/bin/env python3
"""Check workflow status and restart if needed."""
import subprocess
import sys
from pathlib import Path
from datetime import datetime

repo_root = Path(__file__).parent.parent.parent
config_dir = repo_root / 'config' / 'amalgkit'
output_dir = repo_root / 'output' / 'amalgkit'
log_dir = repo_root / 'output'

species_configs = {
    'cfloridanus': 'amalgkit_cfloridanus.yaml',
    'pbarbatus': 'amalgkit_pbarbatus.yaml',
    'mpharaonis': 'amalgkit_mpharaonis.yaml',
    'sinvicta': 'amalgkit_sinvicta.yaml',
}

def count_quantified(species):
    """Count quantified samples."""
    quant_dir = output_dir / species / 'quant'
    if not quant_dir.exists():
        return 0
    return len([d for d in quant_dir.iterdir() 
                if d.is_dir() and (d / 'abundance.tsv').exists()])

def count_total_from_unquantified(species):
    """Estimate total from unquantified list."""
    unquant_file = output_dir / f'{species}_unquantified.txt'
    if not unquant_file.exists():
        return None
    unquantified = len([l.strip() for l in unquant_file.read_text().splitlines() if l.strip()])
    
    # Read metadata to get total
    metadata_file = output_dir / species / 'work' / 'metadata' / 'metadata.tsv'
    if metadata_file.exists():
        with open(metadata_file) as f:
            lines = f.readlines()
            total = len(lines) - 1 if len(lines) > 0 else None
            if total:
                return total
    return None

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
        
        # Count quantified
        quantified = count_quantified(species)
        total = count_total_from_unquantified(species)
        if total is None:
            # Fallback estimates
            total_estimates = {
                'cfloridanus': 307,
                'pbarbatus': 83,
                'mpharaonis': 100,
                'sinvicta': 354,
            }
            total = total_estimates.get(species, 0)
            print(f"  Quantified: {quantified}/{total} (estimated)")
        else:
            pct = (quantified / total * 100) if total > 0 else 0
            print(f"  Quantified: {quantified}/{total} ({pct:.1f}%)")
        
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

