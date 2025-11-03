#!/usr/bin/env python3
"""Restart all amalgkit workflows with nohup."""
import subprocess
import sys
from pathlib import Path
from datetime import datetime

repo_root = Path(__file__).parent.parent.parent
workflow_script = repo_root / 'scripts' / 'rna' / 'workflow_ena_integrated.py'
config_dir = repo_root / 'config' / 'amalgkit'
log_dir = repo_root / 'output'

species_configs = [
    ('cfloridanus', 'amalgkit_cfloridanus.yaml'),
    ('pbarbatus', 'amalgkit_pbarbatus.yaml'),
    ('mpharaonis', 'amalgkit_mpharaonis.yaml'),
    ('sinvicta', 'amalgkit_sinvicta.yaml'),
]

def restart_species(species, config_file):
    """Restart workflow for a species using nohup."""
    config_path = config_dir / config_file
    if not config_path.exists():
        print(f"  ❌ Config not found: {config_path}")
        return False
    
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = log_dir / f'workflow_{species}_restarted_{timestamp}.log'
    
    cmd = [
        'nohup',
        sys.executable,
        str(workflow_script),
        '--config', str(config_path),
        '--batch-size', '12',
        '--threads', '12'
    ]
    
    print(f"  🚀 {species}: Starting workflow...")
    print(f"     Log: {log_file.name}")
    
    try:
        with open(log_file, 'w') as log_f:
            process = subprocess.Popen(
                cmd,
                stdout=log_f,
                stderr=subprocess.STDOUT,
                cwd=str(repo_root),
                start_new_session=True
            )
        print(f"     ✅ Started (PID: {process.pid})")
        return True
    except Exception as e:
        print(f"     ❌ Failed: {e}")
        return False

def main():
    print("=" * 80)
    print("RESTARTING ALL AMALGKIT WORKFLOWS")
    print("=" * 80)
    print(f"\n⏰ {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    print("Starting workflows for all 4 species in parallel...\n")
    
    for species, config_file in species_configs:
        restart_species(species, config_file)
        print()
    
    print("=" * 80)
    print("✅ ALL WORKFLOWS STARTED")
    print("=" * 80)
    print("\nMonitor progress with:")
    print("  python3 scripts/rna/monitor_comprehensive.py")
    print("\nCheck logs:")
    print("  tail -f output/workflow_*_restarted_*.log")

if __name__ == '__main__':
    main()

