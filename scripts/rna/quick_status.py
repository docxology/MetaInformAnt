#!/usr/bin/env python3
"""Quick status check - counts quantified samples and checks log activity."""
from pathlib import Path
from datetime import datetime

repo_root = Path(__file__).parent.parent.parent
output_dir = repo_root / 'output' / 'amalgkit'
log_dir = repo_root / 'output'

species_info = {
    'cfloridanus': {'total': 307, 'name': 'C. floridanus'},
    'pbarbatus': {'total': 83, 'name': 'P. barbatus'},
    'mpharaonis': {'total': 100, 'name': 'M. pharaonis'},
    'sinvicta': {'total': 354, 'name': 'S. invicta'},
}

def count_quantified(species):
    """Count quantified samples."""
    quant_dir = output_dir / species / 'quant'
    if not quant_dir.exists():
        return 0
    return len([d for d in quant_dir.iterdir() 
                if d.is_dir() and (d / 'abundance.tsv').exists()])

def get_latest_log_status(species):
    """Get status from latest log file."""
    logs = sorted(log_dir.glob(f'workflow_{species}_*.log'), 
                  key=lambda p: p.stat().st_mtime, reverse=True)
    if not logs:
        return None, None
    
    latest = logs[0]
    
    # Read last 30 lines
    try:
        content = latest.read_text()
        lines = [l for l in content.split('\n') if l.strip()]
        last_lines = lines[-30:] if len(lines) > 30 else lines
        
        # Find latest activity
        latest_batch = None
        latest_progress = None
        is_downloading = False
        is_quantifying = False
        
        for line in reversed(last_lines):
            if 'BATCH' in line and ':' in line and not latest_batch:
                latest_batch = line.strip()
            if 'Progress:' in line and '/' in line:
                latest_progress = line.strip()
            if 'Downloading batch' in line:
                is_downloading = True
            if 'Quantifying batch' in line:
                is_quantifying = True
            if 'WORKFLOW COMPLETE' in line:
                return 'COMPLETE', latest
            if 'Downloaded:' in line and 'samples' in line:
                is_downloading = False
        
        # Determine status
        if is_downloading:
            status = 'DOWNLOADING'
        elif is_quantifying:
            status = 'QUANTIFYING'
        elif latest_progress:
            status = 'ACTIVE'
        else:
            status = 'IDLE'
        
        return status, latest
    except:
        return 'ERROR', latest

print("=" * 90)
print("AMALGKIT WORKFLOW STATUS CHECK")
print("=" * 90)
print(f"\n⏰ {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

total_quantified = 0
total_samples = 0
active_count = 0

for species_id, info in species_info.items():
    quantified = count_quantified(species_id)
    total_quantified += quantified
    total_samples += info['total']
    
    pct = (quantified / info['total']) * 100 if info['total'] > 0 else 0
    remaining = info['total'] - quantified
    
    status, log_file = get_latest_log_status(species_id)
    
    if status in ['DOWNLOADING', 'QUANTIFYING', 'ACTIVE']:
        active_count += 1
        status_icon = '🟢'
    elif status == 'COMPLETE':
        status_icon = '✅'
    else:
        status_icon = '⚪'
    
    print(f"{status_icon} {info['name']} ({species_id})")
    print(f"   Progress: {quantified}/{info['total']} ({pct:.1f}%)")
    print(f"   Remaining: {remaining}")
    if status and log_file:
        log_name = log_file.name if isinstance(log_file, Path) else str(log_file)
        print(f"   Status: {status}")
        print(f"   Log: {log_name}")
    print()

overall_pct = (total_quantified / total_samples) * 100 if total_samples > 0 else 0

print("=" * 90)
print("SUMMARY")
print("=" * 90)
print(f"Total: {total_quantified}/{total_samples} ({overall_pct:.1f}%)")
print(f"Remaining: {total_samples - total_quantified}")
print(f"Active workflows: {active_count}/4")
print("=" * 90)



