#!/usr/bin/env python3
"""Check current workflow progress from logs and file system."""
from pathlib import Path
from datetime import datetime

repo_root = Path(__file__).parent.parent.parent
output_dir = repo_root / 'output' / 'amalgkit'
log_dir = repo_root / 'output'

species_info = {
    'cfloridanus': {'total': 307, 'name': 'Camponotus floridanus'},
    'pbarbatus': {'total': 83, 'name': 'Pogonomyrmex barbatus'},
    'mpharaonis': {'total': 100, 'name': 'Monomorium pharaonis'},
    'sinvicta': {'total': 354, 'name': 'Solenopsis invicta'},
}

def count_quantified(species):
    """Count quantified samples."""
    quant_dir = output_dir / species / 'quant'
    if not quant_dir.exists():
        return 0
    return len([d for d in quant_dir.iterdir() 
                if d.is_dir() and (d / 'abundance.tsv').exists()])

def read_log_tail(log_path, n_lines=50):
    """Read last N lines of log file."""
    try:
        with open(log_path, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()
            return lines[-n_lines:] if len(lines) > n_lines else lines
    except:
        return []

def analyze_latest_activity(species):
    """Analyze latest log file for current activity."""
    logs = sorted(log_dir.glob(f'workflow_{species}_*.log'), 
                  key=lambda p: p.stat().st_mtime, reverse=True)
    
    if not logs:
        return None
    
    latest = logs[0]
    mtime = datetime.fromtimestamp(latest.stat().st_mtime)
    
    # Read last 50 lines
    tail_lines = read_log_tail(latest, 50)
    
    # Parse activity
    current_batch = None
    current_activity = None
    latest_progress = None
    downloaded_count = None
    quantified_count = None
    is_complete = False
    
    for line in reversed(tail_lines):
        line_clean = line.strip()
        
        if 'WORKFLOW COMPLETE' in line_clean:
            is_complete = True
            break
        
        if 'BATCH' in line_clean and ':' in line_clean:
            if not current_batch:
                # Extract batch info
                if 'BATCH' in line_clean:
                    parts = line_clean.split('BATCH')
                    if len(parts) > 1:
                        current_batch = parts[1].strip()
        
        if 'Downloading batch' in line_clean:
            current_activity = 'DOWNLOADING'
        
        if 'Quantifying batch' in line_clean:
            current_activity = 'QUANTIFYING'
        
        if 'Progress:' in line_clean and '/' in line_clean:
            latest_progress = line_clean
        
        if 'Downloaded:' in line_clean and 'samples' in line_clean:
            # Extract count
            parts = line_clean.split('Downloaded:')
            if len(parts) > 1:
                try:
                    downloaded_count = int(parts[1].split('/')[0].strip())
                except:
                    pass
            if current_activity == 'DOWNLOADING':
                current_activity = 'QUANTIFYING'  # Transition
        
        if 'Quantified:' in line_clean and 'samples' in line_clean:
            parts = line_clean.split('Quantified:')
            if len(parts) > 1:
                try:
                    quantified_count = int(parts[1].split('/')[0].strip())
                except:
                    pass
    
    time_ago = datetime.now() - mtime
    
    return {
        'file': latest.name,
        'modified': mtime,
        'time_ago': time_ago,
        'batch': current_batch,
        'activity': current_activity,
        'progress': latest_progress,
        'downloaded': downloaded_count,
        'quantified': quantified_count,
        'complete': is_complete
    }

def main():
    print("=" * 100)
    print("COMPREHENSIVE WORKFLOW STATUS CHECK")
    print("=" * 100)
    print(f"\n⏰ {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    total_quantified = 0
    total_samples = 0
    active_count = 0
    complete_count = 0
    
    for species_id, info in species_info.items():
        quantified = count_quantified(species_id)
        total_quantified += quantified
        total_samples += info['total']
        
        pct = (quantified / info['total']) * 100 if info['total'] > 0 else 0
        remaining = info['total'] - quantified
        
        log_info = analyze_latest_activity(species_id)
        
        print("=" * 100)
        print(f"{info['name']} ({species_id.upper()})")
        print("=" * 100)
        print(f"📊 Progress: {quantified}/{info['total']} ({pct:.1f}%) | Remaining: {remaining}")
        
        if log_info:
            if log_info['complete']:
                status_icon = '✅'
                status_text = 'COMPLETE'
                complete_count += 1
            elif log_info['activity'] == 'DOWNLOADING':
                status_icon = '🟢'
                status_text = f'DOWNLOADING (batch {log_info.get("batch", "?")})'
                active_count += 1
            elif log_info['activity'] == 'QUANTIFYING':
                status_icon = '🟡'
                status_text = f'QUANTIFYING (batch {log_info.get("batch", "?")})'
                active_count += 1
            else:
                status_icon = '⚪'
                status_text = 'IDLE'
            
            print(f"Status: {status_icon} {status_text}")
            print(f"Log: {log_info['file']}")
            
            time_ago_str = str(log_info['time_ago']).split('.')[0]
            print(f"Last update: {time_ago_str} ago")
            
            if log_info.get('batch'):
                print(f"Current batch: {log_info['batch']}")
            
            if log_info.get('progress'):
                print(f"Latest: {log_info['progress']}")
            
            if log_info.get('downloaded') is not None:
                print(f"Last batch downloaded: {log_info['downloaded']} samples")
            
            if log_info.get('quantified') is not None:
                print(f"Last batch quantified: {log_info['quantified']} samples")
        else:
            print("Status: ⚪ NO LOG FOUND")
        
        print()
    
    overall_pct = (total_quantified / total_samples) * 100 if total_samples > 0 else 0
    
    print("=" * 100)
    print("OVERALL SUMMARY")
    print("=" * 100)
    print(f"Total Progress: {total_quantified}/{total_samples} ({overall_pct:.1f}%)")
    print(f"Remaining: {total_samples - total_quantified} samples")
    print(f"Active workflows: {active_count}/4")
    print(f"Complete workflows: {complete_count}/4")
    print("=" * 100)
    
    if active_count == 0 and complete_count < 4:
        print("\n💡 All workflows appear idle. Consider restarting:")
        print("   python3 scripts/rna/restart_all_workflows.py")
    elif active_count > 0:
        print(f"\n✅ {active_count} workflow(s) actively processing")
        print("   Monitor with: python3 scripts/rna/monitor_comprehensive.py")
    else:
        print("\n✅ All workflows complete!")
    
    print()

if __name__ == '__main__':
    main()



