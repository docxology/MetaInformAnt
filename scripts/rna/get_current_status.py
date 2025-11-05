#!/usr/bin/env python3
"""Get current status directly from file system and latest logs."""
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

def get_log_info(species):
    """Get info from latest log."""
    logs = sorted(log_dir.glob(f'workflow_{species}_*.log'), 
                  key=lambda p: p.stat().st_mtime, reverse=True)
    
    if not logs:
        return None
    
    latest = logs[0]
    
    try:
        content = latest.read_text(encoding='utf-8', errors='ignore')
        lines = [l for l in content.split('\n') if l.strip()]
        
        # Get last 20 lines for analysis
        tail = lines[-20:] if len(lines) > 20 else lines
        
        # Find key info
        already_quantified = None
        to_process = None
        current_batch = None
        is_downloading = False
        is_quantifying = False
        is_complete = False
        latest_progress = None
        
        # Check header for initial counts
        for line in lines[:15]:
            if 'Already quantified:' in line:
                try:
                    already_quantified = int(line.split(':')[1].strip())
                except:
                    pass
            if 'To process:' in line:
                try:
                    to_process = int(line.split(':')[1].strip())
                except:
                    pass
        
        # Check tail for current activity
        for line in reversed(tail):
            if 'WORKFLOW COMPLETE' in line:
                is_complete = True
                break
            
            if 'BATCH' in line and ':' in line:
                if not current_batch:
                    parts = line.split('BATCH')
                    if len(parts) > 1:
                        current_batch = parts[1].split(':')[0].strip()
            
            if 'Downloading batch' in line:
                is_downloading = True
            
            if 'Quantifying batch' in line:
                is_quantifying = True
            
            if 'Progress:' in line and '/' in line:
                latest_progress = line.strip()
        
        mtime = datetime.fromtimestamp(latest.stat().st_mtime)
        time_ago = datetime.now() - mtime
        
        # Determine status
        if is_complete:
            status = 'COMPLETE'
        elif is_downloading:
            status = 'DOWNLOADING'
        elif is_quantifying:
            status = 'QUANTIFYING'
        else:
            status = 'IDLE' if time_ago.total_seconds() > 3600 else 'ACTIVE'
        
        return {
            'file': latest.name,
            'status': status,
            'modified': mtime,
            'time_ago': time_ago,
            'already_quantified': already_quantified,
            'to_process': to_process,
            'batch': current_batch,
            'progress': latest_progress,
            'size_kb': latest.stat().st_size // 1024
        }
    except Exception as e:
        return {'file': latest.name, 'error': str(e)}

print("=" * 100)
print("CURRENT AMALGKIT WORKFLOW STATUS")
print("=" * 100)
print(f"\n⏰ {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

total_quantified = 0
total_samples = 0
active_workflows = 0
complete_workflows = 0

for species_id, info in species_info.items():
    # Count from file system (most accurate)
    quantified = count_quantified(species_id)
    total_quantified += quantified
    total_samples += info['total']
    
    pct = (quantified / info['total']) * 100 if info['total'] > 0 else 0
    remaining = info['total'] - quantified
    
    # Get log info
    log_info = get_log_info(species_id)
    
    print(f"{info['name']} ({species_id.upper()})")
    print(f"  📊 Quantified: {quantified}/{info['total']} ({pct:.1f}%)")
    print(f"  📋 Remaining: {remaining}")
    
    if log_info:
        if 'error' in log_info:
            print(f"  ⚠️  Log error: {log_info['error']}")
        else:
            status = log_info['status']
            time_ago_str = str(log_info['time_ago']).split('.')[0]
            
            if status == 'COMPLETE':
                print(f"  ✅ Status: COMPLETE")
                complete_workflows += 1
            elif status == 'DOWNLOADING':
                print(f"  🟢 Status: DOWNLOADING (batch {log_info.get('batch', '?')})")
                active_workflows += 1
            elif status == 'QUANTIFYING':
                print(f"  🟡 Status: QUANTIFYING (batch {log_info.get('batch', '?')})")
                active_workflows += 1
            elif status == 'ACTIVE':
                print(f"  🟢 Status: ACTIVE ({time_ago_str} ago)")
                active_workflows += 1
            else:
                print(f"  ⚪ Status: IDLE ({time_ago_str} ago)")
            
            print(f"  📝 Log: {log_info['file']} ({log_info.get('size_kb', 0)} KB, {time_ago_str} ago)")
            
            if log_info.get('already_quantified'):
                print(f"  📌 Log shows: {log_info['already_quantified']} already quantified, {log_info.get('to_process', '?')} to process")
            
            if log_info.get('progress'):
                print(f"  📈 {log_info['progress']}")
    else:
        print(f"  ⚠️  No log file found")
    
    print()

overall_pct = (total_quantified / total_samples) * 100 if total_samples > 0 else 0

print("=" * 100)
print("OVERALL SUMMARY")
print("=" * 100)
progress_bar = '█' * int(50 * total_quantified / total_samples) + '░' * (50 - int(50 * total_quantified / total_samples))
print(f"Progress: [{progress_bar}] {overall_pct:.1f}%")
print(f"Quantified: {total_quantified}/{total_samples}")
print(f"Remaining: {total_samples - total_quantified}")
print(f"Active: {active_workflows}/4")
print(f"Complete: {complete_workflows}/4")
print("=" * 100)

# Assessment
print("\n💡 ASSESSMENT:")
if active_workflows > 0:
    print(f"  ✅ {active_workflows} workflow(s) actively processing")
    print(f"  ⏳ Estimated remaining: {(total_samples - total_quantified) * 7.5 / 60:.1f} hours")
elif complete_workflows == 4:
    print("  ✅ All workflows complete!")
else:
    print("  ⚠️  Some workflows idle - may need restart")

print()






