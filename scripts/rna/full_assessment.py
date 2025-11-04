#!/usr/bin/env python3
"""Comprehensive assessment of all amalgkit workflows."""
from pathlib import Path
from datetime import datetime, timedelta
import sys

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

def count_downloaded_not_quantified(species):
    """Count samples with FASTQ/SRA files but not quantified."""
    fastq_dir = output_dir / species / 'fastq'
    if not fastq_dir.exists():
        return 0
    
    count = 0
    for sample_dir in fastq_dir.iterdir():
        if not sample_dir.is_dir():
            continue
        
        # Check for data files
        has_data = bool(
            list(sample_dir.rglob('*.fastq.gz')) or 
            list(sample_dir.rglob('*.sra'))
        )
        
        if has_data:
            # Check if quantified
            quant_file = output_dir / species / 'quant' / sample_dir.name / 'abundance.tsv'
            if not quant_file.exists():
                count += 1
    
    return count

def analyze_log(species):
    """Analyze workflow log for activity."""
    logs = sorted(log_dir.glob(f'workflow_{species}_*.log'), 
                  key=lambda p: p.stat().st_mtime, reverse=True)
    
    if not logs:
        return None, None, None, None, None
    
    latest_log = logs[0]
    
    # Read log content
    try:
        content = latest_log.read_text()
        lines = content.split('\n')
        
        # Find latest batch info
        latest_batch = None
        latest_download = None
        latest_quant = None
        latest_progress = None
        
        for line in reversed(lines[-100:]):  # Check last 100 lines
            if 'BATCH' in line and ':' in line:
                if not latest_batch:
                    # Extract batch number
                    parts = line.split('BATCH')
                    if len(parts) > 1:
                        latest_batch = parts[1].strip()
            
            if 'Downloaded:' in line:
                latest_download = line.strip()
            
            if 'Quantified:' in line and 'samples' in line:
                latest_quant = line.strip()
            
            if 'Progress:' in line and '/' in line:
                latest_progress = line.strip()
        
        # Get file modification time
        mtime = datetime.fromtimestamp(latest_log.stat().st_mtime)
        time_ago = datetime.now() - mtime
        
        return {
            'file': latest_log.name,
            'modified': mtime,
            'time_ago': time_ago,
            'batch': latest_batch,
            'download': latest_download,
            'quant': latest_quant,
            'progress': latest_progress,
            'size_kb': latest_log.stat().st_size // 1024
        }
    except Exception as e:
        return {
            'file': latest_log.name,
            'error': str(e)
        }

def get_unquantified_list(species):
    """Get unquantified sample list."""
    unquant_file = output_dir / f'{species}_unquantified.txt'
    if not unquant_file.exists():
        return []
    return [l.strip() for l in unquant_file.read_text().splitlines() if l.strip()]

def calculate_directory_size(path):
    """Calculate directory size in bytes."""
    total = 0
    try:
        for file_path in path.rglob('*'):
            if file_path.is_file():
                try:
                    total += file_path.stat().st_size
                except:
                    pass
    except:
        pass
    return total

def format_size(size_bytes):
    """Format bytes to human-readable size."""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f} PB"

def main():
    print("=" * 100)
    print("COMPREHENSIVE AMALGKIT WORKFLOW ASSESSMENT")
    print("=" * 100)
    print(f"\n⏰ {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    total_quantified = 0
    total_samples = 0
    total_downloaded = 0
    total_unquantified = 0
    active_workflows = 0
    
    for species_id, info in species_info.items():
        print("=" * 100)
        print(f"{info['name']} ({species_id.upper()})")
        print("=" * 100)
        
        # Quantified count
        quantified = count_quantified(species_id)
        total_quantified += quantified
        total_samples += info['total']
        
        # Downloaded but not quantified
        downloaded = count_downloaded_not_quantified(species_id)
        total_downloaded += downloaded
        
        # Unquantified list
        unquantified = get_unquantified_list(species_id)
        total_unquantified += len(unquantified)
        
        # Progress
        pct = (quantified / info['total']) * 100 if info['total'] > 0 else 0
        remaining = info['total'] - quantified
        progress_bar = '█' * int(50 * quantified / info['total']) + '░' * (50 - int(50 * quantified / info['total']))
        
        print(f"\n📊 PROGRESS")
        print(f"   [{progress_bar}] {pct:.1f}%")
        print(f"   Quantified: {quantified}/{info['total']}")
        print(f"   Remaining: {remaining}")
        
        # Disk usage
        fastq_dir = output_dir / species_id / 'fastq'
        quant_dir = output_dir / species_id / 'quant'
        
        fastq_size = calculate_directory_size(fastq_dir) if fastq_dir.exists() else 0
        quant_size = calculate_directory_size(quant_dir) if quant_dir.exists() else 0
        
        print(f"\n💾 STORAGE")
        print(f"   FASTQ (temporary): {format_size(fastq_size)}")
        print(f"   Quantification (permanent): {format_size(quant_size)}")
        print(f"   Downloaded (not quantified): {downloaded} samples")
        
        # Log analysis
        log_info = analyze_log(species_id)
        if log_info:
            print(f"\n📝 WORKFLOW STATUS")
            if 'error' in log_info:
                print(f"   ⚠️  Error reading log: {log_info['error']}")
            else:
                time_ago_str = str(log_info['time_ago']).split('.')[0] if log_info.get('time_ago') else 'unknown'
                is_recent = log_info.get('time_ago') and log_info['time_ago'] < timedelta(hours=1)
                status = "🟢 ACTIVE" if is_recent else f"🔴 Inactive ({time_ago_str})"
                
                print(f"   Log: {log_info['file']} ({log_info.get('size_kb', 0)} KB)")
                print(f"   Status: {status}")
                print(f"   Last modified: {log_info.get('modified', 'unknown')}")
                
                if log_info.get('batch'):
                    print(f"   Latest batch: {log_info['batch']}")
                if log_info.get('progress'):
                    print(f"   {log_info['progress']}")
                if log_info.get('download'):
                    print(f"   {log_info['download']}")
                if log_info.get('quant'):
                    print(f"   {log_info['quant']}")
                
                if is_recent:
                    active_workflows += 1
        
        print(f"\n📋 UNQUANTIFIED")
        print(f"   From list file: {len(unquantified)} samples")
        if len(unquantified) > 0 and len(unquantified) <= 5:
            print(f"   Samples: {', '.join(unquantified)}")
        elif len(unquantified) > 5:
            print(f"   Samples: {', '.join(unquantified[:3])} ... and {len(unquantified) - 3} more")
        
        print()
    
    # Overall summary
    print("=" * 100)
    print("OVERALL SUMMARY")
    print("=" * 100)
    
    overall_pct = (total_quantified / total_samples) * 100 if total_samples > 0 else 0
    progress_bar = '█' * int(50 * total_quantified / total_samples) + '░' * (50 - int(50 * total_quantified / total_samples))
    
    print(f"\n📊 PROGRESS")
    print(f"   [{progress_bar}] {overall_pct:.1f}%")
    print(f"   Total Samples: {total_samples}")
    print(f"   Quantified: {total_quantified}/{total_samples}")
    print(f"   Remaining: {total_samples - total_quantified}")
    print(f"   Downloaded (not quantified): {total_downloaded}")
    print(f"   Unquantified list: {total_unquantified}")
    print(f"   Active workflows: {active_workflows}/4")
    
    # Time estimate
    remaining = total_samples - total_quantified
    if remaining > 0:
        # Conservative estimate: 10 minutes per sample (includes download failures)
        est_minutes = remaining * 10
        est_hours = est_minutes / 60
        est_days = est_hours / 24
        print(f"\n⏳ ESTIMATED COMPLETION")
        if est_days >= 1:
            print(f"   {est_days:.1f} days ({est_hours:.1f} hours)")
        else:
            print(f"   {est_hours:.1f} hours")
    
    # Health check
    print(f"\n🏥 HEALTH STATUS")
    if active_workflows == 4:
        print("   ✅ All workflows active")
    elif active_workflows > 0:
        print(f"   ⚠️  {active_workflows}/4 workflows active")
    else:
        print("   ❌ No active workflows detected")
    
    if total_downloaded > 100:
        print(f"   ⚠️  Large backlog: {total_downloaded} downloaded samples not quantified")
    else:
        print(f"   ✅ Download/quantify balance healthy ({total_downloaded} pending)")
    
    print("=" * 100)
    print()

if __name__ == '__main__':
    main()





