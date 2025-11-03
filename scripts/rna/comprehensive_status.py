#!/usr/bin/env python3
"""Comprehensive status check for all amalgkit workflows."""
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

def count_downloaded(species):
    """Count downloaded samples (with FASTQ or SRA files)."""
    fastq_dir = output_dir / species / 'fastq'
    if not fastq_dir.exists():
        return 0
    return len([d for d in fastq_dir.iterdir() 
                if d.is_dir() and (
                    list(d.rglob('*.fastq.gz')) or 
                    list(d.rglob('*.sra'))
                )])

def get_unquantified_count(species):
    """Get unquantified count from list file."""
    unquant_file = output_dir / f'{species}_unquantified.txt'
    if not unquant_file.exists():
        return None
    return len([l.strip() for l in unquant_file.read_text().splitlines() if l.strip()])

def check_log_status(species):
    """Check latest log file status."""
    logs = sorted(log_dir.glob(f'workflow_{species}_*.log'), 
                  key=lambda p: p.stat().st_mtime, reverse=True)
    if not logs:
        return None, None, None
    
    latest = logs[0]
    mtime = latest.stat().st_mtime
    hours_ago = (datetime.now().timestamp() - mtime) / 3600
    size_kb = latest.stat().st_size // 1024
    
    return latest.name, hours_ago, size_kb

def main():
    print("=" * 80)
    print("COMPREHENSIVE AMALGKIT WORKFLOW STATUS")
    print("=" * 80)
    print(f"\n⏰ {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    total_quantified = 0
    total_samples = 0
    total_downloaded = 0
    total_unquantified = 0
    
    for species_id, info in species_info.items():
        quantified = count_quantified(species_id)
        downloaded = count_downloaded(species_id)
        unquantified = get_unquantified_count(species_id)
        
        pct = (quantified / info['total']) * 100 if info['total'] > 0 else 0
        remaining = info['total'] - quantified
        
        log_name, hours_ago, log_size = check_log_status(species_id)
        
        print(f"{info['name']} ({species_id})")
        print(f"  📊 Progress: {quantified}/{info['total']} ({pct:.1f}%)")
        print(f"  📥 Downloaded: {downloaded}")
        if unquantified is not None:
            print(f"  ⚠️  Unquantified (per list): {unquantified}")
        print(f"  📋 Remaining: {remaining}")
        if log_name:
            status = "🟢 ACTIVE" if hours_ago < 6 else f"🔴 Inactive ({hours_ago:.1f}h ago)"
            print(f"  📝 Log: {log_name} ({log_size} KB) - {status}")
        else:
            print(f"  📝 Log: None found")
        print()
        
        total_quantified += quantified
        total_samples += info['total']
        total_downloaded += downloaded
        if unquantified is not None:
            total_unquantified += unquantified
    
    overall_pct = (total_quantified / total_samples) * 100 if total_samples > 0 else 0
    print("=" * 80)
    print("OVERALL SUMMARY")
    print("=" * 80)
    print(f"Total Samples: {total_samples}")
    print(f"Quantified: {total_quantified}/{total_samples} ({overall_pct:.1f}%)")
    print(f"Downloaded: {total_downloaded}")
    if total_unquantified > 0:
        print(f"Unquantified: {total_unquantified}")
    print(f"Remaining: {total_samples - total_quantified}")
    print("=" * 80)
    
    # Recommendation
    print("\n💡 RECOMMENDATIONS:")
    if total_quantified < total_samples * 0.9:
        print("  - Workflows should be running to complete remaining samples")
        print("  - Run: python3 scripts/rna/restart_all_workflows.py")
    else:
        print("  - Workflows are making excellent progress!")
    print("  - Monitor: python3 scripts/rna/monitor_comprehensive.py")
    print()

if __name__ == '__main__':
    main()

