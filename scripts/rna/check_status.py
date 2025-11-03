#!/usr/bin/env python3
"""Check current amalgkit workflow status."""
from pathlib import Path

repo_root = Path(__file__).parent.parent.parent
output_dir = repo_root / 'output' / 'amalgkit'

species_info = {
    'cfloridanus': {'total': 307, 'name': 'Camponotus floridanus'},
    'pbarbatus': {'total': 83, 'name': 'Pogonomyrmex barbatus'},
    'mpharaonis': {'total': 100, 'name': 'Monomorium pharaonis'},
    'sinvicta': {'total': 354, 'name': 'Solenopsis invicta'},
}

print("=" * 80)
print("AMALGKIT WORKFLOW STATUS")
print("=" * 80)

total_quantified = 0
total_downloaded = 0
total_unquantified = 0

for species_id, info in species_info.items():
    species_dir = output_dir / species_id
    
    # Count quantified samples
    quant_dir = species_dir / 'quant'
    quantified = 0
    if quant_dir.exists():
        quantified = len([d for d in quant_dir.iterdir() 
                         if d.is_dir() and (d / 'abundance.tsv').exists()])
    
    # Count downloaded samples (with data files)
    fastq_dir = species_dir / 'fastq'
    downloaded = 0
    if fastq_dir.exists():
        downloaded = len([d for d in fastq_dir.iterdir() 
                          if d.is_dir() and (
                              list(d.rglob('*.fastq.gz')) or 
                              list(d.rglob('*.sra'))
                          )])
    
    # Count unquantified
    unquantified_file = output_dir / f'{species_id}_unquantified.txt'
    unquantified = 0
    if unquantified_file.exists():
        unquantified = len([l.strip() for l in unquantified_file.read_text().splitlines() 
                           if l.strip()])
    
    # Progress percentage
    pct = (quantified / info['total']) * 100 if info['total'] > 0 else 0
    
    print(f"\n{info['name']} ({species_id})")
    print(f"  Quantified: {quantified}/{info['total']} ({pct:.1f}%)")
    print(f"  Downloaded (with files): {downloaded}")
    print(f"  Unquantified (per list): {unquantified}")
    
    total_quantified += quantified
    total_downloaded += downloaded
    total_unquantified += unquantified

total_samples = sum(info['total'] for info in species_info.values())
overall_pct = (total_quantified / total_samples) * 100 if total_samples > 0 else 0

print("\n" + "=" * 80)
print("OVERALL SUMMARY")
print("=" * 80)
print(f"Total Samples: {total_samples}")
print(f"Quantified: {total_quantified}/{total_samples} ({overall_pct:.1f}%)")
print(f"Downloaded (with files): {total_downloaded}")
print(f"Remaining Unquantified: {total_unquantified}")
print(f"Progress: {total_quantified}/{total_samples} samples processed")

# Check for workflow logs
print("\n" + "=" * 80)
print("RECENT WORKFLOW LOGS")
print("=" * 80)
log_dir = repo_root / 'output'
workflow_logs = sorted(log_dir.glob('workflow_*.log'), key=lambda p: p.stat().st_mtime, reverse=True)
for log in workflow_logs[:4]:
    size = log.stat().st_size
    mtime = log.stat().st_mtime
    from datetime import datetime
    mtime_str = datetime.fromtimestamp(mtime).strftime('%Y-%m-%d %H:%M:%S')
    print(f"  {log.name} ({size//1024} KB, {mtime_str})")

print("\n" + "=" * 80)

