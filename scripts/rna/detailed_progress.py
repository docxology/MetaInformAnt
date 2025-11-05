#!/usr/bin/env python3
"""Detailed progress check with per-species breakdown."""
from pathlib import Path
from datetime import datetime

repo_root = Path(__file__).parent.parent.parent
output_dir = repo_root / 'output' / 'amalgkit'

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

def count_downloaded_with_files(species):
    """Count downloaded samples with actual data files."""
    fastq_dir = output_dir / species / 'fastq'
    if not fastq_dir.exists():
        return 0
    
    downloaded = []
    for sample_dir in fastq_dir.iterdir():
        if not sample_dir.is_dir():
            continue
        
        # Check for FASTQ or SRA files
        fastq_files = list(sample_dir.rglob('*.fastq.gz'))
        sra_files = list(sample_dir.rglob('*.sra'))
        
        if fastq_files or sra_files:
            # Check if quantified
            quant_file = output_dir / species / 'quant' / sample_dir.name / 'abundance.tsv'
            if not quant_file.exists():
                downloaded.append(sample_dir.name)
    
    return len(downloaded)

def get_unquantified_list(species):
    """Get list of unquantified samples."""
    unquant_file = output_dir / f'{species}_unquantified.txt'
    if not unquant_file.exists():
        return []
    return [l.strip() for l in unquant_file.read_text().splitlines() if l.strip()]

def get_disk_usage(path):
    """Get disk usage in GB."""
    try:
        import subprocess
        result = subprocess.run(
            ['du', '-sh', str(path)],
            capture_output=True,
            text=True,
            timeout=5
        )
        if result.returncode == 0:
            size_str = result.stdout.split()[0]
            # Convert to GB if needed
            if 'G' in size_str:
                return float(size_str.replace('G', '')), size_str
            elif 'M' in size_str:
                return float(size_str.replace('M', '')) / 1024, size_str
            elif 'K' in size_str:
                return float(size_str.replace('K', '')) / (1024 * 1024), size_str
    except:
        pass
    return None, None

def main():
    print("=" * 90)
    print("DETAILED AMALGKIT PROGRESS REPORT")
    print("=" * 90)
    print(f"\n⏰ {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    total_quantified = 0
    total_samples = 0
    total_unquantified = 0
    total_downloaded = 0
    
    for species_id, info in species_info.items():
        quantified = count_quantified(species_id)
        downloaded = count_downloaded_with_files(species_id)
        unquantified = get_unquantified_list(species_id)
        
        pct = (quantified / info['total']) * 100 if info['total'] > 0 else 0
        remaining = info['total'] - quantified
        progress_bar_length = 50
        filled = int(progress_bar_length * quantified / info['total']) if info['total'] > 0 else 0
        bar = '█' * filled + '░' * (progress_bar_length - filled)
        
        # Get disk usage
        fastq_dir = output_dir / species_id / 'fastq'
        fastq_size, fastq_str = get_disk_usage(fastq_dir) if fastq_dir.exists() else (None, None)
        
        quant_dir = output_dir / species_id / 'quant'
        quant_size, quant_str = get_disk_usage(quant_dir) if quant_dir.exists() else (None, None)
        
        print(f"{info['name']} ({species_id})")
        print(f"  Progress: [{bar}] {pct:.1f}%")
        print(f"  Quantified: {quantified}/{info['total']}")
        print(f"  Remaining: {remaining}")
        print(f"  Downloaded (not quantified): {downloaded}")
        print(f"  Unquantified list: {len(unquantified)} samples")
        if fastq_str:
            print(f"  FASTQ size: {fastq_str}")
        if quant_str:
            print(f"  Quant size: {quant_str}")
        
        if len(unquantified) > 0 and len(unquantified) <= 10:
            print(f"  Unquantified samples: {', '.join(unquantified[:10])}")
        elif len(unquantified) > 10:
            print(f"  Unquantified samples: {', '.join(unquantified[:5])} ... and {len(unquantified) - 5} more")
        
        print()
        
        total_quantified += quantified
        total_samples += info['total']
        total_unquantified += len(unquantified)
        total_downloaded += downloaded
    
    overall_pct = (total_quantified / total_samples) * 100 if total_samples > 0 else 0
    progress_bar_length = 50
    filled = int(progress_bar_length * total_quantified / total_samples) if total_samples > 0 else 0
    bar = '█' * filled + '░' * (progress_bar_length - filled)
    
    print("=" * 90)
    print("OVERALL SUMMARY")
    print("=" * 90)
    print(f"Progress: [{bar}] {overall_pct:.1f}%")
    print(f"Total Samples: {total_samples}")
    print(f"Quantified: {total_quantified}/{total_samples}")
    print(f"Remaining: {total_samples - total_quantified}")
    print(f"Downloaded (not quantified): {total_downloaded}")
    print(f"Unquantified list: {total_unquantified}")
    
    # Calculate estimated time
    remaining = total_samples - total_quantified
    if remaining > 0:
        # Estimate: ~7.5 minutes per sample (download + quant)
        est_minutes = remaining * 7.5
        est_hours = est_minutes / 60
        est_days = est_hours / 24
        print(f"\n⏳ Estimated time to complete:")
        if est_days >= 1:
            print(f"   {est_days:.1f} days ({est_hours:.1f} hours)")
        else:
            print(f"   {est_hours:.1f} hours ({est_minutes:.0f} minutes)")
    
    print("=" * 90)

if __name__ == '__main__':
    main()






