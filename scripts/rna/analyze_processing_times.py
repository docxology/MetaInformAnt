#!/usr/bin/env python3
"""
Analyze Apis mellifera processing times vs sample sizes.

Parses the parallel_processing.log file to extract timing data
and creates scatterplots with regression statistics.

Usage:
    python analyze_processing_times.py
"""

import re
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from datetime import datetime

# Configuration
LOG_FILE = Path("output/amalgkit/apis_mellifera_all/work/parallel_processing.log")
OUTPUT_DIR = Path("output/amalgkit/apis_mellifera_all/work/analysis")


def parse_log_file(log_path: Path) -> list[dict]:
    """Parse log file and extract timing data.
    
    Expected format (new):
    [56/4254] ✓ SRR5008299 (96.1s, 22.1MB) [dl:45s ex:40s qt:10s]
    
    Expected format (legacy):
    [56/4254] ✓ SRR5008299 (96.1s, 22.1MB)
    """
    data = []
    
    # Pattern with per-stage timing
    pattern_full = re.compile(
        r'\[(\d+)/\d+\]\s+✓\s+(\S+)\s+'
        r'\((\d+\.?\d*)s,\s+(\d+\.?\d*)MB\)\s*'
        r'\[dl:(\d+)s\s+ex:(\d+)s\s+qt:(\d+)s\]'
    )
    
    # Pattern without per-stage timing (legacy)
    pattern_simple = re.compile(
        r'\[(\d+)/\d+\]\s+✓\s+(\S+)\s+'
        r'\((\d+\.?\d*)s,\s+(\d+\.?\d*)MB\)'
    )
    
    with open(log_path, 'r') as f:
        for line in f:
            # Try full pattern first
            match = pattern_full.search(line)
            if match:
                data.append({
                    'index': int(match.group(1)),
                    'sample_id': match.group(2),
                    'total_time': float(match.group(3)),
                    'size_mb': float(match.group(4)),
                    'download_time': float(match.group(5)),
                    'extract_time': float(match.group(6)),
                    'quant_time': float(match.group(7)),
                })
                continue
            
            # Try simple pattern
            match = pattern_simple.search(line)
            if match:
                data.append({
                    'index': int(match.group(1)),
                    'sample_id': match.group(2),
                    'total_time': float(match.group(3)),
                    'size_mb': float(match.group(4)),
                    'download_time': None,
                    'extract_time': None,
                    'quant_time': None,
                })
    
    return data


def compute_regression(x: np.ndarray, y: np.ndarray) -> dict:
    """Compute linear regression statistics."""
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    
    return {
        'slope': slope,
        'intercept': intercept,
        'r_squared': r_value ** 2,
        'r_value': r_value,
        'p_value': p_value,
        'std_err': std_err,
    }


def create_scatterplot(data: list[dict], output_dir: Path) -> None:
    """Create scatterplot of size vs total processing time."""
    output_dir.mkdir(parents=True, exist_ok=True)
    
    sizes = np.array([d['size_mb'] for d in data])
    times = np.array([d['total_time'] for d in data])
    
    # Compute regression
    reg = compute_regression(sizes, times)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Scatter plot
    ax.scatter(sizes, times, alpha=0.5, s=20, label='Samples')
    
    # Regression line
    x_line = np.linspace(sizes.min(), sizes.max(), 100)
    y_line = reg['slope'] * x_line + reg['intercept']
    ax.plot(x_line, y_line, 'r-', linewidth=2, label='Regression line')
    
    # Labels and title
    ax.set_xlabel('Sample Size (MB)', fontsize=12)
    ax.set_ylabel('Processing Time (seconds)', fontsize=12)
    ax.set_title('Apis mellifera Sample Processing: Size vs Time', fontsize=14)
    
    # Statistics text box
    stats_text = (
        f"n = {len(data)}\n"
        f"R² = {reg['r_squared']:.4f}\n"
        f"slope = {reg['slope']:.3f} s/MB\n"
        f"intercept = {reg['intercept']:.1f} s\n"
        f"p-value = {reg['p_value']:.2e}"
    )
    ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.3)
    
    # Save
    output_path = output_dir / 'size_vs_time_scatter.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")
    
    return reg


def create_stage_breakdown_plot(data: list[dict], output_dir: Path) -> None:
    """Create stacked bar chart showing time breakdown by stage."""
    # Filter to only samples with per-stage timing
    staged_data = [d for d in data if d['download_time'] is not None]
    
    if len(staged_data) < 10:
        print(f"Not enough samples with per-stage timing ({len(staged_data)} < 10)")
        return
    
    # Sort by size
    staged_data.sort(key=lambda x: x['size_mb'])
    
    # Group into size bins
    sizes = np.array([d['size_mb'] for d in staged_data])
    bin_edges = np.percentile(sizes, [0, 25, 50, 75, 100])
    
    bin_data = []
    bin_labels = []
    
    for i in range(len(bin_edges) - 1):
        mask = (sizes >= bin_edges[i]) & (sizes < bin_edges[i + 1])
        if i == len(bin_edges) - 2:
            mask = (sizes >= bin_edges[i]) & (sizes <= bin_edges[i + 1])
        
        bin_samples = [d for d, m in zip(staged_data, mask) if m]
        if bin_samples:
            avg_dl = np.mean([d['download_time'] for d in bin_samples])
            avg_ex = np.mean([d['extract_time'] for d in bin_samples])
            avg_qt = np.mean([d['quant_time'] for d in bin_samples])
            bin_data.append({'download': avg_dl, 'extract': avg_ex, 'quant': avg_qt})
            bin_labels.append(f'{bin_edges[i]:.0f}-{bin_edges[i+1]:.0f} MB')
    
    if not bin_data:
        return
    
    # Create stacked bar chart
    fig, ax = plt.subplots(figsize=(10, 6))
    
    x = np.arange(len(bin_labels))
    width = 0.6
    
    dl_times = [d['download'] for d in bin_data]
    ex_times = [d['extract'] for d in bin_data]
    qt_times = [d['quant'] for d in bin_data]
    
    ax.bar(x, dl_times, width, label='Download', color='#3498db')
    ax.bar(x, ex_times, width, bottom=dl_times, label='Extract', color='#e74c3c')
    ax.bar(x, qt_times, width, bottom=np.array(dl_times) + np.array(ex_times), 
           label='Quantify', color='#2ecc71')
    
    ax.set_xlabel('Size Bin', fontsize=12)
    ax.set_ylabel('Average Time (seconds)', fontsize=12)
    ax.set_title('Processing Time Breakdown by Size Bin', fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(bin_labels, rotation=45)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    output_path = output_dir / 'time_breakdown_by_size.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def print_summary_statistics(data: list[dict], reg: dict) -> None:
    """Print summary statistics to console."""
    sizes = np.array([d['size_mb'] for d in data])
    times = np.array([d['total_time'] for d in data])
    
    print("\n" + "="*60)
    print("PROCESSING TIME ANALYSIS SUMMARY")
    print("="*60)
    print(f"\nSamples analyzed: {len(data)}")
    print(f"\nSize statistics:")
    print(f"  Min:    {sizes.min():.1f} MB")
    print(f"  Max:    {sizes.max():.1f} MB")
    print(f"  Mean:   {sizes.mean():.1f} MB")
    print(f"  Median: {np.median(sizes):.1f} MB")
    print(f"\nTime statistics:")
    print(f"  Min:    {times.min():.1f} s")
    print(f"  Max:    {times.max():.1f} s")
    print(f"  Mean:   {times.mean():.1f} s")
    print(f"  Median: {np.median(times):.1f} s")
    print(f"\nRegression analysis:")
    print(f"  y = {reg['slope']:.3f}x + {reg['intercept']:.1f}")
    print(f"  R² = {reg['r_squared']:.4f}")
    print(f"  p-value = {reg['p_value']:.2e}")
    print(f"\nThroughput:")
    throughput = sizes.sum() / times.sum()
    print(f"  Average: {throughput:.2f} MB/s = {throughput*60:.1f} MB/min")
    print("="*60 + "\n")


def main():
    """Main analysis function."""
    print(f"Parsing log file: {LOG_FILE}")
    
    if not LOG_FILE.exists():
        print(f"Error: Log file not found: {LOG_FILE}")
        return
    
    data = parse_log_file(LOG_FILE)
    
    if len(data) < 5:
        print(f"Not enough data points ({len(data)}). Need at least 5 successful samples.")
        return
    
    print(f"Found {len(data)} successful samples with timing data")
    
    # Create visualizations
    reg = create_scatterplot(data, OUTPUT_DIR)
    create_stage_breakdown_plot(data, OUTPUT_DIR)
    
    # Print summary
    print_summary_statistics(data, reg)
    
    print(f"\nOutput saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
