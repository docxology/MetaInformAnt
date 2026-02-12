#!/usr/bin/env python3
"""Prepare size-sorted, size-capped metadata for full-download mode.

Reads a metadata TSV, sorts by total_bases ascending (smallest first),
filters out samples exceeding a max size, and writes the result.

Usage:
    python3 scripts/rna/prepare_size_sorted_metadata.py \
        --input blue/amalgkit/amellifera/work/metadata/metadata_selected.tsv \
        --output blue/amalgkit/amellifera/work/metadata/metadata.tsv \
        --max-gb 5.0
"""

import argparse
import pandas as pd
from pathlib import Path


def prepare_metadata(input_path: Path, output_path: Path, max_gb: float) -> dict:
    """Sort metadata by total_bases ascending, filter by max size."""
    df = pd.read_csv(input_path, sep='\t')

    original_count = len(df)

    # Convert total_bases to numeric, coercing errors
    df['total_bases'] = pd.to_numeric(df['total_bases'], errors='coerce')

    # Drop rows with no total_bases
    df = df.dropna(subset=['total_bases'])

    # Filter by max size (convert GB to bases)
    max_bases = max_gb * 1e9
    df_filtered = df[df['total_bases'] <= max_bases].copy()

    # Sort ascending by total_bases (smallest first)
    df_sorted = df_filtered.sort_values('total_bases', ascending=True).reset_index(drop=True)

    # Write output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df_sorted.to_csv(output_path, sep='\t', index=False)

    stats = {
        'original': original_count,
        'after_filter': len(df_sorted),
        'removed': original_count - len(df_sorted),
        'max_gb': max_gb,
        'smallest_gb': df_sorted['total_bases'].min() / 1e9 if len(df_sorted) > 0 else 0,
        'largest_gb': df_sorted['total_bases'].max() / 1e9 if len(df_sorted) > 0 else 0,
        'total_gb': df_sorted['total_bases'].sum() / 1e9 if len(df_sorted) > 0 else 0,
    }
    return stats


def main():
    parser = argparse.ArgumentParser(description='Prepare size-sorted metadata')
    parser.add_argument('--input', type=Path, required=True, help='Input metadata TSV')
    parser.add_argument('--output', type=Path, required=True, help='Output metadata TSV')
    parser.add_argument('--max-gb', type=float, default=5.0, help='Max sample size in GB')
    args = parser.parse_args()

    stats = prepare_metadata(args.input, args.output, args.max_gb)

    print(f"Size-Sorted Metadata Prepared:")
    print(f"  Original samples:  {stats['original']}")
    print(f"  After filter (≤{stats['max_gb']}GB): {stats['after_filter']}")
    print(f"  Removed:           {stats['removed']}")
    print(f"  Smallest sample:   {stats['smallest_gb']:.2f} GB")
    print(f"  Largest sample:    {stats['largest_gb']:.2f} GB")
    print(f"  Total download:    {stats['total_gb']:.0f} GB")


if __name__ == '__main__':
    main()
