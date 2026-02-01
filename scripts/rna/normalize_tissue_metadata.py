#!/usr/bin/env python3
"""Normalize tissue metadata values in amalgkit metadata files.

This script applies tissue normalization to metadata tables, creating
standardized tissue values while preserving originals.

Usage:
    uv run python scripts/rna/normalize_tissue_metadata.py \\
        --input output/amalgkit/apis_mellifera_all/work/metadata/metadata.tsv \\
        --mapping config/amalgkit/tissue_mapping.yaml \\
        --patches config/amalgkit/tissue_patches.yaml \\
        --output output/amalgkit/apis_mellifera_all/work/metadata/metadata_normalized.tsv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Add src to path for local development
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

import pandas as pd

from metainformant.rna.amalgkit.tissue_normalizer import (
    apply_tissue_normalization,
    get_unmapped_tissues,
    get_missing_tissue_samples,
    get_canonical_tissues,
)


def main():
    parser = argparse.ArgumentParser(
        description="Normalize tissue metadata values to canonical forms"
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Path to input metadata TSV file"
    )
    parser.add_argument(
        "--mapping", "-m",
        default="config/amalgkit/tissue_mapping.yaml",
        help="Path to tissue mapping YAML (default: config/amalgkit/tissue_mapping.yaml)"
    )
    parser.add_argument(
        "--patches", "-p",
        default="config/amalgkit/tissue_patches.yaml",
        help="Path to tissue patches YAML (default: config/amalgkit/tissue_patches.yaml)"
    )
    parser.add_argument(
        "--output", "-o",
        help="Path for output TSV (default: input with _normalized suffix)"
    )
    parser.add_argument(
        "--report-unmapped", "-u",
        action="store_true",
        help="Print report of unmapped tissue values"
    )
    parser.add_argument(
        "--report-missing", "-s",
        action="store_true",
        help="Print report of samples with missing tissue"
    )
    parser.add_argument(
        "--list-canonical",
        action="store_true",
        help="List canonical tissue names and exit"
    )
    parser.add_argument(
        "--dry-run", "-n",
        action="store_true",
        help="Don't write output file, just report statistics"
    )
    
    args = parser.parse_args()
    
    # Handle list-canonical request
    if args.list_canonical:
        print("Canonical tissue names:")
        for tissue in get_canonical_tissues(args.mapping):
            print(f"  - {tissue}")
        return 0
    
    # Check input file
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        return 1
    
    # Determine output path
    output_path = args.output
    if not output_path:
        output_path = input_path.parent / f"{input_path.stem}_normalized{input_path.suffix}"
    
    # Load metadata
    print(f"Loading metadata from: {args.input}")
    df = pd.read_csv(input_path, sep="\t", low_memory=False)
    print(f"  Loaded {len(df)} samples")
    
    # Check for tissue column
    if "tissue" not in df.columns:
        print("Error: No 'tissue' column found in metadata", file=sys.stderr)
        print(f"  Available columns: {list(df.columns)[:20]}...", file=sys.stderr)
        return 1
    
    # Apply normalization
    print(f"\nApplying tissue normalization...")
    print(f"  Mapping: {args.mapping}")
    print(f"  Patches: {args.patches}")
    
    df_normalized = apply_tissue_normalization(
        df,
        mapping_path=args.mapping,
        patches_path=args.patches if Path(args.patches).exists() else None
    )
    
    # Report unmapped values
    if args.report_unmapped:
        print("\n" + "=" * 60)
        print("UNMAPPED TISSUE VALUES (add to tissue_mapping.yaml)")
        print("=" * 60)
        unmapped = get_unmapped_tissues(df, args.mapping)
        if unmapped:
            for value, count in list(unmapped.items())[:50]:
                print(f"  {count:5d}  {value}")
            if len(unmapped) > 50:
                print(f"  ... and {len(unmapped) - 50} more")
        else:
            print("  All tissue values are mapped!")
    
    # Report missing samples
    if args.report_missing:
        print("\n" + "=" * 60)
        print("SAMPLES WITH MISSING TISSUE (add to tissue_patches.yaml)")
        print("=" * 60)
        missing = get_missing_tissue_samples(df)
        if not missing.empty:
            # Group by bioproject
            by_project = missing.groupby("bioproject").size().sort_values(ascending=False)
            print(f"  {len(missing)} samples missing tissue across {len(by_project)} bioprojects:")
            for project, count in by_project.head(20).items():
                print(f"    {project}: {count} samples")
            if len(by_project) > 20:
                print(f"    ... and {len(by_project) - 20} more projects")
        else:
            print("  All samples have tissue information!")
    
    # Summary statistics
    print("\n" + "=" * 60)
    print("NORMALIZATION SUMMARY")
    print("=" * 60)
    
    normalized_counts = df_normalized["tissue_normalized"].value_counts()
    total = len(df_normalized)
    mapped = (df_normalized["tissue_normalized"] != "").sum()
    unmapped_count = total - mapped
    
    print(f"  Total samples: {total}")
    print(f"  Successfully normalized: {mapped} ({100*mapped/total:.1f}%)")
    print(f"  Unmapped/empty: {unmapped_count} ({100*unmapped_count/total:.1f}%)")
    print("\n  Top normalized tissues:")
    for tissue, count in normalized_counts.head(15).items():
        if tissue:
            print(f"    {tissue}: {count}")
    
    # Write output
    if not args.dry_run:
        print(f"\nWriting normalized metadata to: {output_path}")
        df_normalized.to_csv(output_path, sep="\t", index=False)
        print("Done!")
    else:
        print("\n[DRY RUN] No output file written")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
