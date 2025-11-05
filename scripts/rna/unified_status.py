#!/usr/bin/env python3
"""Unified status checking script for all RNA-seq workflows.

Replaces: check_status.py, comprehensive_status.py, detailed_progress.py,
          full_assessment.py, get_current_status.py, quick_status.py

Usage:
    python3 scripts/rna/unified_status.py [--species SPECIES] [--detailed]
"""

import sys
from pathlib import Path
from datetime import datetime
from glob import glob
from collections import defaultdict

# Add repo root to path
repo_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(repo_root / "src"))

from metainformant.rna import (
    analyze_species_status,
    check_workflow_progress,
    count_quantified_samples,
)


def discover_species_configs(config_dir: Path = None) -> list[tuple[str, Path]]:
    """Discover all species configuration files."""
    if config_dir is None:
        config_dir = repo_root / "config" / "amalgkit"
    
    if not config_dir.exists():
        config_dir = repo_root / "config"
    
    config_pattern = str(config_dir / "amalgkit_*.yaml")
    config_files = sorted(glob(config_pattern))
    
    species_configs = []
    for config_file in config_files:
        path = Path(config_file)
        if "template" in path.stem.lower():
            continue
        
        species_name = path.stem.replace("amalgkit_", "").replace("_", " ").title()
        species_configs.append((species_name, path))
    
    return species_configs


def print_brief_status(species_configs: list[tuple[str, Path]], species_filter: str = None):
    """Print brief status summary."""
    print("=" * 80)
    print("RNA-SEQ WORKFLOW STATUS")
    print("=" * 80)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    total_quantified = 0
    total_samples = 0
    
    for species_name, config_path in species_configs:
        if species_filter and species_filter.lower() not in species_name.lower():
            continue
        
        if not config_path.exists():
            continue
        
        progress = check_workflow_progress(config_path)
        total_quantified += progress["quantified"]
        total_samples += progress["total"]
        
        pct = progress["percentage"]
        print(f"{species_name:<30} {progress['quantified']:>4}/{progress['total']:<4} ({pct:>5.1f}%)")
    
    print("-" * 80)
    overall_pct = (total_quantified / total_samples * 100) if total_samples > 0 else 0.0
    print(f"{'TOTAL':<30} {total_quantified:>4}/{total_samples:<4} ({overall_pct:>5.1f}%)")
    print("=" * 80)


def print_detailed_status(species_configs: list[tuple[str, Path]], species_filter: str = None):
    """Print detailed status with sample categories."""
    print("=" * 80)
    print("COMPREHENSIVE RNA-SEQ WORKFLOW STATUS")
    print("=" * 80)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    all_totals = defaultdict(int)
    
    for species_name, config_path in species_configs:
        if species_filter and species_filter.lower() not in species_name.lower():
            continue
        
        if not config_path.exists():
            continue
        
        status = analyze_species_status(config_path)
        if not status or status["total_in_metadata"] == 0:
            continue
        
        categories = status["categories"]
        
        print(f"\n{species_name}:")
        print(f"  Total in metadata: {status['total_in_metadata']}")
        print(f"  ✅ Quantified and deleted: {len(categories['quantified_and_deleted'])}")
        print(f"  ⚠️  Quantified but not deleted: {len(categories['quantified_not_deleted'])}")
        print(f"  ⬇️  Currently downloading: {len(categories['downloading'])}")
        print(f"  ❌ Failed download: {len(categories['failed_download'])}")
        print(f"  ⏳ Undownloaded: {len(categories['undownloaded'])}")
        
        # Update totals
        all_totals["total_in_metadata"] += status["total_in_metadata"]
        all_totals["quantified_and_deleted"] += len(categories["quantified_and_deleted"])
        all_totals["quantified_not_deleted"] += len(categories["quantified_not_deleted"])
        all_totals["downloading"] += len(categories["downloading"])
        all_totals["failed_download"] += len(categories["failed_download"])
        all_totals["undownloaded"] += len(categories["undownloaded"])
    
    # Print summary
    print("\n" + "=" * 80)
    print("OVERALL SUMMARY")
    print("=" * 80)
    total = all_totals["total_in_metadata"]
    if total > 0:
        print(f"Total samples: {total}")
        print(f"✅ Quantified and deleted: {all_totals['quantified_and_deleted']} ({all_totals['quantified_and_deleted']/total*100:.1f}%)")
        print(f"⚠️  Quantified but not deleted: {all_totals['quantified_not_deleted']} ({all_totals['quantified_not_deleted']/total*100:.1f}%)")
        print(f"⬇️  Currently downloading: {all_totals['downloading']} ({all_totals['downloading']/total*100:.1f}%)")
        print(f"❌ Failed download: {all_totals['failed_download']} ({all_totals['failed_download']/total*100:.1f}%)")
        print(f"⏳ Undownloaded: {all_totals['undownloaded']} ({all_totals['undownloaded']/total*100:.1f}%)")
        
        total_quantified = all_totals["quantified_and_deleted"] + all_totals["quantified_not_deleted"]
        print(f"\nQuantification Progress: {total_quantified}/{total} ({total_quantified/total*100:.1f}%)")
    print("=" * 80)


def main():
    """Main entry point."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Check status of RNA-seq workflows")
    parser.add_argument("--species", help="Filter by species name (partial match)")
    parser.add_argument("--detailed", action="store_true", help="Show detailed status with categories")
    args = parser.parse_args()
    
    species_configs = discover_species_configs()
    
    if not species_configs:
        print("❌ No species configurations found!")
        print("   Looked in: config/amalgkit/amalgkit_*.yaml")
        return 1
    
    if args.detailed:
        print_detailed_status(species_configs, args.species)
    else:
        print_brief_status(species_configs, args.species)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())


