#!/usr/bin/env python3
"""
Analyze sample status across all species.

Categorizes samples into:
- Quantified and deleted (FASTQ/SRA removed)
- Quantified but not deleted (still have FASTQ/SRA files)
- Currently downloading
- Failed download
- Undownloaded (in metadata but no files)
"""

import os
import sys
from pathlib import Path
from datetime import datetime
from collections import defaultdict
from glob import glob

# Ensure virtual environment is activated
def ensure_venv_activated():
    """Automatically activate virtual environment if it exists."""
    repo_root = Path(__file__).parent.parent.parent.resolve()
    venv_python = repo_root / ".venv" / "bin" / "python3"
    venv_dir = repo_root / ".venv"
    
    current_python = Path(sys.executable)
    
    try:
        current_python.relative_to(repo_root / ".venv")
        if "VIRTUAL_ENV" not in os.environ:
            os.environ["VIRTUAL_ENV"] = str(venv_dir)
            venv_bin = str(venv_dir / "bin")
            if venv_bin not in os.environ.get("PATH", ""):
                os.environ["PATH"] = f"{venv_bin}:{os.environ.get('PATH', '')}"
        return
    except ValueError:
        pass
    
    if venv_python.exists():
        new_env = os.environ.copy()
        new_env["VIRTUAL_ENV"] = str(venv_dir)
        venv_bin = str(venv_dir / "bin")
        new_env["PATH"] = f"{venv_bin}:{new_env.get('PATH', '')}"
        new_env.pop("PYTHONHOME", None)
        os.execve(str(venv_python), [str(venv_python)] + sys.argv, new_env)

ensure_venv_activated()

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.rna import analyze_species_status, check_active_downloads
from metainformant.core.logging import get_logger

logger = get_logger("analyze_status")


def main():
    """Analyze sample status across all species."""
    repo_root = Path(__file__).parent.parent.parent.resolve()
    config_dir = repo_root / "config" / "amalgkit"
    
    if not config_dir.exists() or not list(config_dir.glob("amalgkit_*.yaml")):
        config_dir = repo_root / "config"
    
    # Discover all config files
    config_pattern = str(config_dir / "amalgkit_*.yaml")
    config_files = sorted(glob(config_pattern))
    
    species_configs = []
    for config_file in config_files:
        path = Path(config_file)
        if "template" in path.stem.lower():
            continue
        
        species_code = path.stem.replace("amalgkit_", "")
        display_name = species_code.replace("_", " ").title()
        species_configs.append((display_name, path))
    
    if not species_configs:
        logger.error("⚠️  No species configs found")
        return 1
    
    print("\n" + "=" * 80)
    print("SAMPLE STATUS ANALYSIS")
    print("=" * 80)
    print(f"Date: {datetime.now()}")
    print(f"Analyzing {len(species_configs)} species")
    print("=" * 80 + "\n")
    
    # Check for active downloads
    print("Checking for active downloads...")
    active_downloads = check_active_downloads()
    if active_downloads:
        print(f"  Found {len(active_downloads)} active downloads")
    else:
        print("  No active downloads detected")
    print()
    
    # Analyze each species using metainformant functions
    all_results = []
    totals = defaultdict(int)
    
    for species_name, config_path in species_configs:
        if not config_path.exists():
            continue
        
        # Use metainformant function
        status = analyze_species_status(config_path)
        if status and status.get("total_in_metadata", 0) > 0:
            result = {
                "species": species_name,
                "total_in_metadata": status["total_in_metadata"],
                "categories": status["categories"],
            }
            all_results.append(result)
            totals["total_in_metadata"] += result["total_in_metadata"]
            for category, samples in result["categories"].items():
                totals[category] += len(samples)
    
    # Print detailed results
    print("=" * 80)
    print("DETAILED RESULTS BY SPECIES")
    print("=" * 80)
    print()
    
    for result in all_results:
        print(f"\n{result['species']}:")
        print(f"  Total in metadata: {result['total_in_metadata']}")
        
        cats = result["categories"]
        print(f"  ✅ Quantified and deleted: {len(cats['quantified_and_deleted'])}")
        print(f"  ⚠️  Quantified but not deleted: {len(cats['quantified_not_deleted'])}")
        print(f"  ⬇️  Currently downloading: {len(cats['downloading'])}")
        print(f"  ❌ Failed download: {len(cats['failed_download'])}")
        print(f"  ⏳ Undownloaded: {len(cats['undownloaded'])}")
        
        # Show some examples
        if cats["quantified_not_deleted"]:
            print(f"    Examples (quantified but not deleted): {', '.join(cats['quantified_not_deleted'][:5])}")
        if cats["failed_download"]:
            print(f"    Examples (failed): {', '.join(cats['failed_download'][:5])}")
        if cats["downloading"]:
            print(f"    Examples (downloading): {', '.join(cats['downloading'][:5])}")
    
    # Print summary
    print("\n" + "=" * 80)
    print("SUMMARY STATISTICS")
    print("=" * 80)
    print(f"Total species analyzed: {len(all_results)}")
    print(f"Total samples in metadata: {totals['total_in_metadata']}")
    print()
    print("By Category:")
    print(f"  ✅ Quantified and deleted: {totals['quantified_and_deleted']} ({totals['quantified_and_deleted']/max(totals['total_in_metadata'], 1)*100:.1f}%)")
    print(f"  ⚠️  Quantified but not deleted: {totals['quantified_not_deleted']} ({totals['quantified_not_deleted']/max(totals['total_in_metadata'], 1)*100:.1f}%)")
    print(f"  ⬇️  Currently downloading: {totals['downloading']} ({totals['downloading']/max(totals['total_in_metadata'], 1)*100:.1f}%)")
    print(f"  ❌ Failed download: {totals['failed_download']} ({totals['failed_download']/max(totals['total_in_metadata'], 1)*100:.1f}%)")
    print(f"  ⏳ Undownloaded: {totals['undownloaded']} ({totals['undownloaded']/max(totals['total_in_metadata'], 1)*100:.1f}%)")
    print()
    
    # Calculate completion stats
    total_quantified = totals['quantified_and_deleted'] + totals['quantified_not_deleted']
    print(f"Quantification Progress:")
    print(f"  Total quantified: {total_quantified} ({total_quantified/max(totals['total_in_metadata'], 1)*100:.1f}%)")
    print(f"  Cleanup needed: {totals['quantified_not_deleted']} samples still have FASTQ/SRA files")
    print()
    
    print(f"Download Progress:")
    total_downloaded = total_quantified + totals['failed_download']
    print(f"  Total with files: {total_downloaded} ({total_downloaded/max(totals['total_in_metadata'], 1)*100:.1f}%)")
    print(f"  Remaining: {totals['undownloaded']} samples need downloading")
    print()
    
    print("=" * 80)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

