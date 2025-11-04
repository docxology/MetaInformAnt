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

from metainformant.rna.workflow import load_workflow_config
from metainformant.core.io import read_delimited
from metainformant.core.logging import get_logger

logger = get_logger("analyze_status")


def check_active_downloads() -> set:
    """Check for samples currently being downloaded."""
    active_samples = set()
    
    # Check running amalgkit getfastq processes
    try:
        import subprocess
        result = subprocess.run(
            ["ps", "aux"],
            capture_output=True,
            text=True,
            timeout=5
        )
        import re
        for line in result.stdout.split("\n"):
            if "amalgkit getfastq" in line:
                # Extract all SRR IDs from the line
                matches = re.findall(r"SRR\d+", line)
                active_samples.update(matches)
                
                # Also check for --id parameter
                if "--id" in line:
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if part == "--id" and i + 1 < len(parts):
                            active_samples.add(parts[i + 1])
    except Exception:
        pass
    
    # Also check for recent fastq directory activity
    try:
        repo_root = Path(__file__).parent.parent.parent.resolve()
        fastq_base = repo_root / "output" / "amalgkit"
        
        # Find directories modified in last 10 minutes
        for species_dir in fastq_base.glob("*/fastq"):
            if species_dir.exists():
                # Check for recently modified sample directories
                for sample_dir in species_dir.glob("getfastq/SRR*"):
                    if sample_dir.exists():
                        # Check if directory was modified recently
                        import time
                        mtime = sample_dir.stat().st_mtime
                        if time.time() - mtime < 600:  # 10 minutes
                            sample_id = sample_dir.name
                            active_samples.add(sample_id)
    except Exception:
        pass
    
    return active_samples


def analyze_species(species_name: str, config_path: Path, active_downloads: set) -> dict:
    """Analyze sample status for a single species."""
    try:
        cfg = load_workflow_config(config_path)
    except Exception as e:
        logger.warning(f"Failed to load config for {species_name}: {e}")
        return {}
    
    fastq_dir = Path(cfg.per_step.get("getfastq", {}).get("out_dir", cfg.work_dir / "fastq"))
    quant_dir = Path(cfg.per_step.get("quant", {}).get("out_dir", cfg.work_dir / "quant"))
    
    # Get metadata
    metadata_file = cfg.work_dir / "metadata" / "metadata.tsv"
    if not metadata_file.exists():
        metadata_file = cfg.work_dir / "metadata" / "metadata.filtered.tissue.tsv"
    
    if not metadata_file.exists():
        return {}
    
    try:
        rows = list(read_delimited(metadata_file, delimiter="\t"))
        sample_ids = [row.get("run") for row in rows if row.get("run")]
    except Exception as e:
        logger.warning(f"Failed to read metadata for {species_name}: {e}")
        return {}
    
    # Categorize samples
    categories = {
        "quantified_and_deleted": [],
        "quantified_not_deleted": [],
        "downloading": [],
        "failed_download": [],
        "undownloaded": [],
    }
    
    for sample_id in sample_ids:
        if not sample_id:
            continue
        
        # Check if quantified
        abundance_file = quant_dir / sample_id / "abundance.tsv"
        is_quantified = abundance_file.exists()
        
        # Check for FASTQ/SRA files
        has_fastq = False
        has_sra = False
        
        # Check getfastq subdirectory
        sample_dir_getfastq = fastq_dir / "getfastq" / sample_id
        if sample_dir_getfastq.exists():
            has_fastq = any(sample_dir_getfastq.glob("*.fastq*"))
            has_sra = any(sample_dir_getfastq.glob("*.sra"))
        
        # Check direct structure
        sample_dir_direct = fastq_dir / sample_id
        if sample_dir_direct.exists():
            has_fastq = has_fastq or any(sample_dir_direct.glob("*.fastq*"))
            has_sra = has_sra or any(sample_dir_direct.glob("*.sra"))
        
        has_files = has_fastq or has_sra
        
        # Categorize
        if sample_id in active_downloads:
            categories["downloading"].append(sample_id)
        elif is_quantified and not has_files:
            categories["quantified_and_deleted"].append(sample_id)
        elif is_quantified and has_files:
            categories["quantified_not_deleted"].append(sample_id)
        elif has_files:
            # Has files but not quantified - might be failed or in progress
            # Check if there are error logs or if it's very recent
            categories["failed_download"].append(sample_id)
        else:
            categories["undownloaded"].append(sample_id)
    
    return {
        "species": species_name,
        "total_in_metadata": len(sample_ids),
        "categories": categories,
    }


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
    
    # Analyze each species
    all_results = []
    totals = defaultdict(int)
    
    for species_name, config_path in species_configs:
        if not config_path.exists():
            continue
        
        result = analyze_species(species_name, config_path, active_downloads)
        if result:
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

