#!/usr/bin/env python3
"""Download reference transcriptomes and build kallisto indices for all species.

Reads genome.ftp_url + genome.files.transcriptome_fasta from each species
YAML config, downloads from NCBI FTP, then runs `kallisto index`.

Usage:
    python3 scripts/cloud/prep_genomes.py [--config-dir config/amalgkit] [--threads 8]
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

import yaml


def process_species(config_path: Path, threads: int) -> bool:
    """Download genome + build index for one species."""
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    species = config_path.stem.replace("amalgkit_", "")

    # Get genome config
    genome = cfg.get("genome")
    if not genome:
        print(f"  ⚠ No genome config in {config_path.name}, skipping")
        return False

    ftp_url = genome.get("ftp_url", "")
    files = genome.get("files", {})
    transcriptome = files.get("transcriptome_fasta", "")
    dest_dir = Path(genome.get("dest_dir", f"output/amalgkit/shared/genome/{species}"))

    if not ftp_url or not transcriptome:
        print(f"  ⚠ Missing ftp_url or transcriptome_fasta, skipping")
        return False

    # Check if index already exists
    index_dir = dest_dir / "index"
    if index_dir.exists() and list(index_dir.glob("*.idx")):
        print(f"  ✓ Index already exists")
        return True

    # Create directories
    dest_dir.mkdir(parents=True, exist_ok=True)
    index_dir.mkdir(parents=True, exist_ok=True)

    # Download transcriptome FASTA
    fasta_path = dest_dir / transcriptome
    download_url = f"{ftp_url.rstrip('/')}/{transcriptome}"

    if not fasta_path.exists():
        print(f"  ▸ Downloading {transcriptome}...")
        result = subprocess.run(
            ["wget", "-q", "-O", str(fasta_path), download_url],
            timeout=300
        )
        if result.returncode != 0 or not fasta_path.exists():
            print(f"  ✗ Download failed")
            return False
        print(f"  ✓ Downloaded ({fasta_path.stat().st_size / 1e6:.1f} MB)")
    else:
        print(f"  ✓ FASTA already exists ({fasta_path.stat().st_size / 1e6:.1f} MB)")

    # Also download other files (genomic, annotation) if present
    for key in ["genomic_fasta", "annotation_gff"]:
        fname = files.get(key, "")
        if fname:
            fpath = dest_dir / fname
            if not fpath.exists():
                url = f"{ftp_url.rstrip('/')}/{fname}"
                subprocess.run(["wget", "-q", "-O", str(fpath), url], timeout=300)

    # Build kallisto index
    species_name = cfg.get("species_list", [species])[0]
    idx_name = f"{species_name}_transcripts.idx"
    idx_path = index_dir / idx_name

    if not idx_path.exists():
        print(f"  ▸ Building kallisto index...")
        result = subprocess.run(
            ["kallisto", "index", "-i", str(idx_path), str(fasta_path)],
            capture_output=True, text=True, timeout=600
        )
        if result.returncode != 0:
            print(f"  ✗ Index build failed: {result.stderr[:200]}")
            return False
        print(f"  ✓ Index built: {idx_path.name}")
    else:
        print(f"  ✓ Index already exists")

    return True


def main():
    parser = argparse.ArgumentParser(description="Download genomes and build kallisto indices")
    parser.add_argument("--config-dir", default="config/amalgkit")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    config_dir = Path(args.config_dir)
    configs = sorted(config_dir.glob("amalgkit_*.yaml"))

    # Skip template/test configs
    configs = [c for c in configs if "template" not in c.name and "test" not in c.name]
    print(f"Found {len(configs)} species configs")

    results = {}
    for i, cfg_path in enumerate(configs, 1):
        species = cfg_path.stem.replace("amalgkit_", "")
        print(f"\n[{i}/{len(configs)}] {species}")
        results[species] = process_species(cfg_path, args.threads)

    # Summary
    ok = sum(1 for v in results.values() if v)
    print(f"\n{'='*50}")
    print(f"Genomes ready: {ok}/{len(results)}")
    for sp, success in sorted(results.items()):
        print(f"  {'✓' if success else '✗'} {sp}")


if __name__ == "__main__":
    main()
