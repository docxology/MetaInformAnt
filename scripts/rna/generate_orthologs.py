#!/usr/bin/env python3
"""
Generate Ortholog Table from Amalgkit Configs and OrthoDB

This script automatically scans all amalgkit configuration files to extract
the configured `taxon_id`s, downloads the necessary OrthoDB datasets if they
do not exist locally, and runs the ortholog mapping to produce the
cross-species orthogroups.tsv.
"""

import argparse
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path

import yaml


def parse_args():
    parser = argparse.ArgumentParser(description="Automate generation of OrthoDB orthogroup tables.")
    parser.add_argument(
        "--config-dir", default="config/amalgkit", help="Directory containing amalgkit YAML configurations."
    )
    parser.add_argument("--cache-dir", default=".cache/orthodb", help="Directory to cache OrthoDB gzip downloads.")
    parser.add_argument(
        "--out-dir", default="output/amalgkit/cross_species/orthologs", help="Directory to deposit the outcome."
    )
    parser.add_argument("--output-filename", default="orthogroups.tsv", help="Filename of the completed matrix.")
    parser.add_argument("--skip-tests", action="store_true", help="Skip test configs when scanning for species.")
    return parser.parse_args()


def get_taxons_from_configs(config_dir: str, skip_tests: bool = True) -> list:
    taxons = set()
    conf_path = Path(config_dir)
    print(f"Scanning configuration directory: {conf_path}")

    for yaml_file in conf_path.glob("amalgkit_*.yaml"):
        # Exclude cross_species and templates by convention
        if "cross_species" in yaml_file.name or "template" in yaml_file.name:
            continue
        if skip_tests and "test" in yaml_file.name:
            continue

        try:
            with open(yaml_file, "r") as f:
                config = yaml.safe_load(f)

            # Safely get taxon_id at the root level (most common in the template)
            if config and isinstance(config, dict):
                taxon_id = config.get("taxon_id")
                if taxon_id:
                    taxons.add(str(taxon_id))
        except Exception as e:
            print(f"Failed to parse {yaml_file}: {e}")

    print(f"Discovered {len(taxons)} unique taxonomic IDs from configurations.")
    return sorted(list(taxons))


def download_with_resume(url: str, dest: Path):
    """
    Downloads a file securely utilizing urllib natively, supporting HTTP resumes for large files.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    temp_dest = dest.with_suffix(dest.suffix + ".download")

    headers = {}
    if temp_dest.exists():
        existing_size = temp_dest.stat().st_size
        headers["Range"] = f"bytes={existing_size}-"
        print(f"Resuming download from byte {existing_size}...")
    else:
        existing_size = 0

    req = urllib.request.Request(url, headers=headers)
    try:
        with urllib.request.urlopen(req) as response:
            total_size_hdr = response.headers.get("Content-Length")
            total_size = int(total_size_hdr) if total_size_hdr else None

            mode = "ab" if existing_size > 0 else "wb"
            with open(temp_dest, mode) as f:
                downloaded = 0
                chunk_size = 1024 * 1024 * 5  # 5MB chunks
                while True:
                    chunk = response.read(chunk_size)
                    if not chunk:
                        break
                    f.write(chunk)
                    downloaded += len(chunk)

                    if total_size is not None:
                        progress = ((downloaded + existing_size) / (total_size + existing_size)) * 100
                        print(f"\rDownloading {dest.name}... {progress:.1f}%", end="")
                    else:
                        print(f"\rDownloading {dest.name}... {downloaded + existing_size} bytes", end="")

            print()  # new line after progress loop

        # Download strictly succeeded, promote from temporary
        temp_dest.rename(dest)
        print(f"Successfully downloaded {dest.name}")

    except urllib.error.HTTPError as e:
        if e.code == 416:  # Range Not Satisfiable
            print("Download already completely cached based on temp file size.")
            temp_dest.rename(dest)
        else:
            print(f"HTTP Error {e.code}: {url}")
            sys.exit(1)
    except Exception as e:
        print(f"Failed to download {url}: {e}")
        sys.exit(1)


def main():
    args = parse_args()

    taxons = get_taxons_from_configs(args.config_dir, skip_tests=args.skip_tests)
    if not taxons:
        print("No taxonomic IDs found in configurations. Exiting.")
        sys.exit(1)

    cache_dir = Path(args.cache_dir)
    gene_file = cache_dir / "odb12v2_genes.tab.gz"
    mapping_file = cache_dir / "odb12v2_OG2genes.tab.gz"

    # Updated reliable data mirrors (OrthoDB v12.2 data dump endpoint)
    url_base = "https://data.orthodb.org/current/download/odb_data_dump"
    gene_url = f"{url_base}/odb12v2_genes.tab.gz"
    mapping_url = f"{url_base}/odb12v2_OG2genes.tab.gz"

    print("\n--- Phase 1: Obtaining OrthoDB Datasets ---")
    if not gene_file.exists():
        print(f"Missing {gene_file.name}. Commencing robust 4.5GB download from {url_base}...")
        download_with_resume(gene_url, gene_file)
    else:
        print(f"Found cached {gene_file.name}")

    if not mapping_file.exists():
        print(f"Missing {mapping_file.name}. Commencing 500MB download...")
        download_with_resume(mapping_url, mapping_file)
    else:
        print(f"Found cached {mapping_file.name}")

    print("\n--- Phase 2: Generating Custom Ortholog Map ---")
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / args.output_filename

    # Resolve absolute path to sibling map generation script
    script_dir = Path(__file__).resolve().parent
    extractor_script = script_dir / "create_ortholog_table.py"

    if not extractor_script.exists():
        print(f"Error: Required extractor script not found at {extractor_script}")
        sys.exit(1)

    cmd = [
        sys.executable,
        str(extractor_script),
        "--genes",
        str(gene_file),
        "--og2genes",
        str(mapping_file),
        "--output",
        str(out_file),
        "--species-ids",
    ] + taxons

    print(f"Executing: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    print("\nOrtholog Generation Complete!")
    print(f"Output table securely written to: {out_file.absolute()}")


if __name__ == "__main__":
    main()
