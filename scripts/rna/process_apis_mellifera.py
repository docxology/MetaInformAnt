#!/usr/bin/env python3
"""Simple sequential processor for Apis mellifera samples.

This script processes samples one at a time:
1. Download SRA file using prefetch
2. Extract FASTQ using fasterq-dump
3. Quantify using kallisto
4. Clean up FASTQ files
5. Move to next sample

This avoids the streaming mode complexity in the main workflow.
"""

import csv
import shutil
import subprocess
import sys
import time
from pathlib import Path

# Configuration
WORK_DIR = Path("output/amalgkit/apis_mellifera_all/work")
FASTQ_DIR = Path("output/amalgkit/apis_mellifera_all/fastq/getfastq")
QUANT_DIR = WORK_DIR / "quant"
INDEX_FILE = WORK_DIR / "index/Apis_mellifera_transcripts.idx"
METADATA_FILE = WORK_DIR / "metadata/metadata_selected.tsv"
MAX_BP = 50_000_000  # 50M bp downsampling

# Ensure directories exist
FASTQ_DIR.mkdir(parents=True, exist_ok=True)
QUANT_DIR.mkdir(parents=True, exist_ok=True)


def get_processed_samples():
    """Get set of already quantified samples."""
    processed = set()
    for sample_dir in QUANT_DIR.iterdir():
        if sample_dir.is_dir():
            abundance_file = sample_dir / "abundance.tsv"
            if abundance_file.exists():
                processed.add(sample_dir.name)
    return processed


def download_sample(sample_id: str) -> bool:
    """Download SRA file using prefetch."""
    sample_dir = FASTQ_DIR / sample_id
    sample_dir.mkdir(parents=True, exist_ok=True)

    sra_file = sample_dir / f"{sample_id}.sra"
    if sra_file.exists():
        print(f"  SRA file already exists: {sra_file}")
        return True

    cmd = ["prefetch", "--force", "no", "--max-size", "100G", "-O", str(FASTQ_DIR), sample_id]
    print(f"  Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"  ERROR: prefetch failed: {result.stderr[:200]}")
        return False

    return sra_file.exists()


def extract_fastq(sample_id: str) -> bool:
    """Extract FASTQ from SRA file."""
    sample_dir = FASTQ_DIR / sample_id
    sra_file = sample_dir / f"{sample_id}.sra"

    # Check if already extracted
    fastq_files = list(sample_dir.glob("*.fastq.gz"))
    if fastq_files:
        print(f"  FASTQ files already exist: {len(fastq_files)} files")
        return True

    if not sra_file.exists():
        print(f"  ERROR: SRA file not found: {sra_file}")
        return False

    # Use fasterq-dump with output to sample dir
    cmd = [
        "fasterq-dump",
        "--outdir",
        str(sample_dir),
        "--temp",
        str(sample_dir),
        "--threads",
        "4",
        "--split-3",
        str(sra_file),
    ]
    print(f"  Running: fasterq-dump...")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"  ERROR: fasterq-dump failed: {result.stderr[:200]}")
        return False

    # Compress with gzip
    for fq in sample_dir.glob("*.fastq"):
        subprocess.run(["gzip", "-f", str(fq)], capture_output=True)

    # Delete the .sra file to save space
    if sra_file.exists():
        sra_file.unlink()

    return True


def quantify_sample(sample_id: str) -> bool:
    """Run kallisto quantification."""
    sample_dir = FASTQ_DIR / sample_id
    quant_output = QUANT_DIR / sample_id

    # Check if already quantified
    abundance_file = quant_output / "abundance.tsv"
    if abundance_file.exists():
        print(f"  Already quantified: {abundance_file}")
        return True

    quant_output.mkdir(parents=True, exist_ok=True)

    # Find FASTQ files
    fastq_files = sorted(sample_dir.glob("*.fastq.gz"))
    if not fastq_files:
        print(f"  ERROR: No FASTQ files found in {sample_dir}")
        return False

    # Build kallisto command
    cmd = [
        "kallisto",
        "quant",
        "-i",
        str(INDEX_FILE),
        "-o",
        str(quant_output),
        "-t",
        "4",
    ]

    # Check if paired or single end
    if len(fastq_files) >= 2:
        cmd.extend([str(fastq_files[0]), str(fastq_files[1])])
    else:
        cmd.extend(["--single", "-l", "200", "-s", "30", str(fastq_files[0])])

    print(f"  Running: kallisto quant...")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"  ERROR: kallisto failed: {result.stderr[:200]}")
        return False

    return abundance_file.exists()


def cleanup_sample(sample_id: str):
    """Delete FASTQ files after quantification."""
    sample_dir = FASTQ_DIR / sample_id
    if sample_dir.exists():
        shutil.rmtree(sample_dir, ignore_errors=True)
        print(f"  Cleaned up {sample_dir}")


def process_sample(sample_id: str) -> bool:
    """Process a single sample end-to-end."""
    print(f"\n{'='*60}")
    print(f"Processing: {sample_id}")
    print(f"{'='*60}")

    start_time = time.time()

    # Step 1: Download
    print("Step 1: Download SRA")
    if not download_sample(sample_id):
        return False

    # Step 2: Extract
    print("Step 2: Extract FASTQ")
    if not extract_fastq(sample_id):
        return False

    # Step 3: Quantify
    print("Step 3: Quantify with kallisto")
    if not quantify_sample(sample_id):
        return False

    # Step 4: Cleanup
    print("Step 4: Cleanup")
    cleanup_sample(sample_id)

    elapsed = time.time() - start_time
    print(f"✓ Completed in {elapsed:.1f} seconds")
    return True


def main():
    """Main processing loop."""
    # Verify prerequisites
    if not INDEX_FILE.exists():
        print(f"ERROR: Kallisto index not found: {INDEX_FILE}")
        sys.exit(1)

    if not METADATA_FILE.exists():
        print(f"ERROR: Metadata file not found: {METADATA_FILE}")
        sys.exit(1)

    # Load samples
    with open(METADATA_FILE, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        samples = [row.get("run", "") for row in reader if row.get("run")]

    print(f"Loaded {len(samples)} samples from metadata")

    # Get already processed
    processed = get_processed_samples()
    print(f"Already processed: {len(processed)} samples")

    # Filter to unprocessed
    to_process = [s for s in samples if s not in processed]
    print(f"Remaining to process: {len(to_process)} samples")

    if not to_process:
        print("All samples already processed!")
        return

    # Process samples
    success_count = 0
    fail_count = 0

    for i, sample_id in enumerate(to_process):
        print(f"\n[{i+1}/{len(to_process)}] Processing {sample_id}")

        try:
            if process_sample(sample_id):
                success_count += 1
            else:
                fail_count += 1
        except KeyboardInterrupt:
            print("\nInterrupted by user")
            break
        except Exception as e:
            print(f"ERROR: {e}")
            fail_count += 1

    print(f"\n{'='*60}")
    print(f"SUMMARY")
    print(f"{'='*60}")
    print(f"Processed: {success_count + fail_count}")
    print(f"Successful: {success_count}")
    print(f"Failed: {fail_count}")
    print(f"Total quantified: {len(get_processed_samples())}/{len(samples)}")


if __name__ == "__main__":
    main()
