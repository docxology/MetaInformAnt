import argparse
import concurrent.futures
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import List, Tuple

# Ensure src is in path
sys.path.append(str(Path(__file__).parent.parent.parent / "src"))
try:
    from metainformant.rna.retrieval.ena_downloader import (
        download_sra_samples, 
        get_ena_sample_info,
        parse_size,
        format_size
    )
except ImportError:
    # Fallback if running from root
    sys.path.append(str(Path.cwd() / "src"))
    from metainformant.rna.retrieval.ena_downloader import (
        download_sra_samples,
        get_ena_sample_info,
        parse_size,
        format_size
    )


def get_missing_samples(list_file):
    with open(list_file) as f:
        return [l.strip() for l in f if l.strip()]


def create_single_metadata(sra_id, full_metadata_path, work_dir):
    """Create a temporary metadata file containing only the target SRA."""
    out_path = work_dir / f"metadata_{sra_id}.tsv"
    
    # Check if metadata file exists
    if not full_metadata_path.exists():
        print(f"[{sra_id}] Metadata file not found: {full_metadata_path}")
        return None

    with open(full_metadata_path) as f_in, open(out_path, "w") as f_out:
        # Read header
        header = f_in.readline()
        f_out.write(header)

        header_parts = header.strip().split("\t")
        try:
            run_col_idx = header_parts.index("run")
        except ValueError:
            raise ValueError("Metadata file does not contain a 'run' column in the header")

        found = False
        for line in f_in:
            parts = line.strip().split("\t")
            if len(parts) > run_col_idx and parts[run_col_idx] == sra_id:
                f_out.write(line)
                found = True
                break

    if not found:
        out_path.unlink(missing_ok=True)
        return None

    return out_path


def run_quant_single(sra_id, work_dir, fastq_dir, genome_dir, index_dir, full_metadata_file):
    """Run amalgkit quant on a SINGLE sample via per-sample metadata."""
    print(f"[{sra_id}] Preparing metadata...")
    single_metadata = create_single_metadata(sra_id, full_metadata_file, work_dir)

    if not single_metadata:
        print(f"[{sra_id}] [FAIL] Could not find ID in metadata file.")
        return False

    print(f"[{sra_id}] Quantifying...")

    # Try to find amalgkit in PATH
    amalgkit_exe = shutil.which("amalgkit")
    if not amalgkit_exe:
        # Fallback to .venv/bin/amalgkit
        try_path = Path.cwd() / ".venv" / "bin" / "amalgkit"
        if try_path.exists():
            amalgkit_exe = str(try_path)
        else:
            amalgkit_exe = "amalgkit"

    cmd = [
        amalgkit_exe,
        "quant",
        "--out_dir",
        str(work_dir),
        "--metadata",
        str(single_metadata),
        "--threads",
        "2",
        "--index_dir",
        str(index_dir),
        "--fasta_dir",
        str(genome_dir),
        "--redo",
        "yes",
    ]

    env = os.environ.copy()

    try:
        log_file = work_dir / "logs" / f"quant_{sra_id}.log"
        log_file.parent.mkdir(exist_ok=True)
        with open(log_file, "w") as log:
            subprocess.run(cmd, check=True, env=env, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)
        print(f"[{sra_id}] Quant finished.")
        # Clean up temp metadata
        single_metadata.unlink(missing_ok=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"[{sra_id}] [FAIL] Quant failed. See logs: {log_file}")
        single_metadata.unlink(missing_ok=True)
        return False


def cleanup_single(sra_id, work_dir):
    """Delete files for this specific sample."""
    sample_dir = work_dir / "getfastq" / sra_id
    if sample_dir.exists():
        try:
            shutil.rmtree(sample_dir)
        except Exception as e:
            print(f"[{sra_id}] Cleanup failed: {e}")


def process_sample(sra_id, args):
    """Full lifecycle for one sample."""
    work_dir = Path(args.work_dir)
    genome_dir = Path(args.genome_dir)
    index_dir = Path(args.index_dir)
    metadata_file = work_dir / "metadata" / "metadata.tsv"
    
    # Parse constraints
    max_size_bytes = parse_size(args.max_size) if args.max_size else None
    
    print(f"\n[{sra_id}] Starting pipeline...")

    # 1. Download with new best practices (size limit, timeout, fallback)
    try:
        download_sra_samples(
            [sra_id], 
            work_dir,
            max_size_bytes=max_size_bytes,
            timeout=args.timeout,
            sort_by_size=False,  # Sorting handled at batch level
            use_fallback=args.fallback
        )
    except Exception as e:
        print(f"[{sra_id}] Download crashed: {e}")
        return False

    # Check validity after download attempt
    sample_path = work_dir / "getfastq" / sra_id
    if not sample_path.exists() or not list(sample_path.glob("*.fastq.gz")):
        print(f"[{sra_id}] Download seemingly failed (no files). Skipping quant.")
        return False

    # 2. Quantify
    success = run_quant_single(sra_id, work_dir, sample_path, genome_dir, index_dir, metadata_file)

    # 3. Cleanup
    cleanup_single(sra_id, work_dir)

    return success


def get_size_info(sra_id):
    """Helper for parallel size fetching."""
    info = get_ena_sample_info(sra_id)
    return (sra_id, info.size_bytes if info else float('inf'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", required=True, help="List of missing samples")
    parser.add_argument("--max-workers", type=int, default=3)
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--fastq-dir", help="Deprecated.")
    parser.add_argument("--genome-dir", required=True)
    parser.add_argument("--index-dir", required=True)
    
    # New options
    parser.add_argument("--max-size", help="Skip samples larger than this")
    parser.add_argument("--timeout", type=int, default=600, help="Download timeout")
    parser.add_argument("--sort-by-size", action="store_true", help="Process smallest first")
    parser.add_argument("--fallback", action="store_true", help="Use fallback sources")

    args = parser.parse_args()

    all_samples = get_missing_samples(args.file)
    print(f"Total samples to process: {len(all_samples)}")
    print(f"Parallel Workers: {args.max_workers}")

    # Parallel Sorting Phase
    if args.sort_by_size:
        print("Sorting samples by size (checking ENA in parallel)...")
        sample_sizes = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=20) as executor:
            future_to_id = {executor.submit(get_size_info, sra_id): sra_id for sra_id in all_samples}
            for i, future in enumerate(concurrent.futures.as_completed(future_to_id)):
                if (i+1) % 50 == 0:
                    print(f"Checked size for {i+1}/{len(all_samples)} samples...")
                try:
                    sample_sizes.append(future.result())
                except Exception:
                    sample_sizes.append((future_to_id[future], float('inf')))
        
        # Sort smallest first
        sample_sizes.sort(key=lambda x: x[1])
        all_samples = [s[0] for s in sample_sizes]
        print(f"Sorted {len(all_samples)} samples by size.")

    # Pre-create directories to avoid race conditions
    work_dir = Path(args.work_dir)
    (work_dir / "getfastq").mkdir(parents=True, exist_ok=True)
    (work_dir / "logs").mkdir(parents=True, exist_ok=True)
    
    # Execution Phase
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        futures = {executor.submit(process_sample, sra_id, args): sra_id for sra_id in all_samples}

        for future in concurrent.futures.as_completed(futures):
            sra_id = futures[future]
            try:
                result = future.result()
                status = "SUCCESS" if result else "FAILURE"
                print(f"[{sra_id}] Pipeline Finished: {status}")
            except Exception as e:
                print(f"[{sra_id}] Pipeline Exception: {e}")


if __name__ == "__main__":
    main()
