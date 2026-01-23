
import os
import sys
import subprocess
import shutil
import time
import concurrent.futures
from pathlib import Path

# Ensure src is in path
sys.path.append(str(Path(__file__).parent.parent.parent / "src"))
try:
    from metainformant.rna.retrieval.ena_downloader import download_sra_samples
except ImportError:
    # Fallback if running from root
    sys.path.append(str(Path.cwd() / "src"))
    from metainformant.rna.retrieval.ena_downloader import download_sra_samples

def get_missing_samples(list_file):
    with open(list_file) as f:
        return [l.strip() for l in f if l.strip()]

def create_single_metadata(sra_id, full_metadata_path, work_dir):
    """Create a temporary metadata file containing only the target SRA."""
    out_path = work_dir / f"metadata_{sra_id}.tsv"
    
    with open(full_metadata_path) as f_in, open(out_path, "w") as f_out:
        # Read header
        header = f_in.readline()
        f_out.write(header)
        
        header_parts = header.strip().split('\t')
        try:
            run_col_idx = header_parts.index("run")
        except ValueError:
             raise ValueError("Metadata file does not contain a 'run' column in the header")
             
        found = False
        for line in f_in:
            parts = line.strip().split('\t')
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
            # Fallback for when running via python directly? wrapper script?
            # We know it is installed now.
            amalgkit_exe = "amalgkit"
    
    cmd = [
        amalgkit_exe, "quant",
        "--out_dir", str(work_dir),
        "--metadata", str(single_metadata),
        "--threads", "4", 
        "--index_dir", str(index_dir),
        "--fasta_dir", str(genome_dir),
        "--redo", "yes"
    ]
    
    env = os.environ.copy()
    
    try:
        log_file = work_dir / "logs" / f"quant_{sra_id}.log"
        log_file.parent.mkdir(exist_ok=True)
        with open(log_file, "w") as log:
            subprocess.run(cmd, check=True, env=env, stdout=log, stderr=subprocess.STDOUT)
        print(f"[{sra_id}] Quant finished.")
        # Clean up temp metadata
        single_metadata.unlink(missing_ok=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"[{sra_id}] [FAIL] Quant failed. See logs: {log_file}")
        # Clean up temp metadata even on fail
        single_metadata.unlink(missing_ok=True)
        return False

def cleanup_single(sra_id, work_dir):
    """Delete files for this specific sample."""
    sample_dir = work_dir / "getfastq" / sra_id
    print(f"[{sra_id}] Cleaning up {sample_dir}...")
    if sample_dir.exists():
        try:
            shutil.rmtree(sample_dir)
            print(f"[{sra_id}] Cleanup complete.")
        except Exception as e:
            print(f"[{sra_id}] Cleanup failed: {e}")

def process_sample(sra_id, args):
    """Full lifecycle for one sample."""
    work_dir = Path(args.work_dir)
    genome_dir = Path(args.genome_dir)
    index_dir = Path(args.index_dir)
    metadata_file = work_dir / "metadata" / "metadata.tsv"
    
    print(f"\n[{sra_id}] Starting pipeline...")
    
    # 1. Download
    try:
        # Always call downloader. It handles resumption (-C -) and MD5 checks.
        download_sra_samples([sra_id], work_dir)
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

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", required=True, help="List of missing samples")
    parser.add_argument("--max-workers", type=int, default=3)
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--fastq-dir", help="Deprecated.") 
    parser.add_argument("--genome-dir", required=True)
    parser.add_argument("--index-dir", required=True)
    
    args = parser.parse_args()
    
    all_samples = get_missing_samples(args.file)
    print(f"Total samples to process: {len(all_samples)}")
    print(f"Parallel Workers: {args.max_workers}")
    
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
