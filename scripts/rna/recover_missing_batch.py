
import os
import sys
import subprocess
import shutil
import time
from pathlib import Path

# Ensure src is in path
sys.path.append(str(Path(__file__).parent.parent.parent))
from metainformant.rna.retrieval.ena_downloader import download_sra_samples

def get_missing_samples(list_file):
    with open(list_file) as f:
        return [l.strip() for l in f if l.strip()]

def run_quant_batch(samples, work_dir, fastq_dir, genome_dir, index_dir):
    """Run amalgkit quant on a specific batch of samples."""
    print(f"  Quantifying batch: {samples}")
    
    # We can pass IDs via --id (comma separated? No, likely one by one or --id_list)
    # Amalgkit quant with --id_list is best.
    
    batch_list_file = work_dir / "batch_ids.txt"
    with open(batch_list_file, "w") as f:
        for s in samples:
            f.write(f"{s}\n")
            
    cmd = [
        "amalgkit", "quant",
        "--out_dir", str(work_dir),
        "--id_list", str(batch_list_file),
        "--threads", "8",
        "--index_dir", str(index_dir),
        "--fasta_dir", str(genome_dir),
        "--redo", "yes", # Force redo since we just downloaded them
        "--fastq_dir", str(fastq_dir) # Point to where ENA downloaded them
    ]
    
    env = os.environ.copy()
    # Ensure current venv python is used
    
    print(f"  Command: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True, env=env)
        return True
    except subprocess.CalledProcessError as e:
        print(f"  [FAIL] Quant failed for batch: {e}")
        return False

def clean_batch(fastq_dir):
    """Delete all files in the fastq/getfastq directory."""
    getfastq_dir = fastq_dir / "getfastq"
    print(f"  Cleaning up {getfastq_dir}...")
    if getfastq_dir.exists():
        for item in getfastq_dir.iterdir():
            if item.is_dir():
                shutil.rmtree(item)
            else:
                item.unlink()
        print("  Cleanup complete.")

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", required=True, help="List of missing samples")
    parser.add_argument("--batch-size", type=int, default=5)
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--fastq-dir", required=True)
    parser.add_argument("--genome-dir", required=True)
    parser.add_argument("--index-dir", required=True)
    
    args = parser.parse_args()
    
    all_samples = get_missing_samples(args.file)
    print(f"Total samples to process: {len(all_samples)}")
    
    work_dir = Path(args.work_dir)
    fastq_dir = Path(args.fastq_dir)
    
    # Process in batches
    for i in range(0, len(all_samples), args.batch_size):
        batch = all_samples[i:i+args.batch_size]
        print(f"\n=== Processing Batch {i//args.batch_size + 1} ({len(batch)} samples) ===")
        print(f"Samples: {batch}")
        
        # 1. Download
        print("Step 1: Downloading from ENA...")
        download_sra_samples(batch, fastq_dir) # Uses our fixed module
        
        # 2. Quantify
        print("Step 2: Quantifying...")
        # Note: Amalgkit quant expects the fastq dir structure. 
        # ena_downloader creates fastq_dir/getfastq/SRR/SRR.fastq.gz
        # amalgkit quant looks in fastq_dir/getfastq/...
        run_quant_batch(batch, work_dir, fastq_dir, args.genome_dir, args.index_dir)
        
        # 3. Cleanup
        print("Step 3: Cleaning up FASTQs...")
        clean_batch(fastq_dir)
        
        print("Batch complete.\n")

if __name__ == "__main__":
    main()
