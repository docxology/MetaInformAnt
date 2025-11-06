#!/usr/bin/env python3
"""Quality control simulation script.

This script generates synthetic FASTQ quality scores, contamination patterns,
and quality metrics for testing quality control pipelines.

Usage:
    python3 scripts/simulation/simulate_quality.py --type fastq --n-reads 10000 --read-length 100
    python3 scripts/simulation/simulate_quality.py --type contamination --n-reads 5000 --contamination-rate 0.1
    python3 scripts/simulation/simulate_quality.py --type metrics --n-samples 20
"""

import argparse
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core import io, logging, paths, validation
from metainformant.simulation.sequences import generate_random_dna

logger = logging.get_logger(__name__)


def simulate_fastq(
    output_dir: Path,
    n_reads: int,
    read_length: int,
    mean_quality: int,
    quality_degradation: float,
    seed: int,
) -> dict:
    """Simulate FASTQ quality scores."""
    logger.info(f"Generating FASTQ quality data: {n_reads} reads, length {read_length}")
    rng = random.Random(seed)
    np.random.seed(seed)
    
    # Generate quality scores (Phred+33 encoding)
    # Quality scores degrade along the read
    reads = []
    
    for read_idx in range(n_reads):
        # Generate sequence
        seq = generate_random_dna(read_length, rng=rng)
        
        # Generate quality scores
        qualities = []
        for pos in range(read_length):
            # Quality degrades along read
            base_quality = mean_quality - (pos * quality_degradation)
            base_quality = max(2, min(40, base_quality + np.random.normal(0, 2)))
            qualities.append(int(base_quality))
        
        # Convert to Phred+33 ASCII
        quality_string = "".join([chr(q + 33) for q in qualities])
        
        reads.append({
            "read_id": f"read_{read_idx:06d}",
            "sequence": seq,
            "quality": quality_string,
            "mean_quality": np.mean(qualities),
        })
    
    # Save as simplified FASTQ-like format
    fastq_file = output_dir / "simulated_reads.fastq"
    with open(fastq_file, "w") as f:
        for read in reads:
            f.write(f"@{read['read_id']}\n")
            f.write(f"{read['sequence']}\n")
            f.write("+\n")
            f.write(f"{read['quality']}\n")
    
    # Save quality metrics
    metrics_df = pd.DataFrame(reads)
    metrics_file = output_dir / "quality_metrics.csv"
    metrics_df.to_csv(metrics_file, index=False)
    
    logger.info(f"FASTQ data saved to {fastq_file}")
    
    return {
        "type": "fastq",
        "n_reads": n_reads,
        "read_length": read_length,
        "output_file": str(fastq_file),
        "metrics_file": str(metrics_file),
    }


def simulate_contamination(
    output_dir: Path,
    n_reads: int,
    read_length: int,
    contamination_rate: float,
    seed: int,
) -> dict:
    """Simulate contamination patterns."""
    logger.info(f"Generating contamination data: {n_reads} reads, {contamination_rate} contamination rate")
    rng = random.Random(seed)
    np.random.seed(seed)
    
    n_contaminated = int(n_reads * contamination_rate)
    
    reads = []
    contamination_sources = ["adapter", "vector", "rRNA", "cross_species"]
    
    for read_idx in range(n_reads):
        is_contaminated = read_idx < n_contaminated
        
        if is_contaminated:
            source = rng.choice(contamination_sources)
            # Contaminated reads have different characteristics
            if source == "adapter":
                # Adapter sequence at start/end
                seq = "AGATCGGAAG" + generate_random_dna(read_length - 10, rng=rng)
            elif source == "vector":
                # Vector sequence
                seq = generate_random_dna(read_length, gc_content=0.6, rng=rng)
            elif source == "rRNA":
                # High GC content
                seq = generate_random_dna(read_length, gc_content=0.7, rng=rng)
            else:  # cross_species
                # Different composition
                seq = generate_random_dna(read_length, gc_content=0.3, rng=rng)
            
            contamination_type = source
        else:
            seq = generate_random_dna(read_length, gc_content=0.5, rng=rng)
            contamination_type = "none"
        
        reads.append({
            "read_id": f"read_{read_idx:06d}",
            "sequence": seq,
            "is_contaminated": is_contaminated,
            "contamination_type": contamination_type,
        })
    
    # Save sequences
    fasta_file = output_dir / "contamination_reads.fasta"
    with open(fasta_file, "w") as f:
        for read in reads:
            f.write(f">{read['read_id']}\n")
            f.write(f"{read['sequence']}\n")
    
    # Save metadata
    metadata_df = pd.DataFrame(reads)
    metadata_file = output_dir / "contamination_metadata.csv"
    metadata_df.to_csv(metadata_file, index=False)
    
    logger.info(f"Contamination data saved to {fasta_file}")
    
    return {
        "type": "contamination",
        "n_reads": n_reads,
        "contamination_rate": contamination_rate,
        "n_contaminated": n_contaminated,
        "output_file": str(fasta_file),
        "metadata_file": str(metadata_file),
    }


def simulate_metrics(
    output_dir: Path,
    n_samples: int,
    seed: int,
) -> dict:
    """Simulate quality control metrics."""
    logger.info(f"Generating quality metrics: {n_samples} samples")
    rng = random.Random(seed)
    np.random.seed(seed)
    
    metrics = []
    for sample_idx in range(n_samples):
        metrics.append({
            "sample_id": f"sample_{sample_idx:03d}",
            "total_reads": rng.randint(1000000, 5000000),
            "total_bases": rng.randint(100000000, 500000000),
            "mean_read_length": rng.uniform(80, 150),
            "mean_quality": rng.uniform(25, 35),
            "gc_content": rng.uniform(0.4, 0.6),
            "duplication_rate": rng.uniform(0.1, 0.5),
            "adapter_contamination": rng.uniform(0.0, 0.05),
            "n_contaminated_reads": rng.randint(0, 10000),
        })
    
    df = pd.DataFrame(metrics)
    csv_file = output_dir / "quality_metrics.csv"
    df.to_csv(csv_file, index=False)
    
    logger.info(f"Quality metrics saved to {csv_file}")
    
    return {
        "type": "metrics",
        "n_samples": n_samples,
        "output_file": str(csv_file),
    }


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Quality control simulation for testing and validation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Simulate FASTQ quality scores
  %(prog)s --type fastq --n-reads 10000 --read-length 100

  # Simulate contamination
  %(prog)s --type contamination --n-reads 5000 --contamination-rate 0.1

  # Simulate quality metrics
  %(prog)s --type metrics --n-samples 20
        """,
    )
    
    parser.add_argument(
        "--type",
        required=True,
        choices=["fastq", "contamination", "metrics"],
        help="Type of simulation to run",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/simulation/quality"),
        help="Output directory (default: output/simulation/quality)",
    )
    parser.add_argument("--n-reads", type=int, default=10000, help="Number of reads (fastq/contamination types)")
    parser.add_argument("--read-length", type=int, default=100, help="Read length (fastq/contamination types)")
    parser.add_argument("--mean-quality", type=int, default=30, help="Mean quality score (fastq type)")
    parser.add_argument("--quality-degradation", type=float, default=0.1, help="Quality degradation per position (fastq type)")
    parser.add_argument("--contamination-rate", type=float, default=0.1, help="Contamination rate (contamination type)")
    parser.add_argument("--n-samples", type=int, default=20, help="Number of samples (metrics type)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")
    
    args = parser.parse_args()
    
    if args.verbose:
        import logging as std_logging
        logger.setLevel(std_logging.DEBUG)
    
    output_dir = paths.ensure_directory(args.output)
    
    try:
        if args.type == "fastq":
            results = simulate_fastq(
                output_dir,
                args.n_reads,
                args.read_length,
                args.mean_quality,
                args.quality_degradation,
                args.seed,
            )
        elif args.type == "contamination":
            results = simulate_contamination(
                output_dir, args.n_reads, args.read_length, args.contamination_rate, args.seed
            )
        elif args.type == "metrics":
            results = simulate_metrics(output_dir, args.n_samples, args.seed)
        
        # Save summary
        summary_file = output_dir / "simulation_summary.json"
        io.dump_json(results, summary_file, indent=2)
        logger.info(f"Simulation complete. Summary saved to {summary_file}")
        
        return 0
    except Exception as e:
        logger.error(f"Simulation failed: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())

