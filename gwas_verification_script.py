"""GWAS verification script for Apis mellifera pipeline features.

This script:
1. Generates a synthetic Apis mellifera dataset (16 chromosomes, 200 samples)
2. Runs the full GWAS pipeline including:
   - QC (MAF, HWE)
   - PCA (population structure)
   - Association testing (Linear, Mixed Model)
   - Fine-mapping (SuSiE)
   - Heritability estimation (LDSC)
   - Benchmarking (new feature)
3. Verifies outputs and prints summary.
"""

import sys
import shutil
from pathlib import Path
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger("gwas_verify")

from metainformant.gwas.workflow.workflow_execution import run_gwas, run_multi_trait_gwas
from metainformant.gwas.data.config import estimate_runtime
from metainformant.gwas.analysis.benchmarking import benchmark_subset_run, extrapolate_full_genome_time
from tests.fixtures.gwas.generate_test_data import generate_complete_test_dataset

def main():
    work_dir = Path("gwas_verification")
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir()
    
    data_dir = work_dir / "data"
    results_dir = work_dir / "results"
    
    logger.info("Generating synthetic Apis mellifera dataset...")
    # 200 samples, 16 chromosomes (Amel), 1000 variants/chrom = 16k total
    # Use NCBI names (NC_037638.1 etc)
    dataset = generate_complete_test_dataset(
        output_dir=data_dir,
        n_variants_per_chrom=1000,
        n_samples=500,
        n_chroms=16,
        n_causal=5,
        effect_size=3.0,
        use_ncbi_names=True,
        seed=42
    )
    
    logger.info(f"Generated data at {data_dir}")
    logger.info(f"VCF: {dataset['vcf_path']}")
    logger.info(f"Phenotype: {dataset['phenotype_path']}")

    # --- 1. Benchmarking ---
    logger.info("\n--- Step 1: Compute-Time Benchmarking ---")
    config_bench = {
        "threads": 4,
        "model": "linear",
        "create_plots": False 
    }
    
    # Run benchmark on subset
    timings = benchmark_subset_run(
        dataset["vcf_path"], 
        dataset["phenotype_path"], 
        config_bench, 
        max_samples=50, 
        max_variants=2000
    )
    
    # Extrapolate to full genome (pretend full genome is 10M variants)
    est = extrapolate_full_genome_time(
        timings, 
        target_n_samples=200, 
        target_n_variants=10_000_000
    )
    logger.info("\nBenchmarking Report:")
    logger.info(est.summary())

    # --- 2. Full GWAS Run ---
    logger.info("\n--- Step 2: Full GWAS Workflow Execution ---")
    
    config_full = {
        "study_name": "Apis_Verification",
        "description": "Verification run for Apis mellifera synthetic data",
        "vcf_file": str(dataset["vcf_path"]),
        "phenotype_file": str(dataset["phenotype_path"]),
        "output_dir": str(results_dir),
        
        "maf_threshold": 0.05,
        "hwe_p_threshold": 1e-4,
        
        "association_method": "mixed", # Use mixed model for full test
        "kinship_method": "vanraden",
        "pca_components": 5,
        
        "multiple_testing_correction": "fdr",
        
        "create_plots": True,
        "plot_types": ["manhattan", "qq", "pca"],
        
        "threads": 4,
        "memory_gb": 8
    }
    
    # Runtime estimation using new data-driven function
    runtime_est = estimate_runtime(config_full, n_samples=500, n_variants=16000)
    logger.info(f"Pre-run estimate: {runtime_est['estimated_seconds']:.2f}s")
    
    result = run_gwas(
        dataset["vcf_path"],
        dataset["phenotype_path"],
        config_full,
        output_dir=results_dir
    )
    
    if result["status"] == "success":
        logger.info("GWAS Pipeline completed successfully!")
    else:
        logger.error(f"GWAS Pipeline failed: {result.get('error')}")
        sys.exit(1)
        
    # --- 3. Output Verification ---
    logger.info("\n--- Step 3: Output Verification ---")
    
    # Check key files
    expected_files = [
        "summary_statistics.tsv",
        "significant_hits.tsv",
        "Manhattan_Plot.png",
        "QQ_Plot.png",
        "pca_plot_PC1_PC2.png",
        "results_summary.json"
    ]
    
    for fname in expected_files:
        path = results_dir / fname
        if path.exists():
            size = path.stat().st_size
            logger.info(f"[OK] {fname} ({size} bytes)")
        else:
            logger.error(f"[MISSING] {fname}")
            
    # Check significant hits
    hits = result["results"].get("significant_hits", [])
    logger.info(f"Found {len(hits)} significant associations (FDR < 0.05)")
    if hits:
        top_hit = hits[0]
        logger.info(f"Top hit: {top_hit['chrom']}:{top_hit['pos']} p={top_hit['p_value']:.2e}")

    # Check heritability
    h2 = result["results"].get("heritability", {}).get("h2", "N/A")
    logger.info(f"Heritability estimate: {h2}")

if __name__ == "__main__":
    main()
