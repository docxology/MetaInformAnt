"""End-to-end transcriptome SNP-calling pipeline orchestration.

High-level :func:`run_pipeline` composes the variant-calling and stats
primitives into a full per-species workflow: HISAT2 index build, per-sample
FASTQ acquisition/alignment/variant calling/filtering, and population merge
with summaries. CLI entry points live in ``scripts/eqtl/``.

Example:
    >>> from metainformant.eqtl.pipeline import resolve_run_parameters
    >>> params = resolve_run_parameters(species="amellifera")
    >>> params["species"]
    'amellifera'
"""

from __future__ import annotations

import json
import logging
import os
import sys
import time
from pathlib import Path
from typing import Any

from metainformant.eqtl.variant_stats import (
    compute_allele_frequencies,
    compute_popgen_summary,
    compute_sample_stats,
)
from metainformant.eqtl.workflow.variant_calling import (
    align_reads,
    build_hisat2_index,
    call_variants,
    check_tools,
    download_fastq,
    filter_variants,
    find_completed_samples,
    find_reference_genome,
    merge_vcfs,
)

logger = logging.getLogger(__name__)

DEFAULT_AMALGKIT_OUTPUT = (
    Path(os.environ.get("AMALGKIT_DATA_ROOT", "")).expanduser() if os.environ.get("AMALGKIT_DATA_ROOT") else None
)


def resolve_run_parameters(
    species: str | None = None,
    sample_ids: list[str] | None = None,
    n_samples: int = 3,
    threads: int = 4,
    cleanup_fastq: bool = True,
    config: dict[str, Any] | None = None,
    amalgkit_output: Path | None = None,
    eqtl_output: Path | None = None,
) -> dict[str, Any]:
    """Resolve and validate pipeline parameters from CLI args or a config dict.

    Config keys (when ``config`` is given): ``species``, ``samples.mode``
    ("explicit" with ``samples.explicit_ids``), ``samples.max_samples``,
    ``alignment.threads``, ``output.cleanup_fastq``.

    Args:
        species: Species name (e.g. ``amellifera``).
        sample_ids: Explicit SRR IDs (overrides n_samples).
        n_samples: Max samples when sample_ids is None.
        threads: Alignment threads.
        cleanup_fastq: Delete FASTQs after alignment.
        config: Optional YAML-derived config dict (overrides defaults).
        amalgkit_output: Amalgkit data root override.
        eqtl_output: eQTL output root override.

    Returns:
        Dict with resolved keys: species, sample_ids, n_samples, threads,
        cleanup_fastq, amalgkit_output, eqtl_output.

    Raises:
        ValueError: If species is missing or no sample resolution is possible.
    """
    if config is not None:
        species = config.get("species", species)
        samples_cfg = config.get("samples", {})
        if samples_cfg.get("mode") == "explicit" and "explicit_ids" in samples_cfg:
            sample_ids = samples_cfg["explicit_ids"]
        else:
            sample_ids = sample_ids
        n_samples = samples_cfg.get("max_samples", n_samples)
        threads = config.get("alignment", {}).get("threads", threads)
        cleanup_fastq = config.get("output", {}).get("cleanup_fastq", cleanup_fastq)

    if not species:
        raise ValueError("species is required (argument or config 'species' key)")

    resolved_amalgkit = amalgkit_output or DEFAULT_AMALGKIT_OUTPUT
    return {
        "species": species,
        "sample_ids": sample_ids,
        "n_samples": n_samples,
        "threads": threads,
        "cleanup_fastq": cleanup_fastq,
        "amalgkit_output": resolved_amalgkit,
        "eqtl_output": eqtl_output,
    }


def run_pipeline(
    species: str,
    sample_ids: list[str] | None = None,
    n_samples: int = 3,
    threads: int = 4,
    cleanup_fastq: bool = True,
    amalgkit_output: Path | None = None,
    eqtl_output: Path | None = None,
) -> dict[str, Any]:
    """Run the full transcriptome SNP calling pipeline.

    Args:
        species: Species name (must match amalgkit output dir).
        sample_ids: Optional explicit list of SRR IDs.
        n_samples: Number of samples to process (if sample_ids not given).
        threads: Threads for alignment and indexing.
        cleanup_fastq: Delete FASTQs after alignment to save disk space.
        amalgkit_output: Amalgkit data root (defaults to AMALGKIT_DATA_ROOT).
        eqtl_output: eQTL output root (defaults to output/eqtl under cwd).

    Returns:
        Pipeline run summary dict (species, samples, SNP totals, timing).

    Raises:
        SystemExit: When required tools are missing or no samples/genome
            can be resolved.
    """
    start_time = time.time()
    logger.info("=" * 70)
    logger.info(f"Transcriptome SNP Calling Pipeline - {species}")
    logger.info("=" * 70)

    # Check tools
    if not check_tools():
        sys.exit(1)

    if amalgkit_output is None:
        amalgkit_output = DEFAULT_AMALGKIT_OUTPUT
    if eqtl_output is None:
        eqtl_output = Path("output") / "eqtl"

    # Find or use provided samples
    if sample_ids:
        samples = sample_ids
    else:
        if amalgkit_output is None:
            logger.error("No AMALGKIT_DATA_ROOT set and no explicit samples given")
            sys.exit(1)
        samples = find_completed_samples(amalgkit_output, species, max_samples=n_samples)

    if not samples:
        logger.error(f"No completed samples found for {species}")
        sys.exit(1)

    logger.info(f"Processing {len(samples)} samples: {samples}")

    # Find reference genome
    if amalgkit_output is None:
        logger.error("No AMALGKIT_DATA_ROOT set; cannot locate reference genome")
        sys.exit(1)
    ref_genome = find_reference_genome(amalgkit_output, species)
    if not ref_genome:
        logger.error(f"No reference genome found for {species}")
        sys.exit(1)
    logger.info(f"Reference genome: {ref_genome}")

    # Output directories
    species_out = eqtl_output / species
    index_dir = species_out / "index"
    pop_dir = species_out / "population"
    log_dir = species_out / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    # Set up file logging
    fh = logging.FileHandler(log_dir / "pipeline.log")
    fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logger.addHandler(fh)

    # Step 1: Build HISAT2 index
    logger.info("\n[Step 1/6] Building HISAT2 index...")
    index_prefix = build_hisat2_index(ref_genome, index_dir, threads=threads)

    # Process each sample
    filtered_vcfs: list[Path] = []
    all_sample_stats: list[dict[str, Any]] = []

    for i, srr_id in enumerate(samples, 1):
        logger.info(f"\n{'=' * 50}")
        logger.info(f"Sample {i}/{len(samples)}: {srr_id}")
        logger.info(f"{'=' * 50}")

        sample_dir = species_out / "samples" / srr_id
        fastq_dir = sample_dir / "fastq"
        bam_path = sample_dir / "aligned.bam"
        raw_vcf = sample_dir / "variants_raw.vcf.gz"
        filt_vcf = sample_dir / "variants.vcf.gz"
        stats_json = sample_dir / "variant_stats.json"

        # Skip if already fully processed
        if filt_vcf.exists() and filt_vcf.stat().st_size > 0 and stats_json.exists():
            logger.info(f"Sample {srr_id} already processed, loading cached stats")
            with open(stats_json) as f:
                stats = json.load(f)
            stats["sample_id"] = srr_id
            all_sample_stats.append(stats)
            filtered_vcfs.append(filt_vcf)
            logger.info(f"  {srr_id}: {stats.get('n_snps', '?')} SNPs (cached)")
            continue

        # Step 2: Download FASTQs (checks local amalgkit copies first)
        logger.info(f"\n[Step 2/6] Downloading FASTQs for {srr_id}...")
        fastqs = download_fastq(srr_id, fastq_dir, species=species, amalgkit_output=amalgkit_output)
        if not fastqs:
            logger.warning(f"Skipping {srr_id}: no FASTQs downloaded")
            continue

        # Step 3: Align reads
        logger.info(f"\n[Step 3/6] Aligning {srr_id}...")
        ok = align_reads(fastqs, index_prefix, bam_path, threads=threads)
        if not ok:
            logger.warning(f"Skipping {srr_id}: alignment failed")
            continue

        # Step 4: Call variants
        logger.info(f"\n[Step 4/6] Calling variants for {srr_id}...")
        ok = call_variants(bam_path, ref_genome, raw_vcf)
        if not ok:
            logger.warning(f"Skipping {srr_id}: variant calling failed")
            continue

        # Step 5: Filter variants
        logger.info(f"\n[Step 5/6] Filtering variants for {srr_id}...")
        ok = filter_variants(raw_vcf, filt_vcf)
        if not ok:
            logger.warning(f"Skipping {srr_id}: filtering failed")
            continue

        # Compute per-sample stats
        stats = compute_sample_stats(filt_vcf, stats_json)
        stats["sample_id"] = srr_id
        all_sample_stats.append(stats)
        filtered_vcfs.append(filt_vcf)

        logger.info(f"  {srr_id}: {stats['n_snps']} SNPs, {stats['n_indels']} indels")

        # Cleanup FASTQs to save disk
        if cleanup_fastq and fastq_dir.exists():
            import shutil

            shutil.rmtree(fastq_dir)
            logger.info(f"  Cleaned up FASTQs for {srr_id}")

    # Step 6: Merge and population analysis
    logger.info(f"\n[Step 6/6] Merging {len(filtered_vcfs)} sample VCFs...")
    if filtered_vcfs:
        merged_vcf = pop_dir / "merged.vcf.gz"
        merge_vcfs(filtered_vcfs, merged_vcf)

        # Allele frequencies
        compute_allele_frequencies(merged_vcf, pop_dir / "allele_freqs.tsv")

        # Population summary
        pop_summary = compute_popgen_summary(merged_vcf, all_sample_stats, pop_dir / "popgen_summary.json")
    else:
        pop_summary = {"n_samples": 0, "error": "No samples completed successfully"}

    elapsed = time.time() - start_time

    run_summary = {
        "species": species,
        "n_samples_requested": len(samples),
        "n_samples_completed": len(filtered_vcfs),
        "samples": [s["sample_id"] for s in all_sample_stats],
        "total_snps": pop_summary.get("total_snps", 0),
        "total_indels": pop_summary.get("total_indels", 0),
        "ts_tv_ratio": pop_summary.get("ts_tv_ratio", 0.0),
        "elapsed_seconds": round(elapsed, 1),
        "output_dir": str(species_out),
    }

    with open(species_out / "run_summary.json", "w") as f:
        json.dump(run_summary, f, indent=2)

    logger.info("\n" + "=" * 70)
    logger.info("Pipeline Complete!")
    logger.info(f"  Species: {species}")
    logger.info(f"  Samples: {run_summary['n_samples_completed']}/{run_summary['n_samples_requested']}")
    logger.info(f"  SNPs: {run_summary['total_snps']}")
    logger.info(f"  Ti/Tv: {run_summary['ts_tv_ratio']}")
    logger.info(f"  Time: {elapsed:.0f}s")
    logger.info(f"  Output: {species_out}")
    logger.info("=" * 70)

    return run_summary
