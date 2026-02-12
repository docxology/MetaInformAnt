#!/usr/bin/env python3
"""ENA-first sample-by-sample RNA-seq pipeline orchestrator.

Strategy:
  1. Process species one at a time (smallest sample-count first)
  2. Within each species, sort samples smallest → largest
  3. For EACH sample: download FASTQ from ENA → quant immediately
  4. After all samples: run merge/curate/sanity
  5. Configurable workers, threads, and max sample size

Usage:
    python3 scripts/rna/run_all_species.py [--max-gb 5.0] [--workers 4] [--threads 12]
"""

import argparse
import concurrent.futures
import os
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

# ── Configuration ──────────────────────────────────────────────────

CONFIG_DIR = "config/amalgkit"
LOG_DIR = "blue/amalgkit"

# Species processing order: all ants (smallest first), then bees
SPECIES_ORDER = [
    # ── Ants (21 species, sorted by sample count ascending) ──
    "amalgkit_anoplolepis_gracilipes.yaml",      #    2 samples
    "amalgkit_acromyrmex_echinatior.yaml",        #    8 samples (filtered ≤5GB)
    "amalgkit_dinoponera_quadriceps.yaml",        #   13 samples
    "amalgkit_vollenhovia_emeryi.yaml",           #   15 samples
    "amalgkit_odontomachus_brunneus.yaml",        #   19 samples
    "amalgkit_formica_exsecta.yaml",              #   23 samples
    "amalgkit_temnothorax_americanus.yaml",       #   32 samples
    "amalgkit_wasmannia_auropunctata.yaml",       #   33 samples
    "amalgkit_nylanderia_fulva.yaml",             #   40 samples
    "amalgkit_temnothorax_curvispinosus.yaml",    #   43 samples
    "amalgkit_pbarbatus.yaml",                    #   95 samples
    "amalgkit_cardiocondyla_obscurior.yaml",      #  162 samples
    "amalgkit_temnothorax_nylanderi.yaml",        #  166 samples
    "amalgkit_linepithema_humile.yaml",           #  173 samples
    "amalgkit_atta_cephalotes.yaml",              #  220 samples
    "amalgkit_ooceraea_biroi.yaml",               #  237 samples
    "amalgkit_camponotus_floridanus.yaml",        #  304 samples
    "amalgkit_solenopsis_invicta.yaml",           #  349 samples
    "amalgkit_monomorium_pharaonis.yaml",         #  370 samples
    "amalgkit_temnothorax_longispinosus.yaml",    #  508 samples
    "amalgkit_harpegnathos_saltator.yaml",        #  689 samples
    # ── Bees ──
    "amalgkit_amellifera.yaml",                   # 3154 samples (filtered ≤5GB)
]

DEFAULTS = {
    "max_gb": 5.0,
    "workers": 16,
    "threads": 12,
}


# ── ENA Download Logic ────────────────────────────────────────────

def query_ena_fastq_urls(srr_id: str) -> list[str]:
    """Query ENA API for direct FASTQ download URLs."""
    url = (
        f"https://www.ebi.ac.uk/ena/portal/api/filereport"
        f"?accession={srr_id}&result=read_run"
        f"&fields=run_accession,fastq_ftp&format=tsv"
    )
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            lines = resp.read().decode().strip().split('\n')
            if len(lines) < 2:
                return []
            ftp_field = lines[1].split('\t')[1] if '\t' in lines[1] else ''
            if not ftp_field:
                return []
            return [f"https://{p}" for p in ftp_field.split(';') if p]
    except Exception:
        return []


def download_fastq_from_ena(srr_id: str, out_dir: Path) -> bool:
    """Download FASTQ files directly from ENA using curl."""
    urls = query_ena_fastq_urls(srr_id)
    if not urls:
        return False

    sample_dir = out_dir / srr_id
    sample_dir.mkdir(parents=True, exist_ok=True)

    # Check if already downloaded
    existing_fq = list(sample_dir.glob("*.fastq.gz"))
    if len(existing_fq) >= len(urls):
        return True  # Already done

    success = True
    for url in urls:
        fname = url.split('/')[-1]
        fpath = sample_dir / fname
        if fpath.exists() and fpath.stat().st_size > 0:
            continue
        try:
            result = subprocess.run(
                ["curl", "-L", "-f", "-o", str(fpath), "--retry", "3",
                 "--retry-delay", "5", "-s", "--show-error", url],
                timeout=7200,
                capture_output=True, text=True
            )
            if result.returncode != 0:
                fpath.unlink(missing_ok=True)
                success = False
            else:
                sz = fpath.stat().st_size / (1024**3)
                print(f"    ✓ {fname}: {sz:.2f} GB", flush=True)
        except (subprocess.TimeoutExpired, Exception):
            fpath.unlink(missing_ok=True)
            success = False

    return success


# ── Quantification Logic ──────────────────────────────────────────

def quant_sample(config_path: str, batch_index: int, species_name: str, threads: int) -> bool:
    """Run amalgkit quant for a single sample using --batch."""
    # Read config to get quant out_dir
    import yaml
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    steps_cfg = cfg.get("steps", {})
    quant_cfg = steps_cfg.get("quant", {})
    quant_out = quant_cfg.get("out_dir", f"blue/amalgkit/{species_name}/work")
    getfastq_cfg = steps_cfg.get("getfastq", {})
    fastq_out = getfastq_cfg.get("out_dir", f"blue/amalgkit/{species_name}/fastq")

    # Find metadata
    meta_path = None
    for mp in [
        Path(f"blue/amalgkit/{species_name}/work/metadata/metadata.tsv"),
        Path(f"blue/amalgkit/{species_name}/work/metadata/metadata_selected.tsv"),
    ]:
        if mp.exists():
            meta_path = str(mp)
            break

    if not meta_path:
        return False

    cmd = [
        "amalgkit", "quant",
        "--out_dir", quant_out,
        "--metadata", meta_path,
        "--threads", str(threads),
        "--batch", str(batch_index),
    ]
    cmd.extend(["--clean_fastq", "no"])

    # Check for genome index
    genome_cfg = cfg.get("genome", {})
    index_dir = genome_cfg.get("index_dir", "")
    if index_dir and Path(index_dir).exists():
        cmd.extend(["--index_dir", index_dir])
    fasta_dir = genome_cfg.get("fasta_dir", "")
    if fasta_dir and Path(fasta_dir).exists():
        cmd.extend(["--fasta_dir", fasta_dir])

    log_path = Path(LOG_DIR) / f"{species_name}_quant.log"
    
    # Use capture_output to avoid concurrent file handle issues (Bad file descriptor)
    # Write to log before and after
    try:
        with open(log_path, 'a') as log_f:
            log_f.write(f"\n[{time.strftime('%Y-%m-%d %H:%M:%S')}] Quant batch {batch_index} START: {' '.join(cmd)}\n")
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
        
        with open(log_path, 'a') as log_f:
            log_f.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Quant batch {batch_index} END (Exit {result.returncode})\n")
            if result.stdout:
                log_f.write(f"--- STDOUT ---\n{result.stdout}\n")
            if result.stderr:
                log_f.write(f"--- STDERR ---\n{result.stderr}\n")
    except Exception as e:
        with open(log_path, 'a') as log_f:
            log_f.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Quant batch {batch_index} EXCEPTION: {e}\n")
        return False

    return result.returncode == 0


# ── Metadata Helpers ──────────────────────────────────────────────

def load_metadata(tsv_path: Path) -> tuple[list[dict], list[str]]:
    """Load metadata TSV into list of dicts + headers."""
    with open(tsv_path) as f:
        headers = f.readline().strip().split('\t')
        rows = []
        for line in f:
            vals = line.strip().split('\t')
            row = dict(zip(headers, vals))
            rows.append(row)
    return rows, headers


def get_srr_column(rows: list[dict]) -> str:
    """Find the SRR/run accession column name."""
    for col in ['run_accession', 'run', 'Run', 'accession']:
        if rows and col in rows[0]:
            return col
    return 'run'


def filter_and_sort_metadata(rows: list[dict], max_gb: float) -> list[dict]:
    """Filter by max size, sort ascending by total_bases, enforce required fields."""
    max_bases = max_gb * 1e9
    filtered = []
    for r in rows:
        try:
            tb = float(r.get('total_bases', 0))
        except (ValueError, TypeError):
            continue
        if 0 < tb <= max_bases:
            # Ensure amalgkit-required fields are set (prevents strtobool crash)
            r['is_sampled'] = 'yes'
            r['is_qualified'] = 'yes'
            if not r.get('exclusion') or r['exclusion'] in ('', 'nan', 'NaN'):
                r['exclusion'] = 'no'
            filtered.append(r)
    filtered.sort(key=lambda r: float(r.get('total_bases', 0)))
    return filtered


def write_metadata(rows: list[dict], out_path: Path, headers: list[str]):
    """Write rows back to TSV."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        f.write('\t'.join(headers) + '\n')
        for r in rows:
            f.write('\t'.join(r.get(h, '') for h in headers) + '\n')


# ── Per-Sample Worker ─────────────────────────────────────────────

def is_sample_quantified(species_name: str, srr_id: str) -> bool:
    """Check if a sample already has quantification output ({SRR}_abundance.tsv)."""
    quant_dir = Path(f"blue/amalgkit/{species_name}/work/quant/{srr_id}")
    return any(quant_dir.glob("*_abundance.tsv")) if quant_dir.exists() else False


def process_single_sample(
    srr_id: str, batch_index: int, fastq_dir: Path,
    config_path: str, species_name: str, threads: int
) -> dict:
    """Download + quant a single sample. Returns status dict."""
    result = {
        "srr": srr_id,
        "batch": batch_index,
        "downloaded": False,
        "quantified": False,
        "skipped": False,
        "error": None,
    }

    # Skip if already quantified
    if is_sample_quantified(species_name, srr_id):
        result["downloaded"] = True
        result["quantified"] = True
        result["skipped"] = True
        return result

    try:
        # Step 1: Download FASTQ from ENA
        ok = download_fastq_from_ena(srr_id, fastq_dir)
        if not ok:
            result["error"] = "ENA download failed"
            return result
        result["downloaded"] = True

        # Step 2: Immediately quantify this sample
        ok = quant_sample(config_path, batch_index, species_name, threads)
        result["quantified"] = ok
        if not ok:
            result["error"] = "quant failed (may need genome index)"

    except Exception as e:
        result["error"] = str(e)

    return result


# ── Genome Index Pre-Check ────────────────────────────────────────

def verify_genome_index(config_path: str, species_name: str) -> bool:
    """Verify Kallisto index exists BEFORE downloading any samples.
    
    Returns True if index found, False if missing.
    """
    import yaml
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    # Check config-specified index_dir
    steps_cfg = cfg.get('steps', {})
    quant_cfg = steps_cfg.get('quant', {})
    index_dir = quant_cfg.get('index_dir', '')
    genome_cfg = cfg.get('genome', {})
    if not index_dir:
        index_dir = genome_cfg.get('index_dir', '')

    # Also check default locations
    search_dirs = [
        index_dir,
        f'blue/amalgkit/{species_name}/work/index',
        f'blue/amalgkit/shared/genome/{species_name.replace("_", "_").title().replace("_", "_")}/index',
    ]

    # Derive scientific name from config for index filename matching
    species_list = cfg.get('species_list', [])
    sci_name = species_list[0] if species_list else species_name.replace('_', ' ').title().replace(' ', '_')

    for d in search_dirs:
        if not d or not os.path.isdir(d):
            continue
        idx_files = [f for f in os.listdir(d) if f.endswith('.idx')]
        if idx_files:
            print(f"  ✓ Genome index found: {os.path.join(d, idx_files[0])}", flush=True)
            return True

    # Also check via amalgkit's quant out_dir/index path
    quant_out = quant_cfg.get('out_dir', f'blue/amalgkit/{species_name}/work')
    fallback_index_dir = os.path.join(quant_out, 'index')
    if os.path.isdir(fallback_index_dir):
        idx_files = [f for f in os.listdir(fallback_index_dir) if f.endswith('.idx')]
        if idx_files:
            print(f"  ✓ Genome index found: {os.path.join(fallback_index_dir, idx_files[0])}", flush=True)
            return True

    print(f"  ✗ NO GENOME INDEX found for {species_name}!", flush=True)
    print(f"    Searched: {search_dirs}", flush=True)
    return False


# ── Per-Species Workflow ──────────────────────────────────────────

def process_species(config_name: str, max_gb: float, workers: int, threads: int):
    """Process a single species: sample-by-sample download+quant."""
    config_path = os.path.join(CONFIG_DIR, config_name)
    if not os.path.exists(config_path):
        print(f"  ✗ Config not found: {config_path}", flush=True)
        return False

    species_name = config_name.replace("amalgkit_", "").replace(".yaml", "")
    print(f"\n{'='*60}", flush=True)
    print(f"  Processing: {species_name}", flush=True)
    print(f"  Config: {config_path}", flush=True)
    print(f"  Max size: {max_gb} GB | Workers: {workers} | Threads: {threads}", flush=True)
    print(f"{'='*60}", flush=True)

    # Find metadata
    metadata_path = None
    for mp in [
        Path(f"blue/amalgkit/{species_name}/work/metadata/metadata_selected.tsv"),
        Path(f"blue/amalgkit/{species_name}/work/metadata/metadata.tsv"),
    ]:
        if mp.exists():
            metadata_path = mp
            break

    if not metadata_path:
        print(f"  ✗ No metadata found for {species_name}", flush=True)
        return False

    # Load, filter, sort metadata
    rows, headers = load_metadata(metadata_path)
    srr_col = get_srr_column(rows)
    filtered = filter_and_sort_metadata(rows, max_gb)
    print(f"  Samples: {len(rows)} total → {len(filtered)} after filter (≤{max_gb}GB)", flush=True)

    if not filtered:
        print(f"  ⏩ No samples within size limit for {species_name}", flush=True)
        return True

    # ── Verify genome index BEFORE downloading anything ──
    if not verify_genome_index(config_path, species_name):
        print(f"  ✗ Skipping {species_name}: No genome index. Build index first!", flush=True)
        return True

    print(f"  Size range: {float(filtered[0]['total_bases'])/1e9:.2f} – {float(filtered[-1]['total_bases'])/1e9:.2f} GB", flush=True)

    # Write sorted+filtered metadata (amalgkit quant --batch reads from this)
    sorted_meta = Path(f"blue/amalgkit/{species_name}/work/metadata/metadata.tsv")
    write_metadata(filtered, sorted_meta, headers)
    # Also sync metadata_selected.tsv
    selected_meta = Path(f"blue/amalgkit/{species_name}/work/metadata/metadata_selected.tsv")
    write_metadata(filtered, selected_meta, headers)

    fastq_dir = Path(f"blue/amalgkit/{species_name}/fastq/getfastq")
    fastq_dir.mkdir(parents=True, exist_ok=True)

    # ── Process samples: download + quant per sample ──────────
    print(f"\n  ── Sample-by-Sample Processing ({len(filtered)} samples) ──", flush=True)

    downloaded = 0
    quantified = 0
    skipped = 0
    failed = 0

    # Use thread pool for concurrent download+quant
    threads_per_worker = max(1, threads // workers)

    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {}
        for i, row in enumerate(filtered):
            srr = row[srr_col]
            batch_idx = i + 1  # amalgkit uses 1-based batch index
            future = executor.submit(
                process_single_sample,
                srr, batch_idx, fastq_dir,
                config_path, species_name, threads_per_worker
            )
            futures[future] = (srr, batch_idx)

        for future in concurrent.futures.as_completed(futures):
            srr, batch_idx = futures[future]
            try:
                result = future.result()
                if result.get("skipped"):
                    skipped += 1
                    quantified += 1
                    downloaded += 1
                elif result["quantified"]:
                    quantified += 1
                    downloaded += 1
                    print(f"  ✓ {srr} (#{batch_idx}): downloaded + quantified", flush=True)
                elif result["downloaded"]:
                    downloaded += 1
                    print(f"  ⚠ {srr} (#{batch_idx}): downloaded, quant failed: {result['error']}", flush=True)
                else:
                    failed += 1
                    print(f"  ✗ {srr} (#{batch_idx}): {result['error']}", flush=True)
            except Exception as e:
                failed += 1
                print(f"  ✗ {srr} (#{batch_idx}): exception: {e}", flush=True)

    new_quants = quantified - skipped
    print(f"\n  Summary: {quantified} quantified ({skipped} already done, {new_quants} new), {failed} failed", flush=True)

    # ── Always run downstream steps if we have any quantified samples ──
    # Check if downstream completion flags already exist
    curate_flag = Path(f"blue/amalgkit/{species_name}/work/curate")
    sanity_dir = Path(f"blue/amalgkit/{species_name}/work/sanity")
    has_curate_flag = any(curate_flag.rglob("curate_completion_flag.txt")) if curate_flag.exists() else False
    has_sanity = (sanity_dir.exists() and any(sanity_dir.iterdir())) if sanity_dir.exists() else False

    if quantified > 0 and not (has_curate_flag and has_sanity):
        print(f"\n  ── Running downstream steps (merge → sanity) ──", flush=True)
        cmd = [
            "python3", "scripts/rna/run_workflow.py",
            "--config", config_path,
            "--no-progress",
            "--steps", "merge", "curate", "sanity"
        ]
        log_path = os.path.join(LOG_DIR, os.path.basename(config_path) + ".log")
        with open(log_path, 'a') as log_f:
            log_f.write(f"\n[{time.strftime('%Y-%m-%d %H:%M:%S')}] Downstream: {' '.join(cmd)}\n")
            log_f.flush()
            result = subprocess.run(cmd, stdout=log_f, stderr=subprocess.STDOUT)
        if result.returncode == 0:
            print(f"  ✓ Downstream steps complete!", flush=True)
        else:
            print(f"  ⚠ Downstream steps had errors (check logs)", flush=True)
    elif has_curate_flag and has_sanity:
        print(f"  ✓ Already complete (curate + sanity done)", flush=True)

    return quantified > 0 or downloaded > 0


# ── Main ──────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="ENA-first sample-by-sample RNA-seq pipeline"
    )
    parser.add_argument("--max-gb", type=float, default=DEFAULTS["max_gb"],
                        help=f"Max sample size in GB (default: {DEFAULTS['max_gb']})")
    parser.add_argument("--workers", type=int, default=DEFAULTS["workers"],
                        help=f"Parallel workers (default: {DEFAULTS['workers']})")
    parser.add_argument("--threads", type=int, default=DEFAULTS["threads"],
                        help=f"Total threads (split across workers) (default: {DEFAULTS['threads']})")
    args = parser.parse_args()

    print(f"╔══════════════════════════════════════════════════════════╗", flush=True)
    print(f"║  ENA-First Sample-by-Sample Pipeline                    ║", flush=True)
    print(f"║  Species: {len(SPECIES_ORDER)} | Max: {args.max_gb}GB | Workers: {args.workers} | Threads: {args.threads}  ║", flush=True)
    print(f"╚══════════════════════════════════════════════════════════╝", flush=True)

    os.makedirs(LOG_DIR, exist_ok=True)

    results = {}
    for config_name in SPECIES_ORDER:
        species = config_name.replace("amalgkit_", "").replace(".yaml", "")
        print(f"\n🚀 Starting {species}...", flush=True)
        t0 = time.time()
        try:
            ok = process_species(config_name, args.max_gb, args.workers, args.threads)
            elapsed = time.time() - t0
            results[species] = f"{'✓' if ok else '✗'} ({elapsed/60:.0f} min)"
        except Exception as e:
            elapsed = time.time() - t0
            print(f"  ✗ {species}: Fatal error: {e}", flush=True)
            results[species] = f"✗ Error ({elapsed/60:.0f} min): {e}"

    print(f"\n{'='*60}", flush=True)
    print(f"  Final Summary", flush=True)
    print(f"{'='*60}", flush=True)
    for species, result in results.items():
        print(f"  {species}: {result}", flush=True)


if __name__ == "__main__":
    main()
