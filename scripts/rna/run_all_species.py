#!/usr/bin/env python3
"""ENA-first sample-by-sample RNA-seq pipeline orchestrator.

This script is a thin wrapper around the StreamingPipelineOrchestrator.
It manages current metadata, acquisition, and quantification for multiple
species with robust error handling, retries, and resume capability. Downstream
merge, wsfilter, finalize, and sanity are owned by the external-volume
checkpoint runner after this producer stops.

With ``--cohort DIR`` the species and their tasks come from a frozen,
hash-verified cohort (see ``metainformant.rna.engine.cohort``): Amalgkit
metadata/select never runs and ``metadata_selected.tsv`` is never rewritten.
Without a cohort, a missing or stale metadata table is an error unless
``--refresh-metadata`` explicitly allows the NCBI query.

The producer holds the campaign writer lock of the data root for its whole
life (``metainformant.rna.engine.campaign_lock``), logs INFO to
``<data-root>/logs/run_all_species_<utc>_<pid>.log`` and exits with:

* ``0``  every scheduled task is quantified, quarantined, held by the retry
  bookkeeping (layout mismatch, kallisto exit 1, spent attempts) or
  permanently excluded;
* ``75`` resumable work remains, including tasks in retry backoff (the summary
  names ``next_attempt_due_at``) and species whose discovery failed this run
  (rerun the same command);
* ``73`` another writer holds the data root;
* ``69`` environment or preflight failure (including the environment circuit
  breaker and a refused reuse dry-run);
* ``70`` internal error; ``2`` invalid arguments or cohort.

Resources: without ``--workers``/``--threads`` (or ``PIPELINE_WORKERS``/
``PIPELINE_THREADS``) the profile is derived from this host and the measured
c19 campaign (``build_host_resource_profile``); ``--print-profile`` prints the
resolved profile as JSON and exits without taking the lock or opening the
progress database.

Usage:
    python3 scripts/rna/run_all_species.py --cohort DIR [--cohort-definition-sha256 HEX]
        [--max-gb 1000] [--workers N] [--threads N] [--discovery-workers N]
        [--quant-slots N] [--fastq-threads N] [--fastq-slots N]
        [--compression-threads N] [--validation-slots N] [--max-in-flight N]
        [--acquisition-slots N] [--ena-file-workers 1|2]
        [--scheduling-order size_bands|smallest_first|species_order] [--fast-fail-ena]
        [--requantification-policy preserve|version-drift|all] [--print-profile]
    python3 scripts/rna/run_all_species.py --species solenopsis_invicta --dry-run
"""

import argparse
import json
import os
import shutil
import sys
from dataclasses import replace
from pathlib import Path

# Add src to python path to allow importing metainformant modules
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

from metainformant.rna.engine.cohort import (  # noqa: E402
    CohortError,
    add_producer_cohort_arguments,
    default_max_sample_gb,
    producer_cohort_from_args,
    select_cohort_species,
    verify_cohort_configs,
)
from metainformant.rna.engine.campaign_lock import (  # noqa: E402
    EXIT_ENVIRONMENT,
    EXIT_LOCK_BUSY,
    CampaignLockBusy,
)
from metainformant.rna.engine.campaign_lock import acquire as acquire_campaign_lock  # noqa: E402
from metainformant.rna.engine.species import configured_data_root, discover_species_config_names  # noqa: E402
from metainformant.rna.engine.streaming_orchestrator import (  # noqa: E402
    SCHEDULING_ORDERS,
    PipelineResourceProfile,
    StreamingPipelineOrchestrator,
    apply_fast_fail_ena_environment,
    build_host_resource_profile,
    build_pipeline_resource_profile,
    configure_producer_logging,
    logger,
    producer_exit_code_for_exception,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
CONFIG_DIR = REPO_ROOT / "projects" / "hymenoptera_amalgkit" / "config" / "amalgkit"
if not CONFIG_DIR.is_dir():
    CONFIG_DIR = REPO_ROOT / "config" / "amalgkit"


def _optional_env_int(name: str) -> int | None:
    value = os.environ.get(name, "").strip()
    return int(value) if value else None


DEFAULTS = {
    # Same budget as the campaign launcher (AMALGKIT_MAX_SAMPLE_GB, default 1000 GB).
    "max_gb": default_max_sample_gb(),
    # None: derive both from this host (build_host_resource_profile).
    "workers": _optional_env_int("PIPELINE_WORKERS"),
    "threads": _optional_env_int("PIPELINE_THREADS"),
}


def main() -> int:
    parser = argparse.ArgumentParser(description="ENA-first sample-by-sample RNA-seq pipeline")
    parser.add_argument(
        "--config-dir",
        type=Path,
        default=CONFIG_DIR,
        help=f"Directory containing amalgkit_*.yaml species configs (default: {CONFIG_DIR})",
    )
    parser.add_argument(
        "--species",
        action="append",
        help="Run only this species identifier; repeat for multiple species",
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        help="Data root for work files and the progress database (also sets AMALGKIT_DATA_ROOT)",
    )
    parser.add_argument("--dry-run", action="store_true", help="List the resolved species configs without running")
    parser.add_argument(
        "--max-gb",
        type=float,
        default=DEFAULTS["max_gb"],
        help=f"Max sample size in GB (default: {DEFAULTS['max_gb']})",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=DEFAULTS["workers"],
        help="Parallel sample workers (default: PIPELINE_WORKERS, else acquisition + quant slots of the host profile)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=DEFAULTS["threads"],
        help="Total kallisto threads (default: PIPELINE_THREADS, else one per quant slot of the host profile)",
    )
    parser.add_argument(
        "--acquisition-slots",
        type=int,
        default=None,
        help="Concurrent whole-sample acquisitions, 0 unbounded (host default 8, the measured c19 setting)",
    )
    parser.add_argument(
        "--ena-file-workers",
        type=int,
        choices=(1, 2),
        default=None,
        help="Concurrent ENA files per paired run (host default 2)",
    )
    parser.add_argument(
        "--scheduling-order",
        choices=SCHEDULING_ORDERS,
        default=None,
        help="Fresh-acquisition order across species (default: size_bands, global size bands with species "
        "round-robin; AMALGKIT_PIPELINE_SCHEDULING_ORDER)",
    )
    parser.add_argument(
        "--fast-fail-ena",
        action="store_true",
        help="Opt in to fast-fail ENA: no ENA retries, NCBI prefetch first, 64 KiB/s floor over 20 s "
        "(only for variables not already set)",
    )
    parser.add_argument(
        "--print-profile",
        action="store_true",
        help="Print the resolved resource profile as JSON and exit (no lock, no progress database)",
    )
    parser.add_argument(
        "--discovery-workers",
        type=int,
        default=None,
        help="Concurrent species metadata/reference discovery workers (default: 4; set 1 for serial diagnostics)",
    )
    parser.add_argument(
        "--quant-slots",
        type=int,
        default=None,
        help="Maximum concurrent Kallisto quantifications (default: bounded by workers and total threads)",
    )
    parser.add_argument(
        "--fastq-threads",
        type=int,
        default=None,
        help="Threads for fasterq-dump fallback (default: bounded from the quant budget)",
    )
    parser.add_argument(
        "--fastq-slots",
        type=int,
        default=None,
        help="Maximum concurrent fasterq-dump fallbacks (default: up to four; limits external-disk contention)",
    )
    parser.add_argument(
        "--compression-threads",
        type=int,
        default=None,
        help="Threads per pigz compression process (default: bounded from the quant budget)",
    )
    parser.add_argument(
        "--validation-slots",
        type=int,
        default=None,
        help="Concurrent full local FASTQ validations (default: up to four; limits volume contention)",
    )
    parser.add_argument(
        "--max-in-flight",
        type=int,
        default=None,
        help="Maximum submitted sample tasks (default: two times active workers)",
    )
    parser.add_argument(
        "--requantification-policy",
        choices=("preserve", "version-drift", "all"),
        default=os.environ.get("AMALGKIT_REQUANTIFICATION_POLICY", "preserve"),
        help="Handle contract-compatible version drift (default: preserve)",
    )
    parser.add_argument(
        "--reuse-dry-run",
        type=float,
        default=None,
        metavar="FRACTION",
        help="Before discovery, classify this fraction of each species' recorded-quantified outputs and refuse "
        "(exit 69) if more than 1%% would be quarantined (default: AMALGKIT_PREFLIGHT_REUSE_DRY_RUN, else off)",
    )
    add_producer_cohort_arguments(parser)
    args = parser.parse_args()

    config_dir = args.config_dir.expanduser().resolve()
    if args.data_root:
        os.environ["AMALGKIT_DATA_ROOT"] = str(args.data_root.expanduser().resolve())
    data_root = configured_data_root()
    if args.fast_fail_ena:
        applied = apply_fast_fail_ena_environment()
        print(f"Fast-fail ENA settings applied: {applied or 'none (all already set)'}")
    if args.print_profile:
        try:
            profile = _resolve_profile(args, data_root)
        except ValueError as error:
            parser.error(str(error))
        print(json.dumps(profile.as_dict(), indent=2, sort_keys=True))
        return 0

    lock = None
    if not args.dry_run:
        # Claim the data root before any expensive verification so a second producer fails fast.
        try:
            lock = acquire_campaign_lock(data_root, "producer")
        except CampaignLockBusy as busy:
            print(f"Another writer holds the data root: {busy}", file=sys.stderr)
            return EXIT_LOCK_BUSY
        except (OSError, ValueError) as error:
            print(f"Cannot take the campaign writer lock of {data_root}: {error}", file=sys.stderr)
            return EXIT_ENVIRONMENT
    try:
        return _run(parser, args, config_dir, data_root)
    finally:
        if lock is not None:
            lock.release()


def _cohort_index_bytes(args: argparse.Namespace) -> int | None:
    """The largest frozen kallisto index of the selected cohort species (the quant memory guard input)."""

    try:
        cohort = producer_cohort_from_args(args)
    except CohortError:
        return None
    if cohort is None:
        return None
    sizes = []
    for name in select_cohort_species(cohort, args.species):
        try:
            sizes.append((Path(cohort.root) / cohort.entry(name).index_relpath).stat().st_size)
        except OSError:
            continue
    return max(sizes) if sizes else None


def _host_overrides(args: argparse.Namespace) -> bool:
    """Acquisition budgets exist only in the host profile."""

    return args.acquisition_slots is not None or args.ena_file_workers is not None


def _resolve_profile(args: argparse.Namespace, data_root: Path) -> PipelineResourceProfile:
    """The profile ``run_all`` uses for these arguments: explicit workers and threads keep the historical
    profile unless acquisition budgets are given; everything else is derived from the host."""

    if args.workers is not None and args.threads is not None and not _host_overrides(args):
        profile = build_pipeline_resource_profile(
            args.workers,
            args.threads,
            quant_slots=args.quant_slots,
            fasterq_threads=args.fastq_threads,
            fasterq_slots=args.fastq_slots,
            compression_threads=args.compression_threads,
            validation_slots=args.validation_slots,
            max_in_flight=args.max_in_flight,
        )
        if args.scheduling_order is not None:
            profile = replace(profile, scheduling_order=args.scheduling_order)
        return profile
    try:
        free_bytes: int | None = shutil.disk_usage(data_root).free
    except OSError:
        free_bytes = None
    return build_host_resource_profile(
        kallisto_index_bytes=_cohort_index_bytes(args),
        data_root_free_bytes=free_bytes,
        workers=args.workers,
        threads=args.threads,
        quant_slots=args.quant_slots,
        fasterq_threads=args.fastq_threads,
        fasterq_slots=args.fastq_slots,
        compression_threads=args.compression_threads,
        validation_slots=args.validation_slots,
        max_in_flight=args.max_in_flight,
        acquisition_slots=args.acquisition_slots,
        ena_file_workers=args.ena_file_workers,
        scheduling_order=args.scheduling_order,
    )


def _run(parser: argparse.ArgumentParser, args: argparse.Namespace, config_dir: Path, data_root: Path) -> int:
    try:
        cohort = producer_cohort_from_args(args)
        if cohort is not None:
            selected = select_cohort_species(cohort, args.species)
            verify_cohort_configs(cohort, config_dir, selected)
            config_names = [cohort.entry(name).config_name for name in selected]
        else:
            config_names = discover_species_config_names(config_dir)
            if args.species:
                wanted = set(args.species)
                config_names = [
                    name for name in config_names if name.removeprefix("amalgkit_").removesuffix(".yaml") in wanted
                ]
    except CohortError as error:
        parser.error(str(error))
    if not config_names:
        parser.error(f"no runnable species configurations found under {config_dir}")

    print("╔══════════════════════════════════════════════════════════╗")
    print("║  ENA-First Sample-by-Sample Pipeline (Current)          ║")
    workers_text = args.workers if args.workers is not None else "host"
    threads_text = args.threads if args.threads is not None else "host"
    print(
        f"║  Species: {len(config_names)} | Max: {args.max_gb}GB | Workers: {workers_text} | Threads: {threads_text}  ║"
    )
    print("╚══════════════════════════════════════════════════════════╝")

    print(f"Config directory: {config_dir}")
    print(f"Data root: {data_root}")
    if cohort is not None:
        print(f"Frozen cohort: {cohort.cohort_id} (definition sha256 {cohort.definition_sha256})")
    else:
        print(f"Metadata refresh (NCBI): {'allowed' if args.refresh_metadata else 'refused'}")
    print("Species configs:")
    for config_name in config_names:
        print(f"  - {config_name}")

    if args.dry_run:
        return 0

    log_path = configure_producer_logging(data_root / "logs", prefix="run_all_species")
    if log_path is not None:
        print(f"Run log: {log_path}")
    try:
        orchestrator = StreamingPipelineOrchestrator(
            config_dir=config_dir,
            log_dir=data_root / "logs",
            db_path=data_root / "pipeline_progress.db",
            cohort=cohort,
            allow_metadata_refresh=args.refresh_metadata,
        )
        summary = orchestrator.run_all(
            config_names,
            args.max_gb,
            args.workers,
            args.threads,
            quant_slots=args.quant_slots,
            fasterq_threads=args.fastq_threads,
            fasterq_slots=args.fastq_slots,
            compression_threads=args.compression_threads,
            validation_slots=args.validation_slots,
            max_in_flight=args.max_in_flight,
            requantification_policy=args.requantification_policy,
            discovery_workers=args.discovery_workers,
            reuse_dry_run_fraction=args.reuse_dry_run,
            # Acquisition overrides only exist in the host profile; resolve it here so they apply.
            resource_profile=_resolve_profile(args, data_root) if _host_overrides(args) else None,
            scheduling_order=args.scheduling_order,
        )
    except (Exception, KeyboardInterrupt) as error:  # noqa: BLE001 - every failure maps to an exit code
        code = producer_exit_code_for_exception(error)
        logger.exception("Producer stopped (exit %d): %s", code, error)
        print(f"Producer stopped (exit {code}): {type(error).__name__}: {error}", file=sys.stderr)
        return code
    print(
        f"Producer finished: {summary['status']} (exit {summary['exit_code']}); "
        f"remaining={summary['remaining_by_state']} quantified_this_run={summary['quantified_this_run']} "
        f"held={summary.get('held_by_reason', {})} next_attempt_due_at={summary.get('next_attempt_due_at')} "
        f"discovery_failures={sorted(summary.get('discovery_failures', {}))}"
    )
    return int(summary["exit_code"])


if __name__ == "__main__":
    raise SystemExit(main())
