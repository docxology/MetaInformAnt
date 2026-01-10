#!/usr/bin/env python3
"""Main workflow orchestrator for amalgkit RNA-seq analysis.

This is a thin wrapper that calls methods from metainformant.rna.orchestration
to run complete workflows for single species.

Usage:
    # Run full workflow for a species
    python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml

    # Run specific steps
    python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps getfastq quant merge

    # Check status
    python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --status

    # Cleanup unquantified samples
    python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --cleanup-unquantified
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Import setup utilities (must be before other imports)
sys.path.insert(0, str(Path(__file__).parent))
from _setup_utils import ensure_venv_activated, check_environment_or_exit

# Suppress optional dependency warnings until venv is ready
from metainformant.core.utils.optional_deps import suppress_optional_warnings, enable_optional_warnings
suppress_optional_warnings()

# Auto-setup and activate venv
ensure_venv_activated(auto_setup=True)
check_environment_or_exit(auto_setup=True)

# Now enable warnings since venv should be active
enable_optional_warnings()

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.rna.core.cleanup import cleanup_partial_downloads, fix_abundance_naming_for_species
from metainformant.rna.engine.monitoring import analyze_species_status, check_workflow_progress
from metainformant.rna.engine.orchestration import (
    cleanup_unquantified_samples,
    run_workflow_for_species,
)
from metainformant.core.utils.logging import get_logger

logger = get_logger("run_workflow")

REPO_ROOT = Path(__file__).parent.parent.parent.resolve()


def _invocation_hint() -> str:
    """Return a short command prefix users can copy/paste based on how this script is invoked."""
    argv0 = Path(sys.argv[0]).name
    if argv0 == "run_workflow.py":
        return "python run_workflow.py"
    return "python3 scripts/rna/run_workflow.py"


def _discover_species_config_files() -> list[Path]:
    config_dir = REPO_ROOT / "config" / "amalgkit"
    if not config_dir.exists():
        return []
    out: list[Path] = []
    for p in sorted(config_dir.glob("amalgkit_*.yaml")):
        stem = p.stem.lower()
        if "template" in stem or "test" in stem:
            continue
        out.append(p)
    return out


def _resolve_config_path(raw: Path) -> Path:
    """Resolve a config path provided on CLI.

    Supports:
    - absolute paths
    - paths relative to current working directory
    - paths relative to repo root
    - bare filenames resolved under repo_root/config/amalgkit/
    """
    p = Path(raw).expanduser()
    if p.is_absolute() and p.exists():
        return p.resolve()

    cwd_candidate = (Path.cwd() / p).expanduser()
    if cwd_candidate.exists():
        return cwd_candidate.resolve()

    repo_candidate = (REPO_ROOT / p).expanduser()
    if repo_candidate.exists():
        return repo_candidate.resolve()

    if p.parent == Path("."):
        config_candidate = (REPO_ROOT / "config" / "amalgkit" / p.name).expanduser()
        if config_candidate.exists():
            return config_candidate.resolve()

    try:
        return p.resolve()
    except Exception:
        return p


def main() -> int:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Run amalgkit workflow for a single species",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("config_pos", nargs="?", type=Path, help="Path to species workflow config file")
    parser.add_argument(
        "--config",
        type=Path,
        help="Path to species workflow config file",
    )
    parser.add_argument(
        "--steps",
        nargs="+",
        help="Specific steps to run (default: all steps)",
    )
    parser.add_argument(
        "--status",
        action="store_true",
        help="Check workflow status instead of running",
    )
    parser.add_argument(
        "--detailed",
        action="store_true",
        help="Show detailed status (use with --status)",
    )
    parser.add_argument(
        "--cleanup-unquantified",
        action="store_true",
        help="Quantify downloaded samples and cleanup FASTQs",
    )
    parser.add_argument(
        "--cleanup-partial",
        action="store_true",
        help="Clean up partial downloads",
    )
    parser.add_argument(
        "--fix-abundance-naming",
        action="store_true",
        help="Fix abundance file naming for merge compatibility",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Stop on first failure",
    )
    parser.add_argument(
        "--list-configs",
        action="store_true",
        help="List available amalgkit species configs and exit",
    )
    parser.add_argument(
        "--plan",
        action="store_true",
        help="Print planned step order and the exact amalgkit commands that would run, then exit",
    )
    parser.add_argument(
        "--walk",
        action="store_true",
        help="Pause before each stage and wait for Enter (TTY only)",
    )
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable progress bars (falls back to plain logs)",
    )
    parser.add_argument(
        "--show-commands",
        action="store_true",
        help="Log the exact amalgkit command for each stage before running it",
    )
    parser.add_argument(
        "--validate",
        action="store_true",
        help="Run validation on workflow samples without executing workflow",
    )
    parser.add_argument(
        "--validate-stage",
        choices=["download", "extraction", "quantification", "merge"],
        help="Validate specific pipeline stage (use with --validate)",
    )

    args = parser.parse_args()

    # Process steps argument - split comma-separated values
    if args.steps:
        args.steps = [step.strip() for s in args.steps for step in s.split(',') if step.strip()]

    if args.list_configs:
        files = _discover_species_config_files()
        if not files:
            logger.error(f"No configs found under: {REPO_ROOT / 'config' / 'amalgkit'}")
            return 2
        logger.info("Available configs:")
        for p in files:
            logger.info(f"  - {p.relative_to(REPO_ROOT)}")
        logger.info("Example:")
        logger.info(f"  {_invocation_hint()} config/amalgkit/<file>.yaml --plan")
        return 0

    raw_config = args.config or args.config_pos
    if raw_config is None:
        # Default no-args behavior: show available configs and how to run them.
        # This keeps the command useful without accidentally starting a long workflow.
        files = _discover_species_config_files()
        if not files:
            logger.error("No config provided, and no configs were found under `config/amalgkit/`.")
            logger.info(f"Tip: create one under `config/amalgkit/` or pass `--config <path>`.")
            return 2
        logger.info("No config provided. Available configs:")
        for p in files:
            logger.info(f"  - {p.relative_to(REPO_ROOT)}")
        logger.info("Next:")
        logger.info(f"  {_invocation_hint()} config/amalgkit/<file>.yaml --plan")
        logger.info(f"  {_invocation_hint()} config/amalgkit/<file>.yaml --check")
        return 0

    config_path = _resolve_config_path(raw_config)
    if not config_path.exists():
        logger.error(f"Config file not found: {config_path}")
        return 1

    if args.plan:
        from metainformant.rna.amalgkit.amalgkit import build_amalgkit_command
        from metainformant.rna.engine.workflow import apply_step_defaults, load_workflow_config, plan_workflow, sanitize_params_for_cli

        cfg = load_workflow_config(config_path)
        apply_step_defaults(cfg)
        planned = plan_workflow(cfg)
        logger.info(f"Planned {len(planned)} steps for {config_path.name}:")
        for idx, (step, params) in enumerate(planned, start=1):
            cmd = build_amalgkit_command(step, sanitize_params_for_cli(step, params))
            logger.info(f"  {idx:02d}. {step}: {' '.join(cmd)}")
        return 0

    # Status check
    if args.status:
        if args.detailed:
            status = analyze_species_status(config_path)
        else:
            status = check_workflow_progress(config_path)
        logger.info(f"Status: {status}")
        return 0

    # Cleanup operations
    if args.cleanup_unquantified:
        logger.info("Cleaning up unquantified samples...")
        quantified, failed = cleanup_unquantified_samples(config_path)
        logger.info(f"Quantified: {quantified}, Failed: {failed}")
        return 0 if failed == 0 else 1

    if args.cleanup_partial:
        logger.info("Cleaning up partial downloads...")
        result = cleanup_partial_downloads(config_path, dry_run=False)
        logger.info(f"Deleted: {result['deleted']}, Freed: {result['freed_mb']}MB")
        return 0 if result["errors"] == 0 else 1

    if args.fix_abundance_naming:
        logger.info("Fixing abundance file naming...")
        created, exists = fix_abundance_naming_for_species(config_path)
        logger.info(f"Created: {created}, Already exists: {exists}")
        return 0

    # Validation
    if args.validate:
        from metainformant.rna.engine.workflow import load_workflow_config
        from metainformant.rna.analysis.validation import validate_all_samples, save_validation_report
        
        logger.info(f"Running validation for {config_path.name}")
        config = load_workflow_config(config_path)
        
        validation_result = validate_all_samples(config, stage=args.validate_stage)
        
        # Save validation report
        validation_dir = config.work_dir / "validation"
        validation_dir.mkdir(parents=True, exist_ok=True)
        report_file = validation_dir / "validation_report.json"
        save_validation_report(validation_result, report_file)
        
        # Print summary
        total = validation_result.get('total_samples', 0)
        validated = validation_result.get('validated', 0)
        failed = validation_result.get('failed', 0)
        stage = args.validate_stage or 'all'
        
        logger.info(f"Validation results ({stage}):")
        logger.info(f"  Total samples: {total}")
        logger.info(f"  Validated: {validated}")
        logger.info(f"  Failed: {failed}")
        
        if validation_result.get('summary'):
            logger.info("  Stage breakdown:")
            for stage_name, stage_info in validation_result['summary'].items():
                complete = stage_info.get('complete', 0)
                missing = stage_info.get('missing', 0)
                logger.info(f"    {stage_name}: {complete} complete, {missing} missing")
        
        if failed > 0:
            logger.warning(f"  {failed} samples failed validation. Check {report_file} for details.")
            return 1
        else:
            logger.info("  All samples passed validation")
            return 0

    # Run robust download if getfastq is in steps (or steps is None/all)
    should_run_download = args.steps is None or 'getfastq' in args.steps
    
    if should_run_download and not args.validate: # Skip if validating
        # Determine metadata path and fastq dir
        # We need to load config to get work_dir
        from metainformant.rna.engine.workflow import load_workflow_config
        cfg = load_workflow_config(config_path)
        
        metadata_path = cfg.work_dir / "metadata" / "metadata_selected.tsv"
        fastq_dir = cfg.work_dir.parent / "fastq" / "getfastq" # Default layout assumption
        
        steps_config = cfg.extra_config.get('steps', {})
        if 'getfastq' in steps_config and 'out_dir' in steps_config['getfastq']:
             fastq_dir = Path(steps_config['getfastq']['out_dir'])
             if fastq_dir.name != "getfastq":
                 fastq_dir = fastq_dir / "getfastq"

        if metadata_path.exists():
            print(f"DEBUG: Metadata found at {metadata_path}, running robust download...")
            logger.info(f"Running robust pre-download for SRA files to {fastq_dir}...")
            from metainformant.core.io.download_robust import download_sra_files_from_metadata
            download_sra_files_from_metadata(metadata_path, fastq_dir)
        else:
            print(f"DEBUG: Metadata NOT found at {metadata_path}")
            logger.info("Metadata selected file not found yet, skipping robust pre-download (will run in normal flow)")
    else:
        print(f"DEBUG: Skipping robust download. should_run_download={should_run_download}")

    # Run workflow
    logger.info(f"Running workflow for {config_path.name}")
    results = run_workflow_for_species(
        config_path,
        steps=args.steps,
        check=args.check,
        walk=args.walk,
        progress=not args.no_progress,
        show_commands=args.show_commands,
    )

    if results["success"]:
        logger.info(f"✅ Workflow completed: {len(results.get('completed', []))} steps")
        return 0
    else:
        logger.error(f"❌ Workflow failed: {len(results.get('failed', []))} steps failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())

