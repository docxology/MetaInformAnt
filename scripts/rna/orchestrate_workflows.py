#!/usr/bin/env python3
"""
Unified Amalgkit Workflow Orchestrator

Intelligently manages all amalgkit workflows across all species with flexible step execution.

Usage:
    # Full status assessment
    python3 scripts/rna/orchestrate_workflows.py --assess

    # Cleanup: quantify downloaded samples and delete FASTQs
    python3 scripts/rna/orchestrate_workflows.py --cleanup-unquantified
    
    # Run specific steps for specific species
    python3 scripts/rna/orchestrate_workflows.py --species cfloridanus sinvicta --steps merge curate sanity
    
    # Run downstream analysis for all complete species
    python3 scripts/rna/orchestrate_workflows.py --steps merge curate sanity --auto-species
    
    # Resume incomplete downloads
    python3 scripts/rna/orchestrate_workflows.py --resume-downloads
    
    # Full pipeline
    python3 scripts/rna/orchestrate_workflows.py --species pbarbatus --steps metadata select getfastq quant merge curate sanity
"""

import argparse
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional

import yaml

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.rna import count_quantified_samples, find_unquantified_samples
from metainformant.rna.steps import quantify_sample, delete_sample_fastqs
from metainformant.rna.workflow import load_workflow_config
from metainformant.core.logging import get_logger
from metainformant.core.io import read_delimited

logger = get_logger("orchestrate")

# Repository paths
REPO_ROOT = Path(__file__).parent.parent.parent
CONFIG_DIR = REPO_ROOT / "config" / "amalgkit"
OUTPUT_DIR = REPO_ROOT / "output" / "amalgkit"
LOG_DIR = REPO_ROOT / "output"

# Species metadata
SPECIES_INFO = {
    'cfloridanus': {'name': 'Camponotus floridanus', 'total': 307},
    'pbarbatus': {'name': 'Pogonomyrmex barbatus', 'total': 83},
    'mpharaonis': {'name': 'Monomorium pharaonis', 'total': 100},
    'sinvicta': {'name': 'Solenopsis invicta', 'total': 354},
}

# Valid amalgkit steps in order
VALID_STEPS = [
    'metadata', 'integrate', 'config', 'select', 
    'getfastq', 'quant', 'merge', 
    'cstmm', 'curate', 'csca', 'sanity'
]


def get_config_path(species: str) -> Path | None:
    """Get config path for a species."""
    config_path = CONFIG_DIR / f'amalgkit_{species}.yaml'
    if config_path.exists():
        return config_path
    return None


def count_quantified(species: str) -> int:
    """Count successfully quantified samples using metainformant."""
    config_path = get_config_path(species)
    if not config_path:
        # Fallback to direct counting if no config
        quant_dir = OUTPUT_DIR / species / 'quant'
        if not quant_dir.exists():
            return 0
        return len([d for d in quant_dir.iterdir() 
                    if d.is_dir() and (d / 'abundance.tsv').exists()])
    
    quantified, _ = count_quantified_samples(config_path)
    return quantified


def count_downloaded_unquantified(species: str) -> list[str]:
    """Find samples with FASTQs but no quantification using metainformant."""
    config_path = get_config_path(species)
    if not config_path:
        # Fallback to direct checking if no config
        fastq_dir = OUTPUT_DIR / species / 'fastq'
        quant_dir = OUTPUT_DIR / species / 'quant'
        
        if not fastq_dir.exists():
            return []
        
        unquantified = []
        for sample_dir in fastq_dir.iterdir():
            if not sample_dir.is_dir():
                continue
            
            sample_id = sample_dir.name
            
            # Check if quantified
            if (quant_dir / sample_id / 'abundance.tsv').exists():
                continue
            
            # Check for FASTQ files
            fastq_files = list(sample_dir.glob('*.fastq.gz'))
            if fastq_files and all(f.stat().st_size > 1000000 for f in fastq_files):
                unquantified.append(sample_id)
        
        return unquantified
    
    return find_unquantified_samples(config_path)


def assess_species_status(species: str) -> dict:
    """Comprehensive status assessment for a species."""
    info = SPECIES_INFO.get(species, {})
    species_dir = OUTPUT_DIR / species
    
    status = {
        'species': species,
        'name': info.get('name', species),
        'total': info.get('total', 0),
        'quantified': count_quantified(species),
        'unquantified_downloaded': count_downloaded_unquantified(species),
        'exists': species_dir.exists(),
        'metadata_exists': (species_dir / 'work' / 'metadata' / 'metadata.tsv').exists(),
        'index_exists': len(list((species_dir / 'work' / 'index').glob('*.idx'))) > 0 if (species_dir / 'work' / 'index').exists() else False,
        'merge_complete': (species_dir / 'merged' / 'merge' / 'merge.tpm.tsv').exists(),
        'curate_complete': len(list((species_dir / 'curate').glob('*.tsv'))) > 5 if (species_dir / 'curate').exists() else False,
        'sanity_complete': (species_dir / 'work' / 'sanity').exists(),
    }
    
    status['percent'] = (status['quantified'] / status['total'] * 100) if status['total'] > 0 else 0
    status['remaining'] = status['total'] - status['quantified']
    
    # Only ready for merge if quantification is complete (100% or downloads finished with failures)
    # Use 95% as threshold assuming some samples may fail to download
    status['ready_for_merge'] = status['quantified'] >= status['total'] * 0.95 and not status['merge_complete']
    status['ready_for_curate'] = status['merge_complete'] and not status['curate_complete']
    status['ready_for_sanity'] = status['curate_complete'] and not status['sanity_complete']
    
    return status


def print_assessment(species_list: Optional[list[str]] = None):
    """Print comprehensive assessment of all species."""
    if species_list is None:
        species_list = list(SPECIES_INFO.keys())
    
    logger.info("=" * 100)
    logger.info("WORKFLOW STATUS ASSESSMENT")
    logger.info("=" * 100)
    logger.info(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    total_quantified = 0
    total_samples = 0
    
    for species in species_list:
        status = assess_species_status(species)
        total_quantified += status['quantified']
        total_samples += status['total']
        
        logger.info(f"{'='*100}")
        logger.info(f"{status['name']} ({species})")
        logger.info(f"{'='*100}")
        logger.info(f"  Progress: {status['quantified']}/{status['total']} ({status['percent']:.1f}%)")
        logger.info(f"  Remaining: {status['remaining']}")
        
        if status['unquantified_downloaded']:
            logger.info(f"  ⚠️  Downloaded but unquantified: {len(status['unquantified_downloaded'])} samples")
        
        logger.info(f"\n  Status:")
        logger.info(f"    Metadata: {'✅' if status['metadata_exists'] else '❌'}")
        logger.info(f"    Index: {'✅' if status['index_exists'] else '❌'}")
        logger.info(f"    Quantification: {status['quantified']}/{status['total']} ({'✅' if status['quantified'] == status['total'] else '⏳'})")
        logger.info(f"    Merge: {'✅' if status['merge_complete'] else '⏳'}")
        logger.info(f"    Curate: {'✅' if status['curate_complete'] else '⏳'}")
        logger.info(f"    Sanity: {'✅' if status['sanity_complete'] else '⏳'}")
        
        logger.info(f"\n  Ready for:")
        if status['ready_for_merge']:
            logger.info(f"    🔹 merge (≥80% quantified)")
        if status['ready_for_curate']:
            logger.info(f"    🔹 curate (merge complete)")
        if status['ready_for_sanity']:
            logger.info(f"    🔹 sanity (curate complete)")
        if not any([status['ready_for_merge'], status['ready_for_curate'], status['ready_for_sanity']]):
            logger.info(f"    ✅ All steps complete or not ready")
        
        logger.info("")
    
    overall_pct = (total_quantified / total_samples * 100) if total_samples > 0 else 0
    logger.info(f"{'='*100}")
    logger.info(f"OVERALL: {total_quantified}/{total_samples} ({overall_pct:.1f}%)")
    logger.info(f"{'='*100}\n")


def cleanup_unquantified_for_species(species: str) -> tuple[int, int]:
    """Quantify downloaded samples and cleanup FASTQs using metainformant functions."""
    logger.info(f"\n{'='*80}")
    logger.info(f"CLEANUP: {species}")
    logger.info(f"{'='*80}")
    
    config_path = get_config_path(species)
    if not config_path:
        logger.error(f"  ❌ Config not found for {species}")
        return 0, 0
    
    # Use metainformant function to find unquantified samples
    unquantified = find_unquantified_samples(config_path)
    if not unquantified:
        logger.info(f"  ✅ No downloaded unquantified samples")
        return 0, 0
    
    logger.info(f"  Found {len(unquantified)} downloaded but unquantified samples")
    
    # Load config and get paths
    cfg = load_workflow_config(config_path)
    fastq_dir = Path(cfg.per_step.get("getfastq", {}).get("out_dir", cfg.work_dir / "fastq"))
    quant_dir = Path(cfg.per_step.get("quant", {}).get("out_dir", cfg.work_dir / "quant"))
    
    # Read metadata
    metadata_file = cfg.work_dir / "metadata" / "metadata.tsv"
    if not metadata_file.exists():
        metadata_file = cfg.work_dir / "metadata" / "metadata.filtered.tissue.tsv"
    
    if not metadata_file.exists():
        logger.error(f"  ❌ No metadata file found for {species}")
        return 0, len(unquantified)
    
    rows = list(read_delimited(metadata_file, delimiter="\t"))
    
    # Get quant params from config
    quant_params = dict(cfg.per_step.get("quant", {}))
    quant_params["out_dir"] = str(quant_dir.absolute())
    quant_params["threads"] = cfg.threads or 12
    
    # Inject index_dir if needed
    if "index_dir" not in quant_params and "index-dir" not in quant_params:
        index_dir = quant_dir.parent / "work" / "index"
        if index_dir.exists():
            quant_params["index_dir"] = str(index_dir.absolute())
    
    quantified_count = 0
    failed_count = 0
    cleaned_count = 0
    
    # Process each sample using metainformant functions
    for sample_id in unquantified:
        try:
            # Get metadata rows for this sample
            sample_rows = [row for row in rows if row.get("run") == sample_id]
            
            if not sample_rows:
                logger.warning(f"    ⚠️  {sample_id} not in metadata, skipping")
                failed_count += 1
                continue
            
            # Use metainformant quantify_sample function
            success, message, abundance_file = quantify_sample(
                sample_id=sample_id,
                metadata_rows=sample_rows,
                quant_params=quant_params,
                log_dir=cfg.log_dir or (cfg.work_dir / "logs"),
                step_name=f"quant_{species}_{sample_id}",
            )
            
            if success:
                logger.info(f"    ✅ Quantified: {sample_id}")
                quantified_count += 1
            else:
                logger.warning(f"    ❌ Failed: {sample_id} - {message}")
                failed_count += 1
            
            # Always cleanup FASTQs to free space
            delete_sample_fastqs(sample_id, fastq_dir)
            cleaned_count += 1
            logger.info(f"    🗑️  Cleaned: {sample_id}")
        
        except Exception as e:
            logger.warning(f"    ❌ Error: {sample_id} - {e}")
            failed_count += 1
            # Attempt cleanup
            try:
                delete_sample_fastqs(sample_id, fastq_dir)
                cleaned_count += 1
            except Exception:
                pass
    
    logger.info(f"\n  Summary: {quantified_count} quantified, {cleaned_count} cleaned, {failed_count} failed")
    return quantified_count, failed_count


def run_amalgkit_step(species: str, step: str, config_path: Path) -> bool:
    """Run a single amalgkit step for a species."""
    logger.info(f"\n  📦 Running {step} for {species}...")
    
    config = yaml.safe_load(config_path.read_text())
    work_dir = Path(config['work_dir'])
    log_dir = Path(config.get('log_dir', work_dir.parent / 'logs'))
    
    # Build command
    cmd = ['amalgkit', step]
    
    # Add common parameters
    if 'threads' in config:
        cmd.extend(['--threads', str(config['threads'])])
    
    # Add step-specific parameters
    step_config = config.get('steps', {}).get(step, {})
    for key, value in step_config.items():
        if value is None or value is False:
            continue
        if value is True:
            cmd.append(f'--{key}')
        else:
            cmd.extend([f'--{key}', str(value)])
    
    # Run
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f'{step}_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'
    
    try:
        with open(log_file, 'w') as log_f:
            result = subprocess.run(
                cmd,
                stdout=log_f,
                stderr=subprocess.STDOUT,
                cwd=str(work_dir),
                timeout=7200  # 2 hour timeout
            )
        
        if result.returncode == 0:
            logger.info(f"    ✅ {step} completed successfully")
            return True
        else:
            logger.warning(f"    ⚠️  {step} failed with code {result.returncode}")
            logger.warning(f"    Log: {log_file}")
            return False
    except subprocess.TimeoutExpired:
        logger.error(f"    ❌ {step} timed out after 2 hours")
        return False
    except Exception as e:
        logger.error(f"    ❌ {step} error: {e}")
        return False


def run_steps_for_species(species: str, steps: list[str]) -> dict:
    """Run specified steps for a species."""
    logger.info(f"\n{'='*80}")
    logger.info(f"RUNNING STEPS: {species}")
    logger.info(f"{'='*80}")
    logger.info(f"Steps: {', '.join(steps)}")
    
    config_path = CONFIG_DIR / f'amalgkit_{species}.yaml'
    if not config_path.exists():
        logger.error(f"  ❌ Config not found: {config_path}")
        return {'success': False, 'completed': []}
    
    results = {'success': True, 'completed': [], 'failed': []}
    
    for step in steps:
        if step not in VALID_STEPS:
            logger.warning(f"  ⚠️  Invalid step: {step}")
            continue
        
        success = run_amalgkit_step(species, step, config_path)
        if success:
            results['completed'].append(step)
        else:
            results['failed'].append(step)
            results['success'] = False
    
    logger.info(f"\n  Results: {len(results['completed'])} completed, {len(results['failed'])} failed")
    return results


def resume_downloads(species_list: Optional[list[str]] = None):
    """Resume incomplete downloads by restarting workflows."""
    if species_list is None:
        species_list = list(SPECIES_INFO.keys())
    
    logger.info(f"\n{'='*80}")
    logger.info("RESUMING DOWNLOADS")
    logger.info(f"{'='*80}")
    
    workflow_script = REPO_ROOT / 'scripts' / 'rna' / 'workflow_ena_integrated.py'
    
    for species in species_list:
        status = assess_species_status(species)
        
        if status['remaining'] == 0:
            logger.info(f"  ✅ {species}: All samples quantified")
            continue
        
        config_path = CONFIG_DIR / f'amalgkit_{species}.yaml'
        if not config_path.exists():
            logger.warning(f"  ⚠️  {species}: Config not found")
            continue
        
        logger.info(f"  🚀 {species}: Restarting workflow ({status['remaining']} remaining)...")
        
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        log_file = LOG_DIR / f'workflow_{species}_resumed_{timestamp}.log'
        
        cmd = [
            'nohup',
            sys.executable,
            str(workflow_script),
            '--config', str(config_path),
            '--batch-size', '12',
            '--threads', '12'
        ]
        
        try:
            with open(log_file, 'w') as log_f:
                process = subprocess.Popen(
                    cmd,
                    stdout=log_f,
                    stderr=subprocess.STDOUT,
                    cwd=str(REPO_ROOT),
                    start_new_session=True
                )
            logger.info(f"    ✅ Started (PID: {process.pid}, log: {log_file.name})")
        except Exception as e:
            logger.error(f"    ❌ Failed to start: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Unified Amalgkit Workflow Orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Full assessment
  %(prog)s --assess
  
  # Cleanup unquantified samples for all species
  %(prog)s --cleanup-unquantified
  
  # Run merge+curate+sanity for specific species
  %(prog)s --species cfloridanus sinvicta --steps merge curate sanity
  
  # Run downstream for all species ready
  %(prog)s --steps merge curate sanity --auto-species
  
  # Resume downloads
  %(prog)s --resume-downloads --species pbarbatus
        """
    )
    
    parser.add_argument('--assess', action='store_true',
                        help='Print comprehensive status assessment')
    parser.add_argument('--cleanup-unquantified', action='store_true',
                        help='Quantify downloaded samples and cleanup FASTQs')
    parser.add_argument('--resume-downloads', action='store_true',
                        help='Resume incomplete downloads')
    parser.add_argument('--steps', nargs='+', choices=VALID_STEPS,
                        help='Amalgkit steps to run')
    parser.add_argument('--species', nargs='+', choices=list(SPECIES_INFO.keys()),
                        help='Species to process (default: all)')
    parser.add_argument('--auto-species', action='store_true',
                        help='Auto-select species based on readiness for steps')
    
    args = parser.parse_args()
    
    # Determine species list
    if args.species:
        species_list = args.species
    elif args.auto_species and args.steps:
        # Auto-select species ready for the requested steps
        species_list = []
        for species in SPECIES_INFO.keys():
            status = assess_species_status(species)
            if any([
                'merge' in args.steps and status['ready_for_merge'],
                'curate' in args.steps and status['ready_for_curate'],
                'sanity' in args.steps and status['ready_for_sanity'],
            ]):
                species_list.append(species)
        logger.info(f"Auto-selected species: {', '.join(species_list)}")
    else:
        species_list = list(SPECIES_INFO.keys())
    
    # Execute requested actions
    if args.assess or not any([args.cleanup_unquantified, args.resume_downloads, args.steps]):
        print_assessment(species_list)
        return 0
    
    exit_code = 0
    
    if args.cleanup_unquantified:
        logger.info("\n" + "="*100)
        logger.info("CLEANUP UNQUANTIFIED SAMPLES")
        logger.info("="*100)
        
        total_quantified = 0
        total_failed = 0
        
        for species in species_list:
            quantified, failed = cleanup_unquantified_for_species(species)
            total_quantified += quantified
            total_failed += failed
        
        logger.info(f"\n{'='*100}")
        logger.info(f"CLEANUP COMPLETE: {total_quantified} quantified, {total_failed} failed")
        logger.info(f"{'='*100}")
        
        if total_failed > 0:
            exit_code = 1
    
    if args.resume_downloads:
        resume_downloads(species_list)
    
    if args.steps:
        logger.info("\n" + "="*100)
        logger.info("RUNNING AMALGKIT STEPS")
        logger.info("="*100)
        
        all_success = True
        for species in species_list:
            results = run_steps_for_species(species, args.steps)
            if not results['success']:
                all_success = False
                exit_code = 1
        
        logger.info(f"\n{'='*100}")
        if all_success:
            logger.info("✅ ALL STEPS COMPLETED SUCCESSFULLY")
        else:
            logger.info("⚠️  SOME STEPS FAILED - CHECK LOGS")
        logger.info(f"{'='*100}")
    
    return exit_code


if __name__ == '__main__':
    sys.exit(main())

