#!/usr/bin/env python3
"""Comprehensive progress assessment and dashboard generation.

This script assesses the current state of all RNA-seq workflows and generates
a comprehensive progress dashboard, even if workflows were started before
the progress tracker was integrated.
"""

import sys
from pathlib import Path
from datetime import datetime

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.rna.progress_tracker import get_tracker
from metainformant.rna.monitoring import analyze_species_status
from metainformant.rna.workflow import load_workflow_config
from metainformant.core.io import read_delimited
from metainformant.core.logging import get_logger

logger = get_logger("assess_progress")


def discover_all_species() -> list[tuple[str, Path]]:
    """Discover all species from config files and output directories."""
    config_dir = Path("config/amalgkit")
    output_dir = Path("output/amalgkit")
    
    species_configs = []
    
    # Find configs
    for config_file in sorted(config_dir.glob("amalgkit_*.yaml")):
        if "template" in config_file.stem.lower() or "test" in config_file.stem.lower():
            continue
        
        species_id = config_file.stem.replace("amalgkit_", "")
        species_configs.append((species_id, config_file))
    
    return species_configs


def assess_and_update_all_species(tracker):
    """Assess current state and update tracker for all species."""
    species_configs = discover_all_species()
    
    if not species_configs:
        logger.warning("No species configurations found")
        return
    
    logger.info(f"Assessing {len(species_configs)} species...")
    
    for species_id, config_path in species_configs:
        try:
            cfg = load_workflow_config(config_path)
            metadata_file = cfg.work_dir / "metadata" / "metadata.tsv"
            if not metadata_file.exists():
                metadata_file = cfg.work_dir / "metadata" / "metadata.filtered.tissue.tsv"
            
            if not metadata_file.exists():
                logger.debug(f"  {species_id}: No metadata file")
                continue
            
            rows = list(read_delimited(metadata_file, delimiter="\t"))
            sample_ids = [row.get("run") for row in rows if row.get("run")]
            
            if not sample_ids:
                logger.debug(f"  {species_id}: No samples in metadata")
                continue
            
            # Initialize if not already initialized
            if species_id not in tracker.state:
                tracker.initialize_species(species_id, len(sample_ids), sample_ids)
            
            # Analyze current status from filesystem
            status = analyze_species_status(config_path)
            categories = status.get("categories", {})
            
            # Sync tracker state with filesystem reality
            # Mark completed samples
            for sample_id in categories.get("quantified_and_deleted", []):
                if sample_id not in tracker.state[species_id]["completed"]:
                    # Move through states if needed
                    tracker.state[species_id]["need_download"].discard(sample_id)
                    tracker.state[species_id]["ongoing_download"].discard(sample_id)
                    tracker.state[species_id]["failed_download"].discard(sample_id)
                    tracker.state[species_id]["needs_quant"].discard(sample_id)
                    tracker.state[species_id]["needs_delete"].discard(sample_id)
                    tracker.state[species_id]["completed"].add(sample_id)
            
            # Mark samples that need delete
            for sample_id in categories.get("quantified_not_deleted", []):
                if sample_id not in tracker.state[species_id]["needs_delete"]:
                    tracker.state[species_id]["need_download"].discard(sample_id)
                    tracker.state[species_id]["ongoing_download"].discard(sample_id)
                    tracker.state[species_id]["failed_download"].discard(sample_id)
                    tracker.state[species_id]["needs_quant"].discard(sample_id)
                    tracker.state[species_id]["needs_delete"].add(sample_id)
            
            # Mark downloading samples
            for sample_id in categories.get("downloading", []):
                if sample_id not in tracker.state[species_id]["ongoing_download"]:
                    tracker.state[species_id]["need_download"].discard(sample_id)
                    tracker.state[species_id]["ongoing_download"].add(sample_id)
            
            # Mark failed downloads
            for sample_id in categories.get("failed_download", []):
                if sample_id not in tracker.state[species_id]["failed_download"]:
                    tracker.state[species_id]["need_download"].discard(sample_id)
                    tracker.state[species_id]["ongoing_download"].discard(sample_id)
                    tracker.state[species_id]["failed_download"].add(sample_id)
            
            # Samples that have FASTQ but aren't quantified go to needs_quant
            for sample_id in categories.get("failed_download", []):
                # Check if it's actually quantified (might be mis-categorized)
                quant_dir = Path(cfg.per_step.get("quant", {}).get("out_dir", cfg.work_dir / "quant"))
                if (quant_dir / sample_id / "abundance.tsv").exists():
                    # Actually quantified, move to needs_delete
                    tracker.state[species_id]["failed_download"].discard(sample_id)
                    tracker.state[species_id]["needs_delete"].add(sample_id)
            
            tracker.state[species_id]["last_update"] = datetime.now()
            
            state = tracker.get_species_state(species_id)
            logger.info(f"  ✅ {species_id}: {state['completed']}/{state['total']} completed ({state['completed']/max(1,state['total'])*100:.1f}%)")
            
        except Exception as e:
            logger.error(f"  ❌ {species_id}: Failed to assess - {e}")
    
    # Save state and update dashboard
    tracker._save_state()
    tracker.update_dashboard()
    logger.info(f"Dashboard updated: {tracker.dashboard_file}")


def main():
    """Main entry point."""
    from datetime import datetime
    
    tracker = get_tracker()
    
    logger.info("=" * 100)
    logger.info("COMPREHENSIVE PROGRESS ASSESSMENT")
    logger.info("=" * 100)
    logger.info(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info("")
    
    assess_and_update_all_species(tracker)
    
    # Show summary
    all_state = tracker.get_all_species_state()
    total_samples = sum(s.get("total", 0) for s in all_state.values())
    total_completed = sum(s.get("completed", 0) for s in all_state.values())
    
    logger.info("")
    logger.info("=" * 100)
    logger.info("SUMMARY")
    logger.info("=" * 100)
    logger.info(f"Species tracked: {len(all_state)}")
    logger.info(f"Total samples: {total_samples}")
    logger.info(f"Completed: {total_completed} ({total_completed/max(1,total_samples)*100:.1f}%)")
    logger.info(f"Dashboard: {tracker.dashboard_file}")
    logger.info("=" * 100)


if __name__ == "__main__":
    main()

