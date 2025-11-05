#!/usr/bin/env python3
"""Initialize progress tracking for all species from existing output directories.

This script scans the output/amalgkit directory and initializes the progress tracker
with the current state of all species, so tracking can begin even if workflows
were started before the tracker was integrated.
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.rna.progress_tracker import get_tracker
from metainformant.rna.monitoring import analyze_species_status
from metainformant.rna.workflow import load_workflow_config
from metainformant.core.io import read_delimited
from metainformant.core.logging import get_logger

logger = get_logger("init_tracking")

def discover_species_from_output(output_dir: Path = Path("output/amalgkit")) -> list[tuple[str, Path]]:
    """Discover species from output directories and find their configs."""
    species_configs = []
    
    config_dir = Path("config/amalgkit")
    
    for species_dir in sorted(output_dir.iterdir()):
        if not species_dir.is_dir():
            continue
        
        species_id = species_dir.name
        
        # Try to find matching config
        config_path = config_dir / f"amalgkit_{species_id}.yaml"
        if not config_path.exists():
            # Try alternative naming
            config_path = config_dir / f"amalgkit_{species_id.replace('_', '').lower()}.yaml"
            if not config_path.exists():
                continue
        
        species_configs.append((species_id, config_path))
    
    return species_configs


def initialize_all_species(tracker):
    """Initialize tracker for all species found in output directory."""
    species_configs = discover_species_from_output()
    
    if not species_configs:
        logger.warning("No species found in output/amalgkit/")
        return
    
    logger.info(f"Found {len(species_configs)} species to initialize")
    
    for species_id, config_path in species_configs:
        try:
            cfg = load_workflow_config(config_path)
            metadata_file = cfg.work_dir / "metadata" / "metadata.tsv"
            if not metadata_file.exists():
                metadata_file = cfg.work_dir / "metadata" / "metadata.filtered.tissue.tsv"
            
            if not metadata_file.exists():
                logger.warning(f"  {species_id}: No metadata file found")
                continue
            
            rows = list(read_delimited(metadata_file, delimiter="\t"))
            sample_ids = [row.get("run") for row in rows if row.get("run")]
            
            if not sample_ids:
                logger.warning(f"  {species_id}: No samples in metadata")
                continue
            
            # Initialize tracker
            tracker.initialize_species(species_id, len(sample_ids), sample_ids)
            
            # Analyze current status and populate tracker
            status = analyze_species_status(config_path)
            categories = status.get("categories", {})
            
            # Mark completed samples
            for sample_id in categories.get("quantified_and_deleted", []):
                tracker.on_quant_complete(species_id, sample_id)
                tracker.on_delete_complete(species_id, sample_id)
            
            # Mark samples that need delete
            for sample_id in categories.get("quantified_not_deleted", []):
                tracker.on_quant_complete(species_id, sample_id)
            
            # Mark downloading samples
            for sample_id in categories.get("downloading", []):
                tracker.on_download_start(species_id, sample_id)
            
            # Mark failed downloads
            for sample_id in categories.get("failed_download", []):
                tracker.on_download_failed(species_id, sample_id)
            
            # Mark samples that have FASTQ but need quant (they're in needs_quant)
            # These will be detected by analyze_species_status as "failed_download" if they have files
            # but we need to check if they're actually quantified
            
            logger.info(f"  ✅ {species_id}: {len(sample_ids)} samples initialized")
            
        except Exception as e:
            logger.error(f"  ❌ {species_id}: Failed to initialize - {e}")
    
    # Update dashboard
    tracker.update_dashboard()
    logger.info(f"Dashboard updated: {tracker.dashboard_file}")


def main():
    """Main entry point."""
    tracker = get_tracker()
    
    logger.info("=" * 80)
    logger.info("INITIALIZING PROGRESS TRACKING FOR ALL SPECIES")
    logger.info("=" * 80)
    logger.info("")
    
    initialize_all_species(tracker)
    
    logger.info("")
    logger.info("=" * 80)
    logger.info("INITIALIZATION COMPLETE")
    logger.info("=" * 80)
    logger.info(f"Dashboard: {tracker.dashboard_file}")
    logger.info(f"State: {tracker.state_file}")
    logger.info("")


if __name__ == "__main__":
    main()

