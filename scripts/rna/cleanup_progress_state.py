"""Clean up progress_state.json to only include current species from config files.

This script:
1. Discovers current species configs (only 3: camponotus_floridanus, harpegnathos_saltator, pogonomyrmex_barbatus)
2. For each species, loads actual metadata to get correct sample IDs
3. Rebuilds progress_state.json with only those 3 species
4. Preserves actual progress (completed samples) from filesystem analysis
"""

from __future__ import annotations

import sys
from pathlib import Path

# Add src to path
repo_root = Path(__file__).parent.parent.parent.resolve()
src_path = repo_root / "src"
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

from metainformant.rna.workflow import load_workflow_config
from metainformant.rna.progress_tracker import get_tracker, ProgressTracker
from metainformant.core.io import read_delimited
from metainformant.core.logging import get_logger

logger = get_logger("cleanup_progress_state")


def get_species_id_from_config(config_path: Path) -> str:
    """Extract species ID from config file.
    
    Args:
        config_path: Path to config file
        
    Returns:
        Species ID (e.g., "camponotus_floridanus")
    """
    try:
        cfg = load_workflow_config(config_path)
        # Try to get from species_list first
        if cfg.species_list:
            # Use first species, convert to lowercase with underscores
            species_id = cfg.species_list[0].lower().replace(" ", "_")
            return species_id
    except Exception:
        pass
    
    # Fallback: extract from filename
    stem = config_path.stem
    if stem.startswith("amalgkit_"):
        return stem.replace("amalgkit_", "")
    return stem


def discover_current_species() -> list[tuple[str, Path]]:
    """Discover current species config files.
    
    Returns:
        List of (species_name, config_path) tuples
    """
    config_dir = repo_root / "config" / "amalgkit"
    if not config_dir.exists():
        logger.error(f"Config directory not found: {config_dir}")
        return []
    
    species_configs = []
    for config_file in sorted(config_dir.glob("amalgkit_*.yaml")):
        if "template" in config_file.stem.lower() or "test" in config_file.stem.lower():
            continue
        
        species_id = get_species_id_from_config(config_file)
        species_configs.append((species_id, config_file))
        logger.info(f"Found species config: {species_id} -> {config_file.name}")
    
    return species_configs


def get_sample_ids_from_metadata(config_path: Path) -> list[str]:
    """Get all sample IDs from metadata file.
    
    Args:
        config_path: Path to species config file
        
    Returns:
        List of sample IDs (run accessions)
    """
    try:
        cfg = load_workflow_config(config_path)
        metadata_file = cfg.work_dir / "metadata" / "metadata.tsv"
        if not metadata_file.exists():
            metadata_file = cfg.work_dir / "metadata" / "metadata.filtered.tissue.tsv"
        
        if not metadata_file.exists():
            logger.warning(f"No metadata file found for {config_path}")
            return []
        
        rows = list(read_delimited(metadata_file, delimiter="\t"))
        sample_ids = [row.get("run") for row in rows if row.get("run") if row.get("run")]
        return sample_ids
    except Exception as e:
        logger.error(f"Error reading metadata for {config_path}: {e}")
        return []


def check_sample_completed(species_id: str, sample_id: str, config_path: Path) -> bool:
    """Check if a sample has been completed (quantified).
    
    Args:
        species_id: Species ID
        sample_id: Sample ID
        config_path: Path to config file
        
    Returns:
        True if sample is completed (has abundance.tsv)
    """
    try:
        cfg = load_workflow_config(config_path)
        quant_params = dict(cfg.per_step.get("quant", {}))
        quant_dir = Path(quant_params.get("out_dir", cfg.work_dir / "quant"))
        abundance_file = quant_dir / sample_id / "abundance.tsv"
        return abundance_file.exists()
    except Exception:
        return False


def main() -> int:
    """Main entry point."""
    logger.info("=" * 80)
    logger.info("Cleaning up progress_state.json")
    logger.info("=" * 80)
    
    # Discover current species
    species_configs = discover_current_species()
    if not species_configs:
        logger.error("No species configs found!")
        return 1
    
    logger.info(f"\nFound {len(species_configs)} current species:")
    for species_id, config_path in species_configs:
        logger.info(f"  - {species_id}")
    
    # Backup old state file if it exists
    state_file = repo_root / "output" / "amalgkit" / "progress_state.json"
    if state_file.exists():
        backup_file = state_file.with_suffix(".json.backup")
        logger.info(f"\nBacking up existing state file to {backup_file.name}")
        import shutil
        shutil.copy2(state_file, backup_file)
    
    # Create new tracker instance (fresh state, not loading old one)
    # Use a temporary state file first, then move it
    temp_state_file = state_file.with_suffix(".json.tmp")
    tracker = ProgressTracker(state_file=temp_state_file)
    
    # Initialize tracker for each species with correct sample IDs
    total_samples = 0
    completed_samples = 0
    
    for species_id, config_path in species_configs:
        logger.info(f"\nProcessing {species_id}...")
        
        # Get sample IDs from metadata
        sample_ids = get_sample_ids_from_metadata(config_path)
        if not sample_ids:
            logger.warning(f"  No sample IDs found for {species_id}, skipping")
            continue
        
        logger.info(f"  Found {len(sample_ids)} samples in metadata")
        
        # Initialize species in tracker
        tracker.initialize_species(species_id, len(sample_ids), sample_ids)
        
        # Check which samples are already completed
        species_completed = 0
        for sample_id in sample_ids:
            if check_sample_completed(species_id, sample_id, config_path):
                # Mark as completed
                tracker.on_quant_complete(species_id, sample_id)
                tracker.on_delete_complete(species_id, sample_id)
                species_completed += 1
                completed_samples += 1
        
        total_samples += len(sample_ids)
        logger.info(f"  {species_completed} samples already completed for {species_id}")
    
    # State is already saved by tracker methods (auto-save on state changes)
    # Move temp file to final location
    if temp_state_file.exists():
        import shutil
        shutil.move(temp_state_file, state_file)
        logger.info(f"\n✓ Progress state cleaned and saved to {state_file.name}")
    else:
        logger.warning("Warning: State file was not created")
    
    logger.info(f"  Total species: {len(species_configs)}")
    logger.info(f"  Total samples: {total_samples}")
    logger.info(f"  Completed samples: {completed_samples}")
    
    # Skip dashboard update and visualization generation to avoid hanging
    # These will be generated automatically when batch_download_species.py runs
    logger.info("\nNote: Dashboard and visualization will be generated automatically")
    logger.info("  when batch_download_species.py runs with tracker integration")
    
    logger.info("\n" + "=" * 80)
    logger.info("Cleanup complete!")
    logger.info("=" * 80)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

