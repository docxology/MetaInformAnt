#!/usr/bin/env python3
"""Species discovery and configuration generator.

This is a thin wrapper that calls methods from metainformant.rna.discovery
to discover species with RNA-seq data and generate amalgkit configurations.

Usage:
    # Discover ant species with RNA-seq data
    python3 scripts/rna/discover_species.py --output config/amalgkit/

    # Generate config for specific species
    python3 scripts/rna/discover_species.py --species "Camponotus floridanus" --output config/amalgkit/
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Import setup utilities (must be before other imports)
sys.path.insert(0, str(Path(__file__).parent))
from _setup_utils import check_environment_or_exit, ensure_venv_activated

# Auto-setup and activate venv
ensure_venv_activated(auto_setup=True)
check_environment_or_exit(auto_setup=True)

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.core.utils.logging import get_logger
from metainformant.rna.discovery import (
    generate_config_yaml,
    get_genome_info,
    search_species_with_rnaseq,
)

logger = get_logger("discover_species")


def main() -> int:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Discover species with RNA-seq data and generate configurations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output directory for generated config files",
    )
    parser.add_argument(
        "--species",
        type=str,
        help="Specific species to generate config for (scientific name)",
    )
    parser.add_argument(
        "--search-query",
        type=str,
        help="Custom NCBI Entrez search query (default: ant species with RNA-seq)",
    )
    parser.add_argument(
        "--min-samples",
        type=int,
        default=1,
        help="Minimum number of RNA-seq samples required (default: 1)",
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        help="Repository root directory for paths in config (default: auto-detect)",
    )

    args = parser.parse_args()

    output_dir = args.output.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    # Determine repo root
    if args.repo_root:
        repo_root = args.repo_root.resolve()
    else:
        repo_root = Path(__file__).parent.parent.parent.resolve()

    # Single species mode
    if args.species:
        logger.info(f"Generating config for {args.species}...")

        # Search for this species
        search_query = f'"{args.species}"[Organism] AND RNA-Seq[Strategy] AND Illumina[Platform]'
        species_data = search_species_with_rnaseq(search_query)

        if not species_data:
            logger.error(f"No RNA-seq data found for {args.species}")
            return 1

        species_info = species_data.get(args.species)
        if not species_info or species_info["sample_count"] < args.min_samples:
            logger.error(
                f"Insufficient samples for {args.species}: {species_info.get('sample_count', 0)} < {args.min_samples}"
            )
            return 1

        # Get genome info
        taxonomy_id = species_info.get("taxonomy_id")
        genome_info = None
        if taxonomy_id:
            genome_info = get_genome_info(taxonomy_id, args.species)

        # Generate config
        yaml_content = generate_config_yaml(
            args.species,
            species_info,
            genome_info,
            repo_root=repo_root,
        )

        # Write config file
        species_slug = args.species.replace(" ", "_").lower()
        config_file = output_dir / f"amalgkit_{species_slug}.yaml"
        config_file.write_text(yaml_content)
        logger.info(f"✅ Generated config: {config_file}")

        return 0

    # Discovery mode
    logger.info("Discovering species with RNA-seq data...")

    # Default search query for ants
    if not args.search_query:
        args.search_query = (
            'txid7389[Organism:exp] AND "RNA-Seq"[Strategy] AND "Illumina"[Platform] AND "public"[Access]'
        )

    species_data = search_species_with_rnaseq(args.search_query)

    if not species_data:
        logger.warning("No species found matching search criteria")
        return 1

    logger.info(f"Found {len(species_data)} species")

    # Filter by minimum samples
    filtered_species = {name: data for name, data in species_data.items() if data["sample_count"] >= args.min_samples}

    logger.info(f"Filtered to {len(filtered_species)} species with ≥{args.min_samples} samples")

    # Generate configs
    generated = 0
    for species_name, data in filtered_species.items():
        try:
            # Get genome info
            taxonomy_id = data.get("taxonomy_id")
            genome_info = None
            if taxonomy_id:
                genome_info = get_genome_info(taxonomy_id, species_name)

            # Generate config
            yaml_content = generate_config_yaml(
                species_name,
                data,
                genome_info,
                repo_root=repo_root,
            )

            # Write config file
            species_slug = species_name.replace(" ", "_").lower()
            config_file = output_dir / f"amalgkit_{species_slug}.yaml"
            config_file.write_text(yaml_content)
            logger.info(f"✅ Generated: {config_file.name}")
            generated += 1

        except Exception as e:
            logger.error(f"❌ Failed to generate config for {species_name}: {e}")

    logger.info(f"Generated {generated}/{len(filtered_species)} config files")
    return 0 if generated > 0 else 1


if __name__ == "__main__":
    sys.exit(main())
