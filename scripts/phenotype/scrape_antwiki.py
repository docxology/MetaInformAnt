#!/usr/bin/env python3
"""Command-line script for scraping AntWiki species pages.

Scrapes AntWiki (antwiki.org) to extract comprehensive species phenotype data
from all sections of species pages.

Usage:
    python3 scripts/phenotype/scrape_antwiki.py --species Camponotus_pennsylvanicus
    python3 scripts/phenotype/scrape_antwiki.py --all --limit 10
    python3 scripts/phenotype/scrape_antwiki.py --all --delay 3.0 --output output/phenotype/antwiki/
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from metainformant.phenotype.scraper import AntWikiScraper, AntWikiScraperConfig, load_scraper_config


def main() -> int:
    """Main entry point for AntWiki scraper CLI."""
    parser = argparse.ArgumentParser(
        description="Scrape AntWiki species pages for phenotype data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument(
        "--species",
        type=str,
        help="Scrape single species (URL-safe format, e.g., Camponotus_pennsylvanicus)",
    )

    parser.add_argument(
        "--all",
        action="store_true",
        help="Scrape all species pages",
    )

    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output directory (default: output/phenotype/antwiki/)",
    )

    parser.add_argument(
        "--delay",
        type=float,
        default=None,
        help="Delay between requests in seconds (overrides config)",
    )

    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume from last checkpoint",
    )

    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Limit number of species to scrape (for testing)",
    )

    parser.add_argument(
        "--check-robots",
        action="store_true",
        default=None,
        help="Check robots.txt before scraping (default: True)",
    )

    parser.add_argument(
        "--config",
        type=str,
        default=None,
        help="Path to config file (default: config/phenotype/antwiki_scraper.yaml)",
    )

    args = parser.parse_args()

    # Validate arguments
    if not args.species and not args.all:
        parser.error("Must specify either --species or --all")

    if args.species and args.all:
        parser.error("Cannot specify both --species and --all")

    # Load configuration
    config_path = Path(args.config) if args.config else None
    config = load_scraper_config(config_path)

    # Override config with CLI arguments
    if args.delay is not None:
        config.delay_seconds = args.delay

    if args.output is not None:
        config.output_dir = Path(args.output)

    if args.check_robots is not None:
        config.check_robots = args.check_robots

    # Initialize scraper
    scraper = AntWikiScraper(config)

    try:
        if args.species:
            # Scrape single species
            print(f"Scraping species: {args.species}")
            data = scraper.scrape_species_page(args.species)

            # Save to output directory
            output_dir = config.output_dir
            output_dir.mkdir(parents=True, exist_ok=True)
            species_dir = output_dir / "species"
            species_dir.mkdir(parents=True, exist_ok=True)

            from metainformant.core.paths import sanitize_filename
            from metainformant.core.io import dump_json

            safe_name = sanitize_filename(args.species)
            output_file = species_dir / f"{safe_name}.json"
            dump_json(data, output_file, indent=2)

            print(f"Successfully scraped {args.species}")
            print(f"Output saved to: {output_file}")
            print(f"Found {len(data.get('traits', []))} traits, {len(data.get('measurements', {}))} measurements")

            return 0

        elif args.all:
            # Scrape all species
            print(f"Scraping all species (limit: {args.limit or 'none'})")
            print(f"Output directory: {config.output_dir}")
            print(f"Delay between requests: {config.delay_seconds} seconds")
            print(f"Check robots.txt: {config.check_robots}")

            stats = scraper.scrape_all_species(
                output_dir=config.output_dir,
                limit=args.limit,
                resume=args.resume,
            )

            print("\nScraping complete!")
            print(f"Total species: {stats['total']}")
            print(f"Completed: {stats['completed']}")
            print(f"Failed: {stats['failed']}")

            if stats["failed"] > 0:
                print(f"\nFailed species: {', '.join(stats['failed_species'][:10])}")
                if len(stats["failed_species"]) > 10:
                    print(f"... and {len(stats['failed_species']) - 10} more")

            return 0 if stats["failed"] == 0 else 1

    except KeyboardInterrupt:
        print("\nScraping interrupted by user")
        return 130
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback

        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())

