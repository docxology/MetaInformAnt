#!/usr/bin/env python3
"""Command-line script for scraping AntWiki species pages.

Scrapes AntWiki (antwiki.org) to extract comprehensive species phenotype data
from all sections of species pages.

Usage:
    python3 scripts/phenotype/scrape_antwiki.py --species Camponotus_pennsylvanicus
    python3 scripts/phenotype/scrape_antwiki.py --species Formica_rufa --delay 3.0 --output output/phenotype/antwiki/

Notes:
    robots.txt is honoured by default: pages it disallows are refused, and a Crawl-delay larger than --delay
    replaces --delay. Pages are fetched with requests (not cloudscraper), so a Cloudflare challenge is not bypassed.
    --all needs the AntWiki species index, whose retrieval is not implemented, so it exits with an error; scrape
    named species one at a time with --species.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from metainformant.core.utils.errors import ValidationError
from metainformant.phenotype.data.scraper import AntWikiScraper, RobotsDisallowedError, load_scraper_config


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
        help="Scrape every species in the AntWiki species index (index retrieval is not implemented: "
        "this exits with an error; use --species)",
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
        help="Delay between requests in seconds (overrides config; a larger robots.txt Crawl-delay still wins)",
    )

    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume from last checkpoint (not supported by the current scraper; exits with an error)",
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
        help="Check robots.txt before scraping (default: True; opt out with check_robots=false in a --config file)",
    )

    parser.add_argument(
        "--config",
        type=str,
        default=None,
        help="Path to a JSON file of AntWikiScraperConfig settings (default: built-in defaults)",
    )

    args = parser.parse_args()

    # Validate arguments
    if not args.species and not args.all:
        parser.error("Must specify either --species or --all")

    if args.species and args.all:
        parser.error("Cannot specify both --species and --all")

    if args.resume:
        parser.error("--resume is not supported: AntWikiScraper.scrape_all_species has no checkpointing")

    # Load configuration
    config_path = Path(args.config) if args.config else None
    config = load_scraper_config(config_path)

    # Override config with CLI arguments (delay_seconds and delay_between_requests are the same setting)
    if args.delay is not None:
        try:
            config.delay_between_requests = args.delay
        except ValidationError as e:
            parser.error(f"--delay: {e}")

    # Output directory precedence: --output, then the config file, then the documented default
    config.output_dir = Path(args.output or config.output_dir or "output/phenotype/antwiki")

    if args.check_robots is not None:
        config.check_robots = args.check_robots

    # Initialize scraper
    scraper = AntWikiScraper(config)

    try:
        if args.species:
            # Scrape single species
            print(f"Scraping species: {args.species}")
            data = scraper.scrape_species_page(args.species)
            if data is None:
                print(f"Failed to scrape {args.species}: page could not be parsed", file=sys.stderr)
                return 1

            # Save to output directory
            output_dir = config.output_dir
            output_dir.mkdir(parents=True, exist_ok=True)
            species_dir = output_dir / "species"
            species_dir.mkdir(parents=True, exist_ok=True)

            from metainformant.core.io import dump_json
            from metainformant.core.io.paths import sanitize_filename

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
            print(f"Delay between requests: {config.delay_between_requests} seconds")
            if config.check_robots:
                print("Check robots.txt: True (disallowed pages are refused; a larger Crawl-delay replaces the delay)")
            else:
                print("Check robots.txt: False (opted out: robots.txt is not consulted)")

            try:
                stats = scraper.scrape_all_species(
                    output_dir=config.output_dir,
                    limit=args.limit,
                )
            except NotImplementedError as e:
                print(f"Error: {e}", file=sys.stderr)
                print("Scrape named species one at a time with --species.", file=sys.stderr)
                return 2

            print("\nScraping complete!")
            print(f"Total species: {stats['total']}")
            print(f"Completed: {stats['completed']}")
            print(f"Failed: {stats['failed']}")

            failed_species = stats.get("failed_species", [])
            if failed_species:
                print(f"\nFailed species: {', '.join(failed_species[:10])}")
                if len(failed_species) > 10:
                    print(f"... and {len(failed_species) - 10} more")

            return 0 if stats["failed"] == 0 else 1

    except KeyboardInterrupt:
        print("\nScraping interrupted by user")
        return 130
    except RobotsDisallowedError as e:
        print(f"Refused: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback

        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
