#!/usr/bin/env python3
"""Test AntWiki scraper with cloudscraper.

This script tests the scraper's ability to fetch pages from AntWiki
using cloudscraper to bypass Cloudflare protection.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.phenotype.scraper import AntWikiScraper, AntWikiScraperConfig


def test_scraper() -> None:
    """Test scraper with a single species page."""
    print("Testing AntWiki scraper with cloudscraper...")

    # Create config with reasonable delays
    config = AntWikiScraperConfig(
        base_url="https://www.antwiki.org/wiki/",
        delay_seconds=3.0,  # Longer delay to be respectful
        user_agent="Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36",
        timeout_seconds=30,
        max_retries=2,
        check_robots=True,
        output_dir=Path("output/phenotype/test_scraper"),
    )

    scraper = AntWikiScraper(config)

    # Test with a known species
    species_name = "Camponotus_pennsylvanicus"
    print(f"\nAttempting to scrape: {species_name}")

    try:
        data = scraper.scrape_species_page(species_name)
        print(f"\n✓ Successfully scraped {species_name}!")
        print(f"  - Species: {data.get('species', 'N/A')}")
        print(f"  - Measurements: {len(data.get('measurements', {}))} found")
        print(f"  - Traits: {len(data.get('traits', []))} found")
        print(f"  - Description length: {len(data.get('description', ''))} chars")
        print(f"  - Taxonomy: {data.get('taxonomy', {})}")
        return True
    except Exception as e:
        print(f"\n✗ Failed to scrape: {e}")
        import traceback

        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = test_scraper()
    sys.exit(0 if success else 1)
