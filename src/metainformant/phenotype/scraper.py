"""AntWiki web scraping utilities.

This module provides comprehensive web scraping functionality for AntWiki (antwiki.org)
to extract species phenotype data from all sections of species pages.
"""

from __future__ import annotations

import json
import os
import random
import time
import urllib.parse
import urllib.robotparser
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import requests
from bs4 import BeautifulSoup

from ..core.config import apply_env_overrides, load_mapping_from_file
from ..core.errors import NetworkError, retry_with_backoff
from ..core.io import dump_json, ensure_directory, load_json, read_jsonl
from ..core.logging import get_logger
from ..core.paths import expand_and_resolve, sanitize_filename

logger = get_logger(__name__)

# Try to import cloudscraper for Cloudflare bypass
try:
    import cloudscraper

    CLOUDSCRAPER_AVAILABLE = True
except ImportError:
    CLOUDSCRAPER_AVAILABLE = False
    cloudscraper = None  # type: ignore


@dataclass
class AntWikiScraperConfig:
    """Configuration for AntWiki scraper."""

    base_url: str
    delay_seconds: float
    user_agent: str
    timeout_seconds: int
    max_retries: int
    check_robots: bool
    output_dir: Path


def load_scraper_config(config_path: Path | None = None) -> AntWikiScraperConfig:
    """Load scraper configuration from file with environment overrides.

    Args:
        config_path: Path to config file (default: config/phenotype/antwiki_scraper.yaml)

    Returns:
        AntWikiScraperConfig instance
    """
    if config_path is None:
        repo_root = Path(__file__).parent.parent.parent.parent
        config_path = repo_root / "config" / "phenotype" / "antwiki_scraper.yaml"

    if config_path.exists():
        raw_config = load_mapping_from_file(config_path)
        raw_config = apply_env_overrides(raw_config, prefix="PHEN")
    else:
        raw_config = {}

    # Apply environment variable overrides
    delay = float(raw_config.get("delay_seconds", os.getenv("PHEN_SCRAPE_DELAY", "2.0")))
    base_url = raw_config.get("base_url", "https://www.antwiki.org/wiki/")
    user_agent = raw_config.get(
        "user_agent", "METAINFORMANT/0.2.0 (https://github.com/metainformant; research use)"
    )
    timeout = int(raw_config.get("timeout_seconds", 30))
    max_retries = int(raw_config.get("max_retries", 3))
    check_robots = raw_config.get("check_robots", True)
    output_dir = Path(raw_config.get("output_dir", "output/phenotype/antwiki/"))

    return AntWikiScraperConfig(
        base_url=base_url,
        delay_seconds=delay,
        user_agent=user_agent,
        timeout_seconds=timeout,
        max_retries=max_retries,
        check_robots=check_robots,
        output_dir=expand_and_resolve(output_dir),
    )


class AntWikiScraper:
    """Scraper for AntWiki species pages.

    Provides comprehensive scraping of all sections from AntWiki species pages,
    including measurements, traits, taxonomy, distribution, and descriptions.
    """

    def __init__(self, config: AntWikiScraperConfig | None = None):
        """Initialize AntWiki scraper.

        Args:
            config: Scraper configuration (loads from file if None)
        """
        if config is None:
            config = load_scraper_config()
        self.config = config
        
        # Use cloudscraper if available (handles Cloudflare), otherwise fallback to requests
        if CLOUDSCRAPER_AVAILABLE:
            logger.info("Using cloudscraper for Cloudflare bypass")
            self.session = cloudscraper.create_scraper(
                browser={
                    'browser': 'chrome',
                    'platform': 'linux',
                    'desktop': True
                }
            )
        else:
            logger.warning("cloudscraper not available, using requests (may fail on Cloudflare-protected sites)")
            self.session = requests.Session()
        
        # Use more realistic browser-like headers to avoid 403 errors
        self.session.headers.update({
            "User-Agent": self.config.user_agent,
            "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8",
            "Accept-Language": "en-US,en;q=0.9",
            "Accept-Encoding": "gzip, deflate, br",
            "Connection": "keep-alive",
            "Upgrade-Insecure-Requests": "1",
            "Sec-Fetch-Dest": "document",
            "Sec-Fetch-Mode": "navigate",
            "Sec-Fetch-Site": "none",
            "Sec-Fetch-User": "?1",
            "Cache-Control": "max-age=0",
        })
        self.robot_parser: urllib.robotparser.RobotFileParser | None = None
        self.last_request_time = 0.0
        self._last_url: str | None = None  # Track last URL for Referer header
        self._init_robots_txt()

    def _init_robots_txt(self) -> None:
        """Initialize robots.txt parser if enabled."""
        if not self.config.check_robots:
            return

        try:
            robots_url = urllib.parse.urljoin(self.config.base_url, "/robots.txt")
            self.robot_parser = urllib.robotparser.RobotFileParser()
            self.robot_parser.set_url(robots_url)
            self.robot_parser.read()
            logger.info("Loaded robots.txt successfully")
        except Exception as e:
            logger.warning(f"Failed to load robots.txt: {e}. Continuing without robots.txt check.")
            self.robot_parser = None

    def _check_robots_allowed(self, url: str) -> bool:
        """Check if URL is allowed by robots.txt.

        Args:
            url: URL to check

        Returns:
            True if allowed, False otherwise
        """
        if not self.config.check_robots or self.robot_parser is None:
            return True

        try:
            return self.robot_parser.can_fetch(self.config.user_agent, url)
        except Exception as e:
            logger.warning(f"Error checking robots.txt for {url}: {e}")
            return True  # Default to allowing if check fails

    def _rate_limit(self) -> None:
        """Apply rate limiting delay between requests with randomization.
        
        Adds randomized delay (50-150% of base delay) to mimic human behavior
        and avoid detection patterns.
        """
        elapsed = time.time() - self.last_request_time
        
        # Randomize delay: 50-150% of base delay to avoid detection patterns
        randomized_delay = self.config.delay_seconds * random.uniform(0.5, 1.5)
        
        if elapsed < randomized_delay:
            sleep_time = randomized_delay - elapsed
            time.sleep(sleep_time)
        self.last_request_time = time.time()

    @retry_with_backoff(max_attempts=3, initial_delay=1.0, backoff_factor=2.0)
    def _fetch_page(self, url: str) -> str:
        """Fetch HTML content from URL with retry logic.

        Args:
            url: URL to fetch

        Returns:
            HTML content as string

        Raises:
            NetworkError: If request fails after retries
        """
        if not self._check_robots_allowed(url):
            raise NetworkError(f"URL not allowed by robots.txt: {url}")

        self._rate_limit()

        # Add Referer header to make requests look more legitimate
        # Set referer to AntWiki homepage for first request, then previous page
        headers = {}
        if self._last_url:
            headers['Referer'] = self._last_url
        else:
            headers['Referer'] = self.config.base_url

        try:
            response = self.session.get(url, timeout=self.config.timeout_seconds, headers=headers)
            response.raise_for_status()
            # Store this URL as referer for next request
            self._last_url = url
            return response.text
        except requests.RequestException as e:
            raise NetworkError(f"Failed to fetch {url}: {e}") from e

    def get_species_list(self) -> list[str]:
        """Discover all species pages from AntWiki.

        Attempts to find species pages by scraping category pages and index pages.

        Returns:
            List of species page names (URL-safe format, e.g., "Camponotus_pennsylvanicus")
        """
        logger.info("Discovering species list from AntWiki...")
        species: set[str] = set()

        # Try to find species from category pages
        category_urls = [
            self.config.base_url.rstrip("/") + "/Category:Ants_by_Genus",
            self.config.base_url.rstrip("/") + "/Category:Ants_by_Species",
        ]

        for category_url in category_urls:
            try:
                html = self._fetch_page(category_url)
                soup = BeautifulSoup(html, "html.parser")

                # Find links to species pages (typically in category listing)
                for link in soup.find_all("a", href=True):
                    href = link["href"]
                    if "/wiki/" in href:
                        # Extract page name from URL
                        page_name = href.split("/wiki/")[-1]
                        # Filter out non-species pages (categories, templates, etc.)
                        if (
                            page_name
                            and not page_name.startswith("Category:")
                            and not page_name.startswith("Template:")
                            and not page_name.startswith("File:")
                            and not page_name.startswith("User:")
                            and not page_name.startswith("Help:")
                            and ":" not in page_name
                        ):
                            species.add(page_name)

            except Exception as e:
                logger.warning(f"Failed to scrape category {category_url}: {e}")

        # Also try the main species index if available
        index_url = self.config.base_url.rstrip("/") + "/Category:Ants"
        try:
            html = self._fetch_page(index_url)
            soup = BeautifulSoup(html, "html.parser")
            for link in soup.find_all("a", href=True):
                href = link["href"]
                if "/wiki/" in href:
                    page_name = href.split("/wiki/")[-1]
                    if page_name and ":" not in page_name and not page_name.startswith("Category"):
                        species.add(page_name)
        except Exception as e:
            logger.warning(f"Failed to scrape index {index_url}: {e}")

        species_list = sorted(species)
        logger.info(f"Discovered {len(species_list)} species pages")
        return species_list

    def extract_measurements(self, html: str) -> dict[str, Any]:
        """Extract morphological measurements from HTML.

        Parses tables and infoboxes to find measurement data.

        Args:
            html: HTML content of species page

        Returns:
            Dictionary of measurements (e.g., {"worker_length_mm": [6.0, 13.0]})
        """
        soup = BeautifulSoup(html, "html.parser")
        measurements: dict[str, Any] = {}

        # Look for measurement tables
        for table in soup.find_all("table", class_=lambda x: x and ("infobox" in x.lower() or "measurement" in x.lower())):
            rows = table.find_all("tr")
            for row in rows:
                cells = row.find_all(["td", "th"])
                if len(cells) >= 2:
                    key_cell = cells[0].get_text(strip=True).lower()
                    value_cell = cells[1].get_text(strip=True)

                    # Look for common measurement patterns
                    if any(term in key_cell for term in ["length", "width", "height", "size", "measurement"]):
                        # Try to extract numeric values
                        import re

                        numbers = re.findall(r"\d+\.?\d*", value_cell)
                        if numbers:
                            key = key_cell.replace(" ", "_").replace("(", "").replace(")", "")
                            measurements[key] = [float(n) for n in numbers]

        # Also check infobox divs
        for infobox in soup.find_all("div", class_=lambda x: x and "infobox" in str(x).lower()):
            for item in infobox.find_all(["div", "span"], class_=lambda x: x and ("measurement" in str(x).lower() or "size" in str(x).lower())):
                text = item.get_text(strip=True)
                if "mm" in text or "cm" in text:
                    import re

                    match = re.search(r"(\d+\.?\d*)\s*-\s*(\d+\.?\d*)\s*(mm|cm)", text)
                    if match:
                        key = "measurement_" + match.group(3)
                        measurements[key] = [float(match.group(1)), float(match.group(2))]

        return measurements

    def extract_traits(self, html: str) -> list[str]:
        """Extract behavioral and ecological traits from HTML.

        Searches multiple sections for trait information.

        Args:
            html: HTML content of species page

        Returns:
            List of trait strings
        """
        soup = BeautifulSoup(html, "html.parser")
        traits: set[str] = set()

        # Look for trait lists in various sections
        section_keywords = ["behavior", "ecology", "biology", "habitat", "diet", "nesting"]

        for keyword in section_keywords:
            # Find sections with these keywords
            for heading in soup.find_all(["h2", "h3", "h4"]):
                heading_text = heading.get_text(strip=True).lower()
                if keyword in heading_text:
                    # Extract content from this section
                    section = heading.find_next_sibling(["div", "p", "ul"])
                    if section:
                        # Look for list items or text with trait keywords
                        for item in section.find_all(["li", "p"]):
                            text = item.get_text(strip=True).lower()
                            # Common trait patterns
                            trait_keywords = [
                                "arboreal",
                                "terrestrial",
                                "carnivorous",
                                "herbivorous",
                                "omnivorous",
                                "polygynous",
                                "monogynous",
                                "nocturnal",
                                "diurnal",
                                "social",
                                "solitary",
                            ]
                            for trait_kw in trait_keywords:
                                if trait_kw in text:
                                    traits.add(trait_kw)

        # Also check infobox for traits
        for infobox in soup.find_all("table", class_=lambda x: x and "infobox" in str(x).lower()):
            rows = infobox.find_all("tr")
            for row in rows:
                text = row.get_text(strip=True).lower()
                trait_keywords = [
                    "arboreal",
                    "terrestrial",
                    "carnivorous",
                    "herbivorous",
                    "omnivorous",
                    "polygynous",
                    "monogynous",
                ]
                for trait_kw in trait_keywords:
                    if trait_kw in text:
                        traits.add(trait_kw)

        return sorted(traits)

    def extract_description(self, html: str) -> str:
        """Extract species description from main content.

        Args:
            html: HTML content of species page

        Returns:
            Description text
        """
        soup = BeautifulSoup(html, "html.parser")

        # Find main content area (usually after infobox)
        content = soup.find("div", {"id": "mw-content-text"})
        if not content:
            content = soup.find("div", class_=lambda x: x and "content" in str(x).lower())

        if content:
            # Get first few paragraphs as description
            paragraphs = content.find_all("p", limit=5)
            description_parts = [p.get_text(strip=True) for p in paragraphs if p.get_text(strip=True)]
            return " ".join(description_parts)

        return ""

    def extract_taxonomy(self, html: str) -> dict[str, str]:
        """Extract taxonomic information from HTML.

        Args:
            html: HTML content of species page

        Returns:
            Dictionary with taxonomic ranks (e.g., {"genus": "Camponotus", "species": "pennsylvanicus"})
        """
        soup = BeautifulSoup(html, "html.parser")
        taxonomy: dict[str, str] = {}

        # Look for taxonomy in infobox
        for infobox in soup.find_all("table", class_=lambda x: x and "infobox" in str(x).lower()):
            rows = infobox.find_all("tr")
            for row in rows:
                cells = row.find_all(["td", "th"])
                if len(cells) >= 2:
                    key = cells[0].get_text(strip=True).lower()
                    value = cells[1].get_text(strip=True)

                    taxonomy_ranks = ["kingdom", "phylum", "class", "order", "family", "subfamily", "genus", "species"]
                    for rank in taxonomy_ranks:
                        if rank in key:
                            taxonomy[rank] = value
                            break

        # Also try to extract from page title
        title_tag = soup.find("h1", {"id": "firstHeading"})
        if title_tag:
            title = title_tag.get_text(strip=True)
            # Try to parse binomial name
            parts = title.split()
            if len(parts) >= 2:
                taxonomy["genus"] = parts[0]
                taxonomy["species"] = " ".join(parts[1:])

        return taxonomy

    def extract_distribution(self, html: str) -> dict[str, Any]:
        """Extract geographic distribution data.

        Args:
            html: HTML content of species page

        Returns:
            Dictionary with distribution information
        """
        soup = BeautifulSoup(html, "html.parser")
        distribution: dict[str, Any] = {}

        # Look for distribution section
        for heading in soup.find_all(["h2", "h3", "h4"]):
            heading_text = heading.get_text(strip=True).lower()
            if "distribution" in heading_text or "range" in heading_text:
                section = heading.find_next_sibling(["div", "p", "ul"])
                if section:
                    text = section.get_text(strip=True)
                    distribution["description"] = text

                    # Try to extract countries/regions
                    import re

                    # Common geographic terms
                    geo_patterns = [
                        r"([A-Z][a-z]+(?:\s+[A-Z][a-z]+)*)\s+(?:state|province|region|country)",
                        r"in\s+([A-Z][a-z]+(?:\s+[A-Z][a-z]+)*)",
                    ]
                    locations = []
                    for pattern in geo_patterns:
                        matches = re.findall(pattern, text)
                        locations.extend(matches)
                    if locations:
                        distribution["locations"] = list(set(locations))

        return distribution

    def extract_images(self, html: str) -> list[str]:
        """Extract image URLs from page.

        Args:
            html: HTML content of species page

        Returns:
            List of image URLs
        """
        soup = BeautifulSoup(html, "html.parser")
        images: list[str] = []

        for img in soup.find_all("img"):
            src = img.get("src", "")
            if src:
                # Convert relative URLs to absolute
                if src.startswith("//"):
                    src = "https:" + src
                elif src.startswith("/"):
                    src = urllib.parse.urljoin(self.config.base_url, src)
                images.append(src)

        return images

    def scrape_species_page(self, species_name: str) -> dict[str, Any]:
        """Scrape a single species page and extract all data.

        Args:
            species_name: Species page name (URL-safe format, e.g., "Camponotus_pennsylvanicus")

        Returns:
            Dictionary with all extracted data:
            - species: Species name
            - measurements: Morphological measurements
            - traits: Behavioral/ecological traits
            - description: Species description
            - taxonomy: Taxonomic information
            - distribution: Geographic distribution
            - images: Image URLs
            - url: Source URL
        """
        logger.debug(f"Scraping species page: {species_name}")

        # Construct URL
        page_url = urllib.parse.urljoin(self.config.base_url, species_name)

        try:
            html = self._fetch_page(page_url)

            # Extract all sections
            data: dict[str, Any] = {
                "species": species_name.replace("_", " "),
                "measurements": self.extract_measurements(html),
                "traits": self.extract_traits(html),
                "description": self.extract_description(html),
                "taxonomy": self.extract_taxonomy(html),
                "distribution": self.extract_distribution(html),
                "images": self.extract_images(html),
                "url": page_url,
            }

            logger.debug(f"Successfully scraped {species_name}: {len(data.get('traits', []))} traits, {len(data.get('measurements', {}))} measurements")
            return data

        except Exception as e:
            logger.error(f"Failed to scrape {species_name}: {e}", exc_info=True)
            raise

    def scrape_all_species(
        self,
        output_dir: Path | None = None,
        limit: int | None = None,
        resume: bool = False,
    ) -> dict[str, Any]:
        """Scrape all species pages with progress tracking.

        Args:
            output_dir: Output directory (default: config output_dir)
            limit: Limit number of species to scrape (for testing)
            resume: Resume from checkpoint if available

        Returns:
            Dictionary with scraping statistics
        """
        if output_dir is None:
            output_dir = self.config.output_dir

        ensure_directory(output_dir)
        species_dir = output_dir / "species"
        ensure_directory(species_dir)

        log_file = output_dir / "scraping_log.jsonl"
        checkpoint_file = output_dir / "checkpoint.json"

        # Load checkpoint if resuming
        completed_species: set[str] = set()
        if resume and checkpoint_file.exists():
            try:
                checkpoint_data = load_json(checkpoint_file)
                # Check log file for completed species
                if log_file.exists():
                    for entry in read_jsonl(log_file):
                        if entry.get("status") == "completed":
                            completed_species.add(entry.get("species", ""))
                logger.info(f"Resuming: {len(completed_species)} species already completed")
            except Exception as e:
                logger.warning(f"Failed to load checkpoint: {e}")

        # Get species list
        species_list = self.get_species_list()
        if limit:
            species_list = species_list[:limit]

        # Filter out already completed
        species_to_scrape = [s for s in species_list if s not in completed_species]

        logger.info(f"Scraping {len(species_to_scrape)} species pages...")

        stats = {"total": len(species_to_scrape), "completed": 0, "failed": 0, "failed_species": []}
        all_data: list[dict[str, Any]] = []

        for i, species_name in enumerate(species_to_scrape, 1):
            logger.info(f"Scraping {i}/{len(species_to_scrape)}: {species_name}")

            try:
                data = self.scrape_species_page(species_name)

                # Save individual species file
                safe_name = sanitize_filename(species_name)
                species_file = species_dir / f"{safe_name}.json"
                dump_json(data, species_file, indent=2)

                # Add to combined data
                all_data.append(data)

                # Log success
                log_entry = {
                    "timestamp": time.time(),
                    "species": species_name,
                    "status": "completed",
                }
                # Append to log file
                with open(log_file, "a", encoding="utf-8") as f:
                    f.write(json.dumps(log_entry) + "\n")

                stats["completed"] += 1

                # Save checkpoint periodically
                if i % 10 == 0:
                    checkpoint_data = {
                        "last_species": species_name,
                        "completed_count": stats["completed"],
                        "total_count": len(species_to_scrape),
                    }
                    dump_json(checkpoint_data, checkpoint_file, indent=2)

            except Exception as e:
                logger.error(f"Failed to scrape {species_name}: {e}")
                stats["failed"] += 1
                stats["failed_species"].append(species_name)

                # Log failure
                log_entry = {
                    "timestamp": time.time(),
                    "species": species_name,
                    "status": "failed",
                    "error": str(e),
                }
                # Append to log file
                with open(log_file, "a", encoding="utf-8") as f:
                    f.write(json.dumps(log_entry) + "\n")

        # Save combined dataset
        all_species_file = output_dir / "all_species.json"
        dump_json(all_data, all_species_file, indent=2)

        # Final checkpoint
        checkpoint_data = {
            "last_species": species_list[-1] if species_list else None,
            "completed_count": stats["completed"],
            "total_count": len(species_to_scrape),
            "failed_count": stats["failed"],
        }
        dump_json(checkpoint_data, checkpoint_file, indent=2)

        logger.info(f"Scraping complete: {stats['completed']} succeeded, {stats['failed']} failed")
        return stats

