"""AntWiki web scraping utilities.

This module provides functionality to scrape phenotype data from AntWiki.org,
including species pages, morphological traits, behavioral characteristics,
and ecological information.
"""

from __future__ import annotations

import re
import time
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Set, Tuple
from urllib.parse import urljoin, urlparse

import requests
from bs4 import BeautifulSoup

from metainformant.core import errors, io, logging, validation

logger = logging.get_logger(__name__)


class AntWikiScraperConfig:
    """Configuration for AntWiki scraper."""

    def __init__(
        self,
        base_url: str = "https://www.antwiki.org/wiki/",
        delay_between_requests: float | None = None,
        timeout: int | None = None,
        max_retries: int = 3,
        user_agent: str = "MetaInformant/1.0 (research@metainformant.org)",
        delay_seconds: float | None = None,
        timeout_seconds: int | None = None,
        check_robots: bool = True,
        output_dir: Any = None,
    ):
        """Initialize scraper configuration.

        Args:
            base_url: Base URL for AntWiki
            delay_between_requests: Seconds to wait between requests
            timeout: Request timeout in seconds
            max_retries: Maximum number of retries for failed requests
            user_agent: User agent string for requests
            delay_seconds: Alias for delay_between_requests
            timeout_seconds: Alias for timeout
            check_robots: Whether to check robots.txt
            output_dir: Default output directory
        """
        self.base_url = base_url

        # Handle parameter aliases
        if delay_seconds is not None:
            self._delay = float(delay_seconds)
        elif delay_between_requests is not None:
            self._delay = float(delay_between_requests)
        else:
            self._delay = 2.0

        if timeout_seconds is not None:
            self._timeout = int(timeout_seconds)
        elif timeout is not None:
            self._timeout = int(timeout)
        else:
            self._timeout = 30

        self.delay_between_requests = self._delay
        self.timeout = self._timeout
        self.max_retries = max_retries
        self.user_agent = user_agent
        self.check_robots = check_robots
        self.output_dir = output_dir

        # Validate configuration
        validation.validate_type(self._delay, (int, float), "delay_between_requests")
        validation.validate_range(self._delay, min_val=0.1, name="delay_between_requests")
        validation.validate_range(self._timeout, min_val=1, name="timeout")
        validation.validate_range(max_retries, min_val=0, name="max_retries")

    @property
    def delay_seconds(self) -> float:
        """Alias for delay_between_requests."""
        return self._delay

    @property
    def timeout_seconds(self) -> int:
        """Alias for timeout."""
        return self._timeout


class AntWikiScraper:
    """Web scraper for AntWiki phenotype data."""

    def __init__(self, config: Optional[AntWikiScraperConfig] = None):
        """Initialize the scraper.

        Args:
            config: Scraper configuration
        """
        self.config = config or AntWikiScraperConfig()
        self.session = requests.Session()
        self.session.headers.update(
            {
                "User-Agent": self.config.user_agent,
                "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
                "Accept-Language": "en-US,en;q=0.5",
                "Accept-Encoding": "gzip, deflate",
                "Connection": "keep-alive",
                "Upgrade-Insecure-Requests": "1",
            }
        )

        # Track scraped URLs to avoid duplicates
        self.scraped_urls: Set[str] = set()

        logger.info(f"Initialized AntWiki scraper for {self.config.base_url}")

    def _make_request(self, url: str) -> Optional[BeautifulSoup]:
        """Make an HTTP request with retry logic.

        Args:
            url: URL to request

        Returns:
            BeautifulSoup object if successful, None if failed
        """
        for attempt in range(self.config.max_retries + 1):
            try:
                logger.debug(f"Requesting {url} (attempt {attempt + 1})")
                response = self.session.get(url, timeout=self.config.timeout)

                if response.status_code == 200:
                    soup = BeautifulSoup(response.content, "html.parser")
                    return soup
                elif response.status_code == 404:
                    logger.warning(f"Page not found: {url}")
                    return None
                elif response.status_code == 403:
                    logger.warning(f"HTTP 403 Forbidden for {url}")
                    response.raise_for_status()
                elif response.status_code >= 500:
                    logger.warning(f"Server error {response.status_code} for {url}")
                    if attempt < self.config.max_retries:
                        time.sleep(2**attempt)  # Exponential backoff
                        continue
                else:
                    logger.warning(f"HTTP {response.status_code} for {url}")
                    response.raise_for_status()

            except requests.HTTPError:
                raise  # Let HTTP errors propagate immediately (403, etc.)
            except requests.RequestException as e:
                logger.warning(f"Request failed for {url}: {e}")
                if attempt < self.config.max_retries:
                    time.sleep(2**attempt)
                    continue

        return None

    def _respectful_delay(self) -> None:
        """Add delay between requests to be respectful to the server."""
        time.sleep(self.config.delay_between_requests)

    def scrape_species_page(self, species_name_or_url: str) -> Optional[Dict[str, Any]]:
        """Scrape a single species page from AntWiki.

        Args:
            species_name_or_url: Species name (e.g. "Camponotus_pennsylvanicus") or full URL

        Returns:
            Dictionary with scraped data, or None if scraping failed
        """
        # Build full URL if only species name given
        if species_name_or_url.startswith("http"):
            species_url = species_name_or_url
        elif self.config.base_url.endswith("/"):
            species_url = f"{self.config.base_url}{species_name_or_url}"
        else:
            species_url = f"{self.config.base_url}/{species_name_or_url}"

        if species_url in self.scraped_urls:
            logger.debug(f"Already scraped: {species_url}")
            return None

        soup = self._make_request(species_url)
        if not soup:
            raise errors.NetworkError(f"Failed to fetch species page: {species_url}")

        self.scraped_urls.add(species_url)
        self._respectful_delay()

        try:
            # Derive species name from input or title
            if not species_name_or_url.startswith("http"):
                display_name = species_name_or_url.replace("_", " ")
            else:
                display_name = None

            # Extract species information
            data: Dict[str, Any] = {
                "url": species_url,
                "data_source": "antwiki_scraper",
                "scraped_at": time.time(),
                "measurements": {},
                "traits": [],
                "description": "",
                "taxonomy": {},
                "distribution": {},
                "images": [],
            }

            # Extract species name from title
            title_elem = soup.find("title")
            if title_elem:
                title_text = title_elem.get_text().strip()
                if " - AntWiki" in title_text:
                    species_name = title_text.replace(" - AntWiki", "").strip()
                    data["species"] = species_name
                    data["species_name"] = species_name
                    parts = species_name.split()
                    if len(parts) >= 2:
                        data["genus"] = parts[0]

            if display_name and "species" not in data:
                data["species"] = display_name
                data["species_name"] = display_name

            # Extract taxonomic information
            taxonomy = self._extract_taxonomy(soup)
            data["taxonomy"] = taxonomy

            # Extract morphological data
            morphology = self._extract_morphology(soup)
            if morphology:
                data["morphology"] = morphology
                data["measurements"] = morphology

            # Extract behavioral data
            behavior = self._extract_behavior(soup)
            if behavior:
                data["behavior"] = behavior

            # Extract ecological data
            ecology = self._extract_ecology(soup)
            if ecology:
                data["ecology"] = ecology

            # Extract distribution data
            distribution = self._extract_distribution(soup)
            if distribution:
                data["distribution"] = distribution

            # Extract description
            desc_elem = soup.find("div", {"class": "mw-parser-output"})
            if desc_elem:
                paragraphs = desc_elem.find_all("p", limit=3)
                data["description"] = " ".join(p.get_text().strip() for p in paragraphs)

            # Extract images
            img_elems = soup.find_all("img", src=re.compile(r"\.(jpg|jpeg|png|gif)", re.IGNORECASE))
            data["images"] = [img.get("src", "") for img in img_elems[:10]]

            # Extract traits from behavioral and ecological data
            if behavior:
                data["traits"].extend(v for v in behavior.values() if isinstance(v, str))

            # Calculate confidence score based on data completeness
            data["confidence_score"] = self._calculate_confidence_score(data)

            logger.info(f"Successfully scraped: {species_url}")
            return data

        except Exception as e:
            logger.error(f"Error scraping {species_url}: {e}")
            return None

    def _extract_taxonomy(self, soup: BeautifulSoup) -> Dict[str, Any]:
        """Extract taxonomic information from the page."""
        taxonomy = {}

        # Look for taxonomic hierarchy
        tax_elements = soup.find_all(["div", "span", "p"], class_=re.compile(r"tax|classif"))

        for elem in tax_elements:
            text = elem.get_text().strip()

            # Extract subfamily
            subfamily_match = re.search(r"Subfamily:\s*([^\s\n]+)", text, re.IGNORECASE)
            if subfamily_match:
                taxonomy["subfamily"] = subfamily_match.group(1).strip()

            # Extract tribe
            tribe_match = re.search(r"Tribe:\s*([^\s\n]+)", text, re.IGNORECASE)
            if tribe_match:
                taxonomy["tribe"] = tribe_match.group(1).strip()

        return taxonomy

    def _extract_morphology(self, soup: BeautifulSoup) -> Dict[str, Any]:
        """Extract morphological information."""
        morphology = {}

        # Look for morphology sections
        morph_sections = soup.find_all(
            ["h2", "h3", "h4"], string=re.compile(r"morph|size|color|head|body", re.IGNORECASE)
        )

        for section in morph_sections:
            section_text = self._extract_section_content(section)

            # Extract size information
            size_matches = re.findall(r"(\d+(?:\.\d+)?)\s*(mm|cm)", section_text, re.IGNORECASE)
            if size_matches:
                for value, unit in size_matches:
                    if "body" in section_text.lower():
                        morphology["body_length_mm"] = float(value)
                    elif "head" in section_text.lower():
                        morphology["head_width_mm"] = float(value)

            # Extract color information
            color_info = {}
            if "color" in section_text.lower():
                # Simple color extraction - could be enhanced
                colors = ["black", "red", "yellow", "brown", "orange", "white"]
                for color in colors:
                    if color in section_text.lower():
                        color_info["general"] = color
                        break

            if color_info:
                morphology["color"] = color_info

        return morphology

    def _extract_behavior(self, soup: BeautifulSoup) -> Dict[str, Any]:
        """Extract behavioral information."""
        behavior = {}

        # Look for behavior sections
        behavior_sections = soup.find_all(
            ["h2", "h3", "h4"], string=re.compile(r"behav|forag|nest|colony", re.IGNORECASE)
        )

        for section in behavior_sections:
            section_text = self._extract_section_content(section)

            # Extract foraging strategy
            if "forag" in section_text.lower():
                foraging_patterns = ["generalist", "specialist", "predator", "herbivore", "omnivore"]
                for pattern in foraging_patterns:
                    if pattern in section_text.lower():
                        behavior["foraging_strategy"] = pattern
                        break

            # Extract nest type
            if "nest" in section_text.lower():
                nest_patterns = ["arboreal", "ground", "soil", "wood", "leaf"]
                for pattern in nest_patterns:
                    if pattern in section_text.lower():
                        behavior["nest_type"] = pattern
                        break

            # Extract colony size
            colony_matches = re.search(r"colony.*?(small|medium|large|monogyne|polygyne)", section_text, re.IGNORECASE)
            if colony_matches:
                behavior["colony_size_category"] = colony_matches.group(1).lower()

        return behavior

    def _extract_ecology(self, soup: BeautifulSoup) -> Dict[str, Any]:
        """Extract ecological information."""
        ecology = {}

        # Look for ecology sections
        ecology_sections = soup.find_all(["h2", "h3", "h4"], string=re.compile(r"ecol|habit|environ", re.IGNORECASE))

        for section in ecology_sections:
            section_text = self._extract_section_content(section)

            # Extract habitat information
            habitats = ["forest", "grassland", "desert", "urban", "tropical", "temperate"]
            for habitat in habitats:
                if habitat in section_text.lower():
                    ecology["habitat"] = habitat
                    break

        return ecology

    def _extract_distribution(self, soup: BeautifulSoup) -> Dict[str, Any]:
        """Extract distribution information."""
        distribution = {}

        # Look for distribution sections
        dist_sections = soup.find_all(["h2", "h3", "h4"], string=re.compile(r"distrib|range|location", re.IGNORECASE))

        for section in dist_sections:
            section_text = self._extract_section_content(section)

            # Extract geographic regions
            regions = ["africa", "asia", "europe", "north america", "south america", "australia", "antarctica"]
            found_regions = []

            for region in regions:
                if region in section_text.lower():
                    found_regions.append(region)

            if found_regions:
                distribution["regions"] = found_regions

        return distribution

    def _extract_section_content(self, header_element) -> str:
        """Extract content following a header element."""
        content_parts = []

        # Get all sibling elements until next header
        current = header_element.find_next_sibling()
        while current and current.name not in ["h1", "h2", "h3", "h4", "h5", "h6"]:
            if current.name in ["p", "div", "span"]:
                content_parts.append(current.get_text().strip())
            current = current.find_next_sibling()

        return " ".join(content_parts)

    def _calculate_confidence_score(self, data: Dict[str, Any]) -> float:
        """Calculate confidence score based on data completeness."""
        score_components = []

        # Taxonomic information (30%)
        tax_score = 0
        if data.get("genus"):
            tax_score += 0.5
        if data.get("subfamily"):
            tax_score += 0.3
        if data.get("tribe"):
            tax_score += 0.2
        score_components.append(tax_score * 0.3)

        # Morphological data (30%)
        morph_score = 0
        morphology = data.get("morphology", {})
        if morphology.get("body_length_mm"):
            morph_score += 0.4
        if morphology.get("head_width_mm"):
            morph_score += 0.4
        if morphology.get("color"):
            morph_score += 0.2
        score_components.append(morph_score * 0.3)

        # Behavioral data (20%)
        behavior_score = 0
        behavior = data.get("behavior", {})
        if behavior.get("foraging_strategy"):
            behavior_score += 0.5
        if behavior.get("nest_type"):
            behavior_score += 0.3
        if behavior.get("colony_size_category"):
            behavior_score += 0.2
        score_components.append(behavior_score * 0.2)

        # Ecological data (20%)
        ecology_score = 0
        ecology = data.get("ecology", {})
        if ecology.get("habitat"):
            ecology_score += 0.6
        distribution = data.get("distribution", {})
        if distribution.get("regions"):
            ecology_score += 0.4
        score_components.append(ecology_score * 0.2)

        return sum(score_components)

    def scrape_species_list(self, limit: Optional[int] = None) -> Iterator[str]:
        """Scrape list of species URLs from AntWiki.

        Args:
            limit: Maximum number of species URLs to return

        Yields:
            Species page URLs
        """
        # This is a simplified implementation - real AntWiki scraping would need
        # to handle their actual page structure and pagination

        logger.info("Scraping species list from AntWiki (simplified implementation)")

        # For demonstration, we'll use a hardcoded list of common species
        # In a real implementation, this would scrape directory pages
        common_species = [
            "Camponotus_pennsylvanicus",
            "Formica_rufa",
            "Lasius_niger",
            "Myrmica_rubra",
            "Pheidole_pallidula",
        ]

        for species in common_species:
            if limit and len(self.scraped_urls) >= limit:
                break

            if self.config.base_url.endswith("/"):
                url = f"{self.config.base_url}{species}"
            else:
                url = f"{self.config.base_url}/{species}"
            yield url

    def bulk_scrape_species(self, species_urls: List[str], output_dir: str | Path) -> List[Dict[str, Any]]:
        """Scrape multiple species pages and save results.

        Args:
            species_urls: List of species URLs to scrape
            output_dir: Directory to save individual JSON files

        Returns:
            List of scraped data dictionaries
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        results = []
        successful = 0

        for url in species_urls:
            try:
                data = self.scrape_species_page(url)
                if data:
                    # Save individual file
                    species_name = data.get("species_name", "unknown").replace(" ", "_")
                    output_file = output_dir / f"{species_name}.json"
                    io.dump_json(data, output_file, indent=2)

                    results.append(data)
                    successful += 1

                # Progress logging
                if len(results) % 10 == 0:
                    logger.info(f"Scraped {len(results)} species pages")

            except Exception as e:
                logger.error(f"Failed to scrape {url}: {e}")
                continue

        logger.info(f"Successfully scraped {successful}/{len(species_urls)} species pages")
        return results

    def _fetch_page(self, url: str) -> str:
        """Fetch a page and return raw HTML string.

        Args:
            url: URL to fetch

        Returns:
            HTML string

        Raises:
            NetworkError: If request fails
        """
        for attempt in range(self.config.max_retries + 1):
            try:
                logger.debug(f"Fetching {url} (attempt {attempt + 1})")
                response = self.session.get(url, timeout=self.config.timeout)
                response.raise_for_status()
                self._respectful_delay()
                return response.text
            except requests.HTTPError as e:
                raise errors.NetworkError(f"HTTP error fetching {url}: {e}") from e
            except requests.RequestException as e:
                if attempt < self.config.max_retries:
                    time.sleep(2**attempt)
                    continue
                raise errors.NetworkError(f"Failed to fetch {url}: {e}") from e
        raise errors.NetworkError(f"Failed to fetch {url} after {self.config.max_retries + 1} attempts")

    def extract_measurements(self, html: str) -> Dict[str, Any]:
        """Extract measurements from HTML string.

        Args:
            html: Raw HTML string

        Returns:
            Dictionary of measurements
        """
        soup = BeautifulSoup(html, "html.parser")
        return self._extract_morphology(soup)

    def extract_traits(self, html: str) -> List[str]:
        """Extract traits from HTML string.

        Args:
            html: Raw HTML string

        Returns:
            List of trait strings
        """
        soup = BeautifulSoup(html, "html.parser")
        behavior = self._extract_behavior(soup)
        traits = [v for v in behavior.values() if isinstance(v, str)]
        return traits

    def extract_taxonomy(self, html: str) -> Dict[str, Any]:
        """Extract taxonomy from HTML string.

        Args:
            html: Raw HTML string

        Returns:
            Dictionary of taxonomic information
        """
        soup = BeautifulSoup(html, "html.parser")
        return self._extract_taxonomy(soup)

    def get_species_list(self, limit: Optional[int] = None) -> List[str]:
        """Get list of species names from AntWiki.

        Args:
            limit: Maximum number of species to return

        Returns:
            List of species name strings
        """
        return list(self.scrape_species_list(limit=limit))

    def scrape_all_species(self, output_dir: str | Path | None = None, limit: Optional[int] = None) -> Dict[str, Any]:
        """Scrape all species and save to output directory.

        Args:
            output_dir: Directory to save results
            limit: Maximum number of species to scrape

        Returns:
            Statistics dictionary with total, completed, failed counts
        """
        if output_dir is None:
            output_dir = self.config.output_dir
        if output_dir is None:
            raise ValueError("No output directory specified")

        output_dir = Path(output_dir)
        species_dir = output_dir / "species"
        species_dir.mkdir(parents=True, exist_ok=True)

        species_urls = list(self.scrape_species_list(limit=limit))
        results = []
        failed = 0

        for url in species_urls:
            try:
                data = self.scrape_species_page(url)
                if data:
                    species_name = data.get("species_name", data.get("species", "unknown")).replace(" ", "_")
                    output_file = species_dir / f"{species_name}.json"
                    io.dump_json(data, output_file, indent=2)
                    results.append(data)
                else:
                    failed += 1
            except Exception as e:
                logger.error(f"Failed to scrape {url}: {e}")
                failed += 1

        # Save combined file
        all_species_file = output_dir / "all_species.json"
        io.dump_json(results, all_species_file, indent=2)

        stats = {
            "total": len(species_urls),
            "completed": len(results),
            "failed": failed,
        }

        logger.info(f"Scraped {stats['completed']}/{stats['total']} species")
        return stats

    def close(self) -> None:
        """Close the scraper session."""
        self.session.close()
        logger.info("AntWiki scraper session closed")


def load_scraper_config(config_path: str | Path | None = None) -> AntWikiScraperConfig:
    """Load scraper configuration from file, or return defaults.

    Args:
        config_path: Path to configuration file (optional; returns defaults if None)

    Returns:
        AntWikiScraperConfig object
    """
    if config_path is None:
        return AntWikiScraperConfig()

    config_path = validation.validate_path_exists(Path(config_path))

    try:
        config_data = io.load_json(config_path)
        return AntWikiScraperConfig(**config_data)
    except Exception as e:
        logger.error(f"Failed to load scraper config from {config_path}: {e}")
        raise errors.IOError(f"Failed to load scraper config: {e}") from e


def create_default_scraper_config(output_path: str | Path) -> None:
    """Create a default scraper configuration file.

    Args:
        output_path: Path to save the configuration
    """
    config = AntWikiScraperConfig()
    config_dict = {
        "base_url": config.base_url,
        "delay_between_requests": config.delay_between_requests,
        "timeout": config.timeout,
        "max_retries": config.max_retries,
        "user_agent": config.user_agent,
    }

    io.dump_json(config_dict, Path(output_path), indent=2)
    logger.info(f"Created default scraper config at {output_path}")
