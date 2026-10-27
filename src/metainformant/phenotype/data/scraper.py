"""AntWiki web scraping utilities.

This module provides functionality to scrape phenotype data from AntWiki.org,
including species pages, morphological traits, behavioral characteristics,
and ecological information.

HTTP client: every request goes through one ``requests.Session`` per scraper. The optional ``cloudscraper``
package that appears in the project's dependency extras is **not** used here, and nothing in this module solves
Cloudflare or other bot-protection challenges: a challenged request (for example HTTP 403) surfaces as an error.

robots.txt: when ``AntWikiScraperConfig.check_robots`` is true (the default) the scraper consults the site's
robots.txt through :mod:`urllib.robotparser`. It is downloaded once per origin per scraper (session) - or supplied
up front with :meth:`AntWikiScraper.load_robots_txt` - and then:

* a URL it disallows for the configured user agent raises :class:`RobotsDisallowedError`, before any request for
  that URL is sent; redirects are followed one hop at a time so a redirect cannot reach a disallowed URL either;
* the wait between requests is the larger of the configured delay and the site's ``Crawl-delay``;
* a robots.txt that cannot be retrieved (5xx, 408, 429, network error) makes the scraper refuse to proceed, a 401/403
  disallows everything, and any other 4xx is read as "no robots.txt".

``check_robots=False`` is an explicit opt-out: robots.txt is never requested or applied.

Known limits of :mod:`urllib.robotparser`, which apply here: ``*`` and ``$`` wildcards inside a path are not
understood, ``Allow``/``Disallow`` rules are matched first-to-last rather than by longest match (RFC 9309), and
robots.txt groups are matched against the user agent's leading product token (the text before the first ``/``).
Fractional ``Crawl-delay`` values, which the standard library ignores, are rounded up to whole seconds.
"""

from __future__ import annotations

import math
import re
import time
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Optional, Set, Union
from urllib import robotparser
from urllib.parse import urljoin, urlsplit

import requests
from bs4 import BeautifulSoup

from metainformant.core import io
from metainformant.core.data import validation
from metainformant.core.io.paths import sanitize_filename
from metainformant.core.utils import errors, logging

logger = logging.get_logger(__name__)

_MAX_REDIRECTS = 10
_REDIRECT_STATUSES = frozenset({301, 302, 303, 307, 308})
# urllib.robotparser only reads integer Crawl-delay values and silently drops e.g. "Crawl-delay: 2.5".
_FRACTIONAL_CRAWL_DELAY = re.compile(r"^(\s*crawl-delay\s*:\s*)(\d*\.\d+)(\s*(?:#.*)?)$", re.IGNORECASE)


class RobotsDisallowedError(errors.NetworkError):
    """Raised when robots.txt forbids fetching a URL, or cannot be retrieved so access is refused."""


def _origin(url: str) -> str:
    """Return the ``scheme://host[:port]`` origin that robots.txt rules are scoped to."""
    parts = urlsplit(url)
    scheme = parts.scheme.lower()
    host = parts.hostname
    if not scheme or not host:
        raise errors.ValidationError(f"Cannot apply robots.txt to a non-absolute URL: {url!r}")
    host = host.lower()
    if ":" in host:  # IPv6 literal
        host = f"[{host}]"
    port = parts.port
    if port is not None and port != {"http": 80, "https": 443}.get(scheme):
        host = f"{host}:{port}"
    return f"{scheme}://{host}"


def _normalize_robots_lines(lines: Iterable[str]) -> List[str]:
    """Round fractional ``Crawl-delay`` values up so ``urllib.robotparser`` does not drop them."""
    normalized = []
    for line in lines:
        match = _FRACTIONAL_CRAWL_DELAY.match(line)
        if match:
            line = f"{match.group(1)}{math.ceil(float(match.group(2)))}{match.group(3)}"
        normalized.append(line)
    return normalized


def build_robots_parser(
    robots: Union[str, bytes, Iterable[str], robotparser.RobotFileParser],
    robots_url: str = "",
) -> robotparser.RobotFileParser:
    """Build a :class:`urllib.robotparser.RobotFileParser` from robots.txt content.

    Args:
        robots: The file's text, an iterable of its lines, or an already-parsed ``RobotFileParser`` (returned as is).
        robots_url: URL the rules came from (informational).

    Returns:
        A parser ready for ``can_fetch`` / ``crawl_delay`` queries.
    """
    if isinstance(robots, robotparser.RobotFileParser):
        return robots
    if isinstance(robots, (bytes, bytearray)):
        robots = bytes(robots).decode("utf-8", errors="replace")
    lines = robots.splitlines() if isinstance(robots, str) else list(robots)
    parser = robotparser.RobotFileParser()
    if robots_url:
        parser.set_url(robots_url)
    parser.parse(_normalize_robots_lines(lines))
    return parser


def robots_policy_from_response(
    status_code: int, body: str, robots_url: str = ""
) -> Optional[robotparser.RobotFileParser]:
    """Translate the HTTP response for ``/robots.txt`` into a parser, following ``urllib.robotparser.read()``.

    2xx: parse ``body``. 401/403: the site is refusing us, so disallow everything. Other 4xx: there is no robots.txt,
    so allow everything. Timeouts and rate limiting (408/429), 5xx and any unexpected status mean the rules are
    temporarily unknown: ``None``, which callers retry and then refuse.
    """
    if 200 <= status_code < 300:
        return build_robots_parser(body, robots_url)
    if status_code in (401, 403):
        return build_robots_parser(["User-agent: *", "Disallow: /"], robots_url)
    if 400 <= status_code < 500 and status_code not in (408, 429):
        return build_robots_parser([], robots_url)
    return None


def _redirect_target(current_url: str, status_code: int, location: Optional[str]) -> Optional[str]:
    """Absolute URL a response redirects to, or ``None`` when it is not a followable redirect."""
    if status_code not in _REDIRECT_STATUSES or not location:
        return None
    try:  # http.client decodes headers as latin-1, but servers send UTF-8 (same repair requests applies)
        location = location.encode("latin1").decode("utf-8")
    except (UnicodeEncodeError, UnicodeDecodeError):
        pass
    return urljoin(current_url, location)


def _species_output_path(species_dir: Path, data: Dict[str, Any]) -> Path:
    """Return the JSON file for a scraped species inside ``species_dir``.

    The name comes from the page title, which is remote input, so it is sanitised: a title such as ``../../x`` can
    never place the file outside ``species_dir``.
    """
    raw_name = str(data.get("species_name") or data.get("species") or "unknown")
    return species_dir / f"{sanitize_filename(raw_name.replace(' ', '_'))}.json"


class AntWikiScraperConfig:
    """Configuration for AntWiki scraper.

    ``delay_between_requests`` / ``delay_seconds`` and ``timeout`` / ``timeout_seconds`` are two spellings of the
    same two settings: reading or assigning either spelling reads or assigns one shared value, and every assignment is
    validated.
    """

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
            delay_between_requests: Minimum seconds to wait between requests (raised to the site's robots.txt
                ``Crawl-delay`` when that is larger and ``check_robots`` is true)
            timeout: Request timeout in seconds
            max_retries: Maximum number of retries for failed requests
            user_agent: User agent string for requests. robots.txt groups are matched against its leading product
                token (the text before the first ``/``), so keep the token that identifies this scraper first.
            delay_seconds: Alias for delay_between_requests (wins if both are given)
            timeout_seconds: Alias for timeout (wins if both are given)
            check_robots: Consult robots.txt (default True): it is downloaded once per origin per scraper, URLs it
                disallows raise ``RobotsDisallowedError`` and its ``Crawl-delay`` raises the delay between requests.
                ``False`` is an explicit opt-out: robots.txt is neither requested nor applied.
            output_dir: Default output directory
        """
        self.base_url = base_url

        # Handle parameter aliases (the *_seconds spelling wins when both are given)
        if delay_seconds is not None:
            delay = float(delay_seconds)
        elif delay_between_requests is not None:
            delay = float(delay_between_requests)
        else:
            delay = 2.0

        if timeout_seconds is not None:
            timeout_value = int(timeout_seconds)
        elif timeout is not None:
            timeout_value = int(timeout)
        else:
            timeout_value = 30

        self.delay_between_requests = delay  # validated by the property setter
        self.timeout = timeout_value  # validated by the property setter
        self.max_retries = max_retries
        self.user_agent = user_agent
        self.check_robots = check_robots
        self.output_dir = output_dir

        validation.validate_range(max_retries, min_val=0, name="max_retries")

    @property
    def delay_between_requests(self) -> float:
        """Minimum seconds between requests (``delay_seconds`` is the same setting)."""
        return self._delay

    @delay_between_requests.setter
    def delay_between_requests(self, value: float) -> None:
        validation.validate_type(value, (int, float), "delay_between_requests")
        validation.validate_range(value, min_val=0.1, name="delay_between_requests")
        self._delay = float(value)

    @property
    def delay_seconds(self) -> float:
        """Alias for delay_between_requests (readable and assignable)."""
        return self.delay_between_requests

    @delay_seconds.setter
    def delay_seconds(self, value: float) -> None:
        self.delay_between_requests = value

    @property
    def timeout(self) -> int:
        """Request timeout in seconds (``timeout_seconds`` is the same setting)."""
        return self._timeout

    @timeout.setter
    def timeout(self, value: int) -> None:
        validation.validate_type(value, (int, float), "timeout")
        validation.validate_range(value, min_val=1, name="timeout")
        self._timeout = int(value)

    @property
    def timeout_seconds(self) -> int:
        """Alias for timeout (readable and assignable)."""
        return self.timeout

    @timeout_seconds.setter
    def timeout_seconds(self, value: int) -> None:
        self.timeout = value


class AntWikiScraper:
    """Web scraper for AntWiki phenotype data.

    Pages are fetched with ``requests`` (not ``cloudscraper``; see the module docstring). With
    ``config.check_robots`` true, robots.txt is downloaded once per origin per scraper, disallowed URLs raise
    :class:`RobotsDisallowedError`, and the wait between requests is at least the site's ``Crawl-delay``.
    """

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

        # Parsed robots.txt per origin (scheme://host[:port]); filled once per origin per session
        self._robots: Dict[str, robotparser.RobotFileParser] = {}

        logger.info(f"Initialized AntWiki scraper for {self.config.base_url}")

    # ------------------------------------------------------------------
    # robots.txt compliance
    # ------------------------------------------------------------------

    def load_robots_txt(
        self,
        robots: Union[str, bytes, Iterable[str], robotparser.RobotFileParser],
        site_url: Optional[str] = None,
    ) -> robotparser.RobotFileParser:
        """Install robots.txt rules for the origin of ``site_url`` (default: ``config.base_url``).

        This is the same code path the scraper uses after downloading robots.txt itself, so rules can be supplied up
        front (from a saved copy, or in offline tests) and are then never downloaded for that origin.

        Args:
            robots: The file's text, an iterable of its lines, or an already-parsed
                :class:`urllib.robotparser.RobotFileParser`.
            site_url: Any URL on the site the rules belong to; defaults to ``config.base_url``.

        Returns:
            The parser now in force for that origin. It is ignored while ``check_robots`` is false.
        """
        origin = _origin(site_url or self.config.base_url)
        parser = build_robots_parser(robots, f"{origin}/robots.txt")
        self._robots[origin] = parser
        return parser

    def _robots_for(self, url: str) -> robotparser.RobotFileParser:
        """Return the parsed robots.txt for ``url``'s origin, downloading it on first use."""
        origin = _origin(url)
        parser = self._robots.get(origin)
        if parser is None:
            parser = self._download_robots(origin)
            self._robots[origin] = parser
        return parser

    def _download_robots(self, origin: str) -> robotparser.RobotFileParser:
        """Download and parse ``origin``/robots.txt, retrying transient failures.

        Raises:
            RobotsDisallowedError: If robots.txt stays unretrievable; scraping without it is refused.
        """
        robots_url = f"{origin}/robots.txt"
        problem = "no response"
        for attempt in range(self.config.max_retries + 1):
            try:
                logger.debug(f"Requesting {robots_url} (attempt {attempt + 1})")
                response = self.session.get(robots_url, timeout=self.config.timeout)
            except requests.RequestException as e:
                problem = str(e)
            else:
                parser = robots_policy_from_response(
                    response.status_code, response.content.decode("utf-8", errors="replace"), robots_url
                )
                if parser is not None:
                    logger.info(f"Loaded {robots_url} (HTTP {response.status_code})")
                    return parser
                problem = f"HTTP {response.status_code}"
            if attempt < self.config.max_retries:
                time.sleep(2**attempt)
        raise RobotsDisallowedError(
            f"{robots_url} could not be retrieved ({problem}); refusing to scrape without it "
            "(set check_robots=False to opt out)"
        )

    def _require_allowed(self, url: str) -> None:
        """Raise :class:`RobotsDisallowedError` unless robots.txt permits ``url`` (no-op without ``check_robots``)."""
        if not self.config.check_robots:
            return
        parser = self._robots_for(url)
        if not parser.can_fetch(self.config.user_agent, url):
            raise RobotsDisallowedError(
                f"robots.txt of {_origin(url)} disallows {url} for user agent {self.config.user_agent!r}"
            )

    def can_fetch(self, url: str) -> bool:
        """Whether robots.txt permits fetching ``url`` (always True when ``check_robots`` is false).

        Downloads the site's robots.txt on first use. A robots.txt that cannot be retrieved counts as "not permitted".
        """
        try:
            self._require_allowed(url)
        except RobotsDisallowedError:
            return False
        return True

    def crawl_delay(self, url: str) -> Optional[float]:
        """Seconds robots.txt asks to wait between requests, or ``None`` (also when ``check_robots`` is false).

        Downloads the site's robots.txt on first use.
        """
        if not self.config.check_robots:
            return None
        delay = self._robots_for(url).crawl_delay(self.config.user_agent)
        return None if delay is None else float(delay)

    def effective_delay(self, url: str) -> float:
        """Seconds to wait between requests to ``url``'s site: the larger of configured delay and ``Crawl-delay``.

        Downloads the site's robots.txt on first use (unless ``check_robots`` is false).
        """
        return max(self.config.delay_between_requests, self.crawl_delay(url) or 0.0)

    def _get(self, url: str) -> requests.Response:
        """GET ``url`` with robots.txt enforced on it and on every redirect hop.

        Redirects are followed here rather than by ``requests`` so that a permitted URL cannot forward the scraper to a
        disallowed one, on this host or another.

        Raises:
            RobotsDisallowedError: If robots.txt forbids ``url`` or a URL it redirects to.
            NetworkError: If the redirect chain is longer than ``_MAX_REDIRECTS``.
        """
        current = url
        for _ in range(_MAX_REDIRECTS + 1):
            first_contact = self.config.check_robots and _origin(current) not in self._robots
            self._require_allowed(current)
            if first_contact:
                # robots.txt was just downloaded from this host: that was a request too
                self._respectful_delay(current)
            response = self.session.get(current, timeout=self.config.timeout, allow_redirects=False)
            target = _redirect_target(current, response.status_code, response.headers.get("Location"))
            if target is None:
                return response
            logger.debug(f"Following redirect {current} -> {target}")
            self._respectful_delay(current)
            current = target
        raise errors.NetworkError(f"Too many redirects fetching {url}")

    def _make_request(self, url: str) -> Optional[BeautifulSoup]:
        """Make an HTTP request with retry logic.

        Args:
            url: URL to request

        Returns:
            BeautifulSoup object if successful, None if failed

        Raises:
            RobotsDisallowedError: If robots.txt forbids the URL (never retried).
        """
        for attempt in range(self.config.max_retries + 1):
            try:
                logger.debug(f"Requesting {url} (attempt {attempt + 1})")
                response = self._get(url)

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
                        self._backoff(attempt, url)  # Exponential backoff
                        continue
                else:
                    logger.warning(f"HTTP {response.status_code} for {url}")
                    response.raise_for_status()

            except requests.HTTPError:
                raise  # Let HTTP errors propagate immediately (403, etc.)
            except requests.RequestException as e:
                logger.warning(f"Request failed for {url}: {e}")
                if attempt < self.config.max_retries:
                    self._backoff(attempt, url)
                    continue

        return None

    def _backoff_seconds(self, attempt: int, url: str) -> float:
        """Seconds to wait before retry ``attempt + 1``: exponential backoff, never below the site's ``Crawl-delay``."""
        return max(float(2**attempt), self.crawl_delay(url) or 0.0)

    def _backoff(self, attempt: int, url: str) -> None:
        """Sleep before retrying ``url``."""
        time.sleep(self._backoff_seconds(attempt, url))

    def _respectful_delay(self, url: Optional[str] = None) -> None:
        """Wait between requests: the configured delay, raised to the site's ``Crawl-delay`` when ``url`` is given."""
        time.sleep(self.config.delay_between_requests if url is None else self.effective_delay(url))

    def _species_url(self, species_name_or_url: str) -> str:
        """Absolute page URL for a species name (e.g. ``Formica_rufa``) or a full URL."""
        if species_name_or_url.startswith("http"):
            return species_name_or_url
        if self.config.base_url.endswith("/"):
            return f"{self.config.base_url}{species_name_or_url}"
        return f"{self.config.base_url}/{species_name_or_url}"

    def _species_label(self, url: str) -> str:
        """Short identifier for a species URL that can be passed back as ``species=``: its page name, else the URL."""
        base = self.config.base_url if self.config.base_url.endswith("/") else f"{self.config.base_url}/"
        if url.startswith(base) and len(url) > len(base):
            return url[len(base) :]
        return url

    def scrape_species_page(self, species_name_or_url: str) -> Optional[Dict[str, Any]]:
        """Scrape a single species page from AntWiki.

        Args:
            species_name_or_url: Species name (e.g. "Camponotus_pennsylvanicus") or full URL

        Returns:
            Dictionary with scraped data, or None if scraping failed

        Raises:
            RobotsDisallowedError: If ``check_robots`` is true and robots.txt forbids the page.
            NetworkError: If the page could not be fetched.
        """
        species_url = self._species_url(species_name_or_url)

        if species_url in self.scraped_urls:
            logger.debug(f"Already scraped: {species_url}")
            return None

        soup = self._make_request(species_url)
        if not soup:
            raise errors.NetworkError(f"Failed to fetch species page: {species_url}")

        self.scraped_urls.add(species_url)
        self._respectful_delay(species_url)

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
        morphology: Dict[str, Any] = {}

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

    def _extract_section_content(self, header_element: Any) -> str:
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
        taxonomy = data.get("taxonomy", {})
        tax_score = 0.0
        if data.get("genus"):
            tax_score += 0.5
        if taxonomy.get("subfamily"):
            tax_score += 0.3
        if taxonomy.get("tribe"):
            tax_score += 0.2
        score_components.append(tax_score * 0.3)

        # Morphological data (30%)
        morph_score = 0.0
        morphology = data.get("morphology", {})
        if morphology.get("body_length_mm"):
            morph_score += 0.4
        if morphology.get("head_width_mm"):
            morph_score += 0.4
        if morphology.get("color"):
            morph_score += 0.2
        score_components.append(morph_score * 0.3)

        # Behavioral data (20%)
        behavior_score = 0.0
        behavior = data.get("behavior", {})
        if behavior.get("foraging_strategy"):
            behavior_score += 0.5
        if behavior.get("nest_type"):
            behavior_score += 0.3
        if behavior.get("colony_size_category"):
            behavior_score += 0.2
        score_components.append(behavior_score * 0.2)

        # Ecological data (20%)
        ecology_score = 0.0
        ecology = data.get("ecology", {})
        if ecology.get("habitat"):
            ecology_score += 0.6
        distribution = data.get("distribution", {})
        if distribution.get("regions"):
            ecology_score += 0.4
        score_components.append(ecology_score * 0.2)

        return sum(score_components)

    def scrape_species_list(self, limit: Optional[int] = None) -> Iterator[str]:
        """Yield species page URLs from the AntWiki species index.

        Species-index retrieval is **not implemented**. How AntWiki exposes its species directory (category pages,
        ``Special:AllPages``, a sitemap) has not been verified against the live site, and this method no longer
        substitutes a fixed five-species sample for it: that stand-in made "all species" silently mean five species.

        Args:
            limit: Accepted for API compatibility; unused.

        Raises:
            NotImplementedError: Always. Name the species instead: ``scrape_all_species(species=[...])``,
                ``bulk_scrape_species(urls, output_dir)`` or ``scrape_species_page(name)``.
        """
        raise NotImplementedError(
            "AntWiki species-index retrieval is not implemented, so there is no list of all species to scrape. "
            "Pass the species explicitly: scrape_all_species(species=['Formica_rufa', ...]), "
            "bulk_scrape_species(urls, output_dir), or scrape_species_page(name) for a single species."
        )

    def bulk_scrape_species(self, species_urls: List[str], output_dir: str | Path) -> List[Dict[str, Any]]:
        """Scrape multiple species pages and save results.

        Pages that fail (including those robots.txt forbids) are logged and skipped; use
        :meth:`scrape_all_species` when the failures need to be reported.

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
                    io.dump_json(data, _species_output_path(output_dir, data), indent=2)

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
            RobotsDisallowedError: If ``check_robots`` is true and robots.txt forbids the URL
            NetworkError: If request fails
        """
        for attempt in range(self.config.max_retries + 1):
            try:
                logger.debug(f"Fetching {url} (attempt {attempt + 1})")
                response = self._get(url)
                response.raise_for_status()
                self._respectful_delay(url)
                return str(response.text)
            except requests.HTTPError as e:
                raise errors.NetworkError(f"HTTP error fetching {url}: {e}") from e
            except requests.RequestException as e:
                if attempt < self.config.max_retries:
                    self._backoff(attempt, url)
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
        """Get the species page URLs listed in the AntWiki species index.

        Returns page URLs, not bare species names. Species-index retrieval is not implemented, so this raises
        ``NotImplementedError`` (see :meth:`scrape_species_list`).

        Args:
            limit: Accepted for API compatibility; unused.

        Returns:
            List of species page URLs (never returned today: the call raises).
        """
        return list(self.scrape_species_list(limit=limit))

    def _species_urls_to_scrape(self, species: Optional[Iterable[str]], limit: Optional[int]) -> List[str]:
        """Resolve the species to scrape into distinct absolute URLs, in order, capped at ``limit``."""
        if species is None:
            return list(self.scrape_species_list(limit=limit))  # raises NotImplementedError
        if isinstance(species, str):
            species = [species]  # a bare name is one species, not a sequence of characters
        if limit is not None:
            validation.validate_range(limit, min_val=0, name="limit")
        urls = list(dict.fromkeys(self._species_url(name) for name in species))
        return urls if limit is None else urls[:limit]

    def scrape_all_species(
        self,
        output_dir: str | Path | None = None,
        limit: Optional[int] = None,
        species: Optional[Iterable[str]] = None,
    ) -> Dict[str, Any]:
        """Scrape the given species and save them to the output directory.

        Each species is written to ``<output_dir>/species/<name>.json`` and all results are combined in
        ``<output_dir>/all_species.json`` (rewritten on every run with that run's results). There is no checkpoint or
        resume: a re-run fetches every species again. Failures are reported in the result instead of aborting the
        run, so ``failed_species`` can be passed back as ``species=`` to retry only those.

        Args:
            output_dir: Directory to save results (default: ``config.output_dir``)
            limit: Maximum number of species to scrape
            species: Species to scrape, as page names (``"Formica_rufa"``) or full URLs; duplicates are scraped once.
                When omitted the AntWiki species index would be used, but index retrieval is not implemented
                (see :meth:`scrape_species_list`), so omitting it raises ``NotImplementedError``.

        Returns:
            Statistics dictionary: ``total`` species attempted, ``completed`` and ``failed`` counts, and
            ``failed_species``, the page name (or URL, for pages outside ``base_url``) of every species that failed,
            for example because robots.txt forbids it, the fetch failed, or the page could not be parsed.

        Raises:
            ValueError: If no output directory is given or configured.
            NotImplementedError: If ``species`` is omitted.
        """
        if output_dir is None:
            output_dir = self.config.output_dir
        if output_dir is None:
            raise ValueError("No output directory specified")

        output_dir = Path(output_dir)
        species_urls = self._species_urls_to_scrape(species, limit)

        species_dir = output_dir / "species"
        species_dir.mkdir(parents=True, exist_ok=True)

        results: List[Dict[str, Any]] = []
        failed_species: List[str] = []

        for url in species_urls:
            # A batch run always fetches afresh. scrape_species_page() skips URLs this scraper already scraped, which
            # would turn a re-run (or a retry of failed_species) into silent "failures" and an emptied all_species.json.
            self.scraped_urls.discard(url)
            try:
                data = self.scrape_species_page(url)
                if data:
                    io.dump_json(data, _species_output_path(species_dir, data), indent=2)
                    results.append(data)
                else:
                    logger.error(f"No data scraped for {url}: the page could not be parsed")
                    failed_species.append(self._species_label(url))
            except Exception as e:
                logger.error(f"Failed to scrape {url}: {e}")
                failed_species.append(self._species_label(url))

        # Save combined file
        all_species_file = output_dir / "all_species.json"
        io.dump_json(results, all_species_file, indent=2)

        stats: Dict[str, Any] = {
            "total": len(species_urls),
            "completed": len(results),
            "failed": len(failed_species),
            "failed_species": failed_species,
        }

        logger.info(f"Scraped {stats['completed']}/{stats['total']} species")
        if failed_species:
            logger.warning(f"{len(failed_species)} species failed: {', '.join(failed_species[:10])}")
        return stats

    def close(self) -> None:
        """Close the scraper session."""
        self.session.close()
        logger.info("AntWiki scraper session closed")


def load_scraper_config(config_path: str | Path | None = None) -> AntWikiScraperConfig:
    """Load scraper configuration from a JSON file, or return defaults.

    The file holds ``AntWikiScraperConfig`` keyword arguments (for example ``{"delay_seconds": 3, "check_robots":
    false}``). Only JSON is read: the YAML file under ``config/phenotype/`` cannot be loaded by this function.

    Args:
        config_path: Path to a JSON configuration file (optional; returns defaults if None)

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
        "check_robots": config.check_robots,
    }

    io.dump_json(config_dict, Path(output_path), indent=2)
    logger.info(f"Created default scraper config at {output_path}")
