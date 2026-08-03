"""Tests for AntWiki phenotype scraping and parsing.

The scraper still uses real ``requests`` and BeautifulSoup code paths, but the
test pages are served by a local HTTP server so the suite is reproducible and
does not depend on a third-party site's availability or anti-bot policy.
"""

from __future__ import annotations

import functools
import http.server
import threading
import time
from pathlib import Path

import pytest

from metainformant.core.utils.errors import NetworkError
from metainformant.phenotype.data.scraper import AntWikiScraper, AntWikiScraperConfig, load_scraper_config

_SPECIES = [
    "Camponotus_pennsylvanicus",
    "Formica_rufa",
    "Lasius_niger",
    "Myrmica_rubra",
    "Pheidole_pallidula",
]


def _page_html(species_slug: str) -> str:
    species = species_slug.replace("_", " ")
    return f"""<!doctype html>
<html><head><title>{species} - AntWiki</title></head>
<body>
  <div class="mw-parser-output"><p>{species} is a locally served ant species page.</p></div>
  <div class="taxonomy"><p>Subfamily: Formicinae Tribe: Camponotini</p></div>
  <h2>Morphology</h2>
  <h3>Body</h3><p>Body length 12 mm.</p>
  <h3>Head</h3><p>Head width 2 mm.</p>
  <h3>Color</h3><p>General color black.</p>
  <h2>Behavior</h2>
  <h3>Foraging</h3><p>Foraging strategy generalist.</p>
  <h3>Nest</h3><p>Nest type ground.</p>
  <h2>Ecology</h2><p>Temperate forest habitat.</p>
  <h2>Distribution</h2><p>North America.</p>
  <img src="ant.jpg" />
</body></html>"""


@pytest.fixture
def antwiki_scraper_config(tmp_path: Path) -> AntWikiScraperConfig:
    """Provide local AntWiki-shaped pages through a real HTTP server."""
    root = tmp_path / "antwiki"
    wiki_root = root / "wiki"
    wiki_root.mkdir(parents=True)
    for species in _SPECIES:
        (wiki_root / species).write_text(_page_html(species))

    class QuietHandler(http.server.SimpleHTTPRequestHandler):
        def log_message(self, format: str, *args: object) -> None:
            return

    handler = functools.partial(QuietHandler, directory=str(root))
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield AntWikiScraperConfig(
            base_url=f"http://127.0.0.1:{server.server_port}/wiki/",
            delay_seconds=0.1,
            user_agent="METAINFORMANT-test",
            timeout_seconds=5,
            max_retries=1,
            check_robots=False,
            output_dir=tmp_path,
        )
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


@pytest.mark.network
def test_scraper_config_loading(tmp_path: Path) -> None:
    """Test loading the default AntWiki scraper configuration."""
    config = load_scraper_config()
    assert config.base_url == "https://www.antwiki.org/wiki/"
    assert config.delay_seconds == 2.0
    assert config.timeout_seconds == 30
    assert config.check_robots is True


@pytest.mark.network
def test_scraper_initialization() -> None:
    """Test scraper initialization and configuration aliases."""
    config = AntWikiScraperConfig(
        base_url="https://www.antwiki.org/wiki/",
        delay_seconds=1.0,
        user_agent="test-agent",
        timeout_seconds=10,
        max_retries=2,
        check_robots=False,
        output_dir=Path("output/test"),
    )
    scraper = AntWikiScraper(config)
    assert scraper.config.base_url == "https://www.antwiki.org/wiki/"
    assert scraper.config.delay_seconds == 1.0
    scraper.close()


@pytest.mark.network
def test_fetch_page_single_species(antwiki_scraper_config: AntWikiScraperConfig) -> None:
    """Test fetching and structuring one species page."""
    scraper = AntWikiScraper(antwiki_scraper_config)
    data = scraper.scrape_species_page("Camponotus_pennsylvanicus")

    assert data is not None
    assert data["species"] == "Camponotus pennsylvanicus"
    assert data["url"].startswith("http://127.0.0.1:")
    assert data["measurements"]["body_length_mm"] == 12.0
    assert data["measurements"]["head_width_mm"] == 2.0
    assert data["traits"]
    assert data["description"]
    assert data["taxonomy"]["subfamily"] == "Formicinae"
    assert data["distribution"]["regions"] == ["north america"]
    assert data["images"] == ["ant.jpg"]
    scraper.close()


@pytest.mark.network
def test_extract_measurements(antwiki_scraper_config: AntWikiScraperConfig) -> None:
    """Test morphological measurement extraction from fetched HTML."""
    scraper = AntWikiScraper(antwiki_scraper_config)
    html = scraper._fetch_page(f"{antwiki_scraper_config.base_url}Camponotus_pennsylvanicus")
    measurements = scraper.extract_measurements(html)
    assert measurements == {"body_length_mm": 12.0, "head_width_mm": 2.0, "color": {"general": "black"}}
    scraper.close()


@pytest.mark.network
def test_extract_traits(antwiki_scraper_config: AntWikiScraperConfig) -> None:
    """Test behavioral trait extraction from fetched HTML."""
    scraper = AntWikiScraper(antwiki_scraper_config)
    html = scraper._fetch_page(f"{antwiki_scraper_config.base_url}Camponotus_pennsylvanicus")
    traits = scraper.extract_traits(html)
    assert set(traits) == {"generalist", "ground"}
    scraper.close()


@pytest.mark.network
def test_extract_taxonomy(antwiki_scraper_config: AntWikiScraperConfig) -> None:
    """Test taxonomy extraction from fetched HTML."""
    scraper = AntWikiScraper(antwiki_scraper_config)
    html = scraper._fetch_page(f"{antwiki_scraper_config.base_url}Camponotus_pennsylvanicus")
    taxonomy = scraper.extract_taxonomy(html)
    assert taxonomy == {"subfamily": "Formicinae", "tribe": "Camponotini"}
    scraper.close()


@pytest.mark.network
def test_rate_limiting(antwiki_scraper_config: AntWikiScraperConfig) -> None:
    """Test that the configured delay is applied between requests."""
    scraper = AntWikiScraper(antwiki_scraper_config)
    start_time = time.time()
    scraper._fetch_page(f"{antwiki_scraper_config.base_url}Camponotus_pennsylvanicus")
    scraper._fetch_page(f"{antwiki_scraper_config.base_url}Formica_rufa")
    elapsed = time.time() - start_time
    assert elapsed >= antwiki_scraper_config.delay_seconds
    scraper.close()


@pytest.mark.network
def test_error_handling_invalid_species(antwiki_scraper_config: AntWikiScraperConfig) -> None:
    """Test that a missing species page becomes a structured NetworkError."""
    scraper = AntWikiScraper(antwiki_scraper_config)
    with pytest.raises(NetworkError, match="Failed to fetch species page"):
        scraper.scrape_species_page("Nonexistent_Species_12345")
    scraper.close()


@pytest.mark.network
def test_scrape_single_species_saves_file(antwiki_scraper_config: AntWikiScraperConfig) -> None:
    """Test that a single species scrape returns publication-ready fields."""
    scraper = AntWikiScraper(antwiki_scraper_config)
    data = scraper.scrape_species_page("Camponotus_pennsylvanicus")
    assert data is not None
    assert data["species"] == "Camponotus pennsylvanicus"
    assert data["confidence_score"] > 0
    scraper.close()


@pytest.mark.network
@pytest.mark.slow
def test_get_species_list(antwiki_scraper_config: AntWikiScraperConfig) -> None:
    """Test deterministic species-list discovery."""
    scraper = AntWikiScraper(antwiki_scraper_config)
    species_list = scraper.get_species_list()
    assert len(species_list) == len(_SPECIES)
    assert all(species.startswith(antwiki_scraper_config.base_url) for species in species_list)
    scraper.close()


@pytest.mark.network
@pytest.mark.slow
def test_scrape_all_species_with_limit(antwiki_scraper_config: AntWikiScraperConfig) -> None:
    """Test bounded multi-species scraping and output generation."""
    scraper = AntWikiScraper(antwiki_scraper_config)
    stats = scraper.scrape_all_species(output_dir=antwiki_scraper_config.output_dir, limit=3)

    assert stats == {"total": 3, "completed": 3, "failed": 0}
    species_dir = antwiki_scraper_config.output_dir / "species"
    assert sorted(path.stem for path in species_dir.glob("*.json")) == [
        "Camponotus_pennsylvanicus",
        "Formica_rufa",
        "Lasius_niger",
    ]
    assert (antwiki_scraper_config.output_dir / "all_species.json").exists()
    scraper.close()
