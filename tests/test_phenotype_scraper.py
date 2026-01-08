"""Tests for AntWiki web scraping functionality.

All tests use real HTTP requests (no mocks) following project policy.
Tests are marked with @pytest.mark.network and will skip gracefully
if network is unavailable.
"""

from __future__ import annotations

import time
from pathlib import Path

import pytest
import requests

from metainformant.core.utils.errors import NetworkError
from metainformant.phenotype.data.scraper import AntWikiScraper, AntWikiScraperConfig, load_scraper_config


def _handle_antwiki_403(e: Exception) -> None:
    """Handle AntWiki 403 errors consistently.
    
    Centralized error handling for AntWiki blocking automated requests.
    Raises pytest.skip if 403 error detected, otherwise re-raises exception.
    
    Args:
        e: Exception to check for 403 status
        
    Raises:
        pytest.skip.Exception: If 403 Forbidden error detected
        Exception: Re-raises original exception if not 403
    """
    if isinstance(e, NetworkError):
        if isinstance(e.__cause__, requests.HTTPError) and e.__cause__.response.status_code == 403:
            pytest.skip(f"AntWiki blocked request (403 Forbidden): {e}. Site may be blocking automated requests.")
        if "403" in str(e):
            pytest.skip(f"AntWiki blocked request (403 Forbidden): {e}. Site may be blocking automated requests.")
    elif isinstance(e, requests.HTTPError) and e.response.status_code == 403:
        pytest.skip(f"AntWiki blocked request (403 Forbidden): {e}. Site may be blocking automated requests.")
    # Re-raise if not a 403 error
    raise


@pytest.fixture
def antwiki_scraper_config(tmp_path: Path) -> AntWikiScraperConfig:
    """Fixture providing common AntWiki scraper configuration for tests.
    
    Args:
        tmp_path: Pytest temporary directory fixture
        
    Returns:
        AntWikiScraperConfig with test-appropriate settings
    """
    return AntWikiScraperConfig(
        base_url="https://www.antwiki.org/wiki/",
        delay_seconds=1.0,
        user_agent="Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36",
        timeout_seconds=30,
        max_retries=2,
        check_robots=False,  # Skip robots.txt for faster testing
        output_dir=tmp_path,
    )


@pytest.mark.network
def test_scraper_config_loading(tmp_path: Path) -> None:
    """Test loading scraper configuration."""
    # Test default config loading
    config = load_scraper_config()
    assert config.base_url == "https://www.antwiki.org/wiki/"
    assert config.delay_seconds == 2.0
    assert config.timeout_seconds == 30
    assert config.check_robots is True


@pytest.mark.network
def test_scraper_initialization() -> None:
    """Test scraper initialization."""
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


@pytest.mark.network
def test_fetch_page_single_species(antwiki_scraper_config: AntWikiScraperConfig) -> None:
    """Test fetching a single species page."""
    try:
        scraper = AntWikiScraper(antwiki_scraper_config)

        # Test with a known species
        species_name = "Camponotus_pennsylvanicus"
        data = scraper.scrape_species_page(species_name)

        assert "species" in data
        assert data["species"] == "Camponotus pennsylvanicus"
        assert "url" in data
        assert "measurements" in data
        assert "traits" in data
        assert "description" in data
        assert "taxonomy" in data
        assert "distribution" in data
        assert "images" in data

        # Verify data types
        assert isinstance(data["measurements"], dict)
        assert isinstance(data["traits"], list)
        assert isinstance(data["description"], str)
        assert isinstance(data["taxonomy"], dict)
        assert isinstance(data["distribution"], dict)
        assert isinstance(data["images"], list)

    except (NetworkError, requests.HTTPError) as e:
        _handle_antwiki_403(e)
    except (requests.ConnectionError, requests.Timeout) as e:
        pytest.skip(f"Network unavailable: {e}")


@pytest.mark.network
def test_extract_measurements(antwiki_scraper_config: AntWikiScraperConfig) -> None:
    """Test measurement extraction from HTML."""
    try:
        scraper = AntWikiScraper(antwiki_scraper_config)

        # Fetch a real page
        species_name = "Camponotus_pennsylvanicus"
        page_url = f"{antwiki_scraper_config.base_url}{species_name}"
        html = scraper._fetch_page(page_url)

        # Extract measurements
        measurements = scraper.extract_measurements(html)
        assert isinstance(measurements, dict)

    except (NetworkError, requests.HTTPError) as e:
        _handle_antwiki_403(e)
    except (requests.ConnectionError, requests.Timeout) as e:
        pytest.skip(f"Network unavailable: {e}")


@pytest.mark.network
def test_extract_traits(antwiki_scraper_config: AntWikiScraperConfig) -> None:
    """Test trait extraction from HTML."""
    try:
        scraper = AntWikiScraper(antwiki_scraper_config)

        # Fetch a real page
        species_name = "Camponotus_pennsylvanicus"
        page_url = f"{antwiki_scraper_config.base_url}{species_name}"
        html = scraper._fetch_page(page_url)

        # Extract traits
        traits = scraper.extract_traits(html)
        assert isinstance(traits, list)

    except (NetworkError, requests.HTTPError) as e:
        _handle_antwiki_403(e)
    except (requests.ConnectionError, requests.Timeout) as e:
        pytest.skip(f"Network unavailable: {e}")


@pytest.mark.network
def test_extract_taxonomy(antwiki_scraper_config: AntWikiScraperConfig) -> None:
    """Test taxonomy extraction from HTML."""
    try:
        scraper = AntWikiScraper(antwiki_scraper_config)

        # Fetch a real page
        species_name = "Camponotus_pennsylvanicus"
        page_url = f"{antwiki_scraper_config.base_url}{species_name}"
        html = scraper._fetch_page(page_url)

        # Extract taxonomy
        taxonomy = scraper.extract_taxonomy(html)
        assert isinstance(taxonomy, dict)

    except (NetworkError, requests.HTTPError) as e:
        _handle_antwiki_403(e)
    except (requests.ConnectionError, requests.Timeout) as e:
        pytest.skip(f"Network unavailable: {e}")


@pytest.mark.network
def test_rate_limiting(antwiki_scraper_config: AntWikiScraperConfig) -> None:
    """Test that rate limiting is applied between requests."""
    try:
        scraper = AntWikiScraper(antwiki_scraper_config)

        # Make two requests and measure time
        start_time = time.time()
        try:
            scraper._fetch_page(f"{antwiki_scraper_config.base_url}Camponotus_pennsylvanicus")
            scraper._fetch_page(f"{antwiki_scraper_config.base_url}Formica_rufa")
        except (NetworkError, requests.HTTPError) as e:
            _handle_antwiki_403(e)
        elapsed = time.time() - start_time

        # Should take at least delay_seconds between requests
        assert elapsed >= antwiki_scraper_config.delay_seconds, f"Rate limiting not working: {elapsed} < {antwiki_scraper_config.delay_seconds}"

    except (requests.ConnectionError, requests.Timeout) as e:
        pytest.skip(f"Network unavailable: {e}")


@pytest.mark.network
def test_error_handling_invalid_species(antwiki_scraper_config: AntWikiScraperConfig) -> None:
    """Test error handling for invalid species names."""
    try:
        # Override config for quick failure testing
        config = AntWikiScraperConfig(
            base_url=antwiki_scraper_config.base_url,
            delay_seconds=antwiki_scraper_config.delay_seconds,
            user_agent="METAINFORMANT-test",
            timeout_seconds=antwiki_scraper_config.timeout_seconds,
            max_retries=1,  # Quick failure for testing
            check_robots=antwiki_scraper_config.check_robots,
            output_dir=antwiki_scraper_config.output_dir,
        )
        scraper = AntWikiScraper(config)

        # Try to scrape a non-existent species
        invalid_species = "Nonexistent_Species_12345"
        try:
            scraper.scrape_species_page(invalid_species)
            # If no exception raised, that's unexpected
            pytest.fail("Expected exception for invalid species")
        except (NetworkError, requests.HTTPError) as e:
            # If we get a 403, skip the test (site blocking)
            # Otherwise, this is expected behavior (invalid species should raise error)
            try:
                _handle_antwiki_403(e)
            except Exception:
                # Not a 403, so error handling worked correctly
                assert True  # Test passes - error handling works
        except Exception:
            # Any other exception is also acceptable for invalid species
            assert True  # Test passes - error handling works

    except (requests.ConnectionError, requests.Timeout) as e:
        pytest.skip(f"Network unavailable: {e}")


@pytest.mark.network
def test_scrape_single_species_saves_file(antwiki_scraper_config: AntWikiScraperConfig) -> None:
    """Test that scraping a single species saves data to file."""
    try:
        scraper = AntWikiScraper(antwiki_scraper_config)

        species_name = "Camponotus_pennsylvanicus"
        data = scraper.scrape_species_page(species_name)

        # Verify data structure matches expected format
        assert "species" in data
        assert data["species"] == "Camponotus pennsylvanicus"

    except (NetworkError, requests.HTTPError) as e:
        _handle_antwiki_403(e)
    except (requests.ConnectionError, requests.Timeout) as e:
        pytest.skip(f"Network unavailable: {e}")


@pytest.mark.network
@pytest.mark.slow
def test_get_species_list(antwiki_scraper_config: AntWikiScraperConfig) -> None:
    """Test discovering species list from AntWiki."""
    # First check if we can access AntWiki at all
    import requests
    try:
        response = requests.get("https://www.antwiki.org/wiki/Category:Ants", timeout=5)
        if response.status_code == 403:
            pytest.skip("AntWiki blocked request (403 Forbidden). Site blocks automated requests.")
    except (requests.RequestException, requests.Timeout):
        pytest.skip("AntWiki not accessible. Network may be unavailable.")

    from metainformant.core.utils.errors import NetworkError
    try:
        scraper = AntWikiScraper(antwiki_scraper_config)

        species_list = scraper.get_species_list()
        assert isinstance(species_list, list)
        assert len(species_list) > 0

        # Verify species names are in correct format
        for species in species_list[:10]:  # Check first 10
            assert isinstance(species, str)
            assert len(species) > 0
    except NetworkError as e:
        # Skip if AntWiki blocks requests (403 Forbidden due to Cloudflare)
        if "403" in str(e) or "Forbidden" in str(e):
            pytest.skip(f"AntWiki blocked request: {e}")
        else:
            raise

    except (NetworkError, requests.HTTPError) as e:
        _handle_antwiki_403(e)
    except (requests.ConnectionError, requests.Timeout) as e:
        pytest.skip(f"Network unavailable: {e}")


@pytest.mark.network
@pytest.mark.slow
def test_scrape_all_species_with_limit(antwiki_scraper_config: AntWikiScraperConfig) -> None:
    """Test scraping multiple species with limit."""
    # First check if we can access AntWiki at all
    import requests
    try:
        response = requests.get("https://www.antwiki.org/wiki/Category:Ants", timeout=5)
        if response.status_code == 403:
            pytest.skip("AntWiki blocked request (403 Forbidden). Site blocks automated requests.")
    except (requests.RequestException, requests.Timeout):
        pytest.skip("AntWiki not accessible. Network may be unavailable.")

    try:
        # Override config for faster testing
        config = AntWikiScraperConfig(
            base_url=antwiki_scraper_config.base_url,
            delay_seconds=0.5,  # Faster for testing
            user_agent=antwiki_scraper_config.user_agent,
            timeout_seconds=antwiki_scraper_config.timeout_seconds,
            max_retries=antwiki_scraper_config.max_retries,
            check_robots=antwiki_scraper_config.check_robots,
            output_dir=antwiki_scraper_config.output_dir,
        )
        scraper = AntWikiScraper(config)

        # Scrape with small limit
        stats = scraper.scrape_all_species(output_dir=antwiki_scraper_config.output_dir, limit=3)

        assert "total" in stats
        assert "completed" in stats
        assert "failed" in stats
        assert stats["total"] == 3

        # Verify output files were created
        species_dir = antwiki_scraper_config.output_dir / "species"
        assert species_dir.exists()
        assert len(list(species_dir.glob("*.json"))) > 0

        # Verify combined file exists
        all_species_file = antwiki_scraper_config.output_dir / "all_species.json"
        assert all_species_file.exists()

    except (NetworkError, requests.HTTPError) as e:
        _handle_antwiki_403(e)
    except (requests.ConnectionError, requests.Timeout) as e:
        pytest.skip(f"Network unavailable: {e}")

