"""Offline regression tests for the AntWiki scraper's robots.txt compliance and bookkeeping.

Nothing here opens a socket. robots.txt content is fed to the scraper through ``AntWikiScraper.load_robots_txt`` (the
same code path that parses a downloaded robots.txt), and refusals happen before any request is sent. The base URL uses
the reserved ``.invalid`` TLD (RFC 2606), so even a regression that tried to connect could never reach a real site.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Callable
from urllib import robotparser

import pytest

from metainformant.core.utils.errors import IOError as MetainformantIOError, NetworkError, ValidationError
from metainformant.phenotype.data.scraper import (
    AntWikiScraper,
    AntWikiScraperConfig,
    RobotsDisallowedError,
    _redirect_target,
    _species_output_path,
    create_default_scraper_config,
    load_scraper_config,
    robots_policy_from_response,
)

BASE_URL = "https://antwiki.invalid/wiki/"

ROBOTS_LINES = [
    "# robots.txt for a MediaWiki site",
    "User-agent: *",
    "Disallow: /wiki/Special:",
    "Disallow: /wiki/Private_",
    "Crawl-delay: 4",
]


@pytest.fixture
def make_scraper(tmp_path: Path):
    """Build scrapers with robots.txt rules preloaded (``robots=None`` leaves them unloaded); closes them afterwards."""
    scrapers: list[AntWikiScraper] = []

    def factory(*, robots=ROBOTS_LINES, **config_overrides) -> AntWikiScraper:
        settings = {
            "base_url": BASE_URL,
            "delay_seconds": 0.1,
            "max_retries": 0,
            "output_dir": tmp_path,
            **config_overrides,
        }
        scraper = AntWikiScraper(AntWikiScraperConfig(**settings))
        scrapers.append(scraper)
        if robots is not None:
            scraper.load_robots_txt(robots)
        return scraper

    yield factory
    for scraper in scrapers:
        scraper.close()


Factory = Callable[..., AntWikiScraper]


class TestRobotsRefusal:
    def test_disallowed_species_page_is_refused_before_any_connection(self, make_scraper: Factory) -> None:
        scraper = make_scraper()
        with pytest.raises(RobotsDisallowedError, match="Special:Random"):
            scraper.scrape_species_page("Special:Random")
        assert scraper.scraped_urls == set()

    def test_refusal_is_a_network_error_subclass(self, make_scraper: Factory) -> None:
        scraper = make_scraper()
        with pytest.raises(NetworkError):
            scraper.scrape_species_page("Private_notes")

    def test_raw_page_retrieval_is_gated_by_robots_too(self, make_scraper: Factory) -> None:
        scraper = make_scraper()
        with pytest.raises(RobotsDisallowedError):
            scraper._fetch_page(BASE_URL + "Special:Random")

    def test_rules_separate_allowed_from_disallowed_urls(self, make_scraper: Factory) -> None:
        scraper = make_scraper()
        assert scraper.can_fetch(BASE_URL + "Formica_rufa") is True
        assert scraper.can_fetch(BASE_URL + "Special:Random") is False
        assert scraper.can_fetch(BASE_URL + "Private_notes") is False
        assert scraper.can_fetch(BASE_URL + "Privately_named_species") is True  # prefix rule, not substring

    def test_user_agent_group_overrides_wildcard_group(self, make_scraper: Factory) -> None:
        rules = [
            "User-agent: MetaInformant",
            "Disallow: /wiki/Secret",
            "",
            "User-agent: *",
            "Disallow:",
        ]
        ours = make_scraper(robots=rules)  # default user agent starts with the MetaInformant product token
        other = make_scraper(robots=rules, user_agent="OtherBot/2.0")
        assert ours.can_fetch(BASE_URL + "Secret_page") is False
        assert ours.can_fetch(BASE_URL + "Formica_rufa") is True
        assert other.can_fetch(BASE_URL + "Secret_page") is True

    def test_pre_parsed_parser_is_accepted_as_is(self, make_scraper: Factory) -> None:
        parser = robotparser.RobotFileParser()
        parser.parse(["User-agent: *", "Disallow: /wiki/Private_"])
        scraper = make_scraper(robots=None)
        assert scraper.load_robots_txt(parser) is parser
        assert scraper.can_fetch(BASE_URL + "Private_notes") is False
        assert scraper.can_fetch(BASE_URL + "Formica_rufa") is True

    def test_robots_text_and_bytes_are_accepted(self, make_scraper: Factory) -> None:
        text = "User-agent: *\r\nDisallow: /wiki/Private_\r\n"
        for robots in (text, text.encode("utf-8")):
            scraper = make_scraper(robots=robots)
            assert scraper.can_fetch(BASE_URL + "Private_notes") is False
            assert scraper.can_fetch(BASE_URL + "Formica_rufa") is True

    def test_rules_are_scoped_to_each_origin(self, make_scraper: Factory) -> None:
        scraper = make_scraper(robots=["User-agent: *", "Disallow: /wiki/"])
        scraper.load_robots_txt(["User-agent: *", "Disallow:"], site_url="https://mirror.invalid/anything")
        assert scraper.can_fetch("https://antwiki.invalid/wiki/Formica_rufa") is False
        assert scraper.can_fetch("https://mirror.invalid/wiki/Formica_rufa") is True
        # a different scheme or port is a different origin, with its own robots.txt
        scraper.load_robots_txt(["User-agent: *", "Disallow: /"], site_url="https://antwiki.invalid:8443/")
        assert scraper.can_fetch("https://antwiki.invalid:8443/wiki/Formica_rufa") is False
        assert scraper.can_fetch("https://ANTWIKI.invalid:443/wiki/Special:Random") is False  # same origin as default

    def test_opt_out_ignores_robots_entirely(self, make_scraper: Factory) -> None:
        scraper = make_scraper(robots=["User-agent: *", "Disallow: /", "Crawl-delay: 30"], check_robots=False)
        url = BASE_URL + "Special:Random"
        assert scraper.can_fetch(url) is True
        assert scraper.crawl_delay(url) is None
        assert scraper.effective_delay(url) == pytest.approx(0.1)

    def test_opt_out_never_tries_to_load_robots(self, make_scraper: Factory) -> None:
        scraper = make_scraper(robots=None, check_robots=False)
        assert scraper.can_fetch(BASE_URL + "Special:Random") is True
        assert scraper._robots == {}  # nothing was parsed or downloaded

    def test_robots_check_is_on_by_default(self) -> None:
        assert AntWikiScraperConfig().check_robots is True

    def test_redirect_targets_are_judged_by_the_same_rules(self, make_scraper: Factory) -> None:
        scraper = make_scraper()
        target = _redirect_target(BASE_URL + "Old_name", 301, "/wiki/Special:Secret")
        assert target == "https://antwiki.invalid/wiki/Special:Secret"
        assert scraper.can_fetch(target) is False
        assert _redirect_target(BASE_URL + "Old_name", 308, "Formica_rufa") == BASE_URL + "Formica_rufa"
        assert _redirect_target(BASE_URL + "x", 302, "//mirror.invalid/wiki/X") == "https://mirror.invalid/wiki/X"

    @pytest.mark.parametrize(
        ("status", "location"),
        [(200, "/wiki/x"), (301, None), (301, ""), (304, "/wiki/x"), (404, "/wiki/x")],
    )
    def test_non_redirect_responses_have_no_target(self, status: int, location: str | None) -> None:
        assert _redirect_target(BASE_URL + "page", status, location) is None


class TestRobotsResponsePolicy:
    """How a robots.txt HTTP response becomes a policy (the pure half of the download path)."""

    def test_success_parses_the_body(self) -> None:
        parser = robots_policy_from_response(200, "User-agent: *\nDisallow: /wiki/Private_\n")
        assert parser is not None
        assert parser.can_fetch("MetaInformant/1.0", BASE_URL + "Private_notes") is False
        assert parser.can_fetch("MetaInformant/1.0", BASE_URL + "Formica_rufa") is True

    @pytest.mark.parametrize("status", [404, 410, 400])
    def test_missing_robots_means_allow_everything(self, status: int) -> None:
        parser = robots_policy_from_response(status, "")
        assert parser is not None
        assert parser.can_fetch("MetaInformant/1.0", BASE_URL + "Special:Random") is True

    @pytest.mark.parametrize("status", [401, 403])
    def test_refused_robots_means_disallow_everything(self, status: int) -> None:
        parser = robots_policy_from_response(status, "")
        assert parser is not None
        assert parser.can_fetch("MetaInformant/1.0", BASE_URL + "Formica_rufa") is False

    @pytest.mark.parametrize("status", [301, 408, 429, 500, 502, 503, 504])
    def test_unknown_rules_are_not_guessed(self, status: int) -> None:
        assert robots_policy_from_response(status, "User-agent: *\nDisallow:\n") is None


class TestCrawlDelay:
    def test_crawl_delay_raises_the_effective_delay(self, make_scraper: Factory) -> None:
        scraper = make_scraper()  # robots Crawl-delay 4, configured delay 0.1
        assert scraper.crawl_delay(BASE_URL + "Formica_rufa") == 4.0
        assert scraper.effective_delay(BASE_URL + "Formica_rufa") == 4.0

    def test_configured_delay_wins_when_larger(self, make_scraper: Factory) -> None:
        scraper = make_scraper(delay_seconds=6.0)
        assert scraper.effective_delay(BASE_URL + "Formica_rufa") == 6.0

    def test_no_crawl_delay_means_the_configured_delay(self, make_scraper: Factory) -> None:
        scraper = make_scraper(robots=["User-agent: *", "Disallow: /wiki/Private_"], delay_seconds=2.0)
        assert scraper.crawl_delay(BASE_URL + "Formica_rufa") is None
        assert scraper.effective_delay(BASE_URL + "Formica_rufa") == 2.0

    def test_user_agent_specific_crawl_delay_wins(self, make_scraper: Factory) -> None:
        rules = ["User-agent: MetaInformant", "Crawl-delay: 7", "", "User-agent: *", "Crawl-delay: 3"]
        assert make_scraper(robots=rules).effective_delay(BASE_URL + "x") == 7.0
        assert make_scraper(robots=rules, user_agent="OtherBot/2.0").effective_delay(BASE_URL + "x") == 3.0

    def test_fractional_crawl_delay_is_rounded_up_not_dropped(self, make_scraper: Factory) -> None:
        lines = ["User-agent: *", "Crawl-delay: 2.5  # please be gentle"]
        stdlib = robotparser.RobotFileParser()
        stdlib.parse(lines)
        assert stdlib.crawl_delay("MetaInformant") is None  # the stdlib silently ignores it
        scraper = make_scraper(robots=lines)
        assert scraper.crawl_delay(BASE_URL + "x") == 3.0
        assert scraper.effective_delay(BASE_URL + "x") == 3.0

    def test_retry_backoff_never_undercuts_the_crawl_delay(self, make_scraper: Factory) -> None:
        scraper = make_scraper()  # Crawl-delay 4
        url = BASE_URL + "Formica_rufa"
        assert scraper._backoff_seconds(0, url) == 4.0  # plain backoff would be 1 s
        assert scraper._backoff_seconds(1, url) == 4.0
        assert scraper._backoff_seconds(3, url) == 8.0  # longer backoff still applies

    def test_backoff_without_robots_is_plain_exponential(self, make_scraper: Factory) -> None:
        scraper = make_scraper(check_robots=False)
        assert [scraper._backoff_seconds(n, BASE_URL + "x") for n in range(4)] == [1.0, 2.0, 4.0, 8.0]

    def test_respectful_delay_actually_waits_for_the_crawl_delay(self, make_scraper: Factory) -> None:
        scraper = make_scraper(robots=["User-agent: *", "Crawl-delay: 1"], delay_seconds=0.1)
        started = time.monotonic()
        scraper._respectful_delay(BASE_URL + "Formica_rufa")
        assert time.monotonic() - started >= 0.99

    def test_respectful_delay_without_crawl_delay_stays_short(self, make_scraper: Factory) -> None:
        scraper = make_scraper(robots=["User-agent: *", "Disallow: /wiki/Private_"], delay_seconds=0.1)
        started = time.monotonic()
        scraper._respectful_delay(BASE_URL + "Formica_rufa")
        assert time.monotonic() - started < 0.9


class TestConfigAliases:
    def test_delay_alias_is_assignable_and_stays_in_sync(self) -> None:
        config = AntWikiScraperConfig(delay_seconds=1.0)
        config.delay_seconds = 5.0
        assert config.delay_between_requests == 5.0
        config.delay_between_requests = 7.0
        assert config.delay_seconds == 7.0

    def test_timeout_alias_is_assignable_and_stays_in_sync(self) -> None:
        config = AntWikiScraperConfig(timeout_seconds=10)
        config.timeout_seconds = 20
        assert config.timeout == 20
        config.timeout = 40
        assert config.timeout_seconds == 40

    def test_assignment_is_validated_and_a_rejected_value_changes_nothing(self) -> None:
        config = AntWikiScraperConfig(delay_seconds=2.0, timeout_seconds=30)
        for bad_delay in (0, -1, 0.05):
            with pytest.raises(ValidationError):
                config.delay_seconds = bad_delay
            with pytest.raises(ValidationError):
                config.delay_between_requests = bad_delay
        for bad_timeout in (0, -5):
            with pytest.raises(ValidationError):
                config.timeout_seconds = bad_timeout
        assert config.delay_seconds == 2.0
        assert config.timeout_seconds == 30

    def test_constructor_spellings_still_work_and_seconds_wins(self) -> None:
        assert AntWikiScraperConfig(delay_between_requests=3.0).delay_seconds == 3.0
        assert AntWikiScraperConfig(delay_seconds=1.5, delay_between_requests=3.0).delay_between_requests == 1.5
        assert AntWikiScraperConfig(timeout=12).timeout_seconds == 12
        with pytest.raises(ValidationError):
            AntWikiScraperConfig(delay_seconds=0.0)

    def test_scraper_follows_a_delay_changed_after_construction(self, make_scraper: Factory) -> None:
        scraper = make_scraper(check_robots=False)
        scraper.config.delay_seconds = 3.0
        assert scraper.effective_delay(BASE_URL + "Formica_rufa") == 3.0


class TestConfigFile:
    def test_opt_out_is_available_from_a_config_file(self, tmp_path: Path) -> None:
        path = tmp_path / "scraper.json"
        path.write_text(json.dumps({"check_robots": False, "delay_seconds": 3}))
        config = load_scraper_config(path)
        assert config.check_robots is False
        assert config.delay_seconds == 3.0

    def test_generated_default_config_shows_the_robots_setting(self, tmp_path: Path) -> None:
        path = tmp_path / "default.json"
        create_default_scraper_config(path)
        data = json.loads(path.read_text())
        assert data["check_robots"] is True
        assert load_scraper_config(path).check_robots is True
        data["check_robots"] = False
        path.write_text(json.dumps(data))
        assert load_scraper_config(path).check_robots is False

    def test_only_json_is_read_so_a_yaml_file_fails_loudly(self, tmp_path: Path) -> None:
        path = tmp_path / "scraper.yaml"
        path.write_text("check_robots: false\ndelay_seconds: 3.0\n")
        with pytest.raises(MetainformantIOError):
            load_scraper_config(path)

    def test_no_path_returns_the_defaults(self) -> None:
        config = load_scraper_config(None)
        assert config.check_robots is True
        assert config.delay_seconds == 2.0
        assert config.base_url == "https://www.antwiki.org/wiki/"


class TestSpeciesSelection:
    def test_species_index_retrieval_is_an_explicit_error(self, make_scraper: Factory) -> None:
        scraper = make_scraper()
        with pytest.raises(NotImplementedError, match="species-index retrieval is not implemented"):
            scraper.scrape_species_list()
        with pytest.raises(NotImplementedError, match="species-index retrieval is not implemented"):
            scraper.get_species_list(limit=3)

    def test_no_canned_sample_is_ever_returned(self, make_scraper: Factory) -> None:
        scraper = make_scraper()
        for call in (scraper.scrape_species_list, scraper.get_species_list):
            with pytest.raises(NotImplementedError):
                list(call())

    def test_scrape_all_species_without_species_raises_and_writes_nothing(
        self, make_scraper: Factory, tmp_path: Path
    ) -> None:
        scraper = make_scraper()
        with pytest.raises(NotImplementedError, match="species="):
            scraper.scrape_all_species(output_dir=tmp_path / "out")
        assert not (tmp_path / "out").exists()

    def test_names_and_urls_resolve_to_distinct_page_urls_in_order(self, make_scraper: Factory) -> None:
        scraper = make_scraper()
        urls = scraper._species_urls_to_scrape(
            ["Lasius_niger", "Formica_rufa", BASE_URL + "Lasius_niger", BASE_URL + "Formica_rufa"], None
        )
        assert urls == [BASE_URL + "Lasius_niger", BASE_URL + "Formica_rufa"]

    def test_limit_caps_explicit_species(self, make_scraper: Factory) -> None:
        scraper = make_scraper()
        assert scraper._species_urls_to_scrape(["A_a", "B_b", "C_c"], 2) == [BASE_URL + "A_a", BASE_URL + "B_b"]
        assert scraper._species_urls_to_scrape(["A_a", "B_b"], 0) == []

    def test_negative_limit_is_rejected(self, make_scraper: Factory) -> None:
        with pytest.raises(ValidationError):
            make_scraper()._species_urls_to_scrape(["A_a"], -1)

    def test_a_single_species_string_is_not_split_into_characters(self, make_scraper: Factory) -> None:
        assert make_scraper()._species_urls_to_scrape("Formica_rufa", None) == [BASE_URL + "Formica_rufa"]

    def test_base_url_without_trailing_slash_is_handled(self, make_scraper: Factory) -> None:
        scraper = make_scraper(base_url="https://antwiki.invalid/wiki")
        assert scraper._species_url("Formica_rufa") == "https://antwiki.invalid/wiki/Formica_rufa"
        assert scraper._species_label("https://antwiki.invalid/wiki/Formica_rufa") == "Formica_rufa"


class TestFailureTracking:
    """``scrape_all_species`` reports failed species; a robots.txt refusal is a deterministic offline failure."""

    BLOCKED = ["User-agent: *", "Disallow: /wiki/Blocked_"]

    def test_failed_species_names_every_refused_species(self, make_scraper: Factory, tmp_path: Path) -> None:
        scraper = make_scraper(robots=self.BLOCKED)
        stats = scraper.scrape_all_species(species=["Blocked_a", "Blocked_b"])
        assert stats == {"total": 2, "completed": 0, "failed": 2, "failed_species": ["Blocked_a", "Blocked_b"]}
        assert json.loads((tmp_path / "all_species.json").read_text()) == []
        assert list((tmp_path / "species").glob("*.json")) == []

    def test_failed_count_always_matches_failed_species(self, make_scraper: Factory) -> None:
        scraper = make_scraper(robots=self.BLOCKED)
        stats = scraper.scrape_all_species(species=["Blocked_a", "Blocked_b", "Blocked_c"], limit=2)
        assert stats["total"] == 2
        assert stats["failed"] == len(stats["failed_species"]) == 2

    def test_failed_species_uses_the_full_url_outside_base_url(self, make_scraper: Factory) -> None:
        scraper = make_scraper(robots=self.BLOCKED)
        scraper.load_robots_txt(["User-agent: *", "Disallow: /"], site_url="https://mirror.invalid/")
        stats = scraper.scrape_all_species(species=["https://mirror.invalid/wiki/Formica_rufa", "Blocked_a"])
        assert stats["failed_species"] == ["https://mirror.invalid/wiki/Formica_rufa", "Blocked_a"]

    def test_failed_names_can_be_passed_back_as_species(self, make_scraper: Factory) -> None:
        scraper = make_scraper(robots=self.BLOCKED)
        first = scraper.scrape_all_species(species=["Blocked_a", "Blocked_b"])
        second = scraper.scrape_all_species(species=first["failed_species"])
        assert second["failed_species"] == first["failed_species"]
        assert second["total"] == 2

    def test_duplicate_species_are_attempted_once(self, make_scraper: Factory) -> None:
        scraper = make_scraper(robots=self.BLOCKED)
        stats = scraper.scrape_all_species(species=["Blocked_a", BASE_URL + "Blocked_a", "Blocked_a"])
        assert stats["total"] == 1
        assert stats["failed_species"] == ["Blocked_a"]

    def test_batch_run_does_not_treat_a_stale_session_entry_as_done(self, make_scraper: Factory) -> None:
        scraper = make_scraper(robots=self.BLOCKED)
        scraper.scraped_urls.add(BASE_URL + "Blocked_a")  # as if an earlier run had scraped it
        scraper.scrape_all_species(species=["Blocked_a"])
        # the species was attempted afresh (and refused) rather than skipped as "already scraped"
        assert BASE_URL + "Blocked_a" not in scraper.scraped_urls

    def test_no_output_directory_is_still_an_error(self, make_scraper: Factory) -> None:
        scraper = make_scraper(robots=self.BLOCKED, output_dir=None)
        with pytest.raises(ValueError, match="No output directory"):
            scraper.scrape_all_species(species=["Blocked_a"])


class TestOutputFilenames:
    """The file name comes from the scraped page title, which is remote input."""

    @pytest.mark.parametrize(
        "title",
        ["../../escape", "/etc/cron.d/evil", "a/b/c", "..", "", "...", "C:\\Windows\\evil", "x/../../y", "a\x00b/../c"],
    )
    def test_remote_title_cannot_escape_the_species_directory(self, tmp_path: Path, title: str) -> None:
        species_dir = tmp_path / "species"
        species_dir.mkdir()
        path = _species_output_path(species_dir, {"species_name": title, "species": title})
        assert path.parent == species_dir
        assert path.suffix == ".json"
        assert path.resolve().parent == species_dir.resolve()

    def test_ordinary_species_names_are_unchanged(self, tmp_path: Path) -> None:
        path = _species_output_path(tmp_path, {"species_name": "Camponotus pennsylvanicus"})
        assert path == tmp_path / "Camponotus_pennsylvanicus.json"

    def test_missing_name_falls_back_to_species_then_unknown(self, tmp_path: Path) -> None:
        assert _species_output_path(tmp_path, {"species": "Formica rufa"}).name == "Formica_rufa.json"
        assert _species_output_path(tmp_path, {}).name == "unknown.json"
