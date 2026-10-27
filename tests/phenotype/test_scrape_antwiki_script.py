"""The AntWiki scraping script says only what the scraper really does.

The script used to import a module that does not exist, pass ``resume=`` to a method without that parameter, assign a
read-only config alias, read a ``failed_species`` key nothing produced, and print "Check robots.txt: True" while
robots.txt was never consulted. Every case below runs the real script and needs no network: the paths exercised here
stop before any request would be sent.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts" / "phenotype" / "scrape_antwiki.py"


def _run_script(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
        timeout=180,
        check=False,
    )


def test_script_starts_and_documents_its_options() -> None:
    result = _run_script("--help")
    assert result.returncode == 0, result.stderr
    assert "--species" in result.stdout and "--check-robots" in result.stdout


def test_all_flag_exits_with_a_clear_error_about_the_species_index(tmp_path: Path) -> None:
    result = _run_script("--all", "--output", str(tmp_path / "out"))
    assert result.returncode == 2
    assert "species-index retrieval is not implemented" in result.stderr
    assert "--species" in result.stderr
    assert "Traceback" not in result.stderr
    assert not (tmp_path / "out" / "species").exists()  # no empty "results" were produced
    assert not (tmp_path / "out" / "all_species.json").exists()


def test_all_flag_reports_that_robots_txt_is_really_honoured(tmp_path: Path) -> None:
    result = _run_script("--all", "--output", str(tmp_path / "out"))
    assert "Check robots.txt: True (disallowed pages are refused" in result.stdout
    assert "Crawl-delay" in result.stdout


def test_resume_is_rejected_as_unsupported_instead_of_crashing() -> None:
    result = _run_script("--all", "--resume")
    assert result.returncode == 2
    assert "--resume is not supported" in result.stderr
    assert "TypeError" not in result.stderr


def test_help_text_matches_the_implementation() -> None:
    help_text = " ".join(_run_script("--help").stdout.split())
    assert "not implemented" in help_text  # --all
    assert "not supported" in help_text  # --resume
    assert "Crawl-delay" in help_text
    assert "built-in defaults" in help_text  # --config; the YAML named here before cannot be loaded
    assert "antwiki_scraper.yaml" not in help_text


def test_invalid_delay_is_a_usage_error_not_a_traceback() -> None:
    result = _run_script("--all", "--delay", "0")
    assert result.returncode == 2
    assert "--delay" in result.stderr and "must be >= 0.1" in result.stderr
    assert "Traceback" not in result.stderr


def test_delay_override_reaches_the_config(tmp_path: Path) -> None:
    # --delay assigns config.delay_between_requests; the script must echo the value it will really use
    result = _run_script("--all", "--delay", "3.5", "--output", str(tmp_path / "out"))
    assert "Delay between requests: 3.5 seconds" in result.stdout
