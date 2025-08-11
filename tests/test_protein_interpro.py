from __future__ import annotations


def _check_online(url: str) -> bool:
    try:
        import requests
        requests.get(url, timeout=5)
        return True
    except Exception:
        return False


def test_interpro_fetch_parses_minimal():
    if not _check_online("https://www.ebi.ac.uk/interpro"):  # pragma: no cover - skip offline
        import pytest
        pytest.skip("No network access for InterPro")

    from metainformant.protein.interpro import fetch_interpro_domains

    results = fetch_interpro_domains("P69905")
    assert isinstance(results, list)
    if results:
        # Shape may vary; only assert list type when non-empty
        assert isinstance(results[0], dict)


