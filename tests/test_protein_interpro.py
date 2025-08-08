from __future__ import annotations


def test_interpro_fetch_parses_minimal(monkeypatch):
    from metainformant.protein.interpro import fetch_interpro_domains

    class DummyResponse:
        def raise_for_status(self):
            return None

        def json(self):
            return {"results": [{"metadata": {"accession": "IPR000001"}}]}

    import requests

    def fake_get(url, headers=None, timeout=None):
        assert "interpro" in url
        return DummyResponse()

    monkeypatch.setattr(requests, "get", fake_get)

    results = fetch_interpro_domains("P69905")
    assert isinstance(results, list) and results
    assert results[0]["metadata"]["accession"].startswith("IPR")


