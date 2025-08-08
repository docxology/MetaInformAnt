from __future__ import annotations

from pathlib import Path


def test_uniprot_id_mapping_builds_url_and_parses_minimal_response(monkeypatch):
    from metainformant.protein.uniprot import map_ids_uniprot

    # Fake network call
    class DummyResponse:
        def __init__(self, json_obj):
            self._json = json_obj

        def raise_for_status(self):
            return None

        def json(self):
            return self._json

    called = {}

    def fake_post(url, data=None, headers=None, timeout=None):
        called["url"] = url
        assert "idmapping/run" in url
        assert data["from"] == "UniProtKB_AC-ID"
        assert data["to"] == "UniProtKB"
        return DummyResponse({"jobId": "abc123"})

    def fake_get(url, headers=None, timeout=None):
        called["get_url"] = url
        if url.endswith("status/abc123"):
            return DummyResponse({"jobStatus": "FINISHED"})
        else:
            return DummyResponse({"results": [{"from": "P69905", "to": {"primaryAccession": "P69905"}}]})

    import requests

    monkeypatch.setattr(requests, "post", fake_post)
    monkeypatch.setattr(requests, "get", fake_get)

    result = map_ids_uniprot(["P69905"])  # no-op mapping but exercises flow
    assert result == {"P69905": "P69905"}
    assert "idmapping/run" in called["url"]


def test_pdb_download_builds_paths_and_writes(tmp_path: Path, monkeypatch):
    from metainformant.protein.pdb import fetch_pdb_structure

    # Fake network stream write
    class DummyResponse:
        def __init__(self, text):
            self.text = text

        def raise_for_status(self):
            return None

        def iter_content(self, chunk_size=8192):
            yield self.text.encode()

    import requests

    def fake_get(url, stream=False, timeout=None):
        assert "download/pdb" in url or url.endswith(".pdb") or url.endswith(".cif")
        return DummyResponse("MODEL 1\nENDMDL\nEND\n")

    monkeypatch.setattr(requests, "get", fake_get)

    out = fetch_pdb_structure("1CRN", tmp_path, fmt="pdb")
    assert out.exists() and out.suffix == ".pdb"
    assert out.read_text().startswith("MODEL")


