from __future__ import annotations

from pathlib import Path


def test_alphafold_fetch_builds_and_writes(tmp_path: Path, monkeypatch):
    from metainformant.protein.alphafold import build_alphafold_url, fetch_alphafold_model

    url = build_alphafold_url("P69905", version=4, fmt="pdb")
    assert url.endswith("AF-P69905-F1-model_v4.pdb")

    class DummyResponse:
        def raise_for_status(self):
            return None

        def iter_content(self, chunk_size=8192):
            yield b"HEADER\nATOM ...\nEND\n"

    import requests

    def fake_get(url, stream=False, timeout=None):
        return DummyResponse()

    monkeypatch.setattr(requests, "get", fake_get)

    out = fetch_alphafold_model("P69905", tmp_path, version=4, fmt="pdb")
    assert out.exists() and out.suffix == ".pdb"
    assert "HEADER" in out.read_text()


