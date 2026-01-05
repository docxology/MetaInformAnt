from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.core.download import download_with_progress


def test_download_with_progress_file_scheme_creates_file_and_heartbeat(tmp_path: Path) -> None:
    src = tmp_path / "src.bin"
    src.write_bytes(b"x" * (128 * 1024))

    dest_dir = tmp_path / "dest"
    dest_dir.mkdir()
    dest = dest_dir / "out.bin"

    result = download_with_progress(
        f"file://{src}",
        dest,
        show_progress=False,
        heartbeat_interval=1,
        max_retries=1,
        chunk_size=16 * 1024,
        timeout=5,
        resume=False,
    )

    assert result.success
    assert dest.exists()
    assert dest.stat().st_size == src.stat().st_size

    hb = dest_dir / ".downloads" / f"{dest.name}.heartbeat.json"
    assert hb.exists()
    # Must contain at least these keys
    text = hb.read_text(encoding="utf-8")
    assert "\"status\"" in text
    assert "\"bytes_downloaded\"" in text


def test_monitor_subprocess_directory_growth_writes_heartbeat(tmp_path: Path) -> None:
    import subprocess
    import sys

    from metainformant.core.download import monitor_subprocess_directory_growth

    out_dir = tmp_path / "out"
    out_dir.mkdir()
    hb = tmp_path / "hb.json"

    # Real subprocess that writes bytes into out_dir over time.
    code = (
        "import time\n"
        "from pathlib import Path\n"
        "p = Path(r'" + str(out_dir) + "') / 'blob.bin'\n"
        "with open(p, 'wb') as f:\n"
        "  for _ in range(5):\n"
        "    f.write(b'x' * 65536)\n"
        "    f.flush()\n"
        "    time.sleep(0.2)\n"
    )
    proc = subprocess.Popen([sys.executable, "-c", code])
    rc, bytes_done = monitor_subprocess_directory_growth(
        process=proc,
        watch_dir=out_dir,
        heartbeat_path=hb,
        total_bytes=5 * 65536,
        heartbeat_interval=1,
        show_progress=False,
        desc="test-monitor",
    )

    assert rc == 0
    assert bytes_done >= 5 * 65536
    assert hb.exists()
    text = hb.read_text(encoding="utf-8")
    assert "\"status\"" in text


def test_monitor_subprocess_file_count_tracks_expected_files(tmp_path: Path) -> None:
    import subprocess
    import sys

    from metainformant.core.download import monitor_subprocess_file_count

    out_dir = tmp_path / "out_files"
    out_dir.mkdir()
    hb = tmp_path / "hb_files.json"

    # Subprocess creates two files over time
    code = (
        "import time\n"
        "from pathlib import Path\n"
        f"out = Path(r'{out_dir}')\n"
        "time.sleep(0.2)\n"
        "(out / 'a.txt').write_text('a')\n"
        "time.sleep(0.4)\n"
        "(out / 'b.txt').write_text('b')\n"
    )
    proc = subprocess.Popen([sys.executable, "-c", code])
    rc = monitor_subprocess_file_count(
        process=proc,
        watch_dir=out_dir,
        heartbeat_path=hb,
        expected_files=['a.txt', 'b.txt'],
        heartbeat_interval=1,
        show_progress=False,
        desc="test-files",
    )

    assert rc == 0
    assert (out_dir / "a.txt").exists()
    assert (out_dir / "b.txt").exists()
    assert hb.exists()
    text = hb.read_text(encoding="utf-8")
    assert "\"progress\"" in text
    assert "\"file_count\"" in text


def test_monitor_subprocess_sample_progress_counts_completion_glob(tmp_path: Path) -> None:
    import subprocess
    import sys

    from metainformant.core.download import monitor_subprocess_sample_progress

    quant_dir = tmp_path / "quant"
    quant_dir.mkdir()
    hb = tmp_path / "hb_samples.json"

    # Subprocess creates sample completion markers like quant/SRR*/abundance.tsv
    code = (
        "import time\n"
        "from pathlib import Path\n"
        f"q = Path(r'{quant_dir}')\n"
        "for i in range(3):\n"
        "  d = q / f'SRR{i:02d}'\n"
        "  d.mkdir(parents=True, exist_ok=True)\n"
        "  time.sleep(0.2)\n"
        "  (d / 'abundance.tsv').write_text('ok')\n"
    )
    proc = subprocess.Popen([sys.executable, "-c", code])
    rc = monitor_subprocess_sample_progress(
        process=proc,
        watch_dir=quant_dir,
        heartbeat_path=hb,
        completion_glob="*/abundance.tsv",
        total_samples=3,
        heartbeat_interval=1,
        show_progress=False,
        desc="test-samples",
    )

    assert rc == 0
    assert hb.exists()
    assert len(list(quant_dir.glob("*/abundance.tsv"))) == 3
    text = hb.read_text(encoding="utf-8")
    assert "\"sample_count\"" in text


@pytest.mark.network
def test_download_with_progress_http_creates_heartbeat(tmp_path: Path) -> None:
    # Use a small, stable endpoint. If network is unavailable, skip.
    url = "https://www.example.com/"
    dest = tmp_path / "example.html"
    try:
        result = download_with_progress(
            url,
            dest,
            show_progress=False,
            heartbeat_interval=1,
            max_retries=1,
            timeout=10,
            chunk_size=8192,
            resume=False,
        )
    except Exception as e:
        pytest.skip(f"Network unavailable or blocked: {e}")

    assert result.success
    assert dest.exists()
    hb = dest.parent / ".downloads" / f"{dest.name}.heartbeat.json"
    assert hb.exists()


