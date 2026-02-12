"""Robust, modular download utilities with progress + heartbeat.

This module centralizes downloading logic across METAINFORMANT.

Design goals:
- Heartbeat files for long-running downloads (machine-readable JSON).
- Optional progress bars (tqdm if available via metainformant.core.progress).
- Best-effort resume for HTTP(S) via Range requests.
- Retry with exponential backoff.
- Pluggable protocol handlers (http/https, ftp, file, and optional external tools).
"""

from __future__ import annotations

import json
import os
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Protocol
from urllib.parse import urlparse

from metainformant.core.utils import logging
from metainformant.core.utils import progress

logger = logging.get_logger(__name__)

#
# NOTE: This module intentionally keeps protocol handling lightweight and dependency-minimal.
# For HTTP(S) we use requests (already used in core.io and elsewhere).
# For FTP/file we use urllib/pathlib to avoid extra dependencies.
#


@dataclass(frozen=True)
class DownloadResult:
    """Result of a download attempt."""

    url: str
    dest_path: Path
    success: bool
    bytes_downloaded: int = 0
    total_bytes: int | None = None
    elapsed_seconds: float = 0.0
    resumed: bool = False
    error: str | None = None
    attempts: int = 1


@dataclass
class DownloadHeartbeatState:
    """Mutable state written to heartbeat file periodically."""

    url: str
    destination: str
    started_at: str
    last_update: str
    bytes_downloaded: int
    total_bytes: int | None
    progress_percent: float | None
    speed_mbps: float | None
    eta_seconds: float | None
    status: str
    step: str | None = None
    progress: dict[str, Any] | None = None
    errors: list[str] = field(default_factory=list)

    def to_json(self) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "url": self.url,
            "destination": self.destination,
            "started_at": self.started_at,
            "last_update": self.last_update,
            "bytes_downloaded": self.bytes_downloaded,
            "total_bytes": self.total_bytes,
            "progress_percent": self.progress_percent,
            "speed_mbps": self.speed_mbps,
            "eta_seconds": self.eta_seconds,
            "status": self.status,
            "errors": list(self.errors),
        }
        if self.step is not None:
            payload["step"] = self.step
        if self.progress is not None:
            payload["progress"] = self.progress
        return payload


class ProgressCallback(Protocol):
    def __call__(self, *, bytes_written: int, total_bytes: int | None) -> None: ...


def _utc_iso() -> str:
    # Avoid datetime import to keep this lightweight; ISO-like is fine for heartbeat.
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def _heartbeat_path(dest_path: Path) -> Path:
    # Keep heartbeat files near outputs but segregated.
    hb_dir = dest_path.parent / ".downloads"
    hb_dir.mkdir(parents=True, exist_ok=True)
    return hb_dir / f"{dest_path.name}.heartbeat.json"


def _atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.write("\n")
    os.replace(tmp, path)


def _compute_speed_mbps(bytes_downloaded: int, elapsed_seconds: float) -> float | None:
    if elapsed_seconds <= 0:
        return None
    return (bytes_downloaded / (1024 * 1024)) / elapsed_seconds


def _compute_eta_seconds(bytes_downloaded: int, total_bytes: int | None, elapsed_seconds: float) -> float | None:
    if total_bytes is None or total_bytes <= 0:
        return None
    if bytes_downloaded <= 0:
        return None
    speed_bps = bytes_downloaded / max(elapsed_seconds, 1e-9)
    remaining = max(0, total_bytes - bytes_downloaded)
    return remaining / max(speed_bps, 1e-9)


def _compute_progress_percent(bytes_downloaded: int, total_bytes: int | None) -> float | None:
    if total_bytes is None or total_bytes <= 0:
        return None
    return min(100.0, (bytes_downloaded / total_bytes) * 100.0)


def _infer_dest_path(url: str, dest_path: str | Path) -> Path:
    dest = Path(dest_path)
    if dest.is_dir():
        name = Path(urlparse(url).path).name or "downloaded_file"
        return dest / name
    return dest


def _download_file_url(src_url: str, dest: Path, *, chunk_size: int = 1024 * 1024) -> DownloadResult:
    """Download a file:// URL by copying in chunks (supports progress/heartbeat in tests)."""
    parsed = urlparse(src_url)
    src = Path(parsed.path)
    if not src.exists():
        return DownloadResult(url=src_url, dest_path=dest, success=False, error=f"Source not found: {src}")
    total = src.stat().st_size
    dest.parent.mkdir(parents=True, exist_ok=True)
    start = time.time()
    written = 0
    with open(src, "rb") as rfh, open(dest, "wb") as wfh:
        while True:
            chunk = rfh.read(chunk_size)
            if not chunk:
                break
            wfh.write(chunk)
            written += len(chunk)
    return DownloadResult(
        url=src_url,
        dest_path=dest,
        success=True,
        bytes_downloaded=written,
        total_bytes=total,
        elapsed_seconds=time.time() - start,
        resumed=False,
    )


def _http_head_size(url: str, *, timeout: int) -> int | None:
    import requests

    try:
        resp = requests.head(url, allow_redirects=True, timeout=timeout)
        if resp.status_code >= 400:
            return None
        cl = resp.headers.get("Content-Length")
        return int(cl) if cl is not None else None
    except (requests.RequestException, requests.Timeout, ValueError, OSError) as e:
        logger.debug(f"Could not determine file size for {url}: {e}")
        return None


def _download_http(
    url: str,
    dest: Path,
    *,
    chunk_size: int,
    timeout: int,
    resume: bool,
    progress_cb: ProgressCallback | None,
) -> DownloadResult:
    import requests

    dest.parent.mkdir(parents=True, exist_ok=True)

    total = _http_head_size(url, timeout=timeout)
    start = time.time()
    bytes_written = 0
    resumed = False

    headers: dict[str, str] = {}
    mode = "wb"
    existing = dest.stat().st_size if dest.exists() else 0
    if resume and existing > 0:
        headers["Range"] = f"bytes={existing}-"
        mode = "ab"
        resumed = True

    with requests.get(url, stream=True, timeout=timeout, headers=headers, allow_redirects=True) as resp:
        # If Range not supported, server may return 200 and ignore range; restart.
        if resumed and resp.status_code == 200:
            mode = "wb"
            resumed = False
            existing = 0
        resp.raise_for_status()

        # If we resumed and server returned 206, total refers to full file size; ok.
        with open(dest, mode) as fh:
            bytes_written = existing
            for chunk in resp.iter_content(chunk_size=chunk_size):
                if not chunk:
                    continue
                fh.write(chunk)
                bytes_written += len(chunk)
                if progress_cb:
                    progress_cb(bytes_written=bytes_written, total_bytes=total)

    return DownloadResult(
        url=url,
        dest_path=dest,
        success=True,
        bytes_downloaded=bytes_written,
        total_bytes=total,
        elapsed_seconds=time.time() - start,
        resumed=resumed,
    )


def _download_ftp(
    url: str,
    dest: Path,
    *,
    chunk_size: int,
    timeout: int,
    progress_cb: ProgressCallback | None,
) -> DownloadResult:
    # urllib can handle ftp URLs; size may be unknown.
    import urllib.request

    dest.parent.mkdir(parents=True, exist_ok=True)
    start = time.time()
    bytes_written = 0

    with urllib.request.urlopen(url, timeout=timeout) as resp, open(dest, "wb") as fh:
        while True:
            chunk = resp.read(chunk_size)
            if not chunk:
                break
            fh.write(chunk)
            bytes_written += len(chunk)
            if progress_cb:
                progress_cb(bytes_written=bytes_written, total_bytes=None)

    return DownloadResult(
        url=url,
        dest_path=dest,
        success=True,
        bytes_downloaded=bytes_written,
        total_bytes=None,
        elapsed_seconds=time.time() - start,
        resumed=False,
    )


def download_with_progress(
    url: str,
    dest_path: str | Path,
    *,
    protocol: str = "auto",
    heartbeat_interval: int = 5,
    show_progress: bool = True,
    resume: bool = True,
    max_retries: int = 3,
    retry_delay: float = 5.0,
    chunk_size: int = 1024 * 1024,
    timeout: int = 300,
) -> DownloadResult:
    """Download a file with optional progress + heartbeat.

    Heartbeat file is written under `dest.parent/.downloads/<filename>.heartbeat.json`.
    """
    dest = _infer_dest_path(url, dest_path)
    hb_path = _heartbeat_path(dest)

    parsed = urlparse(url)
    scheme = parsed.scheme.lower()
    if protocol != "auto":
        scheme = protocol.lower()

    started = _utc_iso()
    last_hb = 0.0
    errors: list[str] = []

    # Setup progress bar (bytes-based) if requested.
    # For http we can discover total via HEAD. For others may be None.
    total_for_bar: int | None = None
    if show_progress and scheme in ("http", "https"):
        total_for_bar = _http_head_size(url, timeout=timeout)

    bar = None
    if show_progress:
        desc = f"Downloading {dest.name}"
        bar = progress.progress_bar(total=total_for_bar, desc=desc, unit="B", unit_scale=True)

    # Initialise start_time BEFORE defining _progress_cb so the closure
    # captures an already-initialised variable (previously it was set after
    # the closure definition, which worked by accident of late binding but
    # was fragile).
    start_time = time.time()

    def _progress_cb(*, bytes_written: int, total_bytes: int | None) -> None:
        nonlocal last_hb
        now = time.time()
        if heartbeat_interval > 0 and (now - last_hb) >= heartbeat_interval:
            elapsed = max(0.0, now - start_time)
            pct = _compute_progress_percent(bytes_written, total_bytes)
            speed = _compute_speed_mbps(bytes_written, elapsed)
            eta = _compute_eta_seconds(bytes_written, total_bytes, elapsed)
            state = DownloadHeartbeatState(
                url=url,
                destination=str(dest),
                started_at=started,
                last_update=_utc_iso(),
                bytes_downloaded=bytes_written,
                total_bytes=total_bytes,
                progress_percent=pct,
                speed_mbps=speed,
                eta_seconds=eta,
                status="downloading",
                errors=list(errors),
            )
            _atomic_write_json(hb_path, state.to_json())
            last_hb = now

        if bar is not None:
            # tqdm expects delta updates
            delta = bytes_written - bar.n
            if delta > 0:
                bar.update(delta)

    attempt = 0
    final: DownloadResult | None = None
    try:
        while attempt < max_retries:
            attempt += 1
            try:
                # Mark started/attempt
                state = DownloadHeartbeatState(
                    url=url,
                    destination=str(dest),
                    started_at=started,
                    last_update=_utc_iso(),
                    bytes_downloaded=dest.stat().st_size if dest.exists() else 0,
                    total_bytes=total_for_bar,
                    progress_percent=_compute_progress_percent(
                        dest.stat().st_size if dest.exists() else 0, total_for_bar
                    ),
                    speed_mbps=None,
                    eta_seconds=None,
                    status="starting",
                    errors=list(errors),
                )
                _atomic_write_json(hb_path, state.to_json())

                if scheme in ("http", "https"):
                    final = _download_http(
                        url,
                        dest,
                        chunk_size=chunk_size,
                        timeout=timeout,
                        resume=resume,
                        progress_cb=_progress_cb,
                    )
                elif scheme == "ftp":
                    final = _download_ftp(url, dest, chunk_size=chunk_size, timeout=timeout, progress_cb=_progress_cb)
                elif scheme == "file":
                    final = _download_file_url(url, dest, chunk_size=chunk_size)
                else:
                    raise ValueError(f"Unsupported protocol/scheme: {scheme} (url={url})")

                # Success heartbeat
                elapsed = final.elapsed_seconds
                pct = _compute_progress_percent(final.bytes_downloaded, final.total_bytes)
                speed = _compute_speed_mbps(final.bytes_downloaded, elapsed)
                state = DownloadHeartbeatState(
                    url=url,
                    destination=str(dest),
                    started_at=started,
                    last_update=_utc_iso(),
                    bytes_downloaded=final.bytes_downloaded,
                    total_bytes=final.total_bytes,
                    progress_percent=pct,
                    speed_mbps=speed,
                    eta_seconds=0.0,
                    status="completed" if final.success else "failed",
                    errors=list(errors),
                )
                _atomic_write_json(hb_path, state.to_json())
                return DownloadResult(**{**final.__dict__, "attempts": attempt})
            except Exception as e:
                msg = f"Attempt {attempt}/{max_retries} failed: {e}"
                logger.warning(msg)
                errors.append(str(e))
                state = DownloadHeartbeatState(
                    url=url,
                    destination=str(dest),
                    started_at=started,
                    last_update=_utc_iso(),
                    bytes_downloaded=dest.stat().st_size if dest.exists() else 0,
                    total_bytes=total_for_bar,
                    progress_percent=_compute_progress_percent(
                        dest.stat().st_size if dest.exists() else 0, total_for_bar
                    ),
                    speed_mbps=None,
                    eta_seconds=None,
                    status="failed",
                    errors=list(errors),
                )
                _atomic_write_json(hb_path, state.to_json())
                if attempt >= max_retries:
                    break
                time.sleep(retry_delay * (2 ** (attempt - 1)))
    finally:
        if bar is not None:
            try:
                bar.close()
            except Exception:
                pass  # Ignore errors when closing

    elapsed = time.time() - start_time
    return DownloadResult(
        url=url,
        dest_path=dest,
        success=False,
        bytes_downloaded=dest.stat().st_size if dest.exists() else 0,
        total_bytes=total_for_bar,
        elapsed_seconds=elapsed,
        resumed=False,
        error=errors[-1] if errors else "Download failed",
        attempts=attempt,
    )


def _directory_size_bytes(root: Path, *, glob: str = "**/*") -> int:
    if not root.exists():
        return 0
    total = 0
    for p in root.rglob(glob):
        try:
            if p.is_file():
                total += p.stat().st_size
        except OSError:
            continue
    return total


class _CachedDirectorySize:
    """TTL-based cache for directory size calculations.

    Avoids re-walking the entire directory tree on every call when the
    monitoring loops already sleep(1) between iterations.  A 1-second
    TTL means at most one full walk per second per unique (root, glob) key.
    """

    def __init__(self, ttl: float = 1.0) -> None:
        self._cache: dict[str, tuple[float, int]] = {}  # key -> (timestamp, size)
        self._ttl = ttl

    def get(self, root: Path, glob: str = "**/*") -> int:
        """Return directory size, served from cache if within TTL."""
        key = f"{root}:{glob}"
        now = time.time()
        entry = self._cache.get(key)
        if entry is not None:
            ts, size = entry
            if now - ts < self._ttl:
                return size
        size = _directory_size_bytes(root, glob=glob)
        self._cache[key] = (now, size)
        return size

    def clear(self) -> None:
        """Evict all cached entries."""
        self._cache.clear()


def monitor_subprocess_directory_growth(
    *,
    process: "Any",
    watch_dir: str | Path,
    heartbeat_path: str | Path,
    total_bytes: int | None = None,
    heartbeat_interval: int = 5,
    show_progress: bool = True,
    desc: str | None = None,
    errors: list[str] | None = None,
) -> tuple[int, int]:
    """Monitor a subprocess by watching bytes written under a directory.

    This is intended for external tools (e.g., `amalgkit getfastq`) where the tool handles
    downloads internally, but we still want heartbeat + progress.

    Notes:
    - Assumes the process stdout/stderr are NOT piped (avoid deadlocks).
    - Updates heartbeat JSON at `heartbeat_path` every `heartbeat_interval` seconds.

    Returns:
        (returncode, bytes_written)
    """
    watch = Path(watch_dir)
    hb = Path(heartbeat_path)
    hb.parent.mkdir(parents=True, exist_ok=True)
    started = _utc_iso()
    errs = errors or []

    bar = None
    if show_progress:
        bar = progress.progress_bar(
            total=total_bytes, desc=desc or f"Running process in {watch.name}", unit="B", unit_scale=True
        )

    last = 0.0
    last_bytes = 0
    start_time = time.time()
    _size_cache = _CachedDirectorySize(ttl=1.0)

    def _write(status: str, bytes_done: int) -> None:
        elapsed = max(0.0, time.time() - start_time)
        state = DownloadHeartbeatState(
            url=str(getattr(process, "args", "subprocess")),
            destination=str(watch),
            started_at=started,
            last_update=_utc_iso(),
            bytes_downloaded=bytes_done,
            total_bytes=total_bytes,
            progress_percent=_compute_progress_percent(bytes_done, total_bytes),
            speed_mbps=_compute_speed_mbps(bytes_done, elapsed),
            eta_seconds=_compute_eta_seconds(bytes_done, total_bytes, elapsed),
            status=status,
            step=desc,
            progress=(
                None
                if total_bytes is None
                else {
                    "type": "directory_size",
                    "current": bytes_done,
                    "total": total_bytes,
                    "percent": _compute_progress_percent(bytes_done, total_bytes),
                }
            ),
            errors=list(errs),
        )
        _atomic_write_json(hb, state.to_json())

    while True:
        rc = process.poll()
        now = time.time()
        bytes_done = _size_cache.get(watch)

        if bar is not None:
            delta = bytes_done - bar.n
            if delta > 0:
                bar.update(delta)

        if heartbeat_interval > 0 and (now - last) >= heartbeat_interval:
            _write("running" if rc is None else ("completed" if rc == 0 else "failed"), bytes_done)
            last = now

        if rc is not None:
            # Final write
            _write("completed" if rc == 0 else "failed", bytes_done)
            last_bytes = bytes_done
            break

        time.sleep(1)

    if bar is not None:
        try:
            bar.close()
        except Exception:
            pass  # Ignore errors when closing

    return int(rc), int(last_bytes)


def monitor_subprocess_file_count(
    *,
    process: "Any",
    watch_dir: str | Path,
    heartbeat_path: str | Path,
    expected_files: Iterable[str],
    heartbeat_interval: int = 5,
    show_progress: bool = True,
    desc: str | None = None,
    errors: list[str] | None = None,
) -> int:
    """Monitor a subprocess by counting expected files under a directory."""
    watch = Path(watch_dir)
    hb = Path(heartbeat_path)
    hb.parent.mkdir(parents=True, exist_ok=True)
    started = _utc_iso()
    errs = errors or []
    expected = list(expected_files)
    total = len(expected)

    bar = None
    if show_progress:
        bar = progress.progress_bar(total=total, desc=desc or "Files", unit="file")

    last = 0.0
    start_time = time.time()
    _size_cache = _CachedDirectorySize(ttl=1.0)

    def _count_done() -> int:
        done = 0
        for rel in expected:
            if (watch / rel).exists():
                done += 1
        return done

    def _write(status: str, done: int) -> None:
        elapsed = max(0.0, time.time() - start_time)
        pct = (done / total) * 100.0 if total else None
        cached_size = _size_cache.get(watch)
        state = DownloadHeartbeatState(
            url=str(getattr(process, "args", "subprocess")),
            destination=str(watch),
            started_at=started,
            last_update=_utc_iso(),
            bytes_downloaded=cached_size,
            total_bytes=None,
            progress_percent=pct,
            speed_mbps=_compute_speed_mbps(cached_size, elapsed),
            eta_seconds=None,
            status=status,
            step=desc,
            progress={"type": "file_count", "current": done, "total": total, "percent": pct},
            errors=list(errs),
        )
        _atomic_write_json(hb, state.to_json())

    while True:
        rc = process.poll()
        done = _count_done()
        if bar is not None:
            try:
                bar.n = done
                bar.refresh()
            except (AttributeError, TypeError) as e:
                # Handle tqdm compatibility issues gracefully
                # If refresh fails, just update n and continue
                if "refresh" in str(e).lower():
                    try:
                        bar.n = done
                    except Exception:
                        pass  # Ignore if we can't update
        now = time.time()
        if heartbeat_interval > 0 and (now - last) >= heartbeat_interval:
            _write("running" if rc is None else ("completed" if rc == 0 else "failed"), done)
            last = now
        if rc is not None:
            _write("completed" if rc == 0 else "failed", done)
            break
        time.sleep(1)

    if bar is not None:
        try:
            bar.close()
        except Exception:
            pass  # Ignore errors when closing
    return int(rc)


def monitor_subprocess_sample_progress(
    *,
    process: "Any",
    watch_dir: str | Path,
    heartbeat_path: str | Path,
    completion_glob: str,
    total_samples: int | None = None,
    heartbeat_interval: int = 5,
    show_progress: bool = True,
    desc: str | None = None,
    errors: list[str] | None = None,
) -> int:
    """Monitor a subprocess by counting completed samples (files matching a glob)."""
    watch = Path(watch_dir)
    hb = Path(heartbeat_path)
    hb.parent.mkdir(parents=True, exist_ok=True)
    started = _utc_iso()
    errs = errors or []

    bar = None
    if show_progress:
        bar = progress.progress_bar(total=total_samples, desc=desc or "Samples", unit="sample")

    last = 0.0
    start_time = time.time()
    _size_cache = _CachedDirectorySize(ttl=1.0)

    def _count_done() -> int:
        return len(list(watch.glob(completion_glob)))

    def _write(status: str, done: int) -> None:
        elapsed = max(0.0, time.time() - start_time)
        pct = None
        if total_samples:
            pct = min(100.0, (done / total_samples) * 100.0)
        cached_size = _size_cache.get(watch)
        state = DownloadHeartbeatState(
            url=str(getattr(process, "args", "subprocess")),
            destination=str(watch),
            started_at=started,
            last_update=_utc_iso(),
            bytes_downloaded=cached_size,
            total_bytes=None,
            progress_percent=pct,
            speed_mbps=_compute_speed_mbps(cached_size, elapsed),
            eta_seconds=None,
            status=status,
            step=desc,
            progress={"type": "sample_count", "current": done, "total": total_samples, "percent": pct},
            errors=list(errs),
        )
        _atomic_write_json(hb, state.to_json())

    while True:
        rc = process.poll()
        done = _count_done()
        if bar is not None:
            try:
                bar.n = done
                bar.refresh()
            except (AttributeError, TypeError) as e:
                # Handle tqdm compatibility issues gracefully
                # If refresh fails, just update n and continue
                if "refresh" in str(e).lower():
                    try:
                        bar.n = done
                    except Exception:
                        pass  # Ignore if we can't update
        now = time.time()
        if heartbeat_interval > 0 and (now - last) >= heartbeat_interval:
            _write("running" if rc is None else ("completed" if rc == 0 else "failed"), done)
            last = now
        if rc is not None:
            _write("completed" if rc == 0 else "failed", done)
            break
        time.sleep(1)

    if bar is not None:
        try:
            bar.close()
        except Exception:
            pass  # Ignore errors when closing
    return int(rc)


# ---- Optional OO-style handler interface (for modularity/extensibility) ----


class DownloadHandler(Protocol):
    """Protocol for pluggable download handlers."""

    def get_file_size(self, url: str, *, timeout: int = 60) -> int | None: ...

    def download(
        self,
        url: str,
        dest_path: str | Path,
        *,
        heartbeat_interval: int = 5,
        show_progress: bool = True,
        resume: bool = True,
        max_retries: int = 3,
        retry_delay: float = 5.0,
        chunk_size: int = 1024 * 1024,
        timeout: int = 300,
    ) -> DownloadResult: ...


class HTTPDownloadHandler:
    def get_file_size(self, url: str, *, timeout: int = 60) -> int | None:
        return _http_head_size(url, timeout=timeout)

    def download(self, url: str, dest_path: str | Path, **kwargs: Any) -> DownloadResult:
        return download_with_progress(url, dest_path, protocol="https" if url.startswith("https") else "http", **kwargs)


class FTPDownloadHandler:
    def get_file_size(self, url: str, *, timeout: int = 60) -> int | None:
        # FTP size is not reliably available without extra control channel work; treat unknown.
        return None

    def download(self, url: str, dest_path: str | Path, **kwargs: Any) -> DownloadResult:
        return download_with_progress(url, dest_path, protocol="ftp", resume=False, **kwargs)


class FileDownloadHandler:
    def get_file_size(self, url: str, *, timeout: int = 60) -> int | None:
        parsed = urlparse(url)
        src = Path(parsed.path)
        return src.stat().st_size if src.exists() else None

    def download(self, url: str, dest_path: str | Path, **kwargs: Any) -> DownloadResult:
        return download_with_progress(url, dest_path, protocol="file", resume=False, **kwargs)


def get_download_handler(url: str) -> DownloadHandler:
    scheme = urlparse(url).scheme.lower()
    if scheme in ("http", "https"):
        return HTTPDownloadHandler()
    if scheme == "ftp":
        return FTPDownloadHandler()
    if scheme == "file":
        return FileDownloadHandler()
    raise ValueError(f"No handler for scheme: {scheme}")
