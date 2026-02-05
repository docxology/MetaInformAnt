"""Comprehensive tests for the mtime-based discovery cache and download module improvements.

Tests cover:
- _DiscoveryCache: hit, invalidation, clear, thread safety
- discover_functions / find_symbol_usage / build_call_graph integration with cache
- invalidate_discovery_cache public API
- _CachedDirectorySize: TTL behaviour, separate keys, empty dirs
- download_with_progress with file:// URL
- DownloadResult and DownloadHeartbeatState serialisation

All tests use real implementations (NO mocking policy).
"""

from __future__ import annotations

import ast
import json
import os
import textwrap
import threading
import time
from pathlib import Path

import pytest

from metainformant.core.execution.discovery import (
    _DiscoveryCache,
    build_call_graph,
    discover_functions,
    find_symbol_usage,
    get_module_dependencies,
    invalidate_discovery_cache,
)
from metainformant.core.io.download import (
    DownloadHeartbeatState,
    DownloadResult,
    _CachedDirectorySize,
    _compute_eta_seconds,
    _compute_progress_percent,
    _compute_speed_mbps,
    download_with_progress,
)

# ---------------------------------------------------------------------------
# Helper: create a temporary Python module file
# ---------------------------------------------------------------------------


def _write_module(path: Path, source: str) -> Path:
    """Write a Python source string to *path* and return the path."""
    path.write_text(textwrap.dedent(source), encoding="utf-8")
    return path


# ===========================================================================
# _DiscoveryCache unit tests
# ===========================================================================


class TestDiscoveryCacheHit:
    """Cache hit: parsing the same unmodified file twice returns a cached AST."""

    def test_cache_hit_returns_same_tree(self, tmp_path: Path) -> None:
        src = _write_module(
            tmp_path / "mod.py",
            """\
            def greet(name: str) -> str:
                return f"hello {name}"
            """,
        )

        cache = _DiscoveryCache()
        tree1 = ast.parse(src.read_text(encoding="utf-8"), filename=str(src))
        cache.put(src, tree1)

        tree2 = cache.get(src)
        assert tree2 is tree1, "Second get() should return the exact cached object"

    def test_cache_hit_avoids_reparse(self, tmp_path: Path) -> None:
        src = _write_module(
            tmp_path / "mod.py",
            """\
            x = 1
            """,
        )

        cache = _DiscoveryCache()
        tree_a = ast.parse(src.read_text(encoding="utf-8"), filename=str(src))
        cache.put(src, tree_a)

        # Read twice; both should be the same cached object.
        first = cache.get(src)
        second = cache.get(src)
        assert first is second


class TestDiscoveryCacheInvalidationOnModification:
    """Cache invalidation when a file is modified on disk."""

    def test_stale_entry_evicted_after_file_change(self, tmp_path: Path) -> None:
        src = _write_module(tmp_path / "mod.py", "x = 1\n")

        cache = _DiscoveryCache()
        tree1 = ast.parse(src.read_text(encoding="utf-8"), filename=str(src))
        cache.put(src, tree1)

        # Ensure the filesystem mtime actually advances (some filesystems have
        # 1-second granularity).
        time.sleep(0.05)
        _write_module(src, "x = 2\n")
        # Force mtime to be strictly greater than what was cached.
        new_mtime = src.stat().st_mtime + 1.0
        os.utime(src, (new_mtime, new_mtime))

        result = cache.get(src)
        assert result is None, "Stale entry should be evicted after file modification"

    def test_reparse_after_modification_gives_new_tree(self, tmp_path: Path) -> None:
        src = _write_module(tmp_path / "mod.py", "x = 1\n")

        cache = _DiscoveryCache()
        tree1 = ast.parse(src.read_text(encoding="utf-8"), filename=str(src))
        cache.put(src, tree1)

        # Mutate file with a guaranteed mtime bump.
        time.sleep(0.05)
        _write_module(src, "x = 2\ny = 3\n")
        new_mtime = src.stat().st_mtime + 1.0
        os.utime(src, (new_mtime, new_mtime))

        assert cache.get(src) is None

        tree2 = ast.parse(src.read_text(encoding="utf-8"), filename=str(src))
        cache.put(src, tree2)

        assert cache.get(src) is tree2
        assert tree2 is not tree1


class TestDiscoveryCacheInvalidateMethod:
    """Explicit invalidate() removes a single entry."""

    def test_invalidate_removes_entry(self, tmp_path: Path) -> None:
        src = _write_module(tmp_path / "mod.py", "a = 1\n")

        cache = _DiscoveryCache()
        tree = ast.parse(src.read_text(encoding="utf-8"), filename=str(src))
        cache.put(src, tree)
        assert cache.get(src) is tree

        cache.invalidate(src)
        assert cache.get(src) is None

    def test_invalidate_nonexistent_key_is_noop(self, tmp_path: Path) -> None:
        cache = _DiscoveryCache()
        # Should not raise.
        cache.invalidate(tmp_path / "nonexistent.py")


class TestDiscoveryCacheClearMethod:
    """clear() removes all entries."""

    def test_clear_empties_cache(self, tmp_path: Path) -> None:
        a = _write_module(tmp_path / "a.py", "a = 1\n")
        b = _write_module(tmp_path / "b.py", "b = 2\n")

        cache = _DiscoveryCache()
        for p in (a, b):
            tree = ast.parse(p.read_text(encoding="utf-8"), filename=str(p))
            cache.put(p, tree)

        assert cache.get(a) is not None
        assert cache.get(b) is not None

        cache.clear()

        assert cache.get(a) is None
        assert cache.get(b) is None


class TestDiscoveryCacheDeletedFile:
    """Cache evicts entry when the underlying file has been deleted."""

    def test_deleted_file_evicted(self, tmp_path: Path) -> None:
        src = _write_module(tmp_path / "mod.py", "x = 1\n")

        cache = _DiscoveryCache()
        tree = ast.parse(src.read_text(encoding="utf-8"), filename=str(src))
        cache.put(src, tree)
        assert cache.get(src) is tree

        src.unlink()
        assert cache.get(src) is None


class TestDiscoveryCacheThreadSafety:
    """Concurrent reads/writes must not raise or corrupt state."""

    def test_concurrent_put_and_get(self, tmp_path: Path) -> None:
        num_files = 20
        files = []
        for i in range(num_files):
            p = _write_module(tmp_path / f"mod_{i}.py", f"x_{i} = {i}\n")
            files.append(p)

        cache = _DiscoveryCache()
        errors: list[str] = []

        def _put(p: Path) -> None:
            try:
                tree = ast.parse(p.read_text(encoding="utf-8"), filename=str(p))
                cache.put(p, tree)
            except Exception as exc:
                errors.append(f"put({p}): {exc}")

        def _get(p: Path) -> None:
            try:
                cache.get(p)  # May be None or a tree; either is fine.
            except Exception as exc:
                errors.append(f"get({p}): {exc}")

        threads: list[threading.Thread] = []
        for p in files:
            threads.append(threading.Thread(target=_put, args=(p,)))
            threads.append(threading.Thread(target=_get, args=(p,)))

        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=10)

        assert errors == [], f"Thread errors: {errors}"

    def test_concurrent_invalidate_and_get(self, tmp_path: Path) -> None:
        src = _write_module(tmp_path / "mod.py", "val = 42\n")

        cache = _DiscoveryCache()
        tree = ast.parse(src.read_text(encoding="utf-8"), filename=str(src))
        cache.put(src, tree)

        errors: list[str] = []
        iterations = 200

        def _invalidator() -> None:
            try:
                for _ in range(iterations):
                    cache.invalidate(src)
                    cache.put(src, tree)
            except Exception as exc:
                errors.append(f"invalidator: {exc}")

        def _reader() -> None:
            try:
                for _ in range(iterations):
                    cache.get(src)  # None or tree.
            except Exception as exc:
                errors.append(f"reader: {exc}")

        t1 = threading.Thread(target=_invalidator)
        t2 = threading.Thread(target=_reader)
        t1.start()
        t2.start()
        t1.join(timeout=10)
        t2.join(timeout=10)

        assert errors == [], f"Thread errors: {errors}"


# ===========================================================================
# discover_functions integration with cache
# ===========================================================================


class TestDiscoverFunctionsWithCache:
    """Ensure discover_functions uses and respects the cache."""

    def test_discover_functions_basic(self, tmp_path: Path) -> None:
        src = _write_module(
            tmp_path / "funcs.py",
            """\
            def alpha(x: int) -> int:
                \"\"\"Return x doubled.\"\"\"
                return x * 2

            def beta(a: str, b: str = "!") -> str:
                return a + b
            """,
        )

        # Clear the module-level cache so tests are isolated.
        invalidate_discovery_cache()

        funcs = discover_functions(src)
        assert len(funcs) == 2
        names = {f.name for f in funcs}
        assert names == {"alpha", "beta"}

        # Verify fields populated.
        alpha = next(f for f in funcs if f.name == "alpha")
        assert alpha.return_type == "int"
        assert alpha.docstring == "Return x doubled."
        assert "x" in alpha.parameters

    def test_discover_functions_with_pattern_filter(self, tmp_path: Path) -> None:
        src = _write_module(
            tmp_path / "funcs.py",
            """\
            def load_data():
                pass

            def save_data():
                pass

            def load_config():
                pass
            """,
        )

        invalidate_discovery_cache()

        funcs = discover_functions(src, pattern="load")
        assert len(funcs) == 2
        assert all("load" in f.name for f in funcs)

    def test_discover_functions_nonexistent_raises(self) -> None:
        with pytest.raises(FileNotFoundError):
            discover_functions(Path("/absolutely/nonexistent/module.py"))

    def test_discover_functions_cache_hit_on_second_call(self, tmp_path: Path) -> None:
        src = _write_module(
            tmp_path / "cached.py",
            """\
            def one():
                pass
            """,
        )

        invalidate_discovery_cache()

        result1 = discover_functions(src)
        result2 = discover_functions(src)
        assert len(result1) == len(result2) == 1
        assert result1[0].name == result2[0].name == "one"

    def test_discover_functions_invalidation_on_file_change(self, tmp_path: Path) -> None:
        src = _write_module(
            tmp_path / "change.py",
            """\
            def original():
                pass
            """,
        )

        invalidate_discovery_cache()

        funcs1 = discover_functions(src)
        assert [f.name for f in funcs1] == ["original"]

        # Rewrite with different content and bump mtime.
        time.sleep(0.05)
        _write_module(
            src,
            """\
            def replacement():
                pass

            def added():
                pass
            """,
        )
        new_mtime = src.stat().st_mtime + 1.0
        os.utime(src, (new_mtime, new_mtime))

        funcs2 = discover_functions(src)
        names2 = {f.name for f in funcs2}
        assert "original" not in names2
        assert names2 == {"replacement", "added"}


# ===========================================================================
# find_symbol_usage integration with cache
# ===========================================================================


class TestFindSymbolUsageWithCache:
    """find_symbol_usage should use _cached_parse internally."""

    def test_find_symbol_in_tmp_repo(self, tmp_path: Path) -> None:
        _write_module(
            tmp_path / "a.py",
            """\
            import os
            my_var = 10
            print(my_var)
            """,
        )
        _write_module(
            tmp_path / "b.py",
            """\
            my_var = 20
            """,
        )

        invalidate_discovery_cache()

        usages = find_symbol_usage("my_var", tmp_path)
        assert len(usages) >= 2  # at least one in a.py, one in b.py
        files_found = {u.file.name for u in usages}
        assert "a.py" in files_found
        assert "b.py" in files_found

    def test_find_symbol_no_match(self, tmp_path: Path) -> None:
        _write_module(tmp_path / "empty.py", "x = 1\n")

        invalidate_discovery_cache()

        usages = find_symbol_usage("absolutely_nonexistent_symbol_xyz", tmp_path)
        assert usages == []


# ===========================================================================
# build_call_graph integration with cache
# ===========================================================================


class TestBuildCallGraphWithCache:
    """build_call_graph should use _cached_parse."""

    def test_call_graph_basic(self, tmp_path: Path) -> None:
        src = _write_module(
            tmp_path / "calls.py",
            """\
            def helper():
                pass

            def main():
                helper()
                print("done")
            """,
        )

        invalidate_discovery_cache()

        graph = build_call_graph(src)
        assert "main" in graph
        assert "helper" in graph["main"]
        assert "print" in graph["main"]

    def test_call_graph_nonexistent_returns_empty(self) -> None:
        graph = build_call_graph(Path("/nonexistent/path.py"))
        assert graph == {}


# ===========================================================================
# get_module_dependencies integration with cache
# ===========================================================================


class TestGetModuleDependenciesWithCache:
    """get_module_dependencies should use _cached_parse."""

    def test_module_deps_basic(self, tmp_path: Path) -> None:
        src = _write_module(
            tmp_path / "deps.py",
            """\
            import os
            import sys
            from pathlib import Path
            from collections import defaultdict, OrderedDict
            """,
        )

        invalidate_discovery_cache()

        deps = get_module_dependencies(src)
        assert "os" in deps.imports
        assert "sys" in deps.imports
        assert "pathlib" in deps.from_imports
        assert "Path" in deps.from_imports["pathlib"]
        assert "collections" in deps.from_imports
        assert set(deps.from_imports["collections"]) == {"defaultdict", "OrderedDict"}

    def test_module_deps_nonexistent_returns_empty(self) -> None:
        deps = get_module_dependencies(Path("/nonexistent/module.py"))
        assert deps.imports == []
        assert deps.from_imports == {}


# ===========================================================================
# invalidate_discovery_cache public API
# ===========================================================================


class TestInvalidateDiscoveryCachePublicAPI:
    """Public helper invalidate_discovery_cache()."""

    def test_invalidate_specific_file(self, tmp_path: Path) -> None:
        a = _write_module(tmp_path / "a.py", "a = 1\n")
        b = _write_module(tmp_path / "b.py", "b = 2\n")

        invalidate_discovery_cache()

        # Populate cache via discover_functions.
        discover_functions(a)
        discover_functions(b)

        # Invalidate only 'a'.
        invalidate_discovery_cache(a)

        # After invalidation, re-discovering 'a' still works.
        funcs_a = discover_functions(a)
        assert len(funcs_a) == 0  # No function defs in "a = 1".

        # 'b' should still be fine (cache or re-parse).
        funcs_b = discover_functions(b)
        assert len(funcs_b) == 0

    def test_invalidate_all(self, tmp_path: Path) -> None:
        a = _write_module(tmp_path / "a.py", "a = 1\n")

        invalidate_discovery_cache()

        discover_functions(a)

        # Clear everything.
        invalidate_discovery_cache(None)

        # Should still work after full clear.
        funcs = discover_functions(a)
        assert isinstance(funcs, list)

    def test_invalidate_with_string_path(self, tmp_path: Path) -> None:
        src = _write_module(tmp_path / "s.py", "s = 1\n")

        invalidate_discovery_cache()
        discover_functions(src)

        # Pass as string.
        invalidate_discovery_cache(str(src))

        # No errors.
        funcs = discover_functions(src)
        assert isinstance(funcs, list)


# ===========================================================================
# _CachedDirectorySize tests
# ===========================================================================


class TestCachedDirectorySize:
    """TTL-based directory size cache."""

    def test_returns_correct_size(self, tmp_path: Path) -> None:
        (tmp_path / "a.bin").write_bytes(b"x" * 100)
        (tmp_path / "b.bin").write_bytes(b"y" * 200)

        cache = _CachedDirectorySize(ttl=10.0)
        size = cache.get(tmp_path)
        assert size == 300

    def test_cache_expires_after_ttl(self, tmp_path: Path) -> None:
        (tmp_path / "a.bin").write_bytes(b"x" * 100)

        cache = _CachedDirectorySize(ttl=0.1)  # 100ms TTL
        size1 = cache.get(tmp_path)
        assert size1 == 100

        # Add more data.
        (tmp_path / "b.bin").write_bytes(b"y" * 200)

        # Within TTL, stale value returned.
        size_stale = cache.get(tmp_path)
        assert size_stale == 100  # Still cached.

        # Wait for TTL to expire.
        time.sleep(0.15)

        size2 = cache.get(tmp_path)
        assert size2 == 300  # Now refreshed.

    def test_different_directories_separate_keys(self, tmp_path: Path) -> None:
        dir_a = tmp_path / "dir_a"
        dir_b = tmp_path / "dir_b"
        dir_a.mkdir()
        dir_b.mkdir()
        (dir_a / "file.bin").write_bytes(b"a" * 50)
        (dir_b / "file.bin").write_bytes(b"b" * 150)

        cache = _CachedDirectorySize(ttl=10.0)
        assert cache.get(dir_a) == 50
        assert cache.get(dir_b) == 150

    def test_empty_directory_returns_zero(self, tmp_path: Path) -> None:
        empty = tmp_path / "empty"
        empty.mkdir()

        cache = _CachedDirectorySize(ttl=10.0)
        assert cache.get(empty) == 0

    def test_nonexistent_directory_returns_zero(self, tmp_path: Path) -> None:
        cache = _CachedDirectorySize(ttl=10.0)
        assert cache.get(tmp_path / "nope") == 0

    def test_clear_evicts_all(self, tmp_path: Path) -> None:
        (tmp_path / "f.bin").write_bytes(b"z" * 77)

        cache = _CachedDirectorySize(ttl=10.0)
        assert cache.get(tmp_path) == 77

        cache.clear()
        # After clear the next call recalculates (same value here, but proves eviction path).
        assert cache.get(tmp_path) == 77

    def test_different_globs_separate_keys(self, tmp_path: Path) -> None:
        (tmp_path / "data.csv").write_bytes(b"c" * 60)
        (tmp_path / "data.txt").write_bytes(b"t" * 40)

        cache = _CachedDirectorySize(ttl=10.0)
        size_csv = cache.get(tmp_path, glob="*.csv")
        size_txt = cache.get(tmp_path, glob="*.txt")
        size_all = cache.get(tmp_path)

        assert size_csv == 60
        assert size_txt == 40
        assert size_all == 100


# ===========================================================================
# download_with_progress with file:// URL
# ===========================================================================


class TestDownloadWithProgressFileScheme:
    """download_with_progress using file:// protocol (no network needed)."""

    def test_basic_file_download(self, tmp_path: Path) -> None:
        src = tmp_path / "source.bin"
        src.write_bytes(b"A" * 4096)

        dest = tmp_path / "dest" / "output.bin"

        result = download_with_progress(
            f"file://{src}",
            dest,
            show_progress=False,
            heartbeat_interval=0,
            max_retries=1,
            chunk_size=1024,
            timeout=5,
            resume=False,
        )

        assert result.success is True
        assert dest.exists()
        assert dest.read_bytes() == src.read_bytes()
        assert result.bytes_downloaded == 4096
        assert result.total_bytes == 4096
        assert result.url == f"file://{src}"
        assert result.dest_path == dest

    def test_heartbeat_file_created(self, tmp_path: Path) -> None:
        src = tmp_path / "source.bin"
        src.write_bytes(b"B" * 2048)

        dest_dir = tmp_path / "dl"
        dest_dir.mkdir()
        dest = dest_dir / "result.bin"

        result = download_with_progress(
            f"file://{src}",
            dest,
            show_progress=False,
            heartbeat_interval=1,
            max_retries=1,
            chunk_size=512,
            timeout=5,
            resume=False,
        )

        assert result.success is True

        hb = dest_dir / ".downloads" / "result.bin.heartbeat.json"
        assert hb.exists()

        payload = json.loads(hb.read_text(encoding="utf-8"))
        assert payload["status"] in ("completed", "starting")
        assert "bytes_downloaded" in payload
        assert "url" in payload

    def test_download_nonexistent_source_fails(self, tmp_path: Path) -> None:
        dest = tmp_path / "out.bin"
        result = download_with_progress(
            f"file://{tmp_path / 'missing.bin'}",
            dest,
            show_progress=False,
            heartbeat_interval=0,
            max_retries=1,
            chunk_size=1024,
            timeout=5,
            resume=False,
        )

        assert result.success is False
        assert result.error is not None
        assert "not found" in result.error.lower() or "Source not found" in result.error

    def test_download_to_directory_infers_filename(self, tmp_path: Path) -> None:
        src = tmp_path / "named_file.dat"
        src.write_bytes(b"C" * 128)

        dest_dir = tmp_path / "target_dir"
        dest_dir.mkdir()

        result = download_with_progress(
            f"file://{src}",
            dest_dir,
            show_progress=False,
            heartbeat_interval=0,
            max_retries=1,
            chunk_size=64,
            timeout=5,
            resume=False,
        )

        assert result.success is True
        expected_dest = dest_dir / "named_file.dat"
        assert expected_dest.exists()
        assert expected_dest.read_bytes() == src.read_bytes()


# ===========================================================================
# DownloadResult dataclass
# ===========================================================================


class TestDownloadResult:
    """DownloadResult frozen dataclass fields."""

    def test_fields_present(self) -> None:
        r = DownloadResult(
            url="file:///test",
            dest_path=Path("/tmp/test"),
            success=True,
            bytes_downloaded=1024,
            total_bytes=1024,
            elapsed_seconds=0.5,
            resumed=False,
            error=None,
            attempts=1,
        )
        assert r.url == "file:///test"
        assert r.dest_path == Path("/tmp/test")
        assert r.success is True
        assert r.bytes_downloaded == 1024
        assert r.total_bytes == 1024
        assert r.elapsed_seconds == 0.5
        assert r.resumed is False
        assert r.error is None
        assert r.attempts == 1

    def test_defaults(self) -> None:
        r = DownloadResult(url="file:///x", dest_path=Path("/x"), success=False)
        assert r.bytes_downloaded == 0
        assert r.total_bytes is None
        assert r.elapsed_seconds == 0.0
        assert r.resumed is False
        assert r.error is None
        assert r.attempts == 1

    def test_frozen(self) -> None:
        r = DownloadResult(url="file:///x", dest_path=Path("/x"), success=True)
        with pytest.raises(AttributeError):
            r.success = False  # type: ignore[misc]


# ===========================================================================
# DownloadHeartbeatState serialisation
# ===========================================================================


class TestDownloadHeartbeatState:
    """DownloadHeartbeatState.to_json() serialisation."""

    def test_to_json_basic(self) -> None:
        state = DownloadHeartbeatState(
            url="file:///src",
            destination="/dest",
            started_at="2025-01-01T00:00:00Z",
            last_update="2025-01-01T00:01:00Z",
            bytes_downloaded=500,
            total_bytes=1000,
            progress_percent=50.0,
            speed_mbps=1.5,
            eta_seconds=30.0,
            status="downloading",
        )

        payload = state.to_json()
        assert payload["url"] == "file:///src"
        assert payload["destination"] == "/dest"
        assert payload["bytes_downloaded"] == 500
        assert payload["total_bytes"] == 1000
        assert payload["progress_percent"] == 50.0
        assert payload["speed_mbps"] == 1.5
        assert payload["eta_seconds"] == 30.0
        assert payload["status"] == "downloading"
        assert payload["errors"] == []
        assert "step" not in payload  # None -> omitted
        assert "progress" not in payload  # None -> omitted

    def test_to_json_with_optional_fields(self) -> None:
        state = DownloadHeartbeatState(
            url="file:///src",
            destination="/dest",
            started_at="2025-01-01T00:00:00Z",
            last_update="2025-01-01T00:01:00Z",
            bytes_downloaded=0,
            total_bytes=None,
            progress_percent=None,
            speed_mbps=None,
            eta_seconds=None,
            status="starting",
            step="getfastq",
            progress={"type": "sample_count", "current": 0, "total": 10},
            errors=["timeout on attempt 1"],
        )

        payload = state.to_json()
        assert payload["step"] == "getfastq"
        assert payload["progress"]["type"] == "sample_count"
        assert payload["errors"] == ["timeout on attempt 1"]
        assert payload["total_bytes"] is None

    def test_to_json_round_trips_through_json(self) -> None:
        state = DownloadHeartbeatState(
            url="https://example.com/file.gz",
            destination="/data/file.gz",
            started_at="2025-06-15T12:00:00Z",
            last_update="2025-06-15T12:05:00Z",
            bytes_downloaded=1048576,
            total_bytes=10485760,
            progress_percent=10.0,
            speed_mbps=3.5,
            eta_seconds=25.7,
            status="downloading",
            errors=[],
        )

        serialised = json.dumps(state.to_json())
        restored = json.loads(serialised)
        assert restored["bytes_downloaded"] == 1048576
        assert restored["status"] == "downloading"
        assert isinstance(restored["errors"], list)


# ===========================================================================
# Helper compute functions
# ===========================================================================


class TestComputeHelpers:
    """Tests for _compute_speed_mbps, _compute_eta_seconds, _compute_progress_percent."""

    def test_speed_mbps_positive(self) -> None:
        # 1 MiB in 1 second = 1.0 MB/s
        assert _compute_speed_mbps(1024 * 1024, 1.0) == pytest.approx(1.0)

    def test_speed_mbps_zero_elapsed(self) -> None:
        assert _compute_speed_mbps(1024, 0.0) is None

    def test_speed_mbps_negative_elapsed(self) -> None:
        assert _compute_speed_mbps(1024, -1.0) is None

    def test_eta_basic(self) -> None:
        # 50% done in 10s => ~10s remaining.
        eta = _compute_eta_seconds(500, 1000, 10.0)
        assert eta is not None
        assert eta == pytest.approx(10.0, rel=0.01)

    def test_eta_unknown_total(self) -> None:
        assert _compute_eta_seconds(500, None, 10.0) is None

    def test_eta_zero_downloaded(self) -> None:
        assert _compute_eta_seconds(0, 1000, 10.0) is None

    def test_progress_percent_half(self) -> None:
        assert _compute_progress_percent(500, 1000) == pytest.approx(50.0)

    def test_progress_percent_full(self) -> None:
        assert _compute_progress_percent(1000, 1000) == pytest.approx(100.0)

    def test_progress_percent_over(self) -> None:
        # Clamped to 100.
        assert _compute_progress_percent(1500, 1000) == pytest.approx(100.0)

    def test_progress_percent_unknown_total(self) -> None:
        assert _compute_progress_percent(500, None) is None

    def test_progress_percent_zero_total(self) -> None:
        assert _compute_progress_percent(0, 0) is None
