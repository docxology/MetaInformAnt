"""Atomic file operations for METAINFORMANT."""

from __future__ import annotations

import os
import shutil
import tempfile
from collections.abc import Generator
from contextlib import contextmanager
from pathlib import Path

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


@contextmanager
def atomic_write(path: Path | str, mode: str = "w", encoding: str = "utf-8") -> Generator:
    """Context manager that writes to a temp file, then does atomic rename on success.

    The file is written to a temporary location in the same directory as the
    target path. On successful exit, the temp file is atomically renamed to
    the target path. On exception, the temp file is cleaned up.

    Args:
        path: Destination file path.
        mode: File open mode ("w" for text, "wb" for binary).
        encoding: Text encoding (ignored for binary mode).

    Yields:
        File handle for writing.

    Example:
        with atomic_write("/data/results.json") as f:
            json.dump(data, f)
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    is_binary = "b" in mode
    fd, tmp_path_str = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
    )
    tmp_path = Path(tmp_path_str)

    try:
        with os.fdopen(fd, mode, encoding=None if is_binary else encoding) as f:
            yield f
        # Atomic rename on success
        tmp_path.replace(path)
        logger.debug("Atomic write completed: %s", path)
    except BaseException:
        # Clean up temp file on failure
        try:
            tmp_path.unlink(missing_ok=True)
        except OSError:
            pass
        raise


def atomic_replace(src: Path | str, dst: Path | str) -> None:
    """Atomically replace dst with src.

    Uses os.replace which is atomic on POSIX systems when src and dst
    are on the same filesystem. The source file is removed after replacement.

    Args:
        src: Source file path.
        dst: Destination file path (will be replaced).

    Raises:
        FileNotFoundError: If src does not exist.
    """
    src = Path(src)
    dst = Path(dst)
    if not src.is_file():
        raise FileNotFoundError(f"Source file not found: {src}")
    dst.parent.mkdir(parents=True, exist_ok=True)
    os.replace(src, dst)
    logger.debug("Atomic replace: %s -> %s", src, dst)


def safe_write_text(path: Path | str, content: str, encoding: str = "utf-8") -> None:
    """Write text content to a file atomically.

    Uses atomic_write internally to ensure the file is either fully written
    or not modified at all.

    Args:
        path: Destination file path.
        content: Text content to write.
        encoding: Text encoding. Default "utf-8".
    """
    with atomic_write(path, mode="w", encoding=encoding) as f:
        f.write(content)


def safe_write_bytes(path: Path | str, content: bytes) -> None:
    """Write binary content to a file atomically.

    Uses atomic_write internally to ensure the file is either fully written
    or not modified at all.

    Args:
        path: Destination file path.
        content: Binary content to write.
    """
    with atomic_write(path, mode="wb") as f:
        f.write(content)


@contextmanager
def temp_directory(prefix: str = "metainformant_", cleanup: bool = True) -> Generator[Path, None, None]:
    """Context manager yielding a temporary directory.

    Creates a temporary directory that is optionally cleaned up when the
    context exits. Uses the system temp directory as the parent.

    Args:
        prefix: Prefix for the temporary directory name.
        cleanup: Whether to remove the directory on exit. Default True.

    Yields:
        Path to the temporary directory.

    Example:
        with temp_directory() as tmp:
            work_file = tmp / "intermediate.dat"
            work_file.write_text("data")
    """
    tmp_dir = Path(tempfile.mkdtemp(prefix=prefix))
    logger.debug("Created temp directory: %s", tmp_dir)
    try:
        yield tmp_dir
    finally:
        if cleanup:
            try:
                shutil.rmtree(tmp_dir, ignore_errors=True)
                logger.debug("Cleaned up temp directory: %s", tmp_dir)
            except OSError:
                logger.warning("Failed to clean up temp directory: %s", tmp_dir)
        else:
            logger.debug("Temp directory retained: %s", tmp_dir)
