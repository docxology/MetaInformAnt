"""File checksum and integrity verification for METAINFORMANT."""

from __future__ import annotations

import hashlib
from pathlib import Path

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

_SUPPORTED_ALGORITHMS = {"md5", "sha256", "sha1", "sha512"}


def _compute_hash(path: Path, algorithm: str, chunk_size: int = 8192) -> str:
    """Compute a hash digest for a file using the specified algorithm.

    Args:
        path: Path to the file to hash.
        algorithm: Hash algorithm name (e.g. "md5", "sha256").
        chunk_size: Read buffer size in bytes.

    Returns:
        Hex digest string.

    Raises:
        FileNotFoundError: If path does not exist.
        ValueError: If algorithm is not supported.
    """
    if algorithm not in _SUPPORTED_ALGORITHMS:
        raise ValueError(f"Unsupported algorithm '{algorithm}'. Supported: {sorted(_SUPPORTED_ALGORITHMS)}")
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"File not found: {path}")

    h = hashlib.new(algorithm)
    with open(path, "rb") as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def compute_md5(path: Path, chunk_size: int = 8192) -> str:
    """Compute MD5 hex digest of a file.

    Args:
        path: Path to the file.
        chunk_size: Read buffer size in bytes.

    Returns:
        MD5 hex digest string (32 characters).
    """
    return _compute_hash(path, "md5", chunk_size)


def compute_sha256(path: Path, chunk_size: int = 8192) -> str:
    """Compute SHA-256 hex digest of a file.

    Args:
        path: Path to the file.
        chunk_size: Read buffer size in bytes.

    Returns:
        SHA-256 hex digest string (64 characters).
    """
    return _compute_hash(path, "sha256", chunk_size)


def verify_checksum(path: Path, expected: str, algorithm: str = "sha256") -> bool:
    """Verify a file's checksum against an expected value.

    Args:
        path: Path to the file to verify.
        expected: Expected hex digest string.
        algorithm: Hash algorithm to use. Default "sha256".

    Returns:
        True if the computed checksum matches the expected value.
    """
    actual = _compute_hash(path, algorithm)
    match = actual.lower() == expected.lower().strip()
    if not match:
        logger.warning("Checksum mismatch for %s: expected %s, got %s", path, expected, actual)
    return match


def write_checksum_file(path: Path, algorithm: str = "sha256") -> Path:
    """Write a checksum sidecar file next to the target file.

    Creates a file like ``myfile.txt.sha256`` containing the hex digest
    and filename in the standard format: ``<hash>  <filename>``.

    Args:
        path: Path to the file to checksum.
        algorithm: Hash algorithm to use. Default "sha256".

    Returns:
        Path to the created checksum sidecar file.
    """
    path = Path(path)
    digest = _compute_hash(path, algorithm)
    sidecar = path.parent / f"{path.name}.{algorithm}"
    sidecar.write_text(f"{digest}  {path.name}\n", encoding="utf-8")
    logger.debug("Wrote %s checksum to %s", algorithm, sidecar)
    return sidecar


def verify_checksum_file(path: Path) -> bool:
    """Verify a file against its sidecar checksum file.

    Looks for ``<filename>.sha256`` or ``<filename>.md5`` sidecar files
    and verifies the file's integrity against the stored checksum.

    Args:
        path: Path to the file to verify.

    Returns:
        True if the checksum matches.

    Raises:
        FileNotFoundError: If neither .sha256 nor .md5 sidecar file is found.
    """
    path = Path(path)

    for algorithm in ("sha256", "md5"):
        sidecar = path.parent / f"{path.name}.{algorithm}"
        if sidecar.is_file():
            content = sidecar.read_text(encoding="utf-8").strip()
            # Standard format: "<hash>  <filename>" or just "<hash>"
            expected_hash = content.split()[0]
            return verify_checksum(path, expected_hash, algorithm)

    raise FileNotFoundError(f"No checksum sidecar file found for {path} (tried .sha256, .md5)")


def compute_checksums_batch(paths: list[Path], algorithm: str = "sha256") -> dict[str, str]:
    """Compute checksums for multiple files.

    Args:
        paths: List of file paths to checksum.
        algorithm: Hash algorithm to use. Default "sha256".

    Returns:
        Dictionary mapping file path strings to hex digest strings.
        Files that do not exist are logged as warnings and skipped.
    """
    results: dict[str, str] = {}
    for p in paths:
        p = Path(p)
        try:
            results[str(p)] = _compute_hash(p, algorithm)
        except FileNotFoundError:
            logger.warning("Skipping missing file in batch checksum: %s", p)
        except Exception:
            logger.exception("Error computing checksum for %s", p)
    return results
