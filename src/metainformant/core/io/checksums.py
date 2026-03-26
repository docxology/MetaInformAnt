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

    Computes the MD5 hash of a file for integrity verification. MD5 is fast
    but cryptographically broken - suitable for quick corruption detection
    but not for security. For secure hashing, use compute_sha256().

    Args:
        path: Path to the file to hash.
        chunk_size: Read buffer size in bytes. Default 8192 (8KB). Larger values
            are faster but use more memory.

    Returns:
        MD5 hex digest string (32 characters).

    Raises:
        FileNotFoundError: If the file does not exist.

    Example:
        >>> from metainformant.core.io.checksums import compute_md5
        >>> hash_val = compute_md5("data/sample.vcf.gz")
        >>> print(f"MD5: {hash_val}")
        MD5: d41d8cd98f00b204e9800998ecf8427e

    Note:
        MD5 is useful for quick integrity checks of large files (FASTQ, BAM, VCF).
        For secure verification, use SHA-256 instead.
    """
    return _compute_hash(path, "md5", chunk_size)


def compute_sha256(path: Path, chunk_size: int = 8192) -> str:
    """Compute SHA-256 hex digest of a file.

    Computes the SHA-256 hash of a file for integrity verification.
    SHA-256 is cryptographically secure and suitable for both integrity
    verification and security purposes.

    Args:
        path: Path to the file to hash.
        chunk_size: Read buffer size in bytes. Default 8192 (8KB). Larger values
            are faster but use more memory.

    Returns:
        SHA-256 hex digest string (64 characters).

    Raises:
        FileNotFoundError: If the file does not exist.

    Example:
        >>> from metainformant.core.io.checksums import compute_sha256
        >>> hash_val = compute_sha256("data/reference.fasta")
        >>> print(f"SHA-256: {hash_val}")
        SHA-256: e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855

    Use Cases:
        - Verifying reference genome downloads
        - Ensuring VCF/FASTQ file integrity after transfer
        - Detecting data corruption in long-running pipelines
    """
    return _compute_hash(path, "sha256", chunk_size)


def verify_checksum(path: Path, expected: str, algorithm: str = "sha256") -> bool:
    """Verify a file's checksum against an expected value.

    Compares the computed hash of a file against an expected value to
    detect corruption or verify integrity. Case-insensitive comparison
    handles both uppercase and lowercase hex strings.

    Args:
        path: Path to the file to verify.
        expected: Expected hex digest string. Can be in upper or lower case.
        algorithm: Hash algorithm to use. Default "sha256". Supported:
            md5, sha256, sha1, sha512.

    Returns:
        True if the computed checksum matches the expected value.

    Example:
        >>> from metainformant.core.io.checksums import verify_checksum
        >>> # Verify a downloaded reference genome
        >>> expected = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        >>> is_valid = verify_checksum("reference.fasta", expected)
        >>> print(f"Valid: {is_valid}")
        Valid: True

    Use Cases:
        - Verify downloaded reference genomes match expected hashes
        - Check FASTQ files after SRA download
        - Validate VCF files after transfer
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

    Efficiently computes checksums for multiple files, useful for
    verifying entire directories of bioinformatics files. Files that
    cannot be accessed are logged as warnings and skipped rather than
    causing the entire batch to fail.

    Args:
        paths: List of file paths to checksum. Can be strings or Path objects.
        algorithm: Hash algorithm to use. Default "sha256". Supported:
            md5, sha256, sha1, sha512.

    Returns:
        Dictionary mapping file path strings to hex digest strings.
        Files that do not exist or cannot be accessed are logged as warnings
        and skipped.

    Example:
        >>> from pathlib import Path
        >>> from metainformant.core.io.checksums import compute_checksums_batch
        >>> files = list(Path("data").glob("*.vcf.gz"))
        >>> hashes = compute_checksums_batch(files)
        >>> for path, hash_val in hashes.items():
        ...     print(f"{path}: {hash_val}")
        data/sample1.vcf.gz: d41d8cd98f00b204e9800998ecf8427e
        data/sample2.vcf.gz: e3b0c44298fc1c149afbf4c8996fb924

    Use Cases:
        - Verify integrity of all files in a directory after bulk download
        - Generate checksums for all reference genomes
        - Pre-flight validation before starting large workflows
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
