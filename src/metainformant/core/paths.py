from __future__ import annotations

from pathlib import Path


def expand_and_resolve(path: str | Path) -> Path:
    """Expand user (~) and resolve to absolute path without requiring existence."""
    return Path(path).expanduser().resolve(strict=False)


def is_within(path: str | Path, parent: str | Path) -> bool:
    """Return True if path is inside parent directory (after resolving)."""
    p = expand_and_resolve(path)
    r = expand_and_resolve(parent)
    try:
        p.relative_to(r)
        return True
    except ValueError:
        return False


def ensure_directory(path: Path) -> None:
    """Create directory and any missing parent directories.

    Args:
        path: Directory path to create
    """
    path.mkdir(parents=True, exist_ok=True)


def prepare_file_path(file_path: Path) -> None:
    """Ensure parent directories exist for a file path.

    Args:
        file_path: File path to prepare
    """
    if file_path.parent:
        ensure_directory(file_path.parent)


def is_safe_path(path: str) -> bool:
    """Check if path is safe (no path traversal attempts).

    Args:
        path: Path string to check

    Returns:
        True if path appears safe
    """
    # Check for path traversal attempts
    if ".." in path or path.startswith("/etc/") or path.startswith("/root/"):
        return False

    # Check for other potentially dangerous patterns
    dangerous_patterns = ["..//", "..\\", ";", "|", "&", "$"]
    for pattern in dangerous_patterns:
        if pattern in path:
            return False

    return True


def get_file_extension(filename: str) -> str:
    """Get file extension from filename.

    Args:
        filename: Name of file

    Returns:
        File extension including the dot, or empty string
    """
    return Path(filename).suffix


def change_extension(path: str, new_extension: str) -> Path:
    """Change file extension.

    Args:
        path: Original file path
        new_extension: New extension (with or without dot)

    Returns:
        Path with new extension
    """
    p = Path(path)
    if not new_extension.startswith("."):
        new_extension = "." + new_extension

    return p.with_suffix(new_extension)


def find_files_by_extension(directory: str | Path, extension: str) -> list[Path]:
    """Find all files with a specific extension in a directory.

    Args:
        directory: Directory to search
        extension: File extension (with or without dot)

    Returns:
        List of matching file paths
    """
    if not extension.startswith("."):
        extension = "." + extension

    dir_path = Path(directory)
    return list(dir_path.rglob(f"*{extension}"))


def get_file_size(path: str | Path) -> int:
    """Get file size in bytes.

    Args:
        path: Path to file

    Returns:
        File size in bytes, or 0 if file doesn't exist
    """
    try:
        return Path(path).stat().st_size
    except (OSError, FileNotFoundError):
        return 0


def get_directory_size(path: str | Path) -> int:
    """Get total size of all files in a directory.

    Args:
        path: Directory path

    Returns:
        Total size in bytes
    """
    total_size = 0
    dir_path = Path(path)

    if not dir_path.exists():
        return 0

    for file_path in dir_path.rglob("*"):
        if file_path.is_file():
            try:
                total_size += file_path.stat().st_size
            except (OSError, FileNotFoundError):
                continue

    return total_size


def sanitize_filename(filename: str) -> str:
    """Sanitize filename for safe filesystem use.

    Args:
        filename: Original filename

    Returns:
        Sanitized filename safe for filesystem
    """
    import re

    # Remove or replace dangerous characters
    sanitized = re.sub(r'[<>:"/\\|?*]', '_', filename)
    # Remove control characters
    sanitized = re.sub(r'[\x00-\x1f\x7f-\x9f]', '', sanitized)
    # Remove leading/trailing dots and spaces
    sanitized = sanitized.strip(' .')

    # Ensure filename is not empty
    if not sanitized:
        sanitized = "untitled"

    # Truncate if too long (255 chars max for most filesystems)
    if len(sanitized) > 255:
        name_part, ext = Path(sanitized).stem, Path(sanitized).suffix
        max_name = 255 - len(ext)
        sanitized = name_part[:max_name] + ext

    return sanitized


def create_temp_file(suffix: str = "", prefix: str = "tmp", directory: str | Path | None = None) -> Path:
    """Create a temporary file path that doesn't exist yet.

    Args:
        suffix: File extension
        prefix: Filename prefix
        directory: Directory for temp file (uses system temp if None)

    Returns:
        Path to non-existent temporary file
    """
    import tempfile

    if directory is None:
        directory = tempfile.gettempdir()

    dir_path = Path(directory)
    ensure_directory(dir_path)

    while True:
        temp_name = f"{prefix}_{hash(str(Path.home()))}_{suffix}"
        temp_path = dir_path / temp_name
        if not temp_path.exists():
            return temp_path


#
# Output pattern discovery functions
#


def discover_output_patterns(module_name: str) -> dict[str, Any]:
    """Get output directory patterns for a module.

    Args:
        module_name: Name of the module (e.g., 'rna', 'gwas', 'dna')

    Returns:
        Dictionary with output pattern information:
        - base_pattern: Base output pattern
        - subdirs: List of common subdirectories
        - examples: Example output paths
    """
    # Map from .cursorrules patterns
    patterns_map: dict[str, dict[str, Any]] = {
        "rna": {
            "base_pattern": "output/amalgkit/<species>/<step>/",
            "subdirs": ["quant/", "work/", "logs/"],
            "examples": ["output/amalgkit/Apis_mellifera/quant/", "output/amalgkit/Apis_mellifera/work/"],
        },
        "gwas": {
            "base_pattern": "output/gwas/<analysis_type>/",
            "subdirs": ["association/", "plots/", "qc/"],
            "examples": ["output/gwas/association/", "output/gwas/plots/"],
        },
        "life_events": {
            "base_pattern": "output/life_events/<workflow>/",
            "subdirs": ["embeddings/", "models/", "plots/"],
            "examples": ["output/life_events/embeddings/", "output/life_events/models/"],
        },
        "multiomics": {
            "base_pattern": "output/multiomics/<integration>/",
            "subdirs": ["integrated/", "plots/"],
            "examples": ["output/multiomics/integrated/", "output/multiomics/plots/"],
        },
        "singlecell": {
            "base_pattern": "output/singlecell/<analysis>/",
            "subdirs": ["preprocessing/", "clustering/"],
            "examples": ["output/singlecell/preprocessing/", "output/singlecell/clustering/"],
        },
        "networks": {
            "base_pattern": "output/networks/<network_type>/",
            "subdirs": ["ppi/", "regulatory/", "pathways/"],
            "examples": ["output/networks/ppi/", "output/networks/regulatory/"],
        },
        "information": {
            "base_pattern": "output/information/<analysis_type>/",
            "subdirs": ["entropy/", "complexity/"],
            "examples": ["output/information/entropy/", "output/information/complexity/"],
        },
        "dna": {
            "base_pattern": "output/dna/<analysis_type>/",
            "subdirs": ["phylogeny/", "population/", "variants/"],
            "examples": ["output/dna/phylogeny/", "output/dna/population/"],
        },
        "protein": {
            "base_pattern": "output/protein/<analysis_type>/",
            "subdirs": ["structures/", "alignments/"],
            "examples": ["output/protein/structures/", "output/protein/alignments/"],
        },
        "quality": {
            "base_pattern": "output/quality/<dataset>/",
            "subdirs": ["fastq/", "reports/"],
            "examples": ["output/quality/fastq/", "output/quality/reports/"],
        },
        "math": {
            "base_pattern": "output/math/<type>/",
            "subdirs": ["simulations/", "models/", "plots/"],
            "examples": ["output/math/simulations/", "output/math/models/"],
        },
        "simulation": {
            "base_pattern": "output/simulation/<type>/",
            "subdirs": ["sequences/", "ecosystems/"],
            "examples": ["output/simulation/sequences/", "output/simulation/ecosystems/"],
        },
    }

    # Return known pattern or default
    if module_name.lower() in patterns_map:
        return patterns_map[module_name.lower()]

    # Default pattern
    return {
        "base_pattern": f"output/{module_name}/",
        "subdirs": [],
        "examples": [f"output/{module_name}/"],
    }


def find_output_locations(repo_root: str | Path, pattern: str | None = None) -> list[Path]:
    """Find existing output directories.

    Args:
        repo_root: Root directory of repository
        pattern: Optional pattern to match (e.g., 'rna', 'gwas')

    Returns:
        List of Path objects for existing output directories
    """
    repo_root = Path(repo_root)
    output_dir = repo_root / "output"
    if not output_dir.exists():
        return []

    locations: list[Path] = []

    if pattern:
        # Search for directories matching pattern
        for path in output_dir.rglob("*"):
            if path.is_dir() and pattern.lower() in str(path).lower():
                locations.append(path)
    else:
        # Return all output subdirectories
        for path in output_dir.rglob("*"):
            if path.is_dir():
                locations.append(path)

    return sorted(set(locations))


def get_module_output_base(module_name: str) -> str:
    """Get default output base for module.

    Args:
        module_name: Name of the module

    Returns:
        Default output base path string
    """
    patterns = discover_output_patterns(module_name)
    base_pattern = patterns.get("base_pattern", f"output/{module_name}/")
    # Remove template variables for base
    base = base_pattern.split("<")[0].rstrip("/")
    return base if base else f"output/{module_name}"


def list_output_structure(repo_root: str | Path) -> dict[str, Any]:
    """Map entire output directory structure.

    Args:
        repo_root: Root directory of repository

    Returns:
        Dictionary with output structure information:
        - total_dirs: Total number of directories
        - total_files: Total number of files
        - total_size: Total size in bytes
        - structure: Nested structure of directories
    """
    repo_root = Path(repo_root)
    output_dir = repo_root / "output"
    if not output_dir.exists():
        return {
            "total_dirs": 0,
            "total_files": 0,
            "total_size": 0,
            "structure": {},
        }

    total_dirs = 0
    total_files = 0
    total_size = 0
    structure: dict[str, Any] = {}

    def build_structure(path: Path, rel_path: Path) -> dict[str, Any]:
        """Recursively build structure."""
        nonlocal total_dirs, total_files, total_size

        node: dict[str, Any] = {
            "type": "directory" if path.is_dir() else "file",
            "path": str(rel_path),
        }

        if path.is_dir():
            total_dirs += 1
            node["children"] = {}
            try:
                for child in sorted(path.iterdir()):
                    if child.name.startswith("."):
                        continue
                    child_rel = rel_path / child.name
                    node["children"][child.name] = build_structure(child, child_rel)
            except (OSError, PermissionError):
                pass
        else:
            total_files += 1
            try:
                size = path.stat().st_size
                total_size += size
                node["size"] = size
            except OSError:
                node["size"] = 0

        return node

    try:
        structure = build_structure(output_dir, Path("output"))
    except Exception:
        pass

    return {
        "total_dirs": total_dirs,
        "total_files": total_files,
        "total_size": total_size,
        "structure": structure,
    }
