"""Data extraction utilities for GWAS pipeline."""

from __future__ import annotations

import re
import shutil
import zipfile
from pathlib import Path

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Regex to parse BeeWAS sample names
SAMPLE_RE = re.compile(
    r"^(?P<colony>[A-Z]\d+)(?P<biological_group_code>G|ITQ|ITW|IV|WORK)_(?P<read>R[12])\.fastq"
)
STANDALONE_RE = re.compile(
    r"^(?P<colony>[A-Z]\d+)(?P<biological_group_code>G|ITQ|ITW|IV|WORK)_(?P<read>R[12])\.fastq-\d+\.gz$"
)


def parse_fastq_filename(filename: str) -> dict[str, str] | None:
    """Parse a BeeWAS FASTQ filename into colony, biological group, and read components."""
    basename = Path(filename).name

    m = STANDALONE_RE.match(basename)
    if m:
        return {
            "colony": m.group("colony"),
            "biological_group_code": m.group("biological_group_code"),
            "read": m.group("read"),
            "sample_id": f"{m.group('colony')}{m.group('biological_group_code')}",
        }

    m = SAMPLE_RE.match(basename)
    if m:
        return {
            "colony": m.group("colony"),
            "biological_group_code": m.group("biological_group_code"),
            "read": m.group("read"),
            "sample_id": f"{m.group('colony')}{m.group('biological_group_code')}",
        }
    return None


def lineage_from_colony(colony: str) -> str:
    """Infer lineage group from colony prefix letter."""
    prefix = colony[0]
    return {"I": "I", "M": "M", "R": "R"}.get(prefix, "unknown")


def build_sample_file_map(zip_dir: Path, accessions: set[str] | None = None) -> dict[str, dict]:
    """Map each sample_id -> {R1: source_info, R2: source_info}."""
    samples: dict[str, dict] = {}

    if not zip_dir.exists():
        return samples

    for item in sorted(zip_dir.iterdir()):
        if item.name.startswith("."):
            continue

        entries = []
        if item.suffix == ".zip":
            try:
                with zipfile.ZipFile(item, "r") as zf:
                    for info in zf.infolist():
                        if info.is_dir():
                            continue
                        parsed = parse_fastq_filename(info.filename)
                        if parsed and (accessions is None or parsed["sample_id"] in accessions):
                            entries.append(
                                {
                                    "source_type": "zip",
                                    "source_path": str(item),
                                    "inner_filename": info.filename,
                                    "inner_size": info.file_size,
                                    **parsed,
                                }
                            )
            except zipfile.BadZipFile:
                logger.error(f"Bad ZIP: {item}")
        elif item.name.endswith(".gz"):
            parsed = parse_fastq_filename(item.name)
            if parsed and (accessions is None or parsed["sample_id"] in accessions):
                entries.append(
                    {
                        "source_type": "standalone",
                        "source_path": str(item),
                        "inner_filename": item.name,
                        "inner_size": item.stat().st_size,
                        **parsed,
                    }
                )

        for e in entries:
            sid = e["sample_id"]
            if sid not in samples:
                samples[sid] = {"R1": None, "R2": None}
            samples[sid][e["read"]] = e

    return samples


def extract_fastq(entry: dict, dest_path: Path) -> None:
    """Extract a single FASTQ file to dest_path securely and resume-safely."""
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    if dest_path.exists():
        logger.debug(f"Already exists: {dest_path.name}")
        return

    logger.info(f"Extracting -> {dest_path.name}")
    if entry["source_type"] == "zip":
        source = Path(entry["source_path"])
        with zipfile.ZipFile(source, "r") as zf:
            with zf.open(entry["inner_filename"]) as src, open(dest_path, "wb") as dst:
                shutil.copyfileobj(src, dst)
    else:
        shutil.copy2(entry["source_path"], dest_path)
