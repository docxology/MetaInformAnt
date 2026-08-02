from pathlib import Path

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.rna.engine.progress_db import ProgressDB
from metainformant.rna.engine.species import configured_data_root


def main() -> None:
    """Print a summary directly from the current progress database."""

    data_root = configured_data_root()
    db_path = data_root / "pipeline_progress.db"
    if not db_path.is_file():
        raise SystemExit(f"Progress database not found: {db_path}")
    db = ProgressDB(db_path)
    counts = db.get_counts()
    summary = {
        "total_species": len(counts),
        "completed": 0,
        "failed": 0,
        "unquantified_samples": 0,
        "quantified_samples": 0,
    }
    for species_counts in counts.values():
        total = sum(species_counts.values())
        quantified = species_counts.get("quantified", 0)
        summary["quantified_samples"] += quantified
        summary["unquantified_samples"] += total - quantified
        summary["failed"] += species_counts.get("failed", 0)
        if total > 0 and quantified == total:
            summary["completed"] += 1
    db.close()
    print("Sample Status Summary:")
    for k, v in summary.items():
        print(f"{k}: {v}")


if __name__ == "__main__":
    main()
