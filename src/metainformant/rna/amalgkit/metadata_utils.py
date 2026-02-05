"""Utilities for Amalgkit metadata manipulation."""

from pathlib import Path

import pandas as pd

from metainformant.core import logging

logger = logging.get_logger(__name__)


def deduplicate_metadata(file_path: str | Path, output_path: str | Path = None) -> bool:
    """Deduplicate amalgkit metadata table and fix placeholder scientific names.

    Args:
        file_path: Path to the metadata file (TSV)
        output_path: Path to save the cleaned file. If None, overwrites input.

    Returns:
        True if successful, False otherwise.
    """
    path = Path(file_path)
    if not path.exists():
        logger.warning(f"Metadata file not found for deduplication: {file_path}")
        return False

    if output_path is None:
        output_path = file_path
    else:
        output_path = Path(output_path)

    try:
        # Load metadata
        df = pd.read_csv(path, sep="\t")
        initial_count = len(df)

        if "scientific_name" in df.columns:
            # Drop rows with placeholder scientific name if duplicates exist
            placeholder = "Please add in format: Genus species"

            # If we find placeholders, we only keep them if there's no better name for that run
            # But in our case, amalgkit integrate ADDs these placeholders even if the original exists
            # So we prefer non-placeholder names.

            # Create a sort key that puts placeholders last (1) and real names first (0)
            df["_is_placeholder"] = (df["scientific_name"] == placeholder).astype(int)

            # Sort by run first, then by placeholder status (real names first)
            df = df.sort_values(by=["run", "_is_placeholder"], ascending=[True, True])

            # Deduplicate by run, keeping the first (which should be the real name if available)
            df = df.drop_duplicates(subset=["run"], keep="first")

            # Remove the helper column
            df = df.drop(columns=["_is_placeholder"])

            # Final check: if we still have a placeholder as the ONLY name, log it
            if (df["scientific_name"] == placeholder).any():
                logger.info(f"Found {sum(df['scientific_name'] == placeholder)} remaining placeholders in {path.name}")
        else:
            # Fallback deduplication by run only
            if "run" in df.columns:
                df = df.drop_duplicates(subset=["run"], keep="first")

        final_count = len(df)
        if initial_count > final_count:
            logger.info(f"Deduplicated {path.name}: {initial_count} -> {final_count} rows")

        # Save cleaned metadata
        df.to_csv(output_path, sep="\t", index=False)
        return True

    except Exception as e:
        logger.error(f"Error deduplicating metadata {file_path}: {e}")
        return False
