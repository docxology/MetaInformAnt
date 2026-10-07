"""Real-file cohort identity boundaries for the within-species orchestrator."""

from pathlib import Path

import pandas as pd
import pytest

from metainformant.rna.analysis.within_species_orchestrator import (
    WithinSpeciesOrchestrator,
)


@pytest.mark.parametrize("runs", [["a"], ["a", "a", "b"]])
def test_metadata_cannot_drop_or_duplicate_expression_samples(tmp_path: Path, runs: list[str]) -> None:
    # Given a two-sample expression matrix and incomplete or duplicated metadata.
    abundance = tmp_path / "abundance.tsv"
    metadata = tmp_path / "metadata.tsv"
    pd.DataFrame({"a": [1, 2], "b": [3, 4]}, index=["g1", "g2"]).to_csv(abundance, sep="\t")
    pd.DataFrame({"run": runs, "tissue": ["brain"] * len(runs)}).to_csv(metadata, sep="\t", index=False)
    orchestrator = WithinSpeciesOrchestrator("species", abundance, metadata, tmp_path / "out")
    # When loaded; then the denominator cannot silently change.
    with pytest.raises(ValueError):
        orchestrator.load_data()
