from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .amalgkit import AmalgkitParams


@dataclass(slots=True)
class SpeciesProfile:
    """Species profile configuration for RNA-seq workflows.
    
    Attributes:
        name: Species name (e.g., "Homo sapiens")
        taxon_id: NCBI taxonomy ID (optional)
        tissues: List of tissue types to filter (optional)
    """
    name: str
    taxon_id: int | None = None
    tissues: list[str] | None = None


@dataclass(slots=True)
class AmalgkitRunLayout:
    """Directory layout for amalgkit workflow outputs.
    
    Attributes:
        base_dir: Base directory for all workflow outputs
    """
    base_dir: Path

    @property
    def fastq_dir(self) -> Path:
        """Path to FASTQ files directory."""
        return self.base_dir / "fastq"

    @property
    def quant_dir(self) -> Path:
        """Path to quantification results directory."""
        return self.base_dir / "quant"

    @property
    def merge_table(self) -> Path:
        """Path to merged abundance table file."""
        return self.base_dir / "merged_abundance.tsv"

    @property
    def cstmm_dir(self) -> Path:
        """Path to CSTMM statistical test results directory."""
        return self.base_dir / "cstmm"

    @property
    def curate_dir(self) -> Path:
        """Path to curation results directory."""
        return self.base_dir / "curate"

    @property
    def csca_dir(self) -> Path:
        """Path to CSCA statistical test results directory."""
        return self.base_dir / "csca"


def build_step_params(species: SpeciesProfile, layout: AmalgkitRunLayout) -> dict[str, AmalgkitParams]:
    """Construct a per-step parameters map from high-level species/tissue layout.

    This function only builds parameter dictionaries; it does not validate the
    exact amalgkit flags. The values are converted to CLI flags by the wrapper.
    """
    md: dict[str, Any] = {}
    if species.taxon_id is not None:
        md["taxon-id"] = species.taxon_id
    if species.tissues:
        md["tissue"] = list(species.tissues)

    params: dict[str, AmalgkitParams] = {
        "metadata": md,
        "integrate": {},
        "config": {},
        "select": dict(md),
        "getfastq": {"out-dir": layout.fastq_dir},
        "quant": {"out-dir": layout.quant_dir},
        "merge": {"out": layout.merge_table},
        "cstmm": {"out-dir": layout.cstmm_dir},
        "curate": {"out-dir": layout.curate_dir},
        "csca": {"out-dir": layout.csca_dir},
        "sanity": {},
    }
    return params


__all__ = [
    "SpeciesProfile",
    "AmalgkitRunLayout",
    "build_step_params",
]
