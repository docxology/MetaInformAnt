"""GWAS-facing adapter for RNA quantification discovery.

The GWAS loader needs to consume Amalgkit/Kallisto output, but domain modules
must not import one another directly.  Keeping this dependency in an explicit
workflow adapter makes the integration boundary visible to static checks and
to readers of the package layout.
"""

from __future__ import annotations

from pathlib import Path

from metainformant.rna.core.sample_utils import find_quantification_file as _find_quantification_file


def find_quantification_file(
    sample_quant_dir: Path,
    sample_id: str | None = None,
    *,
    require_nonempty: bool = True,
) -> Path | None:
    """Find a recognized per-sample quantification file for GWAS loading."""
    return _find_quantification_file(sample_quant_dir, sample_id, require_nonempty=require_nonempty)
