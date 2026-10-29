"""I/O utilities for single-cell transcriptomics.

This module provides functions for fetching datasets from GEO/SRA, handling fastq
compression/decompression, and running Salmon Alevin.

The ``run_*`` helpers for R and pySCENIC tools (Seurat, Monocle3, clusterProfiler, ShinyCell,
pySCENIC) are not implemented. They raise ``NotImplementedError`` and create nothing at the
requested output path, because an empty file or directory there would be mistaken for a real
result by the next pipeline step. :func:`fetch_atlas_datasets` follows the same rule for the datasets
that live in ArrayExpress, for which no fetcher exists.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import List, NoReturn

from metainformant.core.utils import logging as core_logging

logger = core_logging.get_logger(__name__)


def _not_wrapped(function: str, tool: str, guidance: str) -> NoReturn:
    """Raise ``NotImplementedError`` for a helper whose external tool metainformant does not run.

    Nothing is created at the caller's output path: a placeholder file or directory would look like
    the tool's real output.
    """
    raise NotImplementedError(f"{function}() is not implemented: metainformant does not run {tool}. {guidance}")


def download_geo_supplementary(accession: str, output_dir: str | Path) -> Path:
    """Download supplementary files for a GEO accession using wget."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Standard GEO FTP pattern
    # ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE107nnn/GSE107451/suppl/
    nnn = accession[:-3] + "nnn"
    url = f"ftp://ftp.ncbi.nlm.nih.gov/geo/series/{nnn}/{accession}/suppl/"

    logger.info(f"Checking supplementary files for {accession}...")
    cmd = ["wget", "-r", "-np", "-nd", "-q", "--show-progress", "-P", str(output_dir), url]
    logger.info(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    return output_dir


def download_sra_fastqs(accession: str, output_dir: str | Path, threads: int = 1) -> List[Path]:
    """Download FASTQ files for an SRA accession using fasterq-dump and compress them."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Downloading {accession} using fasterq-dump...")
    cmd = ["fasterq-dump", "--outdir", str(output_dir), "--threads", str(threads), "--progress", accession]
    subprocess.run(cmd, check=True)

    # Identify downloaded FASTQ files
    fastqs = list(output_dir.glob(f"{accession}*.fastq"))
    logger.info(f"Downloaded {len(fastqs)} FASTQ files for {accession}.")

    # Compress files
    for fq in fastqs:
        _compress_fastq(fq)

    return list(output_dir.glob(f"{accession}*.fastq.gz"))


def _compress_fastq(path: Path) -> None:
    """Compress a FASTQ file using pigz or gzip."""
    # Check for pigz availability
    try:
        subprocess.run(["pigz", "--version"], capture_output=True, check=True)
        compressor = "pigz"
    except (subprocess.CalledProcessError, FileNotFoundError):
        logger.warning("pigz not found, using standard gzip...")
        compressor = "gzip"

    logger.info(f"Compressing {path.name} with {compressor}...")

    # USE -f TO PREVENT STALLING ON OVERWRITE
    cmd = [compressor, "-f", str(path)]
    subprocess.run(cmd, check=True)


#: Atlas dataset name -> accession. ``GSE`` accessions are fetched from GEO; ``E-`` accessions are ArrayExpress.
_ATLAS_ACCESSIONS = {
    "Davie-2018": "GSE107451",
    "Baker-2021": "GSE152495",
    "Park-2022": "GSE207799",
    "Dopp-2024": "GSE221239",
    "LeeBenton-2023": "GSE247965",
    "Li-2022": "E-MTAB-10519",
}


def fetch_atlas_datasets(datasets: List[str], output_base: str | Path) -> List[Path]:
    """Fetch the GEO supplementary files of atlas datasets (Davie-2018, Baker-2021, etc.) to the output base.

    Each known dataset is downloaded into ``output_base / <name>`` with :func:`download_geo_supplementary`.
    Unknown names are skipped with a warning. A dataset whose accession is in ArrayExpress (an ``E-``
    accession such as ``Li-2022``) cannot be fetched, because there is no ArrayExpress downloader. That is
    an error, not a skip: an empty directory in the returned list would pass for a downloaded dataset.

    Args:
        datasets: Names of the atlas datasets to fetch.
        output_base: Directory that receives one sub-directory per dataset.

    Returns:
        The directory of each fetched dataset, in the order requested (unknown names left out).

    Raises:
        NotImplementedError: If a requested dataset is an ArrayExpress accession. The check runs before
            anything is created or downloaded, so the call leaves ``output_base`` untouched.
        subprocess.CalledProcessError: If a download fails.
    """
    output_base = Path(output_base)

    unfetchable = [name for name in datasets if _ATLAS_ACCESSIONS.get(name, "").startswith("E-")]
    if unfetchable:
        accessions = ", ".join(f"{name} ({_ATLAS_ACCESSIONS[name]})" for name in unfetchable)
        raise NotImplementedError(
            f"fetch_atlas_datasets() cannot fetch {accessions}: they are ArrayExpress accessions and no "
            "ArrayExpress downloader exists. Download them from https://www.ebi.ac.uk/biostudies/arrayexpress "
            "yourself, or leave them out of the request. Nothing was created or downloaded."
        )

    output_base.mkdir(parents=True, exist_ok=True)
    fetched_paths = []

    for ds_name in datasets:
        acc = _ATLAS_ACCESSIONS.get(ds_name)
        if not acc:
            logger.warning(f"Unknown dataset name: {ds_name}")
            continue

        ds_dir = output_base / ds_name
        ds_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Fetching {ds_name} ({acc}) -> {ds_dir}")

        download_geo_supplementary(acc, ds_dir)

        fetched_paths.append(ds_dir)

    return fetched_paths


def run_salmon_alevin(
    fastq_r1: Path, fastq_r2: Path, index_dir: Path, tg_map: Path, output_dir: Path, threads: int = 1
) -> None:
    """Run Salmon Alevin for single-cell quantification."""
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Running Salmon Alevin quantification into {output_dir}...")

    cmd = [
        "salmon",
        "alevin",
        "-l",
        "ISR",
        "-1",
        str(fastq_r1),
        "-2",
        str(fastq_r2),
        "--chromiumV3",  # Assuming V3, can be parameterized
        "-i",
        str(index_dir),
        "-p",
        str(threads),
        "-o",
        str(output_dir),
        "--tgMap",
        str(tg_map),
    ]
    logger.info(f"Command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def run_seurat_integration(input_dirs: dict[str, Path], output_dir: Path, resolution: float = 0.8) -> NoReturn:
    """Not implemented: Seurat integration is an R workflow that metainformant does not run.

    Raises:
        NotImplementedError: Always. Nothing is written to ``output_dir``.
    """
    _not_wrapped(
        "run_seurat_integration",
        "Seurat",
        "Run Seurat's integration (FindIntegrationAnchors / IntegrateData) yourself in R, or integrate in Python "
        "with metainformant.singlecell.data.integration (harmony_integration, bbknn_integration, "
        "scanorama_integration, mnn_integration, combat_integration).",
    )


def run_differential_expression(input_rds: Path, output_csv: Path, contrast: str) -> NoReturn:
    """Not implemented: differential expression on a Seurat RDS object is an R workflow.

    Raises:
        NotImplementedError: Always. Nothing is written to ``output_csv``.
    """
    _not_wrapped(
        "run_differential_expression",
        "R differential expression on Seurat objects",
        "Run it yourself in R (for example Seurat's FindMarkers), or use "
        "metainformant.singlecell.differential.expression (differential_expression, pseudobulk_de) "
        "on a SingleCellData.",
    )


def run_scenic_grnboost2(expression_loom: Path, tf_txt: Path, output_adjacencies: Path, threads: int = 1) -> NoReturn:
    """Not implemented: pySCENIC GRNBoost2 is not wrapped.

    Raises:
        NotImplementedError: Always. Nothing is written to ``output_adjacencies``.
    """
    _not_wrapped(
        "run_scenic_grnboost2",
        "pySCENIC",
        "Run `pyscenic grn` (GRNBoost2) yourself on the expression loom file and the transcription-factor list.",
    )


def run_scenic_cistarget(
    expression_loom: Path, adjacencies: Path, motif_database: Path, motif_annotations: Path, output_regulons: Path
) -> NoReturn:
    """Not implemented: pySCENIC cisTarget is not wrapped.

    Raises:
        NotImplementedError: Always. Nothing is written to ``output_regulons``.
    """
    _not_wrapped(
        "run_scenic_cistarget",
        "pySCENIC",
        "Run `pyscenic ctx` yourself on the adjacencies with your motif ranking database and motif annotations.",
    )


def run_scenic_aucell(expression_loom: Path, regulons: Path, output_loom: Path) -> NoReturn:
    """Not implemented: pySCENIC AUCell is not wrapped.

    Raises:
        NotImplementedError: Always. Nothing is written to ``output_loom``.
    """
    _not_wrapped(
        "run_scenic_aucell",
        "pySCENIC",
        "Run `pyscenic aucell` yourself on the expression loom file and the regulons.",
    )


def filter_regulons_post_scenic(regulons_csv: Path, filtered_output: Path, raw_expression_csv: Path) -> NoReturn:
    """Not implemented: there is no post-SCENIC regulon filter.

    Raises:
        NotImplementedError: Always. Nothing is written to ``filtered_output``.
    """
    _not_wrapped(
        "filter_regulons_post_scenic",
        "a post-SCENIC regulon filter",
        "Filter the pySCENIC regulon table yourself (for example with pandas) using your own criteria.",
    )


def run_pseudotime_trajectory(input_rds: Path, output_rds: Path, root_selector: str) -> NoReturn:
    """Not implemented: Monocle3 trajectory inference is an R workflow that metainformant does not run.

    Raises:
        NotImplementedError: Always. Nothing is written to ``output_rds``.
    """
    _not_wrapped(
        "run_pseudotime_trajectory",
        "Monocle3",
        "Run Monocle3 yourself in R, or use metainformant.singlecell.analysis.trajectory "
        "(compute_diffusion_pseudotime, paga_trajectory, slingshot_trajectory) on a SingleCellData.",
    )


def run_go_enrichment(input_deseq_csv: Path, output_dir: Path, organism: str) -> NoReturn:
    """Not implemented: clusterProfiler GO enrichment is an R workflow that metainformant does not run.

    Raises:
        NotImplementedError: Always. ``output_dir`` is not created.
    """
    _not_wrapped(
        "run_go_enrichment",
        "clusterProfiler",
        "Run clusterProfiler yourself in R on the differential-expression table, or use the Python "
        "enrichment helpers in metainformant.ontology (for example metainformant.ontology.core.go.enrich_genes).",
    )


def run_shinycell_export(input_rds: Path, output_dir: Path) -> NoReturn:
    """Not implemented: the ShinyCell export is an R workflow that metainformant does not run.

    Raises:
        NotImplementedError: Always. ``output_dir`` is not created.
    """
    _not_wrapped(
        "run_shinycell_export",
        "ShinyCell",
        "Run the ShinyCell R package yourself on the Seurat object.",
    )
