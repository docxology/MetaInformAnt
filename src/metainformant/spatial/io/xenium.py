"""10x Xenium in situ spatial transcriptomics data loading.

Xenium provides subcellular-resolution spatial gene expression data with
per-transcript coordinates, cell segmentation boundaries, and feature matrices.

Typical Xenium outputs:
  - cell_feature_matrix/ (MEX format) or cell_feature_matrix.h5
  - transcripts.csv.gz: Per-transcript x/y/z coordinates
  - cells.csv.gz: Cell centroid positions and QC metrics
  - cell_boundaries.csv.gz or cell_boundaries.parquet: Cell polygon boundaries
  - nucleus_boundaries.csv.gz: Nucleus polygon boundaries
"""

from __future__ import annotations

import csv
import gzip
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional dependencies
try:
    import numpy as np
    from numpy.typing import NDArray
except ImportError:
    np = None  # type: ignore[assignment]
    NDArray = None  # type: ignore[assignment,misc]

try:
    from scipy import sparse as sp_sparse
    from scipy.io import mmread
except ImportError:
    sp_sparse = None  # type: ignore[assignment]
    mmread = None  # type: ignore[assignment]

try:
    import pandas as pd
except ImportError:
    pd = None  # type: ignore[assignment]


@dataclass
class XeniumTranscript:
    """A single Xenium detected transcript.

    Attributes:
        transcript_id: Unique transcript identifier.
        gene: Gene name.
        x: X-coordinate in microns.
        y: Y-coordinate in microns.
        z: Z-coordinate in microns.
        cell_id: Assigned cell ID (0 or empty if unassigned).
        quality_value: Phred-scaled quality score.
        overlaps_nucleus: Whether the transcript overlaps a segmented nucleus.
    """

    transcript_id: str
    gene: str
    x: float
    y: float
    z: float = 0.0
    cell_id: str = ""
    quality_value: float = 20.0
    overlaps_nucleus: bool = False


@dataclass
class CellBoundary:
    """Polygon boundary for a single cell.

    Attributes:
        cell_id: Unique cell identifier.
        vertices: List of (x, y) coordinate tuples forming the polygon boundary.
    """

    cell_id: str
    vertices: list[tuple[float, float]]


@dataclass
class XeniumDataset:
    """Complete 10x Xenium spatial dataset.

    Attributes:
        expression: Cell-by-gene expression matrix (numpy or scipy sparse).
        coordinates: Cell centroid coordinates (n_cells, 2).
        cell_ids: List of cell identifier strings.
        gene_names: List of gene names.
        gene_ids: List of gene IDs.
        transcripts: Optional list of individual transcript records.
        cell_boundaries: Optional list of cell boundary polygons.
        metadata: Additional dataset metadata.
    """

    expression: Any  # np.ndarray or scipy sparse
    coordinates: Any  # np.ndarray (n_cells, 2)
    cell_ids: list[str]
    gene_names: list[str]
    gene_ids: list[str] = field(default_factory=list)
    transcripts: list[XeniumTranscript] | None = None
    cell_boundaries: list[CellBoundary] | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def n_cells(self) -> int:
        """Number of cells."""
        return len(self.cell_ids)

    @property
    def n_genes(self) -> int:
        """Number of genes."""
        return len(self.gene_names)


def _open_maybe_gzipped(filepath: Path, mode: str = "rt"):
    """Open a file that may or may not be gzipped."""
    if str(filepath).endswith(".gz"):
        return gzip.open(filepath, mode)
    return open(filepath, mode, newline="")


def read_cell_features(feature_matrix_path: str | Path) -> tuple[Any, list[str], list[str], list[str]]:
    """Read a Xenium cell-level feature matrix.

    Supports both MEX directory format and HDF5 format.

    Args:
        feature_matrix_path: Path to cell_feature_matrix/ directory or .h5 file.

    Returns:
        Tuple of (expression_matrix, cell_ids, gene_names, gene_ids).
        Matrix is (n_cells, n_genes) in CSR sparse format.

    Raises:
        FileNotFoundError: If the path does not exist.
        ImportError: If required dependencies are missing.
    """
    fpath = Path(feature_matrix_path)
    if not fpath.exists():
        raise FileNotFoundError(f"Feature matrix path not found: {fpath}")

    if fpath.is_dir():
        # MEX format directory
        return _read_xenium_mex(fpath)
    elif fpath.suffix in (".h5", ".hdf5"):
        return _read_xenium_h5(fpath)
    else:
        raise ValueError(f"Unrecognized feature matrix format: {fpath}")


def _read_xenium_mex(matrix_dir: Path) -> tuple[Any, list[str], list[str], list[str]]:
    """Read Xenium MEX format feature matrix."""
    if sp_sparse is None or mmread is None:
        raise ImportError("scipy is required for MEX format: uv pip install scipy")
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")

    # Find matrix file
    mtx_file = None
    for candidate in ["matrix.mtx.gz", "matrix.mtx"]:
        if (matrix_dir / candidate).exists():
            mtx_file = matrix_dir / candidate
            break
    if mtx_file is None:
        raise FileNotFoundError(f"No matrix.mtx(.gz) found in {matrix_dir}")

    mat = mmread(str(mtx_file))
    mat = sp_sparse.csc_matrix(mat).T.tocsr()  # Transpose: (genes x cells) -> (cells x genes)

    # Read barcodes/cell IDs
    barcodes_file = None
    for candidate in ["barcodes.tsv.gz", "barcodes.tsv"]:
        if (matrix_dir / candidate).exists():
            barcodes_file = matrix_dir / candidate
            break
    if barcodes_file is None:
        raise FileNotFoundError(f"No barcodes.tsv(.gz) in {matrix_dir}")

    cell_ids: list[str] = []
    with _open_maybe_gzipped(barcodes_file) as fh:
        for line in fh:
            cell_ids.append(line.strip())

    # Read features
    features_file = None
    for candidate in ["features.tsv.gz", "features.tsv"]:
        if (matrix_dir / candidate).exists():
            features_file = matrix_dir / candidate
            break
    if features_file is None:
        raise FileNotFoundError(f"No features.tsv(.gz) in {matrix_dir}")

    gene_ids: list[str] = []
    gene_names: list[str] = []
    with _open_maybe_gzipped(features_file) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            gene_ids.append(parts[0] if len(parts) > 0 else "")
            gene_names.append(parts[1] if len(parts) > 1 else parts[0])

    logger.info(f"Read Xenium MEX: {mat.shape[0]} cells x {mat.shape[1]} genes")
    return mat, cell_ids, gene_names, gene_ids


def _read_xenium_h5(h5_path: Path) -> tuple[Any, list[str], list[str], list[str]]:
    """Read Xenium HDF5 feature matrix."""
    try:
        import h5py
    except ImportError:
        raise ImportError("h5py is required for HDF5 format: uv pip install h5py")
    if sp_sparse is None:
        raise ImportError("scipy is required: uv pip install scipy")
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")

    with h5py.File(str(h5_path), "r") as f:
        if "matrix" in f:
            group = f["matrix"]
            data = np.array(group["data"])
            indices = np.array(group["indices"])
            indptr = np.array(group["indptr"])
            shape = tuple(group["shape"])
            cell_ids = [b.decode("utf-8") if isinstance(b, bytes) else str(b) for b in group["barcodes"]]

            features = group["features"]
            gene_ids = [g.decode("utf-8") if isinstance(g, bytes) else str(g) for g in features["id"]]
            gene_names = [g.decode("utf-8") if isinstance(g, bytes) else str(g) for g in features["name"]]

            mat = sp_sparse.csc_matrix((data, indices, indptr), shape=shape).T.tocsr()
        else:
            raise ValueError(f"Unrecognized HDF5 layout in {h5_path}")

    logger.info(f"Read Xenium H5: {mat.shape[0]} cells x {mat.shape[1]} genes")
    return mat, cell_ids, gene_names, gene_ids


def read_transcripts(
    transcripts_path: str | Path,
    *,
    min_quality: float = 20.0,
) -> list[XeniumTranscript]:
    """Read per-transcript coordinates from a Xenium transcripts file.

    Args:
        transcripts_path: Path to transcripts.csv or transcripts.csv.gz.
        min_quality: Minimum Phred quality score to keep (default 20 = Q20).

    Returns:
        List of XeniumTranscript records passing the quality filter.

    Raises:
        FileNotFoundError: If the transcripts file does not exist.
    """
    tpath = Path(transcripts_path)
    if not tpath.exists():
        raise FileNotFoundError(f"Transcripts file not found: {tpath}")

    transcripts: list[XeniumTranscript] = []

    with _open_maybe_gzipped(tpath) as fh:
        reader = csv.DictReader(fh)
        headers = set(reader.fieldnames or [])

        # Detect column names (Xenium uses specific conventions)
        tid_col = _find_col(headers, ["transcript_id", "TranscriptID"])
        gene_col = _find_col(headers, ["feature_name", "gene", "Gene", "target"])
        x_col = _find_col(headers, ["x_location", "x", "X"])
        y_col = _find_col(headers, ["y_location", "y", "Y"])
        z_col = _find_col(headers, ["z_location", "z", "Z"])
        cell_col = _find_col(headers, ["cell_id", "CellID"])
        qv_col = _find_col(headers, ["qv", "quality_value", "quality"])
        nuc_col = _find_col(headers, ["overlaps_nucleus"])

        if gene_col is None or x_col is None or y_col is None:
            raise ValueError(
                f"Required columns not found in {tpath}. Available: {sorted(headers)}"
            )

        total = 0
        for row in reader:
            total += 1
            qv = float(row[qv_col]) if qv_col and row.get(qv_col) else 20.0
            if qv < min_quality:
                continue

            transcripts.append(
                XeniumTranscript(
                    transcript_id=str(row[tid_col]) if tid_col and row.get(tid_col) else str(total),
                    gene=str(row[gene_col]),
                    x=float(row[x_col]),
                    y=float(row[y_col]),
                    z=float(row[z_col]) if z_col and row.get(z_col) else 0.0,
                    cell_id=str(row[cell_col]) if cell_col and row.get(cell_col) else "",
                    quality_value=qv,
                    overlaps_nucleus=row.get(nuc_col, "0") in ("1", "True", "true") if nuc_col else False,
                )
            )

    logger.info(f"Read {len(transcripts)}/{total} transcripts (min_quality={min_quality})")
    return transcripts


def load_cell_boundaries(boundaries_path: str | Path) -> list[CellBoundary]:
    """Load cell segmentation polygon boundaries.

    Supports CSV/CSV.GZ format where each row has cell_id, vertex_x, vertex_y,
    with multiple rows per cell forming a polygon.

    Args:
        boundaries_path: Path to cell_boundaries.csv(.gz) or cell_boundaries.parquet.

    Returns:
        List of CellBoundary records with polygon vertices.

    Raises:
        FileNotFoundError: If the boundaries file does not exist.
    """
    bpath = Path(boundaries_path)
    if not bpath.exists():
        raise FileNotFoundError(f"Cell boundaries file not found: {bpath}")

    # Handle parquet format
    if str(bpath).endswith(".parquet"):
        return _load_boundaries_parquet(bpath)

    # CSV/CSV.GZ format
    cell_vertices: dict[str, list[tuple[float, float]]] = {}

    with _open_maybe_gzipped(bpath) as fh:
        reader = csv.DictReader(fh)
        headers = set(reader.fieldnames or [])

        cell_col = _find_col(headers, ["cell_id", "CellID"])
        x_col = _find_col(headers, ["vertex_x", "x", "X"])
        y_col = _find_col(headers, ["vertex_y", "y", "Y"])

        if cell_col is None or x_col is None or y_col is None:
            raise ValueError(
                f"Required columns (cell_id, vertex_x, vertex_y) not found. "
                f"Available: {sorted(headers)}"
            )

        for row in reader:
            cid = str(row[cell_col])
            vx = float(row[x_col])
            vy = float(row[y_col])
            if cid not in cell_vertices:
                cell_vertices[cid] = []
            cell_vertices[cid].append((vx, vy))

    boundaries = [
        CellBoundary(cell_id=cid, vertices=verts)
        for cid, verts in cell_vertices.items()
    ]

    logger.info(f"Loaded boundaries for {len(boundaries)} cells")
    return boundaries


def _load_boundaries_parquet(parquet_path: Path) -> list[CellBoundary]:
    """Load cell boundaries from a parquet file using pandas."""
    if pd is None:
        raise ImportError("pandas is required for parquet format: uv pip install pandas")

    df = pd.read_parquet(parquet_path)

    # Detect columns
    cell_col = None
    x_col = None
    y_col = None
    for col in df.columns:
        col_lower = col.lower()
        if "cell" in col_lower and "id" in col_lower:
            cell_col = col
        elif "vertex_x" in col_lower or (col_lower == "x"):
            x_col = col
        elif "vertex_y" in col_lower or (col_lower == "y"):
            y_col = col

    if cell_col is None or x_col is None or y_col is None:
        raise ValueError(f"Required columns not found in parquet. Available: {list(df.columns)}")

    cell_vertices: dict[str, list[tuple[float, float]]] = {}
    for _, row in df.iterrows():
        cid = str(row[cell_col])
        if cid not in cell_vertices:
            cell_vertices[cid] = []
        cell_vertices[cid].append((float(row[x_col]), float(row[y_col])))

    return [CellBoundary(cell_id=cid, vertices=verts) for cid, verts in cell_vertices.items()]


def load_xenium(
    path: str | Path,
    *,
    load_transcripts: bool = False,
    load_boundaries: bool = False,
    min_transcript_quality: float = 20.0,
) -> XeniumDataset:
    """Load a complete 10x Xenium dataset from a directory.

    Expects the Xenium Ranger output structure:
        path/
            cell_feature_matrix/ (or cell_feature_matrix.h5)
            cells.csv.gz: Cell centroid positions
            transcripts.csv.gz (optional): Per-transcript coordinates
            cell_boundaries.csv.gz (optional): Cell polygon boundaries

    Args:
        path: Path to the Xenium output directory.
        load_transcripts: If True, load individual transcript coordinates.
        load_boundaries: If True, load cell segmentation boundaries.
        min_transcript_quality: Minimum quality for transcript filtering.

    Returns:
        XeniumDataset with expression matrix, coordinates, and optional extras.

    Raises:
        FileNotFoundError: If required files are missing.
    """
    if np is None:
        raise ImportError("NumPy is required for Xenium loading: uv pip install numpy")

    base_dir = Path(path)
    if not base_dir.exists():
        raise FileNotFoundError(f"Xenium directory not found: {base_dir}")

    # --- Load feature matrix ---
    h5_path = base_dir / "cell_feature_matrix.h5"
    mex_dir = base_dir / "cell_feature_matrix"

    if h5_path.exists():
        expression, cell_ids, gene_names, gene_ids = read_cell_features(h5_path)
    elif mex_dir.exists():
        expression, cell_ids, gene_names, gene_ids = read_cell_features(mex_dir)
    else:
        raise FileNotFoundError(
            f"No cell_feature_matrix found in {base_dir}. "
            "Expected cell_feature_matrix.h5 or cell_feature_matrix/ directory."
        )

    # --- Load cell positions ---
    cells_file = None
    for candidate in ["cells.csv.gz", "cells.csv"]:
        if (base_dir / candidate).exists():
            cells_file = base_dir / candidate
            break

    coordinates: Any
    if cells_file is not None:
        cell_coords: dict[str, tuple[float, float]] = {}
        with _open_maybe_gzipped(cells_file) as fh:
            reader = csv.DictReader(fh)
            headers = set(reader.fieldnames or [])
            cid_col = _find_col(headers, ["cell_id", "CellID"])
            x_col = _find_col(headers, ["x_centroid", "x", "X"])
            y_col = _find_col(headers, ["y_centroid", "y", "Y"])

            if cid_col and x_col and y_col:
                for row in reader:
                    cell_coords[str(row[cid_col])] = (float(row[x_col]), float(row[y_col]))

        # Build coordinates in same order as cell_ids
        coord_list: list[list[float]] = []
        for cid in cell_ids:
            if cid in cell_coords:
                coord_list.append([cell_coords[cid][0], cell_coords[cid][1]])
            else:
                coord_list.append([0.0, 0.0])
        coordinates = np.array(coord_list, dtype=np.float64)
    else:
        logger.warning("No cells.csv found; coordinates will be zeros")
        coordinates = np.zeros((len(cell_ids), 2), dtype=np.float64)

    # --- Optionally load transcripts ---
    transcripts: list[XeniumTranscript] | None = None
    if load_transcripts:
        tx_file = None
        for candidate in ["transcripts.csv.gz", "transcripts.csv"]:
            if (base_dir / candidate).exists():
                tx_file = base_dir / candidate
                break
        if tx_file is not None:
            transcripts = read_transcripts(tx_file, min_quality=min_transcript_quality)

    # --- Optionally load boundaries ---
    boundaries: list[CellBoundary] | None = None
    if load_boundaries:
        bd_file = None
        for candidate in ["cell_boundaries.csv.gz", "cell_boundaries.csv", "cell_boundaries.parquet"]:
            if (base_dir / candidate).exists():
                bd_file = base_dir / candidate
                break
        if bd_file is not None:
            boundaries = load_cell_boundaries(bd_file)

    dataset = XeniumDataset(
        expression=expression,
        coordinates=coordinates,
        cell_ids=cell_ids,
        gene_names=gene_names,
        gene_ids=gene_ids,
        transcripts=transcripts,
        cell_boundaries=boundaries,
        metadata={
            "source_dir": str(base_dir),
            "platform": "xenium",
            "min_transcript_quality": min_transcript_quality,
        },
    )

    logger.info(f"Loaded Xenium dataset: {dataset.n_cells} cells, {dataset.n_genes} genes")
    return dataset


def _find_col(headers: set[str], candidates: list[str]) -> str | None:
    """Find first matching column name from candidates."""
    for c in candidates:
        if c in headers:
            return c
    return None
