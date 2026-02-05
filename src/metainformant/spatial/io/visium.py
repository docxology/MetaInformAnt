"""10x Visium spatial transcriptomics data loading.

Provides functions to load and parse 10x Visium Spatial Gene Expression data,
including the filtered feature-barcode matrix, tissue positions, and H&E images.

Visium outputs typically include:
  - filtered_feature_bc_matrix/ (MEX or HDF5 format)
  - spatial/tissue_positions.csv (or tissue_positions_list.csv)
  - spatial/tissue_hires_image.png / tissue_lowres_image.png
  - spatial/scalefactors_json.json
"""

from __future__ import annotations

import csv
import json
import os
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
    from PIL import Image
except ImportError:
    Image = None  # type: ignore[assignment]

try:
    import pandas as pd
except ImportError:
    pd = None  # type: ignore[assignment]

try:
    import h5py
except ImportError:
    h5py = None  # type: ignore[assignment]


@dataclass
class TissuePosition:
    """A single Visium spot position on the tissue.

    Attributes:
        barcode: Spot barcode string.
        in_tissue: Whether the spot overlaps tissue (1) or not (0).
        array_row: Row index in the Visium array grid.
        array_col: Column index in the Visium array grid.
        pixel_row: Row pixel coordinate in the full-resolution image.
        pixel_col: Column pixel coordinate in the full-resolution image.
    """

    barcode: str
    in_tissue: bool
    array_row: int
    array_col: int
    pixel_row: float
    pixel_col: float


@dataclass
class SpatialDataset:
    """Unified spatial transcriptomics dataset.

    Attributes:
        expression: Gene expression matrix (spots x genes) as numpy array or scipy sparse.
        coordinates: Spot coordinates array of shape (n_spots, 2) [row, col in pixels].
        barcodes: List of spot barcodes.
        gene_names: List of gene names/symbols.
        gene_ids: List of gene IDs (e.g., Ensembl).
        tissue_positions: Full TissuePosition records per spot.
        image: Optional tissue image as numpy array (H, W, C).
        scale_factors: Dictionary of Visium scale factors.
        platform: Originating platform (visium, merfish, xenium).
        metadata: Additional metadata dictionary.
    """

    expression: Any  # np.ndarray or scipy sparse
    coordinates: Any  # np.ndarray (n_spots, 2)
    barcodes: list[str]
    gene_names: list[str]
    gene_ids: list[str] = field(default_factory=list)
    tissue_positions: list[TissuePosition] = field(default_factory=list)
    image: Any | None = None  # np.ndarray (H, W, C) or None
    scale_factors: dict[str, Any] = field(default_factory=dict)
    platform: str = "visium"
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def n_spots(self) -> int:
        """Number of spots/cells."""
        if np is not None and hasattr(self.expression, "shape"):
            return int(self.expression.shape[0])
        return len(self.barcodes)

    @property
    def n_genes(self) -> int:
        """Number of genes."""
        return len(self.gene_names)


def read_tissue_positions(positions_file: str | Path) -> list[TissuePosition]:
    """Parse a Visium tissue_positions.csv file.

    Supports both Space Ranger v1 (tissue_positions_list.csv, no header)
    and v2 (tissue_positions.csv, with header) formats.

    Args:
        positions_file: Path to tissue_positions.csv or tissue_positions_list.csv.

    Returns:
        List of TissuePosition records.

    Raises:
        FileNotFoundError: If the positions file does not exist.
        ValueError: If the file format is unrecognized.
    """
    positions_path = Path(positions_file)
    if not positions_path.exists():
        raise FileNotFoundError(f"Tissue positions file not found: {positions_path}")

    positions: list[TissuePosition] = []

    with open(positions_path, "r", newline="") as fh:
        # Peek at first line to detect header
        first_line = fh.readline().strip()
        fh.seek(0)

        has_header = False
        if first_line.startswith("barcode"):
            has_header = True

        reader = csv.reader(fh)
        if has_header:
            next(reader)  # skip header

        for row in reader:
            if len(row) < 6:
                continue
            barcode = row[0].strip()
            in_tissue = bool(int(row[1].strip()))
            array_row = int(row[2].strip())
            array_col = int(row[3].strip())
            pixel_row = float(row[4].strip())
            pixel_col = float(row[5].strip())

            positions.append(
                TissuePosition(
                    barcode=barcode,
                    in_tissue=in_tissue,
                    array_row=array_row,
                    array_col=array_col,
                    pixel_row=pixel_row,
                    pixel_col=pixel_col,
                )
            )

    logger.info(f"Read {len(positions)} tissue positions from {positions_path.name}")
    return positions


def read_spatial_image(image_path: str | Path) -> Any:
    """Load an H&E tissue image as a numpy array.

    Args:
        image_path: Path to the tissue image (PNG/JPEG).

    Returns:
        Numpy array of shape (height, width, channels) with values in [0, 255] uint8.

    Raises:
        FileNotFoundError: If image file does not exist.
        ImportError: If Pillow (PIL) is not installed.
    """
    image_path = Path(image_path)
    if not image_path.exists():
        raise FileNotFoundError(f"Tissue image not found: {image_path}")

    if Image is None:
        raise ImportError("Pillow (PIL) is required to load tissue images: uv pip install Pillow")
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")

    img = Image.open(image_path)
    img_array = np.array(img)
    logger.info(f"Loaded tissue image {image_path.name}: shape={img_array.shape}, dtype={img_array.dtype}")
    return img_array


def filter_tissue_spots(
    positions: list[TissuePosition],
    in_tissue_only: bool = True,
) -> list[TissuePosition]:
    """Filter tissue position records to keep only spots overlapping tissue.

    Args:
        positions: List of TissuePosition records.
        in_tissue_only: If True, keep only spots with in_tissue=True.
            If False, return all positions unchanged.

    Returns:
        Filtered list of TissuePosition records.
    """
    if not in_tissue_only:
        return positions

    filtered = [p for p in positions if p.in_tissue]
    logger.info(f"Filtered tissue spots: {len(filtered)}/{len(positions)} in tissue")
    return filtered


def _read_mex_matrix(matrix_dir: Path) -> tuple[Any, list[str], list[str], list[str]]:
    """Read a Market Exchange (MEX) format matrix directory.

    Expected files: matrix.mtx.gz or matrix.mtx, features.tsv.gz or genes.tsv.gz,
    barcodes.tsv.gz or barcodes.tsv.

    Returns:
        Tuple of (sparse_matrix, barcodes, gene_names, gene_ids).
    """
    if sp_sparse is None or mmread is None:
        raise ImportError("scipy is required to read MEX matrices: uv pip install scipy")
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")

    import gzip

    # Find matrix file
    mtx_file = None
    for candidate in ["matrix.mtx.gz", "matrix.mtx"]:
        if (matrix_dir / candidate).exists():
            mtx_file = matrix_dir / candidate
            break
    if mtx_file is None:
        raise FileNotFoundError(f"No matrix.mtx(.gz) found in {matrix_dir}")

    # Read sparse matrix (genes x barcodes) then transpose to (barcodes x genes)
    mat = mmread(str(mtx_file))
    mat = sp_sparse.csc_matrix(mat).T.tocsr()

    # Read barcodes
    barcodes_file = None
    for candidate in ["barcodes.tsv.gz", "barcodes.tsv"]:
        if (matrix_dir / candidate).exists():
            barcodes_file = matrix_dir / candidate
            break
    if barcodes_file is None:
        raise FileNotFoundError(f"No barcodes.tsv(.gz) found in {matrix_dir}")

    barcodes: list[str] = []
    opener = gzip.open if str(barcodes_file).endswith(".gz") else open
    with opener(barcodes_file, "rt") as fh:  # type: ignore[call-overload]
        for line in fh:
            barcodes.append(line.strip())

    # Read features/genes
    features_file = None
    for candidate in ["features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"]:
        if (matrix_dir / candidate).exists():
            features_file = matrix_dir / candidate
            break
    if features_file is None:
        raise FileNotFoundError(f"No features/genes TSV found in {matrix_dir}")

    gene_ids: list[str] = []
    gene_names: list[str] = []
    opener = gzip.open if str(features_file).endswith(".gz") else open
    with opener(features_file, "rt") as fh:  # type: ignore[call-overload]
        for line in fh:
            parts = line.strip().split("\t")
            gene_ids.append(parts[0] if len(parts) > 0 else "")
            gene_names.append(parts[1] if len(parts) > 1 else parts[0])

    logger.info(f"Read MEX matrix: {mat.shape[0]} barcodes x {mat.shape[1]} genes")
    return mat, barcodes, gene_names, gene_ids


def _read_h5_matrix(h5_path: Path) -> tuple[Any, list[str], list[str], list[str]]:
    """Read a 10x HDF5 filtered feature-barcode matrix.

    Returns:
        Tuple of (sparse_matrix, barcodes, gene_names, gene_ids).
    """
    if h5py is None:
        raise ImportError("h5py is required to read HDF5 matrices: uv pip install h5py")
    if sp_sparse is None:
        raise ImportError("scipy is required: uv pip install scipy")
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")

    with h5py.File(str(h5_path), "r") as f:
        # 10x HDF5 v3 format
        if "matrix" in f:
            group = f["matrix"]
            data = np.array(group["data"])
            indices = np.array(group["indices"])
            indptr = np.array(group["indptr"])
            shape = tuple(group["shape"])
            barcodes = [b.decode("utf-8") if isinstance(b, bytes) else b for b in group["barcodes"]]

            features_group = group["features"]
            gene_ids_raw = features_group["id"]
            gene_names_raw = features_group["name"]
            gene_ids = [g.decode("utf-8") if isinstance(g, bytes) else g for g in gene_ids_raw]
            gene_names = [g.decode("utf-8") if isinstance(g, bytes) else g for g in gene_names_raw]

            # Matrix is (genes x barcodes) in CSC format, transpose to (barcodes x genes)
            mat = sp_sparse.csc_matrix((data, indices, indptr), shape=shape).T.tocsr()
        else:
            raise ValueError(f"Unrecognized HDF5 format in {h5_path}")

    logger.info(f"Read HDF5 matrix: {mat.shape[0]} barcodes x {mat.shape[1]} genes")
    return mat, barcodes, gene_names, gene_ids


def _read_scale_factors(spatial_dir: Path) -> dict[str, Any]:
    """Read Visium scale factors JSON."""
    sf_path = spatial_dir / "scalefactors_json.json"
    if not sf_path.exists():
        return {}

    with open(sf_path, "r") as fh:
        scale_factors = json.load(fh)

    logger.info(f"Read scale factors: {list(scale_factors.keys())}")
    return scale_factors


def load_visium(
    path: str | Path,
    *,
    in_tissue_only: bool = True,
    load_image: bool = True,
    image_resolution: str = "hires",
) -> SpatialDataset:
    """Load a 10x Visium Spatial Gene Expression dataset.

    Expects the standard Space Ranger output directory structure:
        path/
            filtered_feature_bc_matrix/ (or filtered_feature_bc_matrix.h5)
            spatial/
                tissue_positions_list.csv (or tissue_positions.csv)
                tissue_hires_image.png
                tissue_lowres_image.png
                scalefactors_json.json

    Args:
        path: Path to the Space Ranger output directory.
        in_tissue_only: If True, keep only spots overlapping tissue.
        load_image: If True, load the tissue image.
        image_resolution: Which resolution image to load ("hires" or "lowres").

    Returns:
        SpatialDataset with expression matrix, coordinates, and optionally image.

    Raises:
        FileNotFoundError: If required files are missing.
    """
    if np is None:
        raise ImportError("NumPy is required for Visium loading: uv pip install numpy")

    base_dir = Path(path)
    if not base_dir.exists():
        raise FileNotFoundError(f"Visium directory not found: {base_dir}")

    # --- Load expression matrix ---
    h5_path = base_dir / "filtered_feature_bc_matrix.h5"
    mex_dir = base_dir / "filtered_feature_bc_matrix"

    if h5_path.exists():
        expression, barcodes, gene_names, gene_ids = _read_h5_matrix(h5_path)
    elif mex_dir.exists():
        expression, barcodes, gene_names, gene_ids = _read_mex_matrix(mex_dir)
    else:
        raise FileNotFoundError(
            f"No expression matrix found in {base_dir}. "
            "Expected filtered_feature_bc_matrix.h5 or filtered_feature_bc_matrix/ directory."
        )

    # --- Load spatial positions ---
    spatial_dir = base_dir / "spatial"
    positions_file = None
    for candidate in ["tissue_positions.csv", "tissue_positions_list.csv"]:
        if (spatial_dir / candidate).exists():
            positions_file = spatial_dir / candidate
            break

    if positions_file is None:
        raise FileNotFoundError(f"No tissue_positions file found in {spatial_dir}")

    all_positions = read_tissue_positions(positions_file)
    positions = filter_tissue_spots(all_positions, in_tissue_only=in_tissue_only)

    # Build barcode-to-position mapping
    barcode_to_pos = {p.barcode: p for p in positions}

    # Filter expression matrix to match tissue positions
    keep_indices: list[int] = []
    keep_barcodes: list[str] = []
    keep_positions: list[TissuePosition] = []

    for idx, bc in enumerate(barcodes):
        if bc in barcode_to_pos:
            keep_indices.append(idx)
            keep_barcodes.append(bc)
            keep_positions.append(barcode_to_pos[bc])

    if not keep_indices:
        logger.warning("No barcodes match between expression matrix and tissue positions")
        coordinates = np.empty((0, 2), dtype=np.float64)
    else:
        expression = expression[keep_indices, :]
        coordinates = np.array([[p.pixel_row, p.pixel_col] for p in keep_positions], dtype=np.float64)

    # --- Load scale factors ---
    scale_factors = _read_scale_factors(spatial_dir)

    # --- Load tissue image ---
    tissue_image = None
    if load_image:
        image_name = f"tissue_{image_resolution}_image.png"
        image_path = spatial_dir / image_name
        if image_path.exists():
            try:
                tissue_image = read_spatial_image(image_path)
            except ImportError:
                logger.warning("Pillow not available; skipping tissue image")

    dataset = SpatialDataset(
        expression=expression,
        coordinates=coordinates,
        barcodes=keep_barcodes,
        gene_names=gene_names,
        gene_ids=gene_ids,
        tissue_positions=keep_positions,
        image=tissue_image,
        scale_factors=scale_factors,
        platform="visium",
        metadata={
            "source_dir": str(base_dir),
            "in_tissue_only": in_tissue_only,
            "image_resolution": image_resolution,
        },
    )

    logger.info(f"Loaded Visium dataset: {dataset.n_spots} spots, {dataset.n_genes} genes")
    return dataset


def create_spatial_dataset(
    matrix: Any,
    positions: list[TissuePosition],
    image: Any | None = None,
    gene_names: list[str] | None = None,
    gene_ids: list[str] | None = None,
    scale_factors: dict[str, Any] | None = None,
    platform: str = "visium",
) -> SpatialDataset:
    """Create a unified SpatialDataset from components.

    Args:
        matrix: Expression matrix (spots x genes) as numpy array or scipy sparse.
        positions: List of TissuePosition records (one per row in matrix).
        image: Optional tissue image array.
        gene_names: List of gene names. Defaults to G0, G1, ...
        gene_ids: List of gene IDs.
        scale_factors: Optional scale factor dictionary.
        platform: Originating platform name.

    Returns:
        SpatialDataset instance.
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")

    n_spots = matrix.shape[0] if hasattr(matrix, "shape") else len(positions)
    n_genes = matrix.shape[1] if hasattr(matrix, "shape") else 0

    barcodes = [p.barcode for p in positions]
    coordinates = np.array([[p.pixel_row, p.pixel_col] for p in positions], dtype=np.float64)

    if gene_names is None:
        gene_names = [f"G{i}" for i in range(n_genes)]
    if gene_ids is None:
        gene_ids = []

    return SpatialDataset(
        expression=matrix,
        coordinates=coordinates,
        barcodes=barcodes,
        gene_names=gene_names,
        gene_ids=gene_ids,
        tissue_positions=positions,
        image=image,
        scale_factors=scale_factors or {},
        platform=platform,
    )
