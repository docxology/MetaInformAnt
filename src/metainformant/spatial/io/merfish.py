"""MERFISH spatial transcriptomics data loading.

MERFISH (Multiplexed Error-Robust Fluorescence In Situ Hybridization) produces
subcellular-resolution spatial data with individual transcript coordinates and
cell segmentation boundaries.

Typical MERFISH outputs:
  - cell_by_gene.csv: Cell-level gene expression matrix
  - cell_metadata.csv: Cell positions, volumes, segmentation IDs
  - detected_transcripts.csv: Per-transcript x/y/z coordinates
"""

from __future__ import annotations

import csv
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
    import pandas as pd
except ImportError:
    pd = None  # type: ignore[assignment]


@dataclass
class CellMetadata:
    """Metadata for a single MERFISH cell.

    Attributes:
        cell_id: Unique cell identifier.
        x: Cell centroid x-coordinate.
        y: Cell centroid y-coordinate.
        z: Cell centroid z-coordinate (layer/z-plane).
        volume: Cell volume (in cubic microns).
        fov: Field of view index.
        extra: Additional metadata key-value pairs.
    """

    cell_id: str
    x: float
    y: float
    z: float = 0.0
    volume: float = 0.0
    fov: int = 0
    extra: dict[str, Any] = field(default_factory=dict)


@dataclass
class TranscriptSpot:
    """A single detected transcript location.

    Attributes:
        gene: Gene name.
        x: Transcript x-coordinate.
        y: Transcript y-coordinate.
        z: Transcript z-coordinate.
        cell_id: Cell assignment (if segmented), or empty.
        quality: Detection quality/confidence score.
    """

    gene: str
    x: float
    y: float
    z: float = 0.0
    cell_id: str = ""
    quality: float = 1.0


@dataclass
class MERFISHDataset:
    """Complete MERFISH spatial dataset.

    Attributes:
        expression: Cell-by-gene expression matrix (numpy array, cells x genes).
        coordinates: Cell centroid coordinates array of shape (n_cells, 2).
        cell_ids: List of cell identifier strings.
        gene_names: List of gene names.
        cell_metadata: List of CellMetadata records.
        transcript_spots: Optional list of individual TranscriptSpot records.
        metadata: Additional dataset metadata.
    """

    expression: Any  # np.ndarray (n_cells, n_genes)
    coordinates: Any  # np.ndarray (n_cells, 2)
    cell_ids: list[str]
    gene_names: list[str]
    cell_metadata: list[CellMetadata] = field(default_factory=list)
    transcript_spots: list[TranscriptSpot] | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def n_cells(self) -> int:
        """Number of cells."""
        return len(self.cell_ids)

    @property
    def n_genes(self) -> int:
        """Number of genes."""
        return len(self.gene_names)


def parse_cell_metadata(metadata_file: str | Path) -> list[CellMetadata]:
    """Parse MERFISH cell metadata CSV file.

    Reads cell centroid positions, volumes, FOVs, and any additional columns.
    The file should have columns including: cell_id (or CellID), center_x, center_y
    (or x, y), and optionally center_z, volume, fov.

    Args:
        metadata_file: Path to cell_metadata.csv.

    Returns:
        List of CellMetadata records.

    Raises:
        FileNotFoundError: If the metadata file does not exist.
        ValueError: If required columns are missing.
    """
    metadata_path = Path(metadata_file)
    if not metadata_path.exists():
        raise FileNotFoundError(f"Cell metadata file not found: {metadata_path}")

    cells: list[CellMetadata] = []

    with open(metadata_path, "r", newline="") as fh:
        reader = csv.DictReader(fh)
        headers = set(reader.fieldnames or [])

        # Detect column naming conventions
        id_col = _find_column(headers, ["cell_id", "CellID", "cellID", "cell", "Cell"])
        x_col = _find_column(headers, ["center_x", "x", "X", "centroid_x", "global_x"])
        y_col = _find_column(headers, ["center_y", "y", "Y", "centroid_y", "global_y"])
        z_col = _find_column(headers, ["center_z", "z", "Z", "centroid_z", "global_z"])
        vol_col = _find_column(headers, ["volume", "Volume", "cell_volume"])
        fov_col = _find_column(headers, ["fov", "FOV", "fov_id"])

        if id_col is None:
            raise ValueError(f"No cell ID column found in {metadata_path}. " f"Available columns: {sorted(headers)}")
        if x_col is None or y_col is None:
            raise ValueError(
                f"No coordinate columns (x, y) found in {metadata_path}. " f"Available columns: {sorted(headers)}"
            )

        known_cols = {id_col, x_col, y_col}
        if z_col:
            known_cols.add(z_col)
        if vol_col:
            known_cols.add(vol_col)
        if fov_col:
            known_cols.add(fov_col)

        for row in reader:
            extra = {k: v for k, v in row.items() if k not in known_cols}
            cells.append(
                CellMetadata(
                    cell_id=str(row[id_col]),
                    x=float(row[x_col]),
                    y=float(row[y_col]),
                    z=float(row[z_col]) if z_col and row.get(z_col) else 0.0,
                    volume=float(row[vol_col]) if vol_col and row.get(vol_col) else 0.0,
                    fov=int(float(row[fov_col])) if fov_col and row.get(fov_col) else 0,
                    extra=extra,
                )
            )

    logger.info(f"Parsed {len(cells)} cells from {metadata_path.name}")
    return cells


def load_transcript_spots(spots_file: str | Path) -> list[TranscriptSpot]:
    """Load individual transcript spot coordinates from a MERFISH detected transcripts file.

    The file should have columns: gene, x (or global_x), y (or global_y),
    and optionally z, cell_id, quality.

    Args:
        spots_file: Path to detected_transcripts.csv.

    Returns:
        List of TranscriptSpot records.

    Raises:
        FileNotFoundError: If the spots file does not exist.
    """
    spots_path = Path(spots_file)
    if not spots_path.exists():
        raise FileNotFoundError(f"Transcript spots file not found: {spots_path}")

    spots: list[TranscriptSpot] = []

    with open(spots_path, "r", newline="") as fh:
        reader = csv.DictReader(fh)
        headers = set(reader.fieldnames or [])

        gene_col = _find_column(headers, ["gene", "Gene", "gene_name", "target", "feature_name"])
        x_col = _find_column(headers, ["x", "global_x", "x_location", "X"])
        y_col = _find_column(headers, ["y", "global_y", "y_location", "Y"])
        z_col = _find_column(headers, ["z", "global_z", "z_location", "Z"])
        cell_col = _find_column(headers, ["cell_id", "CellID", "cell", "cell_assignment"])
        qual_col = _find_column(headers, ["quality", "qv", "score", "confidence"])

        if gene_col is None or x_col is None or y_col is None:
            raise ValueError(
                f"Required columns (gene, x, y) not found in {spots_path}. " f"Available: {sorted(headers)}"
            )

        for row in reader:
            spots.append(
                TranscriptSpot(
                    gene=str(row[gene_col]),
                    x=float(row[x_col]),
                    y=float(row[y_col]),
                    z=float(row[z_col]) if z_col and row.get(z_col) else 0.0,
                    cell_id=str(row[cell_col]) if cell_col and row.get(cell_col) else "",
                    quality=float(row[qual_col]) if qual_col and row.get(qual_col) else 1.0,
                )
            )

    logger.info(f"Loaded {len(spots)} transcript spots from {spots_path.name}")
    return spots


def aggregate_to_cells(
    transcript_spots: list[TranscriptSpot],
    cell_boundaries: list[CellMetadata],
    *,
    unassigned_label: str = "",
) -> tuple[Any, list[str], list[str]]:
    """Aggregate individual transcript spots to cell-level expression counts.

    For each cell, counts the number of transcripts per gene that fall within
    the cell's assignment. If transcripts have pre-assigned cell_id, uses that;
    otherwise uses nearest-centroid assignment.

    Args:
        transcript_spots: List of TranscriptSpot records.
        cell_boundaries: List of CellMetadata records (used for cell IDs and positions).
        unassigned_label: Label for unassigned transcripts (excluded from output).

    Returns:
        Tuple of (expression_matrix, cell_ids, gene_names) where expression_matrix
        is a numpy array of shape (n_cells, n_genes) with integer counts.

    Raises:
        ImportError: If numpy is not installed.
    """
    if np is None:
        raise ImportError("NumPy is required for transcript aggregation: uv pip install numpy")

    # Build cell ID list and gene set
    cell_id_list = [c.cell_id for c in cell_boundaries]
    cell_id_to_idx = {cid: idx for idx, cid in enumerate(cell_id_list)}

    gene_set: set[str] = set()
    for spot in transcript_spots:
        gene_set.add(spot.gene)
    gene_list = sorted(gene_set)
    gene_to_idx = {g: idx for idx, g in enumerate(gene_list)}

    n_cells = len(cell_id_list)
    n_genes = len(gene_list)
    counts = np.zeros((n_cells, n_genes), dtype=np.int32)

    # Check if transcripts have pre-assigned cell IDs
    has_assignments = any(spot.cell_id and spot.cell_id != unassigned_label for spot in transcript_spots)

    if has_assignments:
        # Use pre-assigned cell IDs
        for spot in transcript_spots:
            if spot.cell_id == unassigned_label or not spot.cell_id:
                continue
            cell_idx = cell_id_to_idx.get(spot.cell_id)
            if cell_idx is not None:
                gene_idx = gene_to_idx[spot.gene]
                counts[cell_idx, gene_idx] += 1
    else:
        # Nearest-centroid assignment
        cell_coords = np.array([[c.x, c.y] for c in cell_boundaries], dtype=np.float64)
        for spot in transcript_spots:
            spot_coord = np.array([spot.x, spot.y], dtype=np.float64)
            distances = np.sqrt(np.sum((cell_coords - spot_coord) ** 2, axis=1))
            nearest_idx = int(np.argmin(distances))
            gene_idx = gene_to_idx[spot.gene]
            counts[nearest_idx, gene_idx] += 1

    assigned_count = int(np.sum(counts))
    logger.info(f"Aggregated {assigned_count} transcripts across {n_cells} cells and {n_genes} genes")
    return counts, cell_id_list, gene_list


def load_merfish(
    path: str | Path,
    *,
    load_transcripts: bool = False,
) -> MERFISHDataset:
    """Load a complete MERFISH dataset from a directory.

    Expects:
        path/
            cell_by_gene.csv: Cell-level expression matrix (cells as rows, genes as columns).
            cell_metadata.csv: Cell positions and metadata.
            detected_transcripts.csv (optional): Per-transcript coordinates.

    Args:
        path: Path to the MERFISH output directory.
        load_transcripts: If True, also load individual transcript spots.

    Returns:
        MERFISHDataset with expression, coordinates, and metadata.

    Raises:
        FileNotFoundError: If required files are missing.
        ImportError: If numpy is not installed.
    """
    if np is None:
        raise ImportError("NumPy is required for MERFISH loading: uv pip install numpy")

    base_dir = Path(path)
    if not base_dir.exists():
        raise FileNotFoundError(f"MERFISH directory not found: {base_dir}")

    # --- Load cell-by-gene matrix ---
    expr_file = None
    for candidate in ["cell_by_gene.csv", "cell_by_gene.tsv", "counts.csv"]:
        if (base_dir / candidate).exists():
            expr_file = base_dir / candidate
            break

    if expr_file is None:
        raise FileNotFoundError(f"No cell-by-gene expression file found in {base_dir}")

    # Parse expression CSV: first column is cell ID, rest are genes
    cell_ids: list[str] = []
    gene_names: list[str] = []
    expression_rows: list[list[float]] = []

    delimiter = "\t" if str(expr_file).endswith(".tsv") else ","

    with open(expr_file, "r", newline="") as fh:
        reader = csv.reader(fh, delimiter=delimiter)
        header = next(reader)
        # First column is cell ID; remaining are gene names
        gene_names = header[1:]

        for row in reader:
            if len(row) < 2:
                continue
            cell_ids.append(row[0])
            expression_rows.append([float(v) if v else 0.0 for v in row[1:]])

    expression = np.array(expression_rows, dtype=np.float32)

    # --- Load cell metadata ---
    metadata_file = None
    for candidate in ["cell_metadata.csv", "cell_metadata.tsv"]:
        if (base_dir / candidate).exists():
            metadata_file = base_dir / candidate
            break

    cell_metadata: list[CellMetadata] = []
    if metadata_file is not None:
        cell_metadata = parse_cell_metadata(metadata_file)
        # Build coordinate array from metadata
        meta_id_to_meta = {m.cell_id: m for m in cell_metadata}
        coords_list: list[list[float]] = []
        for cid in cell_ids:
            meta = meta_id_to_meta.get(cid)
            if meta:
                coords_list.append([meta.x, meta.y])
            else:
                coords_list.append([0.0, 0.0])
        coordinates = np.array(coords_list, dtype=np.float64)
    else:
        logger.warning("No cell_metadata.csv found; coordinates will be zeros")
        coordinates = np.zeros((len(cell_ids), 2), dtype=np.float64)

    # --- Optionally load transcript spots ---
    transcript_spots: list[TranscriptSpot] | None = None
    if load_transcripts:
        transcripts_file = None
        for candidate in ["detected_transcripts.csv", "transcripts.csv"]:
            if (base_dir / candidate).exists():
                transcripts_file = base_dir / candidate
                break
        if transcripts_file is not None:
            transcript_spots = load_transcript_spots(transcripts_file)
        else:
            logger.warning("No detected_transcripts file found; skipping transcript spots")

    dataset = MERFISHDataset(
        expression=expression,
        coordinates=coordinates,
        cell_ids=cell_ids,
        gene_names=gene_names,
        cell_metadata=cell_metadata,
        transcript_spots=transcript_spots,
        metadata={
            "source_dir": str(base_dir),
            "platform": "merfish",
        },
    )

    logger.info(f"Loaded MERFISH dataset: {dataset.n_cells} cells, {dataset.n_genes} genes")
    return dataset


def _find_column(headers: set[str], candidates: list[str]) -> str | None:
    """Find the first matching column name from a list of candidates."""
    for c in candidates:
        if c in headers:
            return c
    return None
