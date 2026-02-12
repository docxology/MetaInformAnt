# Spatial I/O

The spatial I/O module provides data loading functions for the three major spatial transcriptomics platforms: 10x Visium, MERFISH, and 10x Xenium. Each platform loader returns structured data objects that feed into the downstream analysis pipeline.

## Key Concepts

### Platform Overview

- **10x Visium**: Spot-based (55um diameter, ~1-10 cells per spot). Outputs include filtered feature-barcode matrices (MEX or HDF5), tissue position lists, and H&E images.
- **MERFISH**: Subcellular resolution with individual transcript coordinates. Outputs include cell-by-gene matrices, cell metadata with 3D coordinates, and detected transcript positions.
- **10x Xenium**: Subcellular resolution with per-transcript coordinates, cell segmentation boundaries, and feature matrices. Outputs include transcripts.csv.gz, cells.csv.gz, and cell_boundaries files.

### Unified SpatialDataset

All platform loaders produce a `SpatialDataset` object with a consistent interface for downstream analysis, regardless of the originating platform.

## Data Structures

### SpatialDataset

```python
@dataclass
class SpatialDataset:
    expression: Any         # Gene expression matrix (spots/cells x genes)
    coordinates: Any        # Coordinates array (n, 2) [row, col in pixels/microns]
    barcodes: list[str]     # Spot/cell barcodes
    gene_names: list[str]   # Gene names/symbols
    gene_ids: list[str]     # Gene IDs (e.g., Ensembl)
    tissue_positions: list[TissuePosition]  # Visium-specific positions
    image: Any              # Tissue image as numpy array (H, W, C)
    scale_factors: dict     # Platform-specific scale factors
    platform: str           # "visium", "merfish", or "xenium"
    metadata: dict          # Additional metadata
```

### TissuePosition (Visium)

```python
@dataclass
class TissuePosition:
    barcode: str        # Spot barcode
    in_tissue: bool     # Whether spot overlaps tissue
    array_row: int      # Row index in Visium array grid
    array_col: int      # Column index in Visium array grid
    pixel_row: float    # Row pixel coordinate
    pixel_col: float    # Column pixel coordinate
```

### MERFISHDataset, CellMetadata, TranscriptSpot

```python
@dataclass
class CellMetadata:
    cell_id: str; x: float; y: float; z: float
    volume: float; fov: int; extra: dict[str, Any]

@dataclass
class TranscriptSpot:
    gene: str; x: float; y: float; z: float
    cell_id: str; quality: float
```

### XeniumDataset, XeniumTranscript, CellBoundary

```python
@dataclass
class XeniumTranscript:
    transcript_id: str; gene: str; x: float; y: float; z: float
    cell_id: str; quality_value: float; overlaps_nucleus: bool

@dataclass
class CellBoundary:
    cell_id: str
    vertices: list[tuple[float, float]]  # Polygon boundary
```

## Function Reference

### Visium

```python
def load_visium(spaceranger_dir: str | Path) -> SpatialDataset
def read_tissue_positions(positions_file: str | Path) -> list[TissuePosition]
def read_spatial_image(image_path: str | Path) -> Any  # numpy array
def filter_tissue_spots(dataset: SpatialDataset) -> SpatialDataset
def create_spatial_dataset(...) -> SpatialDataset
```

### MERFISH

```python
def load_merfish(data_dir: str | Path) -> MERFISHDataset
def load_transcript_spots(transcripts_file: str | Path) -> list[TranscriptSpot]
def parse_cell_metadata(metadata_file: str | Path) -> list[CellMetadata]
def aggregate_to_cells(transcripts: list[TranscriptSpot]) -> Any
```

### Xenium

```python
def load_xenium(xenium_dir: str | Path) -> XeniumDataset
def read_cell_features(feature_dir: str | Path) -> tuple[Any, list[str]]
def read_transcripts(transcripts_file: str | Path) -> list[XeniumTranscript]
def load_cell_boundaries(boundaries_file: str | Path) -> list[CellBoundary]
```

## Usage Examples

```python
from metainformant.spatial import io

# Load 10x Visium data
dataset = io.load_visium("path/to/spaceranger_output/")
print(f"Spots: {len(dataset.barcodes)}, Genes: {len(dataset.gene_names)}")

# Filter to tissue-overlapping spots only
filtered = io.filter_tissue_spots(dataset)

# Load MERFISH data
merfish = io.load_merfish("path/to/merfish_output/")
transcripts = io.load_transcript_spots("path/to/detected_transcripts.csv")

# Load Xenium data
xenium = io.load_xenium("path/to/xenium_output/")
boundaries = io.load_cell_boundaries("path/to/cell_boundaries.csv.gz")
```

## Configuration

- **Environment prefix**: `SPATIAL_`
- **Optional dependencies**: numpy, scipy, h5py, pandas, Pillow
- Supports MEX (Market Exchange Format) and HDF5 feature matrices
- Gzip-compressed CSV files are handled transparently

## Related Modules

- `spatial.analysis.clustering` -- Spatial clustering on loaded datasets
- `spatial.analysis.autocorrelation` -- Spatial statistics on expression data
- `spatial.deconvolution` -- Cell type deconvolution for Visium spots
- `spatial.visualization` -- Plotting functions for spatial data
