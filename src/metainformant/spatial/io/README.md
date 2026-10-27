# IO

Spatial transcriptomics data loading for 10x Visium, MERFISH, and 10x Xenium platforms, handling expression matrices, tissue positions, transcript coordinates, and cell boundaries.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports visium, merfish, xenium submodules |
| `visium.py` | 10x Visium loader (MEX/HDF5 matrix, tissue positions, scale factors) |
| `merfish.py` | MERFISH loader (cell-by-gene, transcript spots, cell metadata) |
| `xenium.py` | 10x Xenium loader (feature matrix, transcripts, cell boundaries) |

## Key Functions

| Function | Description |
|----------|-------------|
| `visium.load_visium()` | Load complete Visium dataset (matrix + positions + image) |
| `visium.read_tissue_positions()` | Parse tissue position coordinate file |
| `visium.read_spatial_image()` | Load H&E tissue image |
| `visium.create_spatial_dataset()` | Create unified spatial dataset from components |
| `merfish.load_merfish()` | Load complete MERFISH dataset; `cell_metadata.csv` or `.tsv` is required (the delimiter follows the file extension: tab for `.tsv`, comma otherwise) |
| `merfish.parse_cell_metadata()` | Parse cell centroid positions and metadata (`.csv` or `.tsv`) |
| `merfish.load_transcript_spots()` | Load per-transcript x/y/z coordinates (`.csv` or `.tsv`) |
| `merfish.aggregate_to_cells()` | Aggregate transcript spots to cell-level counts |
| `xenium.load_xenium()` | Load complete Xenium dataset; `cells.csv(.gz)` is required |
| `xenium.read_transcripts()` | Load per-transcript coordinates |
| `xenium.load_cell_boundaries()` | Load cell segmentation boundaries |

## Return values and coordinate order

`load_visium` returns a `SpatialDataset`, `load_xenium` a `XeniumDataset` and `load_merfish` a `MERFISHDataset`. Each has `expression`, `coordinates` (`n × 2`) and `gene_names`, plus a `coordinate_order` class constant that names the column convention of `coordinates`: `"row_col"` for Visium (`[pixel_row, pixel_col]`) and `"xy"` for Xenium and MERFISH (`[x, y]`). The plotting functions in `visualization.plots` read it.

## Errors

`load_xenium` and `load_merfish` never give a cell a placeholder `[0, 0]` position. They raise `FileNotFoundError` when `cells.csv(.gz)` or `cell_metadata.csv`/`.tsv` (or the matrix) is missing, and `ValueError` when the id or coordinate columns are missing or a cell of the matrix has no row. `load_visium` raises `FileNotFoundError` for a missing directory, matrix or tissue-positions file, and `ValueError` when no barcode of the matrix has a tissue position (after the `in_tissue_only` filter); the message gives the barcode counts and example barcodes from both files. On success its expression matrix, barcodes, tissue positions and coordinates are filtered to the same spots, in matrix order.

## Usage

```python
from metainformant.spatial.io import visium, merfish, xenium

vis_data = visium.load_visium("path/to/visium_output/")
mer_data = merfish.load_merfish("path/to/merfish_output/")
xen_data = xenium.load_xenium("path/to/xenium_output/")
print(vis_data.coordinate_order, xen_data.coordinate_order)   # row_col xy
```
