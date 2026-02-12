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
| `merfish.load_merfish()` | Load complete MERFISH dataset |
| `merfish.parse_cell_metadata()` | Parse cell centroid positions and metadata |
| `merfish.load_transcript_spots()` | Load per-transcript x/y/z coordinates |
| `merfish.aggregate_to_cells()` | Aggregate transcript spots to cell-level counts |
| `xenium.load_xenium()` | Load complete Xenium dataset |
| `xenium.read_transcripts()` | Load per-transcript coordinates |
| `xenium.load_cell_boundaries()` | Load cell segmentation boundaries |

## Usage

```python
from metainformant.spatial.io import visium, merfish, xenium

vis_data = visium.load_visium("path/to/visium_output/")
mer_data = merfish.load_merfish("path/to/merfish_output/")
xen_data = xenium.load_xenium("path/to/xenium_output/")
```
