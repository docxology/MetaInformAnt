# Data Loading: Spatial

Loaders for 10x Visium, 10x Xenium and MERFISH output directories. Each platform has its own module and its own loader function; there is no dispatcher that detects the platform, and no AnnData is created.

## Common return value

Each loader returns a dataclass:

| Class | Module | `expression` | `coordinates` (float64, `n × 2`) | `coordinate_order` | Identifiers |
|-------|--------|--------------|-----------------------------------|--------------------|-------------|
| `SpatialDataset` | `spatial.io.visium` | spots × genes, SciPy CSR | `[pixel_row, pixel_col]` in full-resolution image pixels | `"row_col"` | `barcodes` |
| `XeniumDataset` | `spatial.io.xenium` | cells × genes, SciPy CSR | `[x_centroid, y_centroid]` | `"xy"` | `cell_ids` |
| `MERFISHDataset` | `spatial.io.merfish` | cells × genes, dense `float32` array | `[center_x, center_y]` | `"xy"` | `cell_ids` |

All three have `gene_names` and a `metadata` dictionary (it includes `source_dir`); the Visium and Xenium datasets also keep `gene_ids` (the first column of `features.tsv`). `SpatialDataset` adds `tissue_positions`, `image`, `scale_factors` and `platform`, `XeniumDataset` adds `transcripts` and `cell_boundaries`, and `MERFISHDataset` adds `cell_metadata` and `transcript_spots`. `SpatialDataset` has `n_spots` and `n_genes`; the other two have `n_cells` and `n_genes`.

Coordinates are taken as stored in the files, with no unit conversion.

### Coordinate order

`coordinate_order` is a class constant that names the column convention of `coordinates`:

- `"row_col"` (Visium): column 0 is the vertical image axis (the pixel row), column 1 the horizontal one.
- `"xy"` (Xenium and MERFISH): column 0 is x, column 1 is y.

The plotting functions in [Visualization](visualization.md) read it: `plot_gene_expression_map` and `plot_cell_type_map` take it from the dataset object, while the functions that take a bare coordinate array default to `"row_col"`, so pass `coordinate_order=dataset.coordinate_order` (or `"xy"`) for Xenium and MERFISH arrays. With the wrong order the tissue appears transposed.

## 10x Visium (`load_visium`)

```python
from metainformant.spatial.io.visium import load_visium

dataset = load_visium(
    'data/visium_sample/',        # Space Ranger output directory
    in_tissue_only=True,          # keep only spots under tissue
    load_image=True,              # read spatial/tissue_<resolution>_image.png
    image_resolution='hires',     # or 'lowres'
)
print(dataset.n_spots, dataset.n_genes, dataset.coordinate_order)
```

**Files read** from `path`:

```
visium_sample/
├── filtered_feature_bc_matrix.h5      (or the filtered_feature_bc_matrix/ MEX directory)
└── spatial/
    ├── tissue_positions.csv           (or tissue_positions_list.csv)
    ├── tissue_hires_image.png         (optional)
    ├── tissue_lowres_image.png        (optional)
    └── scalefactors_json.json         (optional)
```

- A MEX directory holds `matrix.mtx(.gz)`, `barcodes.tsv(.gz)` and `features.tsv(.gz)` (or `genes.tsv(.gz)`). The `.h5` form needs `h5py`.
- `tissue_positions.csv` may or may not have a header row (Space Ranger v2 and v1 layouts).
- Matrix barcodes are matched to the position records: barcodes without a position record are dropped, and with `in_tissue_only=True` the spots outside the tissue are dropped too. The expression matrix, `barcodes`, `tissue_positions` and `coordinates` are filtered to the same spots, in the order of the matrix.
- If no barcode of the matrix has a position (after the `in_tissue_only` filter), the files do not belong together and the loader raises `ValueError`; the message gives the barcode counts, example barcodes from the matrix and from the positions file, and, with `in_tissue_only=True`, a reminder that every off-tissue spot is dropped.
- `scalefactors_json.json` is returned in `dataset.scale_factors` (for example `tissue_hires_scalef`); the image is read with Pillow, and the loader logs a warning and continues without it when Pillow is missing.

**Errors:**

- `FileNotFoundError` when the directory, the expression matrix or the tissue-positions file is missing.
- `ValueError` when no barcode of the matrix has a tissue position (see above).

Lower-level functions: `read_tissue_positions(positions_file)` returns `TissuePosition` records (`barcode`, `in_tissue`, `array_row`, `array_col`, `pixel_row`, `pixel_col`), `filter_tissue_spots(positions, in_tissue_only=True)`, `read_spatial_image(image_path)` and `create_spatial_dataset(matrix, positions, ...)`, which builds a `SpatialDataset` from components.

## Xenium (`load_xenium`)

```python
from metainformant.spatial.io.xenium import load_xenium

dataset = load_xenium(
    'data/xenium_run/',
    load_transcripts=False,         # read transcripts.csv(.gz)
    load_boundaries=False,          # read cell_boundaries.csv(.gz) or .parquet
    min_transcript_quality=20.0,    # Phred quality filter for the transcripts
)
```

**Files read** from `path`:

```
xenium_run/
├── cell_feature_matrix.h5          (or the cell_feature_matrix/ MEX directory)
├── cells.csv.gz                    (or cells.csv)  required
├── transcripts.csv.gz              (optional; read only with load_transcripts=True)
└── cell_boundaries.csv.gz          (optional; read only with load_boundaries=True)
```

Counts are cell-by-gene. `cells.csv(.gz)` is **required**: it supplies the cell centroids through the columns `cell_id` (or `CellID`), `x_centroid` (or `x`, `X`) and `y_centroid` (or `y`, `Y`).

**Errors:**

- `FileNotFoundError` when the directory, the feature matrix or `cells.csv(.gz)` is missing.
- `ValueError` when `cells.csv(.gz)` lacks the cell id and centroid columns, or has no row for a cell of the feature matrix.

The loader never fills in `[0, 0]` for a cell without a position: placeholder coordinates would be indistinguishable from real ones, so a missing file, column or cell row is an error.

Lower-level functions: `read_cell_features(path)` (returns the CSR matrix, cell ids, gene names and gene ids), `read_transcripts(path, *, min_quality=20.0)` and `load_cell_boundaries(path)`.

## MERFISH (`load_merfish`)

```python
from metainformant.spatial.io.merfish import load_merfish

dataset = load_merfish(
    'data/merfish_run/',
    load_transcripts=False,         # read detected_transcripts.csv
)
```

**Files read** from `path`:

```
merfish_run/
├── cell_by_gene.csv                (or cell_by_gene.tsv, counts.csv)  first column = cell id, other columns = genes
├── cell_metadata.csv               (or cell_metadata.tsv)  required
└── detected_transcripts.csv        (or transcripts.csv; optional, read only with load_transcripts=True)
```

`cell_metadata.csv` (or `cell_metadata.tsv`; the `.csv` is used when both exist) is **required**. It needs a cell id column (`cell_id`, `CellID`, `cellID`, `cell` or `Cell`) and coordinate columns `center_x`/`center_y` (or `x`/`y`, `X`/`Y`, `centroid_x`/`centroid_y`, `global_x`/`global_y`); `center_z`, `volume` and `fov` are optional.

**Delimiters.** The field delimiter follows the file extension: a tab for `.tsv` (and `.tab`), a comma for anything else. This holds for the cell-by-gene matrix (`cell_by_gene.tsv` is tab-separated, `cell_by_gene.csv` and `counts.csv` comma-separated), for `parse_cell_metadata` and for `load_transcript_spots`. `load_merfish` itself only looks for `detected_transcripts.csv` or `transcripts.csv`; with `load_transcripts=True` and neither file it logs a warning and leaves `transcript_spots` as `None`, so read a tab-separated transcripts file by calling `load_transcript_spots(path)` yourself.

**Errors:**

- `FileNotFoundError` when the directory, the cell-by-gene matrix or the cell metadata file is missing.
- `ValueError` when the cell metadata has no id or coordinate columns, or has no row for a cell of the matrix.

As for Xenium, no cell is ever given a placeholder `[0, 0]` position.

Lower-level functions: `parse_cell_metadata(path)` (returns `CellMetadata` records), `load_transcript_spots(path)` and `aggregate_to_cells(transcript_spots, cell_boundaries, *, unassigned_label='')`, which counts transcripts per cell and gene.

## Other platforms

There is no platform registry or `load_spatial` function, and no Slide-seq or CosMx loader. For data from another platform, build the dataset yourself: `create_spatial_dataset(matrix, positions)` makes a `SpatialDataset` from a spots × genes matrix and `TissuePosition` records, and the analysis functions only need plain arrays (`expression` and `coordinates`).

## Tips

- Check `dataset.n_spots` after `load_visium`: it should match what you expect from the slide (`in_tissue_only=True` keeps only the spots under tissue).
- Gene identifiers are kept as stored: `gene_names` holds the symbols (second column of `features.tsv`) and `gene_ids` the identifiers (first column, usually Ensembl) of Visium and Xenium datasets.
- Write derived files under `output/`, not next to the raw data.
