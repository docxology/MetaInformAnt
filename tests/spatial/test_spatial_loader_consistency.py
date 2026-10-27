"""Loader consistency: Visium never returns mismatched outputs; MERFISH reads tab-separated files.

- ``load_visium`` used to return the unfiltered expression matrix next to empty coordinates and barcodes when
  no barcode had a tissue position.
- ``load_merfish`` found ``cell_metadata.tsv`` but parsed it as comma-separated, which always failed.

Real files on disk and real loaders only.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse as sp_sparse
from scipy.io import mmwrite

from metainformant.spatial.io.merfish import load_merfish, load_transcript_spots, parse_cell_metadata
from metainformant.spatial.io.visium import load_visium

POSITIONS_HEADER = "barcode,in_tissue,array_row,array_col,pixel_row,pixel_col\n"


def _write_visium_run(base: Path, position_rows: list[str]) -> Path:
    run = base / "visium_run"
    mex = run / "filtered_feature_bc_matrix"
    mex.mkdir(parents=True)
    mmwrite(mex / "matrix.mtx", sp_sparse.csr_matrix(np.array([[1, 0, 2, 0], [0, 1, 0, 3], [4, 0, 0, 1]])))
    (mex / "barcodes.tsv").write_text("BC0\nBC1\nBC2\nBC3\n")
    (mex / "features.tsv").write_text("id0\tG0\nid1\tG1\nid2\tG2\n")
    spatial_dir = run / "spatial"
    spatial_dir.mkdir()
    (spatial_dir / "tissue_positions.csv").write_text(POSITIONS_HEADER + "".join(position_rows))
    return run


class TestVisiumBarcodeMatching:
    def test_no_matching_barcode_raises_instead_of_returning_unfiltered_expression(self, tmp_path: Path) -> None:
        run = _write_visium_run(
            tmp_path,
            ["OTHER0,1,0,0,10.0,20.0\n", "OTHER1,1,0,1,30.0,40.0\n"],
        )

        with pytest.raises(ValueError, match=r"No barcodes match.*BC0.*OTHER0"):
            load_visium(run, load_image=False)

    def test_every_spot_off_tissue_raises_and_names_the_filter(self, tmp_path: Path) -> None:
        run = _write_visium_run(
            tmp_path,
            [f"BC{i},0,0,{i},{10.0 * i},{20.0 * i}\n" for i in range(4)],
        )

        with pytest.raises(ValueError, match="in_tissue_only"):
            load_visium(run, load_image=False)

    def test_the_same_files_load_when_off_tissue_spots_are_kept(self, tmp_path: Path) -> None:
        run = _write_visium_run(
            tmp_path,
            [f"BC{i},0,0,{i},{10.0 * i},{20.0 * i}\n" for i in range(4)],
        )

        dataset = load_visium(run, in_tissue_only=False, load_image=False)

        assert dataset.n_spots == 4

    def test_a_partial_match_keeps_every_output_aligned(self, tmp_path: Path) -> None:
        run = _write_visium_run(
            tmp_path,
            ["BC3,1,1,1,70.0,80.0\n", "BC1,1,0,1,30.0,40.0\n", "NOPE,1,5,5,1.0,2.0\n"],
        )

        dataset = load_visium(run, load_image=False)

        assert dataset.barcodes == ["BC1", "BC3"]  # matrix order
        assert dataset.expression.shape == (2, 3)
        assert dataset.coordinates.shape == (2, 2)
        assert [p.barcode for p in dataset.tissue_positions] == ["BC1", "BC3"]
        np.testing.assert_allclose(dataset.coordinates, [[30.0, 40.0], [70.0, 80.0]])
        assert float(dataset.expression[1, 1]) == 3.0  # BC3 carries 3 counts of G1


def _write_table(path: Path, rows: list[list[str]], delimiter: str) -> None:
    with open(path, "w", newline="") as fh:
        csv.writer(fh, delimiter=delimiter).writerows(rows)


METADATA_ROWS = [
    ["cell_id", "center_x", "center_y", "volume", "fov"],
    ["c1", "10", "20", "100.5", "3"],
    ["c2", "30", "40", "200.5", "4"],
]
MATRIX_ROWS = [["cell_id", "CD3E", "CD8A"], ["c1", "5.0", "3.0"], ["c2", "1.0", "8.0"]]
SPOT_ROWS = [["gene", "x", "y", "cell_id"], ["CD3E", "1.5", "2.5", "c1"], ["CD8A", "3.5", "4.5", "c2"]]


class TestMerfishDelimiters:
    @pytest.mark.parametrize("suffix, delimiter", [(".csv", ","), (".tsv", "\t")])
    def test_cell_metadata_is_read_with_the_delimiter_of_its_extension(
        self, tmp_path: Path, suffix: str, delimiter: str
    ) -> None:
        path = tmp_path / f"cell_metadata{suffix}"
        _write_table(path, METADATA_ROWS, delimiter)

        cells = parse_cell_metadata(path)

        assert [(c.cell_id, c.x, c.y, c.volume, c.fov) for c in cells] == [
            ("c1", 10.0, 20.0, 100.5, 3),
            ("c2", 30.0, 40.0, 200.5, 4),
        ]

    def test_csv_and_tsv_give_identical_records(self, tmp_path: Path) -> None:
        _write_table(tmp_path / "m.csv", METADATA_ROWS, ",")
        _write_table(tmp_path / "m.tsv", METADATA_ROWS, "\t")

        assert parse_cell_metadata(tmp_path / "m.csv") == parse_cell_metadata(tmp_path / "m.tsv")

    @pytest.mark.parametrize("suffix, delimiter", [(".csv", ","), (".tsv", "\t")])
    def test_load_merfish_finds_and_parses_either_metadata_file(
        self, tmp_path: Path, suffix: str, delimiter: str
    ) -> None:
        _write_table(tmp_path / "cell_by_gene.csv", MATRIX_ROWS, ",")
        _write_table(tmp_path / f"cell_metadata{suffix}", METADATA_ROWS, delimiter)

        dataset = load_merfish(tmp_path)

        np.testing.assert_allclose(dataset.coordinates, [[10.0, 20.0], [30.0, 40.0]])
        assert dataset.cell_ids == ["c1", "c2"]
        assert len(dataset.cell_metadata) == 2

    def test_a_fully_tab_separated_run_loads(self, tmp_path: Path) -> None:
        _write_table(tmp_path / "cell_by_gene.tsv", MATRIX_ROWS, "\t")
        _write_table(tmp_path / "cell_metadata.tsv", METADATA_ROWS, "\t")

        dataset = load_merfish(tmp_path)

        assert dataset.gene_names == ["CD3E", "CD8A"]
        np.testing.assert_allclose(dataset.expression, [[5.0, 3.0], [1.0, 8.0]])
        np.testing.assert_allclose(dataset.coordinates, [[10.0, 20.0], [30.0, 40.0]])

    @pytest.mark.parametrize("suffix, delimiter", [(".csv", ","), (".tsv", "\t")])
    def test_transcript_spots_follow_the_same_rule(self, tmp_path: Path, suffix: str, delimiter: str) -> None:
        path = tmp_path / f"detected_transcripts{suffix}"
        _write_table(path, SPOT_ROWS, delimiter)

        spots = load_transcript_spots(path)

        assert [(s.gene, s.x, s.y, s.cell_id) for s in spots] == [("CD3E", 1.5, 2.5, "c1"), ("CD8A", 3.5, 4.5, "c2")]

    def test_an_unknown_extension_keeps_the_comma_default(self, tmp_path: Path) -> None:
        path = tmp_path / "cell_metadata.txt"
        _write_table(path, METADATA_ROWS, ",")

        assert len(parse_cell_metadata(path)) == 2
