"""Xenium and MERFISH loaders refuse to invent [0, 0] coordinates.

A cell without a position used to be placed at the origin (silently when cells.csv / cell_metadata.csv was
missing, and per cell when an id had no row), which downstream spatial statistics read as a real
location. Real files on disk and real loaders only.
"""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse as sp_sparse
from scipy.io import mmwrite

from metainformant.spatial.io.merfish import load_merfish
from metainformant.spatial.io.xenium import load_xenium

CELLS = ["c1", "c2", "c3"]


def _write_xenium_matrix(base: Path) -> None:
    mex = base / "cell_feature_matrix"
    mex.mkdir(parents=True)
    mmwrite(mex / "matrix.mtx", sp_sparse.csr_matrix(np.array([[5, 0, 1], [0, 7, 2]])))  # genes x cells
    (mex / "barcodes.tsv").write_text("\n".join(CELLS) + "\n")
    (mex / "features.tsv").write_text("gid0\tCD3E\ngid1\tCD8A\n")


def _write_merfish_matrix(base: Path) -> None:
    base.mkdir(parents=True, exist_ok=True)
    with open(base / "cell_by_gene.csv", "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["cell_id", "CD3E", "CD8A"])
        for i, cell in enumerate(CELLS):
            writer.writerow([cell, 5.0 + i, 3.0 - i])


class TestXenium:
    def test_missing_cells_file_raises(self, tmp_path: Path) -> None:
        _write_xenium_matrix(tmp_path)

        with pytest.raises(FileNotFoundError, match=r"cells\.csv"):
            load_xenium(tmp_path)

    def test_a_cell_without_a_row_raises_and_names_it(self, tmp_path: Path) -> None:
        _write_xenium_matrix(tmp_path)
        (tmp_path / "cells.csv").write_text("cell_id,x_centroid,y_centroid\nc1,1.5,2.5\nc3,5.5,6.5\n")

        with pytest.raises(ValueError, match=r"1 of 3 cells.*c2.*\[0, 0\]"):
            load_xenium(tmp_path)

    def test_a_cells_file_without_the_position_columns_raises(self, tmp_path: Path) -> None:
        _write_xenium_matrix(tmp_path)
        (tmp_path / "cells.csv").write_text("cell_id,area\nc1,10\nc2,11\nc3,12\n")

        with pytest.raises(ValueError, match="x_centroid"):
            load_xenium(tmp_path)

    def test_a_gzipped_cells_file_loads_real_positions(self, tmp_path: Path) -> None:
        _write_xenium_matrix(tmp_path)
        with gzip.open(tmp_path / "cells.csv.gz", "wt") as fh:
            fh.write("cell_id,x_centroid,y_centroid\nc3,5.5,6.5\nc1,1.5,2.5\nc2,3.5,4.5\n")

        dataset = load_xenium(tmp_path)

        np.testing.assert_allclose(dataset.coordinates, [[1.5, 2.5], [3.5, 4.5], [5.5, 6.5]])  # in matrix order

    def test_a_cell_exactly_at_the_origin_is_still_accepted(self, tmp_path: Path) -> None:
        _write_xenium_matrix(tmp_path)
        (tmp_path / "cells.csv").write_text("cell_id,x_centroid,y_centroid\nc1,0.0,0.0\nc2,3.5,4.5\nc3,5.5,6.5\n")

        dataset = load_xenium(tmp_path)

        np.testing.assert_allclose(dataset.coordinates[0], [0.0, 0.0])


class TestMerfish:
    def test_missing_metadata_file_raises(self, tmp_path: Path) -> None:
        _write_merfish_matrix(tmp_path)

        with pytest.raises(FileNotFoundError, match="cell_metadata"):
            load_merfish(tmp_path)

    def test_a_cell_without_a_row_raises_and_names_it(self, tmp_path: Path) -> None:
        _write_merfish_matrix(tmp_path)
        (tmp_path / "cell_metadata.csv").write_text("cell_id,center_x,center_y\nc1,1,2\nc2,3,4\n")

        with pytest.raises(ValueError, match=r"1 of 3 cells.*c3.*\[0, 0\]"):
            load_merfish(tmp_path)

    def test_coordinates_follow_the_matrix_order_not_the_metadata_order(self, tmp_path: Path) -> None:
        _write_merfish_matrix(tmp_path)
        (tmp_path / "cell_metadata.csv").write_text("cell_id,center_x,center_y\nc3,50,60\nc1,10,20\nc2,30,40\n")

        dataset = load_merfish(tmp_path)

        np.testing.assert_allclose(dataset.coordinates, [[10, 20], [30, 40], [50, 60]])
        assert dataset.cell_ids == CELLS
        assert len(dataset.cell_metadata) == 3
