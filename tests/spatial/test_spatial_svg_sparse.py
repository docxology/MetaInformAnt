"""``detect_spatially_variable_genes`` must accept SciPy sparse input and give the dense result.

``load_visium`` returns the expression matrix as a CSR matrix (spots x genes) and ``spatial_domains`` accepts it,
but the SVG detector began with ``np.asarray(expression, dtype=np.float64)``, which cannot read a sparse matrix
and raised ``ValueError: setting an array element with a sequence``. A sparse matrix is now densified one gene at a
time, so every statistic comes from the same values as the dense path and is bit-for-bit identical to it.

All inputs are real matrices and the Visium case goes through the real ``load_visium`` reader; no test doubles
(real-implementation policy).
"""

from __future__ import annotations

import tracemalloc
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from scipy import sparse
from scipy.io import mmwrite

from metainformant.spatial.analysis.spatially_variable_genes import SVGResult, detect_spatially_variable_genes
from metainformant.spatial.io.visium import load_visium

GENES = ["GRADIENT", "HOTSPOT", "NOISE_A", "NOISE_B", "SPARSE_A", "SPARSE_B", "CONSTANT", "ZERO", "RARE"]
SIDE = 8


def _grid_coordinates() -> np.ndarray:
    xs, ys = np.meshgrid(np.arange(SIDE) / SIDE, np.arange(SIDE) / SIDE)
    return np.column_stack([xs.ravel(), ys.ravel()])


def _counts(coords: np.ndarray) -> np.ndarray:
    """Integer counts, spots x genes: spatial genes, noise, sparse genes, and genes the filters drop."""
    rng = np.random.default_rng(2024)
    x, y = coords[:, 0], coords[:, 1]
    n_spots = len(coords)
    counts = np.zeros((n_spots, len(GENES)), dtype=np.int64)
    counts[:, 0] = np.round(10 * (x + y)).astype(int)  # smooth gradient
    counts[:, 1] = rng.poisson(6 * np.exp(-((x - 0.5) ** 2 + (y - 0.4) ** 2) / 0.05))  # a hotspot
    counts[:, 2] = rng.poisson(4.0, n_spots)
    counts[:, 3] = rng.poisson(1.5, n_spots)
    counts[:, 4] = rng.poisson(0.3, n_spots)
    counts[:, 5] = (rng.random(n_spots) < 0.1) * rng.poisson(5.0, n_spots)
    counts[:, 6] = 7  # constant: zero variance
    counts[:, 7] = 0  # never expressed: zero variance
    counts[10, 8] = 3  # expressed at a single spot
    return counts


def _assert_identical(sparse_result: SVGResult, dense_result: SVGResult) -> None:
    assert sparse_result.genes == dense_result.genes
    assert sparse_result.n_spots == dense_result.n_spots
    for name in ("morans_i", "z_scores", "p_values"):
        np.testing.assert_array_equal(getattr(sparse_result, name), getattr(dense_result, name), err_msg=name)
    if dense_result.q_values is None:
        assert sparse_result.q_values is None
    else:
        np.testing.assert_array_equal(sparse_result.q_values, dense_result.q_values, err_msg="q_values")


@pytest.fixture(scope="module")
def coords() -> np.ndarray:
    return _grid_coordinates()


@pytest.fixture(scope="module")
def counts(coords: np.ndarray) -> np.ndarray:
    return _counts(coords)


@pytest.mark.parametrize("fmt", ["csr", "csc", "coo", "lil", "dok", "bsr"])
@pytest.mark.parametrize("orientation", ["spots_by_genes", "genes_by_spots"])
def test_sparse_input_gives_the_dense_result(
    coords: np.ndarray, counts: np.ndarray, fmt: str, orientation: str
) -> None:
    matrix = counts if orientation == "spots_by_genes" else counts.T
    dense = detect_spatially_variable_genes(matrix, coords, genes=GENES, k=4)

    result = detect_spatially_variable_genes(sparse.csr_matrix(matrix).asformat(fmt), coords, genes=GENES, k=4)

    _assert_identical(result, dense)
    # the filters dropped exactly the constant and the never-expressed gene, and the ranking is not trivial
    assert sorted(result.genes) == sorted(set(GENES) - {"CONSTANT", "ZERO"})
    assert np.all(np.diff(result.morans_i) <= 0)
    assert result.morans_i[0] > 0.5


@pytest.mark.parametrize("make", [sparse.csr_array, sparse.csc_array, sparse.coo_array])
def test_sparse_arrays_are_accepted_like_sparse_matrices(coords: np.ndarray, counts: np.ndarray, make: Any) -> None:
    dense = detect_spatially_variable_genes(counts, coords, genes=GENES, k=4)

    _assert_identical(detect_spatially_variable_genes(make(counts), coords, genes=GENES, k=4), dense)


@pytest.mark.parametrize("dtype", [np.float32, np.float64, np.int32])
def test_value_type_does_not_change_the_result(coords: np.ndarray, counts: np.ndarray, dtype: Any) -> None:
    typed = counts.astype(dtype)
    dense = detect_spatially_variable_genes(typed, coords, genes=GENES, k=4)

    _assert_identical(detect_spatially_variable_genes(sparse.csr_matrix(typed), coords, genes=GENES, k=4), dense)


def test_expression_and_variance_filters_agree_for_sparse_input(coords: np.ndarray, counts: np.ndarray) -> None:
    options = {"min_mean_expression": 1.0, "min_variance": 0.5}
    dense = detect_spatially_variable_genes(counts, coords, genes=GENES, k=4, **options)

    result = detect_spatially_variable_genes(sparse.csr_matrix(counts), coords, genes=GENES, k=4, **options)

    _assert_identical(result, dense)
    assert 0 < len(result.genes) < len(set(GENES) - {"CONSTANT", "ZERO"})  # the filters really removed some genes


def test_explicit_zeros_and_duplicate_entries_count_as_the_dense_matrix_does(
    coords: np.ndarray, counts: np.ndarray
) -> None:
    rows, cols = np.nonzero(counts)
    values = counts[rows, cols].astype(float)
    # split every value into two stored entries for the same cell, and store an explicit zero as well
    coo = sparse.coo_matrix(
        (
            np.concatenate([values - 1, np.ones_like(values), [0.0]]),
            (np.concatenate([rows, rows, [0]]), np.concatenate([cols, cols, [0]])),
        ),
        shape=counts.shape,
    )
    assert coo.nnz == 2 * len(values) + 1  # duplicates are still stored separately
    dense = detect_spatially_variable_genes(coo.toarray(), coords, genes=GENES, k=4)  # toarray() sums duplicates

    _assert_identical(detect_spatially_variable_genes(coo, coords, genes=GENES, k=4), dense)
    _assert_identical(detect_spatially_variable_genes(coo.tocsr(), coords, genes=GENES, k=4), dense)

    # a CSR matrix may also keep the duplicates inside its rows (it is then not in canonical format)
    order = np.argsort(coo.row, kind="stable")
    indptr = np.concatenate([[0], np.cumsum(np.bincount(coo.row, minlength=counts.shape[0]))])
    raw = sparse.csr_matrix((coo.data[order], coo.col[order], indptr), shape=counts.shape)
    assert not raw.has_canonical_format
    _assert_identical(detect_spatially_variable_genes(raw, coords, genes=GENES, k=4), dense)


def test_the_callers_sparse_matrix_is_left_untouched(coords: np.ndarray, counts: np.ndarray) -> None:
    matrix = sparse.csr_matrix(counts)
    before = (matrix.indptr.copy(), matrix.indices.copy(), matrix.data.copy(), matrix.dtype, matrix.shape)

    detect_spatially_variable_genes(matrix, coords, genes=GENES, k=4)

    assert matrix.format == "csr"
    np.testing.assert_array_equal(matrix.indptr, before[0])
    np.testing.assert_array_equal(matrix.indices, before[1])
    np.testing.assert_array_equal(matrix.data, before[2])
    assert (matrix.dtype, matrix.shape) == before[3:]


def test_sparse_input_is_validated_like_dense_input(coords: np.ndarray, counts: np.ndarray) -> None:
    with pytest.raises(ValueError, match="columns after orientation"):
        detect_spatially_variable_genes(sparse.csr_matrix(counts[:, :5].T[:, :7]), coords, k=4)
    with pytest.raises(ValueError, match="gene names"):
        detect_spatially_variable_genes(sparse.csr_matrix(counts), coords, genes=GENES[:3], k=4)
    with pytest.raises(ValueError, match="2-D"):
        detect_spatially_variable_genes(sparse.coo_array(np.arange(float(len(coords)))), coords, k=4)


def test_default_gene_names_and_top_genes_work_for_sparse_input(coords: np.ndarray, counts: np.ndarray) -> None:
    dense = detect_spatially_variable_genes(counts, coords, k=4)

    result = detect_spatially_variable_genes(sparse.csr_matrix(counts), coords, k=4)

    assert result.genes == dense.genes
    assert result.top_genes(3) == dense.top_genes(3)
    assert set(result.genes) <= {f"Gene{i}" for i in range(len(GENES))}


def _write_visium_run(base: Path, coords: np.ndarray, counts: np.ndarray) -> Path:
    """Write a Visium run in 10x MEX layout (genes x barcodes), as Space Ranger does."""
    run = base / "visium_run"
    mex = run / "filtered_feature_bc_matrix"
    mex.mkdir(parents=True)
    mmwrite(mex / "matrix.mtx", sparse.csr_matrix(counts.T))
    (mex / "barcodes.tsv").write_text("".join(f"BC{i}\n" for i in range(len(coords))))
    (mex / "features.tsv").write_text("".join(f"id{i}\t{name}\n" for i, name in enumerate(GENES)))
    spatial_dir = run / "spatial"
    spatial_dir.mkdir()
    rows = ["barcode,in_tissue,array_row,array_col,pixel_row,pixel_col"]
    for i, (x, y) in enumerate(coords):
        rows.append(f"BC{i},1,{i // SIDE},{i % SIDE},{x * 100 + 5.0},{y * 100 + 5.0}")
    (spatial_dir / "tissue_positions.csv").write_text("\n".join(rows) + "\n")
    (spatial_dir / "scalefactors_json.json").write_text('{"spot_diameter_fullres": 65.0}')
    return run


def test_expression_returned_by_load_visium_can_be_ranked(
    tmp_path: Path, coords: np.ndarray, counts: np.ndarray
) -> None:
    dataset = load_visium(_write_visium_run(tmp_path, coords, counts), load_image=False)
    assert sparse.issparse(dataset.expression)  # this is what the detector used to reject
    assert dataset.expression.shape == (len(coords), len(GENES))

    result = detect_spatially_variable_genes(dataset.expression, dataset.coordinates, genes=dataset.gene_names, k=4)
    dense = detect_spatially_variable_genes(
        dataset.expression.toarray(), dataset.coordinates, genes=dataset.gene_names, k=4
    )

    _assert_identical(result, dense)
    assert "GRADIENT" in result.top_genes(2)
    assert result.morans_i[result.genes.index("GRADIENT")] > 0.7
    assert "CONSTANT" not in result.genes and "ZERO" not in result.genes


def test_a_whole_slide_matrix_is_not_densified() -> None:
    """Rows are densified one at a time: peak memory stays far below the size of the dense matrix."""
    n_side, n_genes = 32, 12000
    coords = np.column_stack([a.ravel() for a in np.meshgrid(np.arange(n_side), np.arange(n_side))]).astype(float)
    n_spots = len(coords)
    rng = np.random.default_rng(5)
    rows = rng.integers(0, n_spots, 4000)
    cols = rng.integers(0, 40, 4000)  # only 40 genes carry counts; all the others are zero
    matrix = sparse.csr_matrix((rng.poisson(3.0, 4000) + 1.0, (rows, cols)), shape=(n_spots, n_genes))
    dense_bytes = n_spots * n_genes * 8
    assert dense_bytes > 60e6

    tracemalloc.start()
    try:
        result = detect_spatially_variable_genes(matrix, coords, k=4)
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()

    assert 0 < len(result.genes) <= 40
    assert peak < dense_bytes / 4, f"peak {peak / 1e6:.1f} MB for a {dense_bytes / 1e6:.0f} MB dense matrix"
