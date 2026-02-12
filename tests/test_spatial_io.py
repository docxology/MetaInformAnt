"""Tests for metainformant.spatial.io -- data loading for Visium, MERFISH, Xenium.

All tests use real implementations (NO mocking policy).
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

from metainformant.spatial.io.merfish import (
    CellMetadata,
    MERFISHDataset,
    TranscriptSpot,
    aggregate_to_cells,
    load_merfish,
    load_transcript_spots,
    parse_cell_metadata,
)
from metainformant.spatial.io.visium import (
    SpatialDataset,
    TissuePosition,
    create_spatial_dataset,
    filter_tissue_spots,
    read_tissue_positions,
)
from metainformant.spatial.io.xenium import (
    CellBoundary,
    XeniumDataset,
    XeniumTranscript,
    load_cell_boundaries,
    load_xenium,
    read_transcripts,
)

# ---------------------------------------------------------------------------
# Visium dataclasses and utilities
# ---------------------------------------------------------------------------


class TestSpatialDataset:
    def test_creation_with_numpy(self) -> None:
        expr = np.random.RandomState(42).standard_normal((10, 5))
        coords = np.random.RandomState(42).uniform(0, 100, (10, 2))
        ds = SpatialDataset(
            expression=expr,
            coordinates=coords,
            barcodes=[f"BC{i}" for i in range(10)],
            gene_names=["G0", "G1", "G2", "G3", "G4"],
        )
        assert ds.n_spots == 10
        assert ds.n_genes == 5
        assert ds.platform == "visium"

    def test_creation_with_sparse(self) -> None:
        from scipy import sparse as sp_sparse

        expr = sp_sparse.random(10, 5, density=0.5, random_state=42)
        coords = np.zeros((10, 2))
        ds = SpatialDataset(
            expression=expr,
            coordinates=coords,
            barcodes=[f"BC{i}" for i in range(10)],
            gene_names=["G0", "G1", "G2", "G3", "G4"],
        )
        assert ds.n_spots == 10
        assert ds.n_genes == 5


class TestTissuePosition:
    def test_dataclass_attributes(self) -> None:
        tp = TissuePosition(
            barcode="AACG-1",
            in_tissue=True,
            array_row=5,
            array_col=10,
            pixel_row=1024.5,
            pixel_col=2048.3,
        )
        assert tp.barcode == "AACG-1"
        assert tp.in_tissue is True
        assert tp.array_row == 5
        assert tp.array_col == 10
        assert tp.pixel_row == 1024.5
        assert tp.pixel_col == 2048.3


class TestFilterTissueSpots:
    def test_filters_non_tissue(self) -> None:
        positions = [
            TissuePosition("A", True, 0, 0, 0.0, 0.0),
            TissuePosition("B", False, 0, 1, 1.0, 1.0),
            TissuePosition("C", True, 1, 0, 2.0, 2.0),
        ]
        filtered = filter_tissue_spots(positions, in_tissue_only=True)
        assert len(filtered) == 2
        assert all(p.in_tissue for p in filtered)

    def test_no_filter(self) -> None:
        positions = [
            TissuePosition("A", True, 0, 0, 0.0, 0.0),
            TissuePosition("B", False, 0, 1, 1.0, 1.0),
        ]
        filtered = filter_tissue_spots(positions, in_tissue_only=False)
        assert len(filtered) == 2


class TestCreateSpatialDataset:
    def test_from_components(self) -> None:
        expr = np.random.RandomState(42).standard_normal((3, 4))
        positions = [TissuePosition(f"BC{i}", True, i, i, float(i * 10), float(i * 20)) for i in range(3)]
        ds = create_spatial_dataset(
            expr,
            positions,
            gene_names=["G0", "G1", "G2", "G3"],
        )
        assert isinstance(ds, SpatialDataset)
        assert ds.n_spots == 3
        assert ds.n_genes == 4
        assert ds.barcodes == ["BC0", "BC1", "BC2"]

    def test_auto_gene_names(self) -> None:
        expr = np.zeros((2, 3))
        positions = [
            TissuePosition("A", True, 0, 0, 0.0, 0.0),
            TissuePosition("B", True, 0, 1, 1.0, 1.0),
        ]
        ds = create_spatial_dataset(expr, positions)
        assert ds.gene_names == ["G0", "G1", "G2"]


class TestReadTissuePositions:
    def test_read_csv_v1_format(self, tmp_path: Path) -> None:
        csv_file = tmp_path / "tissue_positions_list.csv"
        with open(csv_file, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["AACG-1", "1", "0", "0", "100.5", "200.3"])
            writer.writerow(["TTCG-1", "0", "1", "1", "300.0", "400.0"])
        positions = read_tissue_positions(csv_file)
        assert len(positions) == 2
        assert positions[0].barcode == "AACG-1"
        assert positions[0].in_tissue is True
        assert positions[1].in_tissue is False

    def test_read_csv_v2_format_with_header(self, tmp_path: Path) -> None:
        csv_file = tmp_path / "tissue_positions.csv"
        with open(csv_file, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"])
            writer.writerow(["AACG-1", "1", "5", "10", "512.0", "768.0"])
        positions = read_tissue_positions(csv_file)
        assert len(positions) == 1
        assert positions[0].array_row == 5

    def test_missing_file_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            read_tissue_positions(tmp_path / "nonexistent.csv")


# ---------------------------------------------------------------------------
# MERFISH dataclasses and loaders
# ---------------------------------------------------------------------------


class TestMERFISHDataset:
    def test_creation(self) -> None:
        expr = np.zeros((5, 3))
        coords = np.zeros((5, 2))
        ds = MERFISHDataset(
            expression=expr,
            coordinates=coords,
            cell_ids=["c0", "c1", "c2", "c3", "c4"],
            gene_names=["g0", "g1", "g2"],
        )
        assert ds.n_cells == 5
        assert ds.n_genes == 3


class TestTranscriptSpot:
    def test_dataclass(self) -> None:
        spot = TranscriptSpot(gene="CD3E", x=10.5, y=20.3, z=1.0, cell_id="c1", quality=0.95)
        assert spot.gene == "CD3E"
        assert spot.x == 10.5
        assert spot.quality == 0.95


class TestCellMetadata:
    def test_dataclass(self) -> None:
        meta = CellMetadata(cell_id="c1", x=100.0, y=200.0, z=0.0, volume=150.0, fov=3)
        assert meta.cell_id == "c1"
        assert meta.volume == 150.0
        assert meta.fov == 3


class TestParseCellMetadata:
    def test_parse_csv(self, tmp_path: Path) -> None:
        csv_file = tmp_path / "cell_metadata.csv"
        with open(csv_file, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=["cell_id", "center_x", "center_y", "volume"])
            writer.writeheader()
            writer.writerow({"cell_id": "c1", "center_x": "10.0", "center_y": "20.0", "volume": "100"})
            writer.writerow({"cell_id": "c2", "center_x": "30.0", "center_y": "40.0", "volume": "200"})
        cells = parse_cell_metadata(csv_file)
        assert len(cells) == 2
        assert cells[0].cell_id == "c1"
        assert cells[0].x == 10.0
        assert cells[1].volume == 200.0

    def test_missing_file_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            parse_cell_metadata(tmp_path / "nonexistent.csv")


class TestLoadTranscriptSpots:
    def test_load_csv(self, tmp_path: Path) -> None:
        csv_file = tmp_path / "detected_transcripts.csv"
        with open(csv_file, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=["gene", "x", "y", "cell_id"])
            writer.writeheader()
            writer.writerow({"gene": "CD3E", "x": "10.0", "y": "20.0", "cell_id": "c1"})
            writer.writerow({"gene": "CD8A", "x": "30.0", "y": "40.0", "cell_id": "c2"})
        spots = load_transcript_spots(csv_file)
        assert len(spots) == 2
        assert spots[0].gene == "CD3E"

    def test_missing_file_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            load_transcript_spots(tmp_path / "nonexistent.csv")


class TestAggregateToCell:
    def test_aggregation_with_preassigned_cells(self) -> None:
        spots = [
            TranscriptSpot(gene="CD3E", x=1.0, y=1.0, cell_id="c1"),
            TranscriptSpot(gene="CD3E", x=1.0, y=1.0, cell_id="c1"),
            TranscriptSpot(gene="CD8A", x=2.0, y=2.0, cell_id="c2"),
        ]
        cells = [
            CellMetadata(cell_id="c1", x=1.0, y=1.0),
            CellMetadata(cell_id="c2", x=2.0, y=2.0),
        ]
        counts, cell_ids, gene_names = aggregate_to_cells(spots, cells)
        assert counts.shape == (2, 2)  # 2 cells, 2 genes
        assert cell_ids == ["c1", "c2"]
        # c1 should have 2 CD3E transcripts
        cd3e_idx = gene_names.index("CD3E")
        assert counts[0, cd3e_idx] == 2


class TestLoadMerfish:
    def test_missing_directory_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            load_merfish(tmp_path / "nonexistent_dir")

    def test_load_from_directory(self, tmp_path: Path) -> None:
        merfish_dir = tmp_path / "merfish"
        merfish_dir.mkdir()

        # Create cell_by_gene.csv
        with open(merfish_dir / "cell_by_gene.csv", "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["cell_id", "CD3E", "CD8A"])
            writer.writerow(["c1", "5.0", "3.0"])
            writer.writerow(["c2", "1.0", "8.0"])

        # Create cell_metadata.csv
        with open(merfish_dir / "cell_metadata.csv", "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=["cell_id", "center_x", "center_y"])
            writer.writeheader()
            writer.writerow({"cell_id": "c1", "center_x": "10", "center_y": "20"})
            writer.writerow({"cell_id": "c2", "center_x": "30", "center_y": "40"})

        ds = load_merfish(merfish_dir)
        assert isinstance(ds, MERFISHDataset)
        assert ds.n_cells == 2
        assert ds.n_genes == 2


# ---------------------------------------------------------------------------
# Xenium dataclasses and loaders
# ---------------------------------------------------------------------------


class TestXeniumDataset:
    def test_creation(self) -> None:
        ds = XeniumDataset(
            expression=np.zeros((3, 2)),
            coordinates=np.zeros((3, 2)),
            cell_ids=["c1", "c2", "c3"],
            gene_names=["g1", "g2"],
        )
        assert ds.n_cells == 3
        assert ds.n_genes == 2


class TestXeniumTranscript:
    def test_dataclass(self) -> None:
        tx = XeniumTranscript(
            transcript_id="t1",
            gene="CD3E",
            x=10.0,
            y=20.0,
            z=0.5,
            cell_id="c1",
            quality_value=30.0,
            overlaps_nucleus=True,
        )
        assert tx.transcript_id == "t1"
        assert tx.gene == "CD3E"
        assert tx.overlaps_nucleus is True


class TestCellBoundary:
    def test_dataclass(self) -> None:
        cb = CellBoundary(
            cell_id="c1",
            vertices=[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        )
        assert cb.cell_id == "c1"
        assert len(cb.vertices) == 4


class TestReadTranscripts:
    def test_read_csv(self, tmp_path: Path) -> None:
        csv_file = tmp_path / "transcripts.csv"
        with open(csv_file, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=["transcript_id", "feature_name", "x_location", "y_location", "qv"])
            writer.writeheader()
            writer.writerow(
                {"transcript_id": "t1", "feature_name": "CD3E", "x_location": "10", "y_location": "20", "qv": "30"}
            )
            writer.writerow(
                {"transcript_id": "t2", "feature_name": "CD8A", "x_location": "30", "y_location": "40", "qv": "15"}
            )
        # Default min_quality=20, so t2 should be filtered out
        transcripts = read_transcripts(csv_file, min_quality=20.0)
        assert len(transcripts) == 1
        assert transcripts[0].gene == "CD3E"

    def test_missing_file_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            read_transcripts(tmp_path / "nonexistent.csv")


class TestLoadCellBoundaries:
    def test_load_csv(self, tmp_path: Path) -> None:
        csv_file = tmp_path / "cell_boundaries.csv"
        with open(csv_file, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=["cell_id", "vertex_x", "vertex_y"])
            writer.writeheader()
            for v in [(0, 0), (1, 0), (1, 1), (0, 1)]:
                writer.writerow({"cell_id": "c1", "vertex_x": str(v[0]), "vertex_y": str(v[1])})
        boundaries = load_cell_boundaries(csv_file)
        assert len(boundaries) == 1
        assert boundaries[0].cell_id == "c1"
        assert len(boundaries[0].vertices) == 4

    def test_missing_file_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            load_cell_boundaries(tmp_path / "nonexistent.csv")


class TestLoadXenium:
    def test_missing_directory_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            load_xenium(tmp_path / "nonexistent_dir")
