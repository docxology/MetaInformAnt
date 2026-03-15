"""Tests for GWAS SNP-to-gene annotation module."""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.gwas.analysis.annotation import (
    annotate_variants_with_genes,
    classify_variant_location,
)


@pytest.fixture
def sample_genes() -> list:
    """Generate sample gene data."""
    return [
        {"chrom": 1, "start": 1000, "end": 3000, "strand": "+", "gene_id": "gene-LOC1", "gene_name": "LOC1"},
        {"chrom": 1, "start": 5000, "end": 7000, "strand": "-", "gene_id": "gene-LOC2", "gene_name": "LOC2"},
        {"chrom": 1, "start": 10000, "end": 12000, "strand": "+", "gene_id": "gene-LOC3", "gene_name": "LOC3"},
    ]


class TestClassifyVariantLocation:
    """Tests for variant location classification."""

    def test_intragenic(self, sample_genes: list) -> None:
        """Variant inside a gene should be classified as intragenic."""
        result = classify_variant_location(1, 2000, sample_genes)
        assert result["location"] == "intragenic"
        assert result["nearest_gene"] == "gene-LOC1"
        assert result["distance"] == 0

    def test_upstream_of_plus_strand(self, sample_genes: list) -> None:
        """Variant upstream (5') of a + strand gene."""
        result = classify_variant_location(1, 800, sample_genes)
        assert result["location"] == "upstream"
        assert result["nearest_gene"] == "gene-LOC1"
        assert result["distance"] == 200

    def test_downstream_of_plus_strand(self, sample_genes: list) -> None:
        """Variant downstream (3') of a + strand gene."""
        result = classify_variant_location(1, 3500, sample_genes)
        assert result["location"] == "downstream"
        assert result["nearest_gene"] == "gene-LOC1"
        assert result["distance"] == 500

    def test_upstream_of_minus_strand(self, sample_genes: list) -> None:
        """Variant upstream (3' in genomic coords) of a - strand gene."""
        # For - strand gene at 5000-7000, variant at 7500 is upstream
        result = classify_variant_location(1, 7500, sample_genes)
        assert result["location"] == "upstream"
        assert result["nearest_gene"] == "gene-LOC2"

    def test_between_genes(self, sample_genes: list) -> None:
        """Variant between two genes should find nearest."""
        result = classify_variant_location(1, 4000, sample_genes)
        assert result["nearest_gene"] is not None
        assert result["distance"] > 0

    def test_no_genes(self) -> None:
        """No genes should return intergenic with None."""
        result = classify_variant_location(1, 1000, [])
        assert result["location"] == "intergenic"
        assert result["nearest_gene"] is None

    def test_beyond_window(self, sample_genes: list) -> None:
        """Variant far from any gene should be intergenic."""
        result = classify_variant_location(1, 1000000, sample_genes, window_bp=1000)
        assert result["location"] == "intergenic"
        # Still reports nearest gene but distance > window
        assert result["distance"] is not None


class TestAnnotateVariantsWithGenes:
    """Tests for full variant annotation pipeline."""

    def test_annotate_basic(self, tmp_path: Path) -> None:
        """Annotate variants using a GFF3 file."""
        # Create GFF3
        gff_content = (
            "##gff-version 3\n"
            "chr1\tRefSeq\tgene\t1000\t3000\t.\t+\t.\tID=gene-LOC1;Name=LOC1\n"
            "chr1\tRefSeq\tgene\t5000\t7000\t.\t-\t.\tID=gene-LOC2;Name=LOC2\n"
        )
        gff_file = tmp_path / "genes.gff3"
        gff_file.write_text(gff_content)

        results = [
            {"chrom": "chr1", "pos": 2000, "beta": 0.5},
            {"chrom": "chr1", "pos": 4000, "beta": 0.1},
            {"chrom": "chr1", "pos": 6000, "beta": -0.3},
        ]

        annotated = annotate_variants_with_genes(results, gff_file, window_kb=50)

        assert len(annotated) == 3
        assert annotated[0]["annotation"]["location"] == "intragenic"
        assert annotated[0]["annotation"]["nearest_gene"] == "gene-LOC1"
        assert annotated[2]["annotation"]["location"] == "intragenic"
        assert annotated[2]["annotation"]["nearest_gene"] == "gene-LOC2"

    def test_annotate_unknown_chrom(self, tmp_path: Path) -> None:
        """Variants on unknown chromosomes should be annotated as unknown."""
        gff_content = "##gff-version 3\nchr1\tRefSeq\tgene\t1000\t3000\t.\t+\t.\tID=gene-LOC1;Name=LOC1\n"
        gff_file = tmp_path / "genes.gff3"
        gff_file.write_text(gff_content)

        results = [{"chrom": "scaffold_unknown", "pos": 1000}]
        annotated = annotate_variants_with_genes(results, gff_file)

        assert annotated[0]["annotation"]["location"] == "unknown"

    def test_annotate_ncbi_accessions(self, tmp_path: Path) -> None:
        """Variants with NCBI accession chromosomes should be annotated."""
        gff_content = "##gff-version 3\n" "NC_037638.1\tRefSeq\tgene\t1000\t3000\t.\t+\t.\tID=gene-LOC1;Name=LOC1\n"
        gff_file = tmp_path / "ncbi_genes.gff3"
        gff_file.write_text(gff_content)

        results = [{"chrom": "NC_037638.1", "pos": 2000}]
        annotated = annotate_variants_with_genes(results, gff_file)

        assert annotated[0]["annotation"]["location"] == "intragenic"
        assert annotated[0]["annotation"]["nearest_gene"] == "gene-LOC1"

    def test_annotate_gff_not_found(self) -> None:
        """Missing GFF3 should raise FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            annotate_variants_with_genes([{"chrom": "chr1", "pos": 100}], "/nonexistent.gff3")
