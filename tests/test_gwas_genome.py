"""Tests for GWAS genome mapping module (Amel_HAv3.1)."""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.gwas.data.genome import (
    AMEL_HAV3_CHROM_SIZES,
    AMEL_HAV3_CHROMOSOMES,
    get_chromosome_order,
    get_genome_size,
    normalize_chromosome_name,
    parse_gff3_genes,
)


class TestChromosomeNormalization:
    """Tests for chromosome name normalization."""

    def test_ncbi_accession(self) -> None:
        """NCBI accessions should map to correct chromosome numbers."""
        assert normalize_chromosome_name("NC_037638.1") == 1
        assert normalize_chromosome_name("NC_037653.1") == 16

    def test_chr_prefix(self) -> None:
        """chr-prefixed names should normalize correctly."""
        assert normalize_chromosome_name("chr1") == 1
        assert normalize_chromosome_name("chr16") == 16
        assert normalize_chromosome_name("Chr5") == 5
        assert normalize_chromosome_name("CHR10") == 10

    def test_group_prefix(self) -> None:
        """Group-prefixed names should normalize correctly."""
        assert normalize_chromosome_name("Group1") == 1
        assert normalize_chromosome_name("group16") == 16

    def test_lg_prefix(self) -> None:
        """LG-prefixed names should normalize correctly."""
        assert normalize_chromosome_name("LG1") == 1
        assert normalize_chromosome_name("LG16") == 16

    def test_bare_integer(self) -> None:
        """Bare integers should normalize correctly."""
        assert normalize_chromosome_name("1") == 1
        assert normalize_chromosome_name("16") == 16

    def test_roman_numerals(self) -> None:
        """Roman numerals should normalize correctly."""
        assert normalize_chromosome_name("I") == 1
        assert normalize_chromosome_name("XVI") == 16
        assert normalize_chromosome_name("IX") == 9

    def test_out_of_range(self) -> None:
        """Out-of-range values should return None."""
        assert normalize_chromosome_name("0") is None
        assert normalize_chromosome_name("17") is None
        assert normalize_chromosome_name("100") is None

    def test_unknown_format(self) -> None:
        """Unknown formats should return None."""
        assert normalize_chromosome_name("scaffold_1") is None
        assert normalize_chromosome_name("MT") is None

    def test_all_ncbi_accessions(self) -> None:
        """All 16 NCBI accessions should be mapped."""
        assert len(AMEL_HAV3_CHROMOSOMES) == 16
        for accession, chrom_num in AMEL_HAV3_CHROMOSOMES.items():
            assert normalize_chromosome_name(accession) == chrom_num


class TestGenomeConstants:
    """Tests for genome constant data."""

    def test_chromosome_order(self) -> None:
        """Chromosome order should be 1..16."""
        order = get_chromosome_order()
        assert order == list(range(1, 17))

    def test_chromosome_sizes(self) -> None:
        """All 16 chromosomes should have positive sizes."""
        assert len(AMEL_HAV3_CHROM_SIZES) == 16
        for chrom, size in AMEL_HAV3_CHROM_SIZES.items():
            assert 1 <= chrom <= 16
            assert size > 0

    def test_genome_size(self) -> None:
        """Total genome size should be reasonable for A. mellifera (~230 Mbp)."""
        total = get_genome_size()
        assert 200_000_000 < total < 300_000_000


class TestGFF3Parsing:
    """Tests for GFF3 gene annotation parsing."""

    def test_parse_gff3(self, tmp_path: Path) -> None:
        """Parse a minimal GFF3 file."""
        gff_content = (
            "##gff-version 3\n"
            "chr1\tRefSeq\tgene\t1000\t5000\t.\t+\t.\tID=gene-LOC1;Name=LOC1;gene_biotype=protein_coding\n"
            "chr2\tRefSeq\tgene\t2000\t8000\t.\t-\t.\tID=gene-LOC2;Name=LOC2;gene_biotype=protein_coding\n"
        )
        gff_file = tmp_path / "test.gff3"
        gff_file.write_text(gff_content)

        genes = parse_gff3_genes(gff_file)
        assert len(genes) == 2
        assert genes[0]["chrom"] == 1
        assert genes[0]["start"] == 1000
        assert genes[0]["end"] == 5000
        assert genes[0]["strand"] == "+"
        assert genes[0]["gene_name"] == "LOC1"
        assert genes[1]["chrom"] == 2
        assert genes[1]["strand"] == "-"

    def test_parse_gff3_ncbi_accessions(self, tmp_path: Path) -> None:
        """Parse GFF3 with NCBI accession chromosome names."""
        gff_content = "##gff-version 3\n" "NC_037638.1\tRefSeq\tgene\t1000\t5000\t.\t+\t.\tID=gene-LOC1;Name=LOC1\n"
        gff_file = tmp_path / "ncbi.gff3"
        gff_file.write_text(gff_content)

        genes = parse_gff3_genes(gff_file)
        assert len(genes) == 1
        assert genes[0]["chrom"] == 1

    def test_parse_gff3_file_not_found(self) -> None:
        """Should raise FileNotFoundError for missing file."""
        with pytest.raises(FileNotFoundError):
            parse_gff3_genes("/nonexistent/file.gff3")

    def test_parse_gff3_empty(self, tmp_path: Path) -> None:
        """Parse empty GFF3 file."""
        gff_file = tmp_path / "empty.gff3"
        gff_file.write_text("##gff-version 3\n")

        genes = parse_gff3_genes(gff_file)
        assert len(genes) == 0
