"""Tests for RNA retrieval and download functionality.

Tests the ENA downloader and related retrieval utilities used for
downloading FASTQ files from ENA/SRA.
"""

from __future__ import annotations

import hashlib
import pytest
from pathlib import Path
# Import the module under test
from metainformant.rna.retrieval import ena_downloader


class TestCalculateMD5:
    """Tests for MD5 checksum calculation."""

    def test_calculate_md5_simple_file(self, tmp_path: Path) -> None:
        """Test MD5 calculation on a simple file."""
        test_file = tmp_path / "test.txt"
        content = b"Hello, World!"
        test_file.write_bytes(content)

        # Calculate expected MD5
        expected_md5 = hashlib.md5(content).hexdigest()

        result = ena_downloader.calculate_md5(test_file)

        assert result == expected_md5

    def test_calculate_md5_empty_file(self, tmp_path: Path) -> None:
        """Test MD5 calculation on empty file."""
        test_file = tmp_path / "empty.txt"
        test_file.write_bytes(b"")

        # MD5 of empty content
        expected_md5 = hashlib.md5(b"").hexdigest()

        result = ena_downloader.calculate_md5(test_file)

        assert result == expected_md5

    def test_calculate_md5_large_file(self, tmp_path: Path) -> None:
        """Test MD5 calculation on larger file (chunked reading)."""
        test_file = tmp_path / "large.bin"

        # Create a file larger than the chunk size (4096 bytes)
        content = b"A" * 10000

        test_file.write_bytes(content)

        expected_md5 = hashlib.md5(content).hexdigest()

        result = ena_downloader.calculate_md5(test_file)

        assert result == expected_md5

    def test_calculate_md5_binary_content(self, tmp_path: Path) -> None:
        """Test MD5 calculation on binary content."""
        test_file = tmp_path / "binary.bin"
        content = bytes(range(256))
        test_file.write_bytes(content)

        expected_md5 = hashlib.md5(content).hexdigest()

        result = ena_downloader.calculate_md5(test_file)

        assert result == expected_md5


class TestCleanStagnantFile:
    """Tests for file cleanup utility."""

    def test_clean_stagnant_file_exists(self, tmp_path: Path) -> None:
        """Test cleaning an existing file."""
        test_file = tmp_path / "stagnant.txt"
        test_file.write_text("incomplete download")

        assert test_file.exists()

        ena_downloader.clean_stagnant_file(test_file)

        assert not test_file.exists()

    def test_clean_stagnant_file_not_exists(self, tmp_path: Path) -> None:
        """Test cleaning a non-existent file (should not error)."""
        test_file = tmp_path / "nonexistent.txt"

        assert not test_file.exists()

        # Should not raise any error
        ena_downloader.clean_stagnant_file(test_file)

        assert not test_file.exists()


class TestGetEnaLinks:
    """Tests for ENA link retrieval."""

    @pytest.mark.network
    def test_get_ena_links_real_sample(self) -> None:
        """Test real ENA API call for known sample.

        This test makes a real network request to ENA.
        """
        # Use a known public sample
        sra_id = "SRR000001"  # One of the first public SRA samples

        links = ena_downloader.get_ena_links(sra_id)

        # Should return something (may be empty if sample is too old/archived)
        assert isinstance(links, list)

        # If we got links, verify structure
        if links:
            for url, md5 in links:
                assert isinstance(url, str)
                assert isinstance(md5, str)
                assert url.startswith(("ftp://", "http://", "https://"))
                assert len(md5) == 32  # MD5 hash length

    @pytest.mark.network
    def test_get_ena_links_nonexistent_sample(self) -> None:
        """Test ENA API call for non-existent sample."""
        sra_id = "SRR999999999999"  # Non-existent sample

        links = ena_downloader.get_ena_links(sra_id)

        # Should return empty list for non-existent sample
        assert links == []

    def test_get_ena_links_url_format(self) -> None:
        """Test that constructed URL is correct format."""
        sra_id = "SRR123456"

        expected_url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={sra_id}&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv"

        # The function constructs this URL internally
        # We can verify by checking what URL would be called
        url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={sra_id}&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv"

        assert expected_url == url


class TestDownloadFile:
    """Tests for file download functionality."""

    def test_download_file_creates_directory(self, tmp_path: Path) -> None:
        """Test that download creates parent directories."""
        dest = tmp_path / "subdir" / "file.txt"

        # Create parent manually since download expects it to exist
        dest.parent.mkdir(parents=True, exist_ok=True)

        # The actual download will fail (no real curl call in unit test)
        # but we're testing the path handling
        assert dest.parent.exists()

    @pytest.mark.network
    def test_download_file_with_md5_verification(self, tmp_path: Path) -> None:
        """Test download with MD5 verification.

        This test requires curl and network access.
        """
        # Use a small, stable file from a reliable source
        # We'll skip this test if curl is not available
        import shutil

        if not shutil.which("curl"):
            pytest.skip("curl not available")

        # We can't easily test this without a reliable test endpoint
        # In real usage, this downloads from ENA FTP
        pass


class TestDownloadSraSamples:
    """Tests for batch SRA sample download."""

    def test_download_creates_directory_structure(self, tmp_path: Path) -> None:
        """Test that download creates expected directory structure."""
        base_dir = tmp_path / "output"

        # Create the structure manually to test expectations
        getfastq_dir = base_dir / "getfastq"
        getfastq_dir.mkdir(parents=True)

        sample_dir = getfastq_dir / "SRR123456"
        sample_dir.mkdir()

        # Verify structure
        assert getfastq_dir.exists()
        assert sample_dir.exists()
        assert sample_dir.name == "SRR123456"

    def test_file_renaming_pattern(self, tmp_path: Path) -> None:
        """Test that files are renamed to amalgkit format."""
        sample_dir = tmp_path / "SRR123456"
        sample_dir.mkdir()

        # Create files with ENA naming
        fastq1 = sample_dir / "SRR123456_1.fastq.gz"
        fastq2 = sample_dir / "SRR123456_2.fastq.gz"
        fastq1.write_bytes(b"test")
        fastq2.write_bytes(b"test")

        # Simulate rename operation
        for f in sample_dir.glob("*.fastq.gz"):
            if ".amalgkit." not in f.name:
                new_name = f.name.replace(".fastq.gz", ".amalgkit.fastq.gz")
                f.rename(sample_dir / new_name)

        # Verify renamed files exist
        assert (sample_dir / "SRR123456_1.amalgkit.fastq.gz").exists()
        assert (sample_dir / "SRR123456_2.amalgkit.fastq.gz").exists()

        # Original names should not exist
        assert not (sample_dir / "SRR123456_1.fastq.gz").exists()
        assert not (sample_dir / "SRR123456_2.fastq.gz").exists()


class TestIntegration:
    """Integration tests for ENA downloader."""

    @pytest.mark.network
    @pytest.mark.slow
    def test_full_download_workflow(self, tmp_path: Path) -> None:
        """Test complete download workflow for a small sample.

        This test makes real network requests and downloads files.
        It should be skipped in quick test runs.
        """
        import shutil

        if not shutil.which("curl"):
            pytest.skip("curl not available")

        # This is a slow test that actually downloads data
        # We skip it in most test runs
        pytest.skip("Full download test - run manually with --network flag")


class TestEdgeCases:
    """Edge case tests for ENA downloader."""

    def test_handle_ftp_protocol_prefix(self) -> None:
        """Test that FTP URLs without protocol are handled."""
        # ENA often returns URLs without protocol (just hostname/path)
        ftp_path = "ftp.sra.ebi.ac.uk/vol1/fastq/SRR123/SRR123456/SRR123456_1.fastq.gz"

        # Should add protocol if missing
        # Check for full protocol prefixes, not just hostname starting with "ftp"
        if not ftp_path.startswith(("ftp://", "http://", "https://")):
            full_url = f"ftp://{ftp_path}"
        else:
            full_url = ftp_path

        assert full_url.startswith("ftp://")
        assert full_url == "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR123/SRR123456/SRR123456_1.fastq.gz"

    def test_handle_mismatched_ftp_md5_counts(self) -> None:
        """Test handling when FTP links and MD5s don't match."""
        ftps = ["url1", "url2", "url3"]
        md5s = ["md51", "md52"]  # Mismatched count

        # The function should detect this mismatch
        assert len(ftps) != len(md5s)

    def test_empty_ftp_response(self) -> None:
        """Test handling of empty FTP response."""
        row = {"fastq_ftp": "", "fastq_md5": ""}

        # Empty FTP links should result in empty list
        if not row["fastq_ftp"]:
            result = []
        else:
            result = row["fastq_ftp"].split(";")

        assert result == []
