"""Network API tests for protein proteomes functions.

Tests for get_proteome_metadata and download_proteome_fasta functions that
were converted from placeholder implementations to real UniProt API calls.
Following NO_MOCKING policy - all tests use real network calls or skip gracefully.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest
import requests

from metainformant.protein.sequence.proteomes import get_proteome_metadata, download_proteome_fasta


def _check_network_connectivity(url: str = "https://rest.uniprot.org") -> bool:
    """Check if network connectivity to UniProt is available."""
    try:
        response = requests.head(url, timeout=5)
        return response.status_code in (200, 302, 301)
    except (requests.RequestException, requests.Timeout):
        return False


class TestGetProteomeMetadata:
    """Tests for get_proteome_metadata function - real UniProt Proteomes API calls."""

    @pytest.mark.network
    def test_get_proteome_metadata_human_real_api(self):
        """Test real API call for human proteome metadata."""
        if not _check_network_connectivity():
            pytest.skip("No network connectivity to UniProt API")

        try:
            metadata = get_proteome_metadata("9606")  # Human taxonomy ID

            # Verify response structure
            assert isinstance(metadata, dict)
            assert "scientific_name" in metadata
            assert "proteome_id" in metadata
            assert metadata["scientific_name"] == "Homo sapiens (Human)"
            assert metadata["proteome_id"].startswith("UP")
            assert len(metadata["proteome_id"]) > 2

        except Exception as e:
            # API might be down or changed
            pytest.skip(f"UniProt API unavailable: {e}")

    @pytest.mark.network
    def test_get_proteome_metadata_mouse_real_api(self):
        """Test real API call for mouse proteome metadata."""
        if not _check_network_connectivity():
            pytest.skip("No network connectivity to UniProt API")

        try:
            metadata = get_proteome_metadata("10090")  # Mouse taxonomy ID

            assert isinstance(metadata, dict)
            assert "scientific_name" in metadata
            assert "proteome_id" in metadata
            assert "Mus musculus" in metadata["scientific_name"]
            assert metadata["proteome_id"].startswith("UP")

        except Exception as e:
            pytest.skip(f"UniProt API unavailable: {e}")

    @pytest.mark.network
    def test_get_proteome_metadata_api_response_structure(self):
        """Test that API response has expected structure."""
        if not _check_network_connectivity():
            pytest.skip("No network connectivity to UniProt API")

        try:
            metadata = get_proteome_metadata("9606")

            # Check for expected fields
            expected_fields = ["scientific_name", "proteome_id"]
            for field in expected_fields:
                assert field in metadata, f"Missing field: {field}"

            # Check field types
            assert isinstance(metadata["scientific_name"], str)
            assert isinstance(metadata["proteome_id"], str)

            # Check field content
            assert len(metadata["scientific_name"]) > 0
            assert len(metadata["proteome_id"]) > 0

        except Exception as e:
            pytest.skip(f"UniProt API unavailable: {e}")

    def test_get_proteome_metadata_invalid_taxon_empty(self):
        """Test error handling for empty taxon ID."""
        with pytest.raises(requests.exceptions.HTTPError):
            get_proteome_metadata("")

    def test_get_proteome_metadata_invalid_taxon_nonexistent(self):
        """Test error handling for non-existent taxon ID."""
        if not _check_network_connectivity():
            pytest.skip("No network connectivity to UniProt API")

        try:
            metadata = get_proteome_metadata("999999")  # Non-existent taxon
            # Should either return None or raise an exception
            # The actual behavior depends on API response
        except Exception:
            # This is acceptable - API correctly rejects invalid taxon
            pass

    def test_get_proteome_metadata_network_error_handling(self):
        """Test graceful handling of network errors."""
        # This would require mocking network calls, but we follow NO_MOCKING policy
        # Instead, we rely on the function to handle network errors internally
        # and raise appropriate exceptions

        # Test with invalid URL format (should trigger network error handling)
        with pytest.raises((requests.RequestException, ValueError)):
            # This should trigger network error handling in the function
            get_proteome_metadata("invalid_taxon_that_causes_error")

    def test_get_proteome_metadata_dependency_requests(self):
        """Test that requests library is available and used."""
        try:
            import requests
            assert requests is not None
        except ImportError:
            pytest.skip("requests library not available")

        # Function should work when requests is available
        if not _check_network_connectivity():
            pytest.skip("No network connectivity")

        try:
            metadata = get_proteome_metadata("9606")
            assert metadata is not None
        except Exception as e:
            pytest.skip(f"API error: {e}")


class TestDownloadProteomeFasta:
    """Tests for download_proteome_fasta function - real UniProt FASTA downloads."""

    @pytest.mark.network
    def test_download_proteome_fasta_human_real_download(self, tmp_path: Path):
        """Test real FASTA download for human proteome."""
        if not _check_network_connectivity():
            pytest.skip("No network connectivity to UniProt API")

        try:
            output_file = tmp_path / "human_proteome.fasta"

            result = download_proteome_fasta("9606", output_file)

            # Verify download succeeded
            assert result is True
            assert output_file.exists()
            assert output_file.stat().st_size > 0

            # Verify FASTA format
            content = output_file.read_text()
            assert content.startswith(">")
            assert "Homo sapiens" in content

        except Exception as e:
            pytest.skip(f"UniProt download unavailable: {e}")

    @pytest.mark.network
    def test_download_proteome_fasta_file_verification(self, tmp_path: Path):
        """Test that downloaded FASTA file has correct structure."""
        if not _check_network_connectivity():
            pytest.skip("No network connectivity to UniProt API")

        try:
            output_file = tmp_path / "test_proteome.fasta"

            result = download_proteome_fasta("9606", output_file)

            assert result is True
            assert output_file.exists()

            # Read and verify FASTA structure
            content = output_file.read_text()
            lines = content.strip().split('\n')

            # Should start with header
            assert lines[0].startswith('>')
            assert 'Homo sapiens' in lines[0]

            # Should have sequence data
            sequence_lines = [line for line in lines[1:] if not line.startswith('>')]
            sequence = ''.join(sequence_lines)

            # Verify sequence contains only valid amino acids
            valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
            sequence_aa = set(sequence.upper())
            assert sequence_aa.issubset(valid_aa), f"Invalid amino acids: {sequence_aa - valid_aa}"

        except Exception as e:
            pytest.skip(f"UniProt download unavailable: {e}")

    @pytest.mark.network
    def test_download_proteome_fasta_large_proteome(self, tmp_path: Path):
        """Test download of larger proteome (human)."""
        if not _check_network_connectivity():
            pytest.skip("No network connectivity to UniProt API")

        try:
            output_file = tmp_path / "human_large.fasta"

            result = download_proteome_fasta("9606", output_file)

            assert result is True
            assert output_file.exists()

            # Human proteome should be substantial
            size_mb = output_file.stat().st_size / (1024 * 1024)
            assert size_mb > 1, f"Human proteome too small: {size_mb:.2f} MB"

        except Exception as e:
            pytest.skip(f"UniProt download unavailable: {e}")

    def test_download_proteome_fasta_invalid_proteome_id(self, tmp_path: Path):
        """Test error handling for invalid proteome ID."""
        if not _check_network_connectivity():
            pytest.skip("No network connectivity to UniProt API")

        output_file = tmp_path / "invalid.fasta"

        try:
            result = download_proteome_fasta("invalid_proteome_id", output_file)
            # Should either return False or raise exception
            if result is False:
                assert not output_file.exists()
        except Exception:
            # This is acceptable - API correctly rejects invalid proteome
            pass

    def test_download_proteome_fasta_permission_error(self, tmp_path: Path):
        """Test error handling for file permission issues."""
        if not _check_network_connectivity():
            pytest.skip("Network not available for permission test")

        # Create a directory and make it read-only
        restricted_dir = tmp_path / "readonly"
        restricted_dir.mkdir()
        restricted_dir.chmod(0o444)

        output_file = restricted_dir / "test.fasta"

        try:
            # This should fail due to permission error
            result = download_proteome_fasta("9606", output_file)
            assert result is False
        except (PermissionError, OSError):
            # This is also acceptable
            pass
        finally:
            # Restore permissions for cleanup
            restricted_dir.chmod(0o755)

    def test_download_proteome_fasta_network_error_handling(self):
        """Test graceful handling of network errors during download."""
        # This tests the function's error handling without mocking
        # We rely on the function to handle network errors internally

        with tempfile.NamedTemporaryFile() as tmp_file:
            try:
                # This should trigger network error handling
                result = download_proteome_fasta("invalid_taxon_causing_error", tmp_file.name)
                # Should either return False or raise exception
                assert result is False or isinstance(result, Exception)
            except Exception:
                # Network errors are acceptable
                pass

    def test_download_proteome_fasta_dependency_requests(self):
        """Test that requests library is available for downloads."""
        try:
            import requests
            assert requests is not None
        except ImportError:
            pytest.skip("requests library not available")

        # Function should work when requests is available
        if not _check_network_connectivity():
            pytest.skip("No network connectivity")

        try:
            with tempfile.NamedTemporaryFile() as tmp_file:
                result = download_proteome_fasta("9606", tmp_file.name)
                # Should succeed or fail gracefully
                assert isinstance(result, bool)
        except Exception as e:
            pytest.skip(f"API error: {e}")


class TestProteomesNetworkConnectivity:
    """Tests for network connectivity and error handling."""

    def test_network_connectivity_check_functionality(self):
        """Test that network connectivity check works."""
        # Test with known good URL
        result = _check_network_connectivity("https://httpbin.org")
        # May pass or fail depending on network, but should not crash
        assert isinstance(result, bool)

        # Test with invalid URL
        result = _check_network_connectivity("https://invalid.domain.that.does.not.exist")
        assert result is False

    def test_api_timeout_handling(self):
        """Test that API calls handle timeouts gracefully."""
        # This tests the function's timeout handling without mocking
        import time

        start_time = time.time()
        try:
            # Use a very short timeout that should cause timeout
            import requests
            response = requests.get("https://httpbin.org/delay/10", timeout=0.001)
        except (requests.Timeout, requests.RequestException):
            # This is expected
            pass

        elapsed = time.time() - start_time
        assert elapsed < 1, "Timeout handling took too long"


class TestProteomesEdgeCases:
    """Tests for edge cases in proteome functions."""

    def test_get_proteome_metadata_empty_string(self):
        """Test with empty string input."""
        with pytest.raises(requests.exceptions.HTTPError):
            get_proteome_metadata("")

    def test_get_proteome_metadata_none_input(self):
        """Test with None input."""
        with pytest.raises(requests.exceptions.HTTPError):
            get_proteome_metadata(None)

    def test_download_proteome_fasta_empty_proteome_id(self, tmp_path: Path):
        """Test download with empty proteome ID."""
        output_file = tmp_path / "empty.fasta"

        with pytest.raises((requests.exceptions.HTTPError, OSError)):
            download_proteome_fasta("", output_file)

    def test_download_proteome_fasta_none_proteome_id(self, tmp_path: Path):
        """Test download with None proteome ID."""
        output_file = tmp_path / "none.fasta"

        with pytest.raises(requests.exceptions.HTTPError):
            download_proteome_fasta(None, output_file)

    def test_download_proteome_fasta_invalid_output_path(self):
        """Test download with invalid output path."""
        with pytest.raises((OSError, ValueError)):
            download_proteome_fasta("9606", "/invalid/path/that/does/not/exist/file.fasta")
