"""Tests for longread IO functions not covered by test_longread.py.

Tests read_long_read_bam error handling, convert_pod5_to_fast5 error cases,
and fast5_to_fastq error cases. All tests use real implementations -- NO MOCKING.
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

import pytest

from metainformant.longread.io.bam import read_long_read_bam
from metainformant.longread.io.formats import convert_pod5_to_fast5, fast5_to_fastq


class TestReadLongReadBam:
    """Tests for read_long_read_bam function -- error handling paths."""

    def test_missing_file_raises_file_not_found(self) -> None:
        """read_long_read_bam raises FileNotFoundError for nonexistent BAM."""
        with pytest.raises((FileNotFoundError, ImportError)):
            read_long_read_bam("/nonexistent/path/to/file.bam")

    def test_missing_file_path_object(self) -> None:
        """read_long_read_bam accepts Path objects and raises FileNotFoundError."""
        with pytest.raises((FileNotFoundError, ImportError)):
            read_long_read_bam(Path("/nonexistent/path/reads.bam"))


class TestConvertPod5ToFast5:
    """Tests for convert_pod5_to_fast5 -- error and edge cases."""

    def test_missing_pod5_file_raises(self) -> None:
        """convert_pod5_to_fast5 raises FileNotFoundError or ImportError for missing file."""
        with pytest.raises((FileNotFoundError, ImportError)):
            convert_pod5_to_fast5(
                "/nonexistent/path/reads.pod5",
                "/tmp/output.fast5",
            )

    def test_import_error_without_pod5(self) -> None:
        """convert_pod5_to_fast5 raises ImportError when pod5 library is missing."""
        # If pod5 is not installed, we expect ImportError.
        # If it IS installed but file is missing, we expect FileNotFoundError.
        # Both are acceptable.
        with pytest.raises((ImportError, FileNotFoundError)):
            convert_pod5_to_fast5(
                "/nonexistent/reads.pod5",
                "/tmp/output.fast5",
            )


class TestFast5ToFastq:
    """Tests for fast5_to_fastq -- error and edge cases."""

    def test_missing_fast5_raises_file_not_found(self) -> None:
        """fast5_to_fastq raises FileNotFoundError for nonexistent FAST5."""
        with pytest.raises(FileNotFoundError, match="FAST5 file not found"):
            fast5_to_fastq(
                "/nonexistent/path/reads.fast5",
                "/tmp/output.fastq",
            )

    def test_missing_fast5_path_object(self) -> None:
        """fast5_to_fastq raises FileNotFoundError for Path to nonexistent file."""
        with pytest.raises(FileNotFoundError, match="FAST5 file not found"):
            fast5_to_fastq(
                Path("/nonexistent/reads.fast5"),
                Path("/tmp/output.fastq"),
            )

    def test_fast5_to_fastq_creates_parent_dirs(self, tmp_path: Path) -> None:
        """fast5_to_fastq raises FileNotFoundError before trying to create output dirs."""
        output_path = tmp_path / "subdir" / "output.fastq"
        with pytest.raises(FileNotFoundError, match="FAST5 file not found"):
            fast5_to_fastq(
                "/nonexistent/reads.fast5",
                output_path,
            )


class TestFormatsMeanPhredQuality:
    """Tests for the _mean_phred_quality internal helper in formats.py."""

    def test_mean_phred_quality_known_values(self) -> None:
        """Verify mean Phred quality calculation from ASCII quality strings."""
        from metainformant.longread.io.formats import _mean_phred_quality

        # Q20 = chr(20+33) = chr(53) = '5'
        q20_str = "5" * 10
        result = _mean_phred_quality(q20_str)
        assert result == pytest.approx(20.0)

    def test_mean_phred_quality_empty_string(self) -> None:
        """Empty quality string returns 0.0."""
        from metainformant.longread.io.formats import _mean_phred_quality

        assert _mean_phred_quality("") == 0.0

    def test_mean_phred_quality_mixed(self) -> None:
        """Mixed quality string produces correct mean."""
        from metainformant.longread.io.formats import _mean_phred_quality

        # Q10 = chr(43) = '+', Q30 = chr(63) = '?'
        q_str = chr(10 + 33) + chr(30 + 33)  # Q10 + Q30
        result = _mean_phred_quality(q_str)
        assert result == pytest.approx(20.0)
