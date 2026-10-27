"""The metabolomics I/O module must not claim file formats it cannot read.

The module docstring of ``metainformant.metabolomics.io.formats`` used to promise "limited readers for mzML and
mzXML formats", but the module only reads CSV intensity tables and MGF spectra. These tests pin the set of readers
and writers that exist and check that the docstring says the same thing, so the two cannot drift apart again.

Uses the real module; no test doubles (real-implementation policy).
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest

from metainformant.metabolomics.io import formats

NEGATION = re.compile(r"\b(not|no|cannot|unsupported)\b", re.IGNORECASE)


def _sentences(text: str) -> list[str]:
    return re.split(r"(?<=[.;])\s+", " ".join(text.split()))


def test_the_module_reads_and_writes_exactly_csv_and_mgf() -> None:
    readers = sorted(name for name in dir(formats) if name.startswith("read_"))
    writers = sorted(name for name in dir(formats) if name.startswith("write_"))

    # adding a reader (for example for mzML) means updating the docstring and this expectation together
    assert readers == ["read_csv", "read_mgf"]
    assert writers == ["write_csv", "write_mgf"]


def test_docstring_names_the_formats_that_exist() -> None:
    doc = formats.__doc__ or ""

    assert "CSV" in doc
    assert "MGF" in doc
    for function in ("read_csv", "write_csv", "read_mgf", "write_mgf"):
        assert function in doc


def test_docstring_does_not_claim_support_for_mzml_or_mzxml() -> None:
    doc = formats.__doc__ or ""

    mentions = [sentence for sentence in _sentences(doc) if re.search(r"mzml|mzxml", sentence, re.IGNORECASE)]
    # every mention must say the formats are not supported
    assert mentions, "the docstring should tell readers that mzML and mzXML are not supported"
    for sentence in mentions:
        assert NEGATION.search(sentence), f"claims support for a format with no reader: {sentence!r}"


def test_an_mzml_file_is_not_silently_accepted_as_mgf(tmp_path: Path) -> None:
    # the documented limit holds in practice: an mzML document yields no spectra rather than made-up ones
    mzml = tmp_path / "run.mzML"
    mzml.write_text('<?xml version="1.0"?><mzML><run id="r1"><spectrumList count="0"/></run></mzML>\n')

    assert formats.read_mgf(mzml) == []


def test_documented_formats_round_trip(tmp_path: Path) -> None:
    spectrum = formats.MassSpectrum(
        scan_number=1,
        retention_time=12.5,
        ms_level=2,
        mz_array=np.array([100.0, 200.5]),
        intensity_array=np.array([10.0, 20.0]),
        precursor_mz=300.25,
        total_ion_current=30.0,
    )
    path = tmp_path / "spectra.mgf"
    formats.write_mgf([spectrum], path)

    (reloaded,) = formats.read_mgf(path)

    assert reloaded.scan_number == 1
    assert reloaded.precursor_mz == pytest.approx(300.25)
    assert reloaded.mz_array.tolist() == pytest.approx([100.0, 200.5])
