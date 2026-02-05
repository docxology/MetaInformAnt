"""FAST5/POD5 file reading for Oxford Nanopore signal data.

Provides real parsing of FAST5 (HDF5-based) and POD5 files produced by
Oxford Nanopore Technologies sequencers. Extracts raw electrical signal
arrays, basecalled sequences, and per-read metadata.

Optional dependencies:
    - h5py: For reading FAST5 (HDF5) files
    - pod5: For reading POD5 files
    - numpy: For signal array processing
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Sequence

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional dependency imports
try:
    import h5py  # type: ignore[import-untyped]
except ImportError:
    h5py = None  # type: ignore[assignment]

try:
    import pod5 as pod5_lib  # type: ignore[import-untyped]
except ImportError:
    pod5_lib = None  # type: ignore[assignment]

try:
    import numpy as np  # type: ignore[import-untyped]
except ImportError:
    np = None  # type: ignore[assignment]


@dataclass
class Fast5Read:
    """Represents a single read from a FAST5/POD5 file.

    Attributes:
        read_id: Unique read identifier.
        signal: Raw electrical signal array (pA values or raw ADC).
        sequence: Basecalled nucleotide sequence (if available).
        quality_string: Per-base quality string (if available).
        channel_id: Sequencing channel number.
        mux: Mux (pore selection) number.
        start_time: Read start time in samples from run start.
        duration: Read duration in samples.
        sampling_rate: ADC sampling rate in Hz.
        run_id: Sequencing run identifier.
        digitisation: ADC digitisation value for signal conversion.
        offset: Signal offset for pA conversion.
        range_value: Signal range for pA conversion.
        metadata: Additional key-value metadata.
    """

    read_id: str
    signal: Any = None  # numpy array when available
    sequence: str | None = None
    quality_string: str | None = None
    channel_id: int = 0
    mux: int = 0
    start_time: int = 0
    duration: int = 0
    sampling_rate: float = 4000.0
    run_id: str = ""
    digitisation: float = 8192.0
    offset: float = 0.0
    range_value: float = 1467.0
    metadata: dict[str, Any] = field(default_factory=dict)


def read_fast5(filepath: str | Path) -> list[Fast5Read]:
    """Parse a FAST5 or POD5 file and extract all reads.

    Supports both single-read and multi-read FAST5 files, as well as
    POD5 files. Automatically detects format from file extension.

    Args:
        filepath: Path to a .fast5 or .pod5 file.

    Returns:
        List of Fast5Read objects containing signal data and metadata.

    Raises:
        FileNotFoundError: If the file does not exist.
        ImportError: If required library (h5py or pod5) is not installed.
        ValueError: If file format is not recognized.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    suffix = filepath.suffix.lower()
    if suffix == ".pod5":
        return _read_pod5(filepath)
    elif suffix == ".fast5":
        return _read_fast5_hdf5(filepath)
    else:
        raise ValueError(f"Unrecognized file extension '{suffix}'. Expected .fast5 or .pod5")


def _read_fast5_hdf5(filepath: Path) -> list[Fast5Read]:
    """Read a FAST5 (HDF5) file, handling both single-read and multi-read formats."""
    if h5py is None:
        raise ImportError("h5py is required for reading FAST5 files. " "Install it with: uv pip install h5py")

    reads: list[Fast5Read] = []

    with h5py.File(str(filepath), "r") as f5:
        # Detect format: multi-read FAST5 has top-level groups starting with 'read_'
        top_keys = list(f5.keys())
        read_groups = [k for k in top_keys if k.startswith("read_")]

        if read_groups:
            # Multi-read FAST5
            for read_group_name in read_groups:
                read_group = f5[read_group_name]
                read_obj = _parse_fast5_read_group(read_group, read_group_name)
                reads.append(read_obj)
        else:
            # Single-read FAST5 (legacy format)
            read_obj = _parse_single_read_fast5(f5, filepath.stem)
            reads.append(read_obj)

    logger.info("Read %d reads from FAST5 file %s", len(reads), filepath.name)
    return reads


def _parse_fast5_read_group(read_group: Any, group_name: str) -> Fast5Read:
    """Parse a single read group from a multi-read FAST5 file."""
    read_id = group_name.replace("read_", "")

    # Extract read_id from attributes if available
    if "Raw" in read_group:
        raw = read_group["Raw"]
        if hasattr(raw, "attrs"):
            read_id = raw.attrs.get("read_id", read_id)
            if isinstance(read_id, bytes):
                read_id = read_id.decode("utf-8")

    signal = None
    if "Raw/Signal" in read_group:
        signal_data = read_group["Raw/Signal"][:]
        if np is not None:
            signal = np.array(signal_data, dtype=np.float32)
        else:
            signal = list(signal_data)

    # Channel information
    channel_id = 0
    mux = 0
    sampling_rate = 4000.0
    digitisation = 8192.0
    offset = 0.0
    range_value = 1467.0

    if "channel_id" in read_group:
        ch = read_group["channel_id"]
        if hasattr(ch, "attrs"):
            channel_id = int(ch.attrs.get("channel_number", 0))
            sampling_rate = float(ch.attrs.get("sampling_rate", 4000.0))
            digitisation = float(ch.attrs.get("digitisation", 8192.0))
            offset = float(ch.attrs.get("offset", 0.0))
            range_value = float(ch.attrs.get("range", 1467.0))

    # Tracking metadata
    run_id = ""
    if "tracking_id" in read_group:
        tracking = read_group["tracking_id"]
        if hasattr(tracking, "attrs"):
            run_id = tracking.attrs.get("run_id", b"")
            if isinstance(run_id, bytes):
                run_id = run_id.decode("utf-8")

    # Read attributes
    start_time = 0
    duration = 0
    if "Raw" in read_group and hasattr(read_group["Raw"], "attrs"):
        raw_attrs = read_group["Raw"].attrs
        start_time = int(raw_attrs.get("start_time", 0))
        duration = int(raw_attrs.get("duration", 0))

    # Basecalls if available
    sequence = None
    quality_string = None
    basecall_groups = [k for k in read_group.keys() if k.startswith("Analyses")]
    if "Analyses" in read_group:
        analyses = read_group["Analyses"]
        # Look for Basecall_1D groups
        bc_groups = sorted([k for k in analyses.keys() if "Basecall" in k], reverse=True)
        for bc_name in bc_groups:
            bc_group = analyses[bc_name]
            fastq_path = None
            for sub in ["BaseCalled_template/Fastq", "BaseCalled_complement/Fastq"]:
                if sub in bc_group:
                    fastq_path = sub
                    break
            if fastq_path is not None:
                raw_fastq = bc_group[fastq_path][()].decode("utf-8")
                lines = raw_fastq.strip().split("\n")
                if len(lines) >= 4:
                    sequence = lines[1]
                    quality_string = lines[3]
                break

    metadata: dict[str, Any] = {
        "source_format": "fast5_multi",
        "digitisation": digitisation,
        "range": range_value,
    }

    return Fast5Read(
        read_id=str(read_id),
        signal=signal,
        sequence=sequence,
        quality_string=quality_string,
        channel_id=channel_id,
        mux=mux,
        start_time=start_time,
        duration=duration,
        sampling_rate=sampling_rate,
        run_id=run_id,
        digitisation=digitisation,
        offset=offset,
        range_value=range_value,
        metadata=metadata,
    )


def _parse_single_read_fast5(f5: Any, fallback_id: str) -> Fast5Read:
    """Parse a single-read (legacy) FAST5 file."""
    read_id = fallback_id
    signal = None
    sequence = None
    quality_string = None
    channel_id = 0
    mux = 0
    start_time = 0
    duration = 0
    sampling_rate = 4000.0
    run_id = ""
    digitisation = 8192.0
    offset = 0.0
    range_value = 1467.0

    # Traverse HDF5 to find Raw signal
    if "Raw" in f5 and "Reads" in f5["Raw"]:
        reads_grp = f5["Raw"]["Reads"]
        read_keys = list(reads_grp.keys())
        if read_keys:
            first_read = reads_grp[read_keys[0]]
            if "Signal" in first_read:
                signal_data = first_read["Signal"][:]
                if np is not None:
                    signal = np.array(signal_data, dtype=np.float32)
                else:
                    signal = list(signal_data)
            if hasattr(first_read, "attrs"):
                read_id = first_read.attrs.get("read_id", fallback_id)
                if isinstance(read_id, bytes):
                    read_id = read_id.decode("utf-8")
                start_time = int(first_read.attrs.get("start_time", 0))
                duration = int(first_read.attrs.get("duration", 0))
                mux = int(first_read.attrs.get("start_mux", 0))

    # Channel info
    if "UniqueGlobalKey/channel_id" in f5:
        ch = f5["UniqueGlobalKey/channel_id"]
        if hasattr(ch, "attrs"):
            channel_id = int(ch.attrs.get("channel_number", 0))
            sampling_rate = float(ch.attrs.get("sampling_rate", 4000.0))
            digitisation = float(ch.attrs.get("digitisation", 8192.0))
            offset = float(ch.attrs.get("offset", 0.0))
            range_value = float(ch.attrs.get("range", 1467.0))

    # Tracking
    if "UniqueGlobalKey/tracking_id" in f5:
        tracking = f5["UniqueGlobalKey/tracking_id"]
        if hasattr(tracking, "attrs"):
            run_id_val = tracking.attrs.get("run_id", b"")
            run_id = run_id_val.decode("utf-8") if isinstance(run_id_val, bytes) else str(run_id_val)

    # Basecalls
    if "Analyses" in f5:
        analyses = f5["Analyses"]
        bc_groups = sorted([k for k in analyses.keys() if "Basecall" in k], reverse=True)
        for bc_name in bc_groups:
            bc_group = analyses[bc_name]
            for sub in ["BaseCalled_template/Fastq", "BaseCalled_complement/Fastq"]:
                if sub in bc_group:
                    raw_fastq = bc_group[sub][()].decode("utf-8")
                    lines = raw_fastq.strip().split("\n")
                    if len(lines) >= 4:
                        sequence = lines[1]
                        quality_string = lines[3]
                    break
            if sequence is not None:
                break

    return Fast5Read(
        read_id=str(read_id),
        signal=signal,
        sequence=sequence,
        quality_string=quality_string,
        channel_id=channel_id,
        mux=mux,
        start_time=start_time,
        duration=duration,
        sampling_rate=sampling_rate,
        run_id=run_id,
        digitisation=digitisation,
        offset=offset,
        range_value=range_value,
        metadata={"source_format": "fast5_single"},
    )


def _read_pod5(filepath: Path) -> list[Fast5Read]:
    """Read a POD5 file and convert to Fast5Read objects."""
    if pod5_lib is None:
        raise ImportError("pod5 is required for reading POD5 files. " "Install it with: uv pip install pod5")

    reads: list[Fast5Read] = []

    with pod5_lib.Reader(filepath) as reader:
        for record in reader.reads():
            signal = None
            if np is not None:
                signal = np.array(record.signal, dtype=np.float32)
            else:
                signal = list(record.signal)

            calibration = record.calibration
            run_info = record.run_info
            pore = record.pore

            read_obj = Fast5Read(
                read_id=str(record.read_id),
                signal=signal,
                sequence=None,  # POD5 files typically don't contain basecalls
                quality_string=None,
                channel_id=pore.channel if pore else 0,
                mux=pore.well if pore else 0,
                start_time=record.start_sample if hasattr(record, "start_sample") else 0,
                duration=record.num_samples if hasattr(record, "num_samples") else 0,
                sampling_rate=run_info.sample_rate if run_info else 4000.0,
                run_id=str(run_info.run_id) if run_info else "",
                digitisation=calibration.digitisation if hasattr(calibration, "digitisation") else 8192.0,
                offset=calibration.offset if calibration else 0.0,
                range_value=calibration.range if hasattr(calibration, "range") else 1467.0,
                metadata={"source_format": "pod5"},
            )
            reads.append(read_obj)

    logger.info("Read %d reads from POD5 file %s", len(reads), filepath.name)
    return reads


def extract_signal(fast5_data: Fast5Read) -> Any:
    """Extract raw electrical signal array from a parsed read.

    Converts raw ADC values to picoampere (pA) values using the calibration
    parameters stored in the read metadata.

    The conversion formula is:
        pA = (raw_adc + offset) * (range / digitisation)

    Args:
        fast5_data: A parsed Fast5Read object.

    Returns:
        Numpy array of signal values in picoamperes, or a plain list if
        numpy is not available. Returns None if no signal data exists.
    """
    if fast5_data.signal is None:
        return None

    # Convert raw ADC to pA using calibration: pA = (raw + offset) * (range / digitisation)
    if np is not None:
        raw = np.asarray(fast5_data.signal, dtype=np.float64)
        scale = fast5_data.range_value / fast5_data.digitisation if fast5_data.digitisation != 0 else 1.0
        pa_signal = (raw + fast5_data.offset) * scale
        return pa_signal.astype(np.float32)
    else:
        scale = fast5_data.range_value / fast5_data.digitisation if fast5_data.digitisation != 0 else 1.0
        return [(v + fast5_data.offset) * scale for v in fast5_data.signal]


def extract_basecalls(fast5_data: Fast5Read) -> tuple[str | None, str | None]:
    """Extract basecalled sequence and quality string from a parsed read.

    Args:
        fast5_data: A parsed Fast5Read object.

    Returns:
        Tuple of (sequence, quality_string). Either may be None if basecalls
        are not available in the read.
    """
    return fast5_data.sequence, fast5_data.quality_string


def get_read_metadata(fast5_data: Fast5Read) -> dict[str, Any]:
    """Extract comprehensive metadata from a parsed read.

    Returns a dictionary with channel, mux, duration, start time, sampling
    rate, run information, and calibration parameters.

    Args:
        fast5_data: A parsed Fast5Read object.

    Returns:
        Dictionary containing all available read metadata.
    """
    duration_seconds = 0.0
    if fast5_data.sampling_rate > 0 and fast5_data.duration > 0:
        duration_seconds = fast5_data.duration / fast5_data.sampling_rate

    signal_length = 0
    if fast5_data.signal is not None:
        if np is not None and hasattr(fast5_data.signal, "__len__"):
            signal_length = len(fast5_data.signal)
        elif isinstance(fast5_data.signal, (list, tuple)):
            signal_length = len(fast5_data.signal)

    return {
        "read_id": fast5_data.read_id,
        "channel_id": fast5_data.channel_id,
        "mux": fast5_data.mux,
        "start_time_samples": fast5_data.start_time,
        "duration_samples": fast5_data.duration,
        "duration_seconds": duration_seconds,
        "sampling_rate": fast5_data.sampling_rate,
        "run_id": fast5_data.run_id,
        "has_signal": fast5_data.signal is not None,
        "signal_length": signal_length,
        "has_basecalls": fast5_data.sequence is not None,
        "sequence_length": len(fast5_data.sequence) if fast5_data.sequence else 0,
        "calibration": {
            "digitisation": fast5_data.digitisation,
            "offset": fast5_data.offset,
            "range": fast5_data.range_value,
        },
        **fast5_data.metadata,
    }
