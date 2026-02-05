"""Format conversion utilities for long-read sequencing data.

Provides conversion between FAST5, POD5, FASTQ, and PAF formats. All
conversions use real file parsing with no placeholder implementations.

Optional dependencies:
    - h5py: For FAST5 reading/writing
    - pod5: For POD5 reading
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Any, Sequence

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

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


def fast5_to_fastq(
    fast5_path: str | Path,
    output_path: str | Path,
    min_quality: float = 0.0,
    gzip_output: bool = False,
) -> int:
    """Convert a FAST5 file to FASTQ format.

    Extracts basecalled sequences and quality scores from FAST5 files
    and writes them in standard FASTQ format.

    Args:
        fast5_path: Path to input FAST5 file.
        output_path: Path for output FASTQ file.
        min_quality: Minimum mean quality score to include a read.
        gzip_output: Whether to gzip-compress the output.

    Returns:
        Number of reads written to the output file.

    Raises:
        FileNotFoundError: If the FAST5 file does not exist.
        ImportError: If h5py is not installed.
    """
    from .fast5 import read_fast5

    fast5_path = Path(fast5_path)
    output_path = Path(output_path)

    if not fast5_path.exists():
        raise FileNotFoundError(f"FAST5 file not found: {fast5_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    reads = read_fast5(fast5_path)
    written = 0

    opener = gzip.open if gzip_output else open
    mode = "wt" if gzip_output else "w"

    with opener(str(output_path), mode) as fout:
        for read_obj in reads:
            if read_obj.sequence is None:
                logger.debug("Skipping read %s: no basecalls", read_obj.read_id)
                continue

            # Calculate mean quality if filtering
            if min_quality > 0.0 and read_obj.quality_string:
                mean_q = _mean_phred_quality(read_obj.quality_string)
                if mean_q < min_quality:
                    continue

            # Construct quality string
            qual = read_obj.quality_string or "!" * len(read_obj.sequence)

            # Write FASTQ entry
            fout.write(f"@{read_obj.read_id}\n")
            fout.write(f"{read_obj.sequence}\n")
            fout.write("+\n")
            fout.write(f"{qual}\n")
            written += 1

    logger.info("Wrote %d reads to FASTQ %s", written, output_path.name)
    return written


def convert_pod5_to_fast5(
    pod5_path: str | Path,
    output_path: str | Path,
) -> int:
    """Convert a POD5 file to FAST5 format.

    Reads signal data from POD5 format and writes it to multi-read FAST5
    (HDF5) format for compatibility with legacy tools.

    Args:
        pod5_path: Path to input POD5 file.
        output_path: Path for output FAST5 file.

    Returns:
        Number of reads converted.

    Raises:
        FileNotFoundError: If the POD5 file does not exist.
        ImportError: If pod5 and/or h5py are not installed.
    """
    if pod5_lib is None:
        raise ImportError(
            "pod5 is required for reading POD5 files. "
            "Install it with: uv pip install pod5"
        )
    if h5py is None:
        raise ImportError(
            "h5py is required for writing FAST5 files. "
            "Install it with: uv pip install h5py"
        )

    pod5_path = Path(pod5_path)
    output_path = Path(output_path)

    if not pod5_path.exists():
        raise FileNotFoundError(f"POD5 file not found: {pod5_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    count = 0

    with pod5_lib.Reader(pod5_path) as reader, h5py.File(str(output_path), "w") as f5out:
        for record in reader.reads():
            read_id = str(record.read_id)
            group_name = f"read_{read_id}"
            read_group = f5out.create_group(group_name)

            # Write raw signal
            raw_group = read_group.create_group("Raw")
            signal_data = record.signal
            if np is not None:
                signal_data = np.array(signal_data, dtype=np.int16)
            raw_group.create_dataset("Signal", data=signal_data)

            # Write read attributes
            raw_group.attrs["read_id"] = read_id
            raw_group.attrs["duration"] = record.num_samples if hasattr(record, "num_samples") else len(signal_data)
            raw_group.attrs["start_time"] = record.start_sample if hasattr(record, "start_sample") else 0

            # Write channel info
            ch_group = read_group.create_group("channel_id")
            pore = record.pore
            if pore:
                ch_group.attrs["channel_number"] = pore.channel
            calibration = record.calibration
            if calibration:
                ch_group.attrs["offset"] = calibration.offset
            run_info = record.run_info
            if run_info:
                ch_group.attrs["sampling_rate"] = run_info.sample_rate

            # Write tracking info
            tracking = read_group.create_group("tracking_id")
            if run_info:
                tracking.attrs["run_id"] = str(run_info.run_id)

            count += 1

    logger.info("Converted %d reads from POD5 to FAST5: %s", count, output_path.name)
    return count


def write_paf(
    alignments: Sequence[dict[str, Any]],
    output_path: str | Path,
) -> int:
    """Write alignments in PAF (Pairwise Alignment Format).

    PAF is a TAB-delimited text format for storing pairwise alignments,
    commonly used by minimap2. Each line has 12 mandatory fields.

    PAF format columns:
        1.  Query name
        2.  Query length
        3.  Query start (0-based)
        4.  Query end (0-based)
        5.  Strand (+/-)
        6.  Target name
        7.  Target length
        8.  Target start (0-based)
        9.  Target end (0-based)
        10. Number of matching bases
        11. Alignment block length
        12. Mapping quality

    Args:
        alignments: Sequence of alignment dictionaries. Each must contain:
            - query_name (str)
            - query_length (int)
            - query_start (int)
            - query_end (int)
            - strand (str: "+" or "-")
            - target_name (str)
            - target_length (int)
            - target_start (int)
            - target_end (int)
            - matches (int)
            - block_length (int)
            - mapping_quality (int)
            Optionally: tags (dict of key-value pairs for extra columns)
        output_path: Path for the output PAF file.

    Returns:
        Number of alignment records written.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    written = 0

    with open(output_path, "w") as fout:
        for aln in alignments:
            # Mandatory fields
            fields = [
                str(aln["query_name"]),
                str(aln["query_length"]),
                str(aln["query_start"]),
                str(aln["query_end"]),
                str(aln.get("strand", "+")),
                str(aln["target_name"]),
                str(aln["target_length"]),
                str(aln["target_start"]),
                str(aln["target_end"]),
                str(aln.get("matches", 0)),
                str(aln.get("block_length", 0)),
                str(aln.get("mapping_quality", 255)),
            ]

            # Optional SAM-like tags
            tags = aln.get("tags", {})
            for key, value in tags.items():
                if isinstance(value, int):
                    fields.append(f"{key}:i:{value}")
                elif isinstance(value, float):
                    fields.append(f"{key}:f:{value:.6f}")
                else:
                    fields.append(f"{key}:Z:{value}")

            fout.write("\t".join(fields) + "\n")
            written += 1

    logger.info("Wrote %d alignments to PAF %s", written, output_path.name)
    return written


def _mean_phred_quality(quality_string: str) -> float:
    """Calculate mean Phred quality score from an ASCII quality string.

    Phred score = ord(char) - 33 (Sanger/Illumina 1.8+ encoding).

    Args:
        quality_string: ASCII-encoded quality string.

    Returns:
        Mean Phred quality score.
    """
    if not quality_string:
        return 0.0
    scores = [ord(c) - 33 for c in quality_string]
    return sum(scores) / len(scores)
