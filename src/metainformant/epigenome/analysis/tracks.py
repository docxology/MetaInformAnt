"""Genomic track processing and analysis.

This module provides functionality to load, manipulate, and analyze
various genomic track formats including BED, BEDgraph, bigWig, and
custom track files for epigenomic data visualization and analysis.
"""

from __future__ import annotations

import heapq
import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

from metainformant.core import io
from metainformant.core.data import validation
from metainformant.core.utils import errors, logging

if TYPE_CHECKING:
    import pandas as pd

logger = logging.get_logger(__name__)


class GenomicTrack:
    """Represents a genomic track with chromosome-based data."""

    def __init__(self, name: str = "", description: str = ""):
        """Initialize a genomic track.

        Args:
            name: Track name
            description: Track description
        """
        self.name = name
        self.description = description
        self.data: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
        self.metadata: Dict[str, Any] = {}

    def add_feature(self, chromosome: str, start: int, end: int, value: float = 0.0, **kwargs: Any) -> None:
        """Add a feature to the track.

        Args:
            chromosome: Chromosome name
            start: Feature start position
            end: Feature end position
            value: Feature value
            **kwargs: Additional feature attributes
        """
        feature = {"start": start, "end": end, "value": value, **kwargs}
        self.data[chromosome].append(feature)

    def get_features(
        self, chromosome: str, start: Optional[int] = None, end: Optional[int] = None
    ) -> List[Dict[str, Any]]:
        """Get features from a genomic region.

        Args:
            chromosome: Chromosome name
            start: Region start position (optional)
            end: Region end position (optional)

        Returns:
            List of features in the region
        """
        if chromosome not in self.data:
            return []

        features = self.data[chromosome]

        if start is None and end is None:
            return features

        # Filter by region
        filtered = []
        for feature in features:
            if start is not None and feature["end"] <= start:
                continue
            if end is not None and feature["start"] >= end:
                continue
            filtered.append(feature)

        return filtered

    def get_chromosomes(self) -> List[str]:
        """Get list of chromosomes in the track."""
        return list(self.data.keys())

    def get_total_features(self) -> int:
        """Get total number of features across all chromosomes."""
        return sum(len(features) for features in self.data.values())

    def merge_overlapping_features(self, chromosome: str, max_distance: int = 0) -> None:
        """Merge overlapping features on a chromosome.

        Args:
            chromosome: Chromosome name
            max_distance: Maximum distance between features to merge
        """
        if chromosome not in self.data:
            return

        features = sorted(self.data[chromosome], key=lambda x: x["start"])

        if not features:
            return

        merged = [features[0]]

        for feature in features[1:]:
            last = merged[-1]

            # Check if features overlap or are within max_distance
            if feature["start"] <= last["end"] + max_distance:
                # Merge features
                last["end"] = max(last["end"], feature["end"])
                last["value"] = max(last["value"], feature["value"])  # Use max value
                # Merge other attributes if they exist
                for key, value in feature.items():
                    if key not in ["start", "end", "value"] and key not in last:
                        last[key] = value
            else:
                merged.append(feature)

        self.data[chromosome] = merged

    def normalize_values(self, method: str = "minmax") -> None:
        """Normalize feature values.

        Args:
            method: Normalization method ("minmax", "zscore", "robust")

        Raises:
            ValueError: If the method is not one of the three above (an unknown name used to change nothing and still
                log that the values had been normalized).
        """
        if method not in ("minmax", "zscore", "robust"):
            raise ValueError(f"Unsupported normalization method {method!r}; use 'minmax', 'zscore' or 'robust'")

        all_values: List[float] = []
        for features in self.data.values():
            all_values.extend(f["value"] for f in features)

        if not all_values:
            return

        if method == "minmax":
            min_val, max_val = min(all_values), max(all_values)
            if max_val > min_val:
                for features in self.data.values():
                    for feature in features:
                        feature["value"] = (feature["value"] - min_val) / (max_val - min_val)

        elif method == "zscore":
            if len(all_values) > 1:
                mean_val = statistics.mean(all_values)
                std_val = statistics.stdev(all_values)
                if std_val > 0:
                    for features in self.data.values():
                        for feature in features:
                            feature["value"] = (feature["value"] - mean_val) / std_val

        elif method == "robust":
            # Use median and MAD (median absolute deviation)
            median_val = statistics.median(all_values)
            mad = statistics.median(abs(v - median_val) for v in all_values)
            if mad > 0:
                for features in self.data.values():
                    for feature in features:
                        feature["value"] = (feature["value"] - median_val) / mad

        logger.info(f"Normalized track values using {method} method")

    def to_dataframe(self) -> "pd.DataFrame":
        """Convert track data to a pandas DataFrame.

        Returns:
            DataFrame with columns: chrom, start, end, value
        """
        import pandas as pd

        rows = []
        for chromosome, features in self.data.items():
            for feature in features:
                rows.append(
                    {
                        "chrom": chromosome,
                        "start": feature["start"],
                        "end": feature["end"],
                        "value": feature["value"],
                    }
                )
        return pd.DataFrame(rows)


def load_bed_track(path: str | Path, name: str = "", description: str = "") -> GenomicTrack:
    """Load a BED format track.

    Args:
        path: Path to BED file
        name: Track name
        description: Track description

    Returns:
        GenomicTrack object
    """
    path = validation.validate_path_exists(Path(path))

    logger.info(f"Loading BED track from {path}")

    track = GenomicTrack(name=name, description=description)

    try:
        with io.open_text_auto(path, "rt") as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                    continue

                parts = line.split("\t")
                if len(parts) < 3:
                    logger.warning(f"Skipping malformed BED line {line_num}: insufficient columns")
                    continue

                try:
                    chromosome = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    name_field = parts[3] if len(parts) > 3 else ""
                    score = float(parts[4]) if len(parts) > 4 and parts[4] != "." else 0.0
                    strand = parts[5] if len(parts) > 5 else "."

                    track.add_feature(
                        chromosome=chromosome, start=start, end=end, value=score, name=name_field, strand=strand
                    )

                except (ValueError, IndexError) as e:
                    logger.warning(f"Error parsing BED line {line_num}: {e}")
                    continue

    except Exception as e:
        logger.error(f"Error loading BED track from {path}: {e}")
        raise errors.IOError(f"Failed to load BED track: {e}") from e

    logger.info(f"Loaded BED track with {track.get_total_features()} features")
    return track


def load_bedgraph_track(path: str | Path, name: str = "", description: str = "") -> GenomicTrack:
    """Load a BEDgraph format track.

    Args:
        path: Path to BEDgraph file
        name: Track name
        description: Track description

    Returns:
        GenomicTrack object
    """
    path = validation.validate_path_exists(Path(path))

    logger.info(f"Loading BEDgraph track from {path}")

    track = GenomicTrack(name=name, description=description)

    try:
        with io.open_text_auto(path, "rt") as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                    continue

                parts = line.split("\t")
                if len(parts) < 4:
                    logger.warning(f"Skipping malformed BEDgraph line {line_num}: insufficient columns")
                    continue

                try:
                    chromosome = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    value = float(parts[3])

                    track.add_feature(chromosome=chromosome, start=start, end=end, value=value)

                except (ValueError, IndexError) as e:
                    logger.warning(f"Error parsing BEDgraph line {line_num}: {e}")
                    continue

    except Exception as e:
        logger.error(f"Error loading BEDgraph track from {path}: {e}")
        raise errors.IOError(f"Failed to load BEDgraph track: {e}") from e

    logger.info(f"Loaded BEDgraph track with {track.get_total_features()} features")
    return track


def read_bedgraph(path: str | Path, name: str = "", description: str = "") -> "pd.DataFrame":
    """Read a BEDgraph file and return as DataFrame.

    Args:
        path: Path to BEDgraph file
        name: Track name (ignored)
        description: Track description (ignored)

    Returns:
        DataFrame with columns: chrom, start, end, value
    """
    track = load_bedgraph_track(path, name, description)
    return track.to_dataframe()


def load_bigwig_track(
    path: str | Path,
    name: str = "",
    description: str = "",
    region: tuple[str, int, int] | None = None,
) -> GenomicTrack:
    """Load a bigWig format track.

    Args:
        path: Path to bigWig file
        name: Track name
        description: Track description
        region: Optional tuple of (chrom, start, end) to load specific region

    Returns:
        GenomicTrack object with loaded data

    Raises:
        ImportError: If pyBigWig library is not available
        FileNotFoundError: If the file does not exist
    """
    path = validation.validate_path_exists(Path(path))

    # Try to import pyBigWig
    try:
        import pyBigWig
    except ImportError:
        raise ImportError(
            "pyBigWig library is required to load bigWig files. "
            "Install with: pip install pyBigWig (or: uv pip install pyBigWig)"
        )

    logger.info(f"Loading bigWig track from {path}")

    track = GenomicTrack(name=name or path.stem, description=description)
    track.metadata["format"] = "bigwig"
    track.metadata["source_file"] = str(path)

    # Open bigWig file
    bw = pyBigWig.open(str(path))

    try:
        # Get chromosome information
        chroms = bw.chroms()
        track.metadata["chromosomes"] = list(chroms.keys())
        track.metadata["total_length"] = sum(chroms.values())

        if region:
            # Load specific region
            chrom, start, end = region
            if chrom in chroms:
                values = bw.values(chrom, start, end)
                if values:
                    # Store as genomic feature
                    track.add_feature(
                        chromosome=chrom,
                        start=start,
                        end=end,
                        value=(
                            sum(v for v in values if v is not None) / len([v for v in values if v is not None])
                            if any(v is not None for v in values)
                            else 0
                        ),
                        values=values,
                    )
        else:
            # Load summary statistics for each chromosome
            for chrom, length in chroms.items():
                # Get stats for each chromosome
                try:
                    stats = bw.stats(chrom, 0, length)
                    if stats and stats[0] is not None:
                        track.add_feature(chromosome=chrom, start=0, end=length, value=stats[0])
                except Exception as e:
                    logger.warning(f"Could not load stats for {chrom}: {e}")
                    continue

        track.metadata["loaded"] = True
        logger.info(f"Loaded bigWig track with {track.get_total_features()} chromosome entries")

    finally:
        bw.close()

    return track


def save_bed_track(track: GenomicTrack, path: str | Path) -> None:
    """Save a GenomicTrack to BED format.

    Args:
        path: Output path
        track: GenomicTrack to save
    """
    path = Path(path)

    logger.info(f"Saving track to BED format: {path}")

    with open(path, "w") as f:
        # Write track header
        if track.name:
            f.write(f'track name="{track.name}" description="{track.description}"\n')

        for chromosome in sorted(track.data.keys()):
            features = sorted(track.data[chromosome], key=lambda x: x["start"])

            for feature in features:
                name = feature.get("name", f'feature_{feature["start"]}_{feature["end"]}')
                score = feature.get("value", 0)
                strand = feature.get("strand", ".")

                bed_line = "\t".join([chromosome, str(feature["start"]), str(feature["end"]), name, str(score), strand])
                f.write(bed_line + "\n")

    logger.info(f"Saved {track.get_total_features()} features to {path}")


def save_bedgraph_track(track: GenomicTrack, path: str | Path) -> None:
    """Save a GenomicTrack to BEDgraph format.

    Args:
        path: Output path
        track: GenomicTrack to save
    """
    path = Path(path)

    logger.info(f"Saving track to BEDgraph format: {path}")

    with open(path, "w") as f:
        # Write track header
        if track.name:
            f.write(f'track type=bedGraph name="{track.name}" description="{track.description}"\n')

        for chromosome in sorted(track.data.keys()):
            features = sorted(track.data[chromosome], key=lambda x: x["start"])

            for feature in features:
                bedgraph_line = "\t".join(
                    [chromosome, str(feature["start"]), str(feature["end"]), str(feature["value"])]
                )
                f.write(bedgraph_line + "\n")

    logger.info(f"Saved {track.get_total_features()} features to {path}")


def _coverage_profile(features: List[Dict[str, Any]]) -> List[Tuple[int, int, float]]:
    """Piecewise-constant coverage of a set of features.

    Args:
        features: Features of one track on one chromosome (dicts with ``start``, ``end`` and ``value``).
            Features with no extent (``end <= start``) are ignored.

    Returns:
        Sorted, non-overlapping ``(start, end, value)`` segments that tile the union of the features. A segment's
        value is the highest value among the features that cover it, and touching segments with the same value
        are joined.
    """
    spans = sorted((f["start"], f["end"], f["value"]) for f in features if f["end"] > f["start"])
    boundaries = sorted({point for start, end, _ in spans for point in (start, end)})

    segments: List[Tuple[int, int, float]] = []
    active: List[Tuple[float, int]] = []  # heap of (-value, end): the covering feature with the highest value first
    index = 0
    for left, right in zip(boundaries, boundaries[1:]):
        while index < len(spans) and spans[index][0] <= left:
            heapq.heappush(active, (-spans[index][2], spans[index][1]))
            index += 1
        while active and active[0][1] <= left:
            heapq.heappop(active)
        if not active:
            continue  # a gap between features

        # Every feature endpoint is a boundary, so a feature covers all of [left, right) or none of it
        value = -active[0][0]
        if segments and segments[-1][1] == left and segments[-1][2] == value:
            segments[-1] = (segments[-1][0], right, value)
        else:
            segments.append((left, right, value))

    return segments


def _intersect_profiles(
    first: List[Tuple[int, int, float]], second: List[Tuple[int, int, float]]
) -> List[Tuple[int, int, float]]:
    """Intersect two coverage profiles; each shared segment takes the lower of the two values."""
    result: List[Tuple[int, int, float]] = []
    i = j = 0
    while i < len(first) and j < len(second):
        start = max(first[i][0], second[j][0])
        end = min(first[i][1], second[j][1])
        if start < end:
            value = min(first[i][2], second[j][2])
            if result and result[-1][1] == start and result[-1][2] == value:
                result[-1] = (result[-1][0], end, value)
            else:
                result.append((start, end, value))
        if first[i][1] < second[j][1]:
            i += 1
        else:
            j += 1
    return result


def merge_tracks(tracks: List[GenomicTrack], operation: str = "union") -> GenomicTrack:
    """Merge multiple genomic tracks.

    The input tracks are not modified. Chromosomes are processed in sorted order.

    - ``"union"``: overlapping features of any track are merged into one feature spanning them, with the highest
      value.
    - ``"mean"`` / ``"max"``: features that overlap, directly or through a chain of overlapping features, form a
      group that becomes one feature spanning the group, with the mean or the highest value of its features.
    - ``"intersection"``: only the regions covered by every track are kept, so a chromosome missing from any track
      has none. Within a track, the value over a region is the highest value among its features covering it; the
      value of an intersection region is the lowest of those values across the tracks. Touching regions with the
      same value are joined. Features keep only ``start``, ``end`` and ``value``.

    Args:
        tracks: List of GenomicTrack objects
        operation: Merge operation ("union", "intersection", "mean", "max")

    Returns:
        Merged GenomicTrack

    Raises:
        ValueError: If the operation is not one of the four above.
    """
    if not tracks:
        return GenomicTrack()

    if operation not in ("union", "intersection", "mean", "max"):
        raise ValueError(f"Unsupported merge operation {operation!r}; use 'union', 'intersection', 'mean' or 'max'")

    logger.info(f"Merging {len(tracks)} tracks using {operation} operation")

    merged_track = GenomicTrack(name=f"Merged_{operation}", description=f"Merged track using {operation} operation")

    # Get all chromosomes, in a fixed order (set order varies from run to run)
    all_chromosomes = set()
    for track in tracks:
        all_chromosomes.update(track.get_chromosomes())

    for chromosome in sorted(all_chromosomes):
        # Collect all features for this chromosome, as copies: merging must not alter the input tracks
        all_features: List[Dict[str, Any]] = []
        for track in tracks:
            all_features.extend(dict(feature) for feature in track.get_features(chromosome))

        if not all_features:
            continue

        # Sort features by position
        all_features.sort(key=lambda x: x["start"])

        if operation == "union":
            # Merge overlapping features
            merged_track.data[chromosome] = all_features
            merged_track.merge_overlapping_features(chromosome)

        elif operation in ["mean", "max"]:
            # Group overlapping features and aggregate values
            current_group = [all_features[0]]
            group_end = all_features[0]["end"]
            merged_features = []

            for feature in all_features[1:]:
                # A feature joins the group when it starts before the group, not just its last feature, has ended
                if feature["start"] <= group_end:
                    current_group.append(feature)
                    group_end = max(group_end, feature["end"])
                else:
                    # Process current group
                    merged_feature = merge_feature_group(current_group, operation)
                    merged_features.append(merged_feature)
                    current_group = [feature]
                    group_end = feature["end"]

            # Process last group
            merged_features.append(merge_feature_group(current_group, operation))

            merged_track.data[chromosome] = merged_features

        elif operation == "intersection":
            # Only keep regions covered by every track
            profiles = [_coverage_profile(track.get_features(chromosome)) for track in tracks]
            shared = profiles[0]
            for profile in profiles[1:]:
                shared = _intersect_profiles(shared, profile)
            if shared:
                merged_track.data[chromosome] = [
                    {"start": start, "end": end, "value": value} for start, end, value in shared
                ]

    logger.info(f"Merged tracks into {merged_track.get_total_features()} features")
    return merged_track


def merge_feature_group(features: List[Dict[str, Any]], operation: str) -> Dict[str, Any]:
    """Merge a group of overlapping features.

    Args:
        features: List of overlapping features
        operation: Merge operation ("mean", "max")

    Returns:
        Merged feature dictionary
    """
    start = min(f["start"] for f in features)
    end = max(f["end"] for f in features)

    if operation == "mean":
        value = statistics.mean(f["value"] for f in features)
    elif operation == "max":
        value = max(f["value"] for f in features)
    else:
        value = features[0]["value"]

    return {
        "start": start,
        "end": end,
        "value": value,
    }


def calculate_track_statistics(track: GenomicTrack) -> Dict[str, Any]:
    """Calculate comprehensive statistics for a genomic track.

    Args:
        track: GenomicTrack to analyze

    Returns:
        Dictionary with track statistics
    """
    if track.get_total_features() == 0:
        return {}

    all_values = []
    all_lengths = []
    chromosome_counts = {}

    for chromosome, features in track.data.items():
        chromosome_counts[chromosome] = len(features)

        for feature in features:
            all_values.append(feature["value"])
            all_lengths.append(feature["end"] - feature["start"])

    stats = {
        "total_features": len(all_values),
        "total_chromosomes": len(track.data),
        "chromosome_distribution": chromosome_counts,
        "mean_value": statistics.mean(all_values),
        "median_value": statistics.median(all_values),
        "min_value": min(all_values),
        "max_value": max(all_values),
        "value_std": statistics.stdev(all_values) if len(all_values) > 1 else 0,
        "mean_length": statistics.mean(all_lengths),
        "median_length": statistics.median(all_lengths),
        "min_length": min(all_lengths),
        "max_length": max(all_lengths),
    }

    # Value distribution
    value_bins: Dict[str, int] = defaultdict(int)
    for value in all_values:
        if value < 1:
            value_bins["<1"] += 1
        elif value < 10:
            value_bins["1-10"] += 1
        elif value < 100:
            value_bins["10-100"] += 1
        else:
            value_bins["100+"] += 1

    stats["value_distribution"] = dict(value_bins)

    return stats


def extract_track_region(track: GenomicTrack, chromosome: str, start: int, end: int) -> GenomicTrack:
    """Extract a specific genomic region from a track.

    Args:
        track: Source GenomicTrack
        chromosome: Chromosome name
        start: Region start position
        end: Region end position

    Returns:
        New GenomicTrack containing only the specified region
    """
    logger.info(f"Extracting region {chromosome}:{start}-{end} from track")

    region_track = GenomicTrack(
        name=f"{track.name}_region", description=f"Region {chromosome}:{start}-{end} from {track.name}"
    )

    features = track.get_features(chromosome, start, end)

    for feature in features:
        # Clip feature to region boundaries
        feature_start = max(feature["start"], start)
        feature_end = min(feature["end"], end)

        if feature_start < feature_end:
            region_track.add_feature(
                chromosome=chromosome,
                start=feature_start,
                end=feature_end,
                value=feature["value"],
                **{k: v for k, v in feature.items() if k not in ["start", "end", "value"]},
            )

    logger.info(f"Extracted {region_track.get_total_features()} features from region")
    return region_track


def _shared_pieces(
    first: List[Tuple[int, int, float]], second: List[Tuple[int, int, float]]
) -> List[Tuple[int, float, float]]:
    """``(length, value in first, value in second)`` for every stretch that both coverage profiles cover."""
    pieces: List[Tuple[int, float, float]] = []
    i = j = 0
    while i < len(first) and j < len(second):
        start = max(first[i][0], second[j][0])
        end = min(first[i][1], second[j][1])
        if start < end:
            pieces.append((end - start, first[i][2], second[j][2]))
        if first[i][1] < second[j][1]:
            i += 1
        else:
            j += 1
    return pieces


def _weighted_correlation(pieces: List[Tuple[int, float, float]]) -> Optional[float]:
    """Pearson correlation of the two values over the pieces, each piece weighted by its length.

    ``None`` when there is nothing to correlate or one of the two tracks takes a single value over the pieces (the
    correlation is undefined there; comparing the variance with zero would turn rounding noise into a correlation).
    """
    total = sum(length for length, _, _ in pieces)
    if total <= 0 or len({x for _, x, _ in pieces}) < 2 or len({y for _, _, y in pieces}) < 2:
        return None

    mean_x = sum(length * x for length, x, _ in pieces) / total
    mean_y = sum(length * y for length, _, y in pieces) / total
    covariance = sum(length * (x - mean_x) * (y - mean_y) for length, x, y in pieces)
    variance_x = sum(length * (x - mean_x) ** 2 for length, x, _ in pieces)
    variance_y = sum(length * (y - mean_y) ** 2 for length, _, y in pieces)
    return max(-1.0, min(1.0, covariance / math.sqrt(variance_x * variance_y)))


def compare_tracks(track1: GenomicTrack, track2: GenomicTrack) -> Dict[str, Any]:
    """Compare two genomic tracks.

    The tracks are compared base by base on the chromosomes they share: a feature's value applies to every base it
    covers, and where several features of one track overlap the highest value counts (the convention of
    ``merge_tracks``). ``correlation_coefficient`` is the Pearson correlation of the two tracks' values over the
    bases that both cover, so it does not depend on how the features are cut or in which order they were added. (It
    used to pair the first features of the two tracks by list position, whatever their coordinates, which made the
    number a function of insertion order.)

    Args:
        track1: First GenomicTrack
        track2: Second GenomicTrack

    Returns:
        Dictionary with comparison statistics: the feature counts, ``shared_chromosomes``, ``unique_to_track1``,
        ``unique_to_track2``, ``shared_bases`` (bases covered by both tracks) and ``correlation_coefficient`` (``None``
        when there is no shared base, or one of the tracks has a single value over the shared bases).
    """
    logger.info("Comparing genomic tracks")

    comparison: Dict[str, Any] = {
        "track1_features": track1.get_total_features(),
        "track2_features": track2.get_total_features(),
        "shared_chromosomes": 0,
        "shared_bases": 0,
        "correlation_coefficient": None,
    }

    # Compare chromosome coverage
    chr1 = set(track1.get_chromosomes())
    chr2 = set(track2.get_chromosomes())
    comparison["shared_chromosomes"] = len(chr1 & chr2)
    comparison["unique_to_track1"] = len(chr1 - chr2)
    comparison["unique_to_track2"] = len(chr2 - chr1)

    # Compare the values over the bases both tracks cover
    pieces: List[Tuple[int, float, float]] = []
    for chromosome in sorted(chr1 & chr2):
        pieces.extend(
            _shared_pieces(
                _coverage_profile(track1.get_features(chromosome)),
                _coverage_profile(track2.get_features(chromosome)),
            )
        )

    comparison["shared_bases"] = sum(length for length, _, _ in pieces)
    comparison["correlation_coefficient"] = _weighted_correlation(pieces)

    return comparison


def generate_track_report(track: GenomicTrack, output_path: Optional[str | Path] = None) -> str:
    """Generate a comprehensive track analysis report.

    Args:
        track: GenomicTrack to analyze
        output_path: Optional path to save the report

    Returns:
        Formatted track report
    """
    report_lines = []
    report_lines.append("=" * 60)
    report_lines.append("GENOMIC TRACK ANALYSIS REPORT")
    report_lines.append("=" * 60)
    report_lines.append("")

    if track.name:
        report_lines.append(f"Track Name: {track.name}")
    if track.description:
        report_lines.append(f"Description: {track.description}")
    report_lines.append("")

    # Basic statistics
    stats = calculate_track_statistics(track)

    if stats:
        report_lines.append("Track Statistics:")
        report_lines.append(f"  Total Features: {stats.get('total_features', 0):,}")
        report_lines.append(f"  Chromosomes: {stats.get('total_chromosomes', 0)}")
        report_lines.append(f"  Mean Value: {stats.get('mean_value', 0):.2f}")
        report_lines.append(f"  Median Value: {stats.get('median_value', 0):.2f}")
        report_lines.append(f"  Mean Length: {stats.get('mean_length', 0):.0f} bp")
        report_lines.append("")

        # Chromosome distribution (top 10)
        chr_dist = stats.get("chromosome_distribution", {})
        if chr_dist:
            report_lines.append("Chromosome Distribution (Top 10):")
            sorted_chrs = sorted(chr_dist.items(), key=lambda x: x[1], reverse=True)[:10]
            for chr_name, count in sorted_chrs:
                report_lines.append(f"  {chr_name}: {count:,} features")
            report_lines.append("")

    report = "\n".join(report_lines)

    if output_path:
        output_path = Path(output_path)
        with open(output_path, "w") as f:
            f.write(report)
        logger.info(f"Track report saved to {output_path}")

    return report
