"""Geographic visualization for GWAS sample and allele frequency data.

This module provides plots for visualizing sample locations, allele frequency
gradients, and population counts on geographic coordinate systems using
matplotlib scatter and bubble charts.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core import logging

logger = logging.get_logger(__name__)

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    plt = None

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None


def sample_map(
    metadata: Dict[str, Dict[str, Any]],
    output_file: Optional[Union[str, Path]] = None,
    color_by: str = "population",
    size_by: Optional[str] = None,
    title: str = "Sample Geographic Distribution",
) -> Dict[str, Any]:
    """Create a geographic scatter plot of sample locations.

    Plots sample positions on a latitude/longitude coordinate system, colored
    by a categorical metadata field and optionally sized by a numeric field.

    Args:
        metadata: Dictionary of {sample_id: {latitude: float, longitude: float,
            population: str, ...}}.  Samples missing latitude or longitude are
            skipped with a warning.
        output_file: Optional path to save the figure.
        color_by: Metadata field name used for coloring points (e.g. "population").
        size_by: Optional metadata field name for scaling point sizes by a
            numeric value.
        title: Title for the plot.

    Returns:
        Dictionary with keys:
            - status: "success", "failed", or "skipped"
            - output_path: Path to saved file or None
            - n_samples_plotted: Number of samples successfully plotted
            - n_skipped: Number of samples skipped (missing coordinates)
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available for sample map")
        return {
            "status": "skipped",
            "output_path": None,
            "n_samples_plotted": 0,
            "n_skipped": 0,
        }

    if not HAS_NUMPY:
        logger.warning("numpy not available for sample map")
        return {
            "status": "skipped",
            "output_path": None,
            "n_samples_plotted": 0,
            "n_skipped": 0,
        }

    try:
        if not metadata:
            logger.error("Empty metadata provided to sample_map")
            return {
                "status": "failed",
                "output_path": None,
                "n_samples_plotted": 0,
                "n_skipped": 0,
            }

        # Collect valid samples with coordinates
        plotable: List[Dict[str, Any]] = []
        n_skipped = 0

        for sample_id, sample_info in metadata.items():
            lat = sample_info.get("latitude")
            lon = sample_info.get("longitude")
            if lat is None or lon is None:
                n_skipped += 1
                continue
            try:
                entry: Dict[str, Any] = {
                    "sample_id": sample_id,
                    "latitude": float(lat),
                    "longitude": float(lon),
                    "group": str(sample_info.get(color_by, "unknown")),
                }
                if size_by is not None and size_by in sample_info:
                    entry["size_val"] = float(sample_info[size_by])
                plotable.append(entry)
            except (ValueError, TypeError):
                n_skipped += 1

        if n_skipped > 0:
            logger.warning(f"Skipped {n_skipped} samples with missing or invalid coordinates")

        if not plotable:
            logger.error("No samples with valid coordinates to plot")
            return {
                "status": "failed",
                "output_path": None,
                "n_samples_plotted": 0,
                "n_skipped": n_skipped,
            }

        # Organize by group for coloring
        groups: Dict[str, List[Dict[str, Any]]] = {}
        for entry in plotable:
            grp = entry["group"]
            if grp not in groups:
                groups[grp] = []
            groups[grp].append(entry)

        unique_groups = sorted(groups.keys())
        cmap_colors = plt.cm.tab10(np.linspace(0, 1, max(len(unique_groups), 1)))
        group_to_color = {g: cmap_colors[i] for i, g in enumerate(unique_groups)}

        fig, ax = plt.subplots(figsize=(12, 8))

        for grp in unique_groups:
            entries = groups[grp]
            lons = [e["longitude"] for e in entries]
            lats = [e["latitude"] for e in entries]

            if size_by is not None and all("size_val" in e for e in entries):
                raw_sizes = np.array([e["size_val"] for e in entries])
                # Normalize sizes to range [30, 300]
                s_min, s_max = raw_sizes.min(), raw_sizes.max()
                if s_max > s_min:
                    sizes = 30 + 270 * (raw_sizes - s_min) / (s_max - s_min)
                else:
                    sizes = np.full_like(raw_sizes, 80.0)
            else:
                sizes = 80

            ax.scatter(
                lons,
                lats,
                c=[group_to_color[grp]],
                s=sizes,
                label=grp,
                alpha=0.7,
                edgecolors="black",
                linewidths=0.5,
            )

        ax.set_xlabel("Longitude", fontsize=12)
        ax.set_ylabel("Latitude", fontsize=12)
        ax.set_title(title, fontsize=14, pad=20)
        ax.grid(True, alpha=0.3)
        ax.legend(
            title=color_by.capitalize(),
            bbox_to_anchor=(1.05, 1),
            loc="upper left",
            fontsize=9,
        )

        # Frame border styling for geographic context
        for spine in ax.spines.values():
            spine.set_linewidth(1.5)

        plt.tight_layout()

        output_path_str: Optional[str] = None
        if output_file:
            output_file = Path(output_file)
            output_file.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(output_file, dpi=300, bbox_inches="tight")
            output_path_str = str(output_file)
            logger.info(f"Saved sample map to {output_file}")

        plt.close(fig)
        return {
            "status": "success",
            "output_path": output_path_str,
            "n_samples_plotted": len(plotable),
            "n_skipped": n_skipped,
        }

    except Exception as e:
        logger.error(f"Error creating sample map: {e}")
        return {
            "status": "failed",
            "output_path": None,
            "n_samples_plotted": 0,
            "n_skipped": 0,
        }


def allele_frequency_map(
    metadata: Dict[str, Dict[str, Any]],
    allele_freqs: Dict[str, Dict[str, float]],
    variant_id: str,
    output_file: Optional[Union[str, Path]] = None,
    title: Optional[str] = None,
) -> Dict[str, Any]:
    """Create a geographic allele frequency gradient map for a variant.

    For each population, computes the centroid latitude/longitude from sample
    metadata and plots a scatter point colored by alternate allele frequency
    (gradient from blue=0 to red=1) and sized by sample count.

    Args:
        metadata: Dictionary of {sample_id: {latitude: float, longitude: float,
            population: str, ...}}.
        allele_freqs: Dictionary mapping population_name to
            {"ref_freq": float, "alt_freq": float}.
        variant_id: Identifier for the variant being plotted.
        output_file: Optional path to save the figure.
        title: Optional custom title.  Defaults to
            "Allele Frequency Map: {variant_id}".

    Returns:
        Dictionary with keys:
            - status: "success", "failed", or "skipped"
            - output_path: Path to saved file or None
            - n_populations: Number of populations plotted
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available for allele frequency map")
        return {"status": "skipped", "output_path": None, "n_populations": 0}

    if not HAS_NUMPY:
        logger.warning("numpy not available for allele frequency map")
        return {"status": "skipped", "output_path": None, "n_populations": 0}

    try:
        if not metadata or not allele_freqs:
            logger.error("Empty metadata or allele_freqs provided")
            return {"status": "failed", "output_path": None, "n_populations": 0}

        # Group samples by population and compute centroids
        pop_samples: Dict[str, List[Tuple[float, float]]] = {}
        for sample_id, sample_info in metadata.items():
            pop = sample_info.get("population")
            lat = sample_info.get("latitude")
            lon = sample_info.get("longitude")
            if pop is None or lat is None or lon is None:
                continue
            try:
                coord = (float(lat), float(lon))
            except (ValueError, TypeError):
                continue
            if pop not in pop_samples:
                pop_samples[pop] = []
            pop_samples[pop].append(coord)

        # Build plot data: only populations present in both metadata and allele_freqs
        plot_data: List[Dict[str, Any]] = []
        for pop_name, freq_info in allele_freqs.items():
            if pop_name not in pop_samples:
                continue
            coords = pop_samples[pop_name]
            centroid_lat = np.mean([c[0] for c in coords])
            centroid_lon = np.mean([c[1] for c in coords])
            alt_freq = freq_info.get("alt_freq", 0.0)
            plot_data.append(
                {
                    "population": pop_name,
                    "lat": centroid_lat,
                    "lon": centroid_lon,
                    "alt_freq": alt_freq,
                    "n_samples": len(coords),
                }
            )

        if not plot_data:
            logger.error("No populations with both coordinates and allele frequencies")
            return {"status": "failed", "output_path": None, "n_populations": 0}

        fig, ax = plt.subplots(figsize=(12, 8))

        lons = [d["lon"] for d in plot_data]
        lats = [d["lat"] for d in plot_data]
        alt_freqs = [d["alt_freq"] for d in plot_data]
        counts = np.array([d["n_samples"] for d in plot_data], dtype=float)

        # Scale bubble sizes by sample count: range [60, 400]
        c_min, c_max = counts.min(), counts.max()
        if c_max > c_min:
            sizes = 60 + 340 * (counts - c_min) / (c_max - c_min)
        else:
            sizes = np.full_like(counts, 150.0)

        scatter = ax.scatter(
            lons,
            lats,
            c=alt_freqs,
            s=sizes,
            cmap="coolwarm",
            vmin=0.0,
            vmax=1.0,
            alpha=0.8,
            edgecolors="black",
            linewidths=0.8,
        )

        cbar = plt.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Alt Allele Frequency", fontsize=10)

        # Label each population
        for d in plot_data:
            ax.annotate(
                f"{d['population']}\n(n={d['n_samples']})",
                (d["lon"], d["lat"]),
                xytext=(8, 8),
                textcoords="offset points",
                fontsize=8,
                alpha=0.9,
            )

        plot_title = title if title else f"Allele Frequency Map: {variant_id}"
        ax.set_xlabel("Longitude", fontsize=12)
        ax.set_ylabel("Latitude", fontsize=12)
        ax.set_title(plot_title, fontsize=14, pad=20)
        ax.grid(True, alpha=0.3)

        for spine in ax.spines.values():
            spine.set_linewidth(1.5)

        plt.tight_layout()

        output_path_str: Optional[str] = None
        if output_file:
            output_file = Path(output_file)
            output_file.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(output_file, dpi=300, bbox_inches="tight")
            output_path_str = str(output_file)
            logger.info(f"Saved allele frequency map to {output_file}")

        plt.close(fig)
        return {
            "status": "success",
            "output_path": output_path_str,
            "n_populations": len(plot_data),
        }

    except Exception as e:
        logger.error(f"Error creating allele frequency map: {e}")
        return {"status": "failed", "output_path": None, "n_populations": 0}


def population_count_map(
    metadata: Dict[str, Dict[str, Any]],
    output_file: Optional[Union[str, Path]] = None,
    title: str = "Sample Counts by Location",
) -> Dict[str, Any]:
    """Create a geographic bubble chart of sample counts per population.

    Groups samples by population, computes centroid coordinates, and draws
    bubbles sized by sample count at each population's geographic centroid.

    Args:
        metadata: Dictionary of {sample_id: {latitude: float, longitude: float,
            population: str, ...}}.
        output_file: Optional path to save the figure.
        title: Title for the plot.

    Returns:
        Dictionary with keys:
            - status: "success", "failed", or "skipped"
            - output_path: Path to saved file or None
            - n_populations: Number of populations plotted
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available for population count map")
        return {"status": "skipped", "output_path": None, "n_populations": 0}

    if not HAS_NUMPY:
        logger.warning("numpy not available for population count map")
        return {"status": "skipped", "output_path": None, "n_populations": 0}

    try:
        if not metadata:
            logger.error("Empty metadata provided to population_count_map")
            return {"status": "failed", "output_path": None, "n_populations": 0}

        # Group samples by population
        pop_coords: Dict[str, List[Tuple[float, float]]] = {}
        for sample_id, sample_info in metadata.items():
            pop = sample_info.get("population")
            lat = sample_info.get("latitude")
            lon = sample_info.get("longitude")
            if pop is None or lat is None or lon is None:
                continue
            try:
                coord = (float(lat), float(lon))
            except (ValueError, TypeError):
                continue
            if pop not in pop_coords:
                pop_coords[pop] = []
            pop_coords[pop].append(coord)

        if not pop_coords:
            logger.error("No populations with valid coordinates")
            return {"status": "failed", "output_path": None, "n_populations": 0}

        # Compute centroids and counts
        pop_data: List[Dict[str, Any]] = []
        for pop_name, coords in pop_coords.items():
            centroid_lat = np.mean([c[0] for c in coords])
            centroid_lon = np.mean([c[1] for c in coords])
            pop_data.append(
                {
                    "population": pop_name,
                    "lat": centroid_lat,
                    "lon": centroid_lon,
                    "count": len(coords),
                }
            )

        fig, ax = plt.subplots(figsize=(12, 8))

        lons = [d["lon"] for d in pop_data]
        lats = [d["lat"] for d in pop_data]
        counts = np.array([d["count"] for d in pop_data], dtype=float)

        # Scale bubble sizes: range [100, 600]
        c_min, c_max = counts.min(), counts.max()
        if c_max > c_min:
            sizes = 100 + 500 * (counts - c_min) / (c_max - c_min)
        else:
            sizes = np.full_like(counts, 250.0)

        unique_pops = [d["population"] for d in pop_data]
        cmap_colors = plt.cm.tab10(np.linspace(0, 1, max(len(unique_pops), 1)))

        ax.scatter(
            lons,
            lats,
            s=sizes,
            c=cmap_colors[: len(unique_pops)],
            alpha=0.6,
            edgecolors="black",
            linewidths=1.0,
        )

        # Label each bubble with population name and count
        for i, d in enumerate(pop_data):
            ax.annotate(
                f"{d['population']}\n(n={d['count']})",
                (d["lon"], d["lat"]),
                ha="center",
                va="center",
                fontsize=9,
                fontweight="bold",
            )

        ax.set_xlabel("Longitude", fontsize=12)
        ax.set_ylabel("Latitude", fontsize=12)
        ax.set_title(title, fontsize=14, pad=20)
        ax.grid(True, alpha=0.3)

        for spine in ax.spines.values():
            spine.set_linewidth(1.5)

        plt.tight_layout()

        output_path_str: Optional[str] = None
        if output_file:
            output_file = Path(output_file)
            output_file.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(output_file, dpi=300, bbox_inches="tight")
            output_path_str = str(output_file)
            logger.info(f"Saved population count map to {output_file}")

        plt.close(fig)
        return {
            "status": "success",
            "output_path": output_path_str,
            "n_populations": len(pop_data),
        }

    except Exception as e:
        logger.error(f"Error creating population count map: {e}")
        return {"status": "failed", "output_path": None, "n_populations": 0}
