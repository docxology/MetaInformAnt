"""AlphaFold protein structure prediction integration.

This module provides tools for accessing AlphaFold predicted structures
from the AlphaFold Protein Structure Database.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import requests

from metainformant.core.utils import logging
from metainformant.protein._network import get_protein_api_timeout

logger = logging.get_logger(__name__)

ALPHAFOLD_FILES_BASE_URL = "https://alphafold.ebi.ac.uk/files"

#: AlphaFold's pLDDT confidence bands as (name, lowest score in the band); a lower score is "very_low".
_PLDDT_BANDS: Tuple[Tuple[str, float], ...] = (("very_high", 90.0), ("confident", 70.0), ("low", 50.0))


def build_alphafold_url(uniprot_acc: str, *, version: int = 4, fmt: str = "pdb") -> str:
    """Build AlphaFold database URL for a UniProt accession.

    Args:
        uniprot_acc: UniProt accession (e.g., "P12345")
        version: AlphaFold version (default: 4)
        fmt: File format ("pdb" or "cif")

    Returns:
        AlphaFold database URL

    Example:
        >>> url = build_alphafold_url("P12345")
        >>> "alphafold.ebi.ac.uk" in url
        True
    """
    if fmt not in ("pdb", "cif"):
        raise ValueError("fmt must be 'pdb' or 'cif'")

    base_url = f"{ALPHAFOLD_FILES_BASE_URL}/AF-{uniprot_acc}-F1-model_v{version}"

    if fmt == "pdb":
        url = f"{base_url}.pdb"
    else:  # cif
        url = f"{base_url}.cif"

    return url


def fetch_alphafold_model(uniprot_acc: str, out_dir: Path, *, version: int = 4, fmt: str = "pdb") -> Path:
    """Download AlphaFold predicted structure for a UniProt accession.

    Args:
        uniprot_acc: UniProt accession
        out_dir: Output directory
        version: AlphaFold version
        fmt: File format

    Returns:
        Path to downloaded structure file

    Raises:
        requests.RequestException: If download fails

    Example:
        >>> # This would download if the accession exists
        >>> # path = fetch_alphafold_model("P12345", Path("structures/"))
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    url = build_alphafold_url(uniprot_acc, version=version, fmt=fmt)
    filename = f"AF-{uniprot_acc}-F1-model_v{version}.{fmt}"
    output_path = out_dir / filename

    logger.info(f"Downloading AlphaFold model for {uniprot_acc}")

    try:
        # The shared helper validates a positive, bounded campaign timeout.
        response = requests.get(  # nosec B113
            url,
            timeout=get_protein_api_timeout(),
        )
        response.raise_for_status()

        with open(output_path, "wb") as f:
            f.write(response.content)

        logger.info(f"Downloaded AlphaFold model to {output_path}")
        return output_path

    except requests.RequestException as e:
        logger.error(f"Failed to download AlphaFold model for {uniprot_acc}: {e}")
        raise


def _plddt_band(score: float) -> str:
    """Name of the AlphaFold confidence band of a pLDDT score (>= 90, >= 70, >= 50, else very low)."""
    for name, lowest in _PLDDT_BANDS:
        if score >= lowest:
            return name
    return "very_low"


def _atom_records(pdb_path: Path) -> List[Tuple[Tuple[str, str, str], Optional[float]]]:
    """Residue and B-factor of every ATOM record of the first model, in file order.

    The residue is the triple (chain, residue number, insertion code), so residues that share a number in
    different chains, or that differ only by insertion code, stay apart. The B-factor is None when columns 61-66
    do not hold a number. HETATM records are ignored and reading stops at the end of the first model.

    Raises:
        FileNotFoundError: If ``pdb_path`` does not exist
    """
    records: List[Tuple[Tuple[str, str, str], Optional[float]]] = []
    with open(pdb_path, "r") as handle:
        for line in handle:
            if line.startswith("ENDMDL"):
                break
            if not line.startswith("ATOM  "):
                continue
            residue = (line[21], line[22:26].strip(), line[26])
            try:
                b_factor: Optional[float] = float(line[60:66])
            except ValueError:
                b_factor = None
            records.append((residue, b_factor))
    return records


def _mean_per_residue(records: List[Tuple[Tuple[str, str, str], Optional[float]]]) -> List[float]:
    """Mean B-factor of each residue, in order of first appearance, over its atoms that have a B-factor.

    AlphaFold writes a residue's pLDDT into the B-factor of every atom of the residue, so the mean is that
    pLDDT. Counting every atom separately would instead weight a residue by its number of atoms.
    """
    by_residue: Dict[Tuple[str, str, str], List[float]] = {}
    for residue, b_factor in records:
        if b_factor is not None:
            by_residue.setdefault(residue, []).append(b_factor)
    return [sum(values) / len(values) for values in by_residue.values()]


def _read_residue_plddt(pdb_path: Path) -> List[float]:
    """One pLDDT per residue from the B-factor column of the ATOM records of the first model.

    Raises:
        FileNotFoundError: If ``pdb_path`` does not exist
        ValueError: If the file has no ATOM records, an ATOM record has no numeric B-factor, or a residue's
            B-factor is outside the 0-100 range of pLDDT scores (the file is not an AlphaFold model)
    """
    records = _atom_records(pdb_path)
    if not records:
        raise ValueError(f"{pdb_path}: no ATOM records found; is this a PDB-format AlphaFold model?")
    for (chain, number, insertion_code), b_factor in records:
        if b_factor is None:
            raise ValueError(
                f"{pdb_path}: an ATOM record of residue {chain}{number}{insertion_code} has no numeric "
                "B-factor in columns 61-66"
            )

    scores = _mean_per_residue(records)
    if not all(0.0 <= score <= 100.0 for score in scores):
        raise ValueError(
            f"{pdb_path}: B-factors outside 0-100, so they are not AlphaFold pLDDT scores "
            f"(range {min(scores):g} to {max(scores):g})"
        )
    return scores


def parse_alphafold_pae(pae_path: Path) -> Dict[str, Any]:
    """Parse a predicted aligned error (PAE) JSON file of an AlphaFold model.

    Both layouts of the AlphaFold database are read, inside a single-element list or as a bare object: the
    current one with a full ``predicted_aligned_error`` matrix, and the older sparse one listing
    ``residue1``, ``residue2`` and ``distance`` for every pair (residues numbered from 1).
    ``max_predicted_aligned_error`` is the cap AlphaFold reports for the model, not a statistic of the matrix.

    Args:
        pae_path: Path to the PAE JSON file

    Returns:
        Dictionary with ``pae_matrix`` (N x N array, row = aligned residue, column = scored residue),
        ``n_residues``, ``mean_pae`` (mean over the off-diagonal pairs; None for a single residue),
        ``max_pae`` (largest value in the matrix), ``reported_max_pae`` (the file's
        ``max_predicted_aligned_error``, or None when absent) and ``format`` (``"matrix"`` or ``"sparse"``)

    Raises:
        FileNotFoundError: If ``pae_path`` does not exist
        ValueError: If the file is not valid JSON or does not hold a complete, square matrix of finite,
            non-negative errors in one of the two layouts
    """
    with open(pae_path, "r") as handle:
        document = json.load(handle)

    if isinstance(document, list) and len(document) == 1:
        document = document[0]
    if not isinstance(document, dict):
        raise ValueError(f"{pae_path}: expected a PAE object, or a list holding one, got {type(document).__name__}")

    reported_max = document.get("max_predicted_aligned_error")
    if "predicted_aligned_error" in document:
        layout = "matrix"
        matrix = np.asarray(document["predicted_aligned_error"], dtype=float)
        if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or matrix.shape[0] == 0:
            raise ValueError(
                f"{pae_path}: predicted_aligned_error must be a non-empty square matrix, got {matrix.shape}"
            )
    elif {"residue1", "residue2", "distance"} <= document.keys():
        layout = "sparse"
        first = np.asarray(document["residue1"], dtype=int)
        second = np.asarray(document["residue2"], dtype=int)
        errors = np.asarray(document["distance"], dtype=float)
        if not (first.shape == second.shape == errors.shape) or first.ndim != 1 or first.size == 0:
            raise ValueError(f"{pae_path}: residue1, residue2 and distance must be equal-length, non-empty lists")
        n_residues = int(max(first.max(), second.max()))
        if first.min() < 1 or second.min() < 1 or first.size != n_residues**2:
            raise ValueError(
                f"{pae_path}: the sparse layout must list every pair of {n_residues} residues exactly once "
                f"(1-based), got {first.size} entries"
            )
        matrix = np.full((n_residues, n_residues), np.nan)
        matrix[first - 1, second - 1] = errors
        if np.isnan(matrix).any():
            raise ValueError(f"{pae_path}: the sparse layout lists some residue pairs more than once and misses others")
    else:
        raise ValueError(f"{pae_path}: neither predicted_aligned_error nor residue1/residue2/distance found")

    if not np.isfinite(matrix).all() or (matrix < 0).any():
        raise ValueError(f"{pae_path}: predicted aligned errors must be finite and non-negative")

    n_residues = matrix.shape[0]
    off_diagonal = matrix[~np.eye(n_residues, dtype=bool)]
    return {
        "pae_matrix": matrix,
        "n_residues": n_residues,
        "mean_pae": float(off_diagonal.mean()) if off_diagonal.size else None,
        "max_pae": float(matrix.max()),
        "reported_max_pae": None if reported_max is None else float(reported_max),
        "format": layout,
    }


def get_alphafold_metadata(
    uniprot_acc: str,
    *,
    version: int = 4,
    pdb_path: Optional[Path] = None,
    pae_path: Optional[Path] = None,
) -> Dict[str, Any]:
    """Get metadata for an AlphaFold predicted structure.

    Nothing is looked up online. The accession and version give the download URL; everything else is read
    from the local model files that are passed, and is None when the file it comes from is not given. Facts
    that only the AlphaFold service holds (organism, gene, sequence length, release date) are not returned.

    Args:
        uniprot_acc: UniProt accession
        version: AlphaFold database version the model belongs to (used for the URL)
        pdb_path: Optional local AlphaFold model in PDB format. Gives ``confidence_score``, ``n_residues``
            and ``plddt_bands``.
        pae_path: Optional local PAE JSON of the same model (see :func:`parse_alphafold_pae`). Gives
            ``mean_pae`` and ``max_pae``.

    Returns:
        Metadata dictionary with ``uniprot_accession``, ``version``, ``source``, ``method``, ``url`` and:

        - ``confidence_score``: mean pLDDT over the residues of ``pdb_path`` (None without it)
        - ``n_residues``: number of residues of ``pdb_path``, or the size of the PAE matrix when only
          ``pae_path`` is given (None without either)
        - ``plddt_bands``: residues per AlphaFold confidence band (``very_high`` >= 90, ``confident`` >= 70,
          ``low`` >= 50, ``very_low`` below), or None without ``pdb_path``
        - ``mean_pae`` / ``max_pae``: off-diagonal mean and largest predicted aligned error of ``pae_path``
          (None without it)

    Raises:
        FileNotFoundError: If a given file does not exist
        ValueError: If a given file cannot be parsed, or the PAE matrix and the model have different residue
            counts

    Example:
        >>> metadata = get_alphafold_metadata("P12345")
        >>> metadata["confidence_score"] is None  # no model file was given
        True
    """
    metadata: Dict[str, Any] = {
        "uniprot_accession": uniprot_acc,
        "version": version,
        "source": f"AlphaFold v{version}",
        "method": "Deep learning structure prediction",
        "url": build_alphafold_url(uniprot_acc, version=version),
        "confidence_score": None,
        "n_residues": None,
        "plddt_bands": None,
        "mean_pae": None,
        "max_pae": None,
    }

    if pdb_path is not None:
        scores = _read_residue_plddt(Path(pdb_path))
        bands = {name: 0 for name, _ in _PLDDT_BANDS}
        bands["very_low"] = 0
        for score in scores:
            bands[_plddt_band(score)] += 1
        metadata["confidence_score"] = math.fsum(scores) / len(scores)
        metadata["n_residues"] = len(scores)
        metadata["plddt_bands"] = bands

    if pae_path is not None:
        pae = parse_alphafold_pae(Path(pae_path))
        if metadata["n_residues"] is not None and pae["n_residues"] != metadata["n_residues"]:
            raise ValueError(
                f"PAE matrix covers {pae['n_residues']} residues but the model has {metadata['n_residues']}"
            )
        metadata["n_residues"] = pae["n_residues"]
        metadata["mean_pae"] = pae["mean_pae"]
        metadata["max_pae"] = pae["max_pae"]

    return metadata


def parse_alphafold_confidence(pdb_path: Path) -> List[float]:
    """Parse AlphaFold confidence scores (pLDDT) from a PDB file, one per residue.

    The score of a residue is the mean of the B-factors of its ATOM records (AlphaFold stores the residue's pLDDT
    in every atom). Residues are told apart by chain, residue number and insertion code, listed in file order,
    and only the first model of the file is read.

    A missing file or an ATOM record without a numeric B-factor is logged and gives an empty list.

    Args:
        pdb_path: Path to AlphaFold PDB file

    Returns:
        List of confidence scores, one per residue

    Example:
        >>> # Assuming PDB file exists
        >>> # confidence = parse_alphafold_confidence(Path("structure.pdb"))
        >>> # isinstance(confidence, list)
        >>> # True
    """
    try:
        records = _atom_records(Path(pdb_path))
    except FileNotFoundError:
        logger.warning(f"PDB file not found: {pdb_path}")
        return []

    if any(b_factor is None for _, b_factor in records):
        logger.error(f"Error parsing PDB file {pdb_path}: an ATOM record has no numeric B-factor in columns 61-66")
        return []
    return _mean_per_residue(records)


def find_alphafold_models_by_sequence(sequence: str, identity_threshold: float = 0.9) -> List[Dict[str, Any]]:
    """Find AlphaFold models for similar sequences (not implemented).

    Finding models by sequence similarity needs a sequence-search service; the AlphaFold database offers no
    similarity search that this module can call, and nothing can be searched offline. The function therefore
    raises instead of returning an empty list, which would claim that no similar model exists.

    Args:
        sequence: Protein sequence
        identity_threshold: Minimum sequence identity

    Raises:
        NotImplementedError: Always. Map the sequence to a UniProt accession with a sequence-search tool, then
            use :func:`fetch_alphafold_model`.
    """
    raise NotImplementedError(
        "find_alphafold_models_by_sequence is not implemented: it would need a sequence-similarity search "
        "service. Map the sequence to a UniProt accession with a sequence-search tool, then call "
        "fetch_alphafold_model(accession)."
    )


def batch_download_alphafold_models(
    uniprot_accessions: List[str], out_dir: Path, max_workers: int = 4
) -> Dict[str, Path | None]:
    """Download multiple AlphaFold models in parallel.

    Args:
        uniprot_accessions: List of UniProt accessions
        out_dir: Output directory
        max_workers: Maximum concurrent downloads

    Returns:
        Dictionary mapping accessions to downloaded file paths

    Example:
        >>> accessions = ["P12345", "P67890"]
        >>> models = batch_download_alphafold_models(accessions, Path("structures/"))
        >>> isinstance(models, dict)
        True
    """
    results: Dict[str, Path | None] = {}

    for acc in uniprot_accessions:
        try:
            path = fetch_alphafold_model(acc, out_dir)
            results[acc] = path
        except Exception as e:
            logger.error(f"Failed to download model for {acc}: {e}")
            results[acc] = None

    return results


def validate_alphafold_structure(pdb_path: Path) -> Dict[str, Any]:
    """Validate AlphaFold predicted structure.

    Args:
        pdb_path: Path to AlphaFold PDB file

    Returns:
        Validation results: ``is_valid`` (the file has ATOM records), ``n_atoms``, ``n_residues`` (distinct
        chain, residue number and insertion code), ``has_confidence_scores``, ``avg_confidence`` (mean pLDDT over
        the residues, not over the atoms) and ``issues``.

    Example:
        >>> # Assuming PDB file exists
        >>> # validation = validate_alphafold_structure(Path("structure.pdb"))
        >>> # validation['is_valid']
        >>> # True
    """
    validation: Dict[str, Any] = {
        "is_valid": False,
        "n_atoms": 0,
        "n_residues": 0,
        "has_confidence_scores": False,
        "avg_confidence": 0.0,
        "issues": [],
    }

    try:
        records = _atom_records(Path(pdb_path))
        confidence_scores = _mean_per_residue(records)

        validation["n_atoms"] = len(records)
        validation["n_residues"] = len({residue for residue, _ in records})
        validation["has_confidence_scores"] = len(confidence_scores) > 0

        if confidence_scores:
            validation["avg_confidence"] = sum(confidence_scores) / len(confidence_scores)

        # Basic validation
        if validation["n_atoms"] > 0:
            validation["is_valid"] = True

            # Check for reasonable confidence scores (AlphaFold range: 0-100)
            if validation["has_confidence_scores"]:
                if not (0 <= validation["avg_confidence"] <= 100):
                    validation["issues"].append("Confidence scores outside expected range")

        else:
            validation["issues"].append("No ATOM records found")

    except FileNotFoundError:
        validation["issues"].append("PDB file not found")
    except Exception as e:
        validation["issues"].append(f"Error reading PDB file: {e}")

    return validation


def get_alphafold_coverage() -> Dict[str, Any]:
    """Get statistics on AlphaFold database coverage (not implemented).

    The statistics would have to be queried from the AlphaFold Protein Structure Database. This function
    used to return fixed numbers written into the source; those were not data, so it raises instead.

    Raises:
        NotImplementedError: Always.
    """
    raise NotImplementedError(
        "get_alphafold_coverage is not implemented: database coverage statistics would have to be queried "
        "from the AlphaFold Protein Structure Database (the figures this function used to return were fixed "
        "approximations, not queried values)."
    )


def search_alphafold_by_keyword(keyword: str, max_results: int = 100) -> List[Dict[str, Any]]:
    """Search AlphaFold database by keyword (not implemented).

    A keyword search would need the AlphaFold or UniProt search services. The function raises instead of
    returning an empty list, which would claim that no structure matches.

    Args:
        keyword: Search keyword
        max_results: Maximum results to return

    Raises:
        NotImplementedError: Always. Find accessions with a UniProt search, then use
            :func:`fetch_alphafold_model`.
    """
    raise NotImplementedError(
        "search_alphafold_by_keyword is not implemented: it would need the AlphaFold or UniProt search "
        "services. Find accessions with a UniProt search, then call fetch_alphafold_model(accession)."
    )


def get_alphafold_structure_quality(pdb_path: Path) -> Dict[str, float]:
    """Assess quality metrics for AlphaFold structure.

    Every figure is per residue (see :func:`parse_alphafold_confidence`): ``plddt_score`` is the mean pLDDT over
    the residues, the ``*_confidence_residues`` entries count residues (>= 90, >= 70, below 70) and
    ``confidence_distribution`` bins residues into 0-50, 50-70, 70-90 and 90-100.

    Args:
        pdb_path: Path to AlphaFold PDB file

    Returns:
        Quality metrics

    Example:
        >>> # Assuming PDB file exists
        >>> # quality = get_alphafold_structure_quality(Path("structure.pdb"))
        >>> # isinstance(quality, dict)
        >>> # True
    """
    quality: Dict[str, Any] = {
        "plddt_score": 0.0,  # Predicted LDDT score
        "confidence_distribution": {},
        "high_confidence_residues": 0,
        "medium_confidence_residues": 0,
        "low_confidence_residues": 0,
    }

    confidence_scores = parse_alphafold_confidence(pdb_path)

    if confidence_scores:
        quality["plddt_score"] = sum(confidence_scores) / len(confidence_scores)

        # Categorize confidence
        for score in confidence_scores:
            if score >= 90:
                quality["high_confidence_residues"] += 1
            elif score >= 70:
                quality["medium_confidence_residues"] += 1
            else:
                quality["low_confidence_residues"] += 1

        # Distribution
        bins = [0, 50, 70, 90, 100]
        hist, _ = np.histogram(confidence_scores, bins=bins)
        quality["confidence_distribution"] = {
            "very_low": int(hist[0]),  # 0-50
            "low": int(hist[1]),  # 50-70
            "medium": int(hist[2]),  # 70-90
            "high": int(hist[3]),  # 90-100
        }

    return quality
