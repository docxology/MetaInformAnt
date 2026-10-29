"""PDB file parsing and manipulation utilities.

This module provides tools for reading, writing, and analyzing PDB files
containing protein structure data.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from metainformant.core.utils import logging
from metainformant.protein._network import get_protein_api_timeout

logger = logging.get_logger(__name__)

PDB_DOWNLOAD_BASE_URL = "https://files.rcsb.org/download"


def _residue_key(atom: Dict[str, Any]) -> Tuple[Any, Any, str]:
    """(chain, residue number, insertion code) of an atom dict, so that residues 52 and 52A stay apart."""
    return (atom["chain_id"], atom["res_seq"], str(atom.get("ins_code") or atom.get("icode") or "").strip())


def fetch_pdb_structure(pdb_id: str, out_dir: Path, *, fmt: str = "pdb") -> Path:
    """Download PDB structure from RCSB PDB database.

    Args:
        pdb_id: PDB identifier (e.g., "1ABC")
        out_dir: Output directory
        fmt: File format ("pdb" or "cif")

    Returns:
        Path to downloaded file

    Raises:
        requests.RequestException: If download fails

    Example:
        >>> # This would download if PDB ID exists
        >>> # path = fetch_pdb_structure("1ABC", Path("structures/"))
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    if fmt == "pdb":
        url = f"{PDB_DOWNLOAD_BASE_URL}/{pdb_id}.pdb"
        filename = f"{pdb_id}.pdb"
    elif fmt == "cif":
        url = f"{PDB_DOWNLOAD_BASE_URL}/{pdb_id}.cif"
        filename = f"{pdb_id}.cif"
    else:
        raise ValueError(f"Unsupported format: {fmt}")

    output_path = out_dir / filename

    import requests

    logger.info(f"Downloading PDB structure {pdb_id}")

    try:
        # The shared helper validates a positive, bounded campaign timeout.
        response = requests.get(  # nosec B113
            url,
            timeout=get_protein_api_timeout(),
        )
        response.raise_for_status()

        with open(output_path, "wb") as f:
            f.write(response.content)

        logger.info(f"Downloaded PDB structure to {output_path}")
        return output_path

    except requests.RequestException as e:
        logger.error(f"Failed to download PDB structure {pdb_id}: {e}")
        raise


def parse_pdb_atoms(pdb_content: str) -> List[Dict[str, Any]]:
    """Parse ATOM records from PDB file content.

    Args:
        pdb_content: PDB file content as string

    Returns:
        List of atom dictionaries

    Example:
        >>> pdb_text = "ATOM      1  CA  ALA A   1      10.000  20.000  30.000  1.00 20.00           C"
        >>> atoms = parse_pdb_atoms(pdb_text)
        >>> len(atoms) >= 1
        True
    """
    atoms = []

    for line in pdb_content.split("\n"):
        line = line.strip()
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                atom = {
                    "record_type": line[0:6].strip(),
                    "serial": int(line[6:11].strip()),
                    "name": line[12:16].strip(),
                    "alt_loc": line[16:17].strip(),
                    "res_name": line[17:20].strip(),
                    "chain_id": line[21:22].strip(),
                    "res_seq": int(line[22:26].strip()),
                    "ins_code": line[26:27].strip(),
                    "x": float(line[30:38].strip()),
                    "y": float(line[38:46].strip()),
                    "z": float(line[46:54].strip()),
                    "occupancy": float(line[54:60].strip()) if line[54:60].strip() else 1.0,
                    "temp_factor": float(line[60:66].strip()) if line[60:66].strip() else 0.0,
                    "element": line[76:78].strip(),
                    "charge": line[78:80].strip(),
                }
                atoms.append(atom)
            except (ValueError, IndexError) as e:
                logger.warning(f"Skipping malformed PDB line: {e}")
                continue

    return atoms


def load_pdb_file(path: Path) -> Dict[str, Any]:
    """Load and parse a PDB file.

    Args:
        path: Path to PDB file

    Returns:
        Dictionary containing parsed PDB data

    Raises:
        FileNotFoundError: If file doesn't exist

    Example:
        >>> # Assuming PDB file exists
        >>> # pdb_data = load_pdb_file(Path("structure.pdb"))
        >>> # "atoms" in pdb_data
        >>> # True
    """
    if not path.exists():
        raise FileNotFoundError(f"PDB file not found: {path}")

    with open(path, "r") as f:
        content = f.read()

    atoms = parse_pdb_atoms(content)

    # Extract additional information
    pdb_id = path.stem.upper()
    chains = set(atom["chain_id"] for atom in atoms if atom["chain_id"])
    residues = set(_residue_key(atom) for atom in atoms)

    pdb_data = {
        "pdb_id": pdb_id,
        "atoms": atoms,
        "n_atoms": len(atoms),
        "chains": sorted(chains),
        "n_residues": len(residues),
        "coordinates": np.array([[atom["x"], atom["y"], atom["z"]] for atom in atoms]),
    }

    logger.info(f"Loaded PDB file {path} with {len(atoms)} atoms")
    return pdb_data


def save_pdb_file(structure: Dict[str, Any], path: Path) -> None:
    """Save structure data to PDB format.

    Args:
        structure: Structure dictionary from load_pdb_file
        path: Output path

    Example:
        >>> # Assuming structure data exists
        >>> # save_pdb_file(structure, Path("output.pdb"))
    """
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w") as f:
        # Write HEADER
        f.write(f"HEADER    PROTEIN STRUCTURE                          {structure.get('pdb_id', 'XXXX')}\n")

        # Write ATOM records
        for atom in structure.get("atoms", []):
            line = (
                f"{atom['record_type']:<6}"
                f"{atom['serial']:>5d} "
                f"{atom['name']:<4}"
                f"{atom['alt_loc']:<1}"
                f"{atom['res_name']:<3} "
                f"{atom['chain_id']:<1}"
                f"{atom['res_seq']:>4d}"
                f"{atom['ins_code']:<1}   "
                f"{atom['x']:>8.3f}"
                f"{atom['y']:>8.3f}"
                f"{atom['z']:>8.3f}"
                f"{atom['occupancy']:>6.2f}"
                f"{atom['temp_factor']:>6.2f}          "
                f"{atom['element']:<2}"
                f"{atom['charge']:<2}\n"
            )
            f.write(line)

        # Write END
        f.write("END\n")

    logger.info(f"Saved PDB structure to {path}")


def extract_backbone_atoms(atoms: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Extract backbone atoms (N, CA, C, O) from atom list.

    Args:
        atoms: List of atom dictionaries

    Returns:
        List of backbone atoms

    Example:
        >>> atoms = [{"name": "CA", "res_seq": 1}, {"name": "CB", "res_seq": 1}]
        >>> backbone = extract_backbone_atoms(atoms)
        >>> len(backbone) == 1
        True
    """
    backbone_names = {"N", "CA", "C", "O"}
    backbone_atoms = [atom for atom in atoms if atom["name"] in backbone_names]

    return backbone_atoms


def extract_sidechain_atoms(atoms: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Extract side chain atoms from atom list.

    Args:
        atoms: List of atom dictionaries

    Returns:
        List of side chain atoms

    Example:
        >>> atoms = [{"name": "CA", "res_seq": 1}, {"name": "CB", "res_seq": 1}]
        >>> sidechain = extract_sidechain_atoms(atoms)
        >>> len(sidechain) == 1
        True
    """
    backbone_names = {"N", "CA", "C", "O"}
    sidechain_atoms = [atom for atom in atoms if atom["name"] not in backbone_names]

    return sidechain_atoms


def get_residue_atoms(atoms: List[Dict[str, Any]], chain_id: str, res_seq: int) -> List[Dict[str, Any]]:
    """Get all atoms for a specific residue.

    Args:
        atoms: List of atom dictionaries
        chain_id: Chain identifier
        res_seq: Residue sequence number

    Returns:
        List of atoms for the residue

    Example:
        >>> atoms = [{"chain_id": "A", "res_seq": 1, "name": "CA"}]
        >>> residue_atoms = get_residue_atoms(atoms, "A", 1)
        >>> len(residue_atoms) == 1
        True
    """
    residue_atoms = [atom for atom in atoms if atom["chain_id"] == chain_id and atom["res_seq"] == res_seq]

    return residue_atoms


def calculate_pdb_statistics(pdb_data: Dict[str, Any]) -> Dict[str, Any]:
    """Calculate statistics for PDB structure.

    ``residue_types`` counts residues, not atoms: every distinct residue (chain, number and insertion code) is counted
    once under its residue name, whatever its number of atoms, so water and ligand molecules are counted as well.

    Args:
        pdb_data: PDB data from load_pdb_file

    Returns:
        Statistics dictionary

    Example:
        >>> atoms = [
        ...     {"name": "N", "res_name": "GLY", "chain_id": "A", "res_seq": 1},
        ...     {"name": "CA", "res_name": "GLY", "chain_id": "A", "res_seq": 1},
        ...     {"name": "N", "res_name": "GLY", "chain_id": "A", "res_seq": 2},
        ... ]
        >>> calculate_pdb_statistics({"atoms": atoms})["residue_types"]
        {'GLY': 2}
    """
    atoms = pdb_data.get("atoms", [])

    if not atoms:
        return {"n_atoms": 0, "n_residues": 0, "chains": []}

    stats = {
        "n_atoms": len(atoms),
        "n_residues": pdb_data["n_residues"] if "n_residues" in pdb_data else len({_residue_key(a) for a in atoms}),
        "chains": pdb_data.get("chains", []),
        "n_chains": len(pdb_data.get("chains", [])),
    }

    # Atom type distribution
    atom_types: Dict[str, int] = {}
    elements: Dict[str, int] = {}

    for atom in atoms:
        # Atom name
        atom_name = atom["name"]
        atom_types[atom_name] = atom_types.get(atom_name, 0) + 1

        # Element
        element = atom.get("element", "")
        if element:
            elements[element] = elements.get(element, 0) + 1

    stats["atom_types"] = atom_types
    stats["elements"] = elements

    # Residue type distribution: one count per residue, however many atoms it has
    residue_names: Dict[Tuple[Any, Any, str], str] = {}
    for atom in atoms:
        residue_names.setdefault(_residue_key(atom), atom["res_name"])
    res_types: Dict[str, int] = {}
    for res_name in residue_names.values():
        res_types[res_name] = res_types.get(res_name, 0) + 1
    stats["residue_types"] = res_types

    return stats


def validate_pdb_file(path: Path) -> Dict[str, Any]:
    """Validate PDB file format and content.

    The whole file is read, so ``n_atoms`` counts every ATOM and HETATM record of it. A file is valid when it has a
    HEADER record and at least one atom record.

    Args:
        path: Path to PDB file

    Returns:
        Validation results

    Example:
        >>> # Assuming PDB file exists
        >>> # validation = validate_pdb_file(Path("structure.pdb"))
        >>> # validation['is_valid']
        >>> # True
    """
    validation: Dict[str, Any] = {
        "is_valid": False,
        "has_atoms": False,
        "has_header": False,
        "n_atoms": 0,
        "issues": [],
    }

    try:
        with open(path, "r") as f:
            n_atoms = 0
            has_header = False

            for line in f:
                if line.startswith("HEADER"):
                    has_header = True
                elif line.startswith("ATOM") or line.startswith("HETATM"):
                    n_atoms += 1

            validation["has_header"] = has_header
            validation["has_atoms"] = n_atoms > 0
            validation["n_atoms"] = n_atoms
            validation["is_valid"] = has_header and n_atoms > 0

            if not has_header:
                validation["issues"].append("Missing HEADER record")

            if n_atoms == 0:
                validation["issues"].append("No ATOM or HETATM records found")

    except Exception as e:
        validation["issues"].append(f"Error reading file: {e}")

    return validation


def get_pdb_sequence(pdb_data: Dict[str, Any], chain_id: Optional[str] = None) -> str:
    """Extract amino acid sequence from PDB data.

    The sequence is that of the polymer residues in the ATOM records, in file order: HETATM groups (water, ligands,
    modified residues) are not part of it. A residue name that is not one of the twenty standard amino acids is
    ``"X"``. Residues are told apart by chain, number and insertion code.

    Args:
        pdb_data: PDB data from load_pdb_file or parse_pdb_file
        chain_id: Specific chain to extract (optional)

    Returns:
        Amino acid sequence

    Example:
        >>> # Assuming PDB data exists
        >>> # sequence = get_pdb_sequence(pdb_data, "A")
        >>> # isinstance(sequence, str)
        >>> # True
    """
    atoms = pdb_data.get("atoms", [])

    if not atoms:
        return ""

    # Group the ATOM records by residue, in file order
    residues: Dict[Tuple[Any, Any, str], str] = {}
    for atom in atoms:
        if chain_id and atom["chain_id"] != chain_id:
            continue
        if str(atom.get("record_type") or atom.get("record_name") or "ATOM").strip().upper() != "ATOM":
            continue

        residues.setdefault(_residue_key(atom), atom["res_name"])

    # Convert 3-letter codes to 1-letter codes
    aa_codes = {
        "ALA": "A",
        "ARG": "R",
        "ASN": "N",
        "ASP": "D",
        "CYS": "C",
        "GLU": "E",
        "GLN": "Q",
        "GLY": "G",
        "HIS": "H",
        "ILE": "I",
        "LEU": "L",
        "LYS": "K",
        "MET": "M",
        "PHE": "F",
        "PRO": "P",
        "SER": "S",
        "THR": "T",
        "TRP": "W",
        "TYR": "Y",
        "VAL": "V",
    }

    sequence = []
    for res_name in residues.values():
        aa = aa_codes.get(res_name, "X")  # X for unknown
        sequence.append(aa)

    return "".join(sequence)


def find_pdb_contacts(atoms: List[Dict[str, Any]], distance_threshold: float = 5.0) -> List[Tuple[int, int, float]]:
    """Find atomic contacts in PDB structure.

    Args:
        atoms: List of atom dictionaries
        distance_threshold: Distance threshold for contacts (Å)

    Returns:
        List of (atom1_index, atom2_index, distance) tuples

    Example:
        >>> atoms = [
        ...     {"x": 0.0, "y": 0.0, "z": 0.0},
        ...     {"x": 4.0, "y": 0.0, "z": 0.0}
        ... ]
        >>> contacts = find_pdb_contacts(atoms, distance_threshold=5.0)
        >>> len(contacts) >= 1
        True
    """
    contacts = []

    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            atom1 = atoms[i]
            atom2 = atoms[j]

            # Calculate distance
            dx = atom1["x"] - atom2["x"]
            dy = atom1["y"] - atom2["y"]
            dz = atom1["z"] - atom2["z"]
            distance = (dx**2 + dy**2 + dz**2) ** 0.5

            if distance <= distance_threshold:
                contacts.append((i, j, distance))

    return contacts
