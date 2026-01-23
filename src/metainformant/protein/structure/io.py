"""Protein structure file I/O utilities.

This module provides functions for reading and writing protein structure files
in various formats (PDB, MMCIF, etc.).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from metainformant.core import io, logging

logger = logging.get_logger(__name__)


def parse_pdb_file(file_path: str | Path) -> Dict[str, Any]:
    """Parse a PDB file into structured data.

    Args:
        file_path: Path to PDB file

    Returns:
        Dictionary containing parsed structure data

    Raises:
        FileNotFoundError: If PDB file doesn't exist
        ValueError: If PDB file format is invalid

    Example:
        >>> # pdb_data = parse_pdb_file("protein.pdb")
        >>> # isinstance(pdb_data, dict)
        >>> # True
    """
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"PDB file not found: {file_path}")

    atoms = []
    coordinates = []
    residues = {}

    with open(file_path, "r") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue

            # Parse ATOM/HETATM records
            if line.startswith(("ATOM", "HETATM")):
                try:
                    atom_data = _parse_atom_line(line)
                    atoms.append(atom_data)

                    # Extract coordinates
                    coord = [atom_data["x"], atom_data["y"], atom_data["z"]]
                    coordinates.append(coord)

                    # Group by residue
                    res_key = (atom_data["chain_id"], atom_data["res_seq"])
                    if res_key not in residues:
                        residues[res_key] = []
                    residues[res_key].append(atom_data)

                except (ValueError, IndexError) as e:
                    logger.warning(f"Error parsing line {line_num}: {e}")
                    continue

    if not atoms:
        raise ValueError(f"No valid atoms found in PDB file: {file_path}")

    return {
        "atoms": atoms,
        "coordinates": np.array(coordinates),
        "residues": residues,
        "file_path": str(file_path),
        "n_atoms": len(atoms),
        "n_residues": len(residues),
    }


def _parse_atom_line(line: str) -> Dict[str, Any]:
    """Parse a single ATOM/HETATM line from PDB format.

    Args:
        line: PDB atom line

    Returns:
        Parsed atom data dictionary
    """
    # PDB format specification
    # Columns: 1-6: Record name, 7-11: Serial, 13-16: Name, 17: AltLoc,
    # 18-20: ResName, 22: ChainID, 23-26: ResSeq, 27: ICode,
    # 31-38: X, 39-46: Y, 47-54: Z, 55-60: Occupancy, 61-66: TempFactor,
    # 77-78: Element, 79-80: Charge

    return {
        "record_name": line[0:6].strip(),
        "serial": int(line[6:11].strip() or "0"),
        "name": line[12:16].strip(),
        "alt_loc": line[16:17].strip(),
        "res_name": line[17:20].strip(),
        "chain_id": line[21:22].strip() or "A",
        "res_seq": int(line[22:26].strip() or "0"),
        "icode": line[26:27].strip(),
        "x": float(line[30:38].strip() or "0.0"),
        "y": float(line[38:46].strip() or "0.0"),
        "z": float(line[46:54].strip() or "0.0"),
        "occupancy": float(line[54:60].strip() or "1.0"),
        "temp_factor": float(line[60:66].strip() or "0.0"),
        "element": line[76:78].strip(),
        "charge": line[78:80].strip(),
    }


def write_pdb_file(structure_data: Dict[str, Any], file_path: str | Path) -> None:
    """Write structure data to a PDB file.

    Args:
        structure_data: Structure data dictionary
        file_path: Output file path

    Example:
        >>> # write_pdb_file(structure_data, "output.pdb")
    """
    file_path = Path(file_path)
    atoms = structure_data.get("atoms", [])

    with open(file_path, "w") as f:
        # Write HEADER
        f.write("HEADER    PROTEIN STRUCTURE                          \n")

        # Write ATOM records
        for atom in atoms:
            line = _format_atom_line(atom)
            f.write(line + "\n")

        # Write END
        f.write("END\n")

    logger.info(f"Written {len(atoms)} atoms to {file_path}")


def _format_atom_line(atom: Dict[str, Any]) -> str:
    """Format atom data into PDB ATOM line.

    Args:
        atom: Atom data dictionary

    Returns:
        Formatted PDB line
    """
    return (
        "{record_name:<6}{serial:>5} {name:<4}{alt_loc:1}{res_name:>3} {chain_id:1}"
        "{res_seq:>4}{icode:1}   {x:>8.3f}{y:>8.3f}{z:>8.3f}{occupancy:>6.2f}"
        "{temp_factor:>6.2f}          {element:>2}{charge:>2}"
    ).format(
        record_name=atom.get("record_name", "ATOM"),
        serial=atom.get("serial", 0),
        name=atom.get("name", "CA"),
        alt_loc=atom.get("alt_loc", ""),
        res_name=atom.get("res_name", "ALA"),
        chain_id=atom.get("chain_id", "A"),
        res_seq=atom.get("res_seq", 1),
        icode=atom.get("icode", ""),
        x=atom.get("x", 0.0),
        y=atom.get("y", 0.0),
        z=atom.get("z", 0.0),
        occupancy=atom.get("occupancy", 1.0),
        temp_factor=atom.get("temp_factor", 0.0),
        element=atom.get("element", "C"),
        charge=atom.get("charge", ""),
    )


def parse_mmcif_file(file_path: str | Path) -> Dict[str, Any]:
    """Parse an mmCIF file into structured data.

    Args:
        file_path: Path to mmCIF file

    Returns:
        Dictionary containing parsed structure data

    Raises:
        FileNotFoundError: If mmCIF file doesn't exist

    Example:
        >>> # cif_data = parse_mmcif_file("protein.cif")
        >>> # isinstance(cif_data, dict)
        >>> # True
    """
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"mmCIF file not found: {file_path}")

    # Simplified mmCIF parsing
    # Real implementation would use a proper mmCIF parser
    data = {}

    with open(file_path, "r") as f:
        content = f.read()

    # Extract basic information
    atoms = []
    coordinates = []

    # Look for atom_site loop
    lines = content.split("\n")
    in_atom_site = False

    for line in lines:
        line = line.strip()
        if line.startswith("_atom_site."):
            in_atom_site = True
            continue
        elif line.startswith("_") or line.startswith("#"):
            continue
        elif in_atom_site and line:
            # Parse atom site data
            # This is highly simplified - real mmCIF has complex format
            parts = line.split()
            if len(parts) >= 11:  # Minimum fields for coordinates
                try:
                    atom_data = {
                        "group_PDB": parts[0],
                        "id": int(parts[1]),
                        "type_symbol": parts[2],
                        "label_atom_id": parts[3],
                        "label_comp_id": parts[5],
                        "label_asym_id": parts[6],
                        "label_seq_id": int(parts[8]) if parts[8] != "." else 0,
                        "Cartn_x": float(parts[10]),
                        "Cartn_y": float(parts[11]),
                        "Cartn_z": float(parts[12]),
                    }
                    atoms.append(atom_data)
                    coordinates.append([atom_data["Cartn_x"], atom_data["Cartn_y"], atom_data["Cartn_z"]])
                except (ValueError, IndexError):
                    continue

    if not atoms:
        logger.warning(f"No atoms found in mmCIF file: {file_path}")
        return {"atoms": [], "coordinates": np.array([]), "file_path": str(file_path)}

    return {"atoms": atoms, "coordinates": np.array(coordinates), "file_path": str(file_path), "n_atoms": len(atoms)}


def convert_pdb_to_mmcif(pdb_file: str | Path, cif_file: str | Path) -> None:
    """Convert PDB file to mmCIF format.

    Args:
        pdb_file: Input PDB file path
        cif_file: Output mmCIF file path

    Example:
        >>> # convert_pdb_to_mmcif("input.pdb", "output.cif")
    """
    # Load PDB data
    pdb_data = parse_pdb_file(pdb_file)

    # Convert to mmCIF format (simplified)
    cif_file = Path(cif_file)

    with open(cif_file, "w") as f:
        f.write("data_protein\n")
        f.write("#\n")
        f.write("loop_\n")
        f.write("_atom_site.group_PDB\n")
        f.write("_atom_site.id\n")
        f.write("_atom_site.type_symbol\n")
        f.write("_atom_site.label_atom_id\n")
        f.write("_atom_site.label_comp_id\n")
        f.write("_atom_site.label_asym_id\n")
        f.write("_atom_site.label_seq_id\n")
        f.write("_atom_site.Cartn_x\n")
        f.write("_atom_site.Cartn_y\n")
        f.write("_atom_site.Cartn_z\n")

        for i, atom in enumerate(pdb_data["atoms"], 1):
            f.write("ATOM\n")
            f.write(f"{i}\n")
            f.write(f"{atom.get('element', 'C')}\n")
            f.write(f"{atom.get('name', 'CA')}\n")
            f.write(f"{atom.get('res_name', 'ALA')}\n")
            f.write(f"{atom.get('chain_id', 'A')}\n")
            f.write(f"{atom.get('res_seq', 1)}\n")
            f.write(f"{atom.get('x', 0.0):.3f}\n")
            f.write(f"{atom.get('y', 0.0):.3f}\n")
            f.write(f"{atom.get('z', 0.0):.3f}\n")

    logger.info(f"Converted {pdb_file} to {cif_file}")


def extract_chains_from_pdb(pdb_file: str | Path, chain_ids: List[str]) -> Dict[str, Dict[str, Any]]:
    """Extract specific chains from a PDB file.

    Args:
        pdb_file: Input PDB file path
        chain_ids: List of chain IDs to extract

    Returns:
        Dictionary mapping chain IDs to structure data

    Example:
        >>> # chains = extract_chains_from_pdb("protein.pdb", ["A", "B"])
        >>> # isinstance(chains, dict)
        >>> # True
    """
    pdb_data = parse_pdb_file(pdb_file)
    atoms = pdb_data["atoms"]

    chains_data = {}

    for chain_id in chain_ids:
        chain_atoms = [atom for atom in atoms if atom["chain_id"] == chain_id]

        if chain_atoms:
            coordinates = np.array([[atom["x"], atom["y"], atom["z"]] for atom in chain_atoms])

            chains_data[chain_id] = {
                "atoms": chain_atoms,
                "coordinates": coordinates,
                "n_atoms": len(chain_atoms),
                "chain_id": chain_id,
            }

    return chains_data


def merge_pdb_files(pdb_files: List[str | Path], output_file: str | Path) -> None:
    """Merge multiple PDB files into one.

    Args:
        pdb_files: List of input PDB file paths
        output_file: Output merged PDB file path

    Example:
        >>> # merge_pdb_files(["chain1.pdb", "chain2.pdb"], "merged.pdb")
    """
    all_atoms = []
    serial_offset = 0

    for pdb_file in pdb_files:
        pdb_data = parse_pdb_file(pdb_file)

        # Adjust serial numbers
        for atom in pdb_data["atoms"]:
            atom_copy = atom.copy()
            atom_copy["serial"] += serial_offset
            all_atoms.append(atom_copy)

        serial_offset += len(pdb_data["atoms"])

    # Create merged structure
    merged_data = {
        "atoms": all_atoms,
        "coordinates": np.array([[atom["x"], atom["y"], atom["z"]] for atom in all_atoms]),
    }

    write_pdb_file(merged_data, output_file)
    logger.info(f"Merged {len(pdb_files)} PDB files into {output_file}")


def validate_pdb_file(file_path: str | Path) -> Tuple[bool, List[str]]:
    """Validate PDB file format and content.

    Args:
        file_path: Path to PDB file

    Returns:
        Tuple of (is_valid, list_of_issues)

    Example:
        >>> # valid, issues = validate_pdb_file("protein.pdb")
        >>> # isinstance(valid, bool)
        >>> # True
    """
    issues = []

    try:
        pdb_data = parse_pdb_file(file_path)
        atoms = pdb_data["atoms"]

        if not atoms:
            issues.append("No atoms found in file")
            return False, issues

        # Check for required fields
        required_fields = ["serial", "name", "res_name", "chain_id", "res_seq", "x", "y", "z"]

        for i, atom in enumerate(atoms[:10]):  # Check first 10 atoms
            for field in required_fields:
                if field not in atom:
                    issues.append(f"Atom {i+1} missing required field: {field}")

        # Check coordinate ranges (should be reasonable)
        coords = pdb_data["coordinates"]
        if np.any(np.abs(coords) > 1000):  # Unreasonably large coordinates
            issues.append("Coordinates appear to be in wrong units or corrupted")

        # Check for duplicate serial numbers
        serials = [atom["serial"] for atom in atoms]
        if len(serials) != len(set(serials)):
            issues.append("Duplicate atom serial numbers found")

        return len(issues) == 0, issues

    except Exception as e:
        issues.append(f"Parse error: {e}")
        return False, issues


def get_pdb_statistics(file_path: str | Path) -> Dict[str, Any]:
    """Get statistics about a PDB file.

    Args:
        file_path: Path to PDB file

    Returns:
        Dictionary of statistics

    Example:
        >>> # stats = get_pdb_statistics("protein.pdb")
        >>> # isinstance(stats, dict)
        >>> # True
    """
    try:
        pdb_data = parse_pdb_file(file_path)
        atoms = pdb_data["atoms"]

        if not atoms:
            return {"n_atoms": 0, "n_residues": 0, "chains": [], "file_size": Path(file_path).stat().st_size}

        # Count chains
        chains = set(atom["chain_id"] for atom in atoms)

        # Count residues
        residues = set((atom["chain_id"], atom["res_seq"]) for atom in atoms)

        # Coordinate statistics
        coords = pdb_data["coordinates"]
        coord_stats = {
            "min_coords": coords.min(axis=0).tolist(),
            "max_coords": coords.max(axis=0).tolist(),
            "center": coords.mean(axis=0).tolist(),
            "dimensions": (coords.max(axis=0) - coords.min(axis=0)).tolist(),
        }

        return {
            "n_atoms": len(atoms),
            "n_residues": len(residues),
            "chains": sorted(chains),
            "file_size": Path(file_path).stat().st_size,
            "coordinate_stats": coord_stats,
        }

    except Exception as e:
        logger.error(f"Error getting PDB statistics: {e}")
        return {"error": str(e)}
