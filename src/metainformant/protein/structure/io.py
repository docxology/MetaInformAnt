"""Protein structure file I/O utilities.

This module provides functions for reading and writing protein structure files
in various formats (PDB, MMCIF, etc.).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def parse_pdb_file(file_path: str | Path) -> Dict[str, Any]:
    """Parse a PDB file into structured data.

    Args:
        file_path: Path to PDB file

    ATOM and HETATM lines whose residue number or coordinates are blank or not numbers are skipped with a warning;
    nothing is filled in for them. A blank chain identifier stays blank. ``residues`` maps (chain, residue number,
    insertion code) to the atoms of the residue.

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
    residues: Dict[Tuple[str, int, str], List[Dict[str, Any]]] = {}

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

                    # Group by residue: chain, number and insertion code (residues 52 and 52A are two residues)
                    res_key = (atom_data["chain_id"], atom_data["res_seq"], atom_data["icode"])
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

    # The residue number and the coordinates must be present: int() and float() raise ValueError on a blank field,
    # and the caller skips the line. The serial number, occupancy and B-factor are optional columns and keep the
    # conventional defaults of the format (0, 1.00 and 0.00).
    return {
        "record_name": line[0:6].strip(),
        "serial": int(line[6:11].strip() or "0"),
        "name": line[12:16].strip(),
        "alt_loc": line[16:17].strip(),
        "res_name": line[17:20].strip(),
        "chain_id": line[21:22].strip(),
        "res_seq": int(line[22:26]),
        "icode": line[26:27].strip(),
        "x": float(line[30:38]),
        "y": float(line[38:46]),
        "z": float(line[46:54]),
        "occupancy": float(line[54:60].strip() or "1.0"),
        "temp_factor": float(line[60:66].strip() or "0.0"),
        "element": line[76:78].strip(),
        "charge": line[78:80].strip(),
    }


def write_pdb_file(structure_data: Dict[str, Any], file_path: str | Path) -> None:
    """Write structure data to a PDB file.

    Every atom needs ``name``, ``res_name``, ``res_seq``, ``x``, ``y`` and ``z``; a missing one is an error, not an
    alanine C-alpha at the origin. The record name (default ``ATOM``), serial number (default: the position in the
    list), alternate location, chain, insertion code, occupancy (1.00), B-factor (0.00), element and charge are
    optional and blank or conventional when absent.

    Args:
        structure_data: Structure data dictionary
        file_path: Output file path

    Raises:
        ValueError: If an atom lacks a required field

    Example:
        >>> # write_pdb_file(structure_data, "output.pdb")
    """
    file_path = Path(file_path)
    atoms = structure_data.get("atoms", [])

    with open(file_path, "w") as f:
        # Write HEADER
        f.write("HEADER    PROTEIN STRUCTURE                          \n")

        # Write ATOM records
        for position, atom in enumerate(atoms, start=1):
            line = _format_atom_line(atom, default_serial=position)
            f.write(line + "\n")

        # Write END
        f.write("END\n")

    logger.info(f"Written {len(atoms)} atoms to {file_path}")


_REQUIRED_ATOM_FIELDS = ("name", "res_name", "res_seq", "x", "y", "z")


def _format_atom_line(atom: Dict[str, Any], default_serial: int = 0) -> str:
    """Format atom data into PDB ATOM line.

    Args:
        atom: Atom data dictionary
        default_serial: Serial number for an atom without one

    Returns:
        Formatted PDB line

    Raises:
        ValueError: If a required field (name, residue name and number, coordinates) is missing
    """
    missing = [field for field in _REQUIRED_ATOM_FIELDS if field not in atom]
    if missing:
        raise ValueError(f"atom {atom.get('serial', default_serial)} lacks required field(s): {', '.join(missing)}")
    return (
        "{record_name:<6}{serial:>5} {name:<4}{alt_loc:1}{res_name:>3} {chain_id:1}"
        "{res_seq:>4}{icode:1}   {x:>8.3f}{y:>8.3f}{z:>8.3f}{occupancy:>6.2f}"
        "{temp_factor:>6.2f}          {element:>2}{charge:>2}"
    ).format(
        record_name=atom.get("record_name", "ATOM"),
        serial=atom.get("serial", default_serial),
        name=atom["name"],
        alt_loc=atom.get("alt_loc", ""),
        res_name=atom["res_name"],
        chain_id=atom.get("chain_id", ""),
        res_seq=atom["res_seq"],
        icode=atom.get("icode", ""),
        x=atom["x"],
        y=atom["y"],
        z=atom["z"],
        occupancy=atom.get("occupancy", 1.0),
        temp_factor=atom.get("temp_factor", 0.0),
        element=atom.get("element", ""),
        charge=atom.get("charge", ""),
    )


def _split_cif_values(line: str) -> List[str]:
    """Split a line of a CIF loop into its values.

    Values are separated by white space; a value may be quoted with ``'`` or ``"``, and a quote ends the value only
    where it is followed by white space or the end of the line, so ``O5'`` and ``"O5' "`` both keep their quote.
    """
    values: List[str] = []
    position, length = 0, len(line)
    while position < length:
        if line[position].isspace():
            position += 1
            continue
        quote = line[position]
        if quote in "'\"":
            end = position + 1
            while end < length and not (line[end] == quote and (end + 1 == length or line[end + 1].isspace())):
                end += 1
            values.append(line[position + 1 : end])
            position = end + 1
        else:
            end = position
            while end < length and not line[end].isspace():
                end += 1
            values.append(line[position:end])
            position = end
    return values


def _cif_number(text: str, caster: Any, default: Any) -> Any:
    """``caster(text)``, or ``default`` for the CIF markers of an absent value (``?`` and ``.``)."""
    return default if text in ("?", ".") else caster(text)


def parse_mmcif_file(file_path: str | Path) -> Dict[str, Any]:
    """Parse the atom records of an mmCIF file into structured data.

    The ``_atom_site`` loop is read by its column names, so the order and the set of columns do not matter, values
    may be quoted, and a row may span several lines. Only the first model is kept when the file holds several
    (``pdbx_PDB_model_num``). Each atom has the mmCIF fields (``group_PDB``, ``id``, ``type_symbol``, ``label_atom_id``,
    ``label_comp_id``, ``label_asym_id``, ``label_seq_id``, ``Cartn_x``, ``Cartn_y``, ``Cartn_z``; ``label_seq_id`` is
    None where the file has no value) and the PDB-style fields of :func:`parse_pdb_file` (``record_name``, ``serial``,
    ``name``, ``alt_loc``, ``res_name``, ``chain_id``, ``res_seq``, ``icode``, ``x``, ``y``, ``z``, ``occupancy``,
    ``temp_factor``, ``element``), so that the analysis functions work on either. The PDB-style chain and residue
    number are the author's (``auth_asym_id``, ``auth_seq_id``) when the file has them and the label ones otherwise.
    An occupancy or B-factor that the file does not give is NaN.

    Args:
        file_path: Path to mmCIF file

    Returns:
        Dictionary with ``atoms``, ``coordinates`` (n_atoms, 3), ``file_path`` and ``n_atoms``

    Raises:
        FileNotFoundError: If mmCIF file doesn't exist
        ValueError: If the file has no ``_atom_site`` loop with atoms, lacks the coordinate or atom-name columns, a
            row does not have one value per column, or a coordinate is not a number

    Example:
        >>> # cif_data = parse_mmcif_file("protein.cif")
        >>> # isinstance(cif_data, dict)
        >>> # True
    """
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"mmCIF file not found: {file_path}")

    with open(file_path, "r") as f:
        lines = f.read().splitlines()

    columns: List[str] = []
    values: List[str] = []
    index = 0
    while index < len(lines) and not columns:
        if lines[index].strip() == "loop_":
            tags: List[str] = []
            index += 1
            while index < len(lines) and lines[index].strip().startswith("_"):
                tags.append(lines[index].split()[0])
                index += 1
            if tags and tags[0].startswith("_atom_site."):
                columns = [tag[len("_atom_site.") :] for tag in tags]
                while index < len(lines):
                    text = lines[index].strip()
                    if text.startswith(("loop_", "data_", "_")):
                        break
                    if text.startswith(";"):
                        raise ValueError(
                            f"{file_path}: multi-line text values in the _atom_site loop are not supported"
                        )
                    if text and not text.startswith("#"):
                        values.extend(_split_cif_values(text))
                    index += 1
        else:
            index += 1

    if not columns:
        raise ValueError(f"No _atom_site loop found in mmCIF file: {file_path}")
    required = ("Cartn_x", "Cartn_y", "Cartn_z")
    if any(name not in columns for name in required) or not ({"label_atom_id", "auth_atom_id"} & set(columns)):
        raise ValueError(f"{file_path}: the _atom_site loop needs Cartn_x/y/z and an atom name column")
    if len(values) % len(columns) != 0:
        raise ValueError(
            f"{file_path}: the _atom_site loop has {len(values)} values, which do not fill rows of "
            f"{len(columns)} columns"
        )

    atoms: List[Dict[str, Any]] = []
    first_model: Optional[str] = None
    for start in range(0, len(values), len(columns)):
        row = dict(zip(columns, values[start : start + len(columns)]))
        model = row.get("pdbx_PDB_model_num")
        if model is not None:
            first_model = model if first_model is None else first_model
            if model != first_model:
                continue
        x, y, z = (float(row[name]) for name in required)  # a ValueError names a malformed coordinate
        name = row.get("label_atom_id", row.get("auth_atom_id", ""))
        res_name = row.get("label_comp_id", row.get("auth_comp_id", ""))
        chain = row.get("auth_asym_id", row.get("label_asym_id", "."))
        number = row.get("auth_seq_id", row.get("label_seq_id", "."))
        insertion = row.get("pdbx_PDB_ins_code", "?")
        alternate = row.get("label_alt_id", ".")
        symbol = row.get("type_symbol", "?")
        group = row.get("group_PDB", "ATOM")
        serial = _cif_number(row.get("id", "?"), int, 0)
        atoms.append(
            {
                "group_PDB": group,
                "id": serial,
                "type_symbol": symbol,
                "label_atom_id": name,
                "label_comp_id": res_name,
                "label_asym_id": row.get("label_asym_id", chain),
                "label_seq_id": _cif_number(row.get("label_seq_id", "."), int, None),
                "Cartn_x": x,
                "Cartn_y": y,
                "Cartn_z": z,
                "record_name": group,
                "serial": serial,
                "name": name,
                "alt_loc": "" if alternate in ("?", ".") else alternate,
                "res_name": res_name,
                "chain_id": "" if chain in ("?", ".") else chain,
                "res_seq": _cif_number(number, int, 0),
                "icode": "" if insertion in ("?", ".") else insertion,
                "x": x,
                "y": y,
                "z": z,
                "occupancy": _cif_number(row.get("occupancy", "?"), float, float("nan")),
                "temp_factor": _cif_number(row.get("B_iso_or_equiv", "?"), float, float("nan")),
                "element": "" if symbol in ("?", ".") else symbol,
            }
        )

    if not atoms:
        raise ValueError(f"No atoms found in mmCIF file: {file_path}")
    return {
        "atoms": atoms,
        "coordinates": np.array([[atom["x"], atom["y"], atom["z"]] for atom in atoms]),
        "file_path": str(file_path),
        "n_atoms": len(atoms),
    }


def _cif_text(value: Any, empty: str) -> str:
    """A value as written in a CIF loop: ``empty`` when blank, quoted when it has white space or a special start."""
    text = str(value).strip()
    if not text:
        return empty
    if (
        any(character.isspace() for character in text)
        or text[0] in "'\"#$_[];"
        or text.lower().startswith(("data_", "loop_"))
    ):
        quote = '"' if "'" in text else "'"
        return f"{quote}{text}{quote}"
    return text


def convert_pdb_to_mmcif(pdb_file: str | Path, cif_file: str | Path) -> None:
    """Convert PDB file to mmCIF format.

    The output is an ``_atom_site`` loop with one atom per line and the columns ``group_PDB``, ``id``,
    ``type_symbol``, ``label_atom_id``, ``label_alt_id``, ``label_comp_id``, ``label_asym_id``, ``label_seq_id``,
    ``pdbx_PDB_ins_code``, ``Cartn_x``, ``Cartn_y``, ``Cartn_z``, ``occupancy``, ``B_iso_or_equiv``, ``auth_seq_id`` and
    ``auth_asym_id``, which :func:`parse_mmcif_file` reads back. Values the PDB atom lacks are written as the CIF
    markers ``?`` (unknown) and ``.`` (not applicable); no residue, chain or element is invented. Entity, model number
    and the other categories of a full mmCIF file are not written.

    Args:
        pdb_file: Input PDB file path
        cif_file: Output mmCIF file path

    Raises:
        FileNotFoundError: If the PDB file does not exist
        ValueError: If the PDB file holds no valid atoms

    Example:
        >>> # convert_pdb_to_mmcif("input.pdb", "output.cif")
    """
    pdb_file = Path(pdb_file)
    pdb_data = parse_pdb_file(pdb_file)

    cif_file = Path(cif_file)
    columns = (
        "group_PDB",
        "id",
        "type_symbol",
        "label_atom_id",
        "label_alt_id",
        "label_comp_id",
        "label_asym_id",
        "label_seq_id",
        "pdbx_PDB_ins_code",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
        "occupancy",
        "B_iso_or_equiv",
        "auth_seq_id",
        "auth_asym_id",
    )

    block = "".join(character if character.isalnum() or character in "_-" else "_" for character in pdb_file.stem)
    with open(cif_file, "w") as f:
        f.write(f"data_{block or 'structure'}\n#\nloop_\n")
        for column in columns:
            f.write(f"_atom_site.{column}\n")

        for position, atom in enumerate(pdb_data["atoms"], 1):
            chain = _cif_text(atom["chain_id"], ".")
            fields = (
                _cif_text(atom["record_name"], "ATOM"),
                str(atom["serial"] or position),
                _cif_text(atom["element"], "?"),
                _cif_text(atom["name"], "?"),
                _cif_text(atom["alt_loc"], "."),
                _cif_text(atom["res_name"], "?"),
                chain,
                str(atom["res_seq"]),
                _cif_text(atom["icode"], "?"),
                f"{atom['x']:.3f}",
                f"{atom['y']:.3f}",
                f"{atom['z']:.3f}",
                f"{atom['occupancy']:.2f}",
                f"{atom['temp_factor']:.2f}",
                str(atom["res_seq"]),
                chain,
            )
            f.write(" ".join(fields) + "\n")
        f.write("#\n")

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


def read_pdb_ca_coordinates(file_path: str | Path) -> List[Tuple[float, float, float]]:
    """Read CA (alpha carbon) coordinates from a PDB file.

    Args:
        file_path: Path to PDB file

    Returns:
        List of (x, y, z) tuples for each CA atom
    """
    pdb_data = parse_pdb_file(file_path)
    ca_coords = []
    for atom in pdb_data["atoms"]:
        if atom["name"] == "CA":
            ca_coords.append((atom["x"], atom["y"], atom["z"]))
    return ca_coords


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

        # Count residues: chain, number and insertion code
        residues = set((atom["chain_id"], atom["res_seq"], atom["icode"]) for atom in atoms)

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
