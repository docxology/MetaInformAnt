"""Tests: secondary structure prediction is local, uses one published table and never pretends to call a service.

Defects:

- ``predict_secondary_structure`` defaulted to ``method="psipred"``. That method (and ``"jpred"``) posted the sequence
  to hard-coded URLs ("this would need to be updated with actual service"), swallowed the failure and silently returned
  the result of a different, local predictor; one of its ``except`` clauses named ``requests`` which is unbound when the
  import fails. They now raise ``NotImplementedError`` and make no request.
- The "Chou-Fasman" values of ``calculate_ss_propensities`` (``E``: 0.51, ``M``: 0.82, ``K``: 0.81, ...) contradicted
  the table of the window predictor (``E``: 1.53, ...) and both contradicted the published parameters (``E``: 1.51,
  ...); ``"coil": 1.0`` was a made-up constant. There is one table now: Chou & Fasman (1978).
- ``parse_dssp_file`` answered 0.0 for an accessibility or angle that the line does not hold. It answers NaN.
- ``predict_transmembrane_regions`` returned one overlapping "region" per qualifying window.
- Three InterPro functions raise ``NotImplementedError`` but their docstrings promised a return value (and one example
  could not run).

Oracles: the published Chou-Fasman parameters as tabulated in the literature, hand-built DSSP lines whose fields are
placed at the documented columns, and the Kyte-Doolittle hydropathy of a poly-leucine stretch (3.8 per residue).
"""

from __future__ import annotations

import math

import pytest

from metainformant.protein.structure import secondary
from metainformant.protein.structure.secondary import (
    calculate_ss_propensities,
    parse_dssp_file,
    predict_secondary_structure,
    predict_transmembrane_regions,
)

# ---------------------------------------------------------------------------
# Methods that need a service are refused, not substituted
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("method", ["psipred", "jpred"])
def test_methods_that_need_an_external_service_are_refused(method: str) -> None:
    with pytest.raises(NotImplementedError, match=method.upper()):
        predict_secondary_structure("MKTIIALSYIFCLVFA", method=method)


def test_an_unknown_method_is_still_a_value_error() -> None:
    with pytest.raises(ValueError, match="Unsupported method"):
        predict_secondary_structure("ACDE", method="invalid")


def test_the_default_method_runs_locally_and_is_the_window_method() -> None:
    sequence = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV"

    default = predict_secondary_structure(sequence)

    assert default == predict_secondary_structure(sequence, method="window")
    assert len(default) == len(sequence) and set(default) <= {"H", "E", "C"}
    assert default == predict_secondary_structure(sequence)  # deterministic
    assert predict_secondary_structure("") == []  # an empty sequence needs no method


# ---------------------------------------------------------------------------
# One published parameter table
# ---------------------------------------------------------------------------

# Chou & Fasman (1978), Adv Enzymol 47:45-148, as tabulated in the literature: (P_alpha, P_beta, P_turn).
PUBLISHED = {
    "E": (1.51, 0.37, 0.74),
    "M": (1.45, 1.05, 0.60),
    "A": (1.42, 0.83, 0.66),
    "V": (1.06, 1.70, 0.50),
    "I": (1.08, 1.60, 0.47),
    "Y": (0.69, 1.47, 1.14),
    "G": (0.57, 0.75, 1.56),
    "P": (0.57, 0.55, 1.52),
    "N": (0.67, 0.89, 1.56),
    "D": (1.01, 0.54, 1.46),
}


@pytest.mark.parametrize("residue, parameters", sorted(PUBLISHED.items()))
def test_the_propensities_are_the_published_parameters(residue: str, parameters) -> None:
    propensities = calculate_ss_propensities(residue)[residue]

    assert (propensities["helix"], propensities["sheet"], propensities["turn"]) == parameters
    assert "coil" not in propensities  # the constant 1.0 that stood for it was not a parameter


def test_the_table_covers_the_twenty_amino_acids_once() -> None:
    assert sorted(secondary.CHOU_FASMAN_PARAMETERS) == sorted("ACDEFGHIKLMNPQRSTVWY")
    table = secondary.CHOU_FASMAN_PARAMETERS
    assert all(len(values) == 3 and all(0.0 < value < 2.0 for value in values) for values in table.values())
    assert set(calculate_ss_propensities("XBZ*-")) == set()  # letters without parameters are left out


def test_the_preferred_conformation_follows_the_larger_of_helix_and_sheet() -> None:
    propensities = calculate_ss_propensities("EV")

    assert propensities["E"]["preferred"] == "helix"
    assert propensities["V"]["preferred"] == "sheet"


def test_the_residue_by_residue_method_follows_the_largest_parameter() -> None:
    """E, A, M and K favour the helix; V, I, Y, F and T the sheet; G, P, N, D and S the turn (coil)."""
    assert predict_secondary_structure("EAMK", method="simple") == ["H"] * 4
    assert predict_secondary_structure("VIYFT", method="simple") == ["E"] * 5
    assert predict_secondary_structure("GPNDS", method="simple") == ["C"] * 5
    assert predict_secondary_structure("XB*", method="simple") == ["C"] * 3  # no parameters: coil


def test_the_window_method_calls_helix_sheet_and_coil_from_the_window_means() -> None:
    assert predict_secondary_structure("EEEEEEEEEEEE") == ["H"] * 12  # mean helix parameter 1.51
    assert predict_secondary_structure("VVVVVVVVVVVV") == ["E"] * 12  # mean sheet parameter 1.70
    assert predict_secondary_structure("GGGGGGGGGGGG") == ["C"] * 12  # neither mean exceeds 1
    assert predict_secondary_structure("PPPPPPPPPPPP") == ["C"] * 12
    # an alternating stretch of helix and sheet formers has no helix or sheet that lasts three residues
    assert len(predict_secondary_structure("EVEVEVEVEV")) == 10


def test_residues_without_parameters_do_not_count_in_a_window() -> None:
    """The former code gave unknown letters a made-up 0.5 and so dragged a helix window below its threshold."""
    padded = "EEEE" + "X" * 3 + "EEEE"

    assert predict_secondary_structure(padded)[:4] == ["H"] * 4
    assert predict_secondary_structure("XXXXXXXX") == ["C"] * 8  # a window without any known residue is coil


# ---------------------------------------------------------------------------
# DSSP
# ---------------------------------------------------------------------------


def dssp_line(index: int, number: int, chain: str, aa: str, structure: str, acc: int, phi: float, psi: float) -> str:
    """A DSSP residue line of 136 characters with the fields at their documented columns."""
    line = [" "] * 136

    def put(first_column: int, text: str) -> None:  # 1-based column, as in the DSSP documentation
        for offset, character in enumerate(text):
            line[first_column - 1 + offset] = character

    put(1, f"{index:5d}")
    put(6, f"{number:5d}")
    put(12, chain)
    put(14, aa)
    put(17, structure)
    put(35, f"{acc:4d}")
    put(104, f"{phi:6.1f}")
    put(110, f"{psi:6.1f}")
    return "".join(line)


def test_dssp_numbers_are_read_from_their_columns() -> None:
    header = "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC"
    content = "\n".join(
        [
            header,
            dssp_line(1, 1, "A", "M", " ", 140, 360.0, 135.0),
            dssp_line(2, 2, "A", "Q", "E", 53, -63.4, 140.2),
            dssp_line(3, 3, "A", "I", "H", 7, -58.3, -45.1),
        ]
    )

    result = parse_dssp_file(content)

    assert result["sequence"] == "MQI"
    assert result["secondary_structure"] == ["C", "E", "H"]
    assert result["accessibility"] == [140.0, 53.0, 7.0]
    assert result["phi_angles"] == [360.0, -63.4, -58.3]  # 360.0 is DSSP's marker of an undefined angle: kept
    assert result["psi_angles"] == [135.0, 140.2, -45.1]


def test_numbers_that_a_line_does_not_hold_are_nan_never_zero() -> None:
    """The former code answered 0.0, a valid-looking angle and an exposed-looking accessibility of nothing."""
    short = "    1    1 A M  H  <         0   0  100"  # no angle columns at all
    content = (
        "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC\n" + short + "\n" + dssp_line(2, 2, "A", "K", "H", 10, -57.0, -47.0)
    )

    result = parse_dssp_file(content)

    assert math.isnan(result["phi_angles"][0]) and math.isnan(result["psi_angles"][0])
    assert result["phi_angles"][1:] == [-57.0] and result["psi_angles"][1:] == [-47.0]
    full = dssp_line(1, 1, "A", "M", "H", 0, -57.0, -47.0)
    not_a_number = "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC\n" + full[:34] + "   X" + full[38:]
    assert math.isnan(parse_dssp_file(not_a_number)["accessibility"][0])


# ---------------------------------------------------------------------------
# Transmembrane regions
# ---------------------------------------------------------------------------


def test_a_hydrophobic_stretch_is_one_region_not_one_per_window() -> None:
    regions = predict_transmembrane_regions("L" * 25)  # seven windows of 19 qualify, and overlap

    assert len(regions) == 1
    assert (regions[0]["start"], regions[0]["end"], regions[0]["length"]) == (0, 24, 25)
    assert regions[0]["avg_hydrophobicity"] == pytest.approx(3.8)  # Kyte-Doolittle value of leucine
    assert regions[0]["sequence"] == "L" * 25


def test_regions_are_separated_where_the_hydropathy_drops() -> None:
    segment = "L" * 22
    linker = "D" * 25  # hydrophilic: no window over it qualifies
    sequence = segment + linker + segment

    regions = predict_transmembrane_regions(sequence)

    # a window of 19 residues needs 14 leucines to average above 1.5, so each region reaches 4 residues into the linker
    assert [(region["start"], region["end"]) for region in regions] == [(0, 26), (42, 68)]
    assert all(region["length"] == 27 for region in regions)
    assert regions[0]["avg_hydrophobicity"] == pytest.approx((22 * 3.8 - 5 * 3.5) / 27)


def test_no_region_in_a_hydrophilic_or_short_sequence() -> None:
    assert predict_transmembrane_regions("D" * 60) == []
    assert predict_transmembrane_regions("LLLLL") == []
    assert predict_transmembrane_regions("") == []
    assert predict_transmembrane_regions("P" * 40) == []


# ---------------------------------------------------------------------------
# Functions that raise instead of returning must say so
# ---------------------------------------------------------------------------


def test_the_unimplemented_interpro_functions_document_that_they_raise() -> None:
    from metainformant.protein.database import interpro

    for function, arguments in (
        (interpro.get_interpro_hierarchy, ("IPR000001",)),
        (interpro.get_interpro_statistics, ()),
        (interpro.find_similar_interpro_entries, ("IPR000001",)),
    ):
        with pytest.raises(NotImplementedError):
            function(*arguments)
        documentation = function.__doc__ or ""
        assert "not implemented" in documentation.lower()
        assert "NotImplementedError" in documentation
        assert ">>> stats" not in documentation  # the example that called the function and expected a result
