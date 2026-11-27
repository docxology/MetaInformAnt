"""Export the software-methods study and replay its deterministic table benchmark."""

from __future__ import annotations

import argparse
import shutil
from hashlib import sha256
from pathlib import Path

import matplotlib
from matplotlib import pyplot as plt
from metainformant.core import io
from metainformant.core.data import validation
from metainformant.core.utils.errors import ValidationError

PROFILE = Path(__file__).resolve().parent
SCALES = (32, 256, 2048)
PERIOD = 17
COLUMNS = {"record": str, "count": int}
CONTROLS = ("valid", "schema", "lexical", "range", "duplicate")


def check_table(path: Path) -> tuple[str, list[dict[str, str]], int | None]:
    """Return the first violated contract or the successfully parsed table and total."""
    rows = list(io.read_delimited(path, delimiter="\t"))
    seen = set()
    total = 0
    for row in rows:
        if list(row) != list(COLUMNS):
            return "schema", rows, None
        try:
            count = int(row["count"])
        except ValueError:
            return "lexical", rows, None
        validation.validate_schema({"record": row["record"], "count": count}, COLUMNS)
        try:
            validation.validate_range(count, min_val=0, name="count")
        except ValidationError:
            return "range", rows, None
        identifier = row["record"]
        if identifier in seen:
            return "duplicate", rows, None
        seen.add(identifier)
        total += count
    return "valid", rows, total


def benchmark(data_dir: Path) -> list[dict]:
    results = []
    for size in SCALES:
        baseline = [
            {"record": f"record_{index:04d}", "count": str(index % PERIOD)}
            for index in range(size)
        ]
        cycles, remainder = divmod(size, PERIOD)
        expected_total = cycles * PERIOD * (PERIOD - 1) // 2 + remainder * (remainder - 1) // 2
        for expected in CONTROLS:
            rows = [row.copy() for row in baseline]
            if expected == "schema":
                rows = [{**row, "extra": "unrecognized"} for row in rows]
            elif expected == "lexical":
                rows[-1]["count"] = "1.5"
            elif expected == "range":
                rows[-1]["count"] = "-1"
            elif expected == "duplicate":
                rows[-1]["record"] = rows[0]["record"]

            fixture = data_dir / f"{size}_{expected}.tsv"
            io.write_delimited(rows, fixture, delimiter="\t")
            actual, parsed, total = check_table(fixture)
            if actual != expected or (expected == "valid" and (parsed != baseline or total != expected_total)):
                raise ValueError(f"Unexpected outcome for {fixture.name}: {actual}")
            payload = fixture.read_bytes()
            results.append(
                {
                    "rows": size,
                    "artifact": fixture.name,
                    "expected": expected,
                    "observed": actual,
                    "expected_total": expected_total if expected == "valid" else None,
                    "observed_total": total,
                    "bytes": len(payload),
                    "sha256": sha256(payload).hexdigest(),
                }
            )
    return results


def render_evidence(results: list[dict], path: Path) -> None:
    valid = [record for record in results if record["expected"] == "valid"]
    scales = [record["rows"] for record in valid]
    matplotlib.use("Agg")
    fig, (checksums, controls) = plt.subplots(1, 2, figsize=(9.5, 3.6), layout="constrained")
    checksums.plot(scales, [record["expected_total"] for record in valid], marker="o", label="Oracle")
    checksums.scatter(
        scales, [record["observed_total"] for record in valid],
        marker="s", facecolors="none", edgecolors="#bd4b12", s=130, label="Observed",
    )
    checksums.set(xscale="log", yscale="log", xlabel="Rows in valid table", ylabel="Sum of counts")
    checksums.set_xticks(scales, labels=[str(scale) for scale in scales])
    checksums.legend()
    malformed = CONTROLS[1:]
    outcomes = {(record["expected"], record["rows"]): record["observed"] for record in results}
    grid = [[outcomes[control, scale] == control for scale in scales] for control in malformed]
    controls.imshow(grid, cmap="Blues", aspect="auto", vmin=0, vmax=1.3)
    controls.set_xticks(range(len(scales)), labels=[str(scale) for scale in scales])
    controls.set_yticks(range(len(malformed)), labels=malformed)
    controls.set(xlabel="Rows in malformed table", title="Expected rejection by control")
    for row_index, row in enumerate(grid):
        for column_index, matched in enumerate(row):
            controls.text(column_index, row_index, "Matched" if matched else "Failed",
                          color="white", va="center", ha="center", fontsize=8)
    fig.savefig(path, dpi=175)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path, help="New Template project directory")
    destination = parser.parse_args().destination
    if destination.exists():
        parser.error("destination must be a new directory")

    source, data, figures = (
        destination / "docs/manuscript",
        destination / "output/data",
        destination / "output/figures",
    )
    for directory in (source, data, figures):
        directory.mkdir(parents=True)
    for item in sorted(PROFILE.glob("[0-9][0-9]_*.md")) + [PROFILE / "config.yaml", PROFILE / "references.bib"]:
        shutil.copyfile(item, source / item.name)

    results = benchmark(data)
    receipt = data / "table_contract.json"
    io.dump_json(
        {
            "valid_tables": sum(case["observed"] == "valid" for case in results),
            "rejected_tables": sum(case["observed"] != "valid" for case in results),
            "cases": results,
        },
        receipt,
        indent=2,
    )
    checked = io.load_json(receipt)
    if checked["cases"] != results:
        raise ValueError("JSON receipt changed during serialization")
    render_evidence(checked["cases"], figures / "table_contract.png")


if __name__ == "__main__":
    main()
