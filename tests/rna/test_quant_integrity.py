"""Invalid real-file quant facts cannot become durable receipts."""

from __future__ import annotations
import json
from pathlib import Path
import pytest
from metainformant.rna.engine.durable_quant import DirectoryStore, lock_quantification
from metainformant.rna.engine.provenance import write_quant_provenance


@pytest.fixture
def quant_sample(tmp_path: Path) -> tuple[Path, Path]:
    config = tmp_path / "config.yaml"
    config.write_text("species_list: [Test_species]\n")
    sample = tmp_path / "SRR123"
    sample.mkdir()
    (sample / "SRR123_abundance.tsv").write_text(
        "target_id\tlength\teff_length\test_counts\ttpm\ntranscript1\t100\t70.5\t3.25\t1000000\n"
    )
    (sample / "SRR123_run_info.json").write_text(
        json.dumps({"n_processed": 100, "n_pseudoaligned": 10, "n_targets": 1})
    )
    write_quant_provenance(
        sample,
        species="test_species",
        run_accession="SRR123",
        config_path=config,
        command=["amalgkit", "quant"],
    )
    return sample, config


@pytest.mark.parametrize(
    "field,value",
    [
        ("n_processed", True),
        ("n_processed", 100.5),
        ("n_pseudoaligned", True),
        ("n_pseudoaligned", 10.5),
        ("n_pseudoaligned", 101),
        ("n_targets", True),
        ("n_targets", 1.5),
        ("n_targets", 2),
    ],
)
def test_impossible_run_facts_never_publish(
    quant_sample: tuple[Path, Path], tmp_path: Path, field: str, value: float
) -> None:
    sample, _ = quant_sample
    path = sample / "SRR123_run_info.json"
    info = json.loads(path.read_text())
    info[field] = value
    path.write_text(json.dumps(info))
    store = DirectoryStore(tmp_path / "store")
    with pytest.raises(ValueError):
        lock_quantification(store, "cohort", "test_species", "SRR123", sample)
    assert not list(store.root.rglob("*.json"))


@pytest.mark.parametrize("column", ["length", "eff_length"])
def test_zero_feature_lengths_never_publish(
    quant_sample: tuple[Path, Path], tmp_path: Path, column: str
) -> None:
    sample, config = quant_sample
    row = ["transcript1", "100", "70.5", "3.25", "1000000"]
    row[["target_id", "length", "eff_length", "est_counts", "tpm"].index(column)] = "0"
    (sample / "SRR123_abundance.tsv").write_text(
        "target_id\tlength\teff_length\test_counts\ttpm\n" + "\t".join(row) + "\n"
    )
    write_quant_provenance(
        sample,
        species="test_species",
        run_accession="SRR123",
        config_path=config,
        command=["amalgkit", "quant"],
    )
    with pytest.raises(ValueError):
        lock_quantification(
            DirectoryStore(tmp_path / "store"),
            "cohort",
            "test_species",
            "SRR123",
            sample,
        )


@pytest.mark.parametrize("column", ["est_counts", "tpm"])
def test_finite_rows_with_overflowing_totals_never_publish(
    quant_sample: tuple[Path, Path], tmp_path: Path, column: str
) -> None:
    sample, config = quant_sample
    counts, tpm = ("1e308", "500000") if column == "est_counts" else ("1.5", "1e308")
    (sample / "SRR123_abundance.tsv").write_text(
        "target_id\tlength\teff_length\test_counts\ttpm\n"
        + "".join(f"transcript{i}\t100\t70.5\t{counts}\t{tpm}\n" for i in [1, 2])
    )
    (sample / "SRR123_run_info.json").write_text(
        json.dumps({"n_processed": 100, "n_pseudoaligned": 10, "n_targets": 2})
    )
    write_quant_provenance(
        sample,
        species="test_species",
        run_accession="SRR123",
        config_path=config,
        command=["amalgkit", "quant"],
    )
    with pytest.raises(ValueError):
        lock_quantification(
            DirectoryStore(tmp_path / "store"),
            "cohort",
            "test_species",
            "SRR123",
            sample,
        )


def test_fractional_estimates_with_matching_targets_are_valid(
    quant_sample: tuple[Path, Path], tmp_path: Path
) -> None:
    sample, _ = quant_sample
    receipt = lock_quantification(
        DirectoryStore(tmp_path / "store"), "cohort", "test_species", "SRR123", sample
    )
    assert receipt["feature_count"] == 1
