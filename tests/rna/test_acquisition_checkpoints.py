from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.rna.engine.acquisition_worker import run_manifest


@pytest.mark.external_tool
def test_checkpoint_precedes_next_admission_and_failure_preserves_inputs(
    tmp_path: Path, native_acquisition_factory
) -> None:
    fixture = native_acquisition_factory(
        accessions=("SRR123", "SRR124"), sequence_length=1000, vary_read_positions=True
    )
    data, work, config_dir, manifest = (
        fixture.data,
        fixture.work,
        fixture.config_dir,
        fixture.manifest,
    )
    checkpoints: list[str] = []

    def checkpoint(species: str, accession: str) -> None:
        assert species == "test_species"
        assert (work / "quant" / accession / f"{accession}_abundance.tsv").is_file()
        if accession == "SRR123":
            assert not (work / "quant/SRR124/SRR124_abundance.tsv").exists()
        checkpoints.append(accession)

    kwargs = {
        "manifest_path": manifest,
        "data_root": data,
        "config_dir": config_dir,
        "workers": 2,
        "threads": 2,
        "fastq_threads": 1,
        "compression_threads": 1,
        "validation_slots": 1,
        "quant_slots": 1,
        "fasterq_slots": 1,
        "max_in_flight": 1,
    }
    summary = run_manifest(**kwargs, on_quantified=checkpoint)
    assert checkpoints == ["SRR123", "SRR124"]
    assert summary["counts"]["newly_quantified"] == 2
    stop_file = tmp_path / "STOP"

    def stop_after_checkpoint(species: str, accession: str) -> None:
        stop_file.write_text(f"{species}/{accession}")

    stopped = run_manifest(**kwargs, on_quantified=stop_after_checkpoint, should_stop=stop_file.exists)
    assert stopped["counts"] == {"quantified": 1, "reused": 1}
    assert stopped["unsubmitted_task_count"] == 1
    assert stopped["stopped_early"] is True
    retained = work / "getfastq/SRR123/SRR123.fastq.gz"
    before = retained.read_bytes()
    published = (data / "cloud_worker_result.json").read_bytes()

    def unavailable_checkpoint(species: str, accession: str) -> None:
        raise OSError(f"Retention unavailable for {species}/{accession}")

    with pytest.raises(OSError, match="Retention unavailable"):
        run_manifest(**kwargs, on_quantified=unavailable_checkpoint)
    assert retained.read_bytes() == before
    assert (data / "cloud_worker_result.json").read_bytes() == published
