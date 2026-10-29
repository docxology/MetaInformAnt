"""Native qualification of the proposed profile; private fixtures only."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, urlsplit

from metainformant.rna.engine.acquisition_manifest import sha256_file
from metainformant.rna.engine.acquisition_worker import run_manifest


def test_overlay_twelve_acquisitions_twenty_four_transfers_q8_and_native_replay(
    native_acquisition_factory,
    monkeypatch,
    tmp_path,
):
    accessions = tuple(f"SRR{50100 + i}" for i in range(12))
    fixture = native_acquisition_factory(accessions=accessions, paired_metadata=True)
    sequences = (tmp_path / "transcripts.fa").read_text().splitlines()[1::2]
    payloads = {}
    for accession in accessions:
        raw = fixture.work / "getfastq" / accession
        (raw / f"{accession}.fastq.gz").rename(raw / "preserved-original.bin")
        for mate in (1, 2):
            records = []
            for number in range(90):
                gene = sequences[number % 3]
                seq = gene[100:175] if mate == 1 else gene[200:275].translate(str.maketrans("ACGT", "TGCA"))[::-1]
                records.append(f"@{accession}.{number + 1}/{mate}\n{seq}\n+\n{'I' * 75}\n")
            payloads[f"{accession}_{mate}.fastq.gz"] = gzip.compress("".join(records).encode(), mtime=0)
    with fixture.metadata.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields, rows = reader.fieldnames, list(reader)
    for row in rows:
        row.update(total_bases="13500", spot_length="150")
    with fixture.metadata.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    tasks = [json.loads(line) for line in fixture.manifest.read_text().splitlines()]
    overlay_rows = []
    for task in tasks:
        task["fastq_bytes"] = 0
        accession = task["accession"]
        names = [f"{accession}_{mate}.fastq.gz" for mate in (1, 2)]
        overlay_rows.append(
            {
                "task_id": task["task_id"],
                "accession": accession,
                "expected_paired": True,
                "reference_index_sha256": task["reference_index_sha256"],
                "fastq_bytes": sum(len(payloads[name]) for name in names),
                "files": [
                    {
                        "name": name,
                        "url": f"https://ftp.sra.ebi.ac.uk/vol1/fastq/{accession}/{name}",
                        "bytes": len(payloads[name]),
                        "md5": hashlib.md5(payloads[name], usedforsecurity=False).hexdigest(),
                    }
                    for name in names
                ],
            }
        )
    fixture.manifest.write_text("".join(json.dumps(task) + "\n" for task in tasks))
    frozen_before = fixture.manifest.read_bytes()
    overlay_path = tmp_path / "size-overlay.json"
    overlay_path.write_text(
        json.dumps(
            {
                "schema": "metainformant.rna.size_overlay.v1",
                "frozen_manifest_sha256": sha256_file(fixture.manifest),
                "records": overlay_rows,
            }
        )
    )
    overlay_sha = sha256_file(overlay_path)
    monkeypatch.setenv("AMALGKIT_LOCAL_SIZE_OVERLAY_PATH", str(overlay_path))
    monkeypatch.setenv("AMALGKIT_LOCAL_SIZE_OVERLAY_SHA256", overlay_sha)
    snapshot_path = tmp_path / "snapshot.json"
    snapshot = json.loads(snapshot_path.read_text())
    snapshot["manifest_sha256"] = sha256_file(fixture.manifest)
    for row in snapshot["input_files"]:
        row["sha256"] = sha256_file(tmp_path / row["path"])
    snapshot_path.write_text(json.dumps(snapshot))
    for key, value in {
        "AMALGKIT_LOCAL_FRONTIER_SIZE": "48",
        "AMALGKIT_LOCAL_RAW_RESERVATION_BYTES": "180000",
        "AMALGKIT_PIPELINE_ACQUISITION_SLOTS": "12",
        "AMALGKIT_ENA_MAX_CONNECTIONS": "24",
        "AMALGKIT_ENA_API_CONNECTIONS": "2",
        "AMALGKIT_PIPELINE_ENA_FILE_WORKERS": "2",
        "AMALGKIT_PIPELINE_DEDICATED_RECOVERY": "yes",
        "AMALGKIT_PIPELINE_VALIDATE_FASTQ_RECORDS": "yes",
        "AMALGKIT_PIPELINE_DOWNLOAD_FIRST_PASS_SECONDS": "60",
        "AMALGKIT_LOCAL_ENA_SCRATCH_DIR": "",
        "AMALGKIT_LOCAL_QUANT_SCRATCH_DIR": "",
    }.items():
        monkeypatch.setenv(key, value)
    lock, release = threading.Lock(), threading.Event()
    state = {"file_active": 0, "file_peak": 0, "api_active": 0, "api_peak": 0, "gets": 0}
    failures = []

    class Server(ThreadingHTTPServer):
        request_queue_size = 64

    class Handler(BaseHTTPRequestHandler):
        def do_HEAD(self):
            body = payloads[self.path.lstrip("/")]
            self.send_response(200)
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()

        def do_GET(self):
            api = self.path.startswith("/api?")
            kind = "api" if api else "file"
            with lock:
                state[kind + "_active"] += 1
                state[kind + "_peak"] = max(state[kind + "_peak"], state[kind + "_active"])
                if not api:
                    state["gets"] += 1
                    if state["file_active"] == 24:
                        release.set()
            try:
                if api:
                    accession = parse_qs(urlsplit(self.path).query)["accession"][0]
                    names = [f"{accession}_{mate}.fastq.gz" for mate in (1, 2)]
                    urls = ";".join(f"http://127.0.0.1:{self.server.server_port}/{name}" for name in names)
                    sizes = ";".join(str(len(payloads[name])) for name in names)
                    md5s = ";".join(hashlib.md5(payloads[name], usedforsecurity=False).hexdigest() for name in names)
                    body = f"run_accession\tfastq_ftp\tfastq_bytes\tfastq_md5\n{accession}\t{urls}\t{sizes}\t{md5s}\n".encode()
                else:
                    if not release.wait(30):
                        failures.append("24 simultaneous native transfers were not reached")
                    body = payloads[self.path.lstrip("/")]
                self.send_response(200)
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
            finally:
                with lock:
                    state[kind + "_active"] -= 1

        def log_message(self, format, *args):
            return

    server = Server(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever)
    thread.start()
    monkeypatch.setenv("AMALGKIT_ENA_API_BASE_URL", f"http://127.0.0.1:{server.server_port}/api")
    kwargs = dict(
        manifest_path=fixture.manifest,
        data_root=fixture.data,
        config_dir=fixture.config_dir,
        workers=12,
        threads=8,
        fastq_threads=2,
        compression_threads=1,
        validation_slots=2,
        quant_slots=8,
        fasterq_slots=1,
        max_in_flight=48,
    )
    try:
        first = run_manifest(**kwargs)
        assert first["size_overlay_sha256"] == overlay_sha
        assert first["size_overlay_records"] == 12
        assert fixture.manifest.read_bytes() == frozen_before
        assert all(t["fastq_bytes"] == 0 for t in tasks)
        assert not failures, failures
        assert state["file_peak"] == 24 and state["api_peak"] <= 2, state
        assert state["gets"] == 24, state
        assert first["counts"]["newly_quantified"] == 12, first
        files = [*fixture.work.glob("getfastq/*/*.gz"), *fixture.work.glob("quant/*/*abundance.tsv")]
        assert len(files) == 36
        before = {str(path): (sha256_file(path), path.stat().st_mtime_ns) for path in files}
        for path in fixture.work.glob("quant/*/*abundance.tsv"):
            with path.open() as handle:
                abundance = list(csv.DictReader(handle, delimiter="\t"))
            assert len(abundance) == 3
            assert abs(sum(float(row["tpm"]) for row in abundance) - 1e6) < 2
            assert sum(float(row["est_counts"]) for row in abundance) > 80
        second = run_manifest(**kwargs)
        assert second["counts"]["reused"] == 12, second
        assert second["size_overlay_sha256"] == overlay_sha
        assert fixture.manifest.read_bytes() == frozen_before
        quant_commands = [
            json.loads(path.read_text()) for path in fixture.work.glob("quant/*/.metainformant_quant_provenance.json")
        ]
        assert len(quant_commands) == 12
        for proof in quant_commands:
            command = proof["command"]
            assert command[command.index("--threads") + 1] == "1"

        assert state["gets"] == 24, state
        assert before == {str(path): (sha256_file(path), path.stat().st_mtime_ns) for path in files}
        (tmp_path / "profile-proof.json").write_text(
            json.dumps(
                {
                    "state": state,
                    "first_counts": first["counts"],
                    "replay_counts": second["counts"],
                    "raw_and_abundance_preserved": before,
                },
                indent=2,
            )
        )
    finally:
        release.set()
        server.shutdown()
        server.server_close()
        thread.join(5)
        assert not thread.is_alive()
