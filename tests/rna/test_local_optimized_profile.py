"""Native qualification of the proposed profile; private fixtures only."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import threading
import shlex
import time
import pytest
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, urlsplit

from metainformant.rna.engine.acquisition_manifest import sha256_file
from metainformant.rna.engine.acquisition_worker import run_manifest


@pytest.mark.parametrize('acquisitions,quantifiers,threads_per_sample', [(4,4,2),(8,8,1)])
def test_eight_paired_acquisitions_profile_and_native_replay(
    native_acquisition_factory, monkeypatch, tmp_path, acquisitions, quantifiers, threads_per_sample,
):
    accessions = tuple(f"SRR{50100 + i}" for i in range(8))
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
    snapshot_path = tmp_path / "snapshot.json"
    snapshot = json.loads(snapshot_path.read_text())
    for row in snapshot["input_files"]:
        row["sha256"] = sha256_file(tmp_path / row["path"])
    snapshot_path.write_text(json.dumps(snapshot))
    for key, value in {
        "AMALGKIT_LOCAL_FRONTIER_SIZE": "48",
        "AMALGKIT_LOCAL_RAW_RESERVATION_BYTES": "180000",
        "AMALGKIT_PIPELINE_ACQUISITION_SLOTS": str(acquisitions),
        "AMALGKIT_ENA_MAX_CONNECTIONS": str(2 * acquisitions),
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
                    if state["file_active"] == 2 * acquisitions:
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
                        failures.append("configured simultaneous native transfers were not reached")
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
    kwargs = dict(manifest_path=fixture.manifest, data_root=fixture.data, config_dir=fixture.config_dir,
                  workers=12, threads=8, fastq_threads=2, compression_threads=1,
                  validation_slots=2, quant_slots=quantifiers, fasterq_slots=1, max_in_flight=48)
    try:
        started = time.monotonic()
        first = run_manifest(**kwargs)
        elapsed = time.monotonic() - started
        assert not failures, failures
        assert state["file_peak"] == 2 * acquisitions and state["api_peak"] <= 2, state
        assert state["gets"] == 16, state
        assert first["counts"]["newly_quantified"] == 8, first
        files = [*fixture.work.glob("getfastq/*/*.gz"), *fixture.work.glob("quant/*/*abundance.tsv")]
        assert len(files) == 24
        before = {str(path): (sha256_file(path), path.stat().st_mtime_ns) for path in files}
        for path in fixture.work.glob("quant/*/*abundance.tsv"):
            with path.open() as handle:
                abundance = list(csv.DictReader(handle, delimiter="\t"))
            assert len(abundance) == 3
            assert abs(sum(float(row["tpm"]) for row in abundance) - 1e6) < 2
            assert sum(float(row["est_counts"]) for row in abundance) > 80
        infos = list(fixture.work.glob("quant/*/*run_info.json"))
        assert len(infos) == 8
        native_calls = [shlex.split(json.loads(path.read_text())["call"]) for path in infos]
        assert all(int(call[call.index("--threads") + 1]) == threads_per_sample
                   for call in native_calls), native_calls
        second = run_manifest(**kwargs)
        assert second["counts"]["reused"] == 8, second
        assert state["gets"] == 16, state
        assert before == {str(path): (sha256_file(path), path.stat().st_mtime_ns) for path in files}
        (tmp_path / "profile-proof.json").write_text(json.dumps({"acquisitions": acquisitions, "quantifiers": quantifiers, "threads_per_sample": threads_per_sample, "elapsed_seconds": elapsed, "state": state, "first_counts": first["counts"], "replay_counts": second["counts"], "raw_and_abundance_preserved": before}, indent=2))
    finally:
        release.set()
        server.shutdown()
        server.server_close()
        thread.join(5)
        assert not thread.is_alive()
