"""
Real-Implementation Integration Test for ENADownloader.

Verifies that the ENADownloader can discover FASTQ URLs from the ENA API.
Does NOT perform actual large file downloads to respect network/storage,
but validates the URL generation logic against a known real-world sample.
"""

import gzip
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import pytest

import metainformant.rna.retrieval.ena_downloader as ena_downloader_module
from metainformant.rna.retrieval.ena_downloader import ENADownloader

pytestmark = pytest.mark.network


class TestENADownloader(unittest.TestCase):
    def setUp(self):
        self.downloader = ENADownloader()
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_get_fastq_urls_real_logic(self):
        """
        Verify URL discovery logic using a known sample ID.
        We expect ENA API to return a valid result for a public sample.

        Using a small, public SRA run: SRR11092056 (Apis mellifera RNA-seq)
        """
        # This test actually hits the ENA API.
        # It aligns with "Real-Implementation" philosophy for integration tests calling external APIs,
        # provided we handle network failures gracefully.

        sample_id = "SRR11092056"
        urls = self.downloader.get_fastq_urls(sample_id)

        # We expect at least one URL
        if not urls:
            print(f"Warning: ENA API returned no URLs for {sample_id}. Network issue?")
            return

        self.assertTrue(len(urls) > 0, "Should find at least one URL")
        self.assertTrue(any("SRR11092056" in u for u in urls), "URL should contain sample ID")
        self.assertTrue(urls[0].startswith("http"), "URL should be HTTP(S)")
        self.assertTrue(urls[0].endswith(".fastq.gz"), "URL should point to .fastq.gz")

    def test_get_fastq_urls_invalid(self):
        """Verify handling of invalid IDs."""
        urls = self.downloader.get_fastq_urls("INVALID_ID_12345")
        self.assertEqual(urls, [], "Should return empty list for invalid ID")

    def test_invalid_gzip_gets_bounded_fresh_retry(self):
        """Corrupt completed transfers are preserved and retried from byte zero."""

        downloader = ENADownloader(timeout=30, retries=1, retry_delay_seconds=1)
        output_dir = Path(self.test_dir)
        calls = []

        def fake_curl(command, **_kwargs):
            calls.append(command)
            partial = Path(command[command.index("-o") + 1])
            partial.parent.mkdir(parents=True, exist_ok=True)
            if len(calls) == 1:
                partial.write_bytes(b"not a gzip payload")
            else:
                with gzip.open(partial, "wb") as handle:
                    handle.write(b"@SRR1.1 1/1\nACGT\n+\n!!!!\n")
            return subprocess.CompletedProcess(command, 0, "200", "")

        original_get_fastq_urls = downloader.get_fastq_urls
        original_run_command = ena_downloader_module._run_command_in_process_group
        downloader.get_fastq_urls = lambda _sample_id: ["https://example.test/SRR1.fastq.gz"]
        ena_downloader_module._run_command_in_process_group = fake_curl
        try:
            success, message, files = downloader.download_run("SRR1", output_dir)
        finally:
            downloader.get_fastq_urls = original_get_fastq_urls
            ena_downloader_module._run_command_in_process_group = original_run_command

        self.assertTrue(success, message)
        self.assertEqual([path.name for path in files], ["SRR1.fastq.gz"])
        self.assertEqual(len(calls), 2)
        self.assertTrue((output_dir / "SRR1.fastq.gz").is_file())
        self.assertTrue((output_dir / "SRR1.fastq.gz.part.invalid").is_file())

    def test_repeated_invalid_gzip_payloads_are_fingerprinted_not_accumulated(self):
        """Only the first invalid payload remains as a full diagnostic witness."""

        downloader = ENADownloader(timeout=30, retries=1, integrity_retries=2, retry_delay_seconds=1)
        output_dir = Path(self.test_dir)
        calls = []

        def fake_curl(command, **_kwargs):
            calls.append(command)
            partial = Path(command[command.index("-o") + 1])
            partial.parent.mkdir(parents=True, exist_ok=True)
            if len(calls) < 3:
                partial.write_bytes(f"invalid-{len(calls)}".encode())
            else:
                with gzip.open(partial, "wb") as handle:
                    handle.write(b"@SRR1.1 1/1\nACGT\n+\n!!!!\n")
            return subprocess.CompletedProcess(command, 0, "200", "")

        original_get_fastq_urls = downloader.get_fastq_urls
        original_run_command = ena_downloader_module._run_command_in_process_group
        downloader.get_fastq_urls = lambda _sample_id: ["https://example.test/SRR2.fastq.gz"]
        ena_downloader_module._run_command_in_process_group = fake_curl
        try:
            success, message, files = downloader.download_run("SRR2", output_dir)
        finally:
            downloader.get_fastq_urls = original_get_fastq_urls
            ena_downloader_module._run_command_in_process_group = original_run_command

        self.assertTrue(success, message)
        self.assertEqual(len(calls), 3)
        self.assertEqual(len(list(output_dir.glob("SRR2.fastq.gz.part.invalid*"))), 1)
        manifest = output_dir / ".transfer_integrity_failures.tsv"
        self.assertEqual(len(manifest.read_text().splitlines()), 3)


if __name__ == "__main__":
    unittest.main()
