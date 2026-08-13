"""Real subprocess tests for the GWAS SRA campaign environment."""

from __future__ import annotations

import os
from pathlib import Path

from metainformant.gwas.data.sra_download import download_sra_run


def test_gwas_sra_tool_inherits_campaign_local_environment(tmp_path: Path, monkeypatch) -> None:
    """GWAS SRA conversion must inherit isolated settings and cache paths."""

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    capture = tmp_path / "environment.txt"
    tool = bin_dir / "fasterq-dump"
    tool.write_text(
        "#!/bin/sh\n"
        'printf \'%s\\n\' "NCBI_SETTINGS=$NCBI_SETTINGS" "TMPDIR=$TMPDIR" '
        '"VDB_CONFIG=$VDB_CONFIG" > "$CAPTURE_FILE"\n'
        "out=''\n"
        'while [ "$#" -gt 0 ]; do\n'
        '  if [ "$1" = "--outdir" ]; then out=$2; shift 2; else shift; fi\n'
        "done\n"
        'mkdir -p "$out"\n'
        "i=0\n"
        "while [ $i -lt 80 ]; do printf '@read\\nACGT\\n+\\nIIII\\n' >> \"$out/SRR123_1.fastq\"; i=$((i + 1)); done\n",
        encoding="utf-8",
    )
    tool.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("CAPTURE_FILE", str(capture))

    result = download_sra_run("SRR123", tmp_path / "campaign")

    assert result == tmp_path / "campaign" / "SRR123"
    captured = capture.read_text(encoding="utf-8")
    for key in ("NCBI_SETTINGS", "TMPDIR", "VDB_CONFIG"):
        value = next(line.split("=", 1)[1] for line in captured.splitlines() if line.startswith(f"{key}="))
        assert Path(value).is_relative_to((tmp_path / "campaign").resolve())
