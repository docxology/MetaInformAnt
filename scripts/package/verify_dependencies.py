#!/usr/bin/env python3
"""Write a structured dependency verification report for MetaInformAnt."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import socket
import subprocess
import sys
import time
from dataclasses import asdict, dataclass, field
from importlib import util as importlib_util
from pathlib import Path
from typing import Iterable


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT_DIR = REPO_ROOT / "output" / "dependencies"

HOMEBREW_BLAST_CANDIDATES = [
    DEFAULT_OUTPUT_DIR / "tools" / "ncbi-blast-2.17.0+" / "bin" / "blastn",
    Path("/opt/homebrew/opt/blast/bin/blastn"),
    Path("/usr/local/opt/blast/bin/blastn"),
    Path("/opt/homebrew/opt/blast@2.2/bin/blastn"),
    Path("/usr/local/opt/blast@2.2/bin/blastn"),
    Path("/opt/homebrew/Cellar/blast/2.17.0/bin/blastn"),
    Path("/usr/local/Cellar/blast/2.17.0/bin/blastn"),
]

PYTHON_IMPORTS = {
    "numpy": "numpy",
    "pandas": "pandas",
    "pyarrow": "pyarrow",
    "matplotlib": "matplotlib",
    "biopython": "Bio",
    "requests": "requests",
    "pyyaml": "yaml",
    "beautifulsoup4": "bs4",
    "scipy": "scipy",
    "networkx": "networkx",
    "scikit-learn": "sklearn",
    "psutil": "psutil",
    "rich": "rich",
    "statsmodels": "statsmodels",
    "seaborn": "seaborn",
    "numba": "numba",
    "llvmlite": "llvmlite",
    "umap-learn": "umap",
    "scanpy": "scanpy",
    "anndata": "anndata",
    "python-louvain": "community",
    "cdlib": "cdlib",
    "plotly": "plotly",
    "bokeh": "bokeh",
    "altair": "altair",
    "graphviz": "graphviz",
    "amalgkit": "amalgkit",
    "psycopg2-binary": "psycopg2",
    "sqlalchemy": "sqlalchemy",
    "pymongo": "pymongo",
    "redis": "redis",
    "pysam": "pysam",
    "dendropy": "dendropy",
    "cloudscraper": "cloudscraper",
    "dask": "dask",
    "ray": "ray",
    "joblib": "joblib",
    "xgboost": "xgboost",
    "lightgbm": "lightgbm",
    "multiqc": "multiqc",
    "pytest": "pytest",
    "pytest-cov": "pytest_cov",
    "pytest-xdist": "xdist",
    "pytest-benchmark": "pytest_benchmark",
    "pytest-asyncio": "pytest_asyncio",
    "pytest-repeat": "pytest_repeat",
    "pytest-rerunfailures": "pytest_rerunfailures",
    "pytest-clarity": "pytest_clarity",
    "black": "black",
    "isort": "isort",
    "flake8": "flake8",
    "mypy": "mypy",
    "bandit": "bandit",
    "safety": "safety",
    "pre-commit": "pre_commit",
    "sphinx": "sphinx",
    "sphinx-rtd-theme": "sphinx_rtd_theme",
    "sphinx-autodoc-typehints": "sphinx_autodoc_typehints",
    "myst-parser": "myst_parser",
    "jupyter": "jupyter",
    "nbsphinx": "nbsphinx",
    "twine": "twine",
    "build": "build",
}

CLI_CHECKS = {
    "amalgkit": [["amalgkit", "--version"]],
    "kallisto": [["kallisto", "version"]],
    "salmon": [["salmon", "--version"]],
    "fastq-dump": [["fastq-dump", "--version"]],
    "fasterq-dump": [["fasterq-dump", "--version"]],
    "prefetch": [["prefetch", "--version"]],
    "seqkit": [["seqkit", "version"]],
    "samtools": [["samtools", "--version"]],
    "bcftools": [["bcftools", "--version"]],
    "bwa": [["bwa"]],
    "bowtie2": [["bowtie2", "--version"]],
    "hisat2": [["hisat2", "--version"]],
    "fastqc": [["fastqc", "--version"]],
    "multiqc": [["multiqc", "--version"]],
    "parallel-fastq-dump": [["parallel-fastq-dump", "--version"]],
    "muscle": [["muscle", "-version"], ["muscle", "--version"]],
    "blastn": [["blastn", "-version"]],
    "Rscript": [["Rscript", "--version"]],
    "gcloud": [["gcloud", "--version"]],
    "docker": [["docker", "--version"]],
    "psql": [["psql", "--version"]],
    "postgres": [["postgres", "--version"]],
    "redis-server": [["redis-server", "--version"]],
    "mongod": [["mongod", "--version"]],
    "mongosh": [["mongosh", "--version"]],
}

CLI_EXECUTABLE_CANDIDATES = {
    "blastn": HOMEBREW_BLAST_CANDIDATES,
}

R_PACKAGES = [
    "amap",
    "RColorBrewer",
    "colorspace",
    "dendextend",
    "NMF",
    "MASS",
    "pvclust",
    "Rtsne",
    "ggplot2",
    "patchwork",
    "optparse",
    "BiocManager",
    "reshape2",
    "gridExtra",
    "limma",
    "tximport",
    "pcaMethods",
    "edgeR",
    "RUVSeq",
    "sva",
]


@dataclass
class CheckResult:
    name: str
    status: str
    detail: str
    version: str | None = None
    path: str | None = None


@dataclass
class VerificationReport:
    generated_at: str
    python: list[CheckResult] = field(default_factory=list)
    cli: list[CheckResult] = field(default_factory=list)
    r: list[CheckResult] = field(default_factory=list)
    services: list[CheckResult] = field(default_factory=list)

    def counts(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for result in [*self.python, *self.cli, *self.r, *self.services]:
            counts[result.status] = counts.get(result.status, 0) + 1
        return counts


def _short_output(result: subprocess.CompletedProcess[str]) -> str:
    text = "\n".join(part for part in (result.stdout.strip(), result.stderr.strip()) if part)
    return "\n".join(text.splitlines()[:8]).strip()


def run_command(command: list[str], timeout: int = 20, ok_codes: Iterable[int] = (0,)) -> tuple[str, str]:
    try:
        result = subprocess.run(command, capture_output=True, text=True, timeout=timeout)
    except FileNotFoundError as exc:
        return "missing", str(exc)
    except subprocess.TimeoutExpired as exc:
        return "failed", f"Timed out after {exc.timeout}s"
    except OSError as exc:
        return "failed", str(exc)

    output = _short_output(result)
    if result.returncode in set(ok_codes):
        return "ok", output or "command completed"
    return "failed", output or f"return code {result.returncode}"


def check_python_imports() -> list[CheckResult]:
    results: list[CheckResult] = []
    for package_name, import_name in PYTHON_IMPORTS.items():
        spec = importlib_util.find_spec(import_name)
        if spec is None:
            results.append(CheckResult(package_name, "missing", f"cannot import {import_name}"))
        else:
            origin = spec.origin if spec.origin else "namespace package"
            results.append(CheckResult(package_name, "ok", f"imports as {import_name}", path=origin))
    return results


def check_cli_tools() -> list[CheckResult]:
    results: list[CheckResult] = []
    for tool, commands in CLI_CHECKS.items():
        executable = shutil.which(commands[0][0])
        command_candidates = commands
        if executable is None:
            for candidate in CLI_EXECUTABLE_CANDIDATES.get(tool, []):
                if candidate.exists() and os.access(candidate, os.X_OK):
                    executable = str(candidate)
                    command_candidates = [[str(candidate), *command[1:]] for command in commands]
                    break
        if executable is None:
            results.append(CheckResult(tool, "missing", "not found on PATH"))
            continue

        last_status = "failed"
        last_detail = ""
        for command in command_candidates:
            ok_codes = (0, 1) if tool == "bwa" else (0,)
            status, detail = run_command(command, ok_codes=ok_codes)
            last_status, last_detail = status, detail
            if status == "ok":
                break
        results.append(CheckResult(tool, last_status, last_detail, path=executable))

    if shutil.which("docker"):
        status, detail = run_command(["docker", "info"], timeout=10)
        results.append(CheckResult("docker-daemon", status, detail))
    return results


def check_r_packages() -> list[CheckResult]:
    if shutil.which("Rscript") is None:
        return [CheckResult(pkg, "missing", "Rscript not found on PATH") for pkg in R_PACKAGES]

    expression = """
pkgs <- strsplit(Sys.getenv("METAINFORMANT_R_PACKAGES"), ",", fixed = TRUE)[[1]]
for (pkg in pkgs) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(pkg, "\t", as.character(packageVersion(pkg)), "\n", sep = "")
  } else {
    cat(pkg, "\tMISSING\n", sep = "")
  }
}
"""
    env = os.environ.copy()
    env["METAINFORMANT_R_PACKAGES"] = ",".join(R_PACKAGES)
    try:
        completed = subprocess.run(
            ["Rscript", "-e", expression],
            capture_output=True,
            text=True,
            timeout=120,
            env=env,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return [CheckResult(pkg, "failed", str(exc)) for pkg in R_PACKAGES]

    parsed: dict[str, str] = {}
    for line in completed.stdout.splitlines():
        if "\t" not in line:
            continue
        name, value = line.split("\t", 1)
        parsed[name.strip()] = value.strip()

    results: list[CheckResult] = []
    for pkg in R_PACKAGES:
        value = parsed.get(pkg)
        if value and value != "MISSING":
            results.append(CheckResult(pkg, "ok", "R package available", version=value))
        elif value == "MISSING":
            results.append(CheckResult(pkg, "missing", "R package not installed"))
        else:
            detail = completed.stderr.strip() or "R package status not reported"
            results.append(CheckResult(pkg, "failed", detail))
    return results


def wait_for_port(port: int, timeout: float = 15.0) -> bool:
    deadline = time.time() + timeout
    while time.time() < deadline:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.settimeout(0.5)
            if sock.connect_ex(("127.0.0.1", port)) == 0:
                return True
        time.sleep(0.25)
    return False


def verify_postgres(service_dir: Path) -> CheckResult:
    missing = [tool for tool in ("initdb", "postgres", "psql") if shutil.which(tool) is None]
    if missing:
        return CheckResult("postgres-temporary-service", "missing", f"missing tools: {', '.join(missing)}")

    port = 55432
    data_dir = service_dir / "postgres-data"
    socket_dir = service_dir / "postgres-socket"
    log_path = service_dir / "postgres.log"
    data_dir.mkdir(parents=True, exist_ok=True)
    socket_dir.mkdir(parents=True, exist_ok=True)

    if not (data_dir / "PG_VERSION").exists():
        status, detail = run_command(["initdb", "-D", str(data_dir), "-A", "trust", "-U", "metainformant"], timeout=60)
        if status != "ok":
            return CheckResult("postgres-temporary-service", status, f"initdb failed: {detail}")

    with log_path.open("a") as log:
        proc = subprocess.Popen(
            ["postgres", "-D", str(data_dir), "-p", str(port), "-k", str(socket_dir)],
            stdout=log,
            stderr=log,
            text=True,
        )
    try:
        if not wait_for_port(port):
            return CheckResult("postgres-temporary-service", "failed", f"port {port} did not open; see {log_path}")
        status, detail = run_command(
            ["psql", "-h", "127.0.0.1", "-p", str(port), "-U", "metainformant", "-d", "postgres", "-c", "SELECT 1"],
            timeout=15,
        )
        return CheckResult("postgres-temporary-service", status, detail)
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait(timeout=10)


def verify_redis(service_dir: Path) -> CheckResult:
    if shutil.which("redis-server") is None:
        return CheckResult("redis-temporary-service", "missing", "redis-server not found on PATH")

    port = 56379
    log_path = service_dir / "redis.log"
    with log_path.open("a") as log:
        proc = subprocess.Popen(
            [
                "redis-server",
                "--port",
                str(port),
                "--bind",
                "127.0.0.1",
                "--save",
                "",
                "--appendonly",
                "no",
                "--dir",
                str(service_dir),
            ],
            stdout=log,
            stderr=log,
            text=True,
        )
    try:
        if not wait_for_port(port):
            return CheckResult("redis-temporary-service", "failed", f"port {port} did not open; see {log_path}")
        with socket.create_connection(("127.0.0.1", port), timeout=5) as sock:
            sock.sendall(b"*1\r\n$4\r\nPING\r\n")
            reply = sock.recv(64).decode("utf-8", errors="replace")
        status = "ok" if "PONG" in reply else "failed"
        return CheckResult("redis-temporary-service", status, reply.strip())
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait(timeout=10)


def verify_mongodb(service_dir: Path) -> CheckResult:
    if shutil.which("mongod") is None:
        return CheckResult("mongodb-temporary-service", "missing", "mongod not found on PATH")

    port = 57017
    data_dir = service_dir / "mongodb-data"
    data_dir.mkdir(parents=True, exist_ok=True)
    log_path = service_dir / "mongodb.log"
    with log_path.open("a") as log:
        proc = subprocess.Popen(
            ["mongod", "--dbpath", str(data_dir), "--port", str(port), "--bind_ip", "127.0.0.1", "--quiet"],
            stdout=log,
            stderr=log,
            text=True,
        )
    try:
        if not wait_for_port(port):
            return CheckResult("mongodb-temporary-service", "failed", f"port {port} did not open; see {log_path}")
        if shutil.which("mongosh") is not None:
            status, detail = run_command(
                ["mongosh", "--quiet", f"mongodb://127.0.0.1:{port}/admin", "--eval", "db.runCommand({ping: 1}).ok"],
                timeout=20,
            )
        else:
            status, detail = run_command(
                [
                    sys.executable,
                    "-c",
                    (
                        "from pymongo import MongoClient; "
                        f"print(MongoClient('mongodb://127.0.0.1:{port}', serverSelectionTimeoutMS=5000).admin.command('ping')['ok'])"
                    ),
                ],
                timeout=20,
            )
        return CheckResult("mongodb-temporary-service", status, detail)
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait(timeout=10)


def check_services(output_dir: Path) -> list[CheckResult]:
    service_dir = (output_dir / "services").resolve()
    service_dir.mkdir(parents=True, exist_ok=True)
    return [verify_postgres(service_dir), verify_redis(service_dir), verify_mongodb(service_dir)]


def result_to_dict(result: CheckResult) -> dict[str, str | None]:
    return asdict(result)


def markdown_cell(value: str | None) -> str:
    if not value:
        return ""
    return value.replace("|", "\\|").replace("\n", "<br>")


def write_reports(report: VerificationReport, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "verification.json"
    md_path = output_dir / "verification.md"

    payload = {
        "generated_at": report.generated_at,
        "counts": report.counts(),
        "python": [result_to_dict(result) for result in report.python],
        "cli": [result_to_dict(result) for result in report.cli],
        "r": [result_to_dict(result) for result in report.r],
        "services": [result_to_dict(result) for result in report.services],
    }
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    lines = [
        "# MetaInformAnt Dependency Verification",
        "",
        f"- Generated: `{report.generated_at}`",
        f"- Counts: `{json.dumps(report.counts(), sort_keys=True)}`",
        "",
    ]
    for title, results in (
        ("Python Imports", report.python),
        ("CLI Tools", report.cli),
        ("R Packages", report.r),
        ("Temporary Services", report.services),
    ):
        lines.extend([f"## {title}", "", "| Name | Status | Detail | Version | Path |", "|---|---|---|---|---|"])
        for result in results:
            lines.append(
                f"| `{markdown_cell(result.name)}` | `{result.status}` | {markdown_cell(result.detail)} | `{markdown_cell(result.version)}` | `{markdown_cell(result.path)}` |"
            )
        lines.append("")
    md_path.write_text("\n".join(lines))


def build_report(output_dir: Path, include_services: bool) -> VerificationReport:
    return VerificationReport(
        generated_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        python=check_python_imports(),
        cli=check_cli_tools(),
        r=check_r_packages(),
        services=check_services(output_dir) if include_services else [],
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--skip-services", action="store_true", help="Do not start temporary service checks")
    parser.add_argument("--strict", action="store_true", help="Return non-zero if any check is not ok")
    args = parser.parse_args()

    report = build_report(args.output_dir, include_services=not args.skip_services)
    write_reports(report, args.output_dir)
    print(f"Wrote dependency verification reports to {args.output_dir}")
    print(f"Counts: {json.dumps(report.counts(), sort_keys=True)}")

    if args.strict and any(status != "ok" for status in report.counts()):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
