"""Randomness audit for the pharmacogenomics package.

Audit result: nothing in the package draws random numbers. Allele calling, diplotype scoring,
phenotype classification, drug-interaction analysis, reporting and plotting are all deterministic
functions of their inputs (the one non-input value is the report's ``generated_at`` timestamp).
These tests keep it that way: a lexical scan for random-number use, and a check that results do
not depend on string-hash ordering of the sets and dicts the analysis code iterates, which Python
randomises per process unless ``PYTHONHASHSEED`` is fixed.

REAL IMPLEMENTATION -- the hash-seed check runs the real analysis in fresh interpreters.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from pathlib import Path

import metainformant.pharmacogenomics as pharmacogenomics

PACKAGE_DIR = Path(pharmacogenomics.__file__).resolve().parent

# Any use of a random-number or random-id source. If a stochastic method is ever added it must
# take an explicit seed, and this pattern list is the place to allow it.
RANDOMNESS = re.compile(
    r"^\s*(?:import|from)\s+(?:random|secrets)\b"
    r"|\b(?:np|numpy)\.random\b"
    r"|\b(?:default_rng|RandomState|SeedSequence)\b"
    r"|\bos\.urandom\b"
    r"|\buuid\.uuid[14]\b",
    re.MULTILINE,
)

# Exercises every analysis path with inputs that make set/dict ordering matter if it can:
# several overlapping variants, ties in the number of defining variants, many genes and drugs.
ANALYSIS = r"""
import json
import logging

logging.disable(logging.CRITICAL)

from metainformant.pharmacogenomics.alleles.diplotype import (
    determine_diplotype,
    phased_diplotype,
    resolve_ambiguous_diplotypes,
)
from metainformant.pharmacogenomics.alleles.phenotype import classify_phenotype, population_phenotype_frequencies
from metainformant.pharmacogenomics.alleles.star_allele import (
    call_star_alleles,
    detect_novel_alleles,
    handle_cyp2d6_cnv,
    load_allele_definitions,
)
from metainformant.pharmacogenomics.clinical.drug_interaction import polypharmacy_analysis
from metainformant.pharmacogenomics.clinical.reporting import export_report, generate_clinical_report
from metainformant.pharmacogenomics.metabolism.metabolizer_status import predict_metabolizer_status

out = {}
variants = {"rs3892097", "rs16947", "rs1065852", "rs5030655", "CYP2D6_DUP", "rs28371706", "rs35742686", "rs9999999"}
out["call"] = [a.name for a in call_star_alleles(variants, "CYP2D6")]
out["tpmt"] = [a.name for a in call_star_alleles({"rs1800460", "rs1142345", "rs1800462"}, "TPMT")]
out["novel"] = detect_novel_alleles(variants, load_allele_definitions("CYP2D6"))
out["cnv"] = {k: v for k, v in handle_cyp2d6_cnv(variants, 3).items() if k != "alleles"}
phased = phased_diplotype(variants, {"rs3892097": 0, "rs16947": 1, "rs1065852": 1}, "CYP2D6")
out["phased"] = [phased.diplotype_string, phased.activity_score, phased.confidence]
out["resolved"] = resolve_ambiguous_diplotypes(
    [determine_diplotype("*1", "*4", "CYP2D6"), determine_diplotype("*2", "*3", "CYP2D6")]
).diplotype_string
out["classify"] = {
    f"{gene} {dip}": classify_phenotype(dip, gene)["phenotype"].value
    for gene, dip in [("CYP2D6", "*1/*4"), ("CYP2C19", "*9/*17"), ("SLCO1B1", "*1/*5"), ("DPYD", "*1/HapB3")]
}
allele_frequencies = {"*1": 0.5, "*2": 0.1, "*4": 0.2, "*10": 0.1, "*41": 0.1}
out["population"] = population_phenotype_frequencies(allele_frequencies, "CYP2D6")
out["predict"] = predict_metabolizer_status({"diplotype": "*2/*41"}, "CYP2D6")
genotypes = {
    "CYP2D6": {"diplotype": "*4/*10"},
    "CYP2C19": {"diplotype": "*2/*17"},
    "CYP2C9": {"diplotype": "*1/*2"},
    "DPYD": {"diplotype": "*1/c.2846A>T"},
    "TPMT": {"diplotype": "*1/*3A"},
    "SLCO1B1": {"diplotype": "*1/*5"},
}
drugs = ["codeine", "clopidogrel", "warfarin", "fluorouracil", "azathioprine", "simvastatin", "tamoxifen", "tramadol"]
phenotypes = {g: {"phenotype": classify_phenotype(v["diplotype"], g)["phenotype"]} for g, v in genotypes.items()}
poly = polypharmacy_analysis(drugs, phenotypes)
out["poly"] = {k: v for k, v in poly.items() if k != "interactions"}
report = generate_clinical_report({"patient_id": "P"}, genotypes)
report["header"].pop("generated_at")
out["report"] = json.loads(export_report(report, format="json"))
print(json.dumps(out, default=str))
"""


def _run_analysis(hash_seed: str) -> str:
    env = {**os.environ, "PYTHONHASHSEED": hash_seed}
    completed = subprocess.run(
        [sys.executable, "-c", ANALYSIS],
        env=env,
        capture_output=True,
        text=True,
        check=False,
        timeout=300,
    )
    assert completed.returncode == 0, completed.stderr[-2000:]
    return completed.stdout.strip().splitlines()[-1]


class TestNoUnseededRandomness:
    def test_package_sources_draw_no_random_numbers(self) -> None:
        sources = sorted(PACKAGE_DIR.rglob("*.py"))
        assert len(sources) >= 15  # the scan really covers the package
        offenders: dict[str, list[str]] = {}
        for path in sources:
            found = sorted({m.group(0).strip() for m in RANDOMNESS.finditer(path.read_text(encoding="utf-8"))})
            if found:
                offenders[str(path.relative_to(PACKAGE_DIR))] = found
        assert offenders == {}

    def test_scan_pattern_detects_the_calls_it_exists_to_forbid(self) -> None:
        for snippet in [
            "import random",
            "from random import choice",
            "x = np.random.rand(3)",
            "numpy.random.seed(0)",
            "rng = default_rng()",
            "state = RandomState()",
            "os.urandom(8)",
            "token = uuid.uuid4()",
        ]:
            assert RANDOMNESS.search(snippet), snippet
        assert not RANDOMNESS.search("operand = 1  # a randomised design is not used here")

    def test_results_do_not_depend_on_string_hash_ordering(self) -> None:
        outputs = {seed: _run_analysis(seed) for seed in ("0", "1", "4242")}
        reference = outputs["0"]
        assert json.loads(reference)["report"]["genotype_results"]  # a real analysis ran
        for seed, output in outputs.items():
            assert output == reference, f"PYTHONHASHSEED={seed} changed an analysis result"

    def test_repeated_calls_in_one_process_agree(self) -> None:
        from metainformant.pharmacogenomics.alleles.phenotype import population_phenotype_frequencies
        from metainformant.pharmacogenomics.alleles.star_allele import call_star_alleles

        variants = {"rs3892097", "rs16947", "rs1065852", "rs28371706"}
        first = [a.name for a in call_star_alleles(variants, "CYP2D6")]
        frequencies = {"*1": 0.6, "*4": 0.25, "*10": 0.15}
        first_freq = population_phenotype_frequencies(frequencies, "CYP2D6")
        for _ in range(5):
            assert [a.name for a in call_star_alleles(variants, "CYP2D6")] == first
            assert population_phenotype_frequencies(frequencies, "CYP2D6") == first_freq
