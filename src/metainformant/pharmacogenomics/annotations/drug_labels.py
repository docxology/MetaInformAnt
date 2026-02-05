"""FDA drug label parsing for pharmacogenomic biomarker information.

Parses FDA-approved drug labels to extract pharmacogenomic testing requirements,
biomarker information, and genotype-specific recommendations. Classifies labels
by the strength of PGx testing language (required, recommended, actionable, informational).

The FDA categorizes pharmacogenomic biomarkers in drug labels as:
    - Required: Genetic testing mandated before prescribing
    - Recommended: Testing strongly suggested
    - Actionable: Genotype information affects dosing or drug selection
    - Informational: PGx information provided for awareness
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


# ── Built-in FDA drug label data ──────────────────────────────────────────────
# Key drugs with pharmacogenomic labeling information.

_FDA_DRUG_LABELS: list[dict[str, Any]] = [
    {
        "drug": "abacavir",
        "brand_name": "Ziagen",
        "gene_biomarker": "HLA-B",
        "variant": "HLA-B*57:01",
        "label_type": "required",
        "boxed_warning": True,
        "section": "Warnings and Precautions",
        "label_text": "Serious and sometimes fatal hypersensitivity reactions have been associated with "
        "HLA-B*5701 allele. Screen for HLA-B*5701 allele before starting abacavir.",
        "recommendation": "Do not use abacavir in patients positive for HLA-B*5701.",
        "therapeutic_area": "HIV",
        "fda_table_url": "https://www.fda.gov/drugs/science-and-research/table-pharmacogenomic-biomarkers-drug-labeling",
    },
    {
        "drug": "carbamazepine",
        "brand_name": "Tegretol",
        "gene_biomarker": "HLA-B",
        "variant": "HLA-B*15:02",
        "label_type": "required",
        "boxed_warning": True,
        "section": "Warnings and Precautions",
        "label_text": "Serious dermatologic reactions including Stevens-Johnson syndrome and toxic epidermal "
        "necrolysis (SJS/TEN). HLA-B*15:02 allele increases risk. Patients of Southeast Asian "
        "ancestry should be screened before starting carbamazepine.",
        "recommendation": "Screen patients with ancestry in populations with high prevalence of HLA-B*15:02. "
        "Do not use carbamazepine in patients positive for HLA-B*15:02 unless benefit clearly outweighs risk.",
        "therapeutic_area": "Epilepsy",
        "fda_table_url": "https://www.fda.gov/drugs/science-and-research/table-pharmacogenomic-biomarkers-drug-labeling",
    },
    {
        "drug": "clopidogrel",
        "brand_name": "Plavix",
        "gene_biomarker": "CYP2C19",
        "variant": "CYP2C19 poor metabolizer",
        "label_type": "actionable",
        "boxed_warning": True,
        "section": "Boxed Warning",
        "label_text": "Diminished effectiveness in CYP2C19 poor metabolizers. Tests are available to identify "
        "CYP2C19 genotype. Consider use of another platelet P2Y12 inhibitor in poor metabolizers.",
        "recommendation": "Consider alternative antiplatelet therapy (prasugrel or ticagrelor) in CYP2C19 poor metabolizers.",
        "therapeutic_area": "Cardiology",
        "fda_table_url": "https://www.fda.gov/drugs/science-and-research/table-pharmacogenomic-biomarkers-drug-labeling",
    },
    {
        "drug": "codeine",
        "brand_name": "Various",
        "gene_biomarker": "CYP2D6",
        "variant": "CYP2D6 ultrarapid metabolizer",
        "label_type": "actionable",
        "boxed_warning": True,
        "section": "Boxed Warning",
        "label_text": "Respiratory depression and death have occurred in children who received codeine following "
        "tonsillectomy/adenoidectomy and had evidence of being CYP2D6 ultrarapid metabolizers.",
        "recommendation": "Avoid use of codeine in CYP2D6 ultrarapid metabolizers. Use alternative analgesics.",
        "therapeutic_area": "Pain Management",
        "fda_table_url": "https://www.fda.gov/drugs/science-and-research/table-pharmacogenomic-biomarkers-drug-labeling",
    },
    {
        "drug": "warfarin",
        "brand_name": "Coumadin",
        "gene_biomarker": "CYP2C9, VKORC1",
        "variant": "CYP2C9 *2, *3; VKORC1 -1639G>A",
        "label_type": "actionable",
        "boxed_warning": False,
        "section": "Clinical Pharmacology",
        "label_text": "Genetic variations in CYP2C9 and VKORC1 genes affect warfarin pharmacokinetics "
        "and pharmacodynamics. Consider lower initial and maintenance doses in patients "
        "known to have CYP2C9 *2 and *3 alleles or VKORC1 A allele.",
        "recommendation": "Consider pharmacogenomic-guided dosing for patients with CYP2C9 and/or VKORC1 variants.",
        "therapeutic_area": "Cardiology",
        "fda_table_url": "https://www.fda.gov/drugs/science-and-research/table-pharmacogenomic-biomarkers-drug-labeling",
    },
    {
        "drug": "simvastatin",
        "brand_name": "Zocor",
        "gene_biomarker": "SLCO1B1",
        "variant": "SLCO1B1 *5 (rs4149056)",
        "label_type": "recommended",
        "boxed_warning": False,
        "section": "Warnings and Precautions",
        "label_text": "Risk of myopathy/rhabdomyolysis is dose-related and increased by SLCO1B1 *5 allele.",
        "recommendation": "Consider SLCO1B1 genotyping before initiating simvastatin 40mg or 80mg doses.",
        "therapeutic_area": "Cardiology",
        "fda_table_url": "https://www.fda.gov/drugs/science-and-research/table-pharmacogenomic-biomarkers-drug-labeling",
    },
    {
        "drug": "fluorouracil",
        "brand_name": "Adrucil",
        "gene_biomarker": "DPYD",
        "variant": "DPYD *2A, *13, D949V, HapB3",
        "label_type": "recommended",
        "boxed_warning": False,
        "section": "Warnings and Precautions",
        "label_text": "Patients with certain DPYD variants have increased risk of severe or fatal "
        "fluoropyrimidine toxicity. Consider DPYD genotyping or phenotyping before "
        "initiating fluoropyrimidine-based chemotherapy.",
        "recommendation": "Test for DPYD deficiency before starting fluoropyrimidine therapy. "
        "Avoid in patients with complete DPD deficiency. Reduce dose in partial deficiency.",
        "therapeutic_area": "Oncology",
        "fda_table_url": "https://www.fda.gov/drugs/science-and-research/table-pharmacogenomic-biomarkers-drug-labeling",
    },
    {
        "drug": "azathioprine",
        "brand_name": "Imuran",
        "gene_biomarker": "TPMT, NUDT15",
        "variant": "TPMT and NUDT15 poor metabolizer alleles",
        "label_type": "recommended",
        "boxed_warning": False,
        "section": "Warnings and Precautions",
        "label_text": "Consider genotyping or phenotyping for TPMT and NUDT15 before starting azathioprine. "
        "Patients with TPMT or NUDT15 deficiency are at increased risk for severe myelotoxicity.",
        "recommendation": "Consider TPMT and NUDT15 testing before starting thiopurine therapy. "
        "Reduce dose or use alternative in poor metabolizers.",
        "therapeutic_area": "Immunology",
        "fda_table_url": "https://www.fda.gov/drugs/science-and-research/table-pharmacogenomic-biomarkers-drug-labeling",
    },
    {
        "drug": "tamoxifen",
        "brand_name": "Nolvadex",
        "gene_biomarker": "CYP2D6",
        "variant": "CYP2D6 poor metabolizer",
        "label_type": "informational",
        "boxed_warning": False,
        "section": "Clinical Pharmacology",
        "label_text": "CYP2D6 is required for conversion of tamoxifen to the active metabolite endoxifen. "
        "Reduced CYP2D6 activity may result in lower endoxifen concentrations.",
        "recommendation": "CYP2D6 poor metabolizers may have reduced benefit from tamoxifen. "
        "Consider aromatase inhibitors as alternative in postmenopausal patients.",
        "therapeutic_area": "Oncology",
        "fda_table_url": "https://www.fda.gov/drugs/science-and-research/table-pharmacogenomic-biomarkers-drug-labeling",
    },
    {
        "drug": "ivacaftor",
        "brand_name": "Kalydeco",
        "gene_biomarker": "CFTR",
        "variant": "Multiple CFTR mutations",
        "label_type": "required",
        "boxed_warning": False,
        "section": "Indications and Usage",
        "label_text": "Indicated for treatment of cystic fibrosis in patients who have specific CFTR mutations "
        "that are responsive to ivacaftor. CFTR genotyping is required to determine eligibility.",
        "recommendation": "Confirm presence of a responsive CFTR mutation using FDA-cleared CFTR genotyping test.",
        "therapeutic_area": "Pulmonology",
        "fda_table_url": "https://www.fda.gov/drugs/science-and-research/table-pharmacogenomic-biomarkers-drug-labeling",
    },
]


def parse_drug_label(
    filepath: str | Path,
) -> dict[str, Any]:
    """Parse an FDA drug label file for pharmacogenomic information.

    Reads a drug label file (JSON or structured text) and extracts PGx-relevant
    sections including biomarker requirements, dosing modifications, and warnings.

    Args:
        filepath: Path to drug label file (JSON format preferred)

    Returns:
        Dictionary with parsed label data including:
            - "drug": Drug name
            - "brand_name": Brand name
            - "gene_biomarker": Gene/biomarker mentioned
            - "sections": Dict of label sections with PGx content
            - "biomarker_info": Extracted biomarker testing information
            - "label_type": Classification of PGx labeling strength
    """
    from metainformant.core import io

    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"Drug label file not found: {filepath}")

    if path.suffix in (".json", ".gz"):
        data = io.load_json(path)
    else:
        # Read as structured text
        text = path.read_text(encoding="utf-8")
        data = _parse_label_text(text)

    # Extract PGx-relevant fields
    result: dict[str, Any] = {
        "drug": data.get(
            "drug",
            data.get("openfda", {}).get("generic_name", [""])[0] if isinstance(data.get("openfda"), dict) else "",
        ),
        "brand_name": data.get(
            "brand_name",
            data.get("openfda", {}).get("brand_name", [""])[0] if isinstance(data.get("openfda"), dict) else "",
        ),
        "gene_biomarker": data.get("gene_biomarker", ""),
        "sections": {},
        "biomarker_info": {},
        "label_type": data.get("label_type", "informational"),
    }

    # Extract PGx-relevant sections
    pgx_section_keys = [
        "boxed_warning",
        "warnings_and_precautions",
        "indications_and_usage",
        "dosage_and_administration",
        "clinical_pharmacology",
        "warnings",
        "precautions",
        "contraindications",
    ]

    for key in pgx_section_keys:
        content = data.get(key, "")
        if isinstance(content, list):
            content = " ".join(str(c) for c in content)
        if content and _contains_pgx_terms(str(content)):
            result["sections"][key] = str(content)

    result["biomarker_info"] = extract_biomarker_info(result)
    result["label_type"] = classify_label_type(result)["label_type"]

    logger.info(
        "Parsed drug label for %s (biomarker: %s, type: %s)",
        result["drug"],
        result["gene_biomarker"],
        result["label_type"],
    )

    return result


def _parse_label_text(text: str) -> dict[str, Any]:
    """Parse unstructured drug label text into sections.

    Args:
        text: Raw drug label text

    Returns:
        Dictionary with identified sections
    """
    sections: dict[str, str] = {}
    current_section = "general"
    current_text: list[str] = []

    section_headers = [
        "BOXED WARNING",
        "INDICATIONS AND USAGE",
        "DOSAGE AND ADMINISTRATION",
        "WARNINGS AND PRECAUTIONS",
        "CONTRAINDICATIONS",
        "CLINICAL PHARMACOLOGY",
        "WARNINGS",
        "PRECAUTIONS",
    ]

    for line in text.split("\n"):
        stripped = line.strip().upper()
        matched_header = None
        for header in section_headers:
            if header in stripped:
                matched_header = header
                break

        if matched_header:
            if current_text:
                sections[current_section] = "\n".join(current_text)
            current_section = matched_header.lower().replace(" ", "_")
            current_text = []
        else:
            current_text.append(line)

    if current_text:
        sections[current_section] = "\n".join(current_text)

    return sections


def _contains_pgx_terms(text: str) -> bool:
    """Check if text contains pharmacogenomic-relevant terms.

    Args:
        text: Text to check

    Returns:
        True if PGx terms found
    """
    pgx_terms = [
        "pharmacogenomic",
        "pharmacogenetic",
        "genotype",
        "genotyping",
        "allele",
        "metabolizer",
        "CYP",
        "HLA",
        "DPYD",
        "TPMT",
        "NUDT15",
        "SLCO1B1",
        "VKORC1",
        "CFTR",
        "G6PD",
        "UGT1A1",
        "poor metabolizer",
        "intermediate metabolizer",
        "ultrarapid metabolizer",
        "gene",
        "polymorphism",
        "variant",
        "biomarker",
    ]
    text_lower = text.lower()
    return any(term.lower() in text_lower for term in pgx_terms)


def extract_biomarker_info(
    label_data: dict[str, Any],
) -> dict[str, Any]:
    """Extract biomarker testing requirements from parsed label data.

    Args:
        label_data: Parsed drug label dictionary

    Returns:
        Dictionary with:
            - "biomarker_gene": Gene/biomarker name
            - "testing_required": Whether testing is mandated
            - "testing_recommended": Whether testing is recommended
            - "boxed_warning": Whether PGx info is in boxed warning
            - "affected_populations": Populations with specific considerations
            - "test_timing": When testing should be performed
    """
    gene = label_data.get("gene_biomarker", "")
    sections = label_data.get("sections", {})
    label_type = label_data.get("label_type", "informational")

    has_boxed = bool(sections.get("boxed_warning", ""))

    # Check for testing language
    all_text = " ".join(str(v) for v in sections.values()).lower()
    testing_required = any(
        phrase in all_text
        for phrase in ["must be tested", "testing is required", "genotyping is required", "screen for", "test for"]
    )
    testing_recommended = any(
        phrase in all_text
        for phrase in ["consider testing", "testing is recommended", "consider genotyping", "should be tested"]
    )

    # Population detection
    populations: list[str] = []
    pop_terms = {
        "southeast asian": "Southeast Asian",
        "east asian": "East Asian",
        "african": "African",
        "european": "European",
        "pediatric": "Pediatric",
    }
    for term, label in pop_terms.items():
        if term in all_text:
            populations.append(label)

    return {
        "biomarker_gene": gene,
        "testing_required": testing_required or label_type == "required",
        "testing_recommended": testing_recommended or label_type == "recommended",
        "boxed_warning": has_boxed,
        "affected_populations": populations,
        "test_timing": "Before initiation of therapy" if (testing_required or testing_recommended) else "Not specified",
    }


def classify_label_type(
    label_data: dict[str, Any],
) -> dict[str, Any]:
    """Classify the strength of PGx labeling for a drug.

    Categorizes the drug label based on the strength of pharmacogenomic
    testing language:
        - "required": Testing mandated before prescribing
        - "recommended": Testing strongly suggested
        - "actionable": Genotype information affects clinical decisions
        - "informational": PGx information provided for awareness

    Args:
        label_data: Parsed drug label dictionary

    Returns:
        Dictionary with:
            - "label_type": Classification string
            - "rationale": Explanation for classification
    """
    # Check if already classified
    explicit_type = label_data.get("label_type", "")
    if explicit_type in ("required", "recommended", "actionable", "informational"):
        return {"label_type": explicit_type, "rationale": "Explicitly classified in label data"}

    sections = label_data.get("sections", {})
    all_text = " ".join(str(v) for v in sections.values()).lower()

    # Required: testing mandated language
    required_phrases = [
        "must be tested",
        "testing is required",
        "genotyping is required",
        "should not be used without",
        "contraindicated without testing",
    ]
    if any(phrase in all_text for phrase in required_phrases):
        return {
            "label_type": "required",
            "rationale": "Label contains mandatory testing language",
        }

    # Recommended: strong suggestion
    recommended_phrases = [
        "testing is recommended",
        "consider genotyping",
        "should be tested",
        "genotyping should be performed",
        "testing should be considered",
    ]
    if any(phrase in all_text for phrase in recommended_phrases):
        return {
            "label_type": "recommended",
            "rationale": "Label contains recommended testing language",
        }

    # Actionable: dosing modifications based on genotype
    actionable_phrases = [
        "dose adjustment",
        "dose reduction",
        "alternative therapy",
        "consider alternative",
        "poor metabolizer",
        "ultrarapid metabolizer",
        "reduced dose",
    ]
    if any(phrase in all_text for phrase in actionable_phrases):
        return {
            "label_type": "actionable",
            "rationale": "Label contains genotype-based clinical action language",
        }

    return {
        "label_type": "informational",
        "rationale": "Label provides PGx information without specific clinical action directives",
    }


def search_labels_by_gene(
    gene: str,
    labels_db: list[dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    """Find drug labels mentioning a specific gene or biomarker.

    Args:
        gene: Gene symbol (case-insensitive)
        labels_db: Optional pre-loaded labels database. If None, uses built-in data.

    Returns:
        List of matching drug label dictionaries
    """
    if labels_db is None:
        labels_db = _FDA_DRUG_LABELS

    gene_upper = gene.upper().strip()
    results: list[dict[str, Any]] = []

    for label in labels_db:
        biomarker = label.get("gene_biomarker", "").upper()
        # Check if gene appears in the biomarker field (may contain multiple genes)
        if gene_upper in biomarker or gene_upper in biomarker.replace(",", " ").split():
            results.append(label)

    logger.info("Found %d drug labels mentioning %s", len(results), gene)
    return results
