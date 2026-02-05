"""PharmGKB (Pharmacogenomics Knowledge Base) annotation support.

Provides querying and parsing of PharmGKB clinical annotations, variant
annotations, drug pathway data, and evidence level classification.

PharmGKB evidence levels:
    1A: Variant-drug combination in a CPIC or DPWG guideline, or FDA PGx biomarker
    1B: Variant-drug combination with strong evidence, not yet in a guideline
    2A: Variant-drug with moderate evidence (from PharmGKB clinical annotation)
    2B: Variant-drug with moderate evidence (from PharmGKB variant annotation)
    3: Variant-drug with low-level evidence (case reports, small studies)
    4: Variant-drug with preliminary evidence (in vitro, preclinical)
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


# ── Built-in PharmGKB annotation data ──────────────────────────────────────────
# Curated clinical annotations for key pharmacogenomic variants.

_PHARMGKB_CLINICAL_ANNOTATIONS: list[dict[str, Any]] = [
    {
        "annotation_id": "PA166104949",
        "gene": "CYP2D6",
        "drug": "codeine",
        "variant": None,
        "evidence_level": "1A",
        "phenotype_category": "Efficacy",
        "annotation_text": "CYP2D6 poor metabolizers have significantly reduced morphine formation "
        "from codeine. Ultrarapid metabolizers have increased morphine formation "
        "with risk of toxicity.",
        "population": "General",
    },
    {
        "annotation_id": "PA166104944",
        "gene": "CYP2C19",
        "drug": "clopidogrel",
        "variant": None,
        "evidence_level": "1A",
        "phenotype_category": "Efficacy",
        "annotation_text": "CYP2C19 poor and intermediate metabolizers have reduced conversion of "
        "clopidogrel to active metabolite, leading to decreased platelet inhibition.",
        "population": "General",
    },
    {
        "annotation_id": "PA166104938",
        "gene": "CYP2C19",
        "drug": "voriconazole",
        "variant": None,
        "evidence_level": "1A",
        "phenotype_category": "Toxicity",
        "annotation_text": "CYP2C19 poor metabolizers have higher voriconazole trough concentrations, "
        "increasing risk of adverse effects including hepatotoxicity and visual disturbances.",
        "population": "General",
    },
    {
        "annotation_id": "PA166153760",
        "gene": "DPYD",
        "drug": "fluorouracil",
        "variant": None,
        "evidence_level": "1A",
        "phenotype_category": "Toxicity",
        "annotation_text": "DPYD poor metabolizers have dramatically increased exposure to fluoropyrimidines "
        "with severe, potentially fatal toxicity including neutropenia and mucositis.",
        "population": "General",
    },
    {
        "annotation_id": "PA166104951",
        "gene": "TPMT",
        "drug": "azathioprine",
        "variant": None,
        "evidence_level": "1A",
        "phenotype_category": "Toxicity",
        "annotation_text": "TPMT poor metabolizers accumulate thioguanine nucleotides with high risk of "
        "severe myelosuppression. Intermediate metabolizers require dose reduction.",
        "population": "General",
    },
    {
        "annotation_id": "PA166105006",
        "gene": "SLCO1B1",
        "drug": "simvastatin",
        "variant": "rs4149056",
        "evidence_level": "1A",
        "phenotype_category": "Toxicity",
        "annotation_text": "SLCO1B1 rs4149056 C allele carriers have increased simvastatin acid exposure "
        "and higher risk of statin-induced myopathy.",
        "population": "General",
    },
    {
        "annotation_id": "PA166155954",
        "gene": "CYP2C9",
        "drug": "warfarin",
        "variant": None,
        "evidence_level": "1A",
        "phenotype_category": "Dosage",
        "annotation_text": "CYP2C9 poor and intermediate metabolizers require lower warfarin doses "
        "to achieve therapeutic INR. Standard doses increase bleeding risk.",
        "population": "General",
    },
    {
        "annotation_id": "PA166176901",
        "gene": "CYP2D6",
        "drug": "tamoxifen",
        "variant": None,
        "evidence_level": "1A",
        "phenotype_category": "Efficacy",
        "annotation_text": "CYP2D6 poor metabolizers have reduced endoxifen concentrations and may "
        "have worse breast cancer outcomes with tamoxifen therapy.",
        "population": "General",
    },
    {
        "annotation_id": "PA166176950",
        "gene": "NUDT15",
        "drug": "azathioprine",
        "variant": "rs116855232",
        "evidence_level": "1A",
        "phenotype_category": "Toxicity",
        "annotation_text": "NUDT15 poor metabolizers (*3/*3) have extremely high risk of leukopenia "
        "with thiopurine therapy. Intermediate metabolizers require dose reduction.",
        "population": "East Asian",
    },
]

# ── Variant-specific annotations ──────────────────────────────────────────────

_VARIANT_ANNOTATIONS: dict[str, dict[str, Any]] = {
    "rs4244285": {
        "rsid": "rs4244285",
        "gene": "CYP2C19",
        "allele": "*2",
        "chromosome": "10",
        "position": 96541616,
        "ref": "G",
        "alt": "A",
        "consequence": "Splicing defect",
        "clinical_significance": "Loss of function",
        "evidence_level": "1A",
        "allele_frequency": {
            "Global": 0.28,
            "African": 0.17,
            "East_Asian": 0.30,
            "European": 0.15,
            "South_Asian": 0.36,
        },
        "associated_drugs": ["clopidogrel", "voriconazole", "escitalopram", "omeprazole"],
    },
    "rs4986893": {
        "rsid": "rs4986893",
        "gene": "CYP2C19",
        "allele": "*3",
        "chromosome": "10",
        "position": 96540410,
        "ref": "G",
        "alt": "A",
        "consequence": "Premature stop codon",
        "clinical_significance": "Loss of function",
        "evidence_level": "1A",
        "allele_frequency": {
            "Global": 0.01,
            "African": 0.002,
            "East_Asian": 0.07,
            "European": 0.001,
            "South_Asian": 0.01,
        },
        "associated_drugs": ["clopidogrel", "voriconazole"],
    },
    "rs12248560": {
        "rsid": "rs12248560",
        "gene": "CYP2C19",
        "allele": "*17",
        "chromosome": "10",
        "position": 96522463,
        "ref": "C",
        "alt": "T",
        "consequence": "Increased transcription",
        "clinical_significance": "Increased function",
        "evidence_level": "1A",
        "allele_frequency": {
            "Global": 0.21,
            "African": 0.27,
            "East_Asian": 0.01,
            "European": 0.21,
            "South_Asian": 0.15,
        },
        "associated_drugs": ["clopidogrel", "voriconazole", "escitalopram"],
    },
    "rs3892097": {
        "rsid": "rs3892097",
        "gene": "CYP2D6",
        "allele": "*4",
        "chromosome": "22",
        "position": 42128945,
        "ref": "C",
        "alt": "T",
        "consequence": "Splicing defect",
        "clinical_significance": "Loss of function",
        "evidence_level": "1A",
        "allele_frequency": {
            "Global": 0.12,
            "African": 0.07,
            "East_Asian": 0.01,
            "European": 0.19,
            "South_Asian": 0.10,
        },
        "associated_drugs": ["codeine", "tramadol", "tamoxifen", "ondansetron"],
    },
    "rs16947": {
        "rsid": "rs16947",
        "gene": "CYP2D6",
        "allele": "*2",
        "chromosome": "22",
        "position": 42130692,
        "ref": "G",
        "alt": "A",
        "consequence": "Missense (R296C)",
        "clinical_significance": "Normal function",
        "evidence_level": "1A",
        "allele_frequency": {
            "Global": 0.33,
            "African": 0.34,
            "East_Asian": 0.15,
            "European": 0.29,
            "South_Asian": 0.35,
        },
        "associated_drugs": ["codeine", "tramadol", "tamoxifen"],
    },
    "rs1065852": {
        "rsid": "rs1065852",
        "gene": "CYP2D6",
        "allele": "*10",
        "chromosome": "22",
        "position": 42129770,
        "ref": "C",
        "alt": "T",
        "consequence": "Missense (P34S)",
        "clinical_significance": "Decreased function",
        "evidence_level": "1A",
        "allele_frequency": {
            "Global": 0.22,
            "African": 0.08,
            "East_Asian": 0.50,
            "European": 0.05,
            "South_Asian": 0.15,
        },
        "associated_drugs": ["codeine", "tramadol", "tamoxifen"],
    },
    "rs4149056": {
        "rsid": "rs4149056",
        "gene": "SLCO1B1",
        "allele": "*5",
        "chromosome": "12",
        "position": 21331549,
        "ref": "T",
        "alt": "C",
        "consequence": "Missense (V174A)",
        "clinical_significance": "Decreased function",
        "evidence_level": "1A",
        "allele_frequency": {
            "Global": 0.08,
            "African": 0.02,
            "East_Asian": 0.12,
            "European": 0.15,
            "South_Asian": 0.05,
        },
        "associated_drugs": ["simvastatin", "atorvastatin", "methotrexate"],
    },
    "rs1799853": {
        "rsid": "rs1799853",
        "gene": "CYP2C9",
        "allele": "*2",
        "chromosome": "10",
        "position": 96702047,
        "ref": "C",
        "alt": "T",
        "consequence": "Missense (R144C)",
        "clinical_significance": "Decreased function",
        "evidence_level": "1A",
        "allele_frequency": {
            "Global": 0.06,
            "African": 0.02,
            "East_Asian": 0.00,
            "European": 0.13,
            "South_Asian": 0.04,
        },
        "associated_drugs": ["warfarin", "phenytoin"],
    },
    "rs1057910": {
        "rsid": "rs1057910",
        "gene": "CYP2C9",
        "allele": "*3",
        "chromosome": "10",
        "position": 96741053,
        "ref": "A",
        "alt": "C",
        "consequence": "Missense (I359L)",
        "clinical_significance": "Loss of function",
        "evidence_level": "1A",
        "allele_frequency": {
            "Global": 0.04,
            "African": 0.01,
            "East_Asian": 0.04,
            "European": 0.06,
            "South_Asian": 0.08,
        },
        "associated_drugs": ["warfarin", "phenytoin"],
    },
    "rs3918290": {
        "rsid": "rs3918290",
        "gene": "DPYD",
        "allele": "*2A",
        "chromosome": "1",
        "position": 97915614,
        "ref": "C",
        "alt": "T",
        "consequence": "Splicing defect (IVS14+1G>A)",
        "clinical_significance": "Loss of function",
        "evidence_level": "1A",
        "allele_frequency": {
            "Global": 0.01,
            "African": 0.001,
            "East_Asian": 0.001,
            "European": 0.01,
            "South_Asian": 0.001,
        },
        "associated_drugs": ["fluorouracil", "capecitabine"],
    },
}

# ── Drug pathway data ──────────────────────────────────────────────────────────

_DRUG_PATHWAYS: dict[str, dict[str, Any]] = {
    "clopidogrel": {
        "drug": "clopidogrel",
        "pharmgkb_id": "PA449053",
        "pathway_type": "Pharmacokinetics",
        "genes": ["CYP2C19", "CYP1A2", "CYP2B6", "CYP3A4", "CYP2C9", "PON1", "ABCB1"],
        "key_metabolizing_enzyme": "CYP2C19",
        "active_metabolite": "clopidogrel thiol active metabolite",
        "mechanism": "CYP2C19 is the primary enzyme for conversion of clopidogrel prodrug to active thiol metabolite "
        "that irreversibly inhibits the P2Y12 ADP receptor on platelets.",
        "transporter_genes": ["ABCB1"],
    },
    "codeine": {
        "drug": "codeine",
        "pharmgkb_id": "PA449088",
        "pathway_type": "Pharmacokinetics",
        "genes": ["CYP2D6", "CYP3A4", "UGT2B7", "UGT2B4"],
        "key_metabolizing_enzyme": "CYP2D6",
        "active_metabolite": "morphine",
        "mechanism": "CYP2D6 O-demethylates codeine to morphine, the active analgesic metabolite. "
        "CYP3A4 N-demethylates codeine to norcodeine (inactive).",
        "transporter_genes": [],
    },
    "warfarin": {
        "drug": "warfarin",
        "pharmgkb_id": "PA451906",
        "pathway_type": "Pharmacokinetics/Pharmacodynamics",
        "genes": ["CYP2C9", "VKORC1", "CYP3A4", "CYP1A2"],
        "key_metabolizing_enzyme": "CYP2C9",
        "active_metabolite": "S-warfarin (more potent enantiomer)",
        "mechanism": "CYP2C9 is the primary metabolizer of S-warfarin. VKORC1 is the pharmacodynamic target "
        "(vitamin K epoxide reductase). VKORC1 promoter variants affect warfarin sensitivity.",
        "transporter_genes": [],
    },
    "simvastatin": {
        "drug": "simvastatin",
        "pharmgkb_id": "PA451363",
        "pathway_type": "Pharmacokinetics",
        "genes": ["SLCO1B1", "CYP3A4", "CYP3A5", "ABCB1", "ABCG2"],
        "key_metabolizing_enzyme": "CYP3A4",
        "active_metabolite": "simvastatin acid",
        "mechanism": "SLCO1B1 encodes OATP1B1 hepatic uptake transporter. Reduced OATP1B1 function increases "
        "systemic exposure to simvastatin acid, raising myopathy risk.",
        "transporter_genes": ["SLCO1B1", "ABCB1", "ABCG2"],
    },
    "tamoxifen": {
        "drug": "tamoxifen",
        "pharmgkb_id": "PA451581",
        "pathway_type": "Pharmacokinetics",
        "genes": ["CYP2D6", "CYP3A4", "CYP3A5", "CYP2C9", "CYP2C19"],
        "key_metabolizing_enzyme": "CYP2D6",
        "active_metabolite": "endoxifen",
        "mechanism": "CYP2D6 is the primary enzyme for conversion of tamoxifen (via N-desmethyl-tamoxifen) "
        "to the potent active metabolite endoxifen.",
        "transporter_genes": [],
    },
    "fluorouracil": {
        "drug": "fluorouracil",
        "pharmgkb_id": "PA128406956",
        "pathway_type": "Pharmacokinetics",
        "genes": ["DPYD", "TYMS", "MTHFR"],
        "key_metabolizing_enzyme": "DPYD",
        "active_metabolite": "fluorodeoxyuridine monophosphate (FdUMP)",
        "mechanism": "DPYD (dihydropyrimidine dehydrogenase) is the rate-limiting enzyme for fluoropyrimidine "
        "catabolism. DPD deficiency leads to accumulation of 5-FU and severe toxicity.",
        "transporter_genes": [],
    },
    "azathioprine": {
        "drug": "azathioprine",
        "pharmgkb_id": "PA448515",
        "pathway_type": "Pharmacokinetics",
        "genes": ["TPMT", "NUDT15", "HPRT1", "XO"],
        "key_metabolizing_enzyme": "TPMT",
        "active_metabolite": "6-thioguanine nucleotides (6-TGN)",
        "mechanism": "TPMT methylates mercaptopurine (azathioprine metabolite), diverting it from the cytotoxic "
        "6-TGN pathway. Low TPMT activity leads to 6-TGN accumulation and myelosuppression.",
        "transporter_genes": [],
    },
}


def query_pharmgkb_annotations(
    gene: str | None = None,
    drug: str | None = None,
    variant: str | None = None,
) -> list[dict[str, Any]]:
    """Query PharmGKB clinical annotations by gene, drug, and/or variant.

    Filters the annotation database by any combination of gene, drug, and variant.
    At least one filter must be provided.

    Args:
        gene: Gene symbol (case-insensitive)
        drug: Drug name (case-insensitive)
        variant: Variant rsID or identifier

    Returns:
        List of matching annotation dictionaries

    Raises:
        ValueError: If no filter criteria provided
    """
    if gene is None and drug is None and variant is None:
        raise ValueError("At least one of gene, drug, or variant must be provided")

    results: list[dict[str, Any]] = []

    for annotation in _PHARMGKB_CLINICAL_ANNOTATIONS:
        match = True

        if gene is not None and annotation["gene"].upper() != gene.upper():
            match = False
        if drug is not None and annotation["drug"].lower() != drug.lower():
            match = False
        if variant is not None:
            ann_variant = annotation.get("variant")
            if ann_variant is None or ann_variant.lower() != variant.lower():
                match = False

        if match:
            results.append(annotation)

    logger.info(
        "PharmGKB query (gene=%s, drug=%s, variant=%s) returned %d annotations",
        gene,
        drug,
        variant,
        len(results),
    )

    return results


def parse_clinical_annotations(
    filepath: str | Path,
) -> list[dict[str, Any]]:
    """Parse PharmGKB clinical annotation data from file.

    Reads a PharmGKB clinical annotation export file (TSV format) and
    parses it into structured annotation records.

    Expected TSV columns: Clinical Annotation ID, Gene, Level of Evidence,
    Drug(s), Phenotype Category, PMID Count, Annotation Text

    Args:
        filepath: Path to PharmGKB clinical annotations TSV file

    Returns:
        List of parsed annotation dictionaries
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"PharmGKB annotation file not found: {filepath}")

    annotations: list[dict[str, Any]] = []

    with open(path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            # Handle PharmGKB TSV column name variations
            annotation_id = row.get("Clinical Annotation ID", row.get("annotation_id", ""))
            gene = row.get("Gene", row.get("gene", "")).upper()
            level = row.get("Level of Evidence", row.get("evidence_level", ""))
            drugs = row.get("Drug(s)", row.get("drug", ""))
            phenotype_cat = row.get("Phenotype Category", row.get("phenotype_category", ""))
            text = row.get("Annotation Text", row.get("annotation_text", ""))

            # Split multi-drug entries
            drug_list = [d.strip() for d in drugs.split(";") if d.strip()] if drugs else []

            for drug in drug_list or [""]:
                annotations.append(
                    {
                        "annotation_id": annotation_id,
                        "gene": gene,
                        "drug": drug.lower(),
                        "evidence_level": level,
                        "phenotype_category": phenotype_cat,
                        "annotation_text": text,
                    }
                )

    logger.info("Parsed %d clinical annotations from %s", len(annotations), filepath)
    return annotations


def get_evidence_level(
    annotation: dict[str, Any],
) -> dict[str, Any]:
    """Extract and interpret the evidence level from a PharmGKB annotation.

    Args:
        annotation: Annotation dictionary with "evidence_level" key

    Returns:
        Dictionary with:
            - "level": Evidence level string (e.g., "1A")
            - "description": Human-readable description
            - "strength": Qualitative strength ("Strong", "Moderate", "Weak", "Preliminary")
            - "actionable": Whether this evidence level supports clinical action
    """
    level = annotation.get("evidence_level", "")

    level_descriptions = {
        "1A": {
            "description": "Variant in CPIC/DPWG guideline or FDA PGx biomarker",
            "strength": "Strong",
            "actionable": True,
        },
        "1B": {
            "description": "Strong evidence, not yet in a guideline",
            "strength": "Strong",
            "actionable": True,
        },
        "2A": {
            "description": "Moderate evidence from PharmGKB clinical annotation",
            "strength": "Moderate",
            "actionable": True,
        },
        "2B": {
            "description": "Moderate evidence from PharmGKB variant annotation",
            "strength": "Moderate",
            "actionable": False,
        },
        "3": {
            "description": "Low-level evidence (case reports, small studies)",
            "strength": "Weak",
            "actionable": False,
        },
        "4": {
            "description": "Preliminary evidence (in vitro, preclinical)",
            "strength": "Preliminary",
            "actionable": False,
        },
    }

    info = level_descriptions.get(
        level,
        {
            "description": f"Unknown evidence level: {level}",
            "strength": "Unknown",
            "actionable": False,
        },
    )

    return {
        "level": level,
        "description": info["description"],
        "strength": info["strength"],
        "actionable": info["actionable"],
    }


def search_drug_pathways(
    drug: str,
) -> dict[str, Any] | None:
    """Get pharmacokinetic/pharmacodynamic pathway information for a drug.

    Args:
        drug: Drug name (case-insensitive)

    Returns:
        Pathway dictionary with genes involved, mechanism, active metabolite,
        etc. Returns None if drug not found.
    """
    drug_lower = drug.lower().strip()

    pathway = _DRUG_PATHWAYS.get(drug_lower)
    if pathway is not None:
        logger.info("Found pathway data for %s (%d genes)", drug, len(pathway["genes"]))
        return dict(pathway)

    logger.debug("No pathway data found for %s", drug)
    return None


def get_variant_annotations(
    rsid: str,
) -> dict[str, Any] | None:
    """Get variant-specific pharmacogenomic annotations.

    Args:
        rsid: Variant rsID (e.g., "rs4244285")

    Returns:
        Variant annotation dictionary with gene, allele, frequency data,
        clinical significance, etc. Returns None if variant not found.
    """
    rsid_normalized = rsid.lower().strip()

    annotation = _VARIANT_ANNOTATIONS.get(rsid_normalized)
    if annotation is not None:
        logger.info("Found variant annotation for %s (%s %s)", rsid, annotation["gene"], annotation["allele"])
        return dict(annotation)

    logger.debug("No variant annotation found for %s", rsid)
    return None
