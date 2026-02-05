"""Star allele calling for pharmacogenes.

Implements star allele identification from observed variants by matching against
known allele definition tables (CPIC-standard). Supports CYP2D6 copy number
variation handling, novel allele detection, and multi-gene pharmacogene panels.

Star alleles are the standard nomenclature for pharmacogene variants, where *1
typically represents the reference (wild-type) allele, and numbered alleles (*2, *3, etc.)
represent variant haplotypes defined by specific combinations of SNPs and indels.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


@dataclass
class StarAllele:
    """Represents a pharmacogene star allele.

    Attributes:
        name: Star allele designation (e.g., "*1", "*2", "*4")
        gene: Gene symbol (e.g., "CYP2D6", "CYP2C19")
        defining_variants: Set of variant identifiers that define this allele
            Each variant is represented as "chr:pos:ref:alt" or rsID
        function: Functional classification (e.g., "Normal function",
            "No function", "Decreased function", "Increased function")
        activity_value: Numeric activity score (e.g., 1.0 for normal, 0.0 for no function)
        clinical_significance: Clinical significance annotation
        evidence_level: Strength of evidence for this allele definition
        substrate_specificity: Optional notes on substrate-specific effects
    """

    name: str
    gene: str
    defining_variants: frozenset[str] = field(default_factory=frozenset)
    function: str = "Unknown function"
    activity_value: float = 1.0
    clinical_significance: str = ""
    evidence_level: str = ""
    substrate_specificity: str = ""

    def __post_init__(self) -> None:
        """Validate and normalize star allele fields."""
        if not self.name.startswith("*"):
            self.name = f"*{self.name}"
        self.gene = self.gene.upper()
        if isinstance(self.defining_variants, (set, list, tuple)):
            object.__setattr__(self, "defining_variants", frozenset(self.defining_variants))

    @property
    def is_reference(self) -> bool:
        """Return True if this is the reference (*1) allele."""
        return self.name == "*1"

    @property
    def is_no_function(self) -> bool:
        """Return True if this allele has no function."""
        return self.activity_value == 0.0 or "no function" in self.function.lower()

    def matches_variants(self, observed: frozenset[str] | set[str]) -> bool:
        """Check if observed variants match this allele's defining variants.

        Args:
            observed: Set of observed variant identifiers

        Returns:
            True if all defining variants are present in observed set
        """
        if not self.defining_variants:
            return False
        return self.defining_variants.issubset(observed)

    def partial_match_score(self, observed: frozenset[str] | set[str]) -> float:
        """Calculate fraction of defining variants that match.

        Args:
            observed: Set of observed variant identifiers

        Returns:
            Float between 0.0 and 1.0 indicating match fraction
        """
        if not self.defining_variants:
            return 0.0
        matched = self.defining_variants.intersection(observed)
        return len(matched) / len(self.defining_variants)


# ── Built-in allele definition tables ──────────────────────────────────────────
# These represent real CPIC-derived allele definitions for key pharmacogenes.
# Each gene maps allele name -> {defining_variants, function, activity_value}.

_BUILTIN_ALLELE_DEFINITIONS: dict[str, dict[str, dict[str, Any]]] = {
    "CYP2D6": {
        "*1": {
            "defining_variants": frozenset(),
            "function": "Normal function",
            "activity_value": 1.0,
        },
        "*2": {
            "defining_variants": frozenset({"rs16947"}),
            "function": "Normal function",
            "activity_value": 1.0,
        },
        "*3": {
            "defining_variants": frozenset({"rs35742686"}),
            "function": "No function",
            "activity_value": 0.0,
        },
        "*4": {
            "defining_variants": frozenset({"rs3892097"}),
            "function": "No function",
            "activity_value": 0.0,
        },
        "*5": {
            "defining_variants": frozenset({"CYP2D6_DEL"}),
            "function": "No function",
            "activity_value": 0.0,
        },
        "*6": {
            "defining_variants": frozenset({"rs5030655"}),
            "function": "No function",
            "activity_value": 0.0,
        },
        "*9": {
            "defining_variants": frozenset({"rs5030656"}),
            "function": "Decreased function",
            "activity_value": 0.5,
        },
        "*10": {
            "defining_variants": frozenset({"rs1065852"}),
            "function": "Decreased function",
            "activity_value": 0.25,
        },
        "*17": {
            "defining_variants": frozenset({"rs28371706"}),
            "function": "Decreased function",
            "activity_value": 0.5,
        },
        "*29": {
            "defining_variants": frozenset({"rs59421388"}),
            "function": "Decreased function",
            "activity_value": 0.5,
        },
        "*36": {
            "defining_variants": frozenset({"CYP2D6_DUP", "rs16947"}),
            "function": "No function",
            "activity_value": 0.0,
        },
        "*41": {
            "defining_variants": frozenset({"rs28371725"}),
            "function": "Decreased function",
            "activity_value": 0.5,
        },
    },
    "CYP2C19": {
        "*1": {
            "defining_variants": frozenset(),
            "function": "Normal function",
            "activity_value": 1.0,
        },
        "*2": {
            "defining_variants": frozenset({"rs4244285"}),
            "function": "No function",
            "activity_value": 0.0,
        },
        "*3": {
            "defining_variants": frozenset({"rs4986893"}),
            "function": "No function",
            "activity_value": 0.0,
        },
        "*4": {
            "defining_variants": frozenset({"rs28399504"}),
            "function": "No function",
            "activity_value": 0.0,
        },
        "*17": {
            "defining_variants": frozenset({"rs12248560"}),
            "function": "Increased function",
            "activity_value": 1.5,
        },
    },
    "CYP2C9": {
        "*1": {
            "defining_variants": frozenset(),
            "function": "Normal function",
            "activity_value": 1.0,
        },
        "*2": {
            "defining_variants": frozenset({"rs1799853"}),
            "function": "Decreased function",
            "activity_value": 0.5,
        },
        "*3": {
            "defining_variants": frozenset({"rs1057910"}),
            "function": "No function",
            "activity_value": 0.0,
        },
        "*5": {
            "defining_variants": frozenset({"rs28371686"}),
            "function": "No function",
            "activity_value": 0.0,
        },
        "*6": {
            "defining_variants": frozenset({"rs9332131"}),
            "function": "No function",
            "activity_value": 0.0,
        },
        "*8": {
            "defining_variants": frozenset({"rs7900194"}),
            "function": "No function",
            "activity_value": 0.0,
        },
        "*11": {
            "defining_variants": frozenset({"rs28371685"}),
            "function": "Decreased function",
            "activity_value": 0.5,
        },
    },
    "CYP3A5": {
        "*1": {
            "defining_variants": frozenset(),
            "function": "Normal function",
            "activity_value": 1.0,
        },
        "*3": {
            "defining_variants": frozenset({"rs776746"}),
            "function": "No function",
            "activity_value": 0.0,
        },
        "*6": {
            "defining_variants": frozenset({"rs10264272"}),
            "function": "No function",
            "activity_value": 0.0,
        },
        "*7": {
            "defining_variants": frozenset({"rs41303343"}),
            "function": "No function",
            "activity_value": 0.0,
        },
    },
    "DPYD": {
        "*1": {
            "defining_variants": frozenset(),
            "function": "Normal function",
            "activity_value": 1.0,
        },
        "*2A": {
            "defining_variants": frozenset({"rs3918290"}),
            "function": "No function",
            "activity_value": 0.0,
        },
        "*13": {
            "defining_variants": frozenset({"rs55886062"}),
            "function": "No function",
            "activity_value": 0.0,
        },
        "c.2846A>T": {
            "defining_variants": frozenset({"rs67376798"}),
            "function": "Decreased function",
            "activity_value": 0.5,
        },
        "HapB3": {
            "defining_variants": frozenset({"rs56038477"}),
            "function": "Decreased function",
            "activity_value": 0.5,
        },
    },
    "TPMT": {
        "*1": {
            "defining_variants": frozenset(),
            "function": "Normal function",
            "activity_value": 1.0,
        },
        "*2": {
            "defining_variants": frozenset({"rs1800462"}),
            "function": "No function",
            "activity_value": 0.0,
        },
        "*3A": {
            "defining_variants": frozenset({"rs1800460", "rs1142345"}),
            "function": "No function",
            "activity_value": 0.0,
        },
        "*3B": {
            "defining_variants": frozenset({"rs1800460"}),
            "function": "No function",
            "activity_value": 0.0,
        },
        "*3C": {
            "defining_variants": frozenset({"rs1142345"}),
            "function": "No function",
            "activity_value": 0.0,
        },
    },
    "NUDT15": {
        "*1": {
            "defining_variants": frozenset(),
            "function": "Normal function",
            "activity_value": 1.0,
        },
        "*2": {
            "defining_variants": frozenset({"rs746071566"}),
            "function": "No function",
            "activity_value": 0.0,
        },
        "*3": {
            "defining_variants": frozenset({"rs116855232"}),
            "function": "No function",
            "activity_value": 0.0,
        },
    },
    "SLCO1B1": {
        "*1": {
            "defining_variants": frozenset(),
            "function": "Normal function",
            "activity_value": 1.0,
        },
        "*5": {
            "defining_variants": frozenset({"rs4149056"}),
            "function": "Decreased function",
            "activity_value": 0.5,
        },
        "*15": {
            "defining_variants": frozenset({"rs4149056", "rs2306283"}),
            "function": "Decreased function",
            "activity_value": 0.5,
        },
        "*1B": {
            "defining_variants": frozenset({"rs2306283"}),
            "function": "Normal function",
            "activity_value": 1.0,
        },
    },
}


def load_allele_definitions(
    gene: str,
    source: str = "cpic",
    filepath: str | Path | None = None,
) -> list[StarAllele]:
    """Load allele definition tables for a pharmacogene.

    Loads allele definitions from either the built-in CPIC-derived table or
    from an external file (TSV or JSON format). The built-in tables cover
    CYP2D6, CYP2C19, CYP2C9, CYP3A5, DPYD, TPMT, NUDT15, and SLCO1B1.

    Args:
        gene: Gene symbol (e.g., "CYP2D6", "CYP2C19")
        source: Definition source - "cpic" for built-in, "file" for external file
        filepath: Path to external allele definition file (required if source="file")

    Returns:
        List of StarAllele objects for the gene

    Raises:
        ValueError: If gene not found in definitions or invalid source
        FileNotFoundError: If filepath specified but not found
    """
    gene_upper = gene.upper()

    if source == "file" and filepath is not None:
        return _load_allele_definitions_from_file(gene_upper, filepath)

    if source == "cpic" or filepath is None:
        return _load_builtin_definitions(gene_upper)

    raise ValueError(f"Invalid source '{source}'. Use 'cpic' or 'file'.")


def _load_builtin_definitions(gene: str) -> list[StarAllele]:
    """Load allele definitions from built-in CPIC-derived tables.

    Args:
        gene: Uppercased gene symbol

    Returns:
        List of StarAllele objects

    Raises:
        ValueError: If gene not found in built-in definitions
    """
    if gene not in _BUILTIN_ALLELE_DEFINITIONS:
        available = sorted(_BUILTIN_ALLELE_DEFINITIONS.keys())
        raise ValueError(
            f"Gene '{gene}' not found in built-in definitions. "
            f"Available genes: {', '.join(available)}. "
            f"Use source='file' with a filepath for custom definitions."
        )

    alleles: list[StarAllele] = []
    for allele_name, defn in _BUILTIN_ALLELE_DEFINITIONS[gene].items():
        alleles.append(
            StarAllele(
                name=allele_name,
                gene=gene,
                defining_variants=defn["defining_variants"],
                function=defn["function"],
                activity_value=defn["activity_value"],
            )
        )

    logger.info("Loaded %d built-in allele definitions for %s", len(alleles), gene)
    return alleles


def _load_allele_definitions_from_file(gene: str, filepath: str | Path) -> list[StarAllele]:
    """Load allele definitions from an external file.

    Supports TSV format with columns: allele, defining_variants (semicolon-separated),
    function, activity_value. Also supports JSON format with the same structure as
    _BUILTIN_ALLELE_DEFINITIONS.

    Args:
        gene: Uppercased gene symbol
        filepath: Path to the definition file

    Returns:
        List of StarAllele objects
    """
    from metainformant.core import io

    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"Allele definition file not found: {filepath}")

    alleles: list[StarAllele] = []

    if path.suffix == ".json" or path.suffix == ".gz" and path.stem.endswith(".json"):
        data = io.load_json(path)
        gene_data = data.get(gene, data)

        for allele_name, defn in gene_data.items():
            variants_raw = defn.get("defining_variants", [])
            if isinstance(variants_raw, list):
                variants = frozenset(variants_raw)
            elif isinstance(variants_raw, str):
                variants = frozenset(v.strip() for v in variants_raw.split(";") if v.strip())
            else:
                variants = frozenset(variants_raw)

            alleles.append(
                StarAllele(
                    name=allele_name,
                    gene=gene,
                    defining_variants=variants,
                    function=defn.get("function", "Unknown function"),
                    activity_value=float(defn.get("activity_value", 1.0)),
                    clinical_significance=defn.get("clinical_significance", ""),
                    evidence_level=defn.get("evidence_level", ""),
                )
            )
    else:
        # TSV format
        import csv

        with open(path, "r", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                allele_name = row.get("allele", row.get("name", ""))
                variants_str = row.get("defining_variants", "")
                variants = frozenset(v.strip() for v in variants_str.split(";") if v.strip())

                alleles.append(
                    StarAllele(
                        name=allele_name,
                        gene=gene,
                        defining_variants=variants,
                        function=row.get("function", "Unknown function"),
                        activity_value=float(row.get("activity_value", 1.0)),
                        clinical_significance=row.get("clinical_significance", ""),
                        evidence_level=row.get("evidence_level", ""),
                    )
                )

    logger.info("Loaded %d allele definitions for %s from %s", len(alleles), gene, filepath)
    return alleles


def call_star_alleles(
    variants: set[str] | frozenset[str] | list[str],
    gene: str,
    allele_definitions: list[StarAllele] | None = None,
) -> list[StarAllele]:
    """Call star alleles from observed variants.

    Matches observed variant identifiers against known allele definitions to
    determine which star alleles are present. Uses a greedy algorithm that
    prioritizes alleles with the most defining variants matched.

    The algorithm:
    1. Sort allele definitions by number of defining variants (descending)
    2. For each allele, check if all defining variants are in the observed set
    3. Matched alleles are collected; variants consumed by matched alleles are tracked
    4. If no non-reference alleles match, the reference (*1) allele is returned

    Args:
        variants: Set of observed variant identifiers (rsIDs or chr:pos:ref:alt)
        gene: Gene symbol (e.g., "CYP2D6")
        allele_definitions: Pre-loaded allele definitions. If None, loads from built-in.

    Returns:
        List of matched StarAllele objects, sorted by allele name.
        Returns [*1] if no variant alleles match.
    """
    gene_upper = gene.upper()
    observed = frozenset(variants) if not isinstance(variants, frozenset) else variants

    if allele_definitions is None:
        allele_definitions = load_allele_definitions(gene_upper)

    # Separate reference from variant alleles
    reference_allele: StarAllele | None = None
    variant_alleles: list[StarAllele] = []

    for allele_def in allele_definitions:
        if allele_def.is_reference:
            reference_allele = allele_def
        elif allele_def.defining_variants:
            variant_alleles.append(allele_def)

    # Sort by number of defining variants descending (greedy: match most specific first)
    variant_alleles.sort(key=lambda a: len(a.defining_variants), reverse=True)

    matched: list[StarAllele] = []
    consumed_variants: set[str] = set()

    for allele_def in variant_alleles:
        # Check if all defining variants are present and not yet consumed
        remaining_observed = observed - consumed_variants
        if allele_def.defining_variants.issubset(remaining_observed):
            matched.append(allele_def)
            consumed_variants.update(allele_def.defining_variants)
            logger.debug(
                "Matched allele %s for %s (variants: %s)",
                allele_def.name,
                gene_upper,
                allele_def.defining_variants,
            )

    if not matched and reference_allele is not None:
        matched.append(reference_allele)
        logger.debug("No variant alleles matched for %s, returning reference *1", gene_upper)

    # Sort by allele name for consistent ordering
    matched.sort(key=lambda a: a.name)

    logger.info(
        "Called %d star allele(s) for %s from %d observed variants: %s",
        len(matched),
        gene_upper,
        len(observed),
        ", ".join(a.name for a in matched),
    )

    return matched


def match_allele_definition(
    observed_variants: set[str] | frozenset[str],
    allele_definition: StarAllele,
) -> dict[str, Any]:
    """Match observed variants against a single allele definition.

    Performs detailed comparison between observed variants and an allele's
    defining variant set, returning match statistics.

    Args:
        observed_variants: Set of observed variant identifiers
        allele_definition: StarAllele definition to match against

    Returns:
        Dictionary with match results:
            - "allele": Star allele name
            - "gene": Gene symbol
            - "is_match": True if all defining variants are present
            - "match_score": Fraction of defining variants matched (0.0 - 1.0)
            - "matched_variants": Set of matched variant IDs
            - "missing_variants": Set of unmatched defining variants
            - "extra_variants": Observed variants not in the definition
    """
    observed = frozenset(observed_variants) if not isinstance(observed_variants, frozenset) else observed_variants
    defining = allele_definition.defining_variants

    if not defining:
        return {
            "allele": allele_definition.name,
            "gene": allele_definition.gene,
            "is_match": False,
            "match_score": 0.0,
            "matched_variants": set(),
            "missing_variants": set(),
            "extra_variants": set(observed),
        }

    matched = defining.intersection(observed)
    missing = defining - observed
    extra = observed - defining

    return {
        "allele": allele_definition.name,
        "gene": allele_definition.gene,
        "is_match": len(missing) == 0,
        "match_score": len(matched) / len(defining) if defining else 0.0,
        "matched_variants": set(matched),
        "missing_variants": set(missing),
        "extra_variants": set(extra),
    }


def detect_novel_alleles(
    variants: set[str] | frozenset[str],
    known_alleles: list[StarAllele],
) -> dict[str, Any]:
    """Detect potentially novel allele combinations.

    Identifies observed variant combinations that do not exactly match any
    known allele definition, which may represent novel or rare alleles.

    Args:
        variants: Set of observed variant identifiers
        known_alleles: List of known StarAllele definitions

    Returns:
        Dictionary with:
            - "has_novel": True if novel combinations detected
            - "unmatched_variants": Variants not explained by any known allele
            - "partial_matches": Alleles with partial but incomplete matches
            - "closest_allele": Best partial match allele name, or None
            - "closest_score": Match score of closest allele
    """
    observed = frozenset(variants) if not isinstance(variants, frozenset) else variants

    # Find all variants explained by known alleles
    explained_variants: set[str] = set()
    partial_matches: list[dict[str, Any]] = []
    best_score = 0.0
    closest_allele: str | None = None

    for allele in known_alleles:
        if not allele.defining_variants:
            continue

        score = allele.partial_match_score(observed)
        if score > 0.0 and score < 1.0:
            partial_matches.append({
                "allele": allele.name,
                "gene": allele.gene,
                "match_score": score,
                "matched_variants": sorted(allele.defining_variants.intersection(observed)),
                "missing_variants": sorted(allele.defining_variants - observed),
            })
        if score >= 1.0:
            explained_variants.update(allele.defining_variants)

        if score > best_score:
            best_score = score
            closest_allele = allele.name

    unmatched = observed - explained_variants
    # Filter to variants that look like pharmacogene-relevant IDs (rsIDs or positional)
    pharmacogene_unmatched = {v for v in unmatched if v.startswith("rs") or ":" in v}

    has_novel = len(pharmacogene_unmatched) > 0

    if has_novel:
        logger.info(
            "Detected %d unmatched variants that may represent novel alleles: %s",
            len(pharmacogene_unmatched),
            sorted(pharmacogene_unmatched),
        )

    return {
        "has_novel": has_novel,
        "unmatched_variants": sorted(pharmacogene_unmatched),
        "partial_matches": sorted(partial_matches, key=lambda x: x["match_score"], reverse=True),
        "closest_allele": closest_allele,
        "closest_score": best_score,
    }


def handle_cyp2d6_cnv(
    variants: set[str] | frozenset[str],
    copy_number: int,
) -> dict[str, Any]:
    """Handle CYP2D6 copy number variation for star allele calling.

    CYP2D6 is subject to gene deletion (*5), gene duplication/multiplication,
    and hybrid gene arrangements. This function adjusts star allele calls
    based on observed copy number.

    Copy number interpretation:
        - 0: Homozygous gene deletion (*5/*5)
        - 1: Hemizygous (one allele deleted): allele/*5
        - 2: Normal diploid (standard calling)
        - 3+: Gene duplication/multiplication

    For duplications (CN >= 3), the activity score of the duplicated allele
    is multiplied by the additional copy count. Per CPIC guidelines, functional
    allele duplications are denoted with "xN" suffix (e.g., *1x2).

    Args:
        variants: Set of observed variant identifiers
        copy_number: CYP2D6 gene copy number (0, 1, 2, 3+)

    Returns:
        Dictionary with:
            - "copy_number": Input copy number
            - "alleles": List of called StarAllele objects (adjusted for CNV)
            - "diplotype_string": Formatted diplotype string
            - "total_activity_score": Sum of all allele activity values
            - "cnv_type": "deletion", "normal", or "duplication"
            - "notes": Clinical notes about the CNV
    """
    gene = "CYP2D6"
    observed = frozenset(variants) if not isinstance(variants, frozenset) else variants

    allele_definitions = load_allele_definitions(gene)
    deletion_allele = StarAllele(
        name="*5",
        gene=gene,
        defining_variants=frozenset({"CYP2D6_DEL"}),
        function="No function",
        activity_value=0.0,
    )

    if copy_number <= 0:
        # Homozygous deletion
        return {
            "copy_number": 0,
            "alleles": [deletion_allele, deletion_allele],
            "diplotype_string": "*5/*5",
            "total_activity_score": 0.0,
            "cnv_type": "deletion",
            "notes": "Homozygous CYP2D6 gene deletion. No enzyme activity expected.",
        }

    # Call alleles from observed variants (non-CNV markers)
    non_cnv_variants = {v for v in observed if v not in {"CYP2D6_DEL", "CYP2D6_DUP"}}
    called = call_star_alleles(non_cnv_variants, gene, allele_definitions)

    if copy_number == 1:
        # Hemizygous: one called allele + *5 deletion
        primary_allele = called[0] if called else StarAllele(
            name="*1", gene=gene, function="Normal function", activity_value=1.0
        )
        alleles = sorted([primary_allele, deletion_allele], key=lambda a: a.name)
        total_score = primary_allele.activity_value + deletion_allele.activity_value

        return {
            "copy_number": 1,
            "alleles": alleles,
            "diplotype_string": f"{alleles[0].name}/{alleles[1].name}",
            "total_activity_score": total_score,
            "cnv_type": "deletion",
            "notes": f"Hemizygous CYP2D6. One gene copy with {primary_allele.name} ({primary_allele.function}).",
        }

    if copy_number == 2:
        # Normal diploid
        if len(called) >= 2:
            allele1, allele2 = called[0], called[1]
        elif len(called) == 1:
            allele1 = called[0]
            # Second allele defaults to reference
            allele2 = StarAllele(
                name="*1", gene=gene, function="Normal function", activity_value=1.0
            )
        else:
            allele1 = StarAllele(name="*1", gene=gene, function="Normal function", activity_value=1.0)
            allele2 = StarAllele(name="*1", gene=gene, function="Normal function", activity_value=1.0)

        alleles = sorted([allele1, allele2], key=lambda a: a.name)
        total_score = allele1.activity_value + allele2.activity_value

        return {
            "copy_number": 2,
            "alleles": alleles,
            "diplotype_string": f"{alleles[0].name}/{alleles[1].name}",
            "total_activity_score": total_score,
            "cnv_type": "normal",
            "notes": "Normal diploid CYP2D6.",
        }

    # copy_number >= 3: Duplication/Multiplication
    extra_copies = copy_number - 2

    if len(called) >= 2:
        allele1, allele2 = called[0], called[1]
    elif len(called) == 1:
        allele1 = called[0]
        allele2 = StarAllele(name="*1", gene=gene, function="Normal function", activity_value=1.0)
    else:
        allele1 = StarAllele(name="*1", gene=gene, function="Normal function", activity_value=1.0)
        allele2 = StarAllele(name="*1", gene=gene, function="Normal function", activity_value=1.0)

    # Per CPIC: assign extra copies to the allele with higher activity
    # (conservative: assume the functional allele is duplicated)
    if allele1.activity_value >= allele2.activity_value:
        duplicated_allele = allele1
        other_allele = allele2
    else:
        duplicated_allele = allele2
        other_allele = allele1

    # Total activity: normal alleles + extra copies of duplicated allele
    total_score = (
        duplicated_allele.activity_value * (1 + extra_copies)
        + other_allele.activity_value
    )

    dup_name = f"{duplicated_allele.name}x{1 + extra_copies}"
    alleles_sorted = sorted([duplicated_allele, other_allele], key=lambda a: a.name)
    diplotype_str = f"{dup_name}/{other_allele.name}"

    return {
        "copy_number": copy_number,
        "alleles": [duplicated_allele, other_allele],
        "diplotype_string": diplotype_str,
        "total_activity_score": total_score,
        "cnv_type": "duplication",
        "notes": (
            f"CYP2D6 gene duplication detected (CN={copy_number}). "
            f"{duplicated_allele.name} duplicated {extra_copies} extra time(s). "
            f"Total activity score: {total_score:.2f}."
        ),
    }
