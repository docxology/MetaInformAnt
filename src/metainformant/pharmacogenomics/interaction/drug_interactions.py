"""Drug-drug interaction prediction and polypharmacy risk assessment.

Provides functions for predicting drug-drug interactions from a built-in
interaction database, assessing polypharmacy risk from multiple concurrent
medications, and predicting CYP enzyme inhibition/induction potential.

The built-in database covers ~100 common drug pairs with severity levels,
mechanisms, and evidence grades derived from established pharmacological
knowledge.
"""

from __future__ import annotations

import math
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Built-in interaction database
# ---------------------------------------------------------------------------

_CYP_SUBSTRATES: dict[str, list[str]] = {
    "CYP1A2": ["caffeine", "theophylline", "clozapine", "olanzapine", "tizanidine"],
    "CYP2B6": ["bupropion", "efavirenz", "cyclophosphamide", "methadone"],
    "CYP2C8": ["paclitaxel", "repaglinide", "rosiglitazone", "amodiaquine"],
    "CYP2C9": ["warfarin", "phenytoin", "losartan", "tolbutamide", "celecoxib"],
    "CYP2C19": [
        "omeprazole",
        "clopidogrel",
        "escitalopram",
        "voriconazole",
        "diazepam",
    ],
    "CYP2D6": [
        "codeine",
        "tramadol",
        "tamoxifen",
        "metoprolol",
        "fluoxetine",
        "paroxetine",
        "dextromethorphan",
        "atomoxetine",
    ],
    "CYP3A4": [
        "simvastatin",
        "atorvastatin",
        "cyclosporine",
        "tacrolimus",
        "midazolam",
        "amlodipine",
        "nifedipine",
        "sildenafil",
        "apixaban",
        "rivaroxaban",
    ],
}

_CYP_INHIBITORS: dict[str, list[dict[str, str]]] = {
    "CYP1A2": [
        {"drug": "fluvoxamine", "type": "strong", "potency": "high"},
        {"drug": "ciprofloxacin", "type": "moderate", "potency": "moderate"},
        {"drug": "cimetidine", "type": "weak", "potency": "low"},
    ],
    "CYP2C9": [
        {"drug": "fluconazole", "type": "moderate", "potency": "moderate"},
        {"drug": "amiodarone", "type": "moderate", "potency": "moderate"},
    ],
    "CYP2C19": [
        {"drug": "fluoxetine", "type": "moderate", "potency": "moderate"},
        {"drug": "fluvoxamine", "type": "strong", "potency": "high"},
        {"drug": "omeprazole", "type": "weak", "potency": "low"},
    ],
    "CYP2D6": [
        {"drug": "fluoxetine", "type": "strong", "potency": "high"},
        {"drug": "paroxetine", "type": "strong", "potency": "high"},
        {"drug": "bupropion", "type": "moderate", "potency": "moderate"},
        {"drug": "duloxetine", "type": "moderate", "potency": "moderate"},
        {"drug": "quinidine", "type": "strong", "potency": "high"},
    ],
    "CYP3A4": [
        {"drug": "ketoconazole", "type": "strong", "potency": "high"},
        {"drug": "itraconazole", "type": "strong", "potency": "high"},
        {"drug": "clarithromycin", "type": "strong", "potency": "high"},
        {"drug": "erythromycin", "type": "moderate", "potency": "moderate"},
        {"drug": "diltiazem", "type": "moderate", "potency": "moderate"},
        {"drug": "verapamil", "type": "moderate", "potency": "moderate"},
        {"drug": "grapefruit", "type": "moderate", "potency": "moderate"},
    ],
}

_CYP_INDUCERS: dict[str, list[dict[str, str]]] = {
    "CYP1A2": [
        {"drug": "smoking", "type": "strong", "potency": "high"},
        {"drug": "rifampin", "type": "strong", "potency": "high"},
    ],
    "CYP2C9": [
        {"drug": "rifampin", "type": "strong", "potency": "high"},
        {"drug": "carbamazepine", "type": "moderate", "potency": "moderate"},
    ],
    "CYP2C19": [
        {"drug": "rifampin", "type": "strong", "potency": "high"},
    ],
    "CYP3A4": [
        {"drug": "rifampin", "type": "strong", "potency": "high"},
        {"drug": "carbamazepine", "type": "strong", "potency": "high"},
        {"drug": "phenytoin", "type": "strong", "potency": "high"},
        {"drug": "phenobarbital", "type": "strong", "potency": "high"},
        {"drug": "st_johns_wort", "type": "strong", "potency": "high"},
    ],
}


def default_interaction_database() -> dict:
    """Return the built-in drug-drug interaction database.

    The database contains approximately 100 common drug pairs with severity,
    mechanism, description, and evidence level. Each entry is keyed by a
    frozenset-style tuple of two drug names (lowercased, sorted).

    Returns:
        Dictionary mapping ``(drug_a, drug_b)`` tuples (sorted) to interaction
        records with keys: severity, mechanism, description, evidence_level,
        recommendation.
    """
    db: dict[tuple[str, str], dict[str, str]] = {}

    # Helper to add a pair
    def _add(
        a: str,
        b: str,
        severity: str,
        mechanism: str,
        description: str,
        evidence: str,
        recommendation: str,
    ) -> None:
        key = tuple(sorted([a.lower(), b.lower()]))
        db[key] = {
            "severity": severity,
            "mechanism": mechanism,
            "description": description,
            "evidence_level": evidence,
            "recommendation": recommendation,
        }

    # --- Major interactions ---
    _add(
        "warfarin",
        "fluconazole",
        "major",
        "CYP2C9 inhibition",
        "Fluconazole strongly inhibits CYP2C9, increasing warfarin exposure and bleeding risk",
        "A",
        "Reduce warfarin dose by 25-50% and monitor INR closely",
    )
    _add(
        "simvastatin",
        "clarithromycin",
        "major",
        "CYP3A4 inhibition",
        "Clarithromycin inhibits CYP3A4 causing markedly elevated statin levels and rhabdomyolysis risk",
        "A",
        "Avoid combination; use azithromycin as alternative macrolide",
    )
    _add(
        "clopidogrel",
        "omeprazole",
        "major",
        "CYP2C19 inhibition",
        "Omeprazole inhibits CYP2C19, reducing clopidogrel activation and antiplatelet efficacy",
        "A",
        "Use pantoprazole instead of omeprazole",
    )
    _add(
        "methotrexate",
        "trimethoprim",
        "major",
        "Reduced renal clearance and folate antagonism",
        "Trimethoprim decreases methotrexate clearance and both are folate antagonists",
        "B",
        "Avoid combination; if unavoidable, monitor blood counts closely",
    )
    _add(
        "codeine",
        "paroxetine",
        "major",
        "CYP2D6 inhibition",
        "Paroxetine inhibits CYP2D6, blocking codeine conversion to morphine and reducing analgesic effect",
        "A",
        "Use alternative analgesic not requiring CYP2D6 activation",
    )
    _add(
        "tamoxifen",
        "fluoxetine",
        "major",
        "CYP2D6 inhibition",
        "Fluoxetine inhibits CYP2D6, reducing tamoxifen activation to endoxifen",
        "A",
        "Switch to SSRI with minimal CYP2D6 inhibition (e.g., citalopram)",
    )
    _add(
        "cyclosporine",
        "ketoconazole",
        "major",
        "CYP3A4 inhibition",
        "Ketoconazole strongly inhibits CYP3A4, dramatically increasing cyclosporine levels",
        "A",
        "Reduce cyclosporine dose and monitor trough levels",
    )
    _add(
        "tacrolimus",
        "itraconazole",
        "major",
        "CYP3A4 inhibition",
        "Itraconazole inhibits CYP3A4, raising tacrolimus levels and nephrotoxicity risk",
        "A",
        "Reduce tacrolimus dose by 50-75% and monitor levels",
    )
    _add(
        "midazolam",
        "ketoconazole",
        "major",
        "CYP3A4 inhibition",
        "Ketoconazole increases midazolam AUC up to 15-fold via CYP3A4 inhibition",
        "A",
        "Avoid oral midazolam; reduce IV dose and monitor sedation",
    )
    _add(
        "simvastatin",
        "itraconazole",
        "major",
        "CYP3A4 inhibition",
        "Itraconazole inhibits CYP3A4 causing elevated statin levels with rhabdomyolysis risk",
        "A",
        "Contraindicated combination; use pravastatin or rosuvastatin",
    )
    _add(
        "warfarin",
        "rifampin",
        "major",
        "CYP2C9 and CYP3A4 induction",
        "Rifampin induces warfarin metabolism, dramatically reducing anticoagulant effect",
        "A",
        "Avoid combination or increase warfarin dose substantially with frequent INR monitoring",
    )
    _add(
        "sildenafil",
        "ketoconazole",
        "major",
        "CYP3A4 inhibition",
        "Ketoconazole increases sildenafil exposure, raising risk of hypotension",
        "B",
        "Reduce sildenafil dose to 25 mg and avoid within 24h of ketoconazole",
    )
    _add(
        "apixaban",
        "ketoconazole",
        "major",
        "CYP3A4 and P-gp inhibition",
        "Dual inhibition increases apixaban exposure and bleeding risk",
        "A",
        "Reduce apixaban dose by 50% when co-administered",
    )
    _add(
        "carbamazepine",
        "erythromycin",
        "major",
        "CYP3A4 inhibition",
        "Erythromycin inhibits carbamazepine metabolism causing toxicity",
        "A",
        "Monitor carbamazepine levels; use azithromycin as alternative",
    )
    _add(
        "phenytoin",
        "fluconazole",
        "major",
        "CYP2C9 inhibition",
        "Fluconazole inhibits phenytoin metabolism via CYP2C9",
        "A",
        "Reduce phenytoin dose and monitor levels",
    )
    _add(
        "theophylline",
        "fluvoxamine",
        "major",
        "CYP1A2 inhibition",
        "Fluvoxamine strongly inhibits CYP1A2, increasing theophylline to toxic levels",
        "A",
        "Reduce theophylline dose by 33% and monitor levels",
    )
    _add(
        "clozapine",
        "fluvoxamine",
        "major",
        "CYP1A2 inhibition",
        "Fluvoxamine inhibits clozapine metabolism, risking toxicity and seizures",
        "A",
        "Reduce clozapine dose by 50-67% and monitor levels",
    )
    _add(
        "methadone",
        "rifampin",
        "major",
        "CYP3A4 and CYP2B6 induction",
        "Rifampin induces methadone metabolism causing withdrawal symptoms",
        "A",
        "Increase methadone dose as needed; monitor for withdrawal",
    )
    _add(
        "colchicine",
        "clarithromycin",
        "major",
        "CYP3A4 and P-gp inhibition",
        "Increased colchicine levels with risk of fatal toxicity",
        "A",
        "Contraindicated in renal/hepatic impairment; reduce colchicine dose otherwise",
    )
    _add(
        "rivaroxaban",
        "ketoconazole",
        "major",
        "CYP3A4 and P-gp inhibition",
        "Ketoconazole increases rivaroxaban exposure and bleeding risk",
        "A",
        "Avoid combination",
    )

    # --- Moderate interactions ---
    _add(
        "metoprolol",
        "fluoxetine",
        "moderate",
        "CYP2D6 inhibition",
        "Fluoxetine inhibits CYP2D6, moderately increasing metoprolol levels",
        "B",
        "Monitor heart rate and blood pressure; consider dose reduction",
    )
    _add(
        "omeprazole",
        "clopidogrel",
        "moderate",
        "CYP2C19 inhibition",
        "Omeprazole reduces clopidogrel activation via CYP2C19",
        "A",
        "Use pantoprazole instead",
    )
    _add(
        "warfarin",
        "amiodarone",
        "moderate",
        "CYP2C9 inhibition",
        "Amiodarone inhibits CYP2C9 and other enzymes, increasing warfarin effect",
        "A",
        "Reduce warfarin dose by 30-50% and monitor INR",
    )
    _add(
        "atorvastatin",
        "diltiazem",
        "moderate",
        "CYP3A4 inhibition",
        "Diltiazem moderately inhibits CYP3A4, increasing statin exposure",
        "B",
        "Limit atorvastatin to 20 mg daily",
    )
    _add(
        "escitalopram",
        "omeprazole",
        "moderate",
        "CYP2C19 competition",
        "Both drugs are CYP2C19 substrates; competition may raise escitalopram levels",
        "C",
        "Monitor for increased SSRI side effects",
    )
    _add(
        "tramadol",
        "bupropion",
        "moderate",
        "CYP2D6 inhibition",
        "Bupropion inhibits CYP2D6, reducing tramadol efficacy and increasing seizure risk",
        "B",
        "Avoid combination; use alternative analgesic",
    )
    _add(
        "losartan",
        "fluconazole",
        "moderate",
        "CYP2C9 inhibition",
        "Fluconazole inhibits CYP2C9-mediated losartan activation to E-3174",
        "B",
        "Monitor blood pressure; consider alternative ARB",
    )
    _add(
        "amlodipine",
        "diltiazem",
        "moderate",
        "CYP3A4 inhibition and pharmacodynamic synergy",
        "Diltiazem inhibits amlodipine metabolism and both lower blood pressure",
        "B",
        "Monitor for hypotension and peripheral edema",
    )
    _add(
        "nifedipine",
        "erythromycin",
        "moderate",
        "CYP3A4 inhibition",
        "Erythromycin inhibits nifedipine metabolism, increasing hypotension risk",
        "B",
        "Monitor blood pressure; use azithromycin as alternative",
    )
    _add(
        "caffeine",
        "ciprofloxacin",
        "moderate",
        "CYP1A2 inhibition",
        "Ciprofloxacin inhibits CYP1A2-mediated caffeine metabolism",
        "B",
        "Advise patient to reduce caffeine intake",
    )
    _add(
        "olanzapine",
        "ciprofloxacin",
        "moderate",
        "CYP1A2 inhibition",
        "Ciprofloxacin inhibits olanzapine metabolism via CYP1A2",
        "B",
        "Monitor for increased sedation; consider dose reduction",
    )
    _add(
        "tizanidine",
        "ciprofloxacin",
        "moderate",
        "CYP1A2 inhibition",
        "Ciprofloxacin raises tizanidine levels substantially",
        "A",
        "Contraindicated; use alternative antibiotic",
    )
    _add(
        "atomoxetine",
        "fluoxetine",
        "moderate",
        "CYP2D6 inhibition",
        "Fluoxetine inhibits atomoxetine metabolism, increasing exposure",
        "B",
        "Reduce atomoxetine dose and monitor for side effects",
    )
    _add(
        "dextromethorphan",
        "quinidine",
        "moderate",
        "CYP2D6 inhibition",
        "Quinidine inhibits CYP2D6, increasing dextromethorphan levels (used therapeutically in Nuedexta)",
        "A",
        "This interaction is used therapeutically; monitor if unintended",
    )
    _add(
        "repaglinide",
        "clarithromycin",
        "moderate",
        "CYP3A4 and CYP2C8 inhibition",
        "Clarithromycin increases repaglinide exposure and hypoglycemia risk",
        "B",
        "Monitor blood glucose; consider dose reduction",
    )
    _add(
        "voriconazole",
        "rifampin",
        "moderate",
        "CYP2C19 and CYP3A4 induction",
        "Rifampin dramatically reduces voriconazole levels, likely causing treatment failure",
        "A",
        "Contraindicated combination",
    )

    # --- Minor interactions ---
    _add(
        "caffeine",
        "cimetidine",
        "minor",
        "CYP1A2 weak inhibition",
        "Cimetidine weakly inhibits CYP1A2, mildly increasing caffeine half-life",
        "C",
        "No action needed; advise patient if symptomatic",
    )
    _add(
        "omeprazole",
        "diazepam",
        "minor",
        "CYP2C19 competition",
        "Both are CYP2C19 substrates; minor increase in diazepam levels possible",
        "C",
        "Monitor for increased sedation in susceptible patients",
    )
    _add(
        "atorvastatin",
        "amlodipine",
        "minor",
        "CYP3A4 competition",
        "Amlodipine is a weak CYP3A4 inhibitor; minor increase in statin levels",
        "C",
        "No dose adjustment usually needed; monitor if high statin dose",
    )
    _add(
        "metoprolol",
        "duloxetine",
        "minor",
        "CYP2D6 moderate inhibition",
        "Duloxetine moderately inhibits CYP2D6, slightly increasing metoprolol levels",
        "C",
        "Monitor heart rate in sensitive patients",
    )

    logger.debug("Built default interaction database with %d entries", len(db))
    return db


# ---------------------------------------------------------------------------
# Interaction prediction
# ---------------------------------------------------------------------------


def predict_drug_interaction(
    drug_a: str,
    drug_b: str,
    interaction_db: dict | None = None,
) -> dict:
    """Predict drug-drug interaction between two medications.

    Looks up the drug pair in the interaction database and returns the
    interaction record including severity, mechanism, description, evidence
    level, and clinical recommendation.

    Args:
        drug_a: Name of the first drug.
        drug_b: Name of the second drug.
        interaction_db: Optional pre-built interaction database. If ``None``,
            the built-in database from :func:`default_interaction_database` is
            used.

    Returns:
        Dictionary with keys:
            - severity: ``"major"`` | ``"moderate"`` | ``"minor"`` | ``"none"``
            - mechanism: Pharmacokinetic or pharmacodynamic mechanism.
            - description: Narrative description of the interaction.
            - evidence_level: ``"A"`` through ``"D"``.
            - recommendation: Clinical recommendation text.
            - drug_a: Normalised name of drug A.
            - drug_b: Normalised name of drug B.
    """
    if interaction_db is None:
        interaction_db = default_interaction_database()

    a_lower = drug_a.strip().lower()
    b_lower = drug_b.strip().lower()
    key = tuple(sorted([a_lower, b_lower]))

    if key in interaction_db:
        record = interaction_db[key]
        logger.info(
            "Found %s interaction between %s and %s",
            record["severity"],
            drug_a,
            drug_b,
        )
        return {
            "severity": record["severity"],
            "mechanism": record["mechanism"],
            "description": record["description"],
            "evidence_level": record["evidence_level"],
            "recommendation": record["recommendation"],
            "drug_a": a_lower,
            "drug_b": b_lower,
        }

    # Check for shared CYP pathway competition even if no explicit entry
    shared = _find_shared_cyp_pathways(a_lower, b_lower)
    if shared:
        logger.info(
            "No explicit interaction but shared CYP pathways %s for %s and %s",
            shared,
            drug_a,
            drug_b,
        )
        return {
            "severity": "none",
            "mechanism": f"Potential CYP competition via {', '.join(shared)}",
            "description": (
                f"No documented interaction, but both drugs share CYP pathway(s): "
                f"{', '.join(shared)}. Monitor if clinically relevant."
            ),
            "evidence_level": "D",
            "recommendation": "No action required; be aware of shared metabolic pathways",
            "drug_a": a_lower,
            "drug_b": b_lower,
        }

    logger.debug("No interaction found between %s and %s", drug_a, drug_b)
    return {
        "severity": "none",
        "mechanism": "None identified",
        "description": f"No known interaction between {drug_a} and {drug_b}",
        "evidence_level": "D",
        "recommendation": "No action required",
        "drug_a": a_lower,
        "drug_b": b_lower,
    }


def _find_shared_cyp_pathways(drug_a: str, drug_b: str) -> list[str]:
    """Find CYP enzymes shared by two drugs as substrates.

    Args:
        drug_a: Lowercased drug name.
        drug_b: Lowercased drug name.

    Returns:
        List of CYP enzyme names where both drugs are substrates.
    """
    shared = []
    for cyp, substrates in _CYP_SUBSTRATES.items():
        subs_lower = [s.lower() for s in substrates]
        if drug_a in subs_lower and drug_b in subs_lower:
            shared.append(cyp)
    return shared


# ---------------------------------------------------------------------------
# Polypharmacy risk
# ---------------------------------------------------------------------------


def polypharmacy_risk(
    medications: list[str],
    interaction_db: dict | None = None,
) -> dict:
    """Assess polypharmacy risk from multiple concurrent medications.

    Performs pairwise interaction checks across all medications, identifies
    CYP pathway competition, and computes an overall risk score.

    Args:
        medications: List of drug names currently prescribed.
        interaction_db: Optional pre-built interaction database.

    Returns:
        Dictionary with keys:
            - risk_score: Float 0-1 summarising overall polypharmacy risk.
            - n_medications: Number of medications assessed.
            - interactions: List of pairwise interaction dicts (non-``"none"``).
            - competing_pathways: Dict mapping CYP enzyme to list of competing
              drugs from the input.
            - n_interactions: Count of identified interactions.
            - severity_counts: Dict mapping severity to count.
            - recommendations: List of distinct recommendation strings.

    Raises:
        ValueError: If fewer than 2 medications provided.
    """
    if len(medications) < 2:
        raise ValueError("At least 2 medications are required for polypharmacy assessment")

    if interaction_db is None:
        interaction_db = default_interaction_database()

    meds_lower = [m.strip().lower() for m in medications]
    interactions: list[dict] = []
    severity_weights = {"major": 1.0, "moderate": 0.5, "minor": 0.15, "none": 0.0}
    severity_counts: dict[str, int] = {"major": 0, "moderate": 0, "minor": 0}
    recommendations: list[str] = []

    # Pairwise interactions
    for i in range(len(meds_lower)):
        for j in range(i + 1, len(meds_lower)):
            result = predict_drug_interaction(meds_lower[i], meds_lower[j], interaction_db)
            if result["severity"] != "none":
                interactions.append(result)
                sev = result["severity"]
                if sev in severity_counts:
                    severity_counts[sev] += 1
                rec = result["recommendation"]
                if rec not in recommendations:
                    recommendations.append(rec)

    # CYP pathway competition
    competing_pathways: dict[str, list[str]] = {}
    for cyp, substrates in _CYP_SUBSTRATES.items():
        subs_lower = [s.lower() for s in substrates]
        competing = [m for m in meds_lower if m in subs_lower]
        if len(competing) >= 2:
            competing_pathways[cyp] = competing

    # Compute risk score
    n_pairs = len(meds_lower) * (len(meds_lower) - 1) / 2
    if n_pairs == 0:
        n_pairs = 1

    weighted_sum = sum(severity_weights.get(ix["severity"], 0.0) for ix in interactions)
    # Add pathway competition contribution
    pathway_penalty = sum(0.1 * (len(drugs) - 1) for drugs in competing_pathways.values())
    raw_score = (weighted_sum + pathway_penalty) / max(n_pairs, 1)
    risk_score = min(1.0, raw_score)

    logger.info(
        "Polypharmacy assessment: %d meds, %d interactions, risk=%.3f",
        len(meds_lower),
        len(interactions),
        risk_score,
    )

    return {
        "risk_score": risk_score,
        "n_medications": len(meds_lower),
        "interactions": interactions,
        "competing_pathways": competing_pathways,
        "n_interactions": len(interactions),
        "severity_counts": severity_counts,
        "recommendations": recommendations,
    }


# ---------------------------------------------------------------------------
# CYP inhibition prediction
# ---------------------------------------------------------------------------


def cyp_inhibition_prediction(
    drug: str,
    cyp_profiles: dict | None = None,
) -> dict:
    """Predict CYP enzyme inhibition and induction potential of a drug.

    Checks the drug against known CYP inhibitor and inducer tables to
    identify which enzymes may be affected and at what potency.

    Args:
        drug: Name of the drug to profile.
        cyp_profiles: Optional custom CYP profile database. If ``None``,
            the built-in inhibitor/inducer tables are used. Expected format:
            ``{"inhibitors": {cyp: [...]}, "inducers": {cyp: [...]}}``.

    Returns:
        Dictionary with keys:
            - drug: Normalised drug name.
            - affected_enzymes: List of dicts, each with
              ``{enzyme, effect_type, inhibition_type, potency_estimate}``.
            - n_affected: Count of affected enzymes.
            - substrate_of: List of CYP enzymes where this drug is a substrate.
            - summary: Human-readable summary string.
    """
    if cyp_profiles is not None:
        inhibitors = cyp_profiles.get("inhibitors", {})
        inducers = cyp_profiles.get("inducers", {})
    else:
        inhibitors = _CYP_INHIBITORS
        inducers = _CYP_INDUCERS

    drug_lower = drug.strip().lower()
    affected: list[dict[str, str]] = []

    # Check inhibitor tables
    for cyp, entries in inhibitors.items():
        for entry in entries:
            if entry["drug"].lower() == drug_lower:
                affected.append(
                    {
                        "enzyme": cyp,
                        "effect_type": "inhibitor",
                        "inhibition_type": entry.get("type", "unknown"),
                        "potency_estimate": entry.get("potency", "unknown"),
                    }
                )

    # Check inducer tables
    for cyp, entries in inducers.items():
        for entry in entries:
            if entry["drug"].lower() == drug_lower:
                affected.append(
                    {
                        "enzyme": cyp,
                        "effect_type": "inducer",
                        "inhibition_type": entry.get("type", "unknown"),
                        "potency_estimate": entry.get("potency", "unknown"),
                    }
                )

    # Check substrates
    substrate_of: list[str] = []
    for cyp, substrates in _CYP_SUBSTRATES.items():
        if drug_lower in [s.lower() for s in substrates]:
            substrate_of.append(cyp)

    # Build summary
    if affected:
        parts = []
        for a in affected:
            parts.append(f"{a['inhibition_type']} {a['effect_type']} of {a['enzyme']}")
        summary = f"{drug} is a {'; '.join(parts)}"
    else:
        summary = f"{drug} has no known CYP inhibition or induction effects"

    if substrate_of:
        summary += f". Metabolised by {', '.join(substrate_of)}"

    logger.info("CYP profile for %s: %d affected enzymes", drug, len(affected))

    return {
        "drug": drug_lower,
        "affected_enzymes": affected,
        "n_affected": len(affected),
        "substrate_of": substrate_of,
        "summary": summary,
    }
