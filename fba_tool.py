"""
fba_tool.py
===========
Agent-callable Flux Balance Analysis tool for strain-design scoring.

DESIGN GOAL
-----------
A single entry-point function, `score_pathway(...)`, that an LLM agent can call.
EVERYTHING is a parameter: the model, the carbon source, the heterologous
reactions to insert, the target product, the gene edits to apply, and the
objective. Nothing about a specific organism or pathway is hardcoded.

The function:
  1. loads a genome-scale model (default: a path you pass in)
  2. sets the carbon-source uptake (the feedstock, e.g. oil-derived fatty acid)
  3. inserts candidate heterologous reactions (e.g. FAR, WS) from a spec
  4. applies gene edits (knockouts / overexpression bounds)
  5. sets the objective (maximize product, or biomass, or a coupled mix)
  6. runs FBA
  7. returns a structured, JSON-serializable result:
       predicted product flux, growth, yield, bottleneck reactions, status

It is deliberately defensive: every step is wrapped so a bad input from the
agent returns a structured error instead of crashing the orchestration loop.

DEPENDENCIES
------------
    pip install cobra --break-system-packages

QUICK TEST
----------
    python fba_tool.py            # runs the built-in self-test on a bundled model
"""

from __future__ import annotations
import json
import re
import traceback
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import cobra
from cobra import Reaction, Metabolite

from model_loader import load_model_robust

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None


# --------------------------------------------------------------------------- #
# Model loading
# --------------------------------------------------------------------------- #
def _load_any_model(model_ref: str):
    """
    Load a model from either:
      - a bundled cobra test model name ('textbook', 'e_coli_core', ...), or
      - a path to an SBML file (.xml / .sbml / .xml.gz) e.g. iYLI647.

    SBML files are validated and auto-repaired for common id omissions before
    loading. Run find_ids.py --preflight-only on new models first.
    """
    model, _report = load_model_robust(model_ref)
    return model


# --------------------------------------------------------------------------- #
# Reaction insertion (the heterologous pathway, e.g. FAR + WS)
# --------------------------------------------------------------------------- #
def _get_or_create_metabolite(model, met_id: str, compartment: str = "c"):
    """Return an existing metabolite or create a new one on the fly."""
    if met_id in model.metabolites:
        return model.metabolites.get_by_id(met_id)
    met = Metabolite(met_id, name=met_id, compartment=compartment)
    model.add_metabolites([met])
    return met


def _insert_reactions(model, reactions_spec: List[Dict[str, Any]]):
    """
    Insert candidate reactions into the model.

    Each entry in reactions_spec is a dict, e.g.:
        {
          "id": "FAR",
          "name": "fatty acyl-CoA reductase",
          "stoichiometry": {"ocoa_c": -1, "nadph_c": -2, "foh_c": 1, ...},
          "lower_bound": 0, "upper_bound": 1000,
          "gene": "MhFAR"
        }

    Metabolite IDs that don't exist in the model are created automatically so
    the agent can describe a pathway without first checking the model's
    namespace. (In production you'd map these to real model metabolite IDs.)
    """
    added = []
    for spec in reactions_spec:
        rid = spec["id"]
        if rid in model.reactions:           # idempotent: skip if already there
            added.append(rid)
            continue
        rxn = Reaction(rid)
        rxn.name = spec.get("name", rid)
        rxn.lower_bound = spec.get("lower_bound", 0.0)
        rxn.upper_bound = spec.get("upper_bound", 1000.0)
        mets = {
            _get_or_create_metabolite(model, m, spec.get("compartment", "c")): coeff
            for m, coeff in spec["stoichiometry"].items()
        }
        rxn.add_metabolites(mets)
        if spec.get("gene"):
            rxn.gene_reaction_rule = spec["gene"]
        model.add_reactions([rxn])
        added.append(rid)
    return added


def _add_demand(model, product_met_id: str, demand_id: Optional[str] = None):
    """
    Add (or reuse) a demand/sink reaction that lets product leave the system,
    so FBA can carry flux toward it. Returns the demand reaction id.
    """
    demand_id = demand_id or f"DM_{product_met_id}"
    if demand_id in model.reactions:
        return demand_id
    met = _get_or_create_metabolite(model, product_met_id)
    dm = Reaction(demand_id)
    dm.name = f"demand {product_met_id}"
    dm.lower_bound = 0.0
    dm.upper_bound = 1000.0
    dm.add_metabolites({met: -1})
    model.add_reactions([dm])
    return demand_id


# --------------------------------------------------------------------------- #
# Gene / reaction edits
# --------------------------------------------------------------------------- #
def _apply_edits(model, knockouts: List[str], bound_overrides: Dict[str, List[float]]):
    """
    knockouts: list of reaction OR gene ids to disable.
        - if it matches a reaction id, that reaction is set to 0/0
        - else if it matches a gene id, the gene is knocked out (cobra handles
          the gene->reaction mapping)
    bound_overrides: {reaction_id: [lower, upper]} to model overexpression
        (raise upper bound) or attenuation (lower it).
    """
    applied = {"knocked_out": [], "bounds_set": [], "not_found": []}

    for ident in knockouts or []:
        if ident in model.reactions:
            model.reactions.get_by_id(ident).bounds = (0.0, 0.0)
            applied["knocked_out"].append(ident)
        elif ident in model.genes:
            model.genes.get_by_id(ident).knock_out()
            applied["knocked_out"].append(ident)
        else:
            applied["not_found"].append(ident)

    for rid, (lb, ub) in (bound_overrides or {}).items():
        if rid in model.reactions:
            model.reactions.get_by_id(rid).bounds = (lb, ub)
            applied["bounds_set"].append(rid)
        else:
            applied["not_found"].append(rid)

    return applied


# --------------------------------------------------------------------------- #
# Medium / exchange constraints (literature-calibrated simulation context)
# --------------------------------------------------------------------------- #
def _apply_exchange_constraints(
    model,
    exchange_constraints: Optional[Dict[str, List[float]]],
) -> Dict[str, Any]:
    """
    Set bounds on any exchange or boundary reaction, e.g. CER, O2, NH4+, co-substrates.
    Values are [lower_bound, upper_bound] in mmol/gDW/h (negative = export).
    """
    applied: List[str] = []
    missing: List[str] = []
    for rid, bounds in (exchange_constraints or {}).items():
        if rid not in model.reactions:
            missing.append(rid)
            continue
        lb, ub = bounds
        model.reactions.get_by_id(rid).bounds = (lb, ub)
        applied.append(rid)
    return {"exchange_bounds_set": applied, "exchange_not_found": missing}


def _reset_to_minimal_medium(model, keep_open: Optional[List[str]] = None) -> None:
    """
    Close all exchange reactions, then caller re-opens the feedstock/co-substrates.
    Mimics minimal-medium FBA setups used in GSMM validation papers.
    """
    keep = set(keep_open or [])
    for ex in model.exchanges:
        if ex.id not in keep:
            ex.bounds = (0.0, 0.0)


# Known metabolite C counts when SBML formula is missing (iYLI647 lipids).
_CARBON_FALLBACK: Dict[str, int] = {
    "ocdcea_e": 18,
    "ocdcea_c": 18,
    "odecoa_c": 18,
    "oleyl_alcohol_c": 18,
    "wax_ester_c": 36,
}


def _carbon_atoms(met) -> Optional[int]:
    """Return number of carbon atoms in a metabolite (None if unknown)."""
    mid = getattr(met, "id", str(met))
    if mid in _CARBON_FALLBACK:
        return _CARBON_FALLBACK[mid]
    formula = (getattr(met, "formula", None) or "").strip()
    if not formula:
        return None
    m = re.search(r"C(\d+)", formula.replace(" ", ""))
    if m:
        return int(m.group(1))
    if formula == "C":
        return 1
    return None


def audit_exchange_carbon(
    model,
    fluxes,
    *,
    feedstock_rxn: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Steady-state carbon audit over boundary exchanges.

    At FBA steady state, internal pools cannot be net carbon sources. If product
    carbon exceeds named feedstock import, check side_door_carbon_imports first.
    """
    imports: List[Dict[str, Any]] = []
    exports: List[Dict[str, Any]] = []

    for ex in model.exchanges:
        flux = float(fluxes[ex.id])
        if abs(flux) < 1e-9:
            continue
        c_flux = 0.0
        for met, coeff in ex.metabolites.items():
            n = _carbon_atoms(met)
            if n:
                c_flux += flux * coeff * n
        if abs(c_flux) < 1e-9:
            continue
        entry = {
            "exchange": ex.id,
            "flux_mmol_per_h": round(flux, 4),
            "carbon_mmol_per_h": round(c_flux, 4),
        }
        if c_flux > 0:
            imports.append(entry)
        else:
            exports.append(entry)

    imports.sort(key=lambda x: -x["carbon_mmol_per_h"])
    exports.sort(key=lambda x: x["carbon_mmol_per_h"])
    side_doors = [e for e in imports if e["exchange"] != feedstock_rxn]
    total_c_in = sum(e["carbon_mmol_per_h"] for e in imports)
    total_c_out = sum(-e["carbon_mmol_per_h"] for e in exports)

    audit_warnings: List[str] = []
    if side_doors:
        ids = ", ".join(e["exchange"] for e in side_doors[:5])
        audit_warnings.append(
            f"Non-feedstock carbon imports detected ({ids}). "
            "Close exchanges (use_minimal_medium=true) before trusting yield."
        )

    return {
        "carbon_imports": imports,
        "carbon_exports": exports,
        "total_carbon_import_mmol_per_h": round(total_c_in, 4),
        "total_carbon_export_mmol_per_h": round(total_c_out, 4),
        "side_door_carbon_imports": side_doors,
        "feedstock_is_sole_carbon_source": len(side_doors) == 0,
        "feedstock_rxn": feedstock_rxn,
        "warnings": audit_warnings,
    }


def _apply_growth_constraints(
    model,
    biomass_rxn: Optional[str],
    *,
    objective: str,
    min_growth_fraction: float,
    min_growth_rate: Optional[float],
    max_growth_rate: Optional[float],
) -> Dict[str, Any]:
    """Apply literature-derived growth floors/caps before product optimization."""
    info: Dict[str, Any] = {
        "max_growth_unconstrained": None,
        "min_growth_applied": None,
        "max_growth_applied": None,
    }
    if not biomass_rxn or biomass_rxn not in model.reactions:
        return info

    bio = model.reactions.get_by_id(biomass_rxn)
    info["max_growth_unconstrained"] = round(float(model.slim_optimize() or 0.0), 4)

    if objective != "biomass":
        floor = (
            float(min_growth_rate)
            if min_growth_rate is not None
            else min_growth_fraction * (info["max_growth_unconstrained"] or 0.0)
        )
        bio.lower_bound = floor
        info["min_growth_applied"] = round(floor, 4)

    if max_growth_rate is not None:
        cap = float(max_growth_rate)
        bio.upper_bound = min(bio.upper_bound, cap)
        info["max_growth_applied"] = round(bio.upper_bound, 4)

    return info


def load_simulation_scenario(path: str | Path) -> Dict[str, Any]:
    """
    Load a YAML scenario file (organism-agnostic simulation context).

    Scenario files hold literature-derived bounds so agents can swap conditions
    without changing code. See scenarios/README examples.
    """
    if yaml is None:
        raise ImportError("PyYAML is required for load_simulation_scenario (pip install pyyaml)")
    with open(path, encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    if not isinstance(data, dict):
        raise ValueError(f"Scenario file must be a mapping: {path}")
    return data


def _merge_scenario(
    scenario: Optional[Dict[str, Any]],
    *,
    carbon_source_rxn: Optional[str],
    carbon_source_uptake: float,
    exchange_constraints: Optional[Dict[str, List[float]]],
    biomass_rxn: Optional[str],
    knockouts: Optional[List[str]],
    bound_overrides: Optional[Dict[str, List[float]]],
    min_growth_fraction: float,
    min_growth_rate: Optional[float],
    max_growth_rate: Optional[float],
    substrate_moles_per_product: float,
    use_minimal_medium: bool,
    simulation_notes: Optional[List[str]],
) -> Dict[str, Any]:
    """Scenario YAML overrides explicit kwargs where provided."""
    if not scenario:
        return {
            "carbon_source_rxn": carbon_source_rxn,
            "carbon_source_uptake": carbon_source_uptake,
            "exchange_constraints": exchange_constraints or {},
            "biomass_rxn": biomass_rxn,
            "knockouts": knockouts or [],
            "bound_overrides": bound_overrides or {},
            "min_growth_fraction": min_growth_fraction,
            "min_growth_rate": min_growth_rate,
            "max_growth_rate": max_growth_rate,
            "substrate_moles_per_product": substrate_moles_per_product,
            "use_minimal_medium": use_minimal_medium,
            "simulation_notes": simulation_notes or [],
        "literature_refs": [],
        "product_literature_refs": [],
    }

    merged = {
        "carbon_source_rxn": scenario.get("carbon_source_rxn", carbon_source_rxn),
        "carbon_source_uptake": scenario.get("carbon_source_uptake", carbon_source_uptake),
        "exchange_constraints": {
            **(exchange_constraints or {}),
            **(scenario.get("exchange_constraints") or {}),
        },
        "biomass_rxn": scenario.get("biomass_rxn", biomass_rxn),
        "knockouts": list(scenario.get("knockouts") or []) + list(knockouts or []),
        "bound_overrides": {
            **(bound_overrides or {}),
            **(scenario.get("bound_overrides") or {}),
        },
        "min_growth_fraction": scenario.get("min_growth_fraction", min_growth_fraction),
        "min_growth_rate": scenario.get("min_growth_rate", min_growth_rate),
        "max_growth_rate": scenario.get("max_growth_rate", max_growth_rate),
        "substrate_moles_per_product": scenario.get(
            "substrate_moles_per_product", substrate_moles_per_product
        ),
        "use_minimal_medium": scenario.get("use_minimal_medium", use_minimal_medium),
        "simulation_notes": (scenario.get("simulation_notes") or []) + (simulation_notes or []),
        "literature_refs": scenario.get("literature_refs") or [],
        "product_literature_refs": scenario.get("product_literature_refs") or [],
    }
    return merged


def _assess_calibration(
    ctx: Dict[str, Any],
    growth_info: Dict[str, Any],
    *,
    status: str,
    yield_corrected: Optional[float],
    objective: str = "product",
    has_product: bool = False,
    carbon_audit: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Tell upstream agents how much to trust a run when literature is sparse.

    Distinguishes **medium/growth calibration** (uptake, CER, mu caps) from
    **product calibration** (experimental flux/titer for the target product).
    A wax run can be medium-calibrated while product flux remains unvalidated.
    """
    missing: List[str] = []
    warnings: List[str] = []

    if not ctx.get("literature_refs"):
        missing.append("literature_refs (no citation-backed scenario)")
    if not ctx.get("exchange_constraints"):
        missing.append("exchange_constraints (CER/O2/NH4/co-feed not pinned)")
    if ctx.get("max_growth_rate") is None:
        missing.append("max_growth_rate (growth not capped to experimental mu)")
    if carbon_audit and not carbon_audit.get("feedstock_is_sole_carbon_source"):
        missing.append(
            "closed_carbon_exchanges (non-feedstock carbon import at boundary)"
        )
    if has_product and not ctx.get("product_literature_refs"):
        missing.append(
            "product_literature_refs (no experimental anchor for target product flux)"
        )
    if ctx.get("min_growth_rate") is None and ctx.get("min_growth_fraction") == 0.1:
        warnings.append(
            "min_growth uses default fraction 0.1 of model max (often unrealistic)"
        )
    if not ctx.get("use_minimal_medium"):
        warnings.append("use_minimal_medium=false; open medium may inflate flux")
    if ctx.get("substrate_moles_per_product", 1.0) <= 1.0:
        warnings.append(
            "substrate_moles_per_product not set above 1; multi-substrate products "
            "may show misleading raw yield"
        )

    mu_unc = growth_info.get("max_growth_unconstrained")
    if mu_unc is not None and mu_unc > 1.0:
        warnings.append(
            f"model max growth {mu_unc}/h is very high; calibrate medium before trusting mu"
        )
    if yield_corrected is not None and yield_corrected > 1.0:
        warnings.append("corrected yield > 1.0; likely internal pool contribution")
    if has_product and not ctx.get("product_literature_refs"):
        warnings.append(
            "Product flux is not anchored to wax-ester (or target) experimental data; "
            "medium constraints alone do not validate product prediction."
        )
    if carbon_audit and not carbon_audit.get("feedstock_is_sole_carbon_source"):
        for w in carbon_audit.get("warnings") or []:
            warnings.append(w)
        if not ctx.get("use_minimal_medium"):
            warnings.append(
                "use_minimal_medium=false with side-door carbon imports; "
                "yield vs feedstock is not meaningful."
            )
    if status != "optimal":
        warnings.append(f"FBA status is {status}; flux values not usable")

    n_missing = len(missing)
    medium_missing = [m for m in missing if not m.startswith("product_literature")]
    n_medium_missing = len(medium_missing)

    if status != "optimal":
        medium_level = "invalid"
        product_level = "invalid"
        level = "invalid"
        use_for = "debug_constraints_only"
    elif n_medium_missing >= 3:
        medium_level = "exploratory"
        product_level = "unvalidated" if has_product else "not_applicable"
        level = "exploratory"
        use_for = "relative_ranking_of_designs_only"
    elif n_medium_missing >= 1:
        medium_level = "partial"
        product_level = "unvalidated" if has_product else "not_applicable"
        level = "partial"
        use_for = "rank_designs_and_identify_bottlenecks; not experimental titers"
    else:
        medium_level = "literature_calibrated"
        if not has_product or objective == "biomass":
            product_level = "not_applicable"
            level = "literature_calibrated"
            use_for = "quantitative_comparison_after_biomass_validation"
        elif ctx.get("product_literature_refs"):
            product_level = "literature_calibrated"
            level = "literature_calibrated"
            use_for = "quantitative_product_comparison_after_biomass_validation"
        else:
            product_level = "unvalidated"
            level = "medium_calibrated"
            use_for = (
                "medium_and_growth_are_literature_pinned; "
                "product_flux_is_in_silico_only_not_experimentally_anchored"
            )

    return {
        "confidence_level": level,
        "medium_confidence_level": medium_level,
        "product_confidence_level": product_level,
        "recommended_use": use_for,
        "missing_literature_inputs": missing,
        "warnings": warnings,
        "literature_refs": ctx.get("literature_refs") or [],
        "product_literature_refs": ctx.get("product_literature_refs") or [],
        "agent_guidance": (
            "Do not report predicted_product_flux as g/L or guaranteed titer unless "
            "product_confidence_level is literature_calibrated and biomass validation "
            "passed. medium_calibrated means growth/medium only — not product yield."
        ),
    }


# --------------------------------------------------------------------------- #
# Bottleneck identification
# --------------------------------------------------------------------------- #
def _find_bottlenecks(model, product_rxn_id: str, top_n: int = 5):
    """
    A pragmatic bottleneck heuristic suitable for a hackathon demo:
    flux variability analysis (FVA) on the pathway-relevant reactions while the
    product reaction is held at its optimum. Reactions whose flux is pinned
    (min == max == bound) while limiting the product are flagged.

    Returns a ranked list of {reaction, flux, min, max, at_bound}.
    """
    try:
        from cobra.flux_analysis import flux_variability_analysis as fva
        sol = model.optimize()
        if sol.status != "optimal":
            return []
        # restrict FVA to a manageable subset: reactions carrying nonzero flux
        active = [r.id for r in model.reactions if abs(sol.fluxes[r.id]) > 1e-6]
        fr = fva(model, reaction_list=active, fraction_of_optimum=0.99)
        rows = []
        for rid in active:
            f = sol.fluxes[rid]
            lo, hi = fr.loc[rid, "minimum"], fr.loc[rid, "maximum"]
            rxn = model.reactions.get_by_id(rid)
            at_bound = (
                abs(f - rxn.lower_bound) < 1e-6 or abs(f - rxn.upper_bound) < 1e-6
            )
            span = hi - lo
            rows.append(
                {
                    "reaction": rid,
                    "flux": round(float(f), 4),
                    "min": round(float(lo), 4),
                    "max": round(float(hi), 4),
                    "span": round(float(span), 4),
                    "at_bound": bool(at_bound),
                }
            )
        # bottlenecks = pinned (tiny span) AND carrying flux, near a bound
        rows.sort(key=lambda x: (x["span"], -abs(x["flux"])))
        return rows[:top_n]
    except Exception:
        return []


# --------------------------------------------------------------------------- #
# Main agent-callable entry point
# --------------------------------------------------------------------------- #
def score_pathway(
    model_ref: str = "iYLI647.xml",
    *,
    carbon_source_rxn: Optional[str] = None,
    carbon_source_uptake: float = 10.0,
    candidate_reactions: Optional[List[Dict[str, Any]]] = None,
    product_metabolite: Optional[str] = None,
    product_demand_rxn: Optional[str] = None,
    knockouts: Optional[List[str]] = None,
    bound_overrides: Optional[Dict[str, List[float]]] = None,
    objective: str = "product",          # "product" | "biomass" | "coupled"
    min_growth_fraction: float = 0.1,
    min_growth_rate: Optional[float] = None,
    max_growth_rate: Optional[float] = None,
    biomass_rxn: Optional[str] = None,
    exchange_constraints: Optional[Dict[str, List[float]]] = None,
    use_minimal_medium: bool = False,
    substrate_moles_per_product: float = 1.0,
    scenario: Optional[Union[str, Path, Dict[str, Any]]] = None,
    simulation_notes: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """
    Score a candidate strain design with FBA. EVERYTHING is a parameter so an
    LLM agent can call this with a JSON payload and modify any aspect.

    Parameters
    ----------
    model_ref            : path to SBML model OR a bundled cobra model name.
    carbon_source_rxn    : exchange reaction id for the feedstock (e.g. an oil-
                           derived fatty-acid exchange). If None, the model's
                           default medium is used.
    carbon_source_uptake : max uptake (positive number, mmol/gDW/h). Applied as
                           lower_bound = -uptake on the exchange.
    candidate_reactions  : list of reaction specs to insert (the heterologous
                           pathway, e.g. FAR + WS). See _insert_reactions docstring.
    product_metabolite   : metabolite id of the target product (e.g. wax ester).
                           A demand reaction is auto-added so flux can reach it.
    product_demand_rxn   : explicit product-exit reaction id (overrides auto demand).
    knockouts            : reaction or gene ids to disable (e.g. POX1..6, MFE1).
    bound_overrides      : {reaction_id: [lb, ub]} for OE / attenuation
                           (e.g. attenuate DGA1, overexpress LIP2/FAA1/ZWF1).
    objective            : what to maximize.
                           "product" -> maximize product flux (growth pinned >=
                                        min_growth_fraction * max growth).
                           "biomass" -> maximize growth (sanity check).
                           "coupled" -> maximize product with growth constraint.
    min_growth_fraction  : fraction of max growth to require when maximizing product.
    min_growth_rate      : absolute growth floor (1/h); overrides min_growth_fraction
                           when set. Use literature batch μ during production phase.
    max_growth_rate      : absolute growth cap (1/h); use experimental μ from papers.
    biomass_rxn          : biomass reaction id. Auto-detected if None.
    exchange_constraints : {exchange_id: [lb, ub]} for CER, O2, NH4+, co-substrates.
    use_minimal_medium   : if True, close all exchanges then open feedstock/constraints.
    substrate_moles_per_product : biochem ratio for yield correction (e.g. 2.0 for
                           C36 wax from C18 oleate). From pathway stoichiometry.
    scenario             : path to YAML scenario or pre-loaded dict (see scenarios/).
    simulation_notes     : free-text assumptions echoed in the result for traceability.

    Returns
    -------
    dict (JSON-serializable):
        {
          "status": "optimal" | "infeasible" | "error",
          "objective_used": ...,
          "predicted_product_flux": float,
          "growth_rate": float,
          "yield_mol_per_mol_substrate": float | None,
          "yield_corrected_mol_per_mol_substrate": float | None,
          "simulation_context": {...},
          "inserted_reactions": [...],
          "edits_applied": {...},
          "bottlenecks": [ {reaction, flux, span, at_bound}, ... ],
          "message": "...",
        }
    """
    result: Dict[str, Any] = {
        "status": "error",
        "objective_used": objective,
        "predicted_product_flux": None,
        "growth_rate": None,
        "yield_mol_per_mol_substrate": None,
        "yield_corrected_mol_per_mol_substrate": None,
        "simulation_context": {},
        "calibration": {},
        "inserted_reactions": [],
        "edits_applied": {},
        "bottlenecks": [],
        "message": "",
    }

    scenario_data: Optional[Dict[str, Any]] = None
    if isinstance(scenario, (str, Path)):
        scenario_data = load_simulation_scenario(scenario)
    elif isinstance(scenario, dict):
        scenario_data = scenario

    ctx = _merge_scenario(
        scenario_data,
        carbon_source_rxn=carbon_source_rxn,
        carbon_source_uptake=carbon_source_uptake,
        exchange_constraints=exchange_constraints,
        biomass_rxn=biomass_rxn,
        knockouts=knockouts,
        bound_overrides=bound_overrides,
        min_growth_fraction=min_growth_fraction,
        min_growth_rate=min_growth_rate,
        max_growth_rate=max_growth_rate,
        substrate_moles_per_product=substrate_moles_per_product,
        use_minimal_medium=use_minimal_medium,
        simulation_notes=simulation_notes,
    )
    result["simulation_context"] = {
        k: v for k, v in ctx.items()
        if k not in ("knockouts", "bound_overrides")
    }

    try:
        model = _load_any_model(model_ref)
    except Exception as e:
        result["message"] = f"Failed to load model '{model_ref}': {e}"
        return result

    try:
        biomass_rxn = ctx["biomass_rxn"]
        carbon_source_rxn = ctx["carbon_source_rxn"]
        carbon_source_uptake = ctx["carbon_source_uptake"]

        # --- detect biomass reaction --------------------------------------- #
        if biomass_rxn is None:
            cand = [r for r in model.reactions if "biomass" in r.id.lower()]
            biomass_rxn = cand[0].id if cand else None
        if biomass_rxn is None:
            biomass_rxn = str(model.objective.expression).split("*")[-1].strip()
        result["simulation_context"]["biomass_rxn"] = biomass_rxn

        # --- medium / exchange constraints (literature-calibrated) --------- #
        if ctx["use_minimal_medium"]:
            _reset_to_minimal_medium(
                model,
                keep_open=list((ctx["exchange_constraints"] or {}).keys())
                + ([carbon_source_rxn] if carbon_source_rxn else []),
            )
        if carbon_source_rxn:
            if carbon_source_rxn in model.reactions:
                model.reactions.get_by_id(carbon_source_rxn).lower_bound = -abs(
                    carbon_source_uptake
                )
            else:
                result["message"] += (
                    f"[warn] carbon source '{carbon_source_rxn}' not in model; "
                    "using default medium. "
                )
        ex_info = _apply_exchange_constraints(model, ctx["exchange_constraints"])
        result["simulation_context"]["exchange_constraints_applied"] = ex_info

        # --- insert candidate heterologous reactions ----------------------- #
        if candidate_reactions:
            result["inserted_reactions"] = _insert_reactions(
                model, candidate_reactions
            )

        # --- make sure product can leave the system ------------------------ #
        product_rxn_id = product_demand_rxn
        if product_metabolite and not product_rxn_id:
            product_rxn_id = _add_demand(model, product_metabolite)
        elif product_demand_rxn and product_demand_rxn not in model.reactions:
            result["message"] += (
                f"[warn] product_demand_rxn '{product_demand_rxn}' not found. "
            )
            product_rxn_id = None

        # --- apply gene / reaction edits ----------------------------------- #
        result["edits_applied"] = _apply_edits(
            model, ctx["knockouts"], ctx["bound_overrides"]
        )

        # --- growth constraints then objective ----------------------------- #
        growth_info = _apply_growth_constraints(
            model,
            biomass_rxn,
            objective=objective,
            min_growth_fraction=ctx["min_growth_fraction"],
            min_growth_rate=ctx["min_growth_rate"],
            max_growth_rate=ctx["max_growth_rate"],
        )
        result["simulation_context"]["growth_constraints"] = growth_info

        if objective == "biomass":
            model.objective = biomass_rxn
        else:
            if not product_rxn_id:
                result["message"] += (
                    "Cannot maximize product: no product reaction defined. "
                    "Pass product_metabolite or product_demand_rxn. "
                )
                model.objective = biomass_rxn
            else:
                model.objective = product_rxn_id

        # --- solve --------------------------------------------------------- #
        sol = model.optimize()
        result["status"] = sol.status

        yield_corrected: Optional[float] = None

        has_product = bool(
            product_metabolite or product_demand_rxn or candidate_reactions
        )

        if sol.status != "optimal":
            result["message"] += "Model infeasible under these constraints. "
            result["calibration"] = _assess_calibration(
                ctx,
                growth_info,
                status=result["status"],
                yield_corrected=None,
                objective=objective,
                has_product=has_product,
            )
            return result

        # --- read outputs -------------------------------------------------- #
        prod_flux = (
            float(sol.fluxes[product_rxn_id]) if product_rxn_id else None
        )
        growth = (
            float(sol.fluxes[biomass_rxn])
            if biomass_rxn and biomass_rxn in model.reactions
            else None
        )
        result["predicted_product_flux"] = (
            round(prod_flux, 4) if prod_flux is not None else None
        )
        result["growth_rate"] = round(growth, 4) if growth is not None else None

        # crude molar yield vs substrate uptake
        if carbon_source_rxn and carbon_source_rxn in model.reactions and prod_flux:
            uptake = abs(sol.fluxes[carbon_source_rxn])
            if uptake > 1e-9:
                raw_yield = prod_flux / uptake
                result["yield_mol_per_mol_substrate"] = round(raw_yield, 4)
                denom = ctx["substrate_moles_per_product"] or 1.0
                result["yield_corrected_mol_per_mol_substrate"] = round(
                    raw_yield / denom, 4
                )
                yield_corrected = result["yield_corrected_mol_per_mol_substrate"]
                if raw_yield / denom > 1.0 + 1e-6:
                    result["message"] += (
                        "[warn] corrected yield > 1.0; check substrate_moles_per_product "
                        "or block internal lipid mobilization. "
                    )

        # --- bottlenecks --------------------------------------------------- #
        carbon_audit: Optional[Dict[str, Any]] = None
        if has_product or objective == "product":
            carbon_audit = audit_exchange_carbon(
                model, sol.fluxes, feedstock_rxn=carbon_source_rxn
            )
            result["carbon_audit"] = carbon_audit
            if carbon_audit.get("warnings"):
                result["message"] += "[warn] " + " ".join(carbon_audit["warnings"]) + " "

        if product_rxn_id:
            result["bottlenecks"] = _find_bottlenecks(model, product_rxn_id)

        result["calibration"] = _assess_calibration(
            ctx,
            growth_info,
            status=result["status"],
            yield_corrected=yield_corrected,
            objective=objective,
            has_product=has_product,
            carbon_audit=carbon_audit,
        )
        result["message"] += "OK"
        return result

    except Exception:
        result["status"] = "error"
        result["message"] = "Exception during scoring:\n" + traceback.format_exc()
        return result


# --------------------------------------------------------------------------- #
# Self-test (runs on a bundled model so it works with no iYLI647 download)
# --------------------------------------------------------------------------- #
def _self_test():
    """
    Demonstrates the full agent call pattern on the bundled e_coli_core model.
    Inserts a toy 2-step pathway off a real core metabolite, knocks out a
    reaction, and maximizes the product. Proves every code path runs.
    """
    print("=" * 70)
    print("SELF-TEST on bundled 'textbook' E. coli model (stand-in for iYLI647)")
    print("=" * 70)

    # toy pathway: pyruvate (pyr_c, a real core metabolite) -> intA -> product
    candidate = [
        {
            "id": "TOY_step1",
            "name": "toy reductase",
            "stoichiometry": {"pyr_c": -1, "nadph_c": -1, "intA_c": 1},
        },
        {
            "id": "TOY_step2",
            "name": "toy synthase",
            "stoichiometry": {"intA_c": -1, "prod_c": 1},
        },
    ]

    out = score_pathway(
        "textbook",
        carbon_source_rxn="EX_glc__D_e",
        carbon_source_uptake=10.0,
        candidate_reactions=candidate,
        product_metabolite="prod_c",
        knockouts=["LDH_D"],          # disable a real core reaction
        bound_overrides={"PGI": [-50, 50]},   # demo OE/attenuation path
        objective="product",
        min_growth_fraction=0.1,
    )
    print(json.dumps(out, indent=2))

    print("\n--- sanity: biomass-only objective (should match ~0.87) ---")
    out2 = score_pathway("textbook", objective="biomass")
    print("growth:", out2["growth_rate"], "| status:", out2["status"])


if __name__ == "__main__":
    _self_test()
