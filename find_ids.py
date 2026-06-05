"""
find_ids.py
===========
Resolve REAL metabolite / reaction / gene IDs inside a genome-scale model so
upstream agents can call score_pathway() with names that exist in the network.

ALWAYS run this before score_pathway() on a new model or design.

USAGE
-----
    python find_ids.py iYLI647.xml
    python find_ids.py iYLI647.xml --term oleate --term pox
    python find_ids.py iYLI647.xml --json > ids.json    # agent-friendly output
    python find_ids.py iYLI647.xml --preflight-only     # validate/repair SBML only
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from typing import Any, Dict, List, Optional, Tuple

from model_loader import (
    IYLI647_BIOMASS_CARBON,
    ensure_sbml_loadable,
    is_sbml_path,
    load_model_robust,
)

# Upstream agents often propose yeast *gene* names; many SBML models only have
# reaction ids. Map common names to search terms for reaction resolution.
GENE_ALIASES: Dict[str, List[str]] = {
    "POX1": ["acyl coa oxidase", "pox", "peroxisomal oxidase"],
    "POX2": ["acyl coa oxidase", "pox", "peroxisomal oxidase"],
    "POX3": ["acyl coa oxidase", "pox", "peroxisomal oxidase"],
    "POX4": ["acyl coa oxidase", "pox", "peroxisomal oxidase"],
    "POX5": ["acyl coa oxidase", "pox", "peroxisomal oxidase"],
    "POX6": ["acyl coa oxidase", "pox", "peroxisomal oxidase"],
    "MFE1": ["mfe", "multifunctional enzyme"],
    "DGA1": ["diacylglycerol acyltransferase", "dga"],
    "LRO1": ["triacylglycerol lipase", "lro"],
    "ARE1": ["acyl coa cholesterol", "are"],
    "FAA1": ["fatty acid coa ligase", "facoa", "faa"],
    "LIP2": ["lipase", "lip"],
    "ZWF1": ["glucose 6 phosphate dehydrogenase", "zwf"],
}


def search_metabolites(model, term: str, limit: int = 25) -> List[Tuple[str, str, str]]:
    term = term.lower()
    hits = []
    for m in model.metabolites:
        blob = f"{m.id} {m.name}".lower()
        if term in blob:
            hits.append((m.id, m.name, m.compartment))
    return hits[:limit]


def search_reactions(model, term: str, limit: int = 25) -> List[Tuple[str, str, str]]:
    term = term.lower()
    hits = []
    for r in model.reactions:
        blob = f"{r.id} {r.name}".lower()
        if term in blob:
            hits.append((r.id, r.name, r.reaction))
    return hits[:limit]


def search_genes(model, term: str, limit: int = 25) -> List[Tuple[str, str]]:
    term = term.lower()
    hits = []
    for g in model.genes:
        blob = f"{g.id} {g.name}".lower()
        if term in blob:
            hits.append((g.id, g.name))
    return hits[:limit]


def is_exchange_reaction(r) -> bool:
    rid = r.id
    name = (r.name or "").lower()
    if rid.startswith(("EX_", "Ex_", "DM_", "R_EX_", "R_DM_")):
        return True
    if "exchange" in name:
        return True
    if r in r.model.boundary:
        return True
    return False


def list_exchanges(model, limit: int = 80) -> List[Tuple[str, str]]:
    exs = [r for r in model.reactions if is_exchange_reaction(r)]
    exs.sort(key=lambda r: r.id)
    return [(r.id, r.name) for r in exs[:limit]]


def find_biomass(model) -> List[Tuple[str, str]]:
    cands = [
        r
        for r in model.reactions
        if "biomass" in r.id.lower()
        or "biomass" in (r.name or "").lower()
        or "growth" in r.id.lower()
    ]
    if cands:
        return [(r.id, r.name) for r in cands]
    return [("(current objective)", str(model.objective.expression))]


def resolve_gene_aliases(model) -> Dict[str, List[Dict[str, str]]]:
    """Map common gene names to reaction ids when the model has no GPR layer."""
    out: Dict[str, List[Dict[str, str]]] = {}
    for gene, terms in GENE_ALIASES.items():
        hits: List[Dict[str, str]] = []
        if gene in model.genes:
            hits.append({"type": "gene", "id": gene, "name": model.genes.get_by_id(gene).name or ""})
        for term in terms:
            for rid, name, _ in search_reactions(model, term, limit=8):
                entry = {"type": "reaction", "id": rid, "name": name, "matched_term": term}
                if entry not in hits:
                    hits.append(entry)
        if hits:
            out[gene] = hits[:12]
    return out


def pick_recommended_ids(model) -> Dict[str, Optional[str]]:
    """Best-effort defaults for the wax-ester / oleate demo on iYLI647."""
    rec: Dict[str, Optional[str]] = {
        "model_ref": "iYLI647.xml",
        "biomass_rxn": None,
        "carbon_source_rxn": None,
        "carbon_source_metabolite": None,
        "oleoyl_coa_metabolite": None,
        "nadph_metabolite": None,
        "nadp_metabolite": None,
        "coa_metabolite": None,
    }

    biomass = find_biomass(model)
    for rid, name in biomass:
        if "carbon" in (name or "").lower() or rid == IYLI647_BIOMASS_CARBON:
            rec["biomass_rxn"] = rid
            break
    if not rec["biomass_rxn"] and biomass:
        rec["biomass_rxn"] = biomass[0][0]

    oleate_ex = search_reactions(model, "octadecenoate exchange", limit=5)
    if oleate_ex:
        rec["carbon_source_rxn"] = oleate_ex[0][0]
    elif search_reactions(model, "ocdcea", limit=5):
        for rid, name, _ in search_reactions(model, "ocdcea", limit=20):
            if is_exchange_reaction(model.reactions.get_by_id(rid)):
                rec["carbon_source_rxn"] = rid
                break

    for term, key in [
        ("octadecenoate", "carbon_source_metabolite"),
        ("octadecenoyl_coa", "oleoyl_coa_metabolite"),
        ("nadph", "nadph_metabolite"),
        ("nadp", "nadp_metabolite"),
        ("coa", "coa_metabolite"),
    ]:
        hits = search_metabolites(model, term, limit=15)
        if not hits:
            continue
        if key == "coa_metabolite":
            exact = [h for h in hits if h[0].endswith("_c") and h[0] in ("coa_c", "M_coa_c")]
            cytoplasm = [h for h in hits if h[2] == "c" and h[0].endswith("coa_c") and "acc" not in h[0]]
            rec[key] = (exact or cytoplasm or hits)[0][0]
        elif key == "nadp_metabolite":
            exact = [h for h in hits if "nadp_c" in h[0] and "nadph" not in h[0]]
            rec[key] = (exact or hits)[0][0]
        elif key == "carbon_source_metabolite":
            extracellular = [h for h in hits if h[2] == "e"]
            rec[key] = (extracellular or hits)[0][0]
        elif key == "oleoyl_coa_metabolite":
            cytoplasm = [h for h in hits if h[2] == "c" and "odecoa" in h[0]]
            rec[key] = (cytoplasm or hits)[0][0]
        else:
            cytoplasm = [h for h in hits if h[2] == "c" and "nadph" in h[0]]
            rec[key] = (cytoplasm or hits)[0][0]

    return rec


def build_report(model_ref: str, extra_terms: List[str]) -> Dict[str, Any]:
    if is_sbml_path(model_ref):
        preflight = ensure_sbml_loadable(model_ref)
        if not preflight.get("loadable"):
            return {
                "status": "error",
                "model_ref": model_ref,
                "preflight": preflight,
                "message": preflight.get("error", "SBML not loadable"),
            }

    model, load_report = load_model_robust(model_ref)
    default_met_terms = [
        "oleate", "oleic", "octadec", "coa", "nadph", "palmit", "stear",
        "hexadec", "glycerol", "alcohol", "ester", "wax",
    ]
    default_gene_terms = ["pox", "mfe", "dga", "lro", "are", "faa", "lip", "ole", "zwf", "fas"]

    searches: Dict[str, Any] = {"metabolites": {}, "genes": {}, "reactions": {}}
    for term in default_met_terms + extra_terms:
        searches["metabolites"][term] = [
            {"id": a, "name": b, "compartment": c}
            for a, b, c in search_metabolites(model, term)
        ]
    for term in default_gene_terms + extra_terms:
        searches["genes"][term] = [
            {"id": a, "name": b} for a, b in search_genes(model, term)
        ]
    for term in ["oxidase", "acyltransferase", "lipase", "synthetase", "desaturase", "reductase"]:
        searches["reactions"][term] = [
            {"id": a, "name": b} for a, b, _ in search_reactions(model, term, limit=10)
        ]

    peroxisomal_pox = search_reactions(model, "acyl coa oxidase", limit=20)
    peroxisomal_pox = [h for h in peroxisomal_pox if "peroxisomal" in h[1].lower()]

    return {
        "status": "ok",
        "model_ref": model_ref,
        "load": load_report,
        "preflight": load_report.get("preflight"),
        "summary": {
            "model_id": model.id,
            "n_reactions": len(model.reactions),
            "n_metabolites": len(model.metabolites),
            "n_genes": len(model.genes),
            "has_gene_associations": len(model.genes) > 0,
            "id_prefix_note": (
                "COBRApy normalizes ids on load (typically drops M_/R_ prefixes from SBML). "
                "Always use ids from this report, not raw XML tags."
            ),
        },
        "recommended": pick_recommended_ids(model),
        "biomass_reactions": [{"id": a, "name": b} for a, b in find_biomass(model)],
        "exchange_reactions": [{"id": a, "name": b} for a, b in list_exchanges(model)],
        "gene_alias_resolution": resolve_gene_aliases(model),
        "peroxisomal_acyl_coa_oxidases": [
            {"id": a, "name": b} for a, b, _ in peroxisomal_pox
        ],
        "searches": searches,
        "agent_rules": [
            "Run find_ids.py before every score_pathway() call on a new model.",
            "Use reaction ids for knockouts when has_gene_associations is false.",
            "Never invent metabolite ids; copy from recommended or searches.",
            "If predicted_product_flux is 0 with status optimal, re-check ids here.",
        ],
    }


def print_human_report(report: Dict[str, Any]) -> None:
    if report.get("status") == "error":
        print("PREFLIGHT FAILED")
        print(json.dumps(report, indent=2))
        return

    s = report["summary"]
    print(f"Loaded: {s['model_id']} | {s['n_reactions']} rxns, "
          f"{s['n_metabolites']} mets, {s['n_genes']} genes")
    if not s["has_gene_associations"]:
        print("\n*** NO GENE ASSOCIATIONS in this model — use REACTION ids for knockouts ***")

    print("\n" + "=" * 72)
    print("RECOMMENDED IDs (copy into score_pathway)")
    print("=" * 72)
    for k, v in report["recommended"].items():
        print(f"  {k:28s}  {v}")

    print("\n" + "=" * 72)
    print("BIOMASS / GROWTH reactions")
    print("=" * 72)
    for row in report["biomass_reactions"]:
        print(f"  {row['id']:40s}  {row['name']}")

    print("\n" + "=" * 72)
    print("EXCHANGE / BOUNDARY reactions (first 40)")
    print("=" * 72)
    for row in report["exchange_reactions"][:40]:
        print(f"  {row['id']:40s}  {row['name']}")

    print("\n" + "=" * 72)
    print("GENE ALIAS -> REACTION resolution (for POX1-6, DGA1, LIP2, ...)")
    print("=" * 72)
    for gene, hits in report["gene_alias_resolution"].items():
        print(f"\n  [{gene}]")
        for h in hits[:6]:
            print(f"      {h['type']:10s}  {h['id']:32s}  {h.get('name', '')}")

    print("\n" + "=" * 72)
    print("PEROXISOMAL acyl-CoA oxidases (typical POX knockout targets)")
    print("=" * 72)
    for row in report["peroxisomal_acyl_coa_oxidases"]:
        print(f"  {row['id']:40s}  {row['name']}")

    print("\n" + "=" * 72)
    print("METABOLITE search highlights")
    print("=" * 72)
    for term, hits in report["searches"]["metabolites"].items():
        if not hits:
            continue
        print(f"\n  [{term}] -> {len(hits)} match(es)")
        for row in hits[:5]:
            print(f"      {row['id']:24s} ({row['compartment']})  {row['name']}")

    if report.get("preflight", {}).get("repairs_applied"):
        print("\n" + "=" * 72)
        print("SBML REPAIRS APPLIED (saved to disk)")
        print("=" * 72)
        for line in report["preflight"]["repairs_applied"]:
            print(f"  - {line}")

    print("\nDONE. Use --json for machine-readable output for the orchestrator agent.")


def main() -> None:
    ap = argparse.ArgumentParser(description="Resolve model IDs for score_pathway()")
    ap.add_argument("model", help="path to SBML model or bundled name (e.g. iYLI647.xml)")
    ap.add_argument("--term", action="append", default=[], help="extra search term; repeatable")
    ap.add_argument("--json", action="store_true", help="emit JSON for upstream agents")
    ap.add_argument("--preflight-only", action="store_true",
                    help="validate/repair SBML only; do not search ids")
    args = ap.parse_args()

    if args.preflight_only:
        if not is_sbml_path(args.model):
            print(json.dumps({"status": "skipped", "reason": "not an SBML path"}))
            return
        result = ensure_sbml_loadable(args.model)
        print(json.dumps(result, indent=2))
        sys.exit(0 if result.get("loadable") else 1)

    if not args.json:
        print(f"Loading model: {args.model} ...")
    report = build_report(args.model, args.term)

    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print_human_report(report)

    sys.exit(0 if report.get("status") == "ok" else 1)


if __name__ == "__main__":
    main()
