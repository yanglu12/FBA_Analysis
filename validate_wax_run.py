"""Literature sanity check for wax-ester FBA on iYLI647."""
from __future__ import annotations

from fba_tool import _add_demand, _apply_edits, _insert_reactions, score_pathway
from model_loader import load_model_robust

MODEL = "c:/Users/Administrator/Documents/Coding/FBA_analysis/iYLI647.xml"

FAR = {
    "id": "FAR",
    "stoichiometry": {
        "odecoa_c": -1,
        "nadph_c": -2,
        "oleyl_alcohol_c": 1,
        "coa_c": 1,
        "nadp_c": 2,
    },
}
WS = {
    "id": "WS",
    "stoichiometry": {
        "oleyl_alcohol_c": -1,
        "odecoa_c": -1,
        "wax_ester_c": 1,
        "coa_c": 1,
    },
}


def main() -> None:
    print("=== iYLI647 literature sanity checks ===\n")

    model, _ = load_model_robust(MODEL)

    # Paper validates on glucose/glycerol at 10 mmol/gDW/h (Mishra et al. 2018)
    for bio in ["BIOMASS_yarrowia_carbon_limiting", "BIOMASS_yarrowia_nitrogen_limiting"]:
        m = model.copy()
        m.objective = bio
        m.reactions.get_by_id("EX_glc_LPAREN_e_RPAREN_").lower_bound = -10
        sol = m.optimize()
        print(f"Glucose 10 mmol/gDW/h | {bio}")
        print(f"  growth rate (1/h): {sol.fluxes[bio]:.4f}")
        print(f"  glucose uptake:    {sol.fluxes['EX_glc_LPAREN_e_RPAREN_']:.4f}")

    print()
    for bio in ["BIOMASS_yarrowia_carbon_limiting", "BIOMASS_yarrowia_nitrogen_limiting"]:
        m = model.copy()
        m.objective = bio
        m.reactions.get_by_id("EX_ocdcea_LPAREN_e_RPAREN_").lower_bound = -10
        sol = m.optimize()
        print(f"Oleate 10 mmol/gDW/h | {bio}")
        print(f"  growth rate (1/h): {sol.fluxes[bio]:.4f}")
        print(f"  oleate uptake:     {sol.fluxes['EX_ocdcea_LPAREN_e_RPAREN_']:.4f}")

    print("\n=== Current wax demo (contract example) ===")
    out = score_pathway(
        MODEL,
        carbon_source_rxn="EX_ocdcea_LPAREN_e_RPAREN_",
        carbon_source_uptake=10.0,
        candidate_reactions=[FAR, WS],
        product_metabolite="wax_ester_c",
        knockouts=["ACOAO4p", "ACOAO5p", "ACOAO7p", "ACOAO8p", "ACOAO9p"],
        biomass_rxn="BIOMASS_yarrowia_carbon_limiting",
        objective="product",
        min_growth_fraction=0.1,
    )
    print(f"  product flux:  {out['predicted_product_flux']}")
    print(f"  growth rate:   {out['growth_rate']}")
    print(f"  yield mol/mol: {out['yield_mol_per_mol_substrate']}")
    print(f"  theoretical max yield (2 C18 -> 1 WE): 0.500")

    # Carbon audit
    m = model.copy()
    _insert_reactions(m, [FAR, WS])
    _add_demand(m, "wax_ester_c")
    _apply_edits(m, ["ACOAO4p", "ACOAO5p", "ACOAO7p", "ACOAO8p", "ACOAO9p"], {})
    m.reactions.get_by_id("EX_ocdcea_LPAREN_e_RPAREN_").lower_bound = -10
    m.reactions.get_by_id("BIOMASS_yarrowia_carbon_limiting").lower_bound = 0.01
    m.objective = "DM_wax_ester_c"
    sol = m.optimize()
    oleate_in = abs(sol.fluxes["EX_ocdcea_LPAREN_e_RPAREN_"])
    facoal = sol.fluxes["FACOAL181"]
    far = sol.fluxes["FAR"]
    ws = sol.fluxes["WS"]
    odecoa = m.metabolites.get_by_id("odecoa_c")
    odecoa_prod = sum(
        sol.fluxes[r.id] * r.metabolites[odecoa]
        for r in m.reactions
        if odecoa in r.metabolites and sol.fluxes[r.id] * r.metabolites[odecoa] > 1e-6
    )
    odecoa_cons = sum(
        sol.fluxes[r.id] * r.metabolites[odecoa]
        for r in m.reactions
        if odecoa in r.metabolites and sol.fluxes[r.id] * r.metabolites[odecoa] < -1e-6
    )
    print("\n=== Carbon / CoA audit ===")
    print(f"  oleate import (FACOAL181): {facoal:.4f}")
    print(f"  odecoa net production:    {odecoa_prod:.4f}")
    print(f"  odecoa net consumption:   {odecoa_cons:.4f}")
    print(f"  FAR + WS flux:              {far:.4f} / {ws:.4f}")
    print(f"  2*WS flux = odecoa for WE:  {2*ws:.4f} (needs <= oleate_in {oleate_in:.4f})")


if __name__ == "__main__":
    main()
