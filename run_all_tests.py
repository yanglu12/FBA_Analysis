"""Run all FBA_analysis checks and print an agent-readable summary."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from fba_tool import score_pathway

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
WAX_ARGS = dict(
    candidate_reactions=[FAR, WS],
    product_metabolite="wax_ester_c",
    objective="product",
)


def row(name: str, out: dict) -> dict:
    cal = out.get("calibration") or {}
    return {
        "test": name,
        "status": out.get("status"),
        "product": out.get("predicted_product_flux"),
        "growth": out.get("growth_rate"),
        "yield_raw": out.get("yield_mol_per_mol_substrate"),
        "yield_corr": out.get("yield_corrected_mol_per_mol_substrate"),
        "confidence": cal.get("confidence_level"),
        "medium_conf": cal.get("medium_confidence_level"),
        "product_conf": cal.get("product_confidence_level"),
        "use_for": cal.get("recommended_use"),
        "warnings": len(cal.get("warnings") or []),
        "missing_lit": len(cal.get("missing_literature_inputs") or []),
    }


def main() -> None:
    model = str(ROOT / "iYLI647.xml")
    results = []

    print("=" * 72)
    print("1. find_ids preflight + JSON")
    print("=" * 72)
    r = subprocess.run(
        [sys.executable, str(ROOT / "find_ids.py"), model, "--json"],
        capture_output=True,
        text=True,
    )
    if r.returncode == 0:
        data = json.loads(r.stdout)
        print(f"  find_ids status: {data.get('status')}")
        print(f"  recommended biomass: {data.get('recommended', {}).get('biomass_rxn')}")
    else:
        print("  find_ids FAILED", r.stderr[-500:])

    print("\n" + "=" * 72)
    print("2. score_pathway scenarios")
    print("=" * 72)

    scenarios = [
        ("wax_open_scenario", ROOT / "scenarios/wax_ester_oleate_open.yaml", "product"),
        ("wax_n_limited_scenario", ROOT / "scenarios/wax_ester_oleate_n_limited.yaml", "product"),
        ("glucose_validation_biomass", ROOT / "scenarios/iyli647_glucose_validation.yaml", "biomass"),
        ("exploratory_no_literature", ROOT / "scenarios/exploratory_no_literature.yaml", "product"),
    ]

    for name, scen, obj in scenarios:
        if name == "glucose_validation_biomass":
            out = score_pathway(
                model,
                scenario=str(scen),
                objective="biomass",
            )
        else:
            out = score_pathway(model, scenario=str(scen), **WAX_ARGS)
        results.append(row(name, out))

    print("\n" + "=" * 72)
    print("3. textbook self-test (bundled model)")
    print("=" * 72)
    r = subprocess.run([sys.executable, str(ROOT / "fba_tool.py")], capture_output=True, text=True)
    self_ok = r.returncode == 0 and "optimal" in r.stdout
    print(f"  self-test exit={r.returncode} optimal_in_output={self_ok}")

    print("\n" + "=" * 72)
    print("SUMMARY TABLE")
    print("=" * 72)
    for r in results:
        print(
            f"  {r['test']:28s} | {str(r['status']):10s} | "
            f"prod={str(r['product']):6s} | growth={str(r['growth']):6s} | "
            f"ycorr={str(r['yield_corr']):6s} | conf={r['confidence']}"
        )
        if r.get("medium_conf") or r.get("product_conf"):
            print(
                f"    └ medium={r['medium_conf']} product={r['product_conf']}"
            )
        if r["warnings"] or r["missing_lit"]:
            print(f"    └ missing_lit={r['missing_lit']} warnings={r['warnings']}")

    print("\n" + "=" * 72)
    print("CALIBRATION DETAIL (exploratory_no_literature — agent should read this)")
    print("=" * 72)
    exploratory = score_pathway(
        model,
        scenario=str(ROOT / "scenarios/exploratory_no_literature.yaml"),
        **WAX_ARGS,
    )
    print(json.dumps(exploratory.get("calibration"), indent=2))


if __name__ == "__main__":
    main()
