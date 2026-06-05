# FBA_Analysis

Agent-callable **Flux Balance Analysis (FBA)** tooling for scoring strain designs
on genome-scale metabolic models (GSMM). The reference implementation targets
*iYLI647* (*Yarrowia lipolytica*) and wax-ester production on oleate, but the
same workflow applies to other SBML models, substrates, and products.

## Quick start

```bash
# Create environment (once)
conda env create -f environment.yml
conda activate fba

# Clone and enter repo
cd FBA_analysis

# 1. Resolve model IDs + validate SBML
python find_ids.py iYLI647.xml --json

# 2. Run trial suite
python run_all_tests.py

# 3. Run wax-ester verification / carbon audit
python testing_debugging/run_all_debug_checks.py
```

## Agent workflow (mandatory)

```
find_ids.py --json  →  build score_pathway payload  →  read calibration  →  rank designs
```

Full API spec: [`FBA_TOOL_CONTRACT.md`](FBA_TOOL_CONTRACT.md)

## Key entry points

| File | Purpose |
|---|---|
| [`fba_tool.py`](fba_tool.py) | `score_pathway()` — main FBA scoring function |
| [`find_ids.py`](find_ids.py) | Resolve reaction/metabolite IDs; SBML preflight |
| [`model_loader.py`](model_loader.py) | Robust SBML load + auto-repair |
| [`scenarios/`](scenarios/) | Literature-backed simulation YAML templates |
| [`run_all_tests.py`](run_all_tests.py) | One-command integration trial |
| [`testing_debugging/`](testing_debugging/) | Pathway, stoichiometry, carbon verification |

## Scenarios

| Scenario | Use |
|---|---|
| `wax_ester_oleate_n_limited.yaml` | **Default** production-mode oleate → wax |
| `wax_ester_oleate_open.yaml` | Optimistic upper bound (partial calibration) |
| `iyli647_glucose_validation.yaml` | Biomass sanity check vs Workman/Mishra |
| `exploratory_no_literature.yaml` | Novel pathways — ranking only |

## Reading results

Every `score_pathway()` result includes a **`calibration`** block:

| `confidence_level` | Meaning |
|---|---|
| `exploratory` | Missing most literature inputs — relative ranking only |
| `partial` | Some medium/growth pins missing |
| `medium_calibrated` | Medium/growth pinned; **product flux not experimentally validated** |
| `literature_calibrated` | Biomass validation or product + medium fully anchored |
| `invalid` | Infeasible / error — debug constraints |

Also check `medium_confidence_level` and `product_confidence_level` separately.

**Do not** report `predicted_product_flux` as g/L or guaranteed titer unless
`product_confidence_level` is `literature_calibrated`.

Inspect **`carbon_audit.feedstock_is_sole_carbon_source`** before citing yield.
Open medium on iYLI647 typically imports glucose carbon via `EX_glc_LPAREN_e_RPAREN_`.

Prefer **`yield_corrected_mol_per_mol_substrate`** over raw yield for multi-carbon
products (e.g. set `substrate_moles_per_product: 2.0` for C36 wax from 2×C18).

## Porting to another organism

1. Add your SBML model to the repo (or pass path to `model_ref`).
2. `python find_ids.py your_model.xml --json`
3. Create a scenario YAML under `scenarios/` with uptake, CER, growth caps from literature.
4. Run biomass validation (`objective="biomass"`) before product objectives.
5. Add `product_literature_refs` when you have experimental product flux/titer data.

## Documentation

- [`FBA_TOOL_CONTRACT.md`](FBA_TOOL_CONTRACT.md) — agent API contract
- [`EXPERT_REPORT.md`](EXPERT_REPORT.md) — team-facing accuracy review
- [`testing_debugging/README.md`](testing_debugging/README.md) — verification findings

## Model files

- `iYLI647.xml` — canonical patched SBML (Mishra et al. 2018)
- `12918_2018_542_MOESM2_ESM.xml` — original supplementary SBML (patched in place)

## Example (Python)

```python
from fba_tool import score_pathway

FAR = {"id": "FAR", "stoichiometry": {"odecoa_c": -1, "nadph_c": -2, "oleyl_alcohol_c": 1, "coa_c": 1, "nadp_c": 2}}
WS = {"id": "WS", "stoichiometry": {"oleyl_alcohol_c": -1, "odecoa_c": -1, "wax_ester_c": 1, "coa_c": 1}}

out = score_pathway(
    "iYLI647.xml",
    scenario="scenarios/wax_ester_oleate_n_limited.yaml",
    candidate_reactions=[FAR, WS],
    product_metabolite="wax_ester_c",
    objective="product",
)
print(out["predicted_product_flux"], out["calibration"]["confidence_level"])
```

## License / citation

Model: Mishra et al. 2018 BMC Systems Biology ([iYLI647](https://doi.org/10.1186/s12918-018-0542-5)).
