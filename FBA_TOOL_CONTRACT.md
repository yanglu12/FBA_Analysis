# FBA Tool Contract — `score_pathway`

**Owner:** Yang (metabolic modeling)  
**For:** the agent orchestration layer (Mahbub)  
**Purpose:** Score a candidate strain design with flux balance analysis. The agent
calls this once per candidate design and reads the result to rank pathways and
decide what to try next.

---

## 0. Mandatory agent workflow (do not skip)

```
1. python find_ids.py iYLI647.xml --json          # resolve real ids + preflight SBML
2. read recommended.* and gene_alias_resolution   # build score_pathway payload
3. score_pathway(...)                             # run FBA
4. interpret status / flux / bottlenecks / calibration  # rank or revise design
```

**Never call `score_pathway()` on a new SBML model without step 1.**

If step 1 reports `has_gene_associations: false`, use **reaction ids** for
knockouts (not yeast gene names like `POX1`). Even when genes exist, iYLI647
gene ids may not match common names like `POX1` — use `gene_alias_resolution`
from find_ids output.

If step 1 reports SBML repairs, the file was auto-fixed on disk (missing
reaction/model `id` attributes). Re-run step 1 to confirm `loadable: true`.

---

## 1. What this tool does (one sentence for the LLM)

> Given a genome-scale model, a feedstock, a set of heterologous reactions to
> insert, and a set of gene edits, run flux balance analysis and return the
> predicted product flux, growth rate, yield, and the limiting (bottleneck)
> reactions.

Use this when you need a **quantitative** prediction of how much product a
proposed strain design can make — not for literature lookup, not for choosing
enzymes (that happens upstream). This tool evaluates a design that has already
been proposed.

---

## 2. CRITICAL RULE — IDs must be real

Every metabolite and **reaction** ID you pass **must already exist in the model**
(except product/intermediate metabolites you deliberately create in
`candidate_reactions`). If you invent IDs, the tool runs but returns a
meaningless flux of 0.

Resolve real IDs with:

```bash
python find_ids.py iYLI647.xml --json
```

Copy values from `recommended`, `exchange_reactions`, `gene_alias_resolution`,
and `peroxisomal_acyl_coa_oxidases`. **Do not guess.**

### iYLI647 resolved ids (from find_ids.py — use these)

| Role | ID (after COBRApy load) |
|---|---|
| Model file | `iYLI647.xml` (same SBML as `12918_2018_542_MOESM2_ESM.xml`) |
| Biomass (carbon-limited) | `BIOMASS_yarrowia_carbon_limiting` |
| Biomass (nitrogen-limited) | `BIOMASS_yarrowia_nitrogen_limiting` |
| Oleate exchange (feedstock) | `EX_ocdcea_LPAREN_e_RPAREN_` |
| Glucose exchange | `EX_glc_LPAREN_e_RPAREN_` |
| Oleate (extracellular) | `ocdcea_e` |
| Oleoyl-CoA (cytoplasm) | `odecoa_c` |
| NADPH / NADP+ | `nadph_c` / `nadp_c` |
| CoA (cytoplasm) | `coa_c` |

COBRApy **strips** `M_` / `R_` prefixes from the raw SBML. Always copy ids
from `find_ids.py --json` → `recommended`, never from the XML file directly.

### Knockouts on iYLI647 — map POX1–6 to reaction ids

Gene names like `POX1` are **not** in this model. Peroxisomal β-oxidation is
represented by **acyl-CoA oxidase reactions**, e.g.:

- `ACOAO8p` — octadecenoyl-CoA (peroxisomal)
- `ACOAO4p`, `ACOAO5p`, `ACOAO7p`, `ACOAO9p` — other chain lengths

Use `gene_alias_resolution` and `peroxisomal_acyl_coa_oxidases` from
`find_ids.py --json`.

---

## 3. Parameters

| Parameter | Type | Required | Meaning |
|---|---|---|---|
| `model_ref` | string | no (default `iYLI647.xml`) | Path to SBML model or bundled name. |
| `carbon_source_rxn` | string | **yes for a real run** | Exchange reaction for feedstock. From `find_ids`. |
| `carbon_source_uptake` | number | no (default 10.0) | Max uptake mmol/gDW/h. Applied as `lower_bound = -uptake`. |
| `candidate_reactions` | list | **yes** for heterologous pathways | Reactions to insert (FAR, WS, …). See schema. |
| `product_metabolite` | string | **yes** | Target product metabolite id (wax ester — usually new). |
| `product_demand_rxn` | string | no | Existing product-exit reaction; overrides auto demand. |
| `knockouts` | list | no | **Reaction ids** (preferred) or gene ids if model has GPR. |
| `bound_overrides` | object | no | `{reaction_id: [lower, upper]}` for OE / attenuation. |
| `objective` | string | no (default `"product"`) | `"product"`, `"biomass"`, or `"coupled"` (same as product today). |
| `min_growth_fraction` | number | no (default 0.1) | Min growth as fraction of max when maximizing product. |
| `min_growth_rate` | number | no | Absolute growth floor (1/h); overrides fraction when set. |
| `max_growth_rate` | number | no | Absolute growth cap (1/h); use experimental μ from papers. |
| `biomass_rxn` | string | no | Biomass reaction id. Default: auto-detect. |
| `exchange_constraints` | object | no | `{exchange_id: [lb, ub]}` for CER, O₂, NH₄⁺, co-substrates. |
| `use_minimal_medium` | bool | no (default false) | Close all exchanges except feedstock + constraints. |
| `substrate_moles_per_product` | number | no (default 1.0) | Biochem ratio for yield correction (e.g. 2.0 for C36 wax). |
| `scenario` | string or dict | no | Path to YAML under `scenarios/` or pre-loaded dict. |
| `simulation_notes` | list | no | Free-text assumptions echoed in result. |

### `candidate_reactions` object schema

```json
{
  "id": "FAR",
  "name": "fatty acyl-CoA reductase",
  "stoichiometry": {
    "odecoa_c": -1,
    "nadph_c": -2,
    "oleyl_alcohol_c": 1,
    "coa_c": 1,
    "nadp_c": 2
  },
  "lower_bound": 0,
  "upper_bound": 1000,
  "subsystem": "user"
}
```

---

## 4. Wax-ester example call (iYLI647, oleate feedstock)

```python
from fba_tool import score_pathway

FAR = {
    "id": "FAR",
    "stoichiometry": {
        "odecoa_c": -1, "nadph_c": -2,
        "oleyl_alcohol_c": 1, "coa_c": 1, "nadp_c": 2,
    },
}
WS = {
    "id": "WS",
    "stoichiometry": {
        "oleyl_alcohol_c": -1, "odecoa_c": -1,
        "wax_ester_c": 1, "coa_c": 1,
    },
}

out = score_pathway(
    "iYLI647.xml",
    scenario="scenarios/wax_ester_oleate_n_limited.yaml",
    candidate_reactions=[FAR, WS],
    product_metabolite="wax_ester_c",
    objective="product",
)
print(out["predicted_product_flux"])
print(out["calibration"]["confidence_level"])
```

---

## 5. Return schema (what the agent reads)

```json
{
  "status": "optimal",
  "objective_used": "product",
  "predicted_product_flux": 5.746,
  "growth_rate": 0.01,
  "yield_mol_per_mol_substrate": 0.5746,
  "yield_corrected_mol_per_mol_substrate": 0.2873,
  "simulation_context": { "...": "..." },
  "calibration": {
    "confidence_level": "literature_calibrated",
    "recommended_use": "quantitative_comparison_after_biomass_validation",
    "missing_literature_inputs": [],
    "warnings": [],
    "literature_refs": ["..."],
    "agent_guidance": "Do not report predicted_product_flux as g/L ..."
  },
  "inserted_reactions": ["FAR", "WS"],
  "edits_applied": { "knocked_out": ["ACOAO4p", "..."], "not_found": [] },
  "bottlenecks": [ { "reaction": "...", "flux": 1.2, "span": 0.0, "at_bound": true } ],
  "carbon_audit": {
    "feedstock_is_sole_carbon_source": true,
    "side_door_carbon_imports": [],
    "total_carbon_import_mmol_per_h": 59.5084,
    "carbon_imports": [ { "exchange": "EX_ocdcea_...", "carbon_mmol_per_h": 59.5084 } ]
  },
  "message": "OK"
}
```

### Field guide for ranking designs

| Field | Use |
|---|---|
| `status` | Must be `optimal` to compare fluxes. |
| `predicted_product_flux` | mmol/gDW/h — **not** g/L. |
| `yield_corrected_mol_per_mol_substrate` | Prefer over raw yield for multi-carbon products. |
| `growth_rate` | Check against experimental μ; unrealistic μ → tighten scenario. |
| `bottlenecks` | Reactions to OE or unblock in next design iteration. |
| `calibration.confidence_level` | Gate how strongly to cite the number. |
| `carbon_audit` | Boundary carbon imports/exports; check `feedstock_is_sole_carbon_source`. |

---

## 6. Calibration levels (mandatory read)

Every result includes `calibration`. **Do not skip it.**

Reported fields:

- `confidence_level` — conservative overall label
- `medium_confidence_level` — uptake / exchange / growth caps
- `product_confidence_level` — `not_applicable` | `unvalidated` | `literature_calibrated`

| `confidence_level` | Meaning | Agent may |
|---|---|---|
| `exploratory` | ≥3 missing medium/growth inputs | Rank designs relative to each other only |
| `partial` | 1–2 missing medium inputs | Rank + bottlenecks; in silico upper bound |
| `medium_calibrated` | Medium/growth pinned; no `product_literature_refs` | Trust μ/medium; **not** product titer |
| `literature_calibrated` | Biomass fully pinned, or product + product literature | Quantitative comparison after biomass validation |
| `invalid` | Infeasible or error | Debug constraints only |

Medium inputs checked: `literature_refs`, `exchange_constraints`, `max_growth_rate`.  
Product runs also require `product_literature_refs` for `product_confidence_level: literature_calibrated`.

---

## 7. Biomass validation (run before product claims)

Before trusting wax-ester product flux, run biomass validation:

```python
out = score_pathway(
    "iYLI647.xml",
    scenario="scenarios/iyli647_glucose_validation.yaml",
    objective="biomass",
)
# Compare out["growth_rate"] to expected_experimental_growth_rate in scenario (0.24 h-1)
```

If biomass validation fails (`status != optimal`), fix medium constraints before
product runs.

---

## 8. Simulation context — what comes from literature vs. code

Most “accuracy” knobs are **not hardcoded in Python**. They should be filled from
literature review (or your experimental dataset) and passed as parameters or a
**scenario YAML** under `scenarios/`.

| Parameter | Typical literature source | Applies to |
|---|---|---|
| `biomass_rxn` | Model paper methods (C- vs N-limited biomass) | Any organism with dual biomass equations |
| `carbon_source_rxn`, `carbon_source_uptake` | Medium composition + measured uptake (mmol/gDW/h) | Any feedstock |
| `exchange_constraints` | Validation tables: CER, O₂, NH₄⁺, co-substrates | Any GSMM validation paper |
| `max_growth_rate` / `min_growth_rate` | Batch μ from growth curves | Any production phase |
| `use_minimal_medium` | Methods section (“minimal medium”, closed exchanges) | Most validation papers |
| `substrate_moles_per_product` | Pathway stoichiometry | Multi-substrate products (wax, polymers) |
| `product_literature_refs` | Experimental product flux/titer papers | Required for product confidence |
| `knockouts`, `bound_overrides` | Strain design papers | Any engineered strain |

### Bundled scenario files

| File | Purpose |
|---|---|
| `scenarios/wax_ester_oleate_n_limited.yaml` | **Default** N-limited production on oleate |
| `scenarios/wax_ester_oleate_open.yaml` | Optimistic upper bound (partial calibration) |
| `scenarios/iyli647_glucose_validation.yaml` | Biomass sanity check vs Mishra 2018 / Workman 2013 |
| `scenarios/exploratory_no_literature.yaml` | Novel pathways — exploratory calibration only |

When `use_minimal_medium: true`, you **must** list essential cofeed exchanges
(O₂, NH₄⁺, Pi, SO₄, ions, H₂O, H⁺) in `exchange_constraints` or the model
will be infeasible.

**Product runs:** always inspect `carbon_audit` in the result. If
`feedstock_is_sole_carbon_source` is false, non-feedstock exchanges (often
`EX_glc_LPAREN_e_RPAREN_`) are importing carbon and yield vs the named
feedstock is not trustworthy until those exchanges are closed.

---

## 9. Accuracy expectations

- **Open medium allows glucose and other carbon side doors** — use
  `use_minimal_medium: true` for oleate-only runs; read `carbon_audit`.
- **Raw FBA on open medium over-predicts** wax yield when `EX_glc` etc. import carbon.
- **iYLI647 was validated on glucose/glycerol**, not wax esters — no direct
  literature benchmark exists for wax product flux on this model.
- **`literature_calibrated`** still means mmol/gDW/h, not g/L until you add
  bioreactor mass-balance conversion.
- Always report **`yield_corrected_mol_per_mol_substrate`** for wax (set
  `substrate_moles_per_product: 2.0` for C36 from 2×C18).

---

## 10. Porting to another organism

1. Obtain GSMM SBML + validation paper for that organism.
2. Run `find_ids.py <model.xml> --json`.
3. Create a scenario YAML under `scenarios/`.
4. Fill literature values (uptake, CER, μ, biomass equation, knockouts).
5. Run biomass validation scenario first (`objective="biomass"`).
6. Only then run product objectives and read `calibration`.

---

## 11. Trial run

```bash
conda activate fba
python run_all_tests.py
```

See `EXPERT_REPORT.md` for the latest summary table and team-facing interpretation.
