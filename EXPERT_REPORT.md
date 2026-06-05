# FBA Agent Tool — Expert Review Report

**Project:** Yarrowia lipolytica wax-ester strain design scoring (iYLI647 GSMM)  
**Date:** 4 June 2026  
**Audience:** Metabolic modeling, strain engineering, and agent-orchestration leads  
**Repository:** `FBA_analysis` (`find_ids.py`, `fba_tool.py` → `score_pathway()`)

---

## 1. Executive summary

We built an **agent-callable Flux Balance Analysis (FBA) scoring tool** so an upstream LLM agent can quantitatively rank wax-ester strain designs on the published **iYLI647** genome-scale model (Mishra et al. 2018).

**Will the agent now give more accurate, requirement-fitting results?**

**Yes — when the agent follows the mandated workflow and uses calibrated scenarios.**  
**No — if the agent skips preflight, invents reaction IDs, or treats exploratory FBA flux as experimental titer.**

The key improvement is not “smarter FBA math” but **structured guardrails**:

1. **`find_ids.py --json`** resolves real COBRApy reaction/metabolite IDs and validates SBML before any run.
2. **Scenario YAML files** encode literature-backed medium constraints (CER, O₂, NH₄⁺, growth caps).
3. **`calibration` block** in every `score_pathway()` result tells the agent **how much to trust** the numbers (`exploratory` → `partial` → `literature_calibrated` → `invalid`).
4. **`yield_corrected_mol_per_mol_substrate`** applies biochemistry (2 mol C18 oleate → 1 mol C36 wax) so yields are not misread.

Without these, the same model **over-predicts** product flux and growth (open medium, internal lipid pools, unconstrained μ).

---

## 2. What the tool does (and does not do)

### Does

| Capability | Purpose |
|---|---|
| Load & repair iYLI647 SBML | Auto-fix missing biomass/model IDs; defensive error messages |
| Insert heterologous reactions (FAR + WS) | Score wax-ester pathway on oleate feedstock |
| Apply knockouts / bound overrides | Reaction-level edits (e.g. peroxisomal β-oxidation blocks) |
| Constrain medium & growth | `exchange_constraints`, `max_growth_rate`, `use_minimal_medium` |
| Return bottlenecks (FVA heuristic) | Identify flux-limited reactions for design iteration |
| Emit `calibration` metadata | Gate agent language (“rank only” vs “quantitative comparison”) |

### Does not

| Limitation | Implication |
|---|---|
| Predict **g/L titers** directly | Output is mmol/gDW/h; no bioreactor mass-balance to concentration |
| Validate **wax ester** against iYLI647 paper | Mishra 2018 validates glucose/glycerol + CER, not wax esters |
| Replace experimental strain testing | FBA is an *in silico* upper bound / ranking aid |
| Map yeast gene names (POX1) automatically | Model uses reaction IDs (`ACOAO4p`–`ACOAO9p`); agent must use `find_ids` |
| Guarantee global optimum on novel pathways | Heterologous reactions are user-specified; stoichiometry must be correct |

---

## 3. Mandatory agent workflow

```
1. python find_ids.py iYLI647.xml --json     # resolve IDs + SBML preflight
2. Build score_pathway payload from recommended.* + scenario YAML
3. score_pathway(...)                         # run FBA
4. Read status, flux, bottlenecks, calibration # rank or revise design
```

**Never** call `score_pathway()` on a new model without step 1.  
**Never** report `predicted_product_flux` as guaranteed titer unless
`product_confidence_level` is `literature_calibrated` **and** biomass validation
has passed. `medium_calibrated` pins growth/medium only — not product yield.

Full API spec: `FBA_TOOL_CONTRACT.md`.

---

## 4. Trial run results (4 June 2026)

Command: `python run_all_tests.py` (conda env `fba`, model `iYLI647.xml`)

| Scenario | FBA status | Product flux (mmol/gDW/h) | Growth μ (1/h) | Yield corrected (mol/mol) | Calibration |
|---|---:|---:|---:|---:|---|
| `wax_ester_oleate_open` | optimal | 5.36 | 0.94 | 0.27 | **partial** |
| `wax_ester_oleate_n_limited` | optimal | ~1.12 | 0.01 | ~0.06–0.17 | **medium_calibrated** |
| `iyli647_glucose_validation` (biomass) | optimal | — | **0.24** | — | **literature_calibrated** |
| `exploratory_no_literature` | optimal | 5.36 | 0.94 | 0.27 | **exploratory** |

Also passed: `find_ids.py --json` (status ok), `fba_tool.py` textbook self-test (optimal).

### Interpretation

**Glucose biomass validation (`iyli647_glucose_validation.yaml`)**  
- Reproduces Workman et al. 2013 batch glucose **μ = 0.24 h⁻¹** on minimal medium with O₂/NH₄⁺/Pi/SO₄/K⁺/Na⁺ cofeeds.  
- Confirms the tool + model load correctly before wax-ester product runs.  
- iYLI647 was validated on glucose/glycerol in the original paper; this scenario aligns with that literature.

**N-limited wax scenario (`wax_ester_oleate_n_limited.yaml`)**  
- Uses **`use_minimal_medium: true`** so oleate is the sole carbon import (closes `EX_glc` side door).  
- Open-medium runs inflated product flux (~5.7 mmol/gDW/h) because **glucose exchange imported 60 mmol C/h** alongside oleate — not internal pool export at steady state.  
- With closed carbon exchanges, product flux ~**1.1 mmol/gDW/h**; wax C ≤ total boundary carbon import.  
- **`medium_calibrated`** — medium/growth pinned; **`product_confidence_level: unvalidated`**.

**Verification (`testing_debugging/`)**  
- **Pathway attribution PASS:** WS/FAR knockout → product flux 0.  
- **Stoichiometry PASS:** FAR uses 2 NADPH per flux; FAR+WS consume 2× oleoyl-CoA per wax.  
- **Carbon audit:** Side-door `EX_glc_LPAREN_e_RPAREN_` explains >100% yield vs oleate alone on open medium; fixed by minimal medium.

**Open / exploratory wax scenarios**  
- Same raw product flux (~5.4 mmol/gDW/h) but **μ ≈ 0.94 h⁻¹** (unrealistic for N-limited production).  
- Raw yield >0.5 mol/mol before correction; corrected yield ~0.27.  
- Calibration correctly flags **partial** or **exploratory** — suitable **only for relative ranking** of knockouts/pathway variants, not for titer claims.

---

## 5. Accuracy assessment: before vs after calibration layer

| Issue (uncalibrated agent call) | Root cause | Mitigation in current tool |
|---|---|---|
| Product flux ~5.4 on open medium | `EX_glc` etc. import carbon when medium open | `use_minimal_medium: true` + read `carbon_audit` |
| Yield >100% vs oleate alone | Side-door carbon exchanges, not internal pools | Close non-feedstock carbon imports |
| Growth μ ~0.94 h⁻¹ on oleate | No experimental μ cap; model max ~9.35 h⁻¹ | `max_growth_rate` in scenario YAML |
| Wrong knockouts (POX1) | Gene names absent in model | `find_ids.py` → `gene_alias_resolution`, reaction IDs |
| Meaningless flux = 0 | Invented metabolite/reaction IDs | Preflight + copy from `recommended` |
| Agent cites flux as “expected titer” | No uncertainty language | `calibration.confidence_level` + `agent_guidance` |
| Infeasible glucose validation | `use_minimal_medium: true` without opening O₂/NH₄/etc. | Full cofeed list in `iyli647_glucose_validation.yaml` |

**Bottom line:** The agent is **more accurate relative to requirements** when it:

1. Runs biomass validation first (`objective="biomass"`, glucose scenario).
2. Uses **`wax_ester_oleate_n_limited.yaml`** (or a team-customized derivative with your fluxomics).
3. Reports **`yield_corrected_mol_per_mol_substrate`** and **`calibration.recommended_use`**.
4. Uses **`exploratory_no_literature.yaml`** only for early design sweeps.

---

## 6. Calibration levels (agent decision table)

Two axes are reported: **`medium_confidence_level`** and **`product_confidence_level`**.

| `confidence_level` | When | Agent should |
|---|---|---|
| `exploratory` | ≥3 missing medium/growth inputs | Rank designs relative to each other only |
| `partial` | 1–2 missing medium inputs | Rank + bottlenecks; in silico upper bound |
| `medium_calibrated` | Medium/growth pinned; no product literature | Trust μ/medium; **do not** cite product flux as validated titer |
| `literature_calibrated` | Biomass fully pinned, or product + `product_literature_refs` | Quantitative comparison after biomass check |
| `invalid` | FBA infeasible/error | Debug constraints; do not rank |

Add experimental wax data to `product_literature_refs` in scenario YAML to upgrade product confidence.

---

## 7. Known gaps & recommended next steps

### Gaps

1. **No direct literature benchmark for wax esters on iYLI647** — product scenarios are engineering templates, not paper reproductions.
2. **CER for Workman batch glucose** — scenario uses a generous CO₂ cap (25 mmol/gDW/h); tighten when exact batch qCO₂ is transcribed from Workman 2013 tables.
3. **Minimal medium completeness** — full ion/trace-metal exchange set may need expansion for other carbon sources.
4. **Units** — agent must not convert mmol/gDW/h to g/L without explicit bioreactor/QP model.

### Recommended team actions

| Priority | Action | Owner suggestion |
|---|---|---|
| High | Adopt `wax_ester_oleate_n_limited.yaml` as default agent scenario for WE scoring | Strain design / agent team |
| High | Require `find_ids.py --json` in agent orchestration (hard gate) | Agent platform (Mahbub) |
| Medium | Fill CER/O₂/NH₄ from your shake-flask or fluxomics into a `wax_ester_oleate_custom.yaml` | Wet-lab + modeling |
| Medium | Add wax-ester-specific validation when experimental flux data exist | Metabolic modeling |
| Low | Suppress COBRApy SBML parser warnings in CI logs | DevOps |

---

## 8. File reference

| File | Role |
|---|---|
| `fba_tool.py` | `score_pathway()` — main agent entry |
| `find_ids.py` | ID resolution + SBML preflight |
| `model_loader.py` | Robust SBML load/repair |
| `FBA_TOOL_CONTRACT.md` | Agent API contract |
| `scenarios/wax_ester_oleate_n_limited.yaml` | **Default production scenario** |
| `scenarios/wax_ester_oleate_open.yaml` | Optimistic upper bound (partial calibration) |
| `scenarios/iyli647_glucose_validation.yaml` | Biomass sanity check vs Mishra/Workman |
| `scenarios/exploratory_no_literature.yaml` | Novel-pathway template (exploratory only) |
| `run_all_tests.py` | One-command trial run + summary table |
| `testing_debugging/` | Pathway, stoichiometry, carbon verification |
| `iYLI647.xml` | Canonical SBML model |

---

## 9. Conclusion for stakeholders

The FBA tool **does improve agent output quality** when integrated as specified: preflight → scenario-driven constraints → calibration-aware reporting. It transforms raw FBA from a misleading “single number” into a **gated recommendation** (rank vs compare vs invalid).

For wax-ester strain design on oleate, **`medium_calibrated`** N-limited runs with
**minimal medium** (~1.1 mmol/gDW/h product, oleate sole carbon source) are the
most honest current predictions: medium is pinned, product flux is **not**
experimentally anchored. **Exploratory/open
runs** remain useful for comparing knockout sets but must not be presented as
expected experimental titers.

Biomass validation on glucose now **passes** (optimal, μ = 0.24 h⁻¹), establishing trust in model loading and medium setup before product objectives are run.

---

*Generated from trial run `python run_all_tests.py` on 4 June 2026. Re-run after scenario or model updates.*
