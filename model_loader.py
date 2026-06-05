"""
model_loader.py
===============
Robust SBML loading for agent-facing FBA tools.

Many published GSM files (including supplementary SBML) omit required reaction
or model `id` attributes. COBRApy refuses to load them. This module validates,
auto-repairs common issues, and returns a structured preflight report so agents
can fail loudly with actionable messages instead of opaque parse errors.
"""

from __future__ import annotations

import re
import warnings
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from cobra import Model
from cobra.io import load_model, read_sbml_model
from cobra.io.sbml import validate_sbml_model

SBML_NS = "http://www.sbml.org/sbml/level2"
SBML_EXTENSIONS = (".xml", ".sbml", ".xml.gz", ".sbml.gz")

# Known biomass reaction ids after repair of the bundled iYLI647 SBML.
IYLI647_BIOMASS_CARBON = "BIOMASS_yarrowia_carbon_limiting"


def is_sbml_path(model_ref: str) -> bool:
    return model_ref.lower().endswith(SBML_EXTENSIONS)


def sanitize_sbml_id(raw: str, prefix: str = "R") -> str:
    """Turn a human-readable name into a valid SBML SId."""
    text = (raw or "unnamed").strip().lower()
    text = re.sub(r"[^a-z0-9_]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    if not text:
        text = "unnamed"
    if text[0].isdigit():
        text = f"{prefix.lower()}_{text}"
    sid = f"{prefix}_{text}" if prefix else text
    if len(sid) > 240:
        sid = sid[:240]
    return sid


def _local(tag: str) -> str:
    return tag.split("}")[-1] if "}" in tag else tag


def repair_sbml_file(path: str | Path, *, write: bool = True) -> Dict[str, Any]:
    """
    Fix missing model / reaction ids in an SBML file.

    Returns a report dict with keys: path, repairs, ok.
    """
    path = Path(path)
    tree = ET.parse(path)
    root = tree.getroot()
    repairs: List[str] = []

    model_elem = None
    for elem in root.iter():
        if _local(elem.tag) == "model":
            model_elem = elem
            break

    if model_elem is None:
        return {"path": str(path), "repairs": repairs, "ok": False,
                "error": "No <model> element found"}

    if not model_elem.get("id"):
        model_elem.set("id", "model")
        repairs.append("Set missing model id -> 'model'")

    seen_ids: set[str] = set()
    for elem in root.iter():
        if _local(elem.tag) != "reaction":
            continue
        rid = elem.get("id")
        if rid:
            seen_ids.add(rid)
            continue
        name = elem.get("name") or "unnamed_reaction"
        base = sanitize_sbml_id(name, prefix="R")
        candidate = base
        n = 2
        while candidate in seen_ids:
            candidate = f"{base}_{n}"
            n += 1
        elem.set("id", candidate)
        seen_ids.add(candidate)
        repairs.append(f"Set missing reaction id for '{name}' -> '{candidate}'")

    if write and repairs:
        tree.write(path, encoding="UTF-8", xml_declaration=True)

    return {"path": str(path), "repairs": repairs, "ok": True}


def validate_sbml(path: str | Path) -> Dict[str, Any]:
    """Run COBRApy SBML validation and flatten blocking errors for agents."""
    _, errors = validate_sbml_model(str(path))
    blocking_levels = {
        "SBML_FATAL",
        "SBML_ERROR",
        "SBML_SCHEMA_ERROR",
        "COBRA_FATAL",
        "COBRA_ERROR",
    }
    flat_blocking: List[str] = []
    flat_warnings: List[str] = []
    for level, msgs in errors.items():
        for msg in msgs:
            text = f"{level}: {msg.strip()}"
            if level in blocking_levels:
                flat_blocking.append(text)
            else:
                flat_warnings.append(text)
    return {
        "valid": len(flat_blocking) == 0,
        "errors": flat_blocking,
        "warnings": flat_warnings[:20],
        "warning_count": len(flat_warnings),
        "by_level": {k: v for k, v in errors.items() if v},
    }


def ensure_sbml_loadable(path: str | Path) -> Dict[str, Any]:
    """
    Validate SBML; repair if needed; validate again.
    Returns a preflight report suitable for JSON agent logs.
    """
    path = Path(path)
    report: Dict[str, Any] = {
        "path": str(path),
        "exists": path.exists(),
        "repairs_applied": [],
        "validation_before": None,
        "validation_after": None,
        "loadable": False,
    }
    if not path.exists():
        report["error"] = f"File not found: {path}"
        return report

    report["validation_before"] = validate_sbml(path)
    if not report["validation_before"]["valid"]:
        repair = repair_sbml_file(path, write=True)
        report["repairs_applied"] = repair.get("repairs", [])
        if not repair.get("ok"):
            report["error"] = repair.get("error", "Repair failed")
            return report

    report["validation_after"] = validate_sbml(path)
    # Attempt load even if validator emits non-blocking warnings
    if report["validation_after"]["valid"]:
        report["loadable"] = True
    else:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                read_sbml_model(str(path))
            report["loadable"] = True
            report["validation_after"]["valid"] = True
            report["validation_after"]["errors"] = []
            report["loadable_via_fallback"] = True
        except Exception as exc:
            report["loadable"] = False
            report["error"] = f"SBML still invalid after repair: {exc}"
    return report


def load_model_robust(
    model_ref: str,
    *,
    auto_repair: bool = True,
    quiet: bool = True,
) -> Tuple[Model, Dict[str, Any]]:
    """
    Load a bundled cobra model or an SBML file.

    SBML files are validated first; common id omissions are auto-repaired.
    Returns (model, preflight_report).
    """
    report: Dict[str, Any] = {"model_ref": model_ref, "source": "bundled"}

    if is_sbml_path(model_ref):
        report["source"] = "sbml"
        path = Path(model_ref)
        if auto_repair:
            preflight = ensure_sbml_loadable(path)
            report["preflight"] = preflight
            if not preflight.get("loadable"):
                msg = preflight.get("error") or "SBML preflight failed"
                raise ValueError(
                    f"{msg}. Validation: {preflight.get('validation_after')}"
                )
        ctx = warnings.catch_warnings()
        ctx.__enter__()
        warnings.simplefilter("ignore")
        try:
            model = read_sbml_model(str(path))
        finally:
            ctx.__exit__(None, None, None)
    else:
        if quiet:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                model = load_model(model_ref)
        else:
            model = load_model(model_ref)

    report["model_id"] = model.id
    report["n_reactions"] = len(model.reactions)
    report["n_metabolites"] = len(model.metabolites)
    report["n_genes"] = len(model.genes)
    report["has_gene_associations"] = report["n_genes"] > 0
    return model, report
