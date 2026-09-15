"""Shared helpers, constants, and the number_input precision monkeypatch for
pbisim-app. Imported by app.py and the view modules.

This module has NO Streamlit entry side effects (no set_page_config / CSS /
session init) — those stay in app.py so import order is irrelevant here.
"""

from __future__ import annotations

import faulthandler
import copy
import dataclasses as _dc
import contextlib
import io
import json
import os
import re
import time
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st
from pbisim import (
    ModelBuilder,
    PBIModel,
    solve_ode,
    DoseSchedule,
    DoseEvent,
    time_to_clearance,
    time_to_log_reduction,
    stationary_phase_ic,
)
from pbisim.strains import StrainDefinition, StrainSet
from pbisim.strains.genotypes import BinaryResistanceGenotypes, BacterialStrain, PhageStrain, Antibiotic
from pbisim.pk.antibiotic import AntibioticDefinition, AntibioticSensitivity
from pbisim.trial.clinical import TreatmentArm
from pbisim_app.agent import SimulationAgent
from pbisim_app.executor import execute_code
from pbisim_app.viz_helper import (
    plot_axis_controls, apply_axis_plotly,
    build_series, series_selector, plot_series, plot_sweep_traces,
    BACTERIAL_BASES, CFU_BASIS_LABEL, CFU_PREFIXES, bacterial_total,
)
from pbisim_app.fit_helper import (
    OBSERVABLES,
    OBS_COMPARTMENTS,
    obs_prefixes,
    predicted_observable,
    normalize_fit_dataframe,
    arm_covariate_values,
    parse_dose_rows,
    DOSE_TARGETS,
    apply_row_filters,
    aggregate_observations,
    fit_residual,
    residual_vector_log10,
    model_information_criteria,
    compare_fit_models,
    build_fit_spec,
    config_param_snapshot,
)
from pbisim_app.trial_helper import (
    IIV_PARAMETERS,
    run_trial_simulation,
    plot_kaplan_meier_plotly,
    plot_metric_distributions_plotly,
    plot_pkpd_trajectories_plotly,
    build_regimen_doses,
)
from pbisim_app.sweep_helper import (
    get_sweep_parameters,
    apply_sweep_parameter,
    parse_comma_separated_series,
    pad_vectors,
    categorize_sweep_params,
)


# ── Precision fix for float number inputs ─────────────────────────────────────
# Streamlit's number_input infers a display format from the step (default "%.2f"
# for floats), so a value like 0.001 rounds to 0.00 on entry — losing precision
# for small rates (death, dormancy, resuscitation, decay, …). Inject a compact
# full-precision format ("%g") for any *float* input that doesn't set one. This
# auto-skips integer inputs (their value is an int) and any input that already
# specifies `format=` (e.g. scientific "%.1e").
_orig_number_input = st.number_input


def _number_input_precise(label, *args, **kwargs):
    if "format" not in kwargs and isinstance(kwargs.get("value"), float):
        kwargs["format"] = "%g"
    return _orig_number_input(label, *args, **kwargs)


st.number_input = _number_input_precise


def scripting_enabled() -> bool:
    """Whether the notebook-style Scripting page is available (opt-in). It runs arbitrary
    Python in the app's server process (research sandbox, not isolated), so it is OFF by
    default and only appears when ``PBISIM_ENABLE_SCRIPTING`` is set to a truthy value —
    turn it on only where you trust the users (a local run, or an authenticated deploy)."""
    import os as _os
    return _os.environ.get("PBISIM_ENABLE_SCRIPTING", "").strip().lower() in (
        "1", "true", "yes", "on")


# ── Dose defaults per target compartment ───────────────────────────────────────
# 1e8 only makes sense for phage (PFU); antibiotics are dosed in mg and nutrient
# in resource units, so give each target a plausible default and unit label.
DOSE_AMOUNT_DEFAULTS = {"phage": 1e8, "antibiotic": 10.0, "nutrient": 1.0}


DOSE_AMOUNT_LABELS = {
    "phage": "Amount (PFU)",
    "antibiotic": "Amount (mg)",
    "nutrient": "Amount (resource units)",
}


def render_regimen_config(prefix, items, target, default_amount, unit_label,
                          default_on=True, initial=None):
    """Render dose-regimen widgets for one agent and return a config dict.

    ``initial`` (an existing config dict) pre-fills the widgets so an arm can be edited
    in place. Returns ``{"on": False}`` when the agent is not included in this arm, else
    ``{"on": True, "index", "amount", "start", "repeat", "interval", "n"}``.
    """
    init = initial or {}
    on = st.checkbox(f"Include {target}", value=bool(init.get("on", default_on)), key=f"{prefix}_on")
    if not on:
        return {"on": False}
    idx = int(init.get("index", 0))
    if len(items) > 1:
        idx = st.selectbox(
            f"{target.title()} target", range(len(items)),
            index=min(idx, len(items) - 1),
            format_func=lambda i: items[i].get("name", f"{target} {i}"),
            key=f"{prefix}_idx",
        )
    c1, c2 = st.columns(2)
    with c1:
        amount = st.number_input(unit_label, min_value=0.0,
                                 value=float(init.get("amount", default_amount)),
                                 format="%.1e", key=f"{prefix}_amt")
    with c2:
        start = st.number_input("Start (h)", min_value=0.0,
                                value=float(init.get("start", 0.0)), step=1.0,
                                key=f"{prefix}_start")
    repeat = st.checkbox("Repeat regimen (qX h × N)", value=bool(init.get("repeat", False)),
                         key=f"{prefix}_rep")
    interval, n_doses = float(init.get("interval", 8.0)), int(init.get("n", 4))
    if repeat:
        c3, c4 = st.columns(2)
        with c3:
            interval = st.number_input("Interval (h)", min_value=0.5, value=interval,
                                       step=1.0, key=f"{prefix}_int")
        with c4:
            n_doses = st.number_input("Doses", min_value=1, value=n_doses, step=1,
                                      key=f"{prefix}_n")
    return {"on": True, "index": int(idx), "amount": float(amount),
            "start": float(start), "repeat": bool(repeat),
            "interval": float(interval), "n": int(n_doses)}


def render_iiv_config(prefix, initial=None):
    """Render the inter-individual-variability (IIV) form and return an iiv dict
    ``{path, dist_type, params, mode}``. ``initial`` pre-fills the widgets for editing."""
    init = initial or {}
    names = list(IIV_PARAMETERS.keys())
    cur_name = next((n for n, p in IIV_PARAMETERS.items() if p == init.get("path")), names[0])
    param_display = st.selectbox("Select Parameter", names,
                                 index=names.index(cur_name) if cur_name in names else 0,
                                 key=f"{prefix}_param")
    dists = ["LogNormal", "Normal", "Uniform"]
    cur_dist = init.get("dist_type", "LogNormal")
    dist_choice = st.selectbox("Distribution Type", dists,
                               index=dists.index(cur_dist) if cur_dist in dists else 0,
                               key=f"{prefix}_dist")
    ip = init.get("params", {})
    c1, c2 = st.columns(2)
    params = {}
    if dist_choice == "LogNormal":
        with c1:
            params["cv"] = st.number_input("CV (coefficient of variation)",
                                           value=float(ip.get("cv", 0.25)), min_value=0.01,
                                           key=f"{prefix}_cv")
        mode = "multiplicative"
    elif dist_choice == "Normal":
        with c1:
            params["mean"] = st.number_input("Mean", value=float(ip.get("mean", 0.0)),
                                             key=f"{prefix}_mean")
        with c2:
            params["sd"] = st.number_input("SD (standard deviation)",
                                           value=float(ip.get("sd", 0.1)), min_value=0.01,
                                           key=f"{prefix}_sd")
        mode = "additive"
    else:
        with c1:
            params["lo"] = st.number_input("Lower Bound", value=float(ip.get("lo", 0.5)),
                                           key=f"{prefix}_lo")
        with c2:
            params["hi"] = st.number_input("Upper Bound", value=float(ip.get("hi", 1.5)),
                                           key=f"{prefix}_hi")
        mode = "replace"
    return {"path": IIV_PARAMETERS[param_display], "dist_type": dist_choice,
            "params": params, "mode": mode}


def render_rate_graph_editor(names, state_key, key_prefix, *, rate_label="Rate (h⁻¹)",
                             default_rate=1e-7, add_label="+ Add edge"):
    """Edit a named directed rate graph stored as ``st.session_state[state_key]`` =
    ``[{"from": name, "to": name, "rate": r}, ...]``.

    Generic editor behind every generator-matrix input (bacterial mutation and
    phenotypic transitions, phage mutation and transitions): ``names`` are the node
    labels (strain names, BRG genotype labels or phage names). Works for any node count.
    """
    edges = st.session_state.get(state_key, [])
    names = list(names)
    for idx, tr in enumerate(list(edges)):
        c1, c2, c3, c4 = st.columns([3, 3, 3, 1])
        with c1:
            tr["from"] = st.selectbox("From", names,
                                      index=names.index(tr["from"]) if tr.get("from") in names else 0,
                                      key=f"{key_prefix}_src_{idx}")
        with c2:
            tr["to"] = st.selectbox("To", names,
                                    index=names.index(tr["to"]) if tr.get("to") in names else 0,
                                    key=f"{key_prefix}_dest_{idx}")
        with c3:
            tr["rate"] = st.number_input(rate_label, value=float(tr.get("rate", default_rate)),
                                         format="%.2e", key=f"{key_prefix}_rate_{idx}")
        with c4:
            st.markdown("<br>", unsafe_allow_html=True)
            if st.button(":material/delete:", key=f"{key_prefix}_del_{idx}"):
                edges.pop(idx)
                st.session_state[state_key] = edges
                st.rerun()
    if st.button(add_label, key=f"{key_prefix}_add"):
        edges.append({"from": names[0] if names else "", "to": names[0] if names else "",
                      "rate": default_rate})
        st.session_state[state_key] = edges
        st.rerun()


def render_mutation_graph_editor(strains, key_prefix):
    """Edit the named bacterial mutation graph (shared `int_transitions`).

    Each entry is {"from": strain_name, "to": strain_name, "rate": mu}. Works for
    any number of strains, so it lifts the 2^m restriction of the per-locus shortcut.
    """
    render_rate_graph_editor([s["name"] for s in strains], "int_transitions", key_prefix,
                             rate_label="Rate (mu)", add_label="+ Add transition")


def rate_matrix_from_graph(edges, names):
    """Build an (n,n) mass-conserving generator matrix from a named edge list.

    Convention (matches pbisim / ``gen_mat_from_graph``): M[dest, origin] = rate
    origin→dest; the diagonal M[o, o] = -(sum of outflows from o), so each column sums
    to zero and the matrix only redistributes (bacteria or phage) between the named
    states. Edges naming unknown nodes, self-loops and non-positive rates are ignored.
    Returns None if there are no valid edges (→ the caller leaves the engine default).
    """
    names = list(names)
    n = len(names)
    name_to_idx = {nm: i for i, nm in enumerate(names)}
    M = np.zeros((n, n))
    any_edge = False
    for tr in edges or []:
        o = name_to_idx.get(tr.get("from"))
        d = name_to_idx.get(tr.get("to"))
        r = float(tr.get("rate", 0.0) or 0.0)
        if o is not None and d is not None and o != d and r > 0:
            M[d, o] += r
            M[o, o] -= r
            any_edge = True
    return M if any_edge else None


def graph_dict_from_edges(edges, names):
    """Edge list → ``{origin: {dest: rate}}`` (the StrainSet ``set_*_graph`` format),
    keeping only valid origin→dest edges between the given node names."""
    names = set(names)
    graph = {}
    for tr in edges or []:
        src, dest = tr.get("from"), tr.get("to")
        r = float(tr.get("rate", 0.0) or 0.0)
        if src in names and dest in names and src != dest and r > 0:
            graph.setdefault(src, {})[dest] = r
    return graph


def mutation_matrix_from_transitions(transitions, strains):
    """Build the (n,n) mass-conserving mutation matrix from a named transition graph
    (see :func:`rate_matrix_from_graph`). Returns None if there are no valid edges."""
    return rate_matrix_from_graph(transitions, [s["name"] for s in strains])


def brg_genotype_labels(n_phages, n_antibiotics):
    """Genotype labels of a Binary-Genotypes model with ``n_phages`` phage loci and
    ``n_antibiotics`` antibiotic loci, in the engine's ordering (mirrors
    ``BinaryResistanceGenotypes.strain_labels``: ``'00'``… without antibiotics,
    ``'phi00_abx0'``… with, ``'abx0'``… when there are no phages)."""
    import itertools
    labels = []
    for comb in itertools.product([0, 1], repeat=n_phages + n_antibiotics):
        if n_antibiotics == 0:
            labels.append("".join(map(str, comb)))
        else:
            p_lbl = "".join(map(str, comb[:n_phages]))
            a_lbl = "".join(map(str, comb[n_phages:]))
            labels.append(f"phi{p_lbl}_abx{a_lbl}" if n_phages > 0 else f"abx{a_lbl}")
    return labels


# Session keys of the four generator-matrix edge lists (all ``int_*`` → captured by
# scenario snapshots). ``int_transitions`` is the (historically named) bacterial
# MUTATION graph; the other three are new (2026-09).
GENERATOR_GRAPH_KEYS = {
    "mutation_rates": "int_transitions",             # bacteria: division-coupled mutation
    "transition_rates": "int_bact_transitions",      # bacteria: phenotypic switching (h⁻¹)
    "mutation_rates_phage": "int_phage_mutations",   # phage: mutation at lysis (per burst)
    "transition_rates_phage": "int_phage_transitions",  # phage: phenotypic switching (h⁻¹)
}


def edges_from_rate_matrix(M, names):
    """Inverse of :func:`rate_matrix_from_graph`: the positive off-diagonal entries of a
    generator matrix as a named edge list (``M[dest, origin]`` → ``origin→dest``)."""
    names = list(names)
    M = np.asarray(M, dtype=float)
    edges = []
    if M.ndim != 2 or M.shape != (len(names), len(names)):
        return edges
    for o in range(len(names)):
        for d in range(len(names)):
            if o != d and M[d, o] > 0:
                edges.append({"from": names[o], "to": names[d], "rate": float(M[d, o])})
    return edges


def write_back_generator_matrices(cfg, bact_names, phage_names, *, include_mutation=True):
    """Write a (fitted) config's generator matrices back into the session edge lists so
    the builder shows the estimated rates. A graph is rewritten only when it was already
    configured (its edge list is non-empty) — an empty graph stays the engine default.
    ``include_mutation=False`` skips the bacterial mutation graph (BRG derives its
    mutation matrix from the per-locus μ, not from ``int_transitions``)."""
    for fld, key in GENERATOR_GRAPH_KEYS.items():
        if fld == "mutation_rates" and not include_mutation:
            continue
        if not st.session_state.get(key):
            continue
        M = getattr(cfg, fld, None)
        if M is None:
            continue
        names = phage_names if fld.endswith("_phage") else bact_names
        edges = edges_from_rate_matrix(M, names)
        if edges:
            st.session_state[key] = edges


def generator_matrices_from_state(bact_names, phage_names):
    """Resolve the three NON-mutation generator matrices (bacterial phenotypic
    transitions, phage mutation, phage transitions) from session state as
    ``{config_field: matrix-or-None}`` — None where the graph has no valid edge, so the
    caller can leave the engine's zero default untouched."""
    g = st.session_state.get
    return {
        "transition_rates": rate_matrix_from_graph(g("int_bact_transitions", []), bact_names),
        "mutation_rates_phage": rate_matrix_from_graph(g("int_phage_mutations", []), phage_names),
        "transition_rates_phage": rate_matrix_from_graph(g("int_phage_transitions", []), phage_names),
    }


# ── Cached calibration data processing ────────────────────────────────────────
# The Calibration page re-runs on every widget interaction; without caching it would
# re-parse the CSV and re-filter/normalise/aggregate the whole dataset each time.
# @st.cache_data memoises by content, so these recompute only when inputs change.
@st.cache_data(show_spinner=False)
def read_uploaded_csv(file):
    return pd.read_csv(file)


@st.cache_data(show_spinner=False)
def calibration_processed(raw, filters_key, time_col, value_col, observable,
                          group_cols, moi_col, stat, band):
    """Filter → normalise → aggregate, cached. Keys must be hashable (tuples)."""
    filtered = apply_row_filters(raw, {c: list(v) for c, v in filters_key})
    long, conds = normalize_fit_dataframe(filtered, time_col, value_col, observable,
                                          list(group_cols), moi_col)
    agg = aggregate_observations(long, stat=stat, band=band)
    return filtered, long, conds, agg


def arm_dose_events(arm):
    """Build the DoseEvent list for a treatment-arm config from its regimens."""
    doses = []
    for key, target in (("phage", "phage"), ("abx", "antibiotic")):
        c = arm.get(key)
        if c and c.get("on"):
            doses += build_regimen_doses(
                target, c["index"], c["amount"], c["start"],
                c["repeat"], c["interval"], c["n"],
            )
    return doses


def arm_regimen_summary(arm):
    """One-line human summary of an arm's dosing."""
    parts = []
    for key, label in (("phage", "P"), ("abx", "A")):
        c = arm.get(key)
        if c and c.get("on"):
            reg = f"q{c['interval']:g}h×{c['n']}" if c.get("repeat") else "single"
            parts.append(f"{label} {c['amount']:.1e} @ {c['start']:g}h ({reg})")
    return ", ".join(parts) if parts else "no doses (control-like)"


# ── State Initialization ──────────────────────────────────────────────────────
def _init_app_state():
    if "agent" not in st.session_state:
        st.session_state.agent = SimulationAgent()
    if "history" not in st.session_state:
        st.session_state.history = []
    if "current_page" not in st.session_state:
        st.session_state.current_page = "Interactive Simulator"
    if "simulation_result" not in st.session_state:
        st.session_state.simulation_result = None
    if "simulation_config" not in st.session_state:
        st.session_state.simulation_config = None
    if "int_builder_mode" not in st.session_state:
        st.session_state.int_builder_mode = "Direct (ModelBuilder)"
    if "api_models_list" not in st.session_state:
        st.session_state.api_models_list = []
    
    # Custom builders states
    if "int_transitions" not in st.session_state:
        st.session_state.int_transitions = []
    for _gk in ("int_bact_transitions", "int_phage_mutations", "int_phage_transitions"):
        if _gk not in st.session_state:
            st.session_state[_gk] = []
    if "int_brg_initial_B" not in st.session_state:
        st.session_state.int_brg_initial_B = {}
        
    # Clinical trial states
    if "trial_iiv_inputs" not in st.session_state:
        st.session_state.trial_iiv_inputs = []
    if "trial_arms" not in st.session_state:
        st.session_state.trial_arms = []
    if "trial_result" not in st.session_state:
        st.session_state.trial_result = None
    if "user_scenarios" not in st.session_state:
        # name -> {"annotation": str, "schema_version": int, "state": {...}}
        st.session_state.user_scenarios = {}
    if "user_models" not in st.session_state:
        # name -> {"description": str, "source": str, "schema_version": int,
        #          "state": {model data keys}} — organism/kinetics only (see dump_model)
        st.session_state.user_models = {}
    if "active_model" not in st.session_state:
        # Which Model the builder currently reflects (label in the Models list).
        st.session_state.active_model = WORKING_DRAFT_LABEL
    if "parts_library" not in st.session_state:
        # {category: {name: {"source","annotation","reference_host?","params"}}}
        st.session_state.parts_library = {"bacteria": {}, "phages": {}, "antibiotics": {}}
    if "fit_dataset" not in st.session_state:
        # {"long": DataFrame[time,arm,observable,value], "conditions": {arm: {"moi": float}}}
        st.session_state.fit_dataset = None


def _next_uid(prefix: str) -> str:
    """A session-stable unique id so per-row widget keys (dose / trial-arm editors)
    survive reorder/delete instead of being re-bound to a different row by list index."""
    st.session_state["_uid_counter"] = st.session_state.get("_uid_counter", 0) + 1
    return f"{prefix}{st.session_state['_uid_counter']}"


def _carry_prerun_debris(ic, kwargs):
    """Carry the OD/debris that accumulated during the stationary-phase pre-run into the
    treatment model — when the user opts to INHERIT it (default) rather than wash the dead
    cells out before treatment. No-op when debris isn't tracked (ic.Debris is None) or the
    'inherit' checkbox is unticked (washout)."""
    if st.session_state.get("int_prerun_inherit_debris", True) and getattr(ic, "Debris", None) is not None:
        kwargs["initial_Debris"] = ic.Debris


def _safe_od(result, total_bacteria):
    """OD trajectory that never raises. Uses the debris-aware ``get_od()`` when the debris
    ODE is in the config; otherwise falls back to biomass / conversion factor WITHOUT
    touching the Debris state (accessing it raises ``KeyError: Debris state not found``
    when debris isn't enabled). Guards against a session debris-flag vs. result mismatch."""
    factor = st.session_state.get("int_od_to_cfu_conversion_factor", 2e8) or 1.0
    try:
        return result.get_od()
    except Exception:
        return np.asarray(total_bacteria, dtype=float) / factor


def phage_pk_enabled(phages):
    """True when any configured phage has PK enabled (pk_mode ≠ 'None'), i.e. the result
    carries central-compartment ``Pc{j}`` (and, with 2-compartment, ``Pp{j}``) states."""
    return any((p or {}).get("pk_mode", "None") != "None" for p in (phages or []))


def central_phage_total(result, phages):
    """Total CENTRAL-compartment (blood) phage **concentration** = Σⱼ Pc{j}/Vc{j}
    (phages/mL), over phages with PK enabled — directly comparable to the infection-site
    titre ``sum_prefixes('P')``. Zeros if no phage PK / no Pc state in the result."""
    tot = None
    for j, p in enumerate(phages or []):
        if (p or {}).get("pk_mode", "None") == "None":
            continue
        try:
            pc = np.asarray(result.get(f"Pc{j}"), dtype=float)
        except Exception:
            continue
        vc = float(p.get("Vc", 5000.0)) or 1.0
        tot = pc / vc if tot is None else tot + pc / vc
    return tot if tot is not None else np.zeros_like(result.time, dtype=float)


def peripheral_phage_total(result, phages):
    """Total PERIPHERAL-tissue phage **amount** = Σⱼ Pp{j}, over phages using the
    2-compartment model (k12 > 0). The peripheral compartment has no volume in the engine,
    so this is an amount (not a concentration). Zeros if no 2-compartment phage / no Pp."""
    tot = None
    for j, p in enumerate(phages or []):
        if (p or {}).get("pk_mode", "None") == "None" or float((p or {}).get("k12", 0.0)) <= 0:
            continue
        try:
            pp = np.asarray(result.get(f"Pp{j}"), dtype=float)
        except Exception:
            continue
        tot = pp if tot is None else tot + pp
    return tot if tot is not None else np.zeros_like(result.time, dtype=float)


def phage_uses_two_compartment(phages):
    """True when any PK-enabled phage has a peripheral compartment (k12 > 0)."""
    return any((p or {}).get("pk_mode", "None") != "None" and float((p or {}).get("k12", 0.0)) > 0
               for p in (phages or []))


def _phage_pk_k12_k21(phages):
    """(k12, k21) arrays for ``PhagePKConfig`` when any PK-enabled phage has a peripheral
    compartment (k12 > 0), else ``(None, None)`` for a 1-compartment central model.
    ``PhagePKConfig`` requires k12 and k21 both set or both None; a PK-enabled phage that
    leaves k12=0 just stays 1-compartment (k12=k21=0)."""
    if not phage_uses_two_compartment(phages):
        return None, None
    k12 = np.array([float(p.get("k12", 0.0)) if p.get("pk_mode", "None") != "None" else 0.0
                    for p in phages])
    k21 = np.array([float(p.get("k21", 0.0)) if p.get("pk_mode", "None") != "None" else 0.0
                    for p in phages])
    return k12, k21


def _sweep_summary_tiles(df_summary):
    """Render a compact metric-tile row summarising a sweep's runs. Robust to
    missing/non-numeric columns (used by the Dose-Response and Parameter sweeps)."""
    import pandas as _pd
    n_runs = len(df_summary)
    if not n_runs:
        return
    nadir = _pd.to_numeric(df_summary.get("Nadir (cells/mL)"), errors="coerce")
    ct = _pd.to_numeric(df_summary.get("Clearance Time (h)"), errors="coerce")
    best_nadir = float(nadir.min()) if nadir is not None and nadir.notna().any() else None
    cleared = int(ct.notna().sum()) if ct is not None else 0
    fastest = float(ct.min()) if cleared else None
    tiles = [
        ("Runs", f"{n_runs}", "in this sweep"),
        ("Best nadir", f"{best_nadir:.2e}" if best_nadir is not None else "—", "lowest cells/mL"),
        ("Runs cleared", f"{cleared}/{n_runs}", f"below threshold"),
        ("Fastest clearance", f"{fastest:.1f} h" if fastest is not None else "—", "earliest eradication"),
    ]
    cols = st.columns(len(tiles))
    for col, (lbl, val, sub) in zip(cols, tiles):
        col.markdown(
            f"""<div class="metric-container">
                <div class="metric-label">{lbl}</div>
                <div class="metric-value">{val}</div>
                <div class="metric-sub">{sub}</div>
            </div>""",
            unsafe_allow_html=True,
        )
    st.markdown("<br>", unsafe_allow_html=True)


def counted_number_input(label, current_len, widget_key, *, min_value=0, max_value=10, help=None):
    """A "number of items" input that doesn't fight the list it resizes.

    The classic Streamlit bug: `st.number_input(value=len(items))` (no key) reverts
    the user's change because `value=` is recomputed from the very list the widget
    mutates. Here the widget is keyed (so it persists), and a marker reconciles it
    with the list: when the list length changes for an EXTERNAL reason (scenario /
    parts load, reset) the widget re-seeds to the new length; otherwise the widget
    drives the resize. Returns the chosen count (int)."""
    mark = f"_{widget_key}_mark"
    # Re-seed the widget from the list when the length changed externally (marker
    # mismatch) OR the widget key was dropped — Streamlit discards a widget's key
    # when it isn't rendered on a rerun (e.g. while another page is showing), and a
    # keyed number_input with no stored value would otherwise default to min_value
    # and silently trim the list.
    if st.session_state.get(mark) != current_len or widget_key not in st.session_state:
        st.session_state[widget_key] = int(min(max(current_len, min_value), max_value))
        st.session_state[mark] = current_len
    n = int(st.number_input(label, min_value=min_value, max_value=max_value, step=1,
                            key=widget_key, help=help))
    if n != current_len:
        st.session_state[mark] = n   # remember: this change was widget-driven
    return n


def load_preset_to_state(params: dict):
    """Deep copy preset parameters into st.session_state variables."""
    # 0. Clear old simulation results to prevent dimension mismatch crashes
    st.session_state.simulation_result = None
    st.session_state.simulation_config = None
    st.session_state.int_builder_mode = "Direct (ModelBuilder)"

    # 1. Clear old per-mode config + widget keys to prevent collisions and so a reset
    #    actually clears BRG / StrainSet (not just the Direct builder). Covers the Direct
    #    widget keys, the BRG data + widgets (int_brg_*, brg_*), StrainSet (int_transitions,
    #    ss_*, trans_*, direct_*), and the mode/signal selector widgets (widget_*).
    _clear_prefixes = (
        "strain_", "phage_", "abx_", "dose_", "ads_",
        "int_brg_", "brg_", "ss_", "str_", "direct_", "trans_", "widget_", "gen_",
    )
    keys_to_clear = [
        k for k in st.session_state.keys()
        if any(k.startswith(p) for p in _clear_prefixes)
        or k in ("int_transitions", "int_bact_transitions", "int_phage_mutations", "int_phage_transitions")
    ]
    for k in keys_to_clear:
        st.session_state.pop(k, None)

    # 2. Global settings
    st.session_state["int_t_end"] = params.get("t_end", 48.0)
    st.session_state["int_dt"] = params.get("dt", 0.25)
    st.session_state["int_extinction_threshold"] = params.get("extinction_threshold", 1.0)
    st.session_state["int_extinction_check_interval"] = params.get("extinction_check_interval", 0.0)
    st.session_state["int_n_latent"] = params.get("n_latent", 5)
    st.session_state["int_solver_method"] = params.get("solver_method", "BDF")
    st.session_state["int_track_nutrients"] = params.get("track_nutrients", True)
    st.session_state["int_initial_S"] = params.get("initial_S", 1.0)
    st.session_state["int_monod_constant"] = params.get("monod_constant", 0.3)
    st.session_state["int_recycle_fraction"] = params.get("recycle_fraction", 0.0)
    st.session_state["int_s_in"] = params.get("s_in", 0.0)
    st.session_state["int_s_out"] = params.get("s_out", 0.0)
    st.session_state["int_carrying_capacity"] = params.get("carrying_capacity", 1e9)
    # Growth signal: prefer an explicit growth_function; else derive from the legacy
    # track_nutrients flag (True → Monod nutrient growth, False → logistic density).
    _valid_gf = {fn for fn, _ in GROWTH_SIGNALS.values()}
    _gf = params.get("growth_function")
    if _gf not in _valid_gf:
        _gf = "monod_growth" if params.get("track_nutrients", True) else "logistic_growth"
    st.session_state["int_growth_function"] = _gf
    st.session_state["int_track_nutrients"] = _gf in _GROWTH_NUTRIENT
    st.session_state["int_density_growth_constant"] = params.get("density_growth_constant", 1e9) or 1e9
    # Diauxic (sequential_monod) per-phase parameters, when present.
    _seq_rf = params.get("growth_phase_rate_factors")
    _seq_mc = params.get("growth_phase_monod")
    _seq_th = params.get("growth_phase_thresholds")
    if _gf == _GROWTH_SEQUENTIAL and _seq_rf is not None and _seq_mc is not None:
        _rf = list(np.atleast_1d(np.asarray(_seq_rf, dtype=float)))
        _mc = list(np.atleast_1d(np.asarray(_seq_mc, dtype=float)))
        _th = list(np.atleast_1d(np.asarray(_seq_th, dtype=float))) if _seq_th is not None else []
        st.session_state["int_growth_n_phases"] = len(_rf)
        st.session_state["gp_monod_0"] = _mc[0] if _mc else 0.3
        for _i in range(1, len(_rf)):
            st.session_state[f"gp_rate_{_i}"] = _rf[_i]
            st.session_state[f"gp_monod_{_i}"] = _mc[_i] if _i < len(_mc) else 0.3
        for _i in range(len(_th)):
            st.session_state[f"gp_thresh_{_i}"] = _th[_i]
    st.session_state["int_death_function"] = params.get("death_function_name", "constant_death")
    # Lysis signal: frac_lysis (nutrient-coupled) vs constant_lysis (default).
    st.session_state["int_lysis_function"] = params.get("lysis_function_name", "constant_lysis")
    st.session_state["int_monod_constant_lysis"] = params.get("monod_constant_lysis", 0.3) or 0.3
    st.session_state["int_lysis_floor"] = float(params.get("lysis_floor", 0.0) or 0.0)
    st.session_state["int_density_total_cells"] = params.get("density_signal_uses_total_cells", False)
    st.session_state["int_superinfection"] = params.get("allow_superinfection", False)
    st.session_state["int_t_prerun"] = params.get("t_prerun", 0.0)
    st.session_state["int_prerun_inherit_debris"] = params.get("prerun_inherit_debris", True)

    # 3. Immunity settings
    st.session_state["int_immunity_enabled"] = params.get("immunity_enabled", False)
    st.session_state["int_innate_kill_rate"] = params.get("innate_kill_rate", 1e7)
    st.session_state["int_innate_kill50"] = params.get("innate_kill50", 1e5)
    st.session_state["int_innate_max"] = params.get("innate_max", 1e7)
    # backward compat: old presets stored "adaptive_decay_rate"
    st.session_state["int_innate_decay_rate"] = params.get(
        "innate_decay_rate", params.get("adaptive_decay_rate", 0.1)
    )
    st.session_state["int_imm_kill_rate_D"] = params.get("innate_kill_rate_D", 0.0)
    # translate legacy "adaptive" (invented by scaffold, not a pbisim module) → "innate"
    _raw_module = params.get("immune_module", "innate")
    st.session_state["int_immune_module"] = "innate" if _raw_module == "adaptive" else _raw_module
    # new fields (missing from older presets → sensible defaults)
    st.session_state["int_imm_stim_rate"] = params.get("imm_stim_rate", 0.1)
    st.session_state["int_imm_stim50"] = params.get("imm_stim50", 1e6)
    st.session_state["int_imm_initial"] = params.get("imm_initial", 0.0)

    # 4. OD & Debris settings
    st.session_state["int_debris_enabled"] = params.get("debris_enabled", False)
    st.session_state["int_debris_u"] = params.get("debris_u", 0.4)
    st.session_state["int_debris_v"] = params.get("debris_v", 0.2)
    st.session_state["int_debris_kdis"] = params.get("debris_kdis", 0.01)
    st.session_state["int_dormant_od_fraction"] = params.get("dormant_od_fraction", 1.0)
    st.session_state["int_od_to_cfu_conversion_factor"] = params.get("od_to_cfu_conversion_factor", 2e8)
    # Nutrient-dependent OD absorptivity (on when the config carries a function; else off).
    st.session_state["int_od_nutrient_enabled"] = params.get("od_absorptivity_function") is not None
    st.session_state["int_od_abs_exp"] = params.get("od_abs_exp", 1.0) or 1.0
    st.session_state["int_od_abs_stat"] = params.get("od_abs_stat", 0.4) or 0.4
    st.session_state["int_od_abs_S50"] = params.get("od_abs_S50", 0.3) or 0.3
    st.session_state["int_od_abs_hill"] = params.get("od_abs_hill", 2.0) or 2.0

    # 5. Strains list
    strains_list = []
    for i, s in enumerate(params.get("strains", [])):
        strains_list.append(
            {
                "name": s.get("name", f"Strain {i}"),
                "initial_B": s.get("initial_B", 1e7),
                "growth_rate": s.get("growth_rate", 1.2),
                "bacteria_to_resource_ratio": s.get("bacteria_to_resource_ratio", 1e9),
                "death_rate_B": s.get("death_rate_B", 0.0),
                "death_rate_D": s.get("death_rate_D", 0.0),
                "dormancy_enabled": s.get("dormancy_enabled", False),
                "dormancy_depth": s.get("dormancy_depth", 1),
                "dormancy_rate": s.get("dormancy_rate", 0.001),
                "resuscitation_rate": s.get("resuscitation_rate", 0.1),
                "dormancy_diffusion_rate": s.get("dormancy_diffusion_rate", 0.05),
                "dormancy_signal": s.get("dormancy_signal", "nutrient"),
                "resuscitation_signal": s.get("resuscitation_signal", "nutrient"),
                "diffusion_signal": s.get("diffusion_signal", "constant"),
                "initial_D": s.get("initial_D", 0.0),
            }
        )
    st.session_state["int_strains"] = strains_list

    # Dormancy signal FUNCTIONS + Ks/threshold are now MODEL-WIDE. Seed them from a
    # config-level value if present, else from the first strain's stored (legacy per-strain)
    # signal — so old presets/scenarios keep their dormancy-signal choice.
    _pf = next((s for s in strains_list if s.get("dormancy_enabled")),
               strains_list[0] if strains_list else {})
    st.session_state["int_dormancy_signal"] = canonical_signal(
        params.get("dormancy_signal", _pf.get("dormancy_signal", "nutrient")))
    st.session_state["int_resuscitation_signal"] = canonical_signal(
        params.get("resuscitation_signal", _pf.get("resuscitation_signal", "nutrient")))
    st.session_state["int_diffusion_signal"] = canonical_signal(
        params.get("diffusion_signal", _pf.get("diffusion_signal", "constant")))
    st.session_state["int_dormancy_monod_constant"] = float(
        params.get("dormancy_monod_constant", _pf.get("dormancy_monod_constant", 0.0)) or 0.0)
    st.session_state["int_dormancy_carrying_capacity"] = float(
        params.get("dormancy_carrying_capacity", _pf.get("dormancy_carrying_capacity", 1e8)) or 1e8)

    # 6. Phages list
    phages_list = []
    for i, p in enumerate(params.get("phages", [])):
        phages_list.append(
            {
                "name": p.get("name", f"Phage {i}"),
                "initial_P": p.get("initial_P", 1e6),
                "burst_sizes": p.get("burst_sizes", 50.0),
                "latent_periods": p.get("latent_periods", 0.5),
                "phage_decay_rates": p.get("phage_decay_rates", 0.1),
                "infected_nutrient_consumption": p.get("infected_nutrient_consumption", 0.0),
                "pk_mode": p.get("pk_mode", "None"),
                "Vc": p.get("Vc", 5000.0),
                "k_elim": p.get("k_elim", 0.2),
                "k_in": p.get("k_in", 0.1),
                "k_out": p.get("k_out", 0.05),
                "Vi": p.get("Vi", 10.0),
                # MM clearances
                "Km_elim": p.get("Km_elim", 0.0),
                "phage_decay_Km": p.get("phage_decay_Km", 0.0),
                # Pseudolysogeny
                "hibernation_rate_s": p.get("hibernation_rate_s", 0.0),
                "hibernation_rate_r": p.get("hibernation_rate_r", 0.0),
                "lytic_resumption_rate_s": p.get("lytic_resumption_rate_s", 0.0),
                "lytic_resumption_rate_r": p.get("lytic_resumption_rate_r", 0.0),
                # BRG specific
                "adsorption_s": p.get("adsorption_s", p.get("adsorption_rates")[0] if isinstance(p.get("adsorption_rates"), (list, np.ndarray)) else (p.get("adsorption_rates") if p.get("adsorption_rates") is not None else 5e-8)),
                "adsorption_r": p.get("adsorption_r", 0.0),
                "fitness_cost": p.get("fitness_cost", 0.05),
                "mu": p.get("mu", 1e-7),
            }
        )
    st.session_state["int_phages"] = phages_list

    # 7. Adsorption matrices (restore from list in presets if any)
    # Default: WT strain (s_idx=0) gets 1e-8; resistant strains get 0.
    for s_idx in range(len(strains_list)):
        for p_idx in range(len(phages_list)):
            p_orig = params.get("phages", [])[p_idx]
            ads = p_orig.get("adsorption_rates", None)
            ads_dorm = p_orig.get("adsorption_rates_dormant", 0.0)

            # Resolve lists or scalars; None → per-strain default
            if isinstance(ads, list):
                val_ads = ads[s_idx]
            elif ads is not None:
                val_ads = ads
            else:
                val_ads = 1e-8 if s_idx == 0 else 0.0
            val_ads_dorm = ads_dorm[s_idx] if isinstance(ads_dorm, list) else ads_dorm

            st.session_state[f"ads_{s_idx}_{p_idx}"] = float(val_ads)
            st.session_state[f"ads_dorm_{s_idx}_{p_idx}"] = float(val_ads_dorm)

    # 8. Antibiotics list
    abx_list = []
    for i, a in enumerate(params.get("antibiotics", [])):
        abx_list.append(
            {
                "name": a.get("name", f"Drug {i}"),
                "Vc": a.get("Vc", 250.0),
                "k_elim": a.get("k_elim", 0.3),
                "k12": a.get("k12", 0.0),
                "k21": a.get("k21", 0.0),
                "emax": a.get("emax", 3.0),
                "ec50": a.get("ec50", 0.2),
                "hill": a.get("hill", 1.5),
                "f_lyse": a.get("f_lyse", 0.0),
                "inoculum_effect_constant": a.get("inoculum_effect_constant", 0.0),
                "Km_elim": a.get("Km_elim", 0.0),
                # BRG specific
                "emax_r": a.get("emax_r", a.get("emax", 3.0) * 0.1),
                "ec50_r": a.get("ec50_r", a.get("ec50", 0.2) * 10.0),
                "fitness_cost": a.get("fitness_cost", 0.05),
                "mu": a.get("mu", 1e-7),
            }
        )
    st.session_state["int_antibiotics"] = abx_list

    # 9. Doses
    doses_list = []
    for d in params.get("doses", []):
        doses_list.append(
            {
                "time": d.get("time", 0.0),
                "amount": d.get("amount", 1e9),
                "target_type": d.get("target_type", "phage"),
                "target_idx": d.get("target_idx", 0),
                "route": d.get("route", "bolus"),
                "duration": d.get("duration", 0.0),
            }
        )
    st.session_state["int_doses"] = doses_list


def configure_summary(config: dict) -> str:
    """One-line human summary of an AI-supplied simulator configuration (no Streamlit)."""
    strains = config.get("strains") or []
    bits = [f"{len(strains)} strain(s)"]
    for key, label in (("phages", "phage"), ("antibiotics", "antibiotic"), ("doses", "dose event")):
        n = len(config.get(key) or [])
        if n:
            bits.append(f"{n} {label}(s)")
    extras = []
    if any(s.get("dormancy_enabled") for s in strains):
        extras.append("dormancy")
    if config.get("immunity_enabled"):
        extras.append("immunity")
    if config.get("debris_enabled"):
        extras.append("OD/debris")
    tail = f"; t_end={config.get('t_end', 48.0)} h"
    if extras:
        tail += ", " + ", ".join(extras)
    return ", ".join(bits) + tail


def apply_ai_configuration(config: dict) -> str:
    """Handler for the assistant's ``configure_simulator`` tool: populate the Interactive
    Simulator from the model's structured config via load_preset_to_state, then set the
    chosen builder mode (Direct / BRG / StrainSet). Returns a summary for the model, or an
    ``ERROR ...`` string (which the model reads as a signal to fall back to run_pbisim_code)."""
    try:
        strains = config.get("strains") or []
        if not strains:
            return "ERROR: configuration needs at least one strain."
        mode = (config.get("builder_mode") or "direct").lower()

        # Shared entities + globals (this also resets to Direct and clears mode-specific keys).
        load_preset_to_state(config)

        if mode == "brg":
            base = strains[0]
            st.session_state["int_builder_mode"] = "Binary Genotypes (BRG)"
            st.session_state["int_brg_base_growth"] = base.get("growth_rate", 1.2)
            st.session_state["int_brg_base_ratio"] = base.get("bacteria_to_resource_ratio", 1e9)
            _dorm = bool(base.get("dormancy_enabled", False))
            st.session_state["int_brg_dormancy_enabled"] = _dorm
            st.session_state["int_brg_dorm_rate"] = base.get("dormancy_rate", 0.001)
            st.session_state["int_brg_resus_rate"] = base.get("resuscitation_rate", 0.1)
            st.session_state["int_brg_diff_rate"] = base.get("dormancy_diffusion_rate", 0.05)
            st.session_state["int_brg_death_rate_B"] = base.get("death_rate_B", 0.0)
            st.session_state["int_brg_death_rate_D"] = base.get("death_rate_D", 0.0)
            st.session_state["int_brg_use_eq_ic"] = bool(config.get("equilibrium_ic", False))
            st.session_state["int_brg_eq_total_B"] = config.get("total_bacteria", 1e7)
            n_phg, n_abx = len(config.get("phages") or []), len(config.get("antibiotics") or [])
            summary = (f"Binary Genotypes (BRG): base strain '{base.get('name', 'WT')}' with "
                       f"{n_phg} phage locus/loci + {n_abx} antibiotic locus/loci "
                       f"({2 ** (n_phg + n_abx)} genotypes)")
        elif mode == "strainset":
            st.session_state["int_builder_mode"] = "Custom Strains & Graph (StrainSet)"
            edges = [
                {"from": g.get("from", ""), "to": g.get("to", ""), "rate": g.get("rate", 0.0)}
                for g in (config.get("mutation_graph") or [])
                if g.get("from") and g.get("to")
            ]
            st.session_state["int_transitions"] = edges
            summary = (f"Custom Strains (StrainSet): {len(strains)} named strain(s), "
                       f"{len(edges)} mutation edge(s)")
        else:  # direct
            st.session_state["int_builder_mode"] = "Direct (ModelBuilder)"
            summary = configure_summary(config)

        # Generator graphs (any mode): bacterial phenotypic transitions + phage mutation /
        # transitions. Nodes are strain names (BRG: genotype labels) / phage names.
        _gen_bits = []
        for _ck, _sk, _lab in (("transition_graph", "int_bact_transitions", "phenotypic transition"),
                               ("phage_mutation_graph", "int_phage_mutations", "phage mutation"),
                               ("phage_transition_graph", "int_phage_transitions", "phage transition")):
            _edges = [
                {"from": g.get("from", ""), "to": g.get("to", ""), "rate": float(g.get("rate", 0.0) or 0.0)}
                for g in (config.get(_ck) or [])
                if g.get("from") and g.get("to")
            ]
            st.session_state[_sk] = _edges
            if _edges:
                _gen_bits.append(f"{len(_edges)} {_lab} edge(s)")
        if _gen_bits:
            summary += ", " + ", ".join(_gen_bits)

        return f"Configured the Interactive Simulator — {summary}; t_end={config.get('t_end', 48.0)} h."
    except Exception as e:
        return f"ERROR applying configuration: {e}"


def summarize_current_results(_inp=None) -> str:
    """Handler for the assistant's ``get_simulation_summary`` tool: a compact metrics
    summary of the current Interactive Simulator result, for the model to interpret."""
    res = st.session_state.get("simulation_result")
    cfg = st.session_state.get("simulation_config")
    if res is None or cfg is None:
        return ("No simulation results are available yet — ask the user to run a simulation "
                "in the Interactive Simulator (or set one up and run it) first.")
    try:
        t = np.asarray(res.time, dtype=float)
        # Culturable CFU = B + D (infected I / hibernating H cells don't form colonies).
        total = np.asarray(res.sum_prefixes("B", "D"), dtype=float)
        lines = [f"Current simulation ({st.session_state.get('int_builder_mode', '?')}), "
                 f"t = {t[0]:.1f}–{t[-1]:.1f} h."]
        lines.append(
            f"Total bacteria (CFU/mL): start {total[0]:.2e}, end {total[-1]:.2e}, "
            f"nadir {total.min():.2e}, peak {total.max():.2e} "
            f"({'net decline' if total[-1] < total[0] else 'net growth/regrowth'})."
        )
        try:
            from pbisim import time_to_clearance, time_to_log_reduction
            thr = st.session_state.get("int_extinction_threshold", 1.0) or 1.0
            tc = time_to_clearance(res, threshold=thr)
            t2 = time_to_log_reduction(res, n_logs=2.0)
            lines.append(
                f"Time to clearance (<{thr:g}): {('%.1f h' % tc) if tc is not None else 'NOT cleared'}. "
                f"Time to 2-log reduction: {('%.1f h' % t2) if t2 is not None else 'not reached'}."
            )
        except Exception:
            pass
        # per-strain final active bacteria + non-WT fraction
        n = int(getattr(cfg, "n_bacteria", 0) or 0)
        names = [s.get("name", f"Strain {i}") for i, s in enumerate(st.session_state.get("int_strains", []))]
        finals = []
        for i in range(n):
            try:
                bi = float(np.asarray(res.get(f"B{i}"))[-1])
            except Exception:
                bi = float("nan")
            finals.append((names[i] if i < len(names) else f"Strain {i}", bi))
        if finals:
            lines.append("Final active bacteria per strain: "
                         + ", ".join(f"{nm} {bi:.2e}" for nm, bi in finals) + ".")
            tot_active = sum(b for _, b in finals if b == b)
            if n >= 2 and tot_active > 0:
                res_frac = sum(b for _, b in finals[1:] if b == b) / tot_active
                lines.append(f"Non-first-strain fraction of active bacteria at end: {100 * res_frac:.1f}%.")
        try:
            p = np.asarray(res.sum_prefixes("P"), dtype=float)
            lines.append(f"Free phage (PFU/mL): end {p[-1]:.2e}, peak {p.max():.2e}.")
        except Exception:
            pass
        if st.session_state.get("int_immunity_enabled", False):
            try:
                imm = np.asarray(res.get("Imm"), dtype=float)
                lines.append(f"Immune level: end {imm[-1]:.2e}, peak {imm.max():.2e}.")
            except Exception:
                pass
        if st.session_state.get("int_debris_enabled", False):
            try:
                od = np.asarray(res.get_od(), dtype=float)
                lines.append(f"Optical density: end {od[-1]:.3g}, peak {od.max():.3g}.")
            except Exception:
                pass
        return "\n".join(lines)
    except Exception as e:
        return f"Could not summarize the results: {e}"


# ── Scenario snapshots (Tier 1: full-config save / load / export / import) ─────
# A "scenario" is everything needed to reproduce a simulation: the builder mode,
# strains / phages / antibiotics, pairwise adsorption, dosing, nutrient, immune,
# debris, solver, prerun, and the trial design. We snapshot the *input* session
# keys directly (rather than maintaining an inverse of load_preset_to_state), so
# new parameters are captured automatically and every builder mode is covered.
SCENARIO_SCHEMA_VERSION = 1


# Non-int_ session keys that are still scenario *data* (not widgets).
_SCENARIO_EXTRA_DATA_KEYS = ("direct_phg_res_rates", "trial_arms", "trial_iiv_inputs")


# Pairwise adsorption data keys written by the phage-config widgets.
_ADS_DATA_RE = re.compile(r"^ads_(dorm_)?\d+_\d+$")


# Widget-key prefixes to clear on load so widgets re-read the restored values
# (Streamlit keeps a widget's value under its key and would otherwise override
# the freshly-loaded data).
_SCENARIO_WIDGET_PREFIXES = (
    "int_", "str_", "phg_", "ss_", "brg_", "abx_", "ads_", "dose_",
    "rep_dose", "single_dose", "new_arm", "trial_phg", "trial_abx",
    "widget_builder_mode",
    "trans_", "dir_trans_", "gen_",   # rate-graph editor widgets (mutation / transitions)
)


def _is_scenario_data_key(k: str) -> bool:
    return (
        k.startswith("int_")
        or k in _SCENARIO_EXTRA_DATA_KEYS
        or bool(_ADS_DATA_RE.match(k))
    )


def _json_safe(obj):
    """Coerce numpy scalars/arrays to plain Python for JSON serialisation."""
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, (str, int, float, bool)) or obj is None:
        return obj
    return str(obj)   # last resort so export never crashes on an odd object


def dump_state_to_scenario() -> dict:
    """Snapshot the current input configuration as a JSON-safe scenario state."""
    return {
        k: _json_safe(copy.deepcopy(st.session_state[k]))
        for k in list(st.session_state.keys())
        if _is_scenario_data_key(k)
    }


def load_scenario_to_state(state: dict) -> None:
    """Restore a scenario snapshot, clearing stale widget state first."""
    # Clear both the data keys we are about to overwrite and the widget keys that
    # would otherwise shadow them.
    for k in list(st.session_state.keys()):
        if k in _SCENARIO_EXTRA_DATA_KEYS or _ADS_DATA_RE.match(k) or any(
            k.startswith(p) for p in _SCENARIO_WIDGET_PREFIXES
        ):
            st.session_state.pop(k, None)
    for k, v in state.items():
        st.session_state[k] = copy.deepcopy(v)
    # Invalidate stale outputs
    st.session_state.simulation_result = None
    st.session_state.simulation_config = None
    st.session_state.trial_result = None


def export_scenarios_json(scenarios: dict) -> str:
    """Serialise the whole scenario library to a portable JSON string."""
    return json.dumps(
        {"schema_version": SCENARIO_SCHEMA_VERSION, "scenarios": _json_safe(scenarios)},
        indent=2,
    )


def import_scenarios_json(text: str) -> dict:
    """Parse an exported library; returns the {name: scenario} mapping.

    Tolerates both the wrapped ({"scenarios": {...}}) and bare ({name: ...})
    forms. Raises ValueError on malformed input.
    """
    data = json.loads(text)
    if isinstance(data, dict) and "scenarios" in data:
        scenarios = data["scenarios"]
    else:
        scenarios = data
    if not isinstance(scenarios, dict):
        raise ValueError("Not a scenario library (expected a JSON object).")
    for name, sc in scenarios.items():
        if not isinstance(sc, dict) or "state" not in sc:
            raise ValueError(f"Scenario '{name}' is missing its 'state'.")
    return scenarios


# ── Model registry ─────────────────────────────────────────────────────────────
# A **Model** is the organism + kinetics + environment + solver configuration — the
# subset of the input state that `build_nominal_config_from_gui` consumes. It is
# deliberately NARROWER than a Scenario: dosing (`int_doses`), trial arms/IIV, and
# every analysis-page setting (sweep ranges, calibration selections) are properties
# of an *experiment*, not the model, so they are EXCLUDED. This decoupling is the
# whole point — a saved Model is reused verbatim across the simulator, sweeps,
# trials, and fitting without the live builder widgets contaminating those tasks.
MODEL_SCHEMA_VERSION = 1

# Analysis-only keys that live inside the scenario partition but must NOT be part of
# a Model (they belong to the page running the experiment).
_MODEL_EXCLUDE_KEYS = {"int_doses", "trial_arms", "trial_iiv_inputs"}


def _is_model_data_key(k: str) -> bool:
    """A model = the scenario data partition minus the analysis-only keys."""
    return _is_scenario_data_key(k) and k not in _MODEL_EXCLUDE_KEYS


def dump_model() -> dict:
    """Snapshot the current live builder state as a JSON-safe Model (organism +
    kinetics + environment + solver only — no dosing / trial / analysis settings)."""
    return {
        k: _json_safe(copy.deepcopy(st.session_state[k]))
        for k in list(st.session_state.keys())
        if _is_model_data_key(k)
    }


# Model widget-key prefixes to clear when loading a Model into the builder. Excludes
# the analysis widget prefixes (dose_/rep_dose/single_dose/new_arm/trial_*) so that
# swapping the organism model leaves the page's dosing / trial design untouched.
_MODEL_WIDGET_PREFIXES = (
    "int_", "str_", "phg_", "ss_", "brg_", "abx_", "ads_", "widget_builder_mode",
)


def apply_model_to_state(state: dict) -> None:
    """Load a Model into the live builder, replacing the current organism/kinetics
    config but leaving dosing / trial / analysis state intact."""
    for k in list(st.session_state.keys()):
        if k in _MODEL_EXCLUDE_KEYS:
            continue
        if (_is_model_data_key(k) or _ADS_DATA_RE.match(k)
                or any(k.startswith(p) for p in _MODEL_WIDGET_PREFIXES)):
            # never drop the analysis keys that share the int_ namespace
            if k in _MODEL_EXCLUDE_KEYS:
                continue
            st.session_state.pop(k, None)
    for k, v in state.items():
        if k in _MODEL_EXCLUDE_KEYS:
            continue
        st.session_state[k] = copy.deepcopy(v)
    st.session_state.simulation_result = None
    st.session_state.simulation_config = None


@contextlib.contextmanager
def model_config_context(snapshot: dict | None):
    """Temporarily swap the live model DATA keys to a frozen Model snapshot, so a
    downstream build reads that Model regardless of the current builder widgets, then
    restore. Only the model data keys are swapped — dosing / analysis state (owned by
    the page) is left in place. ``snapshot=None`` is a no-op (use the live draft)."""
    if snapshot is None:
        yield
        return
    model_keys = [k for k in list(st.session_state.keys()) if _is_model_data_key(k)]
    saved = {k: copy.deepcopy(st.session_state[k]) for k in model_keys}
    try:
        for k in model_keys:
            st.session_state.pop(k, None)
        for k, v in snapshot.items():
            if k in _MODEL_EXCLUDE_KEYS:
                continue
            st.session_state[k] = copy.deepcopy(v)
        yield
    finally:
        for k in [kk for kk in list(st.session_state.keys()) if _is_model_data_key(kk)]:
            st.session_state.pop(k, None)
        for k, v in saved.items():
            st.session_state[k] = v


def build_config_from_model(snapshot: dict | None = None):
    """Build a pbisim config from a frozen Model snapshot (or the live builder draft
    when ``snapshot is None``). Returns whatever `build_nominal_config_from_gui`
    returns: (config, initial_B, initial_P, initial_S, model_kwargs)."""
    with model_config_context(snapshot):
        return build_nominal_config_from_gui()


def _snap_fmt(v):
    """Compact human formatting for a scalar / 1-D / 2-D value in the snapshot view."""
    if v is None:
        return "—"
    if callable(v) or hasattr(v, "__name__"):
        return getattr(v, "__name__", str(v))
    if isinstance(v, (bool, str)):
        return str(v)
    a = np.asarray(v, dtype=object) if not np.isscalar(v) else None
    if a is None or np.ndim(v) == 0:
        try:
            f = float(v)
        except (TypeError, ValueError):
            return str(v)
        if f == 0:
            return "0"
        return f"{f:g}" if 1e-3 <= abs(f) < 1e5 else f"{f:.3e}"
    arr = np.asarray(v)
    if arr.ndim == 1:
        return "[" + ", ".join(_snap_fmt(x) for x in arr) + "]"
    return "[" + "; ".join(_snap_fmt(row) for row in arr) + "]"


def render_model_snapshot(container=None, *, snapshot: dict | None = None):
    """Render a readable, sectioned summary of the fully-resolved model configuration —
    organism / kinetics / nutrient environment / immunity / OD-debris / initial
    conditions / solver — so the whole config is visible in one place instead of
    hunting through the builder tabs. ``snapshot`` = a `dump_model()` dict for a frozen
    Model; ``None`` = the live builder draft. Reads the *built config* so it's
    mode-agnostic (Direct/BRG/StrainSet all resolve to the same fields) and can't drift
    from what actually runs."""
    c = container if container is not None else st
    g = st.session_state.get
    try:
        with model_config_context(snapshot):
            cfg, iB, iP, iS, mk = build_nominal_config_from_gui()
            mode = g("int_builder_mode", "Direct (ModelBuilder)")
            gfn = g("int_growth_function", "monod_growth")
            debris_on = bool(g("int_debris_enabled", False))
            imm_on = bool(g("int_immunity_enabled", False))
            solver = [
                ("Integrator", g("int_solver_method", "BDF")),
                ("Extinction threshold", g("int_extinction_threshold", 1.0)),
                ("Extinction check interval (h)", g("int_extinction_check_interval", 0.0)),
                ("Stationary pre-run (h)", g("int_t_prerun", 0.0)),
                ("Latent compartments (L)", g("int_n_latent", 5)),
            ]
    except Exception as e:  # noqa: BLE001 — snapshot must never crash the page
        c.warning(f"Could not build the model config for the snapshot: {e}")
        return

    def _ga(*names):
        for n in names:
            v = getattr(cfg, n, None)
            if v is not None:
                return v
        return None

    def _section(title, pairs):
        rows = [(k, _snap_fmt(v)) for k, v in pairs if v is not None]
        if rows:
            c.markdown(f"**{title}**")
            c.dataframe(pd.DataFrame(rows, columns=["parameter", "value"]),
                        hide_index=True, width="stretch")

    _na = getattr(cfg, "n_antibiotics", 0) or 0
    c.markdown(f"**Builder mode:** {mode}  ·  **{getattr(cfg, 'n_bacteria', len(np.atleast_1d(iB)))}** "
               f"strain(s) · **{getattr(cfg, 'n_phages', len(np.atleast_1d(iP)))}** phage(s)"
               + (f" · **{_na}** antibiotic(s)" if _na else ""))
    _section("Growth & nutrient environment", [
        ("Growth signal", gfn),
        ("Growth rates (h⁻¹)", _ga("growth_rates")),
        ("Monod constant Ks", _ga("monod_constant") if gfn != "gompertz_growth" else None),
        ("Carrying capacity K", _ga("carrying_capacity") if gfn not in ("gompertz_growth", "sequential_monod") else None),
        ("Gompertz shape k", _ga("monod_constant") if gfn == "gompertz_growth" else None),
        ("Gompertz inflection S∞", _ga("carrying_capacity") if gfn == "gompertz_growth" else None),
        ("Inefficient Monod K (low S)", _ga("monod_K_low") if gfn == "smooth_efficiency_monod" else None),
        ("Efficiency transition θ", _ga("monod_efficiency_theta") if gfn == "smooth_efficiency_monod" else None),
        ("Efficiency transition Hill", _ga("monod_efficiency_hill") if gfn == "smooth_efficiency_monod" else None),
        ("Density throttle Kd", _ga("density_growth_constant") if gfn == "density_throttled_growth" else None),
        ("Diauxic rate factors", _ga("growth_phase_rate_factors") if gfn == "sequential_monod" else None),
        ("Diauxic Monod Kᵢ", _ga("growth_phase_monod") if gfn == "sequential_monod" else None),
        ("Diauxic thresholds θ", _ga("growth_phase_thresholds") if gfn == "sequential_monod" else None),
        ("Bacteria : resource ratio", _ga("bacteria_to_resource_ratio")),
        ("Track nutrients", g("int_track_nutrients", True)),
        ("Initial nutrient S₀", iS),
        ("Recycle fraction", _ga("recycle_fraction")),
        ("Nutrient inflow s_in", _ga("s_in")),
        ("Nutrient washout s_out", _ga("s_out")),
        ("Infected-cell nutrient draw (× per phage)",
         _ga("infected_nutrient_consumption") if np.any(_ga("infected_nutrient_consumption")) else None),
    ])
    _section("Death & dormancy", [
        ("Natural death dB (h⁻¹)", _ga("death_rate_B")),
        ("Dormant death dD (h⁻¹)", _ga("death_rate_D")),
        ("Dormancy function", _ga("dormancy_function")),
        ("Dormancy rate (h⁻¹)", _ga("dormancy_rate")),
        ("Resuscitation rate (h⁻¹)", _ga("resuscitation_rate")),
        ("Depth diffusion (h⁻¹)", _ga("dormancy_diffusion_rate")),
        ("Dormancy density threshold", _ga("dormancy_carrying_capacity")),
        ("Depth layers (Q)", _ga("n_depth")),
    ])
    _section("Phage", [
        ("Adsorption (mL·h⁻¹)", _ga("adsorption_rates")),
        ("Adsorption → dormant", _ga("adsorption_rates_dormant")),
        ("Burst size", _ga("burst_sizes")),
        ("Latent period (h)", _ga("latent_periods")),
        ("Phage decay (h⁻¹)", _ga("phage_decay_rates")),
        ("Dormant attenuation", _ga("attenuation_rate")),
    ])
    def _nz(name):   # generator matrices: show only when some rate is set
        v = _ga(name)
        return v if v is not None and np.any(np.asarray(v, dtype=float) != 0) else None
    _section("Evolution & phenotypic switching (generator matrices, M[dest, origin])", [
        ("Bacterial mutation (per division)", _nz("mutation_rates")),
        ("Bacterial phenotypic transitions (h⁻¹)", _nz("transition_rates")),
        ("Phage mutation (per burst)", _nz("mutation_rates_phage")),
        ("Phage phenotypic transitions (h⁻¹)", _nz("transition_rates_phage")),
    ])
    if imm_on:
        _section("Host immunity", [
            ("Module", g("int_immune_module", "innate")),
            ("Stim rate", _ga("imm_stim_rate")), ("Stim50", _ga("imm_stim50")),
            ("Kill rate", _ga("imm_kill_rate")), ("Kill50", _ga("imm_kill50")),
            ("Imm max (hill)", _ga("imm_max")), ("Decay rate", _ga("imm_decay_rate")),
            ("Dormant kill rate", _ga("imm_kill_rate_D")),
        ])
    if debris_on:
        _section("OD / debris", [
            ("od_to_cfu (CFU per OD)", _ga("od_to_cfu_conversion_factor")),
            ("Debris · deaths (u)", _ga("debris_u")), ("Debris · lysis (v)", _ga("debris_v")),
            ("Dissolution k_dis", _ga("debris_kdis")),
            ("Dormant OD weight", _ga("dormant_od_fraction")),
        ])
    _section("Initial conditions", [
        ("Initial bacteria B₀", iB), ("Initial phage P₀", iP),
        ("Initial dormant D₀", mk.get("initial_D") if isinstance(mk, dict) else None),
    ])
    _section("Solver & structure", solver)


def resolve_model_snapshot(selection: str) -> dict | None:
    """Map a Model selector value to a frozen snapshot (or None = live Working draft).

    `selection` is a registry key: WORKING_DRAFT_LABEL, a demo name, or a user model
    name. Demos are materialised on the default model + overrides."""
    if selection in (None, WORKING_DRAFT_LABEL):
        return None
    demos = {d["name"]: d for d in DEMO_MODELS}
    if selection in demos:
        base = copy.deepcopy(st.session_state.get("_default_model_state") or {})
        base.update(copy.deepcopy(demos[selection]["overrides"]))
        return base
    saved = st.session_state.get("user_models", {})
    if selection in saved:
        return copy.deepcopy(saved[selection]["state"])
    return None


def model_options() -> list[str]:
    """Ordered Model selector options: live draft, prebuilt demos, then saved models."""
    return ([WORKING_DRAFT_LABEL]
            + [d["name"] for d in DEMO_MODELS]
            + list(st.session_state.get("user_models", {}).keys()))


def page_model_selector(page_key: str, label: str = "Model") -> str:
    """Render a per-page Model selector (defaults to the active model) and return the
    selection. Wrap the page body in ``model_config_context(resolve_model_snapshot(sel))``
    so the whole task runs against the chosen Model, not the live builder widgets."""
    opts = model_options()
    cur = st.session_state.active_model if st.session_state.active_model in opts else WORKING_DRAFT_LABEL
    key = f"{page_key}_model_sel"
    if st.session_state.get(key) not in opts:
        st.session_state.pop(key, None)
    sel = st.selectbox(
        label, opts, index=opts.index(cur), key=key,
        help="Which Model this task runs against. 'Working draft (live)' uses the current "
             "Interactive-Simulator builder state; saved/demo models are frozen snapshots, "
             "immune to builder edits.")
    if sel == WORKING_DRAFT_LABEL:
        st.caption("Running against the **live builder draft** — edits in the Interactive "
                   "Simulator flow through here.")
    else:
        st.caption(f"Running against frozen model **{sel}** — builder edits won't affect this.")
    # Click-to-view full config of the selected Model (no tab-hunting). The toggle gates
    # the (moderately expensive) config build so it only runs when the user asks.
    if st.toggle("📋 Show model config", key=f"{page_key}_show_cfg",
                 help="A full snapshot of the selected model's parameters — strains, phages, "
                      "growth & nutrient environment, immunity, OD/debris, ICs, solver."):
        render_model_snapshot(snapshot=resolve_model_snapshot(sel))
    return sel


WORKING_DRAFT_LABEL = "Working draft (live)"


# ── Prebuilt demo models ───────────────────────────────────────────────────────
# Small, curated set spanning distinct scenarios. Each is defined as overrides on
# top of the default model (captured once at init as `_default_model_state`), so the
# specs stay short and readable. Expand slowly as new representative cases are needed.
DEMO_MODELS = [
    {
        "name": "Growth calibration (Monod)",
        "description": "Single strain, nutrient-limited (Monod) growth to a plateau — "
                       "the base model that fits the demonstration dataset "
                       "(examples/tutorial_synthetic_brg.csv).",
        "overrides": {
            "int_builder_mode": "Direct (ModelBuilder)",
            "int_strains": [{
                "name": "Strain 0 (WT)", "initial_B": 5e6, "growth_rate": 1.2,
                "bacteria_to_resource_ratio": 1e8, "death_rate_B": 0.0,
                "dormancy_enabled": False,
            }],
            "int_phages": [{
                "name": "Phage 0", "initial_P": 0.0, "adsorption_rates": 1e-8,
                "adsorption_rates_dormant": 0.0, "burst_sizes": 50.0,
                "latent_periods": 0.5, "phage_decay_rates": 0.1, "pk_mode": "None",
            }],
            "int_antibiotics": [],
            "int_track_nutrients": True,
            "int_growth_function": "monod_growth",
            "int_lysis_function": "constant_lysis",
            "int_monod_constant_lysis": 0.3,
            "int_lysis_floor": 0.0,  # phi_min: residual lysis efficacy at stationary phase
            "int_dormancy_signal": "nutrient",
            "int_resuscitation_signal": "nutrient",
            "int_diffusion_signal": "constant",
            "int_dormancy_monod_constant": 0.0,
            "int_dormancy_carrying_capacity": 1e8,
            "int_monod_constant": 0.3,
            "int_density_growth_constant": 1e9,
            "int_growth_n_phases": 2,  # diauxic (sequential_monod) phase count
            "int_gompertz_k": 10.0,    # gompertz_growth shape k (nutrient scale)
            "int_gompertz_sinf": 0.5,  # gompertz_growth inflection S∞ (nutrient scale)
            "int_monod_K_low": 3.0,            # smooth_efficiency_monod: inefficient low-S K
            "int_monod_efficiency_theta": 0.5, # smooth_efficiency_monod: transition midpoint θ
            "int_monod_efficiency_hill": 4.0,  # smooth_efficiency_monod: transition sharpness
            "int_infected_nutrient_consumption": 0.0,  # infected-cell nutrient draw (B, task B)
            "int_recycle_fraction": 0.5,
            "int_initial_S": 1.0,
            "int_od_to_cfu_conversion_factor": 2e8,
            "int_od_nutrient_enabled": False,
            "int_od_abs_exp": 1.0,
            "int_od_abs_stat": 0.4,
            "int_od_abs_S50": 0.3,
            "int_od_abs_hill": 2.0,
            "int_debris_enabled": True,
            "int_debris_u": 0.4,
            "int_debris_v": 0.2,
            "int_debris_kdis": 0.01,
            "int_immunity_enabled": False,
        },
    },
    {
        "name": "Two-strain resistance (WT + resistant)",
        "description": "WT plus a phage-resistant strain (Direct mode), one phage. A good "
                       "model to try reparameterization on — e.g. tie the two strains' "
                       "growth rates to a single estimated value.",
        "overrides": {
            "int_builder_mode": "Direct (ModelBuilder)",
            "int_strains": [
                {"name": "WT", "initial_B": 1e7, "growth_rate": 1.2,
                 "bacteria_to_resource_ratio": 1e8, "death_rate_B": 0.0, "dormancy_enabled": False},
                {"name": "Resistant", "initial_B": 10.0, "growth_rate": 1.1,
                 "bacteria_to_resource_ratio": 1e8, "death_rate_B": 0.0, "dormancy_enabled": False},
            ],
            "int_phages": [{
                "name": "Phage 0", "initial_P": 1e6, "adsorption_rates": 1e-8,
                "adsorption_rates_dormant": 0.0, "burst_sizes": 50.0,
                "latent_periods": 0.5, "phage_decay_rates": 0.1, "pk_mode": "None",
            }],
            "int_antibiotics": [],
            "int_track_nutrients": True,
            "int_growth_function": "monod_growth",
            "int_lysis_function": "constant_lysis",
            "int_monod_constant_lysis": 0.3,
            "int_lysis_floor": 0.0,  # phi_min: residual lysis efficacy at stationary phase
            "int_dormancy_signal": "nutrient",
            "int_resuscitation_signal": "nutrient",
            "int_diffusion_signal": "constant",
            "int_dormancy_monod_constant": 0.0,
            "int_dormancy_carrying_capacity": 1e8,
            "int_monod_constant": 0.3,
            "int_density_growth_constant": 1e9,
            "int_growth_n_phases": 2,  # diauxic (sequential_monod) phase count
            "int_gompertz_k": 10.0,    # gompertz_growth shape k (nutrient scale)
            "int_gompertz_sinf": 0.5,  # gompertz_growth inflection S∞ (nutrient scale)
            "int_monod_K_low": 3.0,            # smooth_efficiency_monod: inefficient low-S K
            "int_monod_efficiency_theta": 0.5, # smooth_efficiency_monod: transition midpoint θ
            "int_monod_efficiency_hill": 4.0,  # smooth_efficiency_monod: transition sharpness
            "int_infected_nutrient_consumption": 0.0,  # infected-cell nutrient draw (B, task B)
            "int_recycle_fraction": 0.5,
            "int_initial_S": 1.0,
            "int_od_to_cfu_conversion_factor": 2e8,
            "int_od_nutrient_enabled": False,
            "int_od_abs_exp": 1.0,
            "int_od_abs_stat": 0.4,
            "int_od_abs_S50": 0.3,
            "int_od_abs_hill": 2.0,
            "int_debris_enabled": True,
            "int_debris_u": 0.4,
            "int_debris_v": 0.2,
            "int_debris_kdis": 0.01,
            "int_immunity_enabled": False,
        },
    },
]


# ── Parts library (Tier 2: composable bacteria / phages / antibiotics) ─────────
# A "part" is one reusable entity (a bacterium, a phage, or an antibiotic) — its
# parameter dict plus provenance metadata. Parts compose into scenarios: loading a
# part appends its dict to the shared entity list every module reads
# (int_strains / int_phages / int_antibiotics). Phage kinetics (burst/latent/
# adsorption) are phage×host properties, not phage-intrinsic, so phage parts carry
# a `reference_host` tag and a soft "verify for this strain" flag on mismatch.
PARTS_SCHEMA_VERSION = 1


# category -> (shared session entity-list key, max count in the UI, human label)
PART_CATEGORIES = {
    "bacteria":    {"key": "int_strains",     "max": 10, "label": "Bacteria"},
    "phages":      {"key": "int_phages",      "max": 10, "label": "Phages"},
    "antibiotics": {"key": "int_antibiotics", "max": 6, "label": "Antibiotics"},
}


PART_SOURCES = ["educated guess", "literature", "pbisim-fit", "experimental"]


def empty_parts_library() -> dict:
    return {cat: {} for cat in PART_CATEGORIES}


def export_parts_json(library: dict) -> str:
    return json.dumps(
        {"schema_version": PARTS_SCHEMA_VERSION, "parts": _json_safe(library)},
        indent=2,
    )


def import_parts_json(text: str) -> dict:
    """Parse an exported parts library into the {category: {name: part}} mapping."""
    data = json.loads(text)
    parts = data.get("parts", data) if isinstance(data, dict) else None
    if not isinstance(parts, dict):
        raise ValueError("Not a parts library (expected a JSON object).")
    out = empty_parts_library()
    for cat, entries in parts.items():
        if cat not in PART_CATEGORIES:
            continue  # ignore unknown categories rather than fail the whole import
        if not isinstance(entries, dict):
            raise ValueError(f"Parts category '{cat}' must be an object.")
        for name, part in entries.items():
            if not isinstance(part, dict) or "params" not in part:
                raise ValueError(f"Part '{cat}/{name}' is missing its 'params'.")
            out[cat][name] = part
    return out


# Entity-widget key prefixes (NOT data keys). Cleared when a part is appended so the
# strain/phage/antibiotic widgets re-read from the (updated) int_* data lists rather
# than shadowing them with stale per-index widget values.
_ENTITY_WIDGET_PREFIXES = (
    "str_", "phg_", "ss_", "brg_", "abx_", "ads_input_", "ads_dorm_input_",
)


def clear_entity_widgets() -> None:
    for k in list(st.session_state.keys()):
        if any(k.startswith(p) for p in _ENTITY_WIDGET_PREFIXES):
            st.session_state.pop(k, None)


# Minimal sensible starting configuration (one WT strain + one phage, Monod
# nutrients). Used to populate a fresh session and the "Reset" action. Kept in the
# app (not tied to the pbisim tutorials, which may change independently).
DEFAULT_SCENARIO = {
    "t_end": 48.0,
    "dt": 0.25,
    "extinction_threshold": 1.0,
    "solver_method": "BDF",
    "track_nutrients": True,
    "initial_S": 1.0,
    "monod_constant": 0.3,
    "recycle_fraction": 0.0,
    "s_in": 0.0,
    "s_out": 0.0,
    "immunity_enabled": False,
    "debris_enabled": False,
    "strains": [
        {
            "name": "Strain 0 (WT)",
            "initial_B": 1e7,
            "growth_rate": 1.2,
            "bacteria_to_resource_ratio": 1e9,
            "death_rate_B": 0.0,
            "dormancy_enabled": False,
        }
    ],
    "phages": [
        {
            "name": "Phage 0",
            "initial_P": 1e6,
            "adsorption_rates": 1e-8,
            "adsorption_rates_dormant": 0.0,
            "burst_sizes": 50.0,
            "latent_periods": 0.5,
            "phage_decay_rates": 0.1,
            "pk_mode": "None",
        }
    ],
    "antibiotics": [],
    "doses": [],
}


# Dormancy / resuscitation entry-signal options — must match pbisim's
# _DORMANCY_SIGNALS / _RESUSCITATION_SIGNALS keys exactly.
SIGNAL_OPTIONS = ["constant", "nutrient", "density", "nutrient+density"]


# Growth-signal options → pbisim growth function name + whether nutrients are tracked.
GROWTH_SIGNALS = {
    "nutrient (Monod)":            ("monod_growth", True),
    "nutrient + density":          ("monod_logistic_growth", True),
    "density (logistic)":          ("logistic_growth", False),
    "constant (unlimited)":        ("constant_growth", False),
    # Newer pbisim growth functions (all nutrient-tracking):
    "density-throttled (Monod × 1/(1+ΣB/Kd))": ("density_throttled_growth", True),
    "Gompertz (nutrient)":         ("gompertz_growth", True),
    "sequential / diauxic (Monod)": ("sequential_monod", True),
    "smooth two-efficiency Monod": ("smooth_efficiency_monod", True),
}

# Growth functions that read nutrient S (need track_nutrients) and that need a carrying
# capacity / S∞. gompertz_growth REUSES monod_constant as its shape k and carrying_capacity
# as its inflection S∞ (a pbisim convention), so it appears in both sets.
_GROWTH_NUTRIENT = {"monod_growth", "monod_logistic_growth", "density_throttled_growth",
                    "gompertz_growth", "sequential_monod", "smooth_efficiency_monod"}
_GROWTH_NEEDS_K = {"logistic_growth", "monod_logistic_growth", "gompertz_growth"}
_GROWTH_NEEDS_KD = {"density_throttled_growth"}
# sequential_monod (diauxic) needs the per-phase arrays instead of a single Monod K.
_GROWTH_SEQUENTIAL = "sequential_monod"
# smooth_efficiency_monod: efficient Monod K at high S blends to an inefficient K_low
# at low S (Hill of S). Uses monod_constant as the high-S K plus its own low-S K + θ + hill.
_GROWTH_SMOOTH_EFF = "smooth_efficiency_monod"


def growth_tracks_nutrients() -> bool:
    """True when the selected growth signal reads the nutrient substrate S — i.e. any
    nutrient-tracking growth function, not just Monod. This is the correct gate for the
    nutrient-COUPLED signals (frac_lysis, nutrient dormancy/resuscitation/diffusion): they
    need a live S, which every function in ``_GROWTH_NUTRIENT`` provides (Monod, Monod×
    logistic, density-throttled, Gompertz, diauxic, smooth-efficiency)."""
    return st.session_state.get("int_growth_function", "monod_growth") in _GROWTH_NUTRIENT


# Death-signal options → pbisim death function name. constant_death (default) is the
# flat rate d that the app used all along; nutrient = starvation d·(1−S/(Ks+S));
# density = crowding d·min(1, ΣB/K). (No nutrient+density death function exists.)
DEATH_SIGNALS = {
    "constant":                     "constant_death",
    "nutrient (starvation)":        "nutrient_dependent_death",
    "density (crowding)":           "density_dependent_death",
    "nutrient + density":           "nutrient_and_density_death",
}


# Lysis-signal options → pbisim lysis-progression function name. `constant_lysis` (the
# engine default) makes the latent→lysis chain progress at a fixed rate (phi_lysis = 1.0);
# `frac_lysis` couples it to nutrients (phi_lysis = S/(Ks_lysis + S)), slowing lysis when
# nutrients are scarce. frac_lysis reads S, so it needs a nutrient-tracking growth signal.
LYSIS_SIGNALS = {
    "constant":            "constant_lysis",
    "nutrient (Monod)":    "frac_lysis",
}


def lysis_signal_function(name):
    """pbisim lysis-progression function object for a stored name, or ``None`` (the engine
    default, constant_lysis) — ModelConfig treats ``lysis_progression_function=None`` as
    constant lysis, so we don't import a constant_lysis object."""
    if name == "frac_lysis":
        from pbisim import frac_lysis
        return frac_lysis
    return None


def lysis_kwargs():
    """Lysis-progression config for the selected lysis signal, as ModelConfig fields
    (``lysis_progression_function`` + ``monod_constant_lysis`` + ``lysis_floor``). frac_lysis
    is nutrient-coupled, so it is coerced to constant when nutrients aren't tracked (same
    guard as the dormancy signals). ``lysis_floor`` (phi_min ∈ [0,1]) rescales the signal
    ``phi = floor + (1-floor)*signal`` so phage keep some efficacy at stationary phase; only
    meaningful for frac_lysis (a no-op for constant lysis, whose signal is already 1).
    Returns ``{lysis_progression_function, monod_constant_lysis, lysis_floor, coerced}``."""
    name = st.session_state.get("int_lysis_function", "constant_lysis")
    track = growth_tracks_nutrients()   # any nutrient-tracking growth, not just Monod
    coerced = False
    if name == "frac_lysis" and not track:
        name, coerced = "constant_lysis", True
    fn = lysis_signal_function(name)
    ks = st.session_state.get("int_monod_constant_lysis", 0.3) if name == "frac_lysis" else None
    floor = float(st.session_state.get("int_lysis_floor", 0.0) or 0.0) if name == "frac_lysis" else 0.0
    return {"lysis_progression_function": fn, "monod_constant_lysis": ks,
            "lysis_floor": floor, "coerced": coerced}


def canonical_signal(v):
    """Normalise a stored dormancy/resuscitation signal to a pbisim-recognised key.

    Translates the app's legacy ``'nutrient_and_density'`` to the engine's
    ``'nutrient+density'`` (the mismatch that made that option raise a ValueError).
    """
    v = v or "nutrient"
    if v in ("nutrient_and_density", "nutrient+density"):
        return "nutrient+density"
    return v if v in SIGNAL_OPTIONS else "nutrient"


def compat_dormancy_signal(sig, track_nutrients):
    """Coerce a nutrient-based dormancy/resuscitation signal to a nutrient-independent
    one when nutrients are not tracked (S frozen) — otherwise the engine refuses to
    build (a nutrient signal can't read a frozen S). Returns ``(signal, coerced?)``.
    """
    if track_nutrients:
        return sig, False
    if sig == "nutrient":
        return "constant", True
    if sig == "nutrient+density":
        return "density", True
    return sig, False  # constant / density are already compatible


def sequential_growth_phase_params():
    """Read the diauxic (``sequential_monod``) per-phase parameters from session state
    and return ``(rate_factors, monod_constants, thresholds)`` as float ``np.ndarray``.

    ``rate_factors`` / ``monod_constants`` have length X (= ``int_growth_n_phases``);
    ``thresholds`` has length X-1. Phase 1's rate factor is pinned to 1.0 (the engine
    treats phase 1 as the reference). Returns engine-ready arrays; the caller validates
    (see :func:`validate_sequential_growth`).
    """
    n = int(st.session_state.get("int_growth_n_phases", 2) or 2)
    n = max(2, n)
    rf = [1.0]
    mc = [float(st.session_state.get("gp_monod_0", st.session_state.get("int_monod_constant", 0.3)))]
    for i in range(1, n):
        rf.append(float(st.session_state.get(f"gp_rate_{i}", 0.5)))
        mc.append(float(st.session_state.get(f"gp_monod_{i}", 0.3)))
    th = [float(st.session_state.get(f"gp_thresh_{i}", 0.3 / (i + 1))) for i in range(n - 1)]
    return np.asarray(rf, dtype=float), np.asarray(mc, dtype=float), np.asarray(th, dtype=float)


def validate_sequential_growth(rate_factors, monod_constants, thresholds):
    """Return an error string if the diauxic phase parameters are invalid (mirrors the
    engine's ``with_sequential_growth`` checks), else ``None``."""
    rf, mc, th = np.asarray(rate_factors), np.asarray(monod_constants), np.asarray(thresholds)
    if len(rf) != len(mc):
        return "Rate factors and Monod constants must have the same number of phases."
    if len(th) != len(rf) - 1:
        return "Thresholds must have one fewer entry than the phases (X−1)."
    if len(th) and (not np.all(np.diff(th) < 0) or np.any(th <= 0) or np.any(th >= 1)):
        return "Thresholds must be strictly decreasing and lie in (0, 1)."
    return None


def growth_nutrient_kwargs():
    """Growth function + nutrient config for the selected growth signal.

    Returns a dict of ModelConfig fields (growth_function, track_nutrients, and the
    relevant monod_constant / carrying_capacity / recycle / s_in / s_out). Works both
    for the Direct builder (split into with_growth_function + with_nutrient) and for
    BRG / StrainSet ``to_config(**extra_config_kwargs)`` (forwarded to ModelConfig).
    For ``sequential_monod`` the three per-phase arrays are included as ModelConfig
    fields (BRG/StrainSet forward them verbatim; the Direct build instead calls
    ``with_sequential_growth`` — see ``build_nominal_config_from_gui``).
    """
    from pbisim import (monod_growth, logistic_growth, constant_growth, monod_logistic_growth,
                        density_throttled_growth, gompertz_growth, sequential_monod,
                        smooth_efficiency_monod)
    fns = {"monod_growth": monod_growth, "logistic_growth": logistic_growth,
           "constant_growth": constant_growth, "monod_logistic_growth": monod_logistic_growth,
           "density_throttled_growth": density_throttled_growth, "gompertz_growth": gompertz_growth,
           "sequential_monod": sequential_monod, "smooth_efficiency_monod": smooth_efficiency_monod}
    name = st.session_state.get("int_growth_function", "monod_growth")
    fn = fns.get(name, monod_growth)
    nutrient_based = name in _GROWTH_NUTRIENT
    needs_K = name in _GROWTH_NEEDS_K
    # monod_constant + recycle_fraction are always supplied — StrainSet.to_config
    # requires them (they are simply unused by non-nutrient growth functions).
    kw = {
        "growth_function": fn,
        "track_nutrients": nutrient_based,
        "monod_constant": st.session_state.get("int_monod_constant", 0.3),
        "recycle_fraction": st.session_state.get("int_recycle_fraction", 0.0),
        # density signals (dormancy/resuscitation/death) count active B, or all cell
        # states (B+I+D+H) when this is on. Forwarded to with_nutrient (Direct) and
        # to_config (BRG/StrainSet), and captured by the repro recorder automatically.
        "density_signal_uses_total_cells": st.session_state.get("int_density_total_cells", False),
    }
    if nutrient_based:
        kw["s_in"] = st.session_state.get("int_s_in", 0.0)
        kw["s_out"] = st.session_state.get("int_s_out", 0.0)
        # Infected (latent I) cells draw down the shared substrate while building phage, at
        # a multiple of the uninfected per-capita uptake (0 = legacy off; >1 = a hijacked
        # cell consumes more). It's a PER-PHAGE (phage×host) property in the engine, so
        # assemble the (n_phages,) array from the phage dicts. Always an array when phages
        # exist so config.infected_nutrient_consumption is sweepable per phage; a legacy
        # scalar (int_infected_nutrient_consumption) seeds any phage that hasn't set its own.
        _phages = st.session_state.get("int_phages", [])
        _legacy = float(st.session_state.get("int_infected_nutrient_consumption", 0.0) or 0.0)
        if _phages:
            kw["infected_nutrient_consumption"] = np.array(
                [float(p.get("infected_nutrient_consumption", _legacy)) for p in _phages])
        else:
            kw["infected_nutrient_consumption"] = _legacy
    if needs_K:
        kw["carrying_capacity"] = st.session_state.get("int_carrying_capacity", 1e9)
    if name == "gompertz_growth":
        # gompertz_growth = exp(−exp(−k·(S−S∞))) — a NUTRIENT-scale sigmoid. The engine
        # reads k from monod_constant and S∞ from carrying_capacity, but those otherwise
        # hold the Monod half-saturation and the density ceiling (~1e9). Reusing the
        # density K as S∞ makes the exponent overflow → growth ≈ 0 (a flat curve), so
        # Gompertz gets its OWN nutrient-scale k / S∞ keys, mapped in here.
        kw["monod_constant"] = st.session_state.get("int_gompertz_k", 10.0)
        kw["carrying_capacity"] = st.session_state.get("int_gompertz_sinf", 0.5)
    if name in _GROWTH_NEEDS_KD:  # density_throttled_growth's hyperbolic throttle constant
        kw["density_growth_constant"] = st.session_state.get("int_density_growth_constant", 1e9)
    if name == _GROWTH_SMOOTH_EFF:  # smooth two-efficiency Monod: low-S K + Hill transition
        # monod_constant (above) is the efficient high-S K; add the inefficient low-S K
        # and the S-transition (θ midpoint, hill sharpness).
        kw["monod_K_low"] = st.session_state.get("int_monod_K_low", 3.0)
        kw["monod_efficiency_theta"] = st.session_state.get("int_monod_efficiency_theta", 0.5)
        kw["monod_efficiency_hill"] = st.session_state.get("int_monod_efficiency_hill", 4.0)
    if name == _GROWTH_SEQUENTIAL:  # diauxic per-phase arrays
        rf, mc, th = sequential_growth_phase_params()
        kw["growth_phase_rate_factors"] = rf
        kw["growth_phase_monod"] = mc
        kw["growth_phase_thresholds"] = th
    return kw


def dormancy_signal_functions(dsig, rsig):
    """Map dormancy/resuscitation signal strings to pbisim function objects.

    BRG / StrainSet ``to_config`` take function objects (``dormancy_function`` /
    ``resuscitation_function``) rather than the Direct builder's signal strings.
    (Imported from the transitions module — ``nutrient_and_density_resuscitation``
    isn't re-exported at the pbisim top level.)
    """
    from pbisim.dormancy.transitions import (
        nutrient_dependent_dormancy, constant_dormancy, density_dependent_dormancy,
        nutrient_and_density_dormancy, nutrient_dependent_resuscitation,
        constant_resuscitation, density_dependent_resuscitation,
        nutrient_and_density_resuscitation,
    )
    _D = {"constant": constant_dormancy, "nutrient": nutrient_dependent_dormancy,
          "density": density_dependent_dormancy, "nutrient+density": nutrient_and_density_dormancy}
    _R = {"constant": constant_resuscitation, "nutrient": nutrient_dependent_resuscitation,
          "density": density_dependent_resuscitation, "nutrient+density": nutrient_and_density_resuscitation}
    return _D.get(dsig, nutrient_dependent_dormancy), _R.get(rsig, nutrient_dependent_resuscitation)


def diffusion_signal_functions(sig):
    """Map a dormancy-depth diffusion signal string to the (deeper, shallower) pbisim
    function pair (engine ``_DIFFUSION_SIGNALS``). ``constant`` = legacy symmetric,
    nutrient-independent diffusion."""
    from pbisim.dormancy.transitions import (
        constant_diffusion, nutrient_dependent_diffusion_deeper,
        nutrient_dependent_diffusion_shallower, density_dependent_diffusion_deeper,
        density_dependent_diffusion_shallower, nutrient_and_density_diffusion_deeper,
        nutrient_and_density_diffusion_shallower,
    )
    return {
        "constant": (constant_diffusion, constant_diffusion),
        "nutrient": (nutrient_dependent_diffusion_deeper, nutrient_dependent_diffusion_shallower),
        "density": (density_dependent_diffusion_deeper, density_dependent_diffusion_shallower),
        "nutrient+density": (nutrient_and_density_diffusion_deeper, nutrient_and_density_diffusion_shallower),
    }.get(sig, (constant_diffusion, constant_diffusion))


def set_diffusion_functions(config, sig, rec=None):
    """Set a config's dormancy-depth diffusion functions from a signal string (with the
    same nutrient-tracking coercion as Direct mode). BRG / StrainSet ``to_config`` don't
    expose them, so we set the ModelConfig fields directly, post-build. When a repro
    recorder ``rec`` is given, the assignment is mirrored into the generated script."""
    track = growth_tracks_nutrients()   # any nutrient-tracking growth, not just Monod
    sig, _ = compat_dormancy_signal(canonical_signal(sig), track)
    deeper, shallower = diffusion_signal_functions(sig)
    config.dormancy_diffusion_deeper_function = deeper
    config.dormancy_diffusion_shallower_function = shallower
    if rec is not None:
        rec.diffusion_functions(sig)
    return config


def apply_generator_matrices(config, bact_names, phage_names, rec=None):
    """Set the bacterial phenotypic-transition and phage mutation / transition generator
    matrices directly on a built ``ModelConfig`` (BRG's ``to_config`` hard-codes zeros
    for them). Fields whose graph has no valid edge are left untouched. Mirrored into
    the reproduction script via ``rec`` as plain ``cfg.<field> = np.array(...)`` lines."""
    for fld, M in generator_matrices_from_state(bact_names, phage_names).items():
        if M is None:
            continue
        setattr(config, fld, M)
        if rec is not None:
            rec.assign("cfg", fld, M)
    return config


def apply_diffusion_signal(config, strains, rec=None):
    """StrainSet/Direct: set the depth-diffusion functions from the first dormancy-enabled
    strain's ``diffusion_signal``. No-op when no strain has dormancy on."""
    enabled = [s for s in strains if s.get("dormancy_enabled", False)]
    if not enabled:
        return config
    return set_diffusion_functions(config, enabled[0].get("diffusion_signal", "constant"), rec=rec)


def mode_dormancy_kwargs(dsig="nutrient", rsig="nutrient", ks=0.0, kdorm=0.0):
    """``to_config`` dormancy kwargs (function objects + Ks / K_dorm) for BRG / StrainSet
    from the selected signal strings, applying the same compatibility coercion as Direct
    mode (a nutrient signal needs nutrient-tracking growth; density needs a threshold).
    Also covers the dormancy-disabled case — the coerced default is nutrient-independent
    when the growth signal freezes S, so the engine's nutrient default can't crash it.
    """
    track = growth_tracks_nutrients()   # any nutrient-tracking growth, not just Monod
    ds, _ = compat_dormancy_signal(canonical_signal(dsig), track)
    rs, _ = compat_dormancy_signal(canonical_signal(rsig), track)
    dfn, rfn = dormancy_signal_functions(ds, rs)
    kw = {"dormancy_function": dfn, "resuscitation_function": rfn}
    if ks > 0 and any(s in ("nutrient", "nutrient+density") for s in (ds, rs)):
        kw["dormancy_monod_constant"] = ks
    if any(s in ("density", "nutrient+density") for s in (ds, rs)):
        kw["dormancy_carrying_capacity"] = kdorm if kdorm > 0 else st.session_state.get("int_carrying_capacity", 1e9)
    return kw


def death_signal_function(name):
    """Map a death-function name to the pbisim function object."""
    from pbisim import (constant_death, nutrient_dependent_death,
                        density_dependent_death, nutrient_and_density_death)
    return {"constant_death": constant_death,
            "nutrient_dependent_death": nutrient_dependent_death,
            "density_dependent_death": density_dependent_death,
            "nutrient_and_density_death": nutrient_and_density_death}.get(name, constant_death)


def death_kwargs():
    """death_function (+ carrying_capacity when the density death function needs it) for
    the selected death signal. `constant_death` reproduces the previous behaviour (a flat
    rate d, applied regardless of nutrients)."""
    name = st.session_state.get("int_death_function", "constant_death")
    kw = {"death_function": death_signal_function(name)}
    if name in ("density_dependent_death", "nutrient_and_density_death"):
        kw["carrying_capacity"] = st.session_state.get("int_carrying_capacity", 1e9)
    return kw


class _ReproRecorder:
    """Records the exact builder calls made while constructing the config, so the
    reproduction script is a *byproduct of the real build* rather than a parallel
    re-implementation that can drift.

    Every method executes the real call **and** appends the rendered source line;
    arguments (numpy arrays, pbisim signal functions, PhageStrain / Antibiotic /
    DoseSchedule objects, …) are rendered from the same values the build passes, so
    adding a parameter to the build flows into the script automatically. The rendered
    grammar mirrors the builder the user actually chose (ModelBuilder / BRG / StrainSet).
    """

    def __init__(self):
        self.lines = []            # source lines that construct `cfg`
        self.imports = set()       # names to import `from pbisim import ...`
        self.raw_imports = set()   # full `from X import Y` statements (non-top-level)
        self._render_of = {}       # id(obj) -> source string (inline) or variable name

    def diffusion_functions(self, sig):
        """Record the post-build depth-diffusion function assignments. BRG / StrainSet
        ``to_config`` can't take them, so the app sets them on ``cfg`` directly — mirror
        that in the script. ``constant`` is the config default → nothing to emit."""
        names = {
            "nutrient": ("nutrient_dependent_diffusion_deeper", "nutrient_dependent_diffusion_shallower"),
            "density": ("density_dependent_diffusion_deeper", "density_dependent_diffusion_shallower"),
            "nutrient+density": ("nutrient_and_density_diffusion_deeper", "nutrient_and_density_diffusion_shallower"),
        }.get(sig)
        if not names:
            return
        deeper, shallower = names
        self.raw_imports.add("from pbisim.dormancy.transitions import " + ", ".join(sorted({deeper, shallower})))
        self.lines.append(f"cfg.dormancy_diffusion_deeper_function = {deeper}")
        self.lines.append(f"cfg.dormancy_diffusion_shallower_function = {shallower}")

    # ---- value rendering ----------------------------------------------------
    def _fval(self, x):
        x = float(x)
        if np.isinf(x):
            return "np.inf" if x > 0 else "-np.inf"
        if np.isnan(x):
            return "np.nan"
        return repr(x)

    def _pylist(self, lst):
        if isinstance(lst, list):
            return "[" + ", ".join(self._pylist(v) for v in lst) + "]"
        if isinstance(lst, bool):
            return repr(lst)
        if isinstance(lst, float):
            return self._fval(lst)
        return repr(lst)

    def render(self, v):
        if id(v) in self._render_of:
            return self._render_of[id(v)]
        if v is None or isinstance(v, (bool, int, str)):
            return repr(v)
        if isinstance(v, (float, np.floating)):
            return self._fval(v)
        if isinstance(v, np.integer):
            return repr(int(v))
        if isinstance(v, np.bool_):
            return repr(bool(v))
        if isinstance(v, np.ndarray):
            return f"np.array({self._pylist(v.tolist())})"
        if callable(v) and hasattr(v, "__name__"):
            nm = v.__name__
            if nm.isidentifier():
                self.imports.add(nm)
                return nm
            raise ValueError(f"cannot render callable {v!r} (no importable name)")
        if isinstance(v, tuple):
            body = ", ".join(self.render(x) for x in v)
            return f"({body}{',' if len(v) == 1 else ''})"
        if isinstance(v, list):
            return "[" + ", ".join(self.render(x) for x in v) + "]"
        if isinstance(v, dict):
            return "{" + ", ".join(f"{self.render(k)}: {self.render(val)}" for k, val in v.items()) + "}"
        if _dc.is_dataclass(v) and not isinstance(v, type):
            cls = type(v).__name__
            self.imports.add(cls)
            body = ", ".join(f"{f.name}={self.render(getattr(v, f.name))}"
                             for f in _dc.fields(v) if f.init)
            return f"{cls}({body})"
        raise ValueError(f"cannot render value of type {type(v)}: {v!r}")

    def _args(self, args, kwargs):
        return ", ".join([self.render(a) for a in args]
                         + [f"{k}={self.render(v)}" for k, v in kwargs.items()])

    # ---- construction / calls (execute AND record) --------------------------
    def init(self, var, cls, *args, **kwargs):
        """`var = Cls(...)` — a top-level object assigned to a script variable."""
        obj = cls(*args, **kwargs)
        self.imports.add(cls.__name__)
        self.lines.append(f"{var} = {cls.__name__}({self._args(args, kwargs)})")
        self._render_of[id(obj)] = var
        return obj

    def new(self, cls, *args, **kwargs):
        """A nested object rendered inline (e.g. PhageStrain inside a list literal)."""
        obj = cls(*args, **kwargs)
        self.imports.add(cls.__name__)
        self._render_of[id(obj)] = f"{cls.__name__}({self._args(args, kwargs)})"
        return obj

    def var(self, name, value):
        """`name = <value>` — bind a (possibly nested) value to a script variable."""
        self.lines.append(f"{name} = {self.render(value)}")
        self._render_of[id(value)] = name
        return value

    def call(self, var, obj, method, *args, **kwargs):
        """`var = var.method(...)` — a fluent builder call that returns the builder."""
        result = getattr(obj, method)(*args, **kwargs)
        self.lines.append(f"{var} = {var}.{method}({self._args(args, kwargs)})")
        if result is not None:
            self._render_of[id(result)] = var
        return result

    def mutate(self, var, obj, method, *args, **kwargs):
        """`var.method(...)` — an in-place call that returns None (e.g. add_strain)."""
        getattr(obj, method)(*args, **kwargs)
        self.lines.append(f"{var}.{method}({self._args(args, kwargs)})")

    def classcall(self, var, cls, method, *args, **kwargs):
        """`var = Cls.method(...)` — a classmethod constructor (e.g. from_strains)."""
        obj = getattr(cls, method)(*args, **kwargs)
        self.imports.add(cls.__name__)
        self.lines.append(f"{var} = {cls.__name__}.{method}({self._args(args, kwargs)})")
        self._render_of[id(obj)] = var
        return obj

    def result(self, var, srcvar, obj, method, *args, **kwargs):
        """`var = srcvar.method(...)` — terminal call producing the config."""
        res = getattr(obj, method)(*args, **kwargs)
        self.lines.append(f"{var} = {srcvar}.{method}({self._args(args, kwargs)})")
        self._render_of[id(res)] = var
        return res

    def expr(self, value, source):
        """Register ``value`` so it renders as the given source expression (e.g. a
        ``brg.equilibrium_initial_condition(...)`` call) instead of a literal dump."""
        self._render_of[id(value)] = source
        return value

    def assign(self, var, attr, value):
        """`var.attr = <value>` — a post-build field assignment on a config object
        (for fields a builder's ``to_config`` can't take)."""
        self.lines.append(f"{var}.{attr} = {self.render(value)}")


def build_nominal_config_from_gui():
    """
    Constructs and returns the ModelConfig and corresponding state initial values
    based on the selected builder mode in the GUI.

    As it builds, it records the exact builder calls into a ``_ReproRecorder`` stashed
    at ``st.session_state['_repro_rec']`` — the reproduction script is generated from
    that recording, so it can never drift from what is actually simulated.
    """
    rec = _ReproRecorder()
    st.session_state["_repro_rec"] = rec

    builder_mode = st.session_state.get("int_builder_mode", "Direct (ModelBuilder)")
    strains = st.session_state.get("int_strains", [])
    phages = st.session_state.get("int_phages", [])
    antibiotics = st.session_state.get("int_antibiotics", [])
    doses = st.session_state.get("int_doses", [])

    n_bacteria = len(strains)
    n_phages = len(phages)
    n_latent = int(st.session_state.get("int_n_latent", 5))  # latency compartments (all builders)

    # ── Resolve solver settings ───────────────────────────────────────────────
    # track_nutrients is a property of the chosen growth signal (a non-nutrient growth
    # function freezes S), not an independent flag — derive it from int_growth_function
    # so the dormancy-compat coercion below stays consistent with what with_nutrient()
    # actually sets (int_track_nutrients can drift, e.g. when only the growth function
    # is changed, as in a categorical growth-signal sweep).
    track_nutrients = st.session_state.get("int_growth_function", "monod_growth") in _GROWTH_NUTRIENT
    superinfection = st.session_state.get("int_superinfection", False)
    
    # ── Resolve Debris parameters ─────────────────────────────────────────────
    debris_enabled = st.session_state.get("int_debris_enabled", False)
    extra_kwargs = {}
    if debris_enabled:
        extra_kwargs["debris_u"] = st.session_state.get("int_debris_u", 0.4)
        extra_kwargs["debris_v"] = st.session_state.get("int_debris_v", 0.2)
        extra_kwargs["debris_kdis"] = st.session_state.get("int_debris_kdis", 0.01)
        extra_kwargs["od_to_cfu_conversion_factor"] = st.session_state.get("int_od_to_cfu_conversion_factor", 2e8)
        # Optical weight of a dormant/hibernating cell in OD (D, H). BRG/StrainSet forward
        # this to ModelConfig via to_config(**extra_config_kwargs); Direct passes it
        # explicitly to with_od_debris() below. (Previously omitted here → BRG/StrainSet
        # silently ignored int_dormant_od_fraction and used the engine default 1.0.)
        extra_kwargs["dormant_od_fraction"] = st.session_state.get("int_dormant_od_fraction", 1.0)
        # Nutrient-dependent OD absorptivity eps(S) — get_od scales the cell term by it.
        # BRG/StrainSet forward these ModelConfig fields via to_config; Direct instead calls
        # with_nutrient_dependent_od() below (ModelBuilder.build takes no kwargs).
        if st.session_state.get("int_od_nutrient_enabled", False):
            from pbisim import nutrient_dependent_absorptivity
            extra_kwargs["od_absorptivity_function"] = nutrient_dependent_absorptivity
            extra_kwargs["od_abs_exp"] = st.session_state.get("int_od_abs_exp", 1.0)
            extra_kwargs["od_abs_stat"] = st.session_state.get("int_od_abs_stat", 0.4)
            extra_kwargs["od_abs_S50"] = st.session_state.get("int_od_abs_S50", 0.3)
            extra_kwargs["od_abs_hill"] = st.session_state.get("int_od_abs_hill", 2.0)

    # ── Resolve Dose Schedule ─────────────────────────────────────────────────
    dose_events = []
    for d in doses:
        target_type = d["target_type"]
        if target_type == "phage":
            target = "phage"
            target_idx = d["target_idx"]
        elif target_type == "antibiotic":
            target = "antibiotic"
            target_idx = d["target_idx"]
        else:
            target = "nutrient"
            target_idx = 0
            
        event = rec.new(
            DoseEvent,
            time=d["time"],
            amount=d["amount"],
            target=target,
            index=target_idx,
            route=d["route"],
            duration=d.get("duration", 0.0),
        )
        dose_events.append(event)

    if dose_events:
        rec.var("dose_events", dose_events)
        schedule = rec.var("schedule", rec.new(DoseSchedule, dose_events))
    else:
        schedule = None

    # ── BUILDER MODE: Direct (ModelBuilder) ───────────────────────────────────
    if builder_mode == "Direct (ModelBuilder)":
        max_depth = max([s.get("dormancy_depth", 1) for s in strains] if strains else [1])
        builder = rec.init("builder", ModelBuilder, n_bacteria=n_bacteria, n_phages=n_phages, n_latent=n_latent, n_depth=max_depth)

        # Growth rates
        growth_rates = [s["growth_rate"] for s in strains]
        ratios = [s.get("bacteria_to_resource_ratio", 1e9) for s in strains]
        builder = rec.call("builder", builder, "with_growth_rates", growth_rates, bacteria_to_resource_ratio=ratios)

        # Natural death rates
        death_rates_B = [s.get("death_rate_B", 0.0) for s in strains]
        death_rates_D = [s.get("death_rate_D", 0.0) for s in strains]
        _has_active_death = any(db > 0 for db in death_rates_B)
        _dthk = death_kwargs()
        if _has_active_death and "carrying_capacity" in _dthk:  # density death needs K
            builder = rec.call("builder", builder, "with_nutrient", carrying_capacity=_dthk["carrying_capacity"])
        # Always set the death function so the config reflects the chosen signal; the
        # rates are only overridden when > 0 (a None rate means that pathway is off).
        builder = rec.call(
            "builder", builder, "with_death",
            death_rate_B=np.array(death_rates_B) if _has_active_death else None,
            death_rate_D=np.array(death_rates_D) if any(dd > 0 for dd in death_rates_D) else None,
            death_function=_dthk["death_function"],
        )

        # Lysis signal — nutrient-coupled (frac_lysis) vs the default constant lysis.
        _lyk = lysis_kwargs()
        if _lyk["lysis_progression_function"] is not None:
            builder = rec.call(
                "builder", builder, "with_lysis_function",
                _lyk["lysis_progression_function"], Ks_lysis=_lyk["monod_constant_lysis"],
                phi_min=_lyk["lysis_floor"])

        # Dormancy
        any_dormancy = any(s.get("dormancy_enabled", False) for s in strains)
        if any_dormancy:
            dormancy_rates = [s["dormancy_rate"] if s.get("dormancy_enabled", False) else 0.0 for s in strains]
            resus_rates = [s["resuscitation_rate"] if s.get("dormancy_enabled", False) else 0.0 for s in strains]
            diff_rates = [s["dormancy_diffusion_rate"] if s.get("dormancy_enabled", False) else 0.0 for s in strains]
            enabled_strains = [s for s in strains if s.get("dormancy_enabled", False)]
            # Dormancy signal functions + Ks/threshold are model-wide (single engine fields).
            ds = canonical_signal(st.session_state.get("int_dormancy_signal", "nutrient"))
            rs = canonical_signal(st.session_state.get("int_resuscitation_signal", "nutrient"))
            dfs = canonical_signal(st.session_state.get("int_diffusion_signal", "constant"))
            _dorm_ks = float(st.session_state.get("int_dormancy_monod_constant", 0.0) or 0.0)
            _dorm_kdorm = float(st.session_state.get("int_dormancy_carrying_capacity", 0.0) or 0.0)
            # nutrient dormancy/diffusion signals need S tracked; coerce when it isn't.
            ds, _cd = compat_dormancy_signal(ds, track_nutrients)
            rs, _cr = compat_dormancy_signal(rs, track_nutrients)
            dfs, _cdf = compat_dormancy_signal(dfs, track_nutrients)

            _dorm_kwargs = dict(
                dormancy_rate=np.array(dormancy_rates),
                resuscitation_rate=np.array(resus_rates),
                dormancy_diffusion_rate=np.array(diff_rates),
                dormancy_signal=ds,
                resuscitation_signal=rs,
                diffusion_signal=dfs,
            )
            if _dorm_ks > 0 and any(sig in ("nutrient", "nutrient+density") for sig in (ds, rs, dfs)):
                _dorm_kwargs["dormancy_monod_constant"] = _dorm_ks
            # density-based dormancy/resuscitation/diffusion needs a density threshold.
            # Use the per-strain dormancy_carrying_capacity when set, else the growth
            # carrying capacity (Monod growth doesn't set one, so supply it explicitly).
            if any(sig in ("density", "nutrient+density") for sig in (ds, rs, dfs)):
                _dorm_kwargs["dormancy_carrying_capacity"] = (
                    _dorm_kdorm if _dorm_kdorm > 0 else st.session_state.get("int_carrying_capacity", 1e9))
            builder = rec.call("builder", builder, "with_dormancy", **_dorm_kwargs)
        elif not track_nutrients:
            # No dormancy configured, but a frozen S is incompatible with the default
            # nutrient dormancy function — pin nutrient-independent (constant) functions.
            builder = rec.call(
                "builder", builder, "with_dormancy",
                dormancy_rate=0.0, resuscitation_rate=0.0, dormancy_diffusion_rate=0.0,
                dormancy_signal="constant", resuscitation_signal="constant")
            
        # Phages
        if n_phages > 0:
            adsorption_rates = []
            adsorption_rates_dormant = []
            for s_idx in range(n_bacteria):
                s_ads = []
                s_ads_dorm = []
                for p_idx in range(n_phages):
                    s_ads.append(st.session_state.get(f"ads_{s_idx}_{p_idx}", 1e-8 if s_idx == 0 else 0.0))
                    s_ads_dorm.append(st.session_state.get(f"ads_dorm_{s_idx}_{p_idx}", 0.0))
                adsorption_rates.append(s_ads)
                adsorption_rates_dormant.append(s_ads_dorm)
                
            burst_sizes = np.tile(np.array([p["burst_sizes"] for p in phages]), (n_bacteria, 1))
            latent_periods = np.tile(np.array([p["latent_periods"] for p in phages]), (n_bacteria, 1))
            decay_rates = [p["phage_decay_rates"] for p in phages]
            
            # Phage PK
            has_phage_pk = any(p["pk_mode"] != "None" for p in phages)
            if has_phage_pk:
                from pbisim import PhagePKConfig
                vcs = np.array([p.get("Vc", 5000.0) for p in phages])
                k_elims = np.array([p.get("k_elim", 0.2) if p["pk_mode"] != "None" else 0.0 for p in phages])
                k_ins = np.array([p.get("k_in", 0.1) if p["pk_mode"] != "None" else 0.0 for p in phages])
                k_outs = np.array([p.get("k_out", 0.05) if p["pk_mode"] != "None" else 0.0 for p in phages])
                
                # Nonlinear PK Km
                kms = np.array([p.get("Km_elim", 0.0) if p.get("Km_elim", 0.0) > 0 else np.inf for p in phages])
                
                has_mc = any(p["pk_mode"] == "Mass-Conserving" for p in phages)
                vis = np.array([p.get("Vi", 10.0) if p["pk_mode"] == "Mass-Conserving" else 0.0 for p in phages]) if has_mc else None
                _phage_k12, _phage_k21 = _phage_pk_k12_k21(phages)
                
                pk_config = rec.new(
                    PhagePKConfig,
                    n_phages=n_phages, Vc=vcs, k_elim=k_elims, k_in=k_ins, k_out=k_outs, Vi=vis, Km_elim=kms,
                    k12=_phage_k12, k21=_phage_k21,
                )
                builder = rec.call("builder", builder, "with_phage_pk", pk_config)
                
            # Phage decay nonlinear Km
            phage_decay_Km = np.array([p.get("phage_decay_Km", 0.0) if p.get("phage_decay_Km", 0.0) > 0 else np.inf for p in phages])

            # Dormant-adsorption attenuation with dormancy depth (per phage, broadcast
            # across strains): effective dormant rate = adsorption_dormant * exp(-att * depth).
            attenuation_rate = np.tile(
                np.array([p.get("attenuation_rate", 0.0) for p in phages]), (n_bacteria, 1)
            )

            builder = rec.call(
                "builder", builder, "with_phage_params",
                adsorption_rates=np.array(adsorption_rates),
                adsorption_rates_dormant=np.array(adsorption_rates_dormant),
                burst_sizes=np.array(burst_sizes),
                latent_periods=np.array(latent_periods),
                phage_decay_rates=np.array(decay_rates),
                allow_superinfection=superinfection,
                phage_decay_Km=phage_decay_Km,
                attenuation_rate=attenuation_rate,
            )
            
            # Pseudolysogeny
            any_pseudo = any(p.get("hibernation_rate_s", 0.0) > 0 or p.get("hibernation_rate_r", 0.0) > 0 for p in phages)
            if any_pseudo:
                # build n_bacteria x n_phages matrices
                hib_rates = np.zeros((n_bacteria, n_phages))
                res_rates = np.zeros((n_bacteria, n_phages))
                # For Direct mode, broadcast: WT takes susceptible rates, other strains take resistant rates
                for p_idx in range(n_phages):
                    p = phages[p_idx]
                    hib_rates[0, p_idx] = p.get("hibernation_rate_s", 0.0)
                    res_rates[0, p_idx] = p.get("lytic_resumption_rate_s", 0.0)
                    if n_bacteria > 1:
                        hib_rates[1:, p_idx] = p.get("hibernation_rate_r", 0.0)
                        res_rates[1:, p_idx] = p.get("lytic_resumption_rate_r", 0.0)
                builder = rec.call("builder", builder, "with_pseudolysogeny", hibernation_rate=hib_rates, lytic_resumption_rate=res_rates)
                
        # Mutations. A custom mutation-network graph (any n_bacteria) takes
        # precedence; otherwise fall back to the per-phage-locus shortcut, which
        # pbisim only supports when n_bacteria == 2**n_phages.
        _mut_kwargs = {}
        _mut_M = mutation_matrix_from_transitions(st.session_state.get("int_transitions", []), strains)
        if _mut_M is not None:
            _mut_kwargs["mutation_rates"] = _mut_M
        elif n_phages > 0 and n_bacteria == 2**n_phages:
            phg_res_rates = st.session_state.get("direct_phg_res_rates", [1e-7] * n_phages)
            _mut_kwargs["phage_resistance_rates"] = phg_res_rates
        # Phenotypic transitions (bacteria, growth-independent switching) + the phage
        # mutation / transition generators — the same with_mutations call takes all four.
        for _fld, _M in generator_matrices_from_state(
                [s["name"] for s in strains], [p["name"] for p in phages]).items():
            if _M is not None:
                _mut_kwargs[_fld] = _M
        if _mut_kwargs:
            builder = rec.call("builder", builder, "with_mutations", **_mut_kwargs)
            
        # Antibiotics
        for abx in antibiotics:
            builder = rec.call(
                "builder", builder, "with_antibiotic",
                name=abx["name"],
                k_elim=abx["k_elim"],
                Vc=abx.get("Vc", 1.0),
                k12=abx.get("k12", 0.0),
                k21=abx.get("k21", 0.0),
                emax=abx["emax"],
                ec50=abx["ec50"],
                hill=abx.get("hill", 1.0),
                f_lyse=abx.get("f_lyse", 0.0),
                inoculum_effect_constant=abx.get("inoculum_effect_constant", None) if abx.get("inoculum_effect_constant", 0.0) > 0 else None,
                Km_elim=abx.get("Km_elim", None) if abx.get("Km_elim", 0.0) > 0 else None,
            )

        # Nutrients / growth signal
        _gk = growth_nutrient_kwargs()
        _growth_fn = _gk.pop("growth_function")
        _gfn_name = st.session_state.get("int_growth_function")
        if _gfn_name == _GROWTH_SEQUENTIAL:
            # Diauxic growth: with_sequential_growth sets growth_function + the three
            # per-phase arrays; the remaining kwargs (Ks fallback, recycle, s_in/out,
            # density flag) go through with_nutrient. Pop the phase arrays so they
            # aren't also passed to with_nutrient (which doesn't accept them).
            _rf = _gk.pop("growth_phase_rate_factors")
            _mc = _gk.pop("growth_phase_monod")
            _th = _gk.pop("growth_phase_thresholds")
            builder = rec.call("builder", builder, "with_sequential_growth",
                               rate_factors=_rf, monod_constants=_mc, thresholds=_th)
        elif _gfn_name == _GROWTH_SMOOTH_EFF:
            # Smooth two-efficiency Monod: with_smooth_efficiency_growth sets
            # growth_function + monod_K_low/θ/hill; the high-S K (monod_constant) and
            # the rest go through with_nutrient. Pop the smooth-eff fields (with_nutrient
            # doesn't accept them).
            _klow = _gk.pop("monod_K_low")
            _theta = _gk.pop("monod_efficiency_theta")
            _hill = _gk.pop("monod_efficiency_hill")
            builder = rec.call("builder", builder, "with_smooth_efficiency_growth",
                               _klow, theta=_theta, hill=_hill)
        else:
            builder = rec.call("builder", builder, "with_growth_function", _growth_fn)
        builder = rec.call("builder", builder, "with_nutrient", **_gk)

        # Immunity
        immunity_enabled = st.session_state.get("int_immunity_enabled", False)
        if immunity_enabled:
            kill_rate_D = st.session_state.get("int_imm_kill_rate_D", 0.0)
            builder = rec.call(
                "builder", builder, "with_immunity",
                imm_stim_rate=np.full(n_bacteria, st.session_state.get("int_imm_stim_rate", 0.1)),
                imm_stim50=st.session_state.get("int_imm_stim50", 1e6),
                imm_kill_rate=np.full(n_bacteria, st.session_state.get("int_innate_kill_rate", 1e7)),
                imm_kill50=st.session_state.get("int_innate_kill50", 1e5),
                imm_decay_rate=st.session_state.get("int_innate_decay_rate", 0.1),
                immune_module=st.session_state.get("int_immune_module", "innate"),
                imm_max=st.session_state.get("int_innate_max", 1e7),
                imm_kill_rate_D=np.array([kill_rate_D] * n_bacteria) if kill_rate_D > 0 else None
            )

        if schedule:
            builder = rec.call("builder", builder, "with_dose_schedule", schedule)

        # OD / debris ODE (Direct mode). ModelBuilder.build() takes no kwargs — debris
        # must be configured via with_od_debris(), not passed to build().
        if debris_enabled:
            builder = rec.call(
                "builder", builder, "with_od_debris",
                u=extra_kwargs.get("debris_u", 0.4),
                v=extra_kwargs.get("debris_v", 0.2),
                kdis=extra_kwargs.get("debris_kdis", 0.01),
                od_to_cfu_conversion_factor=extra_kwargs.get("od_to_cfu_conversion_factor", 2e8),
                dormant_od_fraction=st.session_state.get("int_dormant_od_fraction", 1.0),
            )
            # Nutrient-dependent OD absorptivity eps(S) — pairs with with_od_debris.
            if st.session_state.get("int_od_nutrient_enabled", False):
                builder = rec.call(
                    "builder", builder, "with_nutrient_dependent_od",
                    exponential=st.session_state.get("int_od_abs_exp", 1.0),
                    stationary=st.session_state.get("int_od_abs_stat", 0.4),
                    S50=st.session_state.get("int_od_abs_S50", 0.3),
                    hill=st.session_state.get("int_od_abs_hill", 2.0))

        config = rec.result("cfg", "builder", builder, "build")

        initial_B = np.array([s["initial_B"] for s in strains])
        initial_P = np.array([p["initial_P"] for p in phages])
        initial_S = st.session_state.get("int_initial_S", 1.0) if track_nutrients else 1.0

        model_kwargs = {}
        if immunity_enabled:
            model_kwargs["initial_Imm"] = st.session_state.get("int_imm_initial", 0.0)
            
        # Dormant initial conditions — use per-strain initial_D (default 0).
        # PBIModel accepts shape (n_bacteria,) and distributes evenly across Q layers.
        if any_dormancy:
            ic_D = np.array([s.get("initial_D", 0.0) for s in strains])
            if np.any(ic_D > 0):
                model_kwargs["initial_D"] = ic_D
            
        return config, initial_B, initial_P, initial_S, model_kwargs

    # ── BUILDER MODE: Binary Resistance Genotypes (BRG) ──────────────────────
    elif builder_mode == "Binary Genotypes (BRG)":
        # Build BacterialStrain
        base_growth = st.session_state.get("int_brg_base_growth", 1.2)
        base_ratio = st.session_state.get("int_brg_base_ratio", 1e9)
        dormancy_enabled = st.session_state.get("int_brg_dormancy_enabled", False)
        
        dorm_rate = st.session_state.get("int_brg_dorm_rate", 0.001) if dormancy_enabled else 0.0
        resus_rate = st.session_state.get("int_brg_resus_rate", 0.1) if dormancy_enabled else 0.0
        diff_rate = st.session_state.get("int_brg_diff_rate", 0.05) if dormancy_enabled else 0.0
        
        b = rec.init(
            "bacteria", BacterialStrain,
            base_growth_rate=base_growth,
            bacteria_to_resource_ratio=base_ratio,
            dormancy_rate=dorm_rate,
            resuscitation_rate=resus_rate,
            dormancy_diffusion_rate=diff_rate,
            death_rate_B=st.session_state.get("int_brg_death_rate_B", 0.0) if st.session_state.get("int_brg_death_rate_B", 0.0) > 0 else None,
            death_rate_D=st.session_state.get("int_brg_death_rate_D", 0.0) if dormancy_enabled and st.session_state.get("int_brg_death_rate_D", 0.0) > 0 else None,
        )

        # Build PhageStrains
        phage_strains = []
        for p in phages:
            phage_strains.append(
                rec.new(
                    PhageStrain,
                    name=p["name"],
                    adsorption_s=p.get("adsorption_s", 5e-8),
                    adsorption_r=p.get("adsorption_r", 0.0),
                    adsorption_dormant_s=p.get("adsorption_dormant_s", 0.0),
                    adsorption_dormant_r=p.get("adsorption_dormant_r", 0.0),
                    burst_size_s=p["burst_sizes"],
                    latent_period_s=p["latent_periods"],
                    decay_rate=p["phage_decay_rates"],
                    fitness_cost=p.get("fitness_cost", 0.0),
                    mu=p.get("mu", 1e-7),
                    hibernation_rate_s=p.get("hibernation_rate_s", None) if p.get("hibernation_rate_s", 0.0) > 0 else None,
                    hibernation_rate_r=p.get("hibernation_rate_r", None) if p.get("hibernation_rate_r", 0.0) > 0 else None,
                    lytic_resumption_rate_s=p.get("lytic_resumption_rate_s", None) if p.get("lytic_resumption_rate_s", 0.0) > 0 else None,
                    lytic_resumption_rate_r=p.get("lytic_resumption_rate_r", None) if p.get("lytic_resumption_rate_r", 0.0) > 0 else None,
                )
            )
        rec.var("phages", phage_strains)

        # Build Antibiotics
        abx_strains = []
        for abx in antibiotics:
            abx_strains.append(
                rec.new(
                    Antibiotic,
                    name=abx["name"],
                    emax_s=abx["emax"],
                    emax_r=abx.get("emax_r", abx["emax"] * 0.1),
                    ec50_s=abx["ec50"],
                    ec50_r=abx.get("ec50_r", abx["ec50"] * 10.0),
                    k_elim=abx["k_elim"],
                    hill=abx.get("hill", 1.0),
                    fitness_cost=abx.get("fitness_cost", 0.0),
                    mu=abx.get("mu", 1e-7),
                    Vc=abx.get("Vc", 1.0),
                    k12=abx.get("k12", 0.0),
                    k21=abx.get("k21", 0.0),
                    inoculum_effect_constant=abx.get("inoculum_effect_constant", None) if abx.get("inoculum_effect_constant", 0.0) > 0 else None,
                )
            )
        if abx_strains:
            rec.var("antibiotics", abx_strains)

        brg = rec.classcall(
            "brg", BinaryResistanceGenotypes, "from_strains",
            phage_strains,
            bacteria=b,
            antibiotics=abx_strains if abx_strains else None,
        )
        
        # Build config — dormancy depth compartments (configurable; 1 when dormancy off)
        max_depth = int(st.session_state.get("int_brg_n_depth", 1)) if dormancy_enabled else 1

        # Resolve Phage PK config
        phage_pk_config = None
        has_phage_pk = any(p["pk_mode"] != "None" for p in phages)
        if has_phage_pk:
            from pbisim import PhagePKConfig
            vcs = np.array([p.get("Vc", 5000.0) for p in phages])
            k_elims = np.array([p.get("k_elim", 0.2) if p["pk_mode"] != "None" else 0.0 for p in phages])
            k_ins = np.array([p.get("k_in", 0.1) if p["pk_mode"] != "None" else 0.0 for p in phages])
            k_outs = np.array([p.get("k_out", 0.05) if p["pk_mode"] != "None" else 0.0 for p in phages])
            kms = np.array([p.get("Km_elim", 0.0) if p.get("Km_elim", 0.0) > 0 else np.inf for p in phages])
            has_mc = any(p["pk_mode"] == "Mass-Conserving" for p in phages)
            vis = np.array([p.get("Vi", 10.0) if p["pk_mode"] == "Mass-Conserving" else 0.0 for p in phages]) if has_mc else None
            _phage_k12, _phage_k21 = _phage_pk_k12_k21(phages)
            
            phage_pk_config = rec.new(
                PhagePKConfig,
                n_phages=n_phages, Vc=vcs, k_elim=k_elims, k_in=k_ins, k_out=k_outs, Vi=vis, Km_elim=kms,
                k12=_phage_k12, k21=_phage_k21,
            )

        # Expose nonlinear clearances
        extra_kwargs["allow_superinfection"] = superinfection

        # Growth signal + nutrient config (forwarded to ModelConfig via to_config).
        extra_kwargs.update(growth_nutrient_kwargs())
        # Dormancy signal functions (+ Ks / K_dorm) — model-wide (topmost panel).
        extra_kwargs.update(mode_dormancy_kwargs(
            st.session_state.get("int_dormancy_signal", "nutrient"),
            st.session_state.get("int_resuscitation_signal", "nutrient"),
            float(st.session_state.get("int_dormancy_monod_constant", 0.0) or 0.0),
            float(st.session_state.get("int_dormancy_carrying_capacity", 0.0) or 0.0),
        ))
        extra_kwargs.update(death_kwargs())  # death signal function (+ K for density death)
        # Lysis signal (nutrient-coupled frac_lysis vs the default constant lysis).
        _lyk = lysis_kwargs()
        extra_kwargs["lysis_progression_function"] = _lyk["lysis_progression_function"]
        if _lyk["monod_constant_lysis"] is not None:
            extra_kwargs["monod_constant_lysis"] = _lyk["monod_constant_lysis"]
        extra_kwargs["lysis_floor"] = _lyk["lysis_floor"]  # phi_min (no-op for constant lysis)

        # Dose schedule
        if schedule:
            extra_kwargs["dose_schedule"] = schedule
            
        # Immunity
        immunity_enabled = st.session_state.get("int_immunity_enabled", False)
        if immunity_enabled:
            extra_kwargs["imm_stim_rate"] = np.full(brg.n_strains, st.session_state.get("int_imm_stim_rate", 0.1))
            extra_kwargs["imm_stim50"] = st.session_state.get("int_imm_stim50", 1e6)
            extra_kwargs["imm_kill_rate"] = np.full(brg.n_strains, st.session_state.get("int_innate_kill_rate", 1e7))
            extra_kwargs["imm_kill50"] = st.session_state.get("int_innate_kill50", 1e5)
            extra_kwargs["imm_decay_rate"] = st.session_state.get("int_innate_decay_rate", 0.1)
            extra_kwargs["immune_module"] = st.session_state.get("int_immune_module", "innate")
            extra_kwargs["imm_max"] = st.session_state.get("int_innate_max", 1e7)
            kill_rate_D = st.session_state.get("int_imm_kill_rate_D", 0.0)
            if kill_rate_D > 0:
                extra_kwargs["imm_kill_rate_D"] = np.array([kill_rate_D] * brg.n_strains)

        # Per-phage dormant-adsorption attenuation (broadcast across genotypes).
        if phages:
            extra_kwargs["attenuation_rate"] = np.array(
                [p.get("attenuation_rate", 0.0) for p in phages]
            )

        # Per-phage nonlinear (Michaelis-Menten) phage decay saturation.
        if phages and any(p.get("phage_decay_Km", 0.0) > 0 for p in phages):
            extra_kwargs["phage_decay_Km"] = np.array(
                [p.get("phage_decay_Km", 0.0) if p.get("phage_decay_Km", 0.0) > 0 else np.inf for p in phages]
            )

        config = rec.result(
            "cfg", "brg", brg, "to_config",
            n_latent=n_latent,
            n_depth=max_depth,
            phage_pk_config=phage_pk_config,
            **extra_kwargs
        )
        # BRG.to_config doesn't expose the depth-diffusion signal — set it on the config
        # directly from the BRG diffusion selector (post-build).
        if dormancy_enabled:
            config = set_diffusion_functions(
                config, st.session_state.get("int_diffusion_signal", "constant"), rec=rec)
        # Likewise BRG.to_config hard-codes zero phenotypic-transition / phage generator
        # matrices (its mutation matrix comes from the per-locus μ) — set them post-build.
        # Bacterial nodes are the genotype labels, phage nodes the locus names.
        config = apply_generator_matrices(
            config, brg.strain_labels, [p["name"] for p in phages], rec=rec)

        # Resolve initial densities
        if st.session_state.get("int_brg_use_eq_ic", False):
            total_B = st.session_state.get("int_brg_eq_total_B", 1e7)
            # Record as the derivation call so the repro shows the pbisim function, not
            # just the resulting numbers.
            initial_B = rec.expr(
                brg.equilibrium_initial_condition(total_bacteria=total_B),
                f"brg.equilibrium_initial_condition(total_bacteria={rec.render(total_B)})",
            )
        else:
            initial_B = np.zeros(brg.n_strains)
            saved_init_B = st.session_state.get("int_brg_initial_B", {})
            for idx, lbl in enumerate(brg.strain_labels):
                initial_B[idx] = saved_init_B.get(lbl, 1e7 if idx == 0 else 0.0)
            
        initial_P = np.array([p["initial_P"] for p in phages])
        initial_S = st.session_state.get("int_initial_S", 1.0) if track_nutrients else 1.0
        
        model_kwargs = {}
        if immunity_enabled:
            model_kwargs["initial_Imm"] = st.session_state.get("int_imm_initial", 0.0)

        return config, initial_B, initial_P, initial_S, model_kwargs

    # ── BUILDER MODE: Custom Strains & Graph (StrainSet) ──────────────────────
    else:
        ss = rec.init("ss", StrainSet, n_phages=n_phages)

        # Register antibiotics
        for abx in antibiotics:
            rec.mutate(
                "ss", ss, "add_antibiotic",
                rec.new(
                    AntibioticDefinition,
                    name=abx["name"],
                    k_elim=abx["k_elim"],
                    Vc=abx.get("Vc", 1.0),
                    k12=abx.get("k12", 0.0),
                    k21=abx.get("k21", 0.0),
                    inoculum_effect_constant=abx.get("inoculum_effect_constant", None) if abx.get("inoculum_effect_constant", 0.0) > 0 else None,
                    Km_elim=abx.get("Km_elim", None) if abx.get("Km_elim", 0.0) > 0 else None,
                )
            )
            
        # Add Strains
        for i, s in enumerate(strains):
            # Parse phage parameters
            ads_rates = []
            ads_rates_dorm = []
            bursts = []
            latents = []
            latents_dorm = []
            hibernations = []
            resumptions = []
            for p_idx in range(n_phages):
                ads_rates.append(st.session_state.get(f"ads_{i}_{p_idx}", 1e-8 if i == 0 else 0.0))
                ads_rates_dorm.append(st.session_state.get(f"ads_dorm_{i}_{p_idx}", 0.0))
                p = phages[p_idx]
                bursts.append(p["burst_sizes"])
                latents.append(p["latent_periods"])
                latents_dorm.append(p["latent_periods"]) # nominal copy
                hibernations.append(p.get("hibernation_rate_s", 0.0) if i == 0 else p.get("hibernation_rate_r", 0.0))
                resumptions.append(p.get("lytic_resumption_rate_s", 0.0) if i == 0 else p.get("lytic_resumption_rate_r", 0.0))
                
            # Parse antibiotic sensitivities
            sensitivities = {}
            for abx in antibiotics:
                # default: first strain is susceptible, others are resistant
                emax_val = abx["emax"] if i == 0 else abx["emax"] * 0.1
                ec50_val = abx["ec50"] if i == 0 else abx["ec50"] * 10.0
                sensitivities[abx["name"]] = rec.new(AntibioticSensitivity, emax=emax_val, ec50=ec50_val)

            rec.mutate(
                "ss", ss, "add_strain",
                rec.new(
                    StrainDefinition,
                    name=s["name"],
                    growth_rate=s["growth_rate"],
                    adsorption_rates=np.array(ads_rates),
                    adsorption_rates_dormant=np.array(ads_rates_dorm),
                    burst_sizes=np.array(bursts),
                    latent_periods=np.array(latents),
                    latent_periods_dormant=np.array(latents_dorm),
                    bacteria_to_resource_ratio=s.get("bacteria_to_resource_ratio", 1e9),
                    dormancy_rate=s["dormancy_rate"] if s.get("dormancy_enabled", False) else 0.0,
                    resuscitation_rate=s["resuscitation_rate"] if s.get("dormancy_enabled", False) else 0.0,
                    dormancy_diffusion_rate=s["dormancy_diffusion_rate"] if s.get("dormancy_enabled", False) else 0.0,
                    imm_stim_rate=st.session_state.get("int_imm_stim_rate", 0.1) if st.session_state.get("int_immunity_enabled", False) else 0.0,
                    imm_kill_rate=st.session_state.get("int_innate_kill_rate", 1e7) if st.session_state.get("int_immunity_enabled", False) else 0.0,
                    attenuation_rate=np.array([p.get("attenuation_rate", 0.0) for p in phages]),
                    death_rate_B=s.get("death_rate_B", None) if s.get("death_rate_B", 0.0) > 0 else None,
                    death_rate_D=s.get("death_rate_D", None) if s.get("death_rate_D", 0.0) > 0 else None,
                    hibernation_rate=np.array(hibernations) if any(h > 0 for h in hibernations) else None,
                    lytic_resumption_rate=np.array(resumptions) if any(r > 0 for r in resumptions) else None,
                    antibiotic_sensitivity=sensitivities
                )
            )

        # Build mutation graph (+ the phenotypic-transition graph and the phage
        # mutation / transition generators via StrainSet's dedicated setters).
        _ss_names = [s["name"] for s in strains]
        graph_dict = graph_dict_from_edges(st.session_state.get("int_transitions", []), _ss_names)
        if graph_dict:
            rec.mutate("ss", ss, "set_mutation_graph", graph_dict)
        _tr_graph = graph_dict_from_edges(st.session_state.get("int_bact_transitions", []), _ss_names)
        if _tr_graph:
            rec.mutate("ss", ss, "set_transition_graph", _tr_graph)
        _gen = generator_matrices_from_state(_ss_names, [p["name"] for p in phages])
        if _gen["mutation_rates_phage"] is not None:
            rec.mutate("ss", ss, "set_phage_mutation_rates", _gen["mutation_rates_phage"])
        if _gen["transition_rates_phage"] is not None:
            rec.mutate("ss", ss, "set_phage_transition_rates", _gen["transition_rates_phage"])
            
        # Phage PK
        phage_pk_config = None
        has_phage_pk = any(p["pk_mode"] != "None" for p in phages)
        if has_phage_pk:
            from pbisim import PhagePKConfig
            vcs = np.array([p.get("Vc", 5000.0) for p in phages])
            k_elims = np.array([p.get("k_elim", 0.2) if p["pk_mode"] != "None" else 0.0 for p in phages])
            k_ins = np.array([p.get("k_in", 0.1) if p["pk_mode"] != "None" else 0.0 for p in phages])
            k_outs = np.array([p.get("k_out", 0.05) if p["pk_mode"] != "None" else 0.0 for p in phages])
            kms = np.array([p.get("Km_elim", 0.0) if p.get("Km_elim", 0.0) > 0 else np.inf for p in phages])
            has_mc = any(p["pk_mode"] == "Mass-Conserving" for p in phages)
            vis = np.array([p.get("Vi", 10.0) if p["pk_mode"] == "Mass-Conserving" else 0.0 for p in phages]) if has_mc else None
            _phage_k12, _phage_k21 = _phage_pk_k12_k21(phages)
            
            phage_pk_config = rec.new(
                PhagePKConfig,
                n_phages=n_phages, Vc=vcs, k_elim=k_elims, k_in=k_ins, k_out=k_outs, Vi=vis, Km_elim=kms,
                k12=_phage_k12, k21=_phage_k21,
            )

        decay_rates = np.array([p["phage_decay_rates"] for p in phages])
        max_depth = max([s.get("dormancy_depth", 1) for s in strains] if strains else [1])
        
        # Expose extra kwargs
        extra_kwargs["allow_superinfection"] = superinfection

        # Growth signal + nutrient config (forwarded to ModelConfig via to_config).
        extra_kwargs.update(growth_nutrient_kwargs())
        # Dormancy signal functions (+ Ks / K_dorm) — model-wide (topmost panel).
        extra_kwargs.update(mode_dormancy_kwargs(
            st.session_state.get("int_dormancy_signal", "nutrient"),
            st.session_state.get("int_resuscitation_signal", "nutrient"),
            float(st.session_state.get("int_dormancy_monod_constant", 0.0) or 0.0),
            float(st.session_state.get("int_dormancy_carrying_capacity", 0.0) or 0.0),
        ))
        extra_kwargs.update(death_kwargs())  # death signal function (+ K for density death)
        # Lysis signal (nutrient-coupled frac_lysis vs the default constant lysis).
        _lyk = lysis_kwargs()
        extra_kwargs["lysis_progression_function"] = _lyk["lysis_progression_function"]
        if _lyk["monod_constant_lysis"] is not None:
            extra_kwargs["monod_constant_lysis"] = _lyk["monod_constant_lysis"]
        extra_kwargs["lysis_floor"] = _lyk["lysis_floor"]  # phi_min (no-op for constant lysis)

        # Dose schedule
        if schedule:
            extra_kwargs["dose_schedule"] = schedule

        # Per-phage nonlinear (Michaelis-Menten) phage decay saturation.
        if phages and any(p.get("phage_decay_Km", 0.0) > 0 for p in phages):
            extra_kwargs["phage_decay_Km"] = np.array(
                [p.get("phage_decay_Km", 0.0) if p.get("phage_decay_Km", 0.0) > 0 else np.inf for p in phages]
            )

        # Immunity defaults
        immunity_enabled = st.session_state.get("int_immunity_enabled", False)

        config = rec.result(
            "cfg", "ss", ss, "to_config",
            n_latent=n_latent,
            n_depth=max_depth,
            phage_decay_rates=decay_rates,
            imm_decay_rate=st.session_state.get("int_innate_decay_rate", 0.1),
            imm_stim50=st.session_state.get("int_imm_stim50", 1e6),
            imm_kill50=st.session_state.get("int_innate_kill50", 1e5),
            phage_pk_config=phage_pk_config,
            immune_module=st.session_state.get("int_immune_module", "innate"),
            imm_max=st.session_state.get("int_innate_max", 1e7),
            **extra_kwargs  # includes growth_function + monod_constant/recycle/etc.
        )
        # StrainSet.to_config doesn't expose the depth-diffusion signal — set it on the
        # config directly from the model-wide diffusion signal (post-build), when any
        # strain has dormancy enabled.
        if any(s.get("dormancy_enabled", False) for s in strains):
            config = set_diffusion_functions(
                config, st.session_state.get("int_diffusion_signal", "constant"), rec=rec)

        initial_B = np.array([s["initial_B"] for s in strains])
        initial_P = np.array([p["initial_P"] for p in phages])
        initial_S = st.session_state.get("int_initial_S", 1.0) if track_nutrients else 1.0

        model_kwargs = {}
        if immunity_enabled:
            model_kwargs["initial_Imm"] = st.session_state.get("int_imm_initial", 0.0)
            
        return config, initial_B, initial_P, initial_S, model_kwargs


def run_sim_from_gui_params():
    """Builds and solves the ODE for the nominal patient."""
    config, initial_B, initial_P, initial_S, model_kwargs = build_nominal_config_from_gui()
    
    # ── Pretreatment equilibrate prerun ───────────────────────────────────────
    t_prerun = st.session_state.get("int_t_prerun", 0.0)
    if t_prerun > 0:
        # Equilibrate to stationary phase (no treatment) and carry the FULL final
        # state into the treatment sim: active B, the dormant reservoir D, nutrient
        # S, and immune priming Imm. Previously only B and S were kept — discarding
        # ic.D silently dropped the (usually dominant) dormant population, so longer
        # preruns lost more of the culture until the residual active cells fell below
        # the extinction floor and the treatment plotted as a flat 0 CFU curve.
        # initial_S sets the pre-run's nutrient level (Monod growth caps the stationary
        # density by S0). Before the engine gained this arg, the pre-run always grew
        # from S=1.0 regardless of the configured S0.
        ic = stationary_phase_ic(config, t_prerun=t_prerun, B0=initial_B, initial_S=initial_S)
        initial_B = ic.B
        initial_S = max(float(ic.S), 0.0)  # prerun can leave S slightly negative (numerical)
        if ic.D is not None:
            model_kwargs["initial_D"] = ic.D
        if ic.Imm is not None:
            model_kwargs["initial_Imm"] = ic.Imm
        _carry_prerun_debris(ic, model_kwargs)

    model = PBIModel(
        config,
        initial_B=initial_B,
        initial_P=initial_P,
        initial_S=initial_S,
        **model_kwargs
    )
    
    t_end = st.session_state.get("int_t_end", 48.0)
    dt = st.session_state.get("int_dt", 0.25)
    method = st.session_state.get("int_solver_method", "BDF")
    extinction_threshold = st.session_state.get("int_extinction_threshold", 1.0) or None
    extinction_check_interval = st.session_state.get("int_extinction_check_interval", 0.0) or None
    result = solve_ode(model, t_end=t_end, dt=dt, method=method,
                       extinction_threshold=extinction_threshold,
                       extinction_check_interval=extinction_check_interval)
    return result, config


def prerun_collapse_fraction(config, B0, t_prerun, initial_S=None):
    """Fraction of the inoculum surviving a stationary-phase pre-run (B + D total).

    Returns ``None`` when there is no pre-run or the pre-run can't be evaluated.
    A small value flags the trap where a natural death rate with no dormancy
    decimates the culture during the pre-run (death keeps acting once nutrients
    exhaust and growth stops), so the treatment — and its CFU/OD curve — starts
    far below the inoculum. ``initial_S`` sets the pre-run nutrient level so the
    estimate matches the configured S0 (a low S0 caps stationary density).
    """
    if not t_prerun or t_prerun <= 0:
        return None
    try:
        ic = stationary_phase_ic(config, t_prerun=t_prerun, B0=B0, initial_S=initial_S)
    except Exception:
        return None
    total = float(np.sum(ic.B)) + (float(np.sum(ic.D)) if ic.D is not None else 0.0)
    b0 = float(np.sum(B0))
    return (total / b0) if b0 > 0 else None


def warn_if_prerun_collapses(config, B0, initial_S=None):
    """Emit a Streamlit warning if the configured pre-run decimates the culture."""
    t_prerun = st.session_state.get("int_t_prerun", 0.0)
    frac = prerun_collapse_fraction(config, B0, t_prerun, initial_S=initial_S)
    if frac is not None and frac < 0.1:
        st.warning(
            f"The {t_prerun:g} h pre-run leaves only ~{frac*100:.2g}% of the inoculum "
            "(active + dormant) at treatment start, so the CFU/OD curves scale very low. "
            "A natural death rate with no dormancy makes the culture decline during the "
            "pre-run once nutrients exhaust (death keeps acting after growth stops). "
            "Reduce the pre-run duration or the death rate, or enable dormancy so persisters "
            "survive the pre-run."
        )


def reseed_widget_config(store_key, prefixes):
    """Re-seed widget selections from a persistent store BEFORE the widgets render.

    Streamlit drops a widget's key from session_state whenever the widget is not
    rendered on a rerun, so leaving a page (e.g. a sweep page) and coming back would
    otherwise reset all its controls. Shadowing the keys into a plain dict (which
    survives) and re-seeding from it keeps the configuration alive across navigation.
    """
    store = st.session_state.setdefault(store_key, {})
    for k, v in list(store.items()):
        if any(k.startswith(p) for p in prefixes) and k not in st.session_state:
            try:
                st.session_state[k] = v
            except Exception:
                pass  # non-settable widget (e.g. a button) — skip


def save_widget_config(store_key, prefixes, exclude=()):
    """Shadow the current widget selections (keys matching *prefixes*) into the store."""
    store = st.session_state.setdefault(store_key, {})
    for k in list(st.session_state.keys()):
        if k in exclude:
            continue
        if any(k.startswith(p) for p in prefixes):
            store[k] = st.session_state[k]


def generate_reproduction_code() -> str:
    """Generate a standalone script that reproduces the current simulation.

    The configuration-building portion is emitted from the *same* builder calls the app
    actually makes (recorded into a ``_ReproRecorder`` by ``build_nominal_config_from_gui``),
    so the script mirrors the chosen builder — ``ModelBuilder`` / ``BinaryResistanceGenotypes``
    / ``StrainSet`` — call-for-call and can never drift from what is simulated. Only the
    initial-conditions / solve / plot tail is templated here (the Plotly charts in the app
    are rendered with matplotlib in the script, as noted).
    """
    config, iB, iP, iS, mkw = build_nominal_config_from_gui()
    rec = st.session_state["_repro_rec"]

    t_prerun = st.session_state.get("int_t_prerun", 0.0)

    imports = set(rec.imports) | {"PBIModel", "solve_ode"}
    if t_prerun > 0:
        imports.add("stationary_phase_ic")

    code = [
        "# \u2500\u2500 pbisim Auto-Generated Reproduction Script \u2500\u2500",
        "# Config built via the exact builder calls the app made (no re-derivation).",
        "import numpy as np",
        "import matplotlib.pyplot as plt",
        "from pbisim import " + ", ".join(sorted(imports)),
        *sorted(rec.raw_imports),
        "",
        "# 1. Build the model configuration",
        *rec.lines,
        "",
        "# 2. Initial conditions",
        f"initial_B = {rec.render(iB)}",
        f"initial_P = {rec.render(iP)}",
        f"initial_S = {rec.render(iS)}",
    ]

    if t_prerun > 0:
        _inherit_debris = (st.session_state.get("int_prerun_inherit_debris", True)
                           and st.session_state.get("int_debris_enabled", False))
        _debris_arg = ", initial_Debris=ic.Debris" if _inherit_debris else ""
        code.append(f"# Stationary-phase pre-run ({t_prerun} h) starting from the configured inoculum;")
        code.append("# carry the full final state — active B, dormant D, nutrient S, immune Imm"
                     + (", debris." if _inherit_debris else " (dead-cell debris washed out)."))
        code.append(f"ic = stationary_phase_ic(cfg, t_prerun={t_prerun}, B0=initial_B, initial_S=initial_S)")
        code.append(
            "model = PBIModel(cfg, initial_B=ic.B, initial_P=initial_P, "
            f"initial_S=max(float(ic.S), 0.0), initial_D=ic.D, initial_Imm=(ic.Imm or 0.0){_debris_arg})"
        )
    else:
        _extra_ic = "".join(f", {k}={rec.render(v)}" for k, v in mkw.items())
        code.append(
            f"model = PBIModel(cfg, initial_B=initial_B, initial_P=initial_P, initial_S=initial_S{_extra_ic})"
        )

    _method = st.session_state.get("int_solver_method", "BDF")
    _thresh = st.session_state.get("int_extinction_threshold", 1.0) or None
    _check_int = st.session_state.get("int_extinction_check_interval", 0.0) or None
    code += [
        "",
        "# 3. Solve",
        f"result = solve_ode(model, t_end={st.session_state.get('int_t_end', 48.0)}, "
        f"dt={st.session_state.get('int_dt', 0.25)}, method='{_method}', "
        f"extinction_threshold={_thresh}, extinction_check_interval={_check_int})",
        "",
        "# 4. Plot trajectories",
        "fig, ax = plt.subplots(figsize=(8, 4))",
        "# CFU = culturable cells only (active B + dormant D); infected I / hibernating H don't plate.",
        "ax.semilogy(result.time, np.maximum(result.sum_prefixes('B', 'D'), 1.0), label='CFU (B+D)')",
        "ax.set(xlabel='Time (h)', ylabel='Density (cells/mL)', title='Simulation Run')",
        "ax.legend()",
        "plt.show()",
    ]

    _script = "\n".join(code)
    st.session_state["_last_repro_code"] = _script  # for tests / debugging
    return _script


def _repro_base_config_block(rec, iB, iP, iS, mk, extra_imports=(), initial_P_override=None):
    """Shared prefix for sweep reproduction scripts: imports + the recorded base-config
    builder calls + initial_B/P/S/model_kwargs statements. ``initial_P_override`` lets the
    caller substitute a modified initial_P array (e.g. swept phages zeroed)."""
    imports = set(rec.imports) | {"PBIModel", "solve_ode"} | set(extra_imports)
    lines = [
        "import numpy as np",
        "import matplotlib.pyplot as plt",
        "from pbisim import " + ", ".join(sorted(imports)),
        "from pbisim_app.sweep_helper import apply_sweep_parameter",
        "",
        "# 1. Nominal model configuration (exact shadow of the app's builder calls)",
        *rec.lines,
        f"initial_B = {rec.render(iB)}",
        f"initial_P = {rec.render(iP if initial_P_override is None else initial_P_override)}",
        f"initial_S = {rec.render(iS)}",
        f"model_kwargs = {rec.render(mk)}",
    ]
    return imports, lines


def _repro_solve_kwargs():
    _thresh = st.session_state.get("int_extinction_threshold", 1.0) or None
    _check = st.session_state.get("int_extinction_check_interval", 0.0) or None
    return (f"t_end={st.session_state.get('int_t_end', 48.0)}, dt={st.session_state.get('int_dt', 0.25)}, "
            f"method='{st.session_state.get('int_solver_method', 'BDF')}', "
            f"extinction_threshold={_thresh}, extinction_check_interval={_check}")


def _repro_prerun_lines(indent, t_prerun):
    """Emit the per-run pre-run block (carry B/D/S/Imm) at the given indent.

    Pops a sweep's ``_t_prerun_override`` (set by apply_sweep_parameter when the
    pre-run duration is the swept parameter) so the swept value is used AND the key
    never reaches PBIModel; falls back to the fixed base duration otherwise. Returns
    [] only when there is neither a base pre-run nor any chance of a swept one."""
    p = indent
    base = float(t_prerun or 0.0)
    lines = [
        f"{p}_tp = mk_k.pop('_t_prerun_override', None)",
        f"{p}_tp = {base} if _tp is None else _tp",
        f"{p}if _tp and _tp > 0:",
        f"{p}    ic = stationary_phase_ic(c_k, t_prerun=_tp, B0=ib_k, initial_S=is_k)",
        f"{p}    ib_k, is_k = ic.B, max(float(ic.S), 0.0)",
        f"{p}    if ic.D is not None: mk_k['initial_D'] = ic.D",
        f"{p}    if ic.Imm is not None: mk_k['initial_Imm'] = ic.Imm",
    ]
    if (st.session_state.get("int_prerun_inherit_debris", True)
            and st.session_state.get("int_debris_enabled", False)):
        lines.append(f"{p}    if ic.Debris is not None: mk_k['initial_Debris'] = ic.Debris")
    return lines


def generate_param_sweep_reproduction_code() -> str:
    """Reproduction script for a 1D parameter sweep: the recorded base config plus a loop
    calling the app's own ``apply_sweep_parameter`` per value (no re-derivation)."""
    if st.session_state.get("ps_sweep_type", "1D Sweep") != "1D Sweep":
        return ("# Reproduction code is currently available for 1D parameter sweeps only.\n"
                "# Switch 'Sweep Dimension' to '1D Sweep' to export a script.")
    config, iB, iP, iS, mk = build_nominal_config_from_gui()
    rec = st.session_state["_repro_rec"]

    sweep_params = get_sweep_parameters(
        config, st.session_state.get("int_strains", []),
        st.session_state.get("int_phages", []), st.session_state.get("int_antibiotics", []))
    label = st.session_state.get("p1_sweep_label")
    if not label or label not in sweep_params:
        return "# Configure and run a 1D parameter sweep first."
    meta = sweep_params[label]

    if meta.get("type") == "categorical":
        # A categorical (signal-function) sweep changes which growth/death/lysis/dormancy
        # function the model uses. That choice is resolved by the full model builder at
        # BUILD time (it also toggles track_nutrients and coerces incompatible dormancy
        # signals), so — unlike a numeric sweep — it can't be reproduced by mutating a
        # single already-built config in a standalone loop. The app therefore rebuilds
        # the model per option internally. Reproduction code is not exported for it.
        _opts = list(meta.get("options", {}))
        return (
            "# ── Categorical (signal-function) sweep ──\n"
            f"# Swept: {label}\n"
            f"# Options compared: {_opts}\n"
            "#\n"
            "# No standalone reproduction script is exported for categorical sweeps. The\n"
            "# signal function is resolved by the full model builder at build time (it also\n"
            "# sets track_nutrients and coerces incompatible dormancy signals), so it can't\n"
            "# be reproduced by mutating one built config in a loop — the app rebuilds the\n"
            "# model per option. To script this yourself, rebuild the model for each option\n"
            "# (e.g. call ModelBuilder.with_growth_function(...) / .with_dormancy(...) etc.\n"
            "# per option) rather than editing a single config's function fields.\n"
            "#\n"
            "# The 'View Python Reproduction Code' export covers 1D NUMERIC sweeps.")

    t_prerun = st.session_state.get("int_t_prerun", 0.0)
    imports, code = _repro_base_config_block(
        rec, iB, iP, iS, mk, extra_imports={"stationary_phase_ic"} if t_prerun else set())
    code[0:0] = ["# ── pbisim 1D Parameter-Sweep Reproduction Script ──"]

    _min = st.session_state.get("ps_1d_min")
    _max = st.session_state.get("ps_1d_max")
    _steps = int(st.session_state.get("ps_1d_steps", 5))
    _log = st.session_state.get("ps_1d_spacing", "Linear") == "Logarithmic"

    code += ["", "# 2. Sweep definition", f"meta = {rec.render(meta)}"]
    if _log:
        code.append(f"sweep_values = np.logspace(np.log10({_min}), np.log10({_max}), {_steps})")
    else:
        code.append(f"sweep_values = np.linspace({_min}, {_max}, {_steps})")
    if meta["type"] == "dimension":
        code.append("sweep_values = np.unique(np.clip(np.round(sweep_values).astype(int), 1, None))")

    _label_lit = label.replace('"', "'")
    code += [
        "",
        "# 3. Run the sweep (apply_sweep_parameter is the same function the app uses)",
        "fig, ax = plt.subplots(figsize=(8, 4))",
        "for val in sweep_values:",
        "    c_k, ib_k, ip_k, is_k, mk_k = apply_sweep_parameter(",
        "        val, meta, cfg, initial_B, initial_P, initial_S, model_kwargs)",
        *_repro_prerun_lines("    ", t_prerun),
        "    model = PBIModel(c_k, initial_B=ib_k, initial_P=ip_k, initial_S=is_k, **mk_k)",
        f"    result = solve_ode(model, {_repro_solve_kwargs()})",
        "    total = np.maximum(result.sum_prefixes('B', 'D'), 1.0)  # CFU (culturable: B+D)",
        '    ax.semilogy(result.time, total, label="' + _label_lit + ' = %.2e" % val)',
        f'ax.set(xlabel="Time (h)", ylabel="Total viable (cells/mL)", title="1D sweep: {_label_lit}")',
        "ax.legend(fontsize=7)",
        "plt.show()",
    ]
    _script = "\n".join(code)
    st.session_state["_last_sweep_repro_code"] = _script
    return _script


def generate_dose_sweep_reproduction_code() -> str:
    """Reproduction script for the dose-response sweep: the recorded base config plus a
    loop that rebuilds the per-run dose schedule exactly as the app does (MOI scaling,
    repeat regimens, nominal overrides) and re-solves."""
    phages = st.session_state.get("int_phages", [])
    antibiotics = st.session_state.get("int_antibiotics", [])
    strains = st.session_state.get("int_strains", [])

    # Reconstruct the sweep inputs from the dr_sweep_* widgets (mirrors the page).
    swept_inputs, swept_units, swept_repeat_configs = {}, {}, {}
    for j in range(len(phages)):
        if st.session_state.get(f"dr_sweep_phg_en_{j}"):
            swept_inputs[f"phage_{j}"] = st.session_state.get(f"dr_sweep_phg_series_{j}", "0, 1e3, 1e5, 1e7, 1e9")
            swept_units[f"phage_{j}"] = st.session_state.get(f"dr_sweep_phg_unit_{j}", "PFU (absolute)")
            if st.session_state.get(f"dr_sweep_phg_rep_en_{j}"):
                swept_repeat_configs[f"phage_{j}"] = {
                    "interval": st.session_state.get(f"dr_sweep_phg_rep_int_{j}", 12.0),
                    "count": int(st.session_state.get(f"dr_sweep_phg_rep_count_{j}", 4)),
                    "start": st.session_state.get(f"dr_sweep_phg_rep_start_{j}", 0.0),
                    "route": st.session_state.get(f"dr_sweep_phg_rep_route_{j}", "bolus"),
                    "duration": st.session_state.get(f"dr_sweep_phg_rep_dur_{j}", 0.0),
                }
    for j in range(len(antibiotics)):
        if st.session_state.get(f"dr_sweep_abx_en_{j}"):
            swept_inputs[f"abx_{j}"] = st.session_state.get(f"dr_sweep_abx_series_{j}", "0.5, 1.0, 2.0")
            swept_units[f"abx_{j}"] = "absolute"
            if st.session_state.get(f"dr_sweep_abx_rep_en_{j}"):
                swept_repeat_configs[f"abx_{j}"] = {
                    "interval": st.session_state.get(f"dr_sweep_abx_rep_int_{j}", 12.0),
                    "count": int(st.session_state.get(f"dr_sweep_abx_rep_count_{j}", 4)),
                    "start": st.session_state.get(f"dr_sweep_abx_rep_start_{j}", 0.0),
                    "route": st.session_state.get(f"dr_sweep_abx_rep_route_{j}", "bolus"),
                    "duration": st.session_state.get(f"dr_sweep_abx_rep_dur_{j}", 0.0),
                }
    if not swept_inputs:
        return "# Enable at least one phage/antibiotic dose sweep, then reopen this."

    parsed = {}
    for k, s in swept_inputs.items():
        parsed[k] = parse_comma_separated_series(s)
    padded, _ = pad_vectors(parsed)

    config, iB, iP, iS, mk = build_nominal_config_from_gui()
    rec = st.session_state["_repro_rec"]

    # Swept phages start at zero free phage — the dose delivers them (mirrors the app).
    swept_phage_idx = [int(k.split("_")[1]) for k in padded if k.startswith("phage_")]
    iP_mod = np.array(iP, dtype=float).copy()
    for jj in swept_phage_idx:
        if jj < len(iP_mod):
            iP_mod[jj] = 0.0

    original_doses = list(st.session_state.get("int_doses", []))
    sum_initial_B = float(sum(s["initial_B"] for s in strains))
    t_prerun = st.session_state.get("int_t_prerun", 0.0)

    _extra = {"DoseEvent", "DoseSchedule"}
    if t_prerun:
        _extra.add("stationary_phase_ic")
    imports, code = _repro_base_config_block(
        rec, iB, iP, iS, mk, extra_imports=_extra, initial_P_override=iP_mod)
    code[0:0] = ["# ── pbisim Dose-Response Sweep Reproduction Script ──"]

    code += [
        "",
        "# 2. Sweep definition (padded value series per swept agent)",
        f"padded = {padded!r}",
        f"swept_units = {swept_units!r}",
        f"swept_repeat_configs = {swept_repeat_configs!r}",
        f"original_doses = {original_doses!r}",
        f"sum_initial_B = {sum_initial_B!r}",
        "M = len(next(iter(padded.values())))",
        "",
        "# 3. Run one simulation per dose combination",
        "fig, ax = plt.subplots(figsize=(8, 4))",
        "for k in range(M):",
        "    custom = [nd for nd in original_doses if not (",
        "        (nd['target_type'] == 'phage' and f\"phage_{nd['target_idx']}\" in padded) or",
        "        (nd['target_type'] == 'antibiotic' and f\"abx_{nd['target_idx']}\" in padded))]",
        "    for key, vec in padded.items():",
        "        val = vec[k]",
        "        if swept_units.get(key) == 'MOI (relative to B(0))':",
        "            val = val * sum_initial_B",
        "        tt = 'phage' if key.startswith('phage') else 'antibiotic'",
        "        ti = int(key.split('_')[1])",
        "        if key in swept_repeat_configs:",
        "            rc = swept_repeat_configs[key]",
        "            for r in range(rc['count']):",
        "                custom.append({'time': rc['start'] + r * rc['interval'], 'amount': val,",
        "                               'target_type': tt, 'target_idx': ti, 'route': rc['route'], 'duration': rc['duration']})",
        "        else:",
        "            nom = [d for d in original_doses if d['target_type'] == tt and d['target_idx'] == ti]",
        "            if nom:",
        "                for nd in nom:",
        "                    c = dict(nd); c['amount'] = val; custom.append(c)",
        "            else:",
        "                custom.append({'time': 0.0, 'amount': val, 'target_type': tt,",
        "                               'target_idx': ti, 'route': 'bolus', 'duration': 0.0})",
        "    events = [DoseEvent(time=d['time'], amount=d['amount'],",
        "                        target=('phage' if d['target_type'] == 'phage' else",
        "                                'antibiotic' if d['target_type'] == 'antibiotic' else 'nutrient'),",
        "                        index=d['target_idx'], route=d['route'], duration=d.get('duration', 0.0))",
        "              for d in custom]",
        "    cfg.dose_schedule = DoseSchedule(events) if events else None",
        "    c_k, ib_k, ip_k, is_k, mk_k = cfg, np.array(initial_B, float), np.array(initial_P, float), initial_S, dict(model_kwargs)",
        *_repro_prerun_lines("    ", t_prerun),
        "    model = PBIModel(c_k, initial_B=ib_k, initial_P=ip_k, initial_S=is_k, **mk_k)",
        f"    result = solve_ode(model, {_repro_solve_kwargs()})",
        "    total = np.maximum(result.sum_prefixes('B', 'D'), 1.0)  # CFU (culturable: B+D)",
        '    ax.semilogy(result.time, total, label="run %d" % (k + 1))',
        'ax.set(xlabel="Time (h)", ylabel="Total viable (cells/mL)", title="Dose-response sweep")',
        "ax.legend(fontsize=7)",
        "plt.show()",
    ]
    _script = "\n".join(code)
    st.session_state["_last_sweep_repro_code"] = _script
    return _script


__all__ = [
    'scripting_enabled',
    'faulthandler',
    'copy',
    '_dc',
    'io',
    'json',
    'os',
    're',
    'time',
    'matplotlib',
    'plt',
    'np',
    'pd',
    'st',
    'ModelBuilder',
    'PBIModel',
    'solve_ode',
    'DoseSchedule',
    'DoseEvent',
    'time_to_clearance',
    'time_to_log_reduction',
    'stationary_phase_ic',
    'StrainDefinition',
    'StrainSet',
    'BinaryResistanceGenotypes',
    'BacterialStrain',
    'PhageStrain',
    'Antibiotic',
    'AntibioticDefinition',
    'AntibioticSensitivity',
    'TreatmentArm',
    'SimulationAgent',
    'execute_code',
    'plot_axis_controls',
    'BACTERIAL_BASES',
    'CFU_BASIS_LABEL',
    'CFU_PREFIXES',
    'bacterial_total',
    'apply_axis_plotly',
    'build_series',
    'series_selector',
    'plot_series',
    'plot_sweep_traces',
    'OBSERVABLES',
    'OBS_COMPARTMENTS',
    'obs_prefixes',
    'predicted_observable',
    'normalize_fit_dataframe',
    'arm_covariate_values',
    'parse_dose_rows',
    'DOSE_TARGETS',
    'apply_row_filters',
    'aggregate_observations',
    'fit_residual',
    'residual_vector_log10',
    'model_information_criteria',
    'compare_fit_models',
    'build_fit_spec',
    'config_param_snapshot',
    'IIV_PARAMETERS',
    'run_trial_simulation',
    'plot_kaplan_meier_plotly',
    'plot_metric_distributions_plotly',
    'plot_pkpd_trajectories_plotly',
    'build_regimen_doses',
    'get_sweep_parameters',
    'categorize_sweep_params',
    'apply_sweep_parameter',
    'parse_comma_separated_series',
    'pad_vectors',
    '_orig_number_input',
    '_number_input_precise',
    'DOSE_AMOUNT_DEFAULTS',
    'DOSE_AMOUNT_LABELS',
    'render_regimen_config',
    'render_iiv_config',
    'render_mutation_graph_editor',
    'mutation_matrix_from_transitions',
    'render_rate_graph_editor',
    'rate_matrix_from_graph',
    'graph_dict_from_edges',
    'brg_genotype_labels',
    'GENERATOR_GRAPH_KEYS',
    'generator_matrices_from_state',
    'apply_generator_matrices',
    'edges_from_rate_matrix',
    'write_back_generator_matrices',
    'read_uploaded_csv',
    'calibration_processed',
    'arm_dose_events',
    'arm_regimen_summary',
    '_init_app_state',
    '_next_uid',
    '_carry_prerun_debris',
    '_safe_od',
    'phage_pk_enabled',
    'central_phage_total',
    'peripheral_phage_total',
    'phage_uses_two_compartment',
    '_sweep_summary_tiles',
    'counted_number_input',
    'load_preset_to_state',
    'configure_summary',
    'apply_ai_configuration',
    'summarize_current_results',
    'SCENARIO_SCHEMA_VERSION',
    '_SCENARIO_EXTRA_DATA_KEYS',
    '_ADS_DATA_RE',
    '_SCENARIO_WIDGET_PREFIXES',
    '_is_scenario_data_key',
    '_json_safe',
    'dump_state_to_scenario',
    'load_scenario_to_state',
    'export_scenarios_json',
    'import_scenarios_json',
    'MODEL_SCHEMA_VERSION',
    'WORKING_DRAFT_LABEL',
    'DEMO_MODELS',
    '_is_model_data_key',
    'dump_model',
    'apply_model_to_state',
    'model_config_context',
    'build_config_from_model',
    'render_model_snapshot',
    'resolve_model_snapshot',
    'model_options',
    'page_model_selector',
    'PARTS_SCHEMA_VERSION',
    'PART_CATEGORIES',
    'PART_SOURCES',
    'empty_parts_library',
    'export_parts_json',
    'import_parts_json',
    '_ENTITY_WIDGET_PREFIXES',
    'clear_entity_widgets',
    'DEFAULT_SCENARIO',
    'SIGNAL_OPTIONS',
    'GROWTH_SIGNALS',
    'DEATH_SIGNALS',
    'LYSIS_SIGNALS',
    'lysis_signal_function',
    'lysis_kwargs',
    'canonical_signal',
    'compat_dormancy_signal',
    'growth_nutrient_kwargs',
    'growth_tracks_nutrients',
    'sequential_growth_phase_params',
    'validate_sequential_growth',
    'dormancy_signal_functions',
    'diffusion_signal_functions',
    'apply_diffusion_signal',
    'set_diffusion_functions',
    'mode_dormancy_kwargs',
    'death_signal_function',
    'death_kwargs',
    '_ReproRecorder',
    'build_nominal_config_from_gui',
    'run_sim_from_gui_params',
    'prerun_collapse_fraction',
    'warn_if_prerun_collapses',
    'reseed_widget_config',
    'save_widget_config',
    'generate_reproduction_code',
    '_repro_base_config_block',
    '_repro_solve_kwargs',
    '_repro_prerun_lines',
    'generate_param_sweep_reproduction_code',
    'generate_dose_sweep_reproduction_code',
]
