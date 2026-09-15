"""Rendered by app.py when this page is selected."""
import re
import multiprocessing as _mp
import time as _time

from pbisim_app.common import *  # noqa: F401,F403

# Covariate-link forms (mirror pbisim-fit's covariates.FORMS).
_COV_FORMS = ("power", "linear", "exponential")


def _apply_arm_covariates(config, cov):
    """Return a per-arm config with covariate links applied, given the arm's covariate
    values ``cov`` ({name: value}). No-op (returns ``config``) unless the config carries
    ``covariate_effects`` (only the post-fit fitted config does), so the manual overlay is
    unchanged. Deep-copies so arms don't stomp each other's scaled paths."""
    if not getattr(config, "covariate_effects", None):
        return config
    try:
        import copy as _copy
        from pbisim_fit.refinement.covariates import apply_covariate_effects

        class _Arm:  # minimal arm: covariates dict (+ empty dose_events for the moi reader)
            def __init__(self, c):
                self.covariates = dict(c or {})
                self.dose_events = []

        return apply_covariate_effects(_copy.deepcopy(config), _Arm(cov), None)
    except Exception:
        return config

# Widget-key prefixes owned by the shared model builder (render_model_builder). These
# are WIDGET keys only — NOT model-data keys like ``ads_<i>_<j>`` (pairwise adsorption)
# or ``direct_phg_res_rates`` (mutation list), which must survive. Popping these on
# fit-apply forces the builder inputs to re-seed from the updated int_* dicts.
_BUILDER_WIDGET_PREFIXES = (
    "str_", "ss_", "phg_", "brg_", "trans_", "dir_trans_", "gen_", "direct_mu_",
    "ads_input_", "ads_dorm_input_",
    "widget_builder_mode", "widget_density_total_cells", "widget_brg_",
)


def _pf(s):
    """Parse a possibly-scientific-notation string to float. Blank/invalid → None.
    (data_editor NumberColumn won't accept typed `1e9`; TextColumn + this does.)"""
    if s is None:
        return None
    t = str(s).strip()
    if t == "" or t.lower() in ("nan", "none"):
        return None
    try:
        return float(t)
    except (TypeError, ValueError):
        return None


def _start_fit_process(payload):
    """Launch the NLS fit in a SEPARATE PROCESS and return ``(process, queue)``.

    Uses the ``fork`` start method: the child inherits the parent's memory, so the
    (possibly non-picklable) ``payload`` — config, dataset, covariate resolvers — is NOT
    pickled to pass it (only the *result* is pickled back through the queue). ``fork`` is
    the Linux/Render default and, unlike ``spawn``/``forkserver``, does not re-import the
    Streamlit app script in the child. Daemonised so it can't outlive the server; the fit
    passes ``n_arm_jobs=1`` so it never spawns nested workers (a daemon process can't)."""
    from pbisim_app.fit_process import run_fit_in_process
    # Pre-import the fit machinery in the PARENT so the forked child runs NO imports —
    # the import lock is the most common fork-in-a-multithreaded-process deadlock source,
    # and the child inherits already-imported modules. (Best-effort; harmless if it fails.)
    try:
        import pbisim_app.nls_fit  # noqa: F401
        import pbisim_fit.refinement.nls  # noqa: F401
    except Exception:
        pass
    _ctx = _mp.get_context("fork")
    _q = _ctx.Queue()
    _proc = _ctx.Process(target=run_fit_in_process, args=(_q, payload), daemon=True)
    _proc.start()
    return _proc, _q


def _apply_map_to_state(map_dict):
    """Write pbisim-fit MAP path→value results back into the live model dicts /
    session keys, so a fit updates the Interactive-Simulator model in place.

    Handled paths: growth_rates[i], bacteria_to_resource_ratio[i], death_rate_B[i]
    → int_strains[i]; burst_sizes[i,j], latent_periods[i,j], phage_decay_rates[j]
    → int_phages[j]; adsorption_rates[i,j] → ads_{i}_{j} key (+ phage adsorption_s);
    monod_constant → int_monod_constant.
    """
    strains = st.session_state.get("int_strains", [])
    phages = st.session_state.get("int_phages", [])
    for path, val in map_dict.items():
        val = float(val)
        m1 = re.fullmatch(r"([a-zA-Z_]+)\[(\d+)\]", path)
        m2 = re.fullmatch(r"([a-zA-Z_]+)\[(\d+),(\d+)\]", path)
        if path == "monod_constant":
            st.session_state["int_monod_constant"] = val
        elif m1:
            name, i = m1.group(1), int(m1.group(2))
            if name == "growth_rates" and i < len(strains):
                strains[i]["growth_rate"] = val
            elif name in ("bacteria_to_resource_ratio", "death_rate_B") and i < len(strains):
                strains[i][name] = val
            elif name == "phage_decay_rates" and i < len(phages):
                phages[i]["phage_decay_rates"] = val
        elif m2:
            name, i, j = m2.group(1), int(m2.group(2)), int(m2.group(3))
            if name == "adsorption_rates":
                st.session_state[f"ads_{i}_{j}"] = val
                if j < len(phages):
                    phages[j]["adsorption_rates"] = val
                    phages[j]["adsorption_s"] = val
            elif name in ("burst_sizes", "latent_periods") and j < len(phages):
                phages[j][name] = val
    st.session_state["int_strains"] = strains
    st.session_state["int_phages"] = phages


def _apply_config_to_state(cfg):
    """Write a fitted pbisim ModelConfig's resolved biological parameters back into the
    live builder dicts / session keys. Robust for reparameterized fits (reads final
    values off the config rather than parsing theta names)."""
    strains = st.session_state.get("int_strains", [])
    phages = st.session_state.get("int_phages", [])

    def _a1(x):
        return np.atleast_1d(np.asarray(x, dtype=float)) if x is not None else np.array([])

    _gr, _rr = _a1(getattr(cfg, "growth_rates", None)), _a1(getattr(cfg, "bacteria_to_resource_ratio", None))
    _db = _a1(getattr(cfg, "death_rate_B", None))
    for _i, _s in enumerate(strains):
        if _i < len(_gr): _s["growth_rate"] = float(_gr[_i])
        if _i < len(_rr): _s["bacteria_to_resource_ratio"] = float(_rr[_i])
        if _i < len(_db): _s["death_rate_B"] = float(_db[_i])
    _bs = np.atleast_2d(np.asarray(getattr(cfg, "burst_sizes", [[]]), dtype=float))
    _lp = np.atleast_2d(np.asarray(getattr(cfg, "latent_periods", [[]]), dtype=float))
    _ads = np.atleast_2d(np.asarray(getattr(cfg, "adsorption_rates", [[]]), dtype=float))
    _pdr = _a1(getattr(cfg, "phage_decay_rates", None))
    for _j, _p in enumerate(phages):
        if _bs.shape[1] > _j: _p["burst_sizes"] = float(_bs[0, _j])
        if _lp.shape[1] > _j: _p["latent_periods"] = float(_lp[0, _j])
        if _j < len(_pdr): _p["phage_decay_rates"] = float(_pdr[_j])
        if _ads.shape[1] > _j:
            _p["adsorption_rates"] = float(_ads[0, _j]); _p["adsorption_s"] = float(_ads[0, _j])
    for _i in range(_ads.shape[0]):
        for _j in range(_ads.shape[1]):
            st.session_state[f"ads_{_i}_{_j}"] = float(_ads[_i, _j])
    if getattr(cfg, "monod_constant", None) is not None:
        st.session_state["int_monod_constant"] = float(cfg.monod_constant)
    # Estimated initial conditions (fit_initial_cfu / fit_initial_pfu) → the model's
    # inoculum, so a re-run / the Interactive Simulator matches the fit.
    _ic = getattr(cfg, "fit_initial_cfu", None)
    if _ic is not None and np.isfinite(_ic) and _ic > 0 and strains:
        _tot = sum(float(s.get("initial_B", 0.0)) for s in strains)
        for _s in strains:
            _s["initial_B"] = (float(_s.get("initial_B", 0.0)) * float(_ic) / _tot
                               if _tot > 0 else float(_ic) / len(strains))
    _ip = getattr(cfg, "fit_initial_pfu", None)
    if _ip is not None and np.isfinite(_ip) and _ip > 0 and phages:
        phages[0]["initial_P"] = float(_ip)
    st.session_state["int_strains"] = strains
    st.session_state["int_phages"] = phages
    # Generator matrices (mutation / phenotypic-transition graphs) — rewrite the edge
    # lists the builder was already using so estimated rates show up in the editors.
    _brg = st.session_state.get("int_builder_mode") == "Binary Genotypes (BRG)"
    _bn = (brg_genotype_labels(len(phages), len(st.session_state.get("int_antibiotics", [])))
           if _brg else [s["name"] for s in strains])
    write_back_generator_matrices(cfg, _bn, [p["name"] for p in phages], include_mutation=not _brg)


def _broadcast_growth(inits, n_bacteria):
    """Curve stripping estimates ONE growth rate (g_max = the control's growth). Map it to ALL
    bacterial species (growth_rates[i]), not only growth_rates[0], so every strain's growth is
    set/seeded — the resistant genotype gets the same growth absent a per-strain estimate. Other
    stripped params (cap, cap_r, adsorption, …) stay per-index."""
    out = dict(inits)
    g = out.get("growth_rates[0]")
    if g is not None:
        for i in range(int(n_bacteria)):
            out[f"growth_rates[{i}]"] = float(g)
    return out


def _split_initial_B_by_f0(strains, f0):
    """Map a resistant fraction f0 onto the first two strains' ``initial_B`` (strain 0 = WT /
    sensitive, strain 1 = resistant — the 2-genotype stripping convention), preserving their
    combined total. In Direct/StrainSet the resistant fraction IS this split (there's no
    ``init_resistant_fraction`` ModelConfig field). Returns True if applied."""
    if f0 is None or len(strains) < 2:
        return False
    f0 = float(min(max(float(f0), 0.0), 1.0))
    tot = float(strains[0].get("initial_B", 0.0)) + float(strains[1].get("initial_B", 0.0))
    if tot <= 0:
        return False
    strains[0]["initial_B"] = (1.0 - f0) * tot
    strains[1]["initial_B"] = f0 * tot
    return True


def _ensure_resistant_genotype(inits, mode):
    """A curve-strip / amortized estimate is a 2-GENOTYPE model (WT + resistant: the resistant
    strain does NOT adsorb the phage — its defining trait — and carries its own capacity `cap_r`
    and initial fraction `f0`). In Direct/StrainSet, make sure the builder has a resistant strain
    at index 1 — add a WT clone if it has only one — and mark it resistant (pairwise adsorption 0
    to every phage), so `cap_r` / `f0` / the resistant growth have somewhere to land. BRG is
    inherently 2-genotype (its resistant genotype is derived from the fitness cost), so it's left
    alone. Returns a short list of what it did (for the apply message)."""
    if mode == "Binary Genotypes (BRG)":
        return []
    if not ("bacteria_to_resource_ratio[1]" in inits or "init_resistant_fraction" in inits):
        return []
    _strns = st.session_state.get("int_strains", [])
    if not _strns:
        return []
    _did = []
    if len(_strns) < 2:
        import copy as _copy
        _res = _copy.deepcopy(_strns[0])
        _res["name"] = str(_strns[0].get("name", "WT")) + " (resistant)"
        st.session_state["int_strains"] = _strns + [_res]
        _did.append("added a resistant strain")
    # Resistant genotype (strain 1) does not adsorb the phage → set its pairwise adsorption to 0.
    _nph = max(len(st.session_state.get("int_phages", [])), 1)
    for _j in range(_nph):
        st.session_state[f"ads_1_{_j}"] = 0.0
    _did.append("resistant adsorption = 0")
    return _did


def _config_structure(config, extra=None):
    """Structural function CHOICES (not just numbers) to mirror from a fitted config into the
    builder, so an applied model reproduces the fit: the lysis progression function (e.g.
    ``frac_lysis`` — the app default is ``constant_lysis``), the growth-signal function, the
    Monod Ks, and the infected-cell nutrient draw. Returns ``{session_key: value}`` (plus a
    special ``"_inc"`` key applied per-phage). ``extra`` may carry ``{"inc": value}`` from the
    estimate's MAP when the config's own field is a placeholder."""
    _s = {}
    _lf = getattr(config, "lysis_progression_function", None)
    if getattr(_lf, "__name__", None):
        _s["int_lysis_function"] = _lf.__name__                 # e.g. "frac_lysis"
    _gf = getattr(config, "growth_function", None)
    if getattr(_gf, "__name__", None):
        _s["int_growth_function"] = _gf.__name__
    elif getattr(config, "track_nutrients", False):
        _s["int_growth_function"] = "monod_growth"              # None + nutrient → Monod
    _mc = getattr(config, "monod_constant", None)
    if _mc is not None:
        _s["int_monod_constant"] = float(_mc)
    _inc = (extra or {}).get("inc")
    if _inc is None:
        try:
            _inc = float(np.max(np.asarray(getattr(config, "infected_nutrient_consumption", 0.0))))
        except Exception:
            _inc = None
    if _inc is not None and _inc > 0:
        _s["_inc"] = float(_inc)
    return _s


def _apply_initials_to_builder(inits, structure=None):
    """Write estimate initials (``{config_path: value}``) into the Working-draft builder (the
    manual-calibration values), reusing the fit-apply write path. Shared by curve-strip Apply
    and amortized-ONNX Apply. Ensures a resistant genotype exists (2-genotype estimate), mirrors
    the fitted model's structural function choices (``structure``, from :func:`_config_structure`
    — lysis/growth signal, Monod Ks, infected-nutrient draw), broadcasts g_max to all species,
    mirrors BRG base growth/cap, maps f0 to the WT/resistant B0 split, and enables + populates the
    OD/debris module. Builder widgets are popped here; the caller recomputes the overlay. Returns
    ``(applied, skipped)`` path lists."""
    from pbisim_app.nls_fit import set_config_path
    _mode = st.session_state.get("int_builder_mode", "Direct (ModelBuilder)")
    # Ensure the 2-genotype structure BEFORE building the config, so cap_r / f0 have a strain 1.
    _extra = _ensure_resistant_genotype(inits, _mode)
    # Mirror the fitted model's structural function choices BEFORE build_nominal so the rebuilt
    # config uses them (a param that exists in the fit but whose function/module isn't selected
    # in the builder must still take effect).
    _struct = dict(structure or {})
    for _key in ("int_lysis_function", "int_growth_function", "int_monod_constant"):
        if _struct.get(_key) is not None:
            st.session_state[_key] = _struct[_key]
            _extra.append(_key.replace("int_", "") + " → " + str(_struct[_key]))
    if _struct.get("_inc") is not None:
        _incv = float(_struct["_inc"])
        st.session_state["int_infected_nutrient_consumption"] = _incv       # legacy scalar seed
        for _ph in st.session_state.get("int_phages", []):
            _ph["infected_nutrient_consumption"] = _incv                    # per-phage (authoritative)
        _extra.append("infected consumption → %.3g" % _incv)
    _cfg0 = build_nominal_config_from_gui()[0]
    _ai = _broadcast_growth(inits, getattr(_cfg0, "n_bacteria", 1))
    _ap, _sk = [], []
    for _p, _v in _ai.items():
        (_ap if set_config_path(_cfg0, _p, _v) else _sk).append(_p)
    _apply_config_to_state(_cfg0)
    # Debris rates live on session keys, and the OD/debris module is FEATURE-GATED
    # (int_debris_enabled). A fitted param that exists in the estimate but whose module isn't
    # pre-activated in the builder must still take effect — so enable the module when applying
    # its params (the general rule; debris is the gated case in the strip/amortized param set).
    if "debris_v" in _ai or "debris_kdis" in _ai:
        st.session_state["int_debris_enabled"] = True
        if "debris_v" in _ai:
            st.session_state["int_debris_v"] = float(_ai["debris_v"])
        if "debris_kdis" in _ai:
            st.session_state["int_debris_kdis"] = float(_ai["debris_kdis"])
        st.session_state.setdefault("int_debris_u", 0.0)
        _extra.append("enabled OD/debris")
    if _mode == "Binary Genotypes (BRG)":
        if "growth_rates[0]" in _ai:
            st.session_state["int_brg_base_growth"] = float(_ai["growth_rates[0]"])
        if "bacteria_to_resource_ratio[0]" in _ai:
            st.session_state["int_brg_base_ratio"] = float(_ai["bacteria_to_resource_ratio[0]"])
    else:
        _f0v = _ai.get("init_resistant_fraction")
        _strns = st.session_state.get("int_strains", [])
        if _split_initial_B_by_f0(_strns, _f0v):
            st.session_state["int_strains"] = _strns
            if "init_resistant_fraction" in _sk:
                _sk.remove("init_resistant_fraction")
            _ap.append("init_resistant_fraction (B₀ split)")
    for _k in [k for k in list(st.session_state.keys())
               if k.startswith("fit_edit_") or k.startswith(_BUILDER_WIDGET_PREFIXES)]:
        st.session_state.pop(_k, None)
    return _ap + _extra, _sk


def _compute_overlay(config, iB, iP, iS, mk, ctx):
    """Simulate `config` once per arm and project every selected observable into a
    small-multiples overlay-vs-data result dict (also computes per-obs RMSE/R² and the
    pooled log10 combined objective J). Shared by the manual overlay button and the
    post-fit fitted-curve overlay, so both draw the same way. `ctx` bundles the arm/
    observable selection, aggregated data, and plotting options."""
    _B0 = float(np.sum(iB))
    _panels, _metrics = {}, []
    for _arm in ctx["sel_arms"]:
        _cond = ctx["arm_cond"].get(_arm, {})
        _arm_b0 = float(_cond.get("b0", _B0)) or _B0
        _arm_prerun = float(_cond.get("prerun", 0.0) or 0.0)
        _moi = float(_cond.get("moi", ctx["conds"].get(_arm, {}).get("moi", 0.0)))
        _armB = iB * (_arm_b0 / _B0) if _B0 > 0 else iB
        _armP = np.zeros(len(iP))
        if len(iP):
            # dose value is absolute PFU/mL (pfu) or a multiple of B₀ (moi); a per-arm
            # unit (from an imported dose record) overrides the global default.
            _du = _cond.get("moi_unit", ctx.get("dose_unit"))
            _armP[0] = _moi if _du == "pfu" else _moi * _arm_b0
        # Per-arm covariate scaling: a fitted config carries covariate_effects, so each
        # arm's model parameters are modulated by that arm's covariate (θ_i = θ_ref·f(cov)).
        # No-op for the manual overlay (base config has no covariate_effects).
        _cfg_arm = _apply_arm_covariates(config, ctx.get("arm_covariates", {}).get(_arm, {}))
        _mk_arm = dict(mk)
        _iS_arm = iS
        if _arm_prerun > 0:
            _ic = stationary_phase_ic(_cfg_arm, t_prerun=_arm_prerun, B0=_armB, initial_S=_iS_arm)
            _armB = _ic.B
            _iS_arm = max(float(_ic.S), 0.0)
            if _ic.D is not None:
                _mk_arm["initial_D"] = _ic.D
            if _ic.Imm is not None:
                _mk_arm["initial_Imm"] = _ic.Imm
            _carry_prerun_debris(_ic, _mk_arm)
        _m = PBIModel(_cfg_arm, initial_B=_armB, initial_P=_armP, initial_S=_iS_arm, **_mk_arm)
        _r = solve_ode(_m, t_end=ctx["t_end"], dt=0.25, method=ctx["method"],
                       extinction_threshold=ctx["thr"])
        for _ok in ctx["sel_obs"]:
            _sp = OBSERVABLES.get(_ok, {"log": True, "link": None, "label": _ok})
            _umo = (_ok == "od" and ctx["debris_on"])
            _px = obs_prefixes(_ok, ctx.get("obs_compartments"))
            _pred = predicted_observable(_r, _ok, ctx["link_vals"].get(_ok),
                                         use_model_od=_umo, prefixes=_px)
            _d = ctx["agg"][(ctx["agg"]["arm"] == _arm) & (ctx["agg"]["observable"] == _ok)].sort_values("time")
            if not len(_d):
                continue
            _has_band = ctx["band"] is not None and _d["lo"].notna().any()
            _pan = _panels.setdefault(_ok, {"obs": _ok, "label": _sp.get("label", _ok),
                                            "log": bool(_sp.get("log")), "series": []})
            _pan["series"].append({
                "label": _arm, "time": np.asarray(_r.time), "pred": np.asarray(_pred),
                "obs_time": _d["time"].to_numpy(), "obs_value": _d["value"].to_numpy(),
                "obs_lo": _d["lo"].to_numpy() if _has_band else None,
                "obs_hi": _d["hi"].to_numpy() if _has_band else None,
                "is_raw": ctx["stat_key"] == "raw",
            })
            _resid = residual_vector_log10(_r.time, _pred, _d["time"].values,
                                           _d["value"].values, _sp.get("floor_log10", 1.0))
            _metrics.append({"observable": _sp.get("label", _ok), "group": _arm,
                             "B₀": _arm_b0, "pre-run (h)": _arm_prerun,
                             ctx.get("dose_label", "MOI"): _moi,
                             "n_points": len(_d),
                             "RMSE (log₁₀)": float(np.sqrt(np.mean(_resid ** 2))) if _resid.size else float("nan")})

    _panel_list, _all_resid = [], []
    for _ok, _pan in _panels.items():
        _floor = 10.0 ** float(OBSERVABLES.get(_ok, {}).get("floor_log10", 1.0))
        _oa, _pa = [], []
        for _s in _pan["series"]:
            _pt = np.interp(_s["obs_time"], _s["time"], _s["pred"])
            _ov = np.asarray(_s["obs_value"], dtype=float)
            _oa.append(np.log10(np.maximum(_ov, _floor)))
            _pa.append(np.log10(np.maximum(_pt, _floor)))
        _oa = np.concatenate(_oa) if _oa else np.array([])
        _pa = np.concatenate(_pa) if _pa else np.array([])
        _msk = np.isfinite(_oa) & np.isfinite(_pa)
        _oa, _pa = _oa[_msk], _pa[_msk]
        _pan["rmse"] = float(np.sqrt(np.mean((_oa - _pa) ** 2))) if _oa.size else float("nan")
        _sst = float(np.sum((_oa - _oa.mean()) ** 2)) if _oa.size else 0.0
        _pan["r2"] = (1.0 - float(np.sum((_oa - _pa) ** 2)) / _sst) if _sst > 0 else float("nan")
        _pan["n"] = int(_oa.size)
        _panel_list.append(_pan)
        if _oa.size:
            _all_resid.append(_oa - _pa)

    _all = np.concatenate(_all_resid) if _all_resid else np.array([])
    _stat_label = ctx["stat"] if ctx["stat_key"] != "raw" else "raw points"
    return {
        "panels": _panel_list, "metrics": _metrics,
        "combined": float(np.sqrt(np.mean(_all ** 2))) if _all.size else float("nan"),
        # Pooled log10 residual vector (the AIC/BIC model-comparison input) + its size.
        "residuals": _all.tolist(), "n_resid": int(_all.size),
        "stat_label": _stat_label,
        "title": ctx.get("title") or (f"Model vs observations ({_stat_label}"
                 + (f" + {ctx['band_choice']} band)" if ctx["band"] else ")")),
    }


def _harvest_fit_done(_job):
    """Read a finished fit's result (``_job['result']`` = the picklable dict the fit
    subprocess returned: ``{'map','ci','fitted_config'}``) into session_state:
    ``calib_fit_result`` (MAP + CI table), ``calib_fitted_config``, and the fitted-curve
    ``calib_overlay_result``. Returns True on success; on failure records
    ``calib_fit_error`` and returns False (so a fit that converged to an unreadable
    result surfaces a clean message instead of a raw traceback)."""
    _res = _job["result"]
    try:
        _map = _res["map"]
        _ci = _res.get("ci") or {}
        _mapped = {m["path"] for m in _job["mappings"]}
        _pl = _job["path_label"]
        _pmeta = [{"label": _pl.get(t["path"], t["path"]), "key": f"free{k}"}
                  for k, t in enumerate(_job["targets"])
                  if t["free"] and t["path"] not in _mapped]
        _pmeta += [{"label": f"θ {th['name']}", "key": th["name"]}
                   for th in _job["thetas"]]
        st.session_state["calib_fit_result"] = {
            "map": {k: float(v) for k, v in _map.items()},
            "ci": {k: [float(a), float(b)] for k, (a, b) in _ci.items()},
            "params": _pmeta, "model": _job["fit_model"],
        }
        _fcfg = _res["fitted_config"]
        st.session_state["calib_fitted_config"] = _fcfg
        try:
            # If B₀ / initial phage were estimated (fit_initial_cfu/pfu), reflect
            # them in the overlay's inoculum so the fitted curve matches the data.
            _fB2, _fP2 = _job["fB"], _job["fP"]
            _fovl = dict(_job["ovl_ctx"])
            _ic = getattr(_fcfg, "fit_initial_cfu", None)
            if _ic is not None and np.isfinite(_ic) and _ic > 0:
                if float(np.sum(_fB2)) > 0:
                    _fB2 = _fB2 * (float(_ic) / float(np.sum(_fB2)))
                _fovl["arm_cond"] = {a: {**c, "b0": float(_ic)}
                                     for a, c in _job["ovl_ctx"]["arm_cond"].items()}
            _ip = getattr(_fcfg, "fit_initial_pfu", None)
            if _ip is not None and np.isfinite(_ip) and _ip > 0 and len(_fP2):
                _fP2 = np.asarray(_fP2, dtype=float).copy(); _fP2[0] = float(_ip)
            _fovl["title"] = "Fitted model vs observations (NLS MAP)"
            if _job["od_link"]:
                _fovl["link_vals"] = dict(_job["ovl_ctx"]["link_vals"], od=_job["od_link"])
            st.session_state["calib_overlay_result"] = _compute_overlay(
                _fcfg, _fB2, _fP2, _job["fS"], _job["fmk"], _fovl)
        except Exception:
            pass
        return True
    except Exception as _he:  # noqa: BLE001 — surface, don't crash the page
        st.session_state["calib_fit_error"] = (
            f"The fit finished but its result couldn't be read: {type(_he).__name__}: {_he}")
        return False


def _monitor_fit_job():
    """Observe the background NLS fit at the TOP of the page — decoupled from
    section 5c's parameter-table code.

    The fit runs in a separate PROCESS (``_start_fit_process``) that returns its result on
    a ``queue``; the process has its own interpreter/GIL, so the Streamlit app stays fully
    responsive and navigable while it runs. This function only DISPLAYS and HARVESTS it, at
    the top of ``render()`` — outside section 5c's ``try`` — so a completed fit is picked up
    the moment the user returns and the banner can't be hidden by an unrelated error.

    **No auto-poll.** Earlier thread-based versions polled with a whole-app
    ``sleep``+``st.rerun`` loop (starved the script runner → froze + nav desync) and then an
    ``st.fragment(run_every=…)`` (timer kept re-triggering a whole-app rerun after nav →
    froze again). Both fight Streamlit's model, and a *thread* also couldn't keep the app
    responsive (GIL). Now: the fit is a process (responsive), and the running banner is
    STATIC with a manual **Refresh** button. Navigation is a plain rerun (nothing to
    starve); the result is harvested on the next render — when the user clicks Refresh or
    simply returns to this page. The queue is checked here each render."""
    _err = st.session_state.pop("calib_fit_error", None)
    if _err:
        st.error(_err)
    _job = st.session_state.get("fit_job")
    if _job is None:
        return
    if _job.get("status") != "running":
        st.session_state["fit_job"] = None
        return

    # Non-blocking check of the fit process: harvest a result off the queue, or detect a
    # process that died without producing one.
    _q, _proc = _job.get("queue"), _job.get("proc")
    _finished = False
    try:
        if _q is not None and not _q.empty():
            _kind, _data = _q.get_nowait()
            if _kind == "done":
                _job["result"] = _data
                _harvest_fit_done(_job)
            else:
                st.session_state["calib_fit_error"] = f"Fit failed: {_data}"
            _finished = True
        elif _proc is not None and not _proc.is_alive():
            # Exited without a result on the queue (crash / OOM-killed / segfault).
            st.session_state["calib_fit_error"] = (
                "The fit process exited unexpectedly without a result (it may have run out "
                "of memory or crashed). Try fewer restarts / max-evaluations, or a simpler model.")
            _finished = True
    except Exception as _pe:  # noqa: BLE001 — queue/process error, don't crash the page
        st.session_state["calib_fit_error"] = f"Could not read the fit result: {_pe}"
        _finished = True

    if _finished:
        if _proc is not None:
            try:
                _proc.join(timeout=1)
            except Exception:
                pass
        st.session_state["fit_job"] = None
        st.rerun()
        return

    # Still running — static banner + manual Refresh + Stop (no auto-rerun loop).
    _el = _time.time() - _job.get("t0", _time.time())
    _pnames = ", ".join(_job.get("param_preview", []))
    st.info(
        f"⏳ Fitting **{len(_job.get('param_preview', []))} parameter(s)** in a background "
        f"process (started ~{_el:0.0f}s ago). The app stays responsive — **navigate away and "
        f"come back** freely; the fit keeps running and the result will be waiting. Click "
        f"**Refresh status** to check on it here."
        + (f"  \nEstimating: {_pnames}" if _pnames else ""))
    _rc1, _rc2 = st.columns(2)
    if _rc1.button("🔄 Refresh status", key="fit_refresh", width="stretch"):
        st.rerun()   # single rerun (not a loop): re-checks the queue / harvests if done
    if _rc2.button("Stop fit", key="fit_stop", width="stretch"):
        if _proc is not None and _proc.is_alive():
            try:
                _proc.terminate()
            except Exception:
                pass
        st.session_state["fit_job"] = None
        st.warning("Fit stopped — result discarded. (Fix the bounds and re-run.)")
        st.rerun()


def render():
    theme_mode = st.session_state.get("theme_mode", "Light")
    st.title("Calibration — data overlay")
    st.caption(
        "Upload experimental data and overlay the **current model's** prediction (configured in "
        "the Interactive Simulator) on the observations. Tune parameters there to match; a "
        "manual-tuning panel and the pbisim-fit hand-off come next."
    )

    # Re-seed the Calibration widgets from a persistent config BEFORE they render.
    # Streamlit drops a widget's key from session_state whenever the widget isn't
    # rendered on a rerun, so navigating to the Simulator (to change the model) and
    # back would otherwise reset the filters / grouping / statistics. `fit_config`
    # is a plain (non-widget) key, so it survives.
    # Buttons and the file-uploader can't be re-seeded via session_state, so they
    # are never persisted; everything else (filters/grouping/statistics/overlay
    # selections) is.
    # Buttons + the file-uploader must never be shadowed into fit_config: re-seeding a
    # button's value pre-sets it, which makes the later st.button() raise. (Text/number
    # widgets are fine to persist.)
    _FIT_NOPERSIST = {"fit_csv", "fit_config", "fit_dataset", "fit_overlay", "fit_clear",
                      "fit_load", "fit_save_scenario", "fit_run_nls", "fit_apply_map",
                      "fit_model_sel", "fit_job", "fit_stop", "fit_refresh", "fit_tbl_reset",
                      "fit_spec_from_tables", "fit_spec_to_tables", "fit_spec_rev",
                      "fit_share_go", "fit_cmp_add", "fit_cmp_clear",
                      "strip_compute", "strip_seed", "strip_apply"}
    _fcfg = st.session_state.setdefault("fit_config", {})
    for _wk, _wv in list(_fcfg.items()):
        if _wk in _FIT_NOPERSIST:
            _fcfg.pop(_wk, None)  # scrub any stale non-persistable key
            continue
        if _wk not in st.session_state:
            try:
                st.session_state[_wk] = _wv
            except Exception:
                pass  # widget type refuses assignment (e.g. a button) — skip it

    # Background-fit monitor — rendered FIRST, before the (fragile, non-persisted)
    # section-5c code, so a running fit's progress and a finished fit's result are always
    # shown/harvested regardless of the rest of the page's state. The fit runs in a separate
    # PROCESS (own GIL → app stays responsive/navigable); a running fit shows a STATIC banner
    # + manual Refresh (no auto-poll loop). The result is harvested off the process queue on
    # the next render (Refresh / return).
    _monitor_fit_job()

    # ── 1. Upload + column mapping ───────────────────────────────────────────
    st.markdown("### 1 · Upload data")
    _up = st.file_uploader("Experimental data (CSV)", type=["csv"], key="fit_csv")
    if _up is not None:
        try:
            _raw = read_uploaded_csv(_up)
        except Exception as e:
            st.error(f"Could not read CSV: {e}")
            _raw = None
        if _raw is not None:
            st.dataframe(_raw.head(8), width="stretch")
            _cols = list(_raw.columns)
            _low = [c.lower() for c in _cols]

            def _guess(cands, default=0):
                for c in cands:
                    if c in _low:
                        return _low.index(c)
                return default

            _canonical = all(k in _low for k in ("time", "arm", "observable", "value"))
            if _canonical:
                st.success("Detected pbisim-fit long format (time, arm, observable, value).")
            with st.expander("Column mapping", expanded=not _canonical):
                _tc = st.selectbox("Time column", _cols, index=_guess(["time"]))
                _vc = st.selectbox("Value (measurement) column", _cols, index=_guess(["value", "dv"]))
                _obs_from_col = st.checkbox("Observable is in a column", value=("observable" in _low))
                if _obs_from_col:
                    _obs = st.selectbox("Observable column", _cols, index=_guess(["observable"]))
                else:
                    _obs = st.selectbox("Observable type (fixed for all rows)", list(OBSERVABLES),
                                        format_func=lambda k: OBSERVABLES[k]["label"])
                _default_arms = [c for c in _cols if c.lower() in ("phage", "moi", "arm", "experi", "experi_num")]
                _ac = st.multiselect("Arm-defining column(s)", _cols, default=_default_arms or ([_cols[0]] if _cols else []))
                # Auto-detect a dose column by substring (dose/pfu/moi/titre) so e.g.
                # `dose_phage_pfu` is picked up without manual selection.
                _dose_idx = next((i for i, c in enumerate(_low)
                                  if any(k in c for k in ("dose", "pfu", "moi", "titre", "titer"))), None)
                _mc = st.selectbox("Phage-dose column (optional — drives the simulated dose per arm)",
                                   ["(none)"] + _cols,
                                   index=(1 + _dose_idx) if _dose_idx is not None else 0)
                _mc = None if _mc == "(none)" else _mc
                _dunit_lbl = st.radio(
                    "Dose unit", ["MOI (× B₀)", "PFU/mL (absolute)"],
                    index=(1 if (_mc and "pfu" in _mc.lower()) else 0), horizontal=True,
                    help="How to read the dose column: MOI multiplies each arm's B₀; "
                         "PFU/mL is an absolute phage titre (e.g. a column like dose_phage_pfu).")
                # NONMEM/Monolix dose rows (optional): a dose-event/EVID column marks
                # dose rows (=1), routed to the compartment named in the observable
                # column. When present these become the per-arm doses (inoculum + phage),
                # overriding the manual per-arm fields for the targets they specify.
                st.markdown("**Dose records (NONMEM/Monolix, optional)**")
                st.caption("Map an EVID column to import interleaved dose rows. Each dose row's "
                           "*observable* names its target compartment (bacteria / phage / "
                           "antibiotic / nutrient); its amount is read from the AMT column.")
                _evid_i = next((i for i, c in enumerate(_low) if c in ("evid", "dose_event", "mdv")), None)
                _evid = st.selectbox("Dose-event column (EVID: 0 = observation, 1 = dose)",
                                     ["(none)"] + _cols,
                                     index=(1 + _evid_i) if _evid_i is not None else 0)
                _evid = None if _evid == "(none)" else _evid
                _amt = _unit_col = None
                if _evid:
                    _amt_i = next((i for i, c in enumerate(_low) if c in ("amt", "amount", "dose")), None)
                    _amt = st.selectbox("Dose amount column (AMT)", ["(none)"] + _cols,
                                        index=(1 + _amt_i) if _amt_i is not None else 0)
                    _amt = None if _amt == "(none)" else _amt
                    _unit_i = next((i for i, c in enumerate(_low) if "unit" in c), None)
                    _unit_col = st.selectbox("Dose unit column (optional; else a per-target default)",
                                             ["(none)"] + _cols,
                                             index=(1 + _unit_i) if _unit_i is not None else 0)
                    _unit_col = None if _unit_col == "(none)" else _unit_col
            if st.button("Load dataset", key="fit_load", width="stretch"):
                st.session_state.fit_dataset = {
                    "raw": _raw, "time": _tc, "value": _vc, "observable": _obs,
                    "arm_cols": _ac, "moi": _mc,
                    "dose_unit": ("pfu" if str(_dunit_lbl).startswith("PFU") else "moi"),
                    "evid": _evid, "amount": _amt, "unit_col": _unit_col,
                }
                st.success(f"Loaded {len(_raw)} rows. Configure grouping / filters / statistics below.")
                st.rerun()

    # ── 2. Filter · group · statistics · overlay ─────────────────────────────
    _ds = st.session_state.get("fit_dataset")
    if not _ds:
        st.info("Upload a dataset above to begin.")
    else:
        _raw = _ds["raw"]
        _cols = list(_raw.columns)
        _tc, _vc, _obs, _mc = _ds["time"], _ds["value"], _ds["observable"], _ds["moi"]
        _dose_unit = _ds.get("dose_unit", "moi")
        _dose_lbl = "PFU/mL" if _dose_unit == "pfu" else "MOI"

        # -- Filters --------------------------------------------------------
        st.markdown("### 2 · Filter rows")
        _filter_cols = st.multiselect(
            "Filter on column(s) (leave a value list empty = include all)", _cols,
            default=[], key="fit_filter_cols",
        )
        _filters = {}
        for _fc in _filter_cols:
            _uniques = sorted(_raw[_fc].dropna().astype(str).unique().tolist())
            _filters[_fc] = st.multiselect(f"Include {_fc} =", _uniques, default=[], key=f"fit_filter_{_fc}")

        # -- Grouping + statistic -------------------------------------------
        st.markdown("### 3 · Grouping & statistics")
        gc1, gc2, gc3 = st.columns(3)
        with gc1:
            _group_cols = st.multiselect("Grouping variables (define the curves/arms)", _cols,
                                         default=_ds["arm_cols"], key="fit_group_cols")
        with gc2:
            _stat = st.selectbox("Statistic over replicates", ["Raw points", "Mean", "Median"], key="fit_stat")
        with gc3:
            _band_choice = st.selectbox("Percentile band", ["None", "10–90", "25–75", "5–95"],
                                        index=0, disabled=(_stat == "Raw points"), key="fit_band")
        _stat_key = {"Raw points": "raw", "Mean": "mean", "Median": "median"}[_stat]
        _band = None if (_band_choice == "None" or _stat_key == "raw") else tuple(int(x) for x in _band_choice.split("–"))

        # Filter → normalise → aggregate (cached; recomputes only when inputs change)
        try:
            _filters_key = tuple((c, tuple(v)) for c, v in _filters.items())
            _filtered, _long, _conds, _agg = calibration_processed(
                _raw, _filters_key, _tc, _vc, _obs, tuple(_group_cols), _mc, _stat_key, _band)
        except Exception as e:
            st.error(f"Could not build the grouped dataset: {e}")
            _filtered, _long, _agg = _raw, None, None
        st.caption(f"{len(_filtered)} / {len(_raw)} rows after filtering.")

        # NONMEM/Monolix dose rows (EVID=1) → per-arm dose records, keyed by the same
        # grouping as the observations. Used to gate the manual per-arm gap-fillers and
        # to emit the doses into the fit dataset.
        try:
            _arm_doses = parse_dose_rows(_filtered, tuple(_group_cols), _ds.get("evid"),
                                         _obs, _ds.get("amount"), _tc, _ds.get("unit_col"))
        except Exception:
            _arm_doses = {}
        if _arm_doses:
            _n_bd = sum(1 for ds in _arm_doses.values() for d in ds if d["target"] == "bacteria")
            _n_pd = sum(1 for ds in _arm_doses.values() for d in ds if d["target"] == "phage")
            st.caption(f"Imported dose records: {_n_bd} bacteria + {_n_pd} phage across "
                       f"{len(_arm_doses)} arm(s) — these override the matching manual fields.")

        if _long is not None and len(_long):
            # key=str guards the sort even if a stale/exotic dtype slips a non-str in
            # (normalize now coerces arm/observable to str, but be defensive).
            _arms = sorted(_long["arm"].astype(str).unique(), key=str)
            _obs_keys = sorted(_long["observable"].astype(str).unique(), key=str)

            # Per-arm covariate values (for covariate-link fitting): numeric grouping
            # columns (constant within an arm) plus the dose-derived MOI. Available names
            # feed the covariate-effects table in the fit section below.
            _num_cov_cols = [c for c in _group_cols if c in _filtered.columns
                             and pd.api.types.is_numeric_dtype(_filtered[c])]
            _arm_covariates = arm_covariate_values(_filtered, _group_cols, _num_cov_cols)
            _has_moi = any(float(_c.get("moi", 0) or 0) > 0 for _c in _conds.values())
            if _has_moi:
                for _a, _c in _conds.items():
                    if float(_c.get("moi", 0) or 0) > 0:
                        _arm_covariates.setdefault(_a, {})["moi"] = float(_c["moi"])
            _cov_names = list(dict.fromkeys(_num_cov_cols + (["moi"] if _has_moi else [])))

            st.markdown("### 4 · Overlay")
            st.caption(f"{len(_arms)} group(s) · {len(_obs_keys)} observable(s) in the data · "
                       "pick what to overlay against the current model. Each arm is simulated once; "
                       "every observable is projected from that one trajectory.")
            oc1, oc2 = st.columns(2)
            with oc1:
                _sel_arms = st.multiselect("Groups to overlay", _arms, default=_arms[:min(4, len(_arms))], key="fit_arms")
            with oc2:
                _sel_obs = st.multiselect("Observables", _obs_keys, default=_obs_keys,
                                          format_func=lambda k: OBSERVABLES.get(k, {}).get("label", k), key="fit_obs_sel")

            # Per-observable link parameters (OD = biomass / od_to_cfu; luminescence =
            # active biomass × rlu_per_cell). One input per selected observable that needs
            # one; OD instead uses the debris module's get_od() when it is enabled.
            _debris_on = st.session_state.get("int_debris_enabled", False)
            _link_vals = {}
            _link_obs = [k for k in _sel_obs if OBSERVABLES.get(k, {}).get("link")]
            if _link_obs:
                _lcols = st.columns(min(3, len(_link_obs)))
                for _li, _ok in enumerate(_link_obs):
                    _sp = OBSERVABLES[_ok]
                    if _ok == "od":
                        # OD↔CFU is ONE physical factor. Bind this input to the model's
                        # int_od_to_cfu_conversion_factor so the SAME value drives B₀
                        # (OD×factor), the overlay OD (biomass÷factor) and the debris
                        # get_od() — no separate/competing od_to_cfu fields.
                        st.session_state["int_od_to_cfu_conversion_factor"] = \
                            _lcols[_li % len(_lcols)].number_input(
                                "od_to_cfu (CFU per OD unit)",
                                value=float(st.session_state.get(
                                    "int_od_to_cfu_conversion_factor", 2e8) or 2e8),
                                format="%.3e", key="fit_link_od",
                                help="The single OD↔CFU factor: B₀ = first OD × this; overlay "
                                     "OD = biomass ÷ this; the debris module uses it too.")
                        _link_vals["od"] = float(st.session_state["int_od_to_cfu_conversion_factor"])
                    else:
                        _pn, _op, _dflt = _sp["link"]
                        _link_vals[_ok] = _lcols[_li % len(_lcols)].number_input(
                            f"Link · {_sp['label']} ({_pn})", value=float(_dflt), format="%.3e",
                            key=f"fit_link_{_ok}",
                            help="Scales model state → signal. Tunable below / future fit parameter.")
            if "od" in _sel_obs and _debris_on:
                st.caption("OD includes lysed-cell **debris** (`get_od`) — it uses the same "
                           "od_to_cfu set above; the debris rates are in the OD/debris block below.")

            # ── Observation model — which model compartments each signal counts ──
            # obs = f(B, D, I, H, debris). CFU defaults to culturable cells only (B+D):
            # infected (I) / hibernating (H) cells don't form colonies. This SAME set is
            # used for the overlay, the RMSE, and the pbisim-fit residuals (threaded through).
            _obs_comp = {}
            # Only registry observables have a defined compartment mapping. A dataset's
            # observable column can contain arbitrary strings (e.g. "colony_count",
            # "od600") — skip those here (they still overlay via the defensive path).
            _bact_obs = [k for k in _sel_obs if k in OBSERVABLES and k != "pfu"]
            if _bact_obs:
                with st.expander("🧫 Observation model — which compartments each signal reflects"):
                    st.caption("Define what each measurement counts. **CFU = culturable only "
                               "(B + D)** by default; toggle I/H if your assay recovers them. "
                               "OD/turbidity normally includes all cells (+ debris when the OD "
                               "module is on).")
                    _comp_labels = {"B": "B (active)", "D": "D (dormant)",
                                    "I": "I (infected)", "H": "H (hibernating)"}
                    for _ok in _bact_obs:
                        _sp = OBSERVABLES[_ok]
                        _dflt = set(_sp["prefixes"])
                        st.markdown(f"**{_sp['label']}** — obs = sum of:")
                        _ccols = st.columns(len(OBS_COMPARTMENTS))
                        _chosen = [
                            _cp for _ci, _cp in enumerate(OBS_COMPARTMENTS)
                            if _ccols[_ci].checkbox(_comp_labels[_cp], value=(_cp in _dflt),
                                                    key=f"fit_obscomp_{_ok}_{_cp}")
                        ]
                        _obs_comp[_ok] = tuple(_chosen) if _chosen else _sp["prefixes"]
                        if _ok == "od" and _debris_on:
                            st.caption("OD/debris module is on → OD uses the model's `get_od()` "
                                       "(all cells + debris); this choice applies only when debris is off.")

            _t_end_fit = st.number_input("Overlay duration (h)", value=float(np.ceil(_long["time"].max())), step=1.0, key="fit_tend")

            # Per-arm conditions — each group can carry its own growth phase (pre-run),
            # initial density B₀, and dose. Like pbisim-fit, each arm's B₀ defaults to
            # THAT ARM'S FIRST DATA POINT (CFU directly, or OD × od_to_cfu → CFU): the
            # inoculum is anchored to the data, not to the builder. The model's B₀ vector
            # supplies only the strain/genotype *ratio* — its magnitude is renormalised to
            # this B₀ (which is why changing the builder's B₀ magnitude doesn't move the
            # overlay). Falls back to the model total when the first point isn't a bacteria
            # proxy (e.g. PFU-only arm).
            def _model_total_b0():
                if st.session_state.get("int_builder_mode", "").startswith("Binary"):
                    if st.session_state.get("int_brg_use_eq_ic", False):
                        return float(st.session_state.get("int_brg_eq_total_B", 1e7)) or 1e7
                    _sv = st.session_state.get("int_brg_initial_B", {})
                    return float(sum(float(v) for v in _sv.values())) or 1e7
                return float(sum(float(_s.get("initial_B", 0.0))
                                 for _s in st.session_state.get("int_strains", []))) or 1e7
            _fallback_b0 = _model_total_b0()
            # The single OD↔CFU factor (bound to the §4 input above) — B₀ = first OD × this.
            _od2cfu_b0 = float(st.session_state.get("int_od_to_cfu_conversion_factor", 2e8) or 2e8)

            def _arm_first_b0(arm):
                """Arm's earliest bacteria observation as CFU/mL — CFU verbatim, OD×od_to_cfu,
                else the model total (pbisim-fit's data-anchored initial condition)."""
                _dd = _agg[_agg["arm"] == arm]
                if not len(_dd):
                    return _fallback_b0
                _d0 = _dd[_dd["time"] == _dd["time"].min()]
                _c = _d0[_d0["observable"] == "cfu"]
                if len(_c) and float(_c["value"].iloc[0]) > 0:
                    return float(_c["value"].iloc[0])
                _o = _d0[_d0["observable"] == "od"]
                if len(_o) and float(_o["value"].iloc[0]) > 0:
                    return float(_o["value"].iloc[0]) * _od2cfu_b0
                return _fallback_b0
            _arm_cond = {}
            _dose_help = ("Dose seeds the phage inoculum as an absolute PFU/mL titre for that arm."
                          if _dose_unit == "pfu" else
                          "MOI seeds the phage inoculum as MOI × B₀ for that arm.")
            with st.expander(f"Per-arm conditions (B₀ · growth phase · {_dose_lbl})", expanded=False):
                # B₀ source — the SINGLE control for the initial bacterial density (there is
                # deliberately no fit_initial_cfu row in the parameter table that could
                # silently override it). Under pbisim-fit's additive model a known inoculum is
                # a DoseRecord(target="bacteria"), exactly like the phage dose — not the (noisy)
                # first observation. The "Estimate (shared/per-arm)" options ARE how you free B₀
                # (they wire free_initial_conditions); the others fix it.
                _b0_mode = st.radio(
                    "Initial bacterial density B₀",
                    ["First observation", "Shared value", "Per-arm values",
                     "Estimate (shared)", "Estimate (per-arm)"],
                    horizontal=True, key="fit_b0_mode",
                    help="First observation → earliest CFU/OD point (noisy). Shared / Per-arm → "
                         "a known inoculum, recorded as a bacteria dose. Estimate → fit B₀ as an "
                         "additive offset (free_initial_conditions), one value shared across arms "
                         "or one per arm.")
                _estimate_b0 = {"Estimate (shared)": "shared",
                                "Estimate (per-arm)": "per_arm"}.get(_b0_mode, "none")
                if _b0_mode == "First observation":
                    st.warning("B₀ is taken from each arm's **first observation, which includes "
                               "measurement noise**. For a known inoculum choose *Shared*/*Per-arm*, "
                               "or *Estimate* to fit it as an additive offset.")
                elif _estimate_b0 != "none":
                    st.info("B₀ is **estimated** ("
                            + ("one value shared across all arms" if _estimate_b0 == "shared"
                               else "one value per arm")
                            + ") as an additive offset over any bacteria dose, via pbisim-fit's "
                            "`free_initial_conditions`. The per-arm B₀ below is a starting "
                            "placeholder for the overlay; the fitted value replaces it after the fit.")
                _shared_b0 = None
                if _b0_mode == "Shared value":
                    _shared_b0 = st.number_input("Shared B₀ (CFU/mL) — applied to every arm",
                                                 value=_fallback_b0, format="%.2e", key="fit_b0_shared")
                st.caption("Pre-run 0 → log phase (fresh inoculum). Pre-run > 0 → equilibrate "
                           "toward stationary phase before t=0. " + _dose_help)
                if "od" in _sel_obs and "cfu" not in _sel_obs:
                    st.caption(f"↳ Your data is OD, so first-observation B₀ = (first OD) × od_to_cfu "
                               f"= (first OD) × **{_od2cfu_b0:.2e}**. Set od_to_cfu in the **§4 "
                               "Overlay** OD field above — B₀ recomputes live.")

                def _arm_bac_dose(a):
                    """Total bacteria dose (CFU) imported for this arm, else None."""
                    tot, found = 0.0, False
                    for d in _arm_doses.get(a, []):
                        if d["target"] == "bacteria":
                            tot += d["amount"] * (_od2cfu_b0 if d.get("unit") == "od_units" else 1.0)
                            found = True
                    return tot if found else None

                def _arm_phage_dose(a):
                    for d in _arm_doses.get(a, []):
                        if d["target"] == "phage":
                            return d
                    return None

                for _arm in _sel_arms:
                    _cc = st.columns([2, 1, 1, 1])
                    _cc[0].markdown(f"**{_arm}**")
                    _bac, _ph = _arm_bac_dose(_arm), _arm_phage_dose(_arm)
                    # ── B₀ ── data bacteria dose (if any) wins over the manual gap-filler
                    if _bac is not None:
                        _b0v, _is_dose = _bac, False   # emitted from the dataset by build_dataset
                        _cc[1].caption(f"B₀ = {_b0v:.2e}\n(data dose)")
                    elif _b0_mode == "Shared value":
                        _b0v, _is_dose = float(_shared_b0), True
                        _cc[1].caption(f"B₀ = {_b0v:.2e}")
                    elif _b0_mode == "Per-arm values":
                        _b0v = _cc[1].number_input("B₀ (CFU/mL)", value=_arm_first_b0(_arm),
                                                   format="%.2e", key=f"fit_cond_b0_{_arm}")
                        _is_dose = True
                    else:  # First observation / Estimate — COMPUTED live (a caption, not a
                        # keyed widget: a disabled number_input would cache its first value
                        # and never reflect a changed od_to_cfu / mode switch).
                        _b0v, _is_dose = _arm_first_b0(_arm), False
                        _cc[1].caption(f"B₀ = {_b0v:.2e}\n"
                                       + ("(est. start)" if _estimate_b0 != "none" else "(first obs)"))
                    _prv = _cc[2].number_input("Pre-run (h)", value=0.0, step=4.0,
                                               key=f"fit_cond_prerun_{_arm}")
                    # ── phage dose ── data phage dose (if any) wins over the manual MOI/PFU
                    if _ph is not None:
                        _moiv, _moi_unit = float(_ph["amount"]), _ph.get("unit", "pfu")
                        _cc[3].caption(f"{_moiv:.3g} {_moi_unit}\n(data dose)")
                    else:
                        _moiv, _moi_unit = _cc[3].number_input(
                            _dose_lbl, format="%g", value=float(_conds.get(_arm, {}).get("moi", 0.0)),
                            key=f"fit_cond_moi_{_arm}"), _dose_unit
                    _arm_cond[_arm] = {"b0": _b0v, "b0_is_dose": _is_dose, "prerun": _prv,
                                       "moi": _moiv, "moi_unit": _moi_unit}

            # ── 5. Manual parameter tuning (Phase B) ─────────────────────────
            # Edit the model's ACTUAL parameter values via the SAME builder panels as
            # the Interactive Simulator (render_model_builder) — every parameter and
            # every builder mode (Direct / BRG / StrainSet) is available, and the two
            # can never drift. Params that live OUTSIDE the Tab-1 builder (latent
            # stages, nutrient environment, OD/debris) get a compact block here. All
            # edits write the shared int_* state, so they ARE the live model (no
            # separate apply step) and the overlay below reflects them immediately.
            st.markdown("### 5 · Manual parameter tuning")
            st.caption("Edit the model's real parameter values, then re-overlay. These are the "
                       "**same controls as the Interactive Simulator** — all builder modes and "
                       "every parameter — so nothing is missing. Edits update the live model "
                       "directly (no separate apply step) and are savable as a Scenario / Parts.")
            _show_builder = st.toggle(
                "Show model builder & structural parameters",
                value=bool(st.session_state.get("fit_show_builder", False)), key="fit_show_builder")
            if _show_builder:
                # Structural parameters that are NOT part of the Tab-1 model builder
                # (they live in the Simulator's Environment / Solver tabs): latent-stage
                # count, the nutrient environment, and the OD/debris module — including
                # the dormant-cell optical weight `dormant_od_fraction` (default 1.0).
                with st.container(border=True):
                    st.markdown("**Global & structural**")
                    _track_nut = st.session_state.get("int_track_nutrients", True)
                    gk1, gk2 = st.columns([1, 2])
                    with gk1:
                        st.session_state["int_n_latent"] = int(st.number_input(
                            "Latent compartments (L)", min_value=1, max_value=50,
                            value=int(st.session_state.get("int_n_latent", 5)), step=1,
                            key="fit_edit_n_latent",
                            help="Number of phage latent (eclipse) stages — Erlang shape of the latent period."))
                    with gk2:
                        st.caption("Carrying capacity **K** and Monod constant **Ks** are in the "
                                   "builder's *Growth signal* section below (they depend on the growth model).")
                    if _track_nut:
                        nk1, nk2, nk3, nk4 = st.columns(4)
                        with nk1:
                            st.session_state["int_initial_S"] = st.number_input(
                                "Initial nutrient (S₀)", value=float(st.session_state.get("int_initial_S", 1.0)),
                                format="%g", key="fit_edit_S0")
                        with nk2:
                            st.session_state["int_recycle_fraction"] = st.number_input(
                                "Recycle fraction", value=float(st.session_state.get("int_recycle_fraction", 0.0)),
                                format="%g", key="fit_edit_recycle")
                        with nk3:
                            st.session_state["int_s_in"] = st.number_input(
                                "Nutrient inflow (s_in)", value=float(st.session_state.get("int_s_in", 0.0)),
                                format="%g", key="fit_edit_s_in")
                        with nk4:
                            st.session_state["int_s_out"] = st.number_input(
                                "Nutrient washout (s_out)", value=float(st.session_state.get("int_s_out", 0.0)),
                                format="%g", key="fit_edit_s_out")
                        st.caption("Infected-cell nutrient consumption is now **per-phage** — edit it "
                                   "on each phage in the model builder below.")
                    else:
                        st.caption("Nutrient tracking is off (constant/logistic growth) — S₀/recycle/inflow/"
                                   "washout are inactive. Choose a nutrient growth signal in the builder to fit them.")
                    # OD / debris module — enableable HERE (its checkbox lives in the
                    # Simulator's Environment tab, which the embedded builder doesn't
                    # include, so calibration would otherwise be unable to turn on OD
                    # fitting). Required to fit the OD observable via lysed-cell debris.
                    st.session_state["int_debris_enabled"] = st.checkbox(
                        "Track cell debris → model OD (needed to fit the OD observable)",
                        value=bool(st.session_state.get("int_debris_enabled", False)),
                        key="fit_edit_debris_enabled",
                        help="OD = (active B + I + dormant_OD_weight·(D+H) + debris) / od_to_cfu. "
                             "Off → OD is a plain biomass×link scaling.")
                    if st.session_state.get("int_debris_enabled", False):
                        st.markdown("*OD / debris module*")
                        st.caption("od_to_cfu is set once in the **§4 Overlay** OD field above "
                                   "(the same value the debris OD uses) — edit it there.")
                        dk2, dk3, dk4, dk5 = st.columns(4)
                        with dk2:
                            st.session_state["int_debris_u"] = st.number_input(
                                "Debris · deaths (u)", value=float(st.session_state.get("int_debris_u", 0.4)),
                                format="%g", key="fit_edit_debris_u")
                        with dk3:
                            st.session_state["int_debris_v"] = st.number_input(
                                "Debris · lysis (v)", value=float(st.session_state.get("int_debris_v", 0.2)),
                                format="%g", key="fit_edit_debris_v")
                        with dk4:
                            st.session_state["int_debris_kdis"] = st.number_input(
                                "Dissolution (k_dis)", value=float(st.session_state.get("int_debris_kdis", 0.01)),
                                format="%g", key="fit_edit_debris_kdis")
                        with dk5:
                            st.session_state["int_dormant_od_fraction"] = st.number_input(
                                "Dormant OD weight", min_value=0.0, max_value=1.0,
                                value=float(st.session_state.get("int_dormant_od_fraction", 1.0)),
                                format="%g", key="fit_edit_dorm_od",
                                help="Optical weight of a dormant/hibernating cell (D, H) relative to an "
                                     "active cell in OD. 1.0 = counts fully; 0.0 = invisible to OD.")
                    else:
                        st.caption("OD/debris off → the OD observable is a plain biomass × link "
                                   "scaling. Tick the box above to model OD from debris (and expose "
                                   "od_to_cfu, debris rates, and the dormant-cell OD weight).")

                # The full model builder — mode selector + growth/death signals + the
                # selected Direct / BRG / StrainSet body. SAME function the Simulator uses,
                # so every parameter (dormancy kinetics in every mode, pseudolysogeny,
                # signals, PK, …) is present and always in sync. Rendered OUTSIDE any
                # st.expander (the builder has its own expanders; nesting is illegal).
                from pbisim_app.views.simulator import render_model_builder
                st.markdown("**Model builder**")
                render_model_builder(inoculum_mode="ratio")
                st.caption("Tip: B₀ may be overridden by an equilibrium/pre-run initial condition in some "
                           "builder modes; the phage inoculum in the overlay comes from each group's MOI × B₀.")

            # Compute the overlay only when the button is clicked; store the plot
            # data in session_state so the visualization stays alive across page
            # navigation (and reruns) until it is explicitly re-run or the dataset
            # is cleared.
            # Context shared by the manual overlay and the post-fit fitted overlay.
            _ovl_ctx = {
                "sel_arms": _sel_arms, "sel_obs": _sel_obs, "arm_cond": _arm_cond,
                "conds": _conds, "agg": _agg, "link_vals": _link_vals,
                "obs_compartments": _obs_comp, "arm_covariates": _arm_covariates,
                "debris_on": _debris_on, "band": _band, "band_choice": _band_choice,
                "stat_key": _stat_key, "stat": _stat, "t_end": _t_end_fit,
                "dose_unit": _dose_unit, "dose_label": _dose_lbl,
                "method": st.session_state.get("int_solver_method", "BDF"),
                "thr": st.session_state.get("int_extinction_threshold", 1.0) or None,
            }
            if st.button("Overlay model on data", key="fit_overlay", width="stretch", type="primary"):
                try:
                    _config, _iB, _iP, _iS, _mk = build_nominal_config_from_gui()
                    st.session_state["calib_overlay_result"] = _compute_overlay(
                        _config, _iB, _iP, _iS, _mk, _ovl_ctx)
                except Exception as e:
                    st.session_state["calib_overlay_result"] = None
                    st.error(f"Overlay failed: {e}")
                    import traceback
                    st.code(traceback.format_exc())

            # Render the (persisted) overlay result if one exists.
            _ovr = st.session_state.get("calib_overlay_result")
            if _ovr and _ovr.get("panels"):
                import plotly.graph_objects as go
                _palette = ["#0d7a68", "#c1873a", "#5457a6", "#b5487f", "#3b6fb5",
                            "#2e8b57", "#a0522d", "#6a5acd"]

                def _rgba(hexc, a):
                    h = hexc.lstrip("#")
                    return f"rgba({int(h[0:2], 16)},{int(h[2:4], 16)},{int(h[4:6], 16)},{a})"

                st.markdown(f"#### {_ovr['title']}")
                _npan = len(_ovr["panels"])
                # Combined objective (the single number to minimise; matches what
                # pbisim-fit will jointly minimise) + one R²/RMSE chip per observable.
                _cols = st.columns(1 + _npan)
                _cval = _ovr.get("combined", float("nan"))
                _cols[0].markdown(
                    f"""<div class="metric-container">
                        <div class="metric-label">Combined objective J</div>
                        <div class="metric-value">{_cval:.3f}</div>
                        <div class="metric-sub">pooled log₁₀ least-squares — the pbisim-fit NLS objective</div>
                    </div>""", unsafe_allow_html=True)
                for _ci, _pan in enumerate(_ovr["panels"]):
                    _cols[_ci + 1].markdown(
                        f"""<div class="metric-container">
                            <div class="metric-label">{_pan['label']} · R²</div>
                            <div class="metric-value">{_pan['r2']:.3f}</div>
                            <div class="metric-sub">RMSE (log₁₀) {_pan['rmse']:.3f} · n={_pan['n']}</div>
                        </div>""", unsafe_allow_html=True)
                st.markdown("<br>", unsafe_allow_html=True)

                # One panel per observable (small multiples — units differ).
                for _pan in _ovr["panels"]:
                    _pfig = go.Figure()
                    for _i, _s in enumerate(_pan["series"]):
                        _color = _palette[_i % len(_palette)]
                        _yp = np.maximum(_s["pred"], 1e-30) if _pan["log"] else _s["pred"]
                        if (not _s["is_raw"]) and _s["obs_lo"] is not None:
                            _pfig.add_trace(go.Scatter(x=_s["obs_time"], y=_s["obs_hi"], mode="lines",
                                                       line=dict(width=0), showlegend=False, hoverinfo="skip"))
                            _pfig.add_trace(go.Scatter(x=_s["obs_time"], y=_s["obs_lo"], mode="lines",
                                                       line=dict(width=0), fill="tonexty",
                                                       fillcolor=_rgba(_color, 0.15), showlegend=False, hoverinfo="skip"))
                        _pfig.add_trace(go.Scatter(x=_s["time"], y=_yp, mode="lines",
                                                   name=f"{_s['label']} (model)", line=dict(color=_color, width=2)))
                        _mk = dict(color=_color, size=5 if _s["is_raw"] else 7,
                                   opacity=0.45 if _s["is_raw"] else 1.0)
                        if not _s["is_raw"]:
                            _mk["line"] = dict(color="#16211f", width=0.4)
                        _pfig.add_trace(go.Scatter(x=_s["obs_time"], y=_s["obs_value"], mode="markers",
                                                   name=f"{_s['label']} (obs)", marker=_mk))
                    _pfig.update_layout(title=_pan["label"], xaxis_title="Time (h)", yaxis_title=_pan["label"],
                                        template="plotly_white", height=380, margin=dict(t=44, b=40),
                                        legend=dict(orientation="h", yanchor="bottom", y=-0.32, x=0))
                    apply_axis_plotly(_pfig, plot_axis_controls(
                        f"calib_ovl_{_pan['obs']}", default_y="Log" if _pan["log"] else "Linear"))
                    st.plotly_chart(_pfig, width="stretch")

                st.markdown(f"#### Fit quality per group × observable (vs {_ovr['stat_label']})")
                st.dataframe(pd.DataFrame(_ovr["metrics"]), width="stretch", hide_index=True)
                st.caption("Edit the parameter values above and re-overlay to improve the fit "
                           "(minimise the combined objective). Edits update the live model directly "
                           "and can be saved as Parts in the Library.")

                # ── Model comparison (AIC / BIC) ─────────────────────────────
                # RSS alone always favours the bigger model; AIC/BIC penalise parameters.
                # Snapshot the current overlay (its pooled log10 residuals + a free-param
                # count) under a name, then rank candidates via pbisim-fit's compare_models.
                _resid = _ovr.get("residuals") or []
                with st.expander("Compare models (AIC / BIC)", expanded=False):
                    st.caption("Snapshot each candidate model (a growth form, a resistance "
                               "mechanism, …) after overlaying it, then compare. AIC/BIC "
                               "penalise extra parameters so a richer model must *earn* its "
                               "complexity. **Only valid across models overlaid on the same "
                               "data** (same points).")
                    _cmp = st.session_state.setdefault("fit_model_comparison", [])
                    _cc1, _cc2, _cc3 = st.columns([2, 1, 1])
                    with _cc1:
                        _cmp_name = st.text_input("Candidate name", key="fit_cmp_name",
                                                  placeholder=f"model {len(_cmp) + 1}")
                    with _cc2:
                        _cmp_k = int(st.number_input("Free params (k)", min_value=1, step=1,
                                                     value=int(st.session_state.get("fit_cmp_k", 3)),
                                                     key="fit_cmp_k",
                                                     help="Number of estimated parameters for THIS "
                                                          "candidate — the AIC/BIC complexity penalty."))
                    with _cc3:
                        st.markdown("<div style='height:1.8em'></div>", unsafe_allow_html=True)
                        _add = st.button("Add current model", key="fit_cmp_add", width="stretch",
                                         disabled=not _resid)
                    if _add and _resid:
                        _cmp.append({"name": _cmp_name.strip() or f"model {len(_cmp) + 1}",
                                     "residuals": list(_resid), "k": _cmp_k,
                                     "n": int(_ovr.get("n_resid", len(_resid)))})
                        st.session_state["fit_model_comparison"] = _cmp
                        st.rerun()
                    if not _resid:
                        st.info("Overlay a model first — a snapshot needs residuals.")
                    if _cmp:
                        _ns = {c["n"] for c in _cmp}
                        if len(_ns) > 1:
                            st.warning(f"Candidates were overlaid on different numbers of data "
                                       f"points ({sorted(_ns)}) — AIC/BIC are only comparable on "
                                       "identical data. Re-snapshot with the same group/observable "
                                       "selection.")
                        _ranked = compare_fit_models(
                            {c["name"]: (c["residuals"], c["k"]) for c in _cmp})
                        _rows = [{"model": r["name"], "k": r["n_params"], "n": r["n_data"],
                                  "RSS (log₁₀)": r["rss"], "AIC": r["aic"], "BIC": r["bic"],
                                  "AICc": r["aicc"], "ΔAIC": r["delta_aic"], "ΔBIC": r["delta_bic"]}
                                 for r in _ranked]
                        st.dataframe(
                            pd.DataFrame(_rows).style.format({
                                "RSS (log₁₀)": "{:.3f}", "AIC": "{:.1f}", "BIC": "{:.1f}",
                                "AICc": "{:.1f}", "ΔAIC": "{:.2f}", "ΔBIC": "{:.2f}"}),
                            width="stretch", hide_index=True)
                        st.caption("Lower is better; ΔAIC = 0 is the best-supported model. "
                                   "ΔAIC < 2 ≈ comparable support; > 10 ≈ decisively worse.")
                        if st.button("Clear comparison", key="fit_cmp_clear"):
                            st.session_state.pop("fit_model_comparison", None)
                            st.rerun()

            # ── 5b · Export fit specification (pbisim-fit hand-off) ───────────
            st.markdown("### 5b · Export fit specification")
            st.caption("Bundle the data (long format), per-arm conditions (growth phase / B₀ / MOI), "
                       "the selected observables + detection floors, and the current parameter values "
                       "into a **pbisim-fit specification** — exactly the inputs its NLS `refine_nls` / "
                       "`fit()` consume (dataset → ExperimentalDataset, nls_cfg → NLSConfig).")
            try:
                _spec_cfg, _sB, _sP, _sS, _smk = build_nominal_config_from_gui()
                _od_link = _link_vals.get("od")
                if _od_link is None and _debris_on:
                    _od_link = st.session_state.get("int_od_to_cfu_conversion_factor")
                _spec = build_fit_spec(_agg, _sel_arms, _sel_obs, _arm_cond,
                                       od_to_cfu=_od_link, dose_unit=_dose_unit,
                                       model_params=config_param_snapshot(_spec_cfg))
                st.download_button(
                    "Download fit spec (JSON)", data=json.dumps(_spec, indent=2),
                    file_name="pbisim_fit_spec.json", mime="application/json",
                    key="fit_export_spec", width="stretch")
            except Exception as _e:
                st.caption(f"(Select at least one group and observable first — {_e})")

            # ── 5c · Run NLS fit (pbisim-fit, in-app) ─────────────────────────
            st.markdown("### 5c · Run NLS fit (pbisim-fit)")
            st.caption("Fit the selected parameters to the selected arms/observables with "
                       "pbisim-fit's non-linear least squares (`refine_nls`) — the same engine "
                       "the export spec feeds. Runs locally; needs `pbisim-fit` installed. "
                       "Fits use the **BDF** stiff solver with Jacobian variable-scaling "
                       "(`x_scale='jac'`) — pbisim-fit's defaults; LSODA can silently stall on "
                       "these stiff phage–bacteria / sequential-growth ODEs.")
            try:
                from pbisim_app import nls_fit as _nls
                # Fit runs against a CHOSEN Model, not the live builder — a frozen
                # saved/demo model can't be contaminated by unrelated builder edits.
                _mopts = model_options()
                _mdef = st.session_state.active_model if st.session_state.active_model in _mopts else WORKING_DRAFT_LABEL
                _fit_model = st.selectbox(
                    "Model to fit", _mopts, index=_mopts.index(_mdef), key="fit_model_sel",
                    help="The base model whose free parameters are estimated. Saved/demo "
                         "models are frozen snapshots; 'Working draft (live)' uses the "
                         "current builder state.")
                _fit_snap = resolve_model_snapshot(_fit_model)
                _fit_cfg, _fB, _fP, _fS, _fmk = build_config_from_model(_fit_snap)
                # Builder mode of the CHOSEN model (frozen snapshot, or the live draft) —
                # makes the parameter catalog model-aware (BRG → fitness cost, not raw
                # resistant growth; Direct/StrainSet → raw growth, no selection virtuals).
                _fit_mode = ((_fit_snap or {}).get("int_builder_mode")
                             if _fit_snap is not None
                             else st.session_state.get("int_builder_mode", "Direct (ModelBuilder)"))
                # ── Comprehensive parameter table: fix / free each model parameter ─
                _targets_cat = _nls.available_targets(
                    _fit_cfg, initial_cfu=float(np.sum(_fB)) if _fB is not None else None,
                    initial_pfu=(float(_fP[0]) if _fP is not None and len(_fP) else None),
                    builder_mode=_fit_mode)
                _path_label = {p: lab for (lab, p, *_r) in _targets_cat}
                _label_to_path = {lab: p for (lab, p, *_r) in _targets_cat}
                def _bound(s, log, is_lower):
                    """Blank bound → unconstrained on that side: +inf above; -inf (or 0
                    for a log-space parameter) below. A set value (scientific notation
                    accepted) is used verbatim."""
                    v = _pf(s)
                    if v is not None:
                        return float(v)
                    if not is_lower:
                        return float("inf")
                    return 0.0 if log else float("-inf")

                import hashlib as _hl
                # A revision counter (bumped when a text spec is applied) is appended to
                # the editor KEYS only — it forces the keyed data_editors to reload from
                # the freshly-parsed dataframes WITHOUT triggering the defaults rebuild
                # (which keys off _sig).
                _rev = int(st.session_state.get("fit_spec_rev", 0))
                _sig = _hl.md5(("|".join(p for (_l, p, *_r) in _targets_cat) + "|" + _fit_model)
                               .encode()).hexdigest()[:10]
                _reset_tbl = st.button("Reset table to model defaults", key="fit_tbl_reset")
                # Default table on model change / reset: one row per model
                # parameter, role=Fixed, bounds blank. Cells are TEXT (sci notation).
                if st.session_state.get("fit_targets_sig") != _sig or _reset_tbl:
                    _rows = [{"parameter": lab, "path": p, "role": "Fixed", "value": f"{v:g}",
                              "lower": "", "upper": "", "log": bool(log),
                              "prior \u03bc": "", "prior \u03c3": "", "expression": ""}
                             for (lab, p, v, lo_f, hi_f, log) in _targets_cat]
                    st.session_state["fit_targets_df"] = pd.DataFrame(_rows)
                    st.session_state["fit_targets_sig"] = _sig

                # ── Curve stripping — analytic initial estimates (warm-start) ─────
                with st.expander("\U0001F3AF Curve stripping \u2014 analytic initial estimates", expanded=False):
                    st.caption(
                        "Read growth rate, carrying capacity, adsorption, resistant fraction \u2026 "
                        "straight off the **OD** curves (pbisim-fit's `propose_initials` \u2014 no ODE "
                        "fit), then **seed the fit table** with them. A good warm start makes the "
                        "NLS fit converge far more reliably on these stiff phage\u2013bacteria ODEs.")
                    if "_strip_apply_msg" in st.session_state:
                        st.success(st.session_state.pop("_strip_apply_msg"))
                    _od_arms = (_agg[_agg["observable"] == "od"]["arm"].unique().tolist()
                                if "od" in set(_agg["observable"].unique()) else [])

                    def _armmoi(_a):
                        _v = (_conds.get(_a) or {}).get("moi", None)
                        try:
                            return float(_v)
                        except (TypeError, ValueError):
                            return float("nan")
                    _mois = {a: _armmoi(a) for a in _od_arms}
                    _auto_ctrl = next((a for a in _od_arms if _mois[a] <= 0.0), None)   # NaN<=0 False
                    _treat_arms = [a for a in _od_arms if np.isfinite(_mois[a]) and _mois[a] > 0.0]
                    if not _od_arms:
                        st.info("Curve stripping needs an **OD/turbidity** observable \u2014 none in "
                                "the selected data.")
                    elif not _treat_arms:
                        st.info("Curve stripping needs **\u22651 phage arm (MOI > 0)** \u2014 none found. "
                                "Set per-arm MOI in \u00a74.")
                    else:
                        # Control arm: auto-detected MOI=0 by default; overridable for datasets
                        # whose control has a blank/NaN MOI (auto-detect can't find those).
                        if _auto_ctrl is None:
                            st.warning("No MOI = 0 control auto-detected (the control's MOI may be "
                                       "blank/NaN). Pick the control arm below.")
                        _ctrl_choice = st.selectbox(
                            "Control arm (no phage)", _od_arms,
                            index=(_od_arms.index(_auto_ctrl) if _auto_ctrl in _od_arms else 0),
                            key="strip_ctrl",
                            help="The no-phage reference. Auto-set to the MOI = 0 arm; override if "
                                 "your control's MOI is blank/NaN.")
                        try:
                            _b_def = float(np.asarray(getattr(_fit_cfg, "burst_sizes", [[50.0]]))[0, 0])
                        except Exception:
                            _b_def = 50.0
                        _b_fixed = st.number_input(
                            "Fixed burst size (b)", min_value=1.0,
                            value=float(_b_def if _b_def and _b_def > 0 else 50.0),
                            key="strip_bfixed",
                            help="Curve stripping treats burst size as known while it solves "
                                 "adsorption from the lysis curve.")
                        with st.expander("Advanced picker guards (real-data artifacts)", expanded=False):
                            _req_peak = st.checkbox(
                                "Require growth before the lysis peak", value=False,
                                key="strip_reqpeak",
                                help="The lysis peak must clear the inoculum OD, skipping a lag-phase "
                                     "dip at t\u22480 \u2014 avoids an inflated adsorption k on screens that dip "
                                     "before growing.")
                            _min_vir = st.number_input(
                                "Min. virulence to keep a phage arm (0 = off)", min_value=0.0,
                                value=0.0, key="strip_minvir", format="%g",
                                help="A phage arm whose max local virulence vs the control is below "
                                     "this abstains (no k / lysis peak) \u2014 filters resistant / no-kill "
                                     "arms. 0 disables.")
                        _fit_f0 = st.checkbox(
                            "Refine f\u2080 (1-D fit)", value=False, key="strip_fitf0",
                            help="After stripping, refine the pre-existing resistant fraction f\u2080 "
                                 "with pbisim-fit's single-parameter forward-simulation fit "
                                 "(`refine_f0`) \u2014 the one fragile analytic estimate. Slower: a "
                                 "few ODE solves.")
                        _od2cfu = float(st.session_state.get("int_od_to_cfu_conversion_factor", 2e8) or 2e8)
                        _mk = float(getattr(_fit_cfg, "monod_constant", 0.3) or 0.3)

                        # Growth signal + debris of the CHOSEN model (frozen snapshot, or the live
                        # draft) \u2192 so the g-correction and f\u2080 refine use the right forward model, and
                        # f\u2080 refine uses the model's debris (not propose_initials' frozen defaults).
                        def _model_val(_k, _d):
                            return ((_fit_snap or {}).get(_k, _d) if _fit_snap is not None
                                    else st.session_state.get(_k, _d))
                        _gs_fn = _model_val("int_growth_function", "monod_growth") or "monod_growth"
                        _strip_gs = _nls.strip_growth_signal(_gs_fn)
                        _strip_debris = None
                        if _model_val("int_debris_enabled", False):
                            _strip_debris = {
                                "v": float(_model_val("int_debris_v", 0.2) or 0.2),
                                "kdis": float(_model_val("int_debris_kdis", 0.01) or 0.01),
                                "u": float(_model_val("int_debris_u", 0.4) or 0.4)}
                        # f₀-refine forward-model nuisances — editable inputs, defaulting to
                        # the chosen model's values (burst = the Fixed-burst input above). Shown
                        # only when Refine f₀ is on (they affect ONLY the f₀ refine).
                        try:
                            _lat_def = float(np.asarray(getattr(_fit_cfg, "latent_periods", [[0.5]]))[0, 0])
                        except Exception:
                            _lat_def = 0.5
                        _dbg_def = _strip_debris or {"v": 0.3, "kdis": 0.1, "u": 0.0}
                        if _fit_f0:
                            st.caption("f₀-refine forward-model nuisances (default to the chosen "
                                       "model; burst = the fixed-burst input above):")
                            _lat = st.number_input(
                                "f₀ · latent period (h)", min_value=0.0, value=float(_lat_def),
                                key="strip_f0_latent")
                            _dcols = st.columns(3)
                            _strip_debris = {
                                "v": _dcols[0].number_input(
                                    "f₀ · debris v", min_value=0.0,
                                    value=float(_dbg_def["v"]), key="strip_f0_dv"),
                                "kdis": _dcols[1].number_input(
                                    "f₀ · debris k_dis", min_value=0.0,
                                    value=float(_dbg_def["kdis"]), key="strip_f0_kdis"),
                                "u": _dcols[2].number_input(
                                    "f₀ · debris u", min_value=0.0,
                                    value=float(_dbg_def["u"]), key="strip_f0_u"),
                            }
                        else:
                            _lat = _lat_def
                        st.caption(f"Using od_to_cfu = {_od2cfu:.3e} (§4) · monod_constant = {_mk:g} · "
                                   f"growth model = {_strip_gs}. B₀ defaults to the control's first OD.")
                        if _fit_f0 and _gs_fn not in ("monod_growth", "logistic_growth",
                                                      "constant_growth"):
                            st.caption(f"Note: '{_gs_fn}' has no dedicated stripping forward model \u2014 "
                                       f"approximated as '{_strip_gs}' for the g-correction and f\u2080 refine.")
                        if st.button("Compute initial estimates", key="strip_compute"):
                            try:
                                st.session_state["calib_strip_result"] = _nls.strip_curves(
                                    _agg, _conds, b_fixed=_b_fixed, od_to_cfu=_od2cfu,
                                    monod_constant=_mk, control_arm=_ctrl_choice, fit_f0=_fit_f0,
                                    growth_signal=_strip_gs, debris=_strip_debris, f0_latent=_lat,
                                    require_growth_peak=_req_peak, min_virulence=_min_vir)
                            except Exception as _se:
                                st.session_state.pop("calib_strip_result", None)
                                st.error(f"Curve stripping failed: {_se}")
                        _R = st.session_state.get("calib_strip_result")
                        if _R is not None:
                            _erows = [{"estimate": _nm,
                                       "value": (f"{_val:.4g}" if np.isfinite(_val) else "\u2014"),
                                       "reliability": _rel}
                                      for _nm, (_val, _rel) in _R.estimates.items()]
                            st.dataframe(pd.DataFrame(_erows), hide_index=True, width="stretch")
                            with st.expander("Full report", expanded=False):
                                st.code(_R.report())
                            _inits = dict(_R.initials)
                            if st.button("Seed fit table with these estimates", key="strip_seed",
                                         width="stretch"):
                                _df = st.session_state.get("fit_targets_df")
                                # g_max maps to ALL bacterial species' growth_rates[i] rows.
                                _seed_inits = _broadcast_growth(
                                    _inits, getattr(_fit_cfg, "n_bacteria", 1))
                                if _df is not None and _seed_inits:
                                    _paths = set(_df["path"])
                                    _applied = [p for p in _seed_inits if p in _paths]
                                    _skipped = [p for p in _seed_inits if p not in _paths]
                                    for _p in _applied:
                                        _df.loc[_df["path"] == _p, "value"] = f"{_seed_inits[_p]:g}"
                                        _df.loc[_df["path"] == _p, "role"] = "Free"
                                    st.session_state["fit_targets_df"] = _df
                                    st.session_state["fit_spec_rev"] = _rev + 1
                                    if _applied:
                                        st.success("Seeded (set Free): " + ", ".join(_applied))
                                    if _skipped:
                                        st.caption("Not in this model's fit table (skipped): "
                                                   + ", ".join(_skipped))
                                    st.rerun()
                                elif not _inits:
                                    st.warning("No solid/medium estimates to seed.")

                            # Push the estimates into the Working-draft builder (the manual-
                            # calibration values) AND refresh the overlay in one click.
                            if st.button("Apply estimates to model & refresh overlay",
                                         key="strip_apply", width="stretch"):
                                if not _inits:
                                    st.warning("No solid/medium estimates to apply.")
                                else:
                                    _ap, _sk = _apply_initials_to_builder(_inits)
                                    try:
                                        _c1, _b1, _p1, _s1, _m1 = build_nominal_config_from_gui()
                                        st.session_state["calib_overlay_result"] = _compute_overlay(
                                            _c1, _b1, _p1, _s1, _m1, _ovl_ctx)
                                    except Exception:
                                        st.session_state["calib_overlay_result"] = None
                                    # Persist the outcome across the rerun (a bare st.success here
                                    # is discarded by st.rerun) → shown at the panel top next render.
                                    _msg = ("**Applied to the model & refreshed the overlay.** "
                                            + ("Set: " + ", ".join(_ap) + ". " if _ap else "")
                                            + ("Skipped (no builder slot — e.g. BRG resistant "
                                               "fraction): " + ", ".join(_sk) + ". " if _sk else ""))
                                    if st.session_state.get("fit_show_builder"):
                                        _msg += ("\n\n_The model-builder inputs update in state but a "
                                                 "Streamlit repaint quirk can leave an **open** builder "
                                                 "showing the old numbers — collapse & reopen that "
                                                 "section to see the refreshed values._")
                                    st.session_state["_strip_apply_msg"] = _msg
                                    st.rerun()

                # ── Amortized fit (ONNX) — instant calibrated estimate ────────
                _am_avail = _nls.amortized_available()
                _am_nets = _nls.list_amortized_nets() if _am_avail else []
                _am_od = "od" in set(_agg["observable"].unique())
                with st.expander("⚡ Amortized fit (ONNX) — instant calibrated estimate", expanded=False):
                    if "_amort_apply_msg" in st.session_state:
                        st.success(st.session_state.pop("_amort_apply_msg"))
                    st.caption("An instant, torch-free amortized estimate for an OD assay — a trained "
                               "neural net, complementary to the slow NLS fit and the analytic strip. "
                               "Pick the model trained for your phage/host; the estimate + 95% CI are "
                               "immediate. Seed the NLS fit with it, or apply it to the model.")
                    if not _am_avail:
                        st.info("Needs `onnxruntime` (install `pbisim-fit[onnx]`) — not available here.")
                    elif not (_am_od and _am_nets):
                        st.info("Needs an **OD** observable and an available amortized model.")
                    else:
                        try:
                            _am_bdef = float(np.asarray(getattr(_fit_cfg, "burst_sizes", [[50.0]]))[0, 0])
                        except Exception:
                            _am_bdef = 50.0
                        try:
                            _am_latdef = float(np.asarray(getattr(_fit_cfg, "latent_periods", [[0.5]]))[0, 0])
                        except Exception:
                            _am_latdef = 0.5
                        _am_net = st.selectbox(
                            "Amortized model", _am_nets, key="amort_net",
                            help="Each net is trained for a specific phage/host + OD assay; it is only "
                                 "valid inside that trained envelope (single-phage, OD-only here).")
                        _amc = st.columns(2)
                        _am_burst = _amc[0].number_input(
                            "Burst size (context)", min_value=1.0,
                            value=float(_am_bdef if _am_bdef and _am_bdef > 0 else 50.0),
                            key="amort_burst", help="Supplied context (measured/assumed), not inferred.")
                        _am_latent = _amc[1].number_input(
                            "Latent period (h, context)", min_value=0.0, value=float(_am_latdef),
                            key="amort_latent")
                        st.caption("B₀ and MOI are read from the arms automatically.")
                        if st.button("Run amortized fit", key="amort_run"):
                            try:
                                _am_ds = _nls.build_dataset(_agg, _sel_arms, ["od"], _arm_cond,
                                                            od_to_cfu=_od_link, dose_unit=_dose_unit)
                                st.session_state["calib_amortized_result"] = _nls.amortized_fit(
                                    _am_ds, _am_net, burst=_am_burst, latent=_am_latent)
                            except Exception as _ae:
                                st.session_state.pop("calib_amortized_result", None)
                                st.error(f"Amortized fit failed: {_ae}")
                        _amr = st.session_state.get("calib_amortized_result")
                        if _amr is not None:
                            _arows = []
                            for _k, _v in _amr["map"].items():
                                _lo, _hi = _amr["ci"].get(_k, [None, None])
                                _arows.append({"parameter": _k, "MAP": f"{_v:.4g}",
                                               "95% CI low": (f"{_lo:.4g}" if _lo is not None else "—"),
                                               "95% CI high": (f"{_hi:.4g}" if _hi is not None else "—")})
                            st.dataframe(pd.DataFrame(_arows), hide_index=True, width="stretch")
                            # OOD guard: cross-check against the analytic strip when available.
                            _sr = st.session_state.get("calib_strip_result")
                            if _sr is not None:
                                try:
                                    _rep = _nls.strip_sanity_check(_sr, _amr["initials"])
                                    if _rep.rows and not _rep.ok:
                                        st.warning("Possible out-of-envelope: the amortized estimate "
                                                   "disagrees with the analytic strip on — "
                                                   + ", ".join(_rep.flagged) + " — consider the NLS fit.")
                                    elif _rep.rows:
                                        st.caption("Cross-check: agrees with the analytic strip.")
                                except Exception:
                                    pass
                            else:
                                st.caption("Tip: run **curve stripping** too — its sanity check "
                                           "cross-checks the amortized estimate for out-of-distribution assays.")
                            _amcols = st.columns(2)
                            if _amcols[0].button("Seed fit table (warm-start NLS)", key="amort_seed"):
                                _df = st.session_state.get("fit_targets_df")
                                _si = _broadcast_growth(_amr["initials"], getattr(_fit_cfg, "n_bacteria", 1))
                                if _df is not None and _si:
                                    _paths = set(_df["path"])
                                    _appl = [p for p in _si if p in _paths]
                                    for _p in _appl:
                                        _df.loc[_df["path"] == _p, "value"] = f"{_si[_p]:g}"
                                        _df.loc[_df["path"] == _p, "role"] = "Free"
                                    st.session_state["fit_targets_df"] = _df
                                    st.session_state["fit_spec_rev"] = _rev + 1
                                    if _appl:
                                        st.success("Seeded (set Free): " + ", ".join(_appl))
                                    st.rerun()
                            if _amcols[1].button("Apply to model & refresh overlay", key="amort_apply"):
                                _struct = (_config_structure(_amr.get("config"),
                                                             {"inc": _amr["map"].get("inc")})
                                           if _amr.get("config") is not None else None)
                                _ap, _sk = _apply_initials_to_builder(_amr["initials"],
                                                                      structure=_struct)
                                try:
                                    _c1, _b1, _p1, _s1, _m1 = build_nominal_config_from_gui()
                                    st.session_state["calib_overlay_result"] = _compute_overlay(
                                        _c1, _b1, _p1, _s1, _m1, _ovl_ctx)
                                except Exception:
                                    st.session_state["calib_overlay_result"] = None
                                _m = ("**Applied the amortized fit & refreshed the overlay.** "
                                      + ("Set: " + ", ".join(_ap) + ". " if _ap else ""))
                                if st.session_state.get("fit_show_builder"):
                                    _m += ("\n\n_Collapse & reopen an open model builder to see the "
                                           "refreshed values (Streamlit repaint quirk)._")
                                st.session_state["_amort_apply_msg"] = _m
                                st.rerun()

                st.markdown("**Fit parameters** \u2014 set each parameter's **role**: "
                            "**Fixed** (held at value) \u00b7 **Free** (estimated with bounds/prior) "
                            "\u00b7 **Derived** (`= expression` of a \u03b8, defined below).")
                st.caption("Bounds/priors optional (blank = unconstrained / no prior); values accept "
                           "scientific notation (`1e7`); **log** = log-space search. A **Derived** row "
                           "draws its value/bounds/prior from its \u03b8 (those cells blank out \u2014 set them "
                           "on the \u03b8 below). **Sharing** = set several rows to *Derived* with the same "
                           "\u03b8 (use the **Share** helper). It's all in this one table \u2014 no separate sections.")
                _mode_note = (
                    "This table shows the **" + str(_fit_mode) + "** model's own parameters, so each "
                    "quantity has exactly one control. "
                    + ("**BRG:** the resistant genotype's growth is the **Fitness cost** "
                       "(growth\u00b7(1\u2212cost)), so its raw growth rate isn't listed separately. "
                       if _fit_mode == "Binary Genotypes (BRG)"
                       else "Strains are independent here, so each growth rate is listed and the "
                            "BRG-only *fitness cost* / *resistant fraction* aren't. ")
                    + "**Initial bacterial density (B\u2080)** is not in this table \u2014 it's set by the "
                      "**B\u2080 source** control in \u00a74 (its *Estimate* option is how you free B\u2080).")
                st.caption(_mode_note)
                _tg_ed = st.data_editor(
                    st.session_state["fit_targets_df"], key=f"fit_targets_editor_{_sig}_{_rev}",
                    hide_index=True, width="stretch",
                    column_config={
                        "parameter": st.column_config.TextColumn("parameter", disabled=True),
                        "path": None,
                        "role": st.column_config.SelectboxColumn(
                            "role", options=["Fixed", "Free", "Derived"], required=True,
                            help="Fixed: held at value \u00b7 Free: estimated 1:1 \u00b7 Derived: = expression of \u03b8"),
                        "value": st.column_config.TextColumn("value / start",
                                                             help="Fixed value, or start for a Free fit. e.g. 1e8"),
                        "lower": st.column_config.TextColumn("lower", help="Free only. Blank = unbounded. e.g. 1e7"),
                        "upper": st.column_config.TextColumn("upper", help="Free only. Blank = unbounded. e.g. 1e10"),
                        "log": st.column_config.CheckboxColumn("log", help="Free/\u03b8: log-space search."),
                        "prior \u03bc": st.column_config.TextColumn("prior \u03bc", help="Free only. Gaussian prior mean -> MAP."),
                        "prior \u03c3": st.column_config.TextColumn("prior \u03c3", help="Free only. Prior stdev."),
                        "expression": st.column_config.TextColumn("= expression (Derived)",
                                                                  help="For Derived rows: an expression of \u03b8, e.g. g*(1-cost)"),
                    })

                # A Derived row inherits value/bounds/priors from its \u03b8, so blank those
                # cells (they read as inactive) and keep only role + expression live.
                _recs = _tg_ed.to_dict("records")
                _norm = False
                for _r in _recs:
                    if str(_r.get("role")) == "Derived":
                        for _c in ("value", "lower", "upper", "prior \u03bc", "prior \u03c3"):
                            if str(_r.get(_c, "")).strip():
                                _r[_c] = ""; _norm = True
                        if bool(_r.get("log")):
                            _r["log"] = False; _norm = True
                if _norm:
                    st.session_state["fit_targets_df"] = pd.DataFrame(_recs)
                    st.session_state["fit_spec_rev"] = _rev + 1
                    st.rerun()

                # -- Share helper: tie selected rows to one new theta --
                if st.session_state.pop("_share_clear", None):
                    for _k in ("fit_share_pick", "fit_share_name"):
                        st.session_state.pop(_k, None)
                with st.expander("Share \u2014 tie several parameters to one \u03b8", expanded=False):
                    st.caption("Sets the selected parameters to **Derived = one new \u03b8**; define that "
                               "\u03b8's bounds/prior in the \u03b8 table below.")
                    _sh_pick = st.multiselect("Parameters to share",
                                              [l for (l, *_r) in _targets_cat], key="fit_share_pick")
                    _shc = st.columns([3, 1])
                    _sh_name = _shc[0].text_input("\u03b8 name (blank = auto)", key="fit_share_name")
                    if _shc[1].button("Share \u2192", key="fit_share_go", width="stretch"):
                        if len(_sh_pick) < 2:
                            st.warning("Pick at least two parameters to share.")
                        else:
                            _exist = (set(st.session_state["fit_thetas_df"]["name"])
                                      if "fit_thetas_df" in st.session_state else set())
                            _nm = (_sh_name.strip() or next(
                                f"shared{i}" for i in range(1, 999) if f"shared{i}" not in _exist))
                            _picked = {_label_to_path[l] for l in _sh_pick}
                            _df = st.session_state["fit_targets_df"].copy()
                            for _i, _r in _df.iterrows():
                                if _r["path"] in _picked:
                                    _df.at[_i, "role"] = "Derived"
                                    _df.at[_i, "expression"] = _nm
                            st.session_state["fit_targets_df"] = _df
                            _thd = (st.session_state["fit_thetas_df"].copy()
                                    if "fit_thetas_df" in st.session_state else
                                    pd.DataFrame(columns=["name", "lower", "upper", "log",
                                                          "initial", "prior \u03bc", "prior \u03c3"]))
                            _thd = pd.concat([_thd, pd.DataFrame([{
                                "name": _nm, "lower": "", "upper": "", "log": False,
                                "initial": "", "prior \u03bc": "", "prior \u03c3": ""}])], ignore_index=True)
                            st.session_state["fit_thetas_df"] = _thd
                            st.session_state["fit_spec_rev"] = _rev + 1
                            st.session_state["_share_clear"] = True
                            st.rerun()

                # -- Custom parameters (theta) --
                with st.expander("Custom parameters (\u03b8) \u2014 used by Derived rows", expanded=False):
                    st.caption("Estimated quantities referenced by Derived expressions. Name + bounds "
                               "(blank = unconstrained); give large-scale \u03b8 an initial. Optional "
                               "prior \u03bc/\u03c3. Unreferenced \u03b8 are ignored.")
                    if "fit_thetas_df" not in st.session_state:
                        st.session_state["fit_thetas_df"] = pd.DataFrame(
                            [{"name": "", "lower": "", "upper": "", "log": False, "initial": "",
                              "prior \u03bc": "", "prior \u03c3": ""}])
                    _th_ed = st.data_editor(
                        st.session_state["fit_thetas_df"], key=f"fit_thetas_editor_{_rev}",
                        num_rows="dynamic", hide_index=True, width="stretch",
                        column_config={
                            "name": st.column_config.TextColumn("\u03b8 name"),
                            "lower": st.column_config.TextColumn("lower", help="Blank = unbounded below. e.g. 1e7"),
                            "upper": st.column_config.TextColumn("upper", help="Blank = unbounded above. e.g. 1e10"),
                            "log": st.column_config.CheckboxColumn("log", help="Log-space search."),
                            "initial": st.column_config.TextColumn("initial (start)",
                                                                   help="Start value \u2014 important when unbounded. e.g. 1e9"),
                            "prior \u03bc": st.column_config.TextColumn("prior \u03bc", help="Optional prior mean -> MAP."),
                            "prior \u03c3": st.column_config.TextColumn("prior \u03c3", help="Prior stdev."),
                        })
                    _th_names = [str(r["name"]).strip() for r in _th_ed.to_dict("records")
                                 if str(r.get("name", "")).strip()]

                # -- Covariate effects: per-arm parameter links (NONMEM/Monolix style) --
                _covariate_effects = []
                with st.expander("Covariate effects — per-arm parameter links (advanced)", expanded=False):
                    if not _cov_names:
                        st.caption("No covariates available. A covariate is a **numeric grouping "
                                   "column** (e.g. temperature, inoculum) or **MOI** (a phage dose). "
                                   "Add one as a grouping variable above to link a parameter to it.")
                    else:
                        st.caption("Let a parameter vary **per arm** by a covariate: "
                                   "θᵢ = θ_ref · (cov/ref)^β (power), or linear / "
                                   "exponential. The slope **β is estimated**; the parameter's own "
                                   "row above sets θ_ref. Available covariates: "
                                   + ", ".join(f"`{c}`" for c in _cov_names) + ".")
                        if "fit_covariates_df" not in st.session_state:
                            st.session_state["fit_covariates_df"] = pd.DataFrame(
                                [{"parameter": "", "covariate": "", "form": "power",
                                  "ref": "", "β lower": "-3", "β upper": "3", "β init": "0"}])
                        _cov_ed = st.data_editor(
                            st.session_state["fit_covariates_df"], key=f"fit_covariates_editor_{_rev}",
                            num_rows="dynamic", hide_index=True, width="stretch",
                            column_config={
                                "parameter": st.column_config.SelectboxColumn(
                                    "parameter", options=[l for (l, *_r) in _targets_cat],
                                    help="Model parameter to modulate (its row above sets θ_ref)."),
                                "covariate": st.column_config.SelectboxColumn(
                                    "covariate", options=_cov_names, help="Per-arm covariate."),
                                "form": st.column_config.SelectboxColumn(
                                    "form", options=list(_COV_FORMS), required=True),
                                "ref": st.column_config.TextColumn(
                                    "ref (cov₀)", help="Reference covariate value; multiplier=1 here. Blank=1."),
                                "β lower": st.column_config.TextColumn("β lower", help="Blank = -3."),
                                "β upper": st.column_config.TextColumn("β upper", help="Blank = 3."),
                                "β init": st.column_config.TextColumn("β init", help="Start slope. Blank = 0."),
                            })
                        for r in _cov_ed.to_dict("records"):
                            _plab = str(r.get("parameter", "")).strip()
                            _cvn = str(r.get("covariate", "")).strip()
                            if not _plab or not _cvn or _plab not in _label_to_path:
                                continue
                            _blo, _bhi = _pf(r.get("β lower")), _pf(r.get("β upper"))
                            _covariate_effects.append({
                                "path": _label_to_path[_plab], "covariate": _cvn,
                                "form": str(r.get("form", "power") or "power"),
                                "ref": (_pf(r.get("ref")) or 1.0),
                                "beta_lo": (_blo if _blo is not None else -3.0),
                                "beta_hi": (_bhi if _bhi is not None else 3.0),
                                "beta_init": (_pf(r.get("β init")) or 0.0),
                            })

                # -- Assembly: role -> free targets + Derived mappings; thetas --
                _targets, _mappings, _bad = [], [], []
                for r in _tg_ed.to_dict("records"):
                    _role = str(r.get("role", "Fixed"))
                    _log = bool(r["log"])
                    _lo = _bound(r.get("lower"), _log, True)
                    _hi = _bound(r.get("upper"), _log, False)
                    if np.isfinite(_lo) and np.isfinite(_hi) and _lo >= _hi:
                        if _role == "Free":
                            st.warning(f"'{_path_label.get(r['path'], r['path'])}': lower \u2265 upper \u2014 unbounded.")
                        _lo, _hi = (0.0 if _log else float("-inf")), float("inf")
                    _val = _pf(r.get("value"))
                    _targets.append({"path": r["path"], "free": (_role == "Free"),
                                     "value": (_val if _val is not None else 0.0),
                                     "lo": _lo, "hi": _hi, "log": _log,
                                     "prior_mu": _pf(r.get("prior \u03bc")), "prior_sd": _pf(r.get("prior \u03c3"))})
                    if _role == "Derived":
                        _ex = str(r.get("expression", "")).strip()
                        if not _ex:
                            _bad.append(f"{r['path']}: Derived but no expression")
                        else:
                            _ok, _msg = _nls.validate_expr(_ex, _th_names)
                            if _ok:
                                _mappings.append({"path": r["path"], "expr": _ex})
                            else:
                                _bad.append(f"{r['path']} = {_ex}: {_msg}")
                for _b in _bad:
                    st.warning(_b)
                _used = {n for n in _th_names
                         if any(re.search(rf"\b{re.escape(n)}\b", m["expr"]) for m in _mappings)}
                _thetas, _theta_bad = [], []
                for r in _th_ed.to_dict("records"):
                    _nm = str(r.get("name", "")).strip()
                    if _nm not in _used:
                        continue
                    _log = bool(r.get("log"))
                    _lo = _bound(r.get("lower"), _log, True)
                    _hi = _bound(r.get("upper"), _log, False)
                    if np.isfinite(_lo) and np.isfinite(_hi) and _lo >= _hi:
                        _theta_bad.append(_nm)
                        continue
                    _thetas.append({"name": _nm, "lo": _lo, "hi": _hi, "log": _log,
                                    "initial": _pf(r.get("initial")),
                                    "prior_mu": _pf(r.get("prior \u03bc")), "prior_sd": _pf(r.get("prior \u03c3"))})
                if _theta_bad:
                    st.warning("\u03b8 with lower \u2265 upper: " + ", ".join(sorted(set(_theta_bad))))

                # -- Fit spec (text) -- two-way with the table above --
                _sp_pending = st.session_state.pop("_spec_pending", None)
                if _sp_pending is not None:
                    st.session_state["fit_spec_text"] = _sp_pending
                with st.expander("Fit spec (text) \u2014 two-way with the table above", expanded=False):
                    st.caption("Compact text form; round-trips with the table. Grammar (one per line): "
                               "`free <path> [init= bounds=LO..HI prior=MU,SD log]`, `fix <path> = <v>`, "
                               "`theta <name> [...]`, `map <path> = <expr>`. Generate-from-table lists the "
                               "model's parameters as a comment header.")
                    _spec_txt = st.text_area("Spec", key="fit_spec_text", height=220,
                                             placeholder="theta g bounds=0.1..3.0 prior=1.2,0.3\n"
                                                         "map growth_rates[0] = g\nmap growth_rates[1] = g*(1-cost)")
                    _spc = st.columns(2)
                    if _spc[0].button("\u2191 Generate from table", key="fit_spec_from_tables", width="stretch"):
                        st.session_state["_spec_pending"] = _nls.serialize_fit_spec(
                            _targets, _thetas, _mappings, _targets_cat)
                        st.rerun()
                    if _spc[1].button("\u2193 Apply to table", key="fit_spec_to_tables",
                                      type="primary", width="stretch"):
                        _tdf, _thdf, _errs = _nls.parse_fit_spec(_spec_txt or "", _targets_cat)
                        if _errs:
                            for _e in _errs:
                                st.error(_e)
                        else:
                            st.session_state["fit_targets_df"] = _tdf
                            st.session_state["fit_thetas_df"] = _thdf
                            st.session_state["fit_spec_rev"] = _rev + 1
                            st.success("Applied to the table above.")
                            st.rerun()

                # OD link (CFU per OD unit): pbisim-fit reads it off the config, so it is
                # REQUIRED to fit 'od'. Use the ONE od_to_cfu (§4 Overlay) that also drives
                # B₀ and the overlay — no separate fit-only factor. Surface the data's own
                # median CFU/OD ratio as a suggestion the user can adopt in §4.
                _od_link = None
                if "od" in _sel_obs:
                    _od_link = float(st.session_state.get("int_od_to_cfu_conversion_factor", 2e8) or 2e8)
                    _od_suggest = _nls.estimate_od_to_cfu(_agg, _sel_arms)
                    _hint = (f" — data suggests ≈ **{_od_suggest:.2e}** (set it in §4 to adopt)"
                             if _od_suggest else "")
                    st.caption(f"OD fit uses od_to_cfu = **{_od_link:.2e}** (from §4 Overlay){_hint}.")
                _c1, _c2 = st.columns(2)
                with _c1:
                    _restarts = int(st.number_input("Restarts", 1, 10, 3, key="fit_nls_restarts"))
                with _c2:
                    _maxnfev = int(st.number_input("Max evaluations / restart", 50, 2000, 300,
                                                   step=50, key="fit_nls_maxnfev"))
                if _nls.has_unbounded(_targets, _thetas):
                    st.caption("Unbounded parameter(s) present → the fit uses a **single start** "
                               "(multi-start needs finite bounds).")
                _any_free = any(t["free"] for t in _targets) or bool(_thetas)
                _job = st.session_state.get("fit_job")
                _running = bool(_job and _job.get("status") == "running")
                # The fit runs in a separate PROCESS (not a thread): a ~1-min NLS fit is
                # CPU-heavy and a thread would starve Streamlit's main thread via the GIL,
                # freezing navigation. A process has its own interpreter/GIL, so the app
                # stays responsive and navigable while the fit runs; Stop terminates it.
                if st.button("Run NLS fit", key="fit_run_nls", type="primary",
                             width="stretch", disabled=_running):
                    if _theta_bad:
                        st.error("theta(s) have lower ≥ upper: "
                                 + ", ".join(sorted(set(_theta_bad))) + " — fix or blank a bound.")
                    elif not _any_free:
                        st.error("Set at least one parameter's role to **Free** "
                                 "(or add a referenced θ) — nothing is being estimated.")
                    elif not _sel_arms or not _sel_obs:
                        st.error("Select at least one arm and observable above.")
                    else:
                        try:
                            _ds_fit = _nls.build_dataset(_agg, _sel_arms, _sel_obs, _arm_cond,
                                                         od_to_cfu=_od_link, dose_unit=_dose_unit,
                                                         arm_doses=_arm_doses,
                                                         arm_covariates=_arm_covariates)
                            _proc, _q = _start_fit_process({
                                "base_config": _fit_cfg, "targets": _targets, "thetas": _thetas,
                                "mappings": _mappings, "dataset": _ds_fit, "obs_keys": list(_sel_obs),
                                "od_to_cfu": _od_link, "n_restarts": _restarts,
                                "max_nfev": _maxnfev, "estimate_b0": _estimate_b0,
                                "obs_compartments": dict(_obs_comp),
                                "covariate_effects": list(_covariate_effects),
                            })
                            st.session_state["fit_job"] = {
                                "status": "running", "t0": _time.time(),
                                "proc": _proc, "queue": _q, "result": None,
                                "od_link": _od_link,
                                "targets": _targets, "thetas": _thetas, "mappings": _mappings,
                                # post-processing context captured now, so edits during the
                                # fit don't change how the result is interpreted/overlaid.
                                "fit_model": _fit_model, "path_label": dict(_path_label),
                                "fB": _fB, "fP": _fP, "fS": _fS, "fmk": _fmk,
                                "ovl_ctx": dict(_ovl_ctx),
                                # short labels of what's being estimated (for the status line)
                                "param_preview": (
                                    [_path_label.get(t["path"], t["path"]) for t in _targets
                                     if t["free"] and t["path"] not in {m["path"] for m in _mappings}]
                                    + [f"θ {th['name']}" for th in _thetas]),
                            }
                            st.rerun()
                        except Exception as _fe:
                            st.error(f"Could not start fit: {_fe}")

                # The running banner, the Stop button, and the done→result harvest
                # are handled by _monitor_fit_job() at the TOP of render() (so they
                # survive navigation and can't be hidden by an error in this section).
                # Here we only surface the result table once it's ready.
                _fr = st.session_state.get("calib_fit_result")
                if _fr:
                    _rows = []
                    for _pm in _fr.get("params", []):
                        _mv = _fr["map"].get(_pm["key"])
                        _lo, _hi = (_fr["ci"].get(_pm["key"]) or [None, None])
                        _rows.append({
                            "parameter": _pm["label"],
                            "MAP": (f"{_mv:.4g}" if _mv is not None else "—"),
                            "95% CI low": (f"{_lo:.4g}" if _lo is not None else "—"),
                            "95% CI high": (f"{_hi:.4g}" if _hi is not None else "—"),
                        })
                    st.dataframe(pd.DataFrame(_rows), hide_index=True, width="stretch")

                    # Fit vs. curve-stripping agreement (if a strip result exists)
                    _sr = st.session_state.get("calib_strip_result")
                    if _sr is not None and _fr.get("map"):
                        with st.expander("🔎 Fit vs. curve-stripping agreement"):
                            st.caption("Compares each fitted parameter to its analytic curve-"
                                       "stripping estimate. Bands are wide by design — an out-of-band "
                                       "ratio is a caution flag (possible local minimum / model or "
                                       "data mismatch), not necessarily an error.")
                            try:
                                _rep = _nls.strip_sanity_check(_sr, _fr["map"])
                                if not _rep.rows:
                                    st.caption("No parameters overlap between the fit and the "
                                               "stripping estimates.")
                                else:
                                    if _rep.ok:
                                        st.success("All checked parameters agree with the analytic estimates.")
                                    else:
                                        st.warning("Outside the analytic band: " + ", ".join(_rep.flagged))
                                    st.dataframe(pd.DataFrame([
                                        {"parameter": _rw[0], "stripped": f"{_rw[1]:.4g}",
                                         "fitted": f"{_rw[2]:.4g}", "ratio": f"{_rw[3]:.2f}",
                                         "band": f"{_rw[4][0]:g}–{_rw[4][1]:g}",
                                         "ok": ("✓" if _rw[5] else "⚠")}
                                        for _rw in _rep.rows]), hide_index=True, width="stretch")
                            except Exception as _re:
                                st.caption(f"(Agreement check unavailable — {_re})")

                    st.caption("The overlay above already shows the fitted curves. **Apply** writes "
                               "these values into the model the fit ran against, so every task can reuse it.")
                    if st.button("Apply fitted values to model", key="fit_apply_map",
                                 width="stretch"):
                        _tgt = _fr.get("model", WORKING_DRAFT_LABEL)
                        # Make the builder the fit's base model, then write the FITTED
                        # CONFIG's resolved parameters onto it (robust for reparameterized
                        # fits, where the MAP keys are theta names rather than paths).
                        _tsnap = resolve_model_snapshot(_tgt)
                        if _tsnap is not None:
                            apply_model_to_state(_tsnap)
                        _fcfg = st.session_state.get("calib_fitted_config")
                        if _fcfg is not None:
                            _apply_config_to_state(_fcfg)
                        # Drop the manual-tuning + model-builder widgets so they
                        # re-seed from the updated int_* dicts — otherwise their stale
                        # keyed values would overwrite the just-applied fit on the next
                        # rerun (the builder writes each widget back into the dict).
                        for _k in [k for k in list(st.session_state.keys())
                                   if k.startswith("fit_edit_")
                                   or k.startswith(_BUILDER_WIDGET_PREFIXES)]:
                            st.session_state.pop(_k, None)
                        if _tgt in st.session_state.user_models:
                            st.session_state.user_models[_tgt]["state"] = dump_model()
                            st.session_state["_pending_active_model"] = _tgt
                            st.success(f"Applied fitted values to saved model '{_tgt}' and updated it.")
                        else:
                            _reason = ("demo model is read-only" if _tgt != WORKING_DRAFT_LABEL
                                       else "working draft")
                            st.success(f"Applied fitted values to the builder ({_reason}) — "
                                       "save it as a Model (sidebar) to keep.")
                        st.rerun()
            except Exception as _e:
                st.caption(f"(Load data and select parameters to run a fit — {_e})")

        # ── 6. Save the calibrated model ─────────────────────────────────────
        if _ds:
            st.markdown("### 6 · Save the calibrated model")
            st.caption("Parameter edits in the tuning panel are **already applied** to the live "
                       "Interactive-Simulator model. Save the whole calibrated configuration as a "
                       "Scenario to reload it later (from the Library page), or save individual "
                       "strains/phages as Parts in the Library.")
            _cs1, _cs2 = st.columns([3, 2])
            with _cs1:
                _cal_name = st.text_input("Scenario name", value="calibrated", key="fit_save_name")
            with _cs2:
                st.markdown("<br>", unsafe_allow_html=True)
                if st.button("Save calibrated config as Scenario", key="fit_save_scenario", width="stretch"):
                    _nm = (_cal_name or "").strip()
                    if not _nm:
                        st.error("Enter a scenario name.")
                    else:
                        _scen = st.session_state.user_scenarios
                        _scen[_nm] = {
                            "annotation": "Saved from Calibration",
                            "schema_version": SCENARIO_SCHEMA_VERSION,
                            "state": dump_state_to_scenario(),
                        }
                        st.session_state.user_scenarios = _scen
                        st.success(f"Saved scenario '{_nm}'. Reload it from the Library page.")

        if st.button("Clear dataset", key="fit_clear"):
            st.session_state.fit_dataset = None
            st.session_state.fit_config = {}
            st.session_state.calib_overlay_result = None
            st.rerun()

    # Save the current Calibration widget selections to the persistent config so they
    # survive navigation (see the re-seed block at the top of this page). The
    # parameter-tuning widgets (fit_edit_*) are excluded: they mirror the live
    # int_strains/int_phages dicts (authoritative + already persistent), so caching
    # and re-seeding them would let a stale copy override edits made elsewhere.
    for _wk in list(st.session_state.keys()):
        if (_wk.startswith("fit_") and _wk not in _FIT_NOPERSIST
                and not _wk.startswith("fit_edit_")
                and not _wk.startswith("fit_targets")   # parameter table (df + editor + sig)
                and not _wk.startswith("fit_thetas")     # theta table
                and not _wk.startswith("fit_share")      # Share helper (transient; self-clears)
                and not _wk.startswith("fit_map")        # (removed) mapping table
                and not _wk.startswith("fit_shrq_")      # (removed) shared-param quick builder
                and not _wk.startswith("fit_shr_")       # (removed) sharing-builder widgets
                and not _wk.startswith("fit_der_")):     # (removed) derived-link builder widgets
            st.session_state.fit_config[_wk] = st.session_state[_wk]


# ── AI Simulation Assistant Page ──────────────────────────────────────────────
