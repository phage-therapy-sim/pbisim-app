"""Rendered by app.py when this page is selected."""
from pbisim_app.common import *  # noqa: F401,F403


_SWEEP_MISSING = object()  # sentinel: session key absent before the sweep set it


def _short_sweep_label(label):
    """Compact a sweep-parameter label for trajectory legends (strip the ' - strain'
    suffix and any parenthetical)."""
    return label.rsplit(" - ", 1)[0].split("(")[0].strip()


def _fmt_sweep_value(v):
    """Format a coupled-sweep step value: numeric → 2 sig-figs, categorical → verbatim."""
    if isinstance(v, bool):
        return str(v)
    if isinstance(v, (int, float)):
        return f"{v:.2g}"
    return str(v)


def render():
    theme_mode = st.session_state.get("theme_mode", "Light")
    st.title("Model Parameter Sweeps")
    st.caption("Sweep any model parameter in 1D or 2D and visualize cellular trajectories and outcome heatmaps.")
    # Run the whole sweep against a CHOSEN Model (frozen snapshot or live draft), so
    # the base config, parameter picker, and solves are all consistent and immune to
    # unrelated builder edits.
    _sel = page_model_selector("psweep")
    with model_config_context(resolve_model_snapshot(_sel)):
        _render_body(theme_mode)


def _render_body(theme_mode):
    # Keep the sweep controls alive across navigation (re-seed before they render).
    reseed_widget_config("param_sweep_config", ("p1_", "p2_", "ps_", "pc_"))

    # Build nominal configuration
    try:
        nominal_config, initial_B, initial_P, initial_S, model_kwargs = build_nominal_config_from_gui()
    except Exception as e:
        st.error(f"Nominal configuration build failed. Please configure it in the Interactive Simulator first. Error: {e}")
        st.stop()

    strains_gui = st.session_state.get("int_strains", [])
    phages_gui = st.session_state.get("int_phages", [])
    antibiotics_gui = st.session_state.get("int_antibiotics", [])
    theme_mode = st.session_state.get("theme_mode", "Light")

    st.markdown(
        "<div class='info-banner'>Sweeps are based on the <b>Model selected above</b>. "
        "Toggle 1D vs 2D sweep type, select parameters, and run.</div>",
        unsafe_allow_html=True,
    )

    sweep_type = st.radio("Sweep Dimension", ["1D Sweep", "2D Sweep", "Coupled (linked)"], horizontal=True, key="ps_sweep_type")

    sweep_params = get_sweep_parameters(nominal_config, strains_gui, phages_gui, antibiotics_gui)
    param_labels = sorted(list(sweep_params.keys()))
    # Group the (many) sweepable parameters by builder-style category so they're
    # easy to find. `_cat_ordered` = all labels in category order (for the coupled
    # multiselect, so related params cluster).
    _cats = categorize_sweep_params(sweep_params)
    _cat_names = list(_cats)
    _cat_ordered = [lbl for c in _cat_names for lbl in _cats[c]]

    def _pick_param(prefix, cat_label, param_label):
        """A category selectbox + a parameter selectbox filtered to that category.
        Drops a persisted parameter choice that no longer belongs to the category."""
        cat = st.selectbox(cat_label, _cat_names, key=f"{prefix}_cat")
        opts = _cats.get(cat, [])
        if st.session_state.get(f"{prefix}_sweep_label") not in opts:
            st.session_state.pop(f"{prefix}_sweep_label", None)
        return st.selectbox(param_label, opts, key=f"{prefix}_sweep_label")

    col_setup, col_run = st.columns([1, 2])

    with col_setup:
        st.markdown("### Configure Parameters")
        
        if sweep_type == "1D Sweep":
            param1_label = _pick_param("p1", "Category", "Parameter")
            meta1 = sweep_params[param1_label]
            _is_categorical = meta1.get("type") == "categorical"

            if _is_categorical:
                # Categorical (signal-function) sweep — pick which discrete options to
                # compare. Each rebuilds+re-solves the model (the signal choice resolves
                # at BUILD time, so we can't mutate a built config), so keep it modest.
                _cat_opts = list(meta1["options"].keys())
                if st.session_state.get("ps_1d_cat_opts") and any(
                        o not in _cat_opts for o in st.session_state["ps_1d_cat_opts"]):
                    st.session_state.pop("ps_1d_cat_opts", None)  # options changed
                chosen_opts = st.multiselect(
                    "Options to compare", _cat_opts, default=_cat_opts, key="ps_1d_cat_opts")
                st.caption("Each option rebuilds and re-solves the model, then overlays the "
                           "trajectories and compares outcome metrics. A nutrient-dependent "
                           "signal under a signal-frozen environment may coerce to *constant* "
                           "— such runs are labelled.")
                run_sweep = st.button("Run Categorical Sweep", width="stretch", type="primary")

            if not _is_categorical:
                # Re-autoscale the range widgets when the swept parameter changes (their
                # value= default depends on the nominal value; a persisted key would pin it).
                if st.session_state.get("_ps_1d_last_param") != param1_label:
                    for _k in ("ps_1d_min", "ps_1d_max", "ps_1d_steps", "ps_1d_spacing"):
                        st.session_state.pop(_k, None)
                        st.session_state.setdefault("param_sweep_config", {}).pop(_k, None)
                    st.session_state["_ps_1d_last_param"] = param1_label

                # Default values
                default_val = 1e-9
                if meta1["type"] == "scalar":
                    default_val = getattr(nominal_config, meta1["field"])
                    if default_val is None:  # e.g. dormancy_monod_constant when inheriting
                        default_val = meta1.get("default", 1.0)
                elif meta1["type"] == "dimension":
                    default_val = getattr(nominal_config, meta1["field"])
                elif meta1["type"] == "array1d":
                    default_val = getattr(nominal_config, meta1["field"])[meta1["index"]]
                elif meta1["type"] == "array1d_or_none":
                    arr = getattr(nominal_config, meta1["field"])
                    default_val = arr[meta1["index"]] if arr is not None else 0.0
                elif meta1["type"] in ("array1d_broadcast", "array1d_broadcast_or_none"):
                    arr = getattr(nominal_config, meta1["field"])
                    default_val = float(arr[0]) if arr is not None and len(arr) else 0.0
                elif meta1["type"] == "initial_B_broadcast":
                    default_val = float(initial_B[0]) if len(initial_B) else 1e7
                elif meta1["type"] == "array2d":
                    default_val = getattr(nominal_config, meta1["field"])[meta1["index_row"], meta1["index_col"]]
                elif meta1["type"] == "array2d_broadcast":
                    _arr = getattr(nominal_config, meta1["field"])
                    default_val = float(_arr.flat[0]) if _arr is not None and _arr.size else 0.0
                elif meta1["type"] == "pd_array2d_broadcast":
                    _arr = getattr(nominal_config.pd_config, meta1["field"])
                    default_val = float(_arr.flat[0]) if _arr is not None and _arr.size else 0.0
                elif meta1["type"] == "pk_array1d":
                    pk_config = nominal_config.phage_pk_config or nominal_config.pk_config
                    default_val = getattr(pk_config, meta1["field"])[meta1["index"]]
                elif meta1["type"] == "pd_array2d":
                    default_val = getattr(nominal_config.pd_config, meta1["field"])[meta1["index_row"], meta1["index_col"]]
                elif meta1["type"] == "initial_B":
                    default_val = initial_B[meta1["index"]]
                elif meta1["type"] == "initial_P":
                    default_val = initial_P[meta1["index"]]
                elif meta1["type"] == "initial_S":
                    default_val = initial_S
                elif meta1["type"] == "prerun":
                    default_val = st.session_state.get("int_t_prerun", 0.0) or meta1.get("default", 24.0)
                elif meta1["type"] == "mutation":
                    _M = getattr(nominal_config, meta1.get("field", "mutation_rates"), None)
                    default_val = float(_M[meta1["dest"], meta1["origin"]]) if _M is not None else 0.0
                    if not default_val:
                        default_val = meta1.get("default", 1e-7)

                st.caption(f"Nominal Value: `{default_val:.2e}`" if isinstance(default_val, (int, float)) else f"Nominal Value: `{default_val}`")

                _is_dim = meta1["type"] == "dimension"
                c1, c2, c3 = st.columns(3)
                if _is_dim:
                    # integer compartment count: integer bounds ≥ 1
                    _dv = max(1, int(round(default_val)))
                    with c1:
                        min_val = st.number_input("Min Value", min_value=1, value=1, step=1, key="ps_1d_min")
                    with c2:
                        max_val = st.number_input("Max Value", min_value=1, value=max(_dv + 3, 5), step=1, key="ps_1d_max")
                    with c3:
                        steps = st.number_input("Steps", min_value=2, max_value=25, value=5, key="ps_1d_steps")
                else:
                    with c1:
                        min_val = st.number_input("Min Value", value=float(default_val * 0.1) if default_val > 0 else 0.0, format="%.2e", key="ps_1d_min")
                    with c2:
                        max_val = st.number_input("Max Value", value=float(default_val * 10.0) if default_val > 0 else 1.0, format="%.2e", key="ps_1d_max")
                    with c3:
                        steps = st.number_input("Steps", min_value=2, max_value=25, value=5, key="ps_1d_steps")

                spacing = st.selectbox("Spacing", ["Linear", "Logarithmic"], key="ps_1d_spacing")
                run_sweep = st.button("Run 1D Sweep", width="stretch", type="primary")

        elif sweep_type == "2D Sweep":
            param1_label = _pick_param("p1", "Parameter 1 category (X-axis)", "Parameter 1")
            meta1 = sweep_params[param1_label]
            
            c1, c2, c3 = st.columns(3)
            with c1:
                min_val = st.number_input("P1 Min Value", value=1e-10, format="%.2e", key="p1_min")
            with c2:
                max_val = st.number_input("P1 Max Value", value=1e-6, format="%.2e", key="p1_max")
            with c3:
                steps = st.number_input("P1 Steps", min_value=2, max_value=10, value=3, key="p1_steps")
            spacing = st.selectbox("P1 Spacing", ["Linear", "Logarithmic"], key="p1_spacing")
            
            st.markdown("---")
            param2_label = _pick_param("p2", "Parameter 2 category (Y-axis)", "Parameter 2")
            meta2 = sweep_params[param2_label]
            
            c4, c5, c6 = st.columns(3)
            with c4:
                min_val2 = st.number_input("P2 Min Value", value=10.0, format="%.2e", key="p2_min")
            with c5:
                max_val2 = st.number_input("P2 Max Value", value=200.0, format="%.2e", key="p2_max")
            with c6:
                steps2 = st.number_input("P2 Steps", min_value=2, max_value=10, value=3, key="p2_steps")
            spacing2 = st.selectbox("P2 Spacing", ["Linear", "Logarithmic"], key="p2_spacing")

            # 2D axes are both NUMERIC (heatmap grid). A signal-function parameter can't
            # be a 2D axis — steer the user to the modes that do support it.
            _2d_cat = [l for l, m in ((param1_label, meta1), (param2_label, meta2))
                       if m.get("type") == "categorical"]
            if _2d_cat:
                st.warning(
                    "2D sweeps vary two **numeric** parameters (the heatmap axes). Signal-"
                    "function parameter(s) — " + ", ".join(_2d_cat) + " — can't be a 2D axis. "
                    "Use a **1D** sweep to compare signal-function options, or a **Coupled "
                    "(linked)** sweep to pair a signal function with a continuous parameter.")
            run_sweep = st.button("Run 2D Sweep", width="stretch", type="primary",
                                  disabled=bool(_2d_cat))

        else:  # Coupled (linked) sweep
            st.caption(
                "Sweep several parameters **together**: pick the parameters and give each a "
                "value series of the **same length**. At step *k*, value[k] of every parameter "
                "is applied at once. Use the **(ALL strains)** parameters to share one value "
                "across strains. Example — (dormancy rate, resuscitation rate) = "
                "(0, 1), (0.5, 0.5), (1, 0)."
            )
            st.caption(
                "You can mix a **categorical** signal-function parameter (e.g. Lysis "
                "signal function) with **continuous** ones (e.g. Initial Phage Density): "
                "**choose the categorical options from a dropdown** (each = one step, in the "
                "order picked; the run rebuilds the model per step) and give the continuous "
                "parameter(s) an equal-length value list. Example — (Lysis signal function, "
                "Initial Phage Density) = (constant, 1e6), (nutrient, 1e8)."
            )
            coupled_labels = st.multiselect("Parameters to sweep together (grouped by category)", _cat_ordered, key="pc_labels")
            for _ci, _lbl in enumerate(coupled_labels):
                _meta = sweep_params[_lbl]
                if _meta.get("type") == "categorical":
                    # One dropdown PER STEP (not a single multiselect) so an option can be
                    # REPEATED across steps — e.g. constant, constant, constant, nutrient,
                    # nutrient, nutrient to pair with 0, 1e3, 1e5, 0, 1e3, 1e5. The step
                    # count sets how many dropdowns render.
                    _opt_names = list(_meta["options"].keys())
                    _n = int(st.number_input(
                        f"'{_lbl}' — number of steps", min_value=1, max_value=25,
                        value=int(st.session_state.get(f"pc_catn_{_ci}", 3)), step=1,
                        key=f"pc_catn_{_ci}"))
                    st.caption("Choose an option per step (repeats allowed). Give the "
                               f"continuous parameter(s) **{_n}** value(s) to match.")
                    _sbcols = st.columns(min(_n, 4))
                    for _j in range(_n):
                        _sk = f"pc_catopt_{_ci}_{_j}"
                        # Sanitize a stale/invalid persisted value (option set changed) so
                        # the keyed selectbox doesn't raise.
                        if st.session_state.get(_sk) not in _opt_names:
                            st.session_state.pop(_sk, None)
                        _sbcols[_j % len(_sbcols)].selectbox(
                            f"Step {_j + 1}", _opt_names, key=_sk)
                else:
                    st.text_input(
                        f"Values for '{_lbl}' (comma-separated)",
                        key=f"pc_series_{_ci}",
                        placeholder="e.g. 0, 0.5, 1",
                    )
            run_sweep = st.button("Run Coupled Sweep", width="stretch", type="primary")

    with col_run:
        st.markdown("### Sweep Results")
        if run_sweep:
            # Warn once if the pre-run decimates the culture (death w/o dormancy).
            warn_if_prerun_collapses(nominal_config, initial_B, initial_S=initial_S)
            progress_bar = st.progress(0)
            status_text = st.empty()

            if sweep_type == "1D Sweep" and _is_categorical:
                # Categorical sweep: the signal-function choice resolves at BUILD time
                # (via session state + the *_kwargs() helpers), so we can't mutate a
                # built config — instead set the session key and rebuild+re-solve per
                # option (save/restore the key like the dose-response sweep does).
                if not chosen_opts:
                    st.error("Select at least one option to compare.")
                    st.stop()
                _sess_key = meta1["session_key"]
                _opt_map = meta1["options"]
                _saved = st.session_state.get(_sess_key)

                runs_outcomes = []
                trajectories = []           # (time, cfu_b, label) — culturable CFU (B+D)
                total_incl_trajectories = []  # (time, B+D+I+H, label)
                active_trajectories = []    # (time, B, label)
                phage_trajectories = []     # (time, total_free_phage, label)
                central_phage_trajectories = []  # (time, central/blood conc. Σ Pc/Vc, label)
                od_trajectories = []        # (time, od, label) — only if OD/debris enabled
                _od_enabled = st.session_state.get("int_debris_enabled", False)
                _seen = []                  # (label, cfu) for coercion/duplicate flagging
                _trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz

                try:
                    for idx, opt in enumerate(chosen_opts):
                        status_text.text(f"Running option {idx+1} of {len(chosen_opts)}: {opt}...")
                        st.session_state[_sess_key] = _opt_map[opt]
                        # Full rebuild + solve (handles prerun / debris / immunity itself).
                        result, _cfg = run_sim_from_gui_params()

                        _cfu = result.sum_prefixes("B", "D")
                        _total_incl = result.sum_prefixes("B", "D", "I", "H")
                        _active = result.sum_prefixes("B")
                        nadir_val = float(np.min(_cfu))
                        auc_val = float(_trapz(_cfu, result.time))
                        t_clear = time_to_clearance(
                            result, threshold=st.session_state.get("int_extinction_threshold", 1.0))
                        t_log_red = time_to_log_reduction(result, n_logs=2.0)

                        _lbl = str(opt)
                        # A nutrient-dependent signal under a signal-frozen environment
                        # coerces to constant, giving dynamics identical to another option
                        # — flag the collapse rather than silently overlap the lines.
                        for _plbl, _pcfu in _seen:
                            if _pcfu.shape == _cfu.shape and np.allclose(_pcfu, _cfu, rtol=1e-6, atol=0.0):
                                _lbl = f"{opt} (≡ {_plbl})"
                                break
                        _seen.append((_lbl, _cfu))

                        runs_outcomes.append({
                            "Option": _lbl,
                            "Nadir (cells/mL)": nadir_val,
                            "AUC (cells·h/mL)": auc_val,
                            "Clearance Time (h)": t_clear if t_clear is not None else np.nan,
                            "2-Log Red Time (h)": t_log_red if t_log_red is not None else np.nan,
                        })
                        trajectories.append((result.time, _cfu, _lbl))
                        total_incl_trajectories.append((result.time, _total_incl, _lbl))
                        active_trajectories.append((result.time, _active, _lbl))
                        phage_trajectories.append((result.time, result.sum_prefixes("P"), _lbl))
                        if phage_pk_enabled(phages_gui):
                            central_phage_trajectories.append((result.time, central_phage_total(result, phages_gui), _lbl))
                        if _od_enabled:
                            od_trajectories.append((result.time, _safe_od(result, _cfu), _lbl))
                        progress_bar.progress((idx + 1) / len(chosen_opts))
                finally:
                    # Restore the swept key so the sweep leaves no residue on the model.
                    if _saved is None:
                        st.session_state.pop(_sess_key, None)
                    else:
                        st.session_state[_sess_key] = _saved

                status_text.text("Sweep completed!")
                st.session_state.param_sweep_result = {
                    "type": "categorical",
                    "param1_label": param1_label,
                    "summary": runs_outcomes,
                    "trajectories": [(np.asarray(t), np.asarray(b), lbl) for t, b, lbl in trajectories],
                    "total_incl_trajectories": [(np.asarray(t), np.asarray(b), lbl) for t, b, lbl in total_incl_trajectories],
                    "active_trajectories": [(np.asarray(t), np.asarray(b), lbl) for t, b, lbl in active_trajectories],
                    "phage_trajectories": [(np.asarray(t), np.asarray(p), lbl) for t, p, lbl in phage_trajectories],
                    "central_phage_trajectories": [(np.asarray(t), np.asarray(p), lbl) for t, p, lbl in central_phage_trajectories],
                    "od_trajectories": [(np.asarray(t), np.asarray(o), lbl) for t, o, lbl in od_trajectories],
                }

            elif sweep_type == "1D Sweep":
                # Compute sweep values
                if spacing == "Logarithmic":
                    if min_val <= 0 or max_val <= 0:
                        st.error("Logarithmic spacing requires positive bounds.")
                        st.stop()
                    sweep_values = np.logspace(np.log10(min_val), np.log10(max_val), int(steps))
                else:
                    sweep_values = np.linspace(min_val, max_val, int(steps))
                # Compartment counts (n_depth / n_latent) are integers >= 1 — sweep unique
                # integers so the results aren't fractional/duplicated.
                if meta1["type"] == "dimension":
                    sweep_values = np.unique(np.clip(np.round(sweep_values).astype(int), 1, None))
                    st.caption("Integer parameter — swept over unique integer values ≥ 1: "
                               f"{list(sweep_values)}")

                runs_outcomes = []
                trajectories = [] # (time, cfu_b, label) — culturable CFU (B+D)
                total_incl_trajectories = [] # (time, B+D+I+H, label)
                active_trajectories = [] # (time, B, label)
                phage_trajectories = []  # (time, total_free_phage, label)
                central_phage_trajectories = []  # (time, central/blood conc. Σ Pc/Vc, label)
                od_trajectories = [] # (time, od, label) — only if OD/debris enabled
                _od_enabled = st.session_state.get("int_debris_enabled", False)

                for idx, val in enumerate(sweep_values):
                    status_text.text(f"Running simulation {idx+1} of {len(sweep_values)} (Value: {val:.2e})...")
                    
                    # Apply parameter
                    c_k, ib_k, ip_k, is_k, mk_k = apply_sweep_parameter(
                        val, meta1, nominal_config, initial_B, initial_P, initial_S, model_kwargs
                    )

                    # equilibrate pre-treatment prerun — carry the full stationary
                    # state (B, D, S, Imm), not just B/S (see run_sim_from_gui_params).
                    t_prerun = mk_k.pop("_t_prerun_override", None)
                    if t_prerun is None:
                        t_prerun = st.session_state.get("int_t_prerun", 0.0)
                    if t_prerun > 0:
                        ic = stationary_phase_ic(c_k, t_prerun=t_prerun, B0=ib_k, initial_S=is_k)
                        ib_k = ic.B
                        is_k = max(float(ic.S), 0.0)
                        if ic.D is not None:
                            mk_k["initial_D"] = ic.D
                        if ic.Imm is not None:
                            mk_k["initial_Imm"] = ic.Imm
                        _carry_prerun_debris(ic, mk_k)

                    model = PBIModel(c_k, initial_B=ib_k, initial_P=ip_k, initial_S=is_k, **mk_k)
                    result = solve_ode(model, t_end=st.session_state.get("int_t_end", 48.0), dt=st.session_state.get("int_dt", 0.25), method=st.session_state.get("int_solver_method", "BDF"), extinction_threshold=st.session_state.get("int_extinction_threshold", 1.0) or None, extinction_check_interval=st.session_state.get("int_extinction_check_interval", 0.0) or None)

                    # Compute metrics. Summary metrics use culturable CFU (B+D).
                    _cfu = result.sum_prefixes("B", "D")
                    _total_incl = result.sum_prefixes("B", "D", "I", "H")
                    _active = result.sum_prefixes("B")
                    total_bacteria = _cfu
                    nadir_val = np.min(total_bacteria)

                    _trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
                    auc_val = _trapz(total_bacteria, result.time)

                    t_clear = time_to_clearance(
                        result,
                        threshold=st.session_state.get("int_extinction_threshold", 1.0),
                    )
                    t_log_red = time_to_log_reduction(result, n_logs=2.0)

                    runs_outcomes.append({
                        param1_label: val,
                        "Nadir (cells/mL)": nadir_val,
                        "AUC (cells·h/mL)": auc_val,
                        "Clearance Time (h)": t_clear if t_clear is not None else np.nan,
                        "2-Log Red Time (h)": t_log_red if t_log_red is not None else np.nan
                    })
                    _lbl = f"{param1_label} = {val:.2e}"
                    trajectories.append((result.time, _cfu, _lbl))
                    total_incl_trajectories.append((result.time, _total_incl, _lbl))
                    active_trajectories.append((result.time, _active, _lbl))
                    phage_trajectories.append((result.time, result.sum_prefixes("P"), _lbl))
                    if phage_pk_enabled(phages_gui):
                        central_phage_trajectories.append((result.time, central_phage_total(result, phages_gui), _lbl))
                    if _od_enabled:
                        _od = (_safe_od(result, total_bacteria))
                        od_trajectories.append((result.time, _od, _lbl))
                    progress_bar.progress((idx + 1) / len(sweep_values))

                status_text.text("Sweep completed!")
                st.session_state.param_sweep_result = {
                    "type": "1D",
                    "param1_label": param1_label,
                    "spacing": spacing,
                    "summary": runs_outcomes,
                    "trajectories": [(np.asarray(t), np.asarray(b), lbl) for t, b, lbl in trajectories],
                    "total_incl_trajectories": [(np.asarray(t), np.asarray(b), lbl) for t, b, lbl in total_incl_trajectories],
                    "active_trajectories": [(np.asarray(t), np.asarray(b), lbl) for t, b, lbl in active_trajectories],
                    "phage_trajectories": [(np.asarray(t), np.asarray(p), lbl) for t, p, lbl in phage_trajectories],
                    "central_phage_trajectories": [(np.asarray(t), np.asarray(p), lbl) for t, p, lbl in central_phage_trajectories],
                    "od_trajectories": [(np.asarray(t), np.asarray(o), lbl) for t, o, lbl in od_trajectories],
                }

            elif sweep_type == "Coupled (linked)":
                # Parse each selected parameter's value series; all must be equal length.
                labels = st.session_state.get("pc_labels", [])
                if not labels:
                    st.error("Select at least one parameter to sweep together.")
                    st.stop()
                # Categorical (signal-function) params resolve at BUILD time (session key
                # + *_kwargs helpers), so those steps set the key and REBUILD the config;
                # continuous params mutate the built config via apply_sweep_parameter.
                _cat_labels = [l for l in labels if sweep_params[l].get("type") == "categorical"]
                _cont_labels = [l for l in labels if l not in _cat_labels]
                series = {}
                for _ci, lbl in enumerate(labels):
                    if lbl in _cat_labels:
                        # Options come from the per-step dropdowns (pc_catopt_{ci}_{j}),
                        # in step order, repeats allowed.
                        _opts = sweep_params[lbl]["options"]
                        _n = int(st.session_state.get(f"pc_catn_{_ci}", 0))
                        seq = [st.session_state.get(f"pc_catopt_{_ci}_{_j}") for _j in range(_n)]
                        seq = [o for o in seq if o in _opts]
                        if not seq:
                            st.error(f"Choose at least one option for '{lbl}' from its dropdowns.")
                            st.stop()
                        series[lbl] = seq                # canonical option-label strings
                    else:
                        vals = parse_comma_separated_series(st.session_state.get(f"pc_series_{_ci}", ""))
                        if not vals:
                            st.error(f"Provide a value series for '{lbl}'.")
                            st.stop()
                        series[lbl] = vals
                M = len(next(iter(series.values())))
                if any(len(v) != M for v in series.values()):
                    st.error("All value series must have the same number of points.")
                    st.stop()

                # Categorical steps rebuild from GUI, mutating their session keys — save
                # them so the sweep leaves no residue on the model.
                _saved_cat = {sweep_params[l]["session_key"]:
                              st.session_state.get(sweep_params[l]["session_key"], _SWEEP_MISSING)
                              for l in _cat_labels}

                runs_outcomes = []
                trajectories = []               # culturable CFU (B+D)
                total_incl_trajectories = []    # B+D+I+H
                active_trajectories = []        # B
                phage_trajectories = []  # (time, total_free_phage, label)
                central_phage_trajectories = []  # (time, central/blood conc. Σ Pc/Vc, label)
                od_trajectories = []
                _od_enabled = st.session_state.get("int_debris_enabled", False)
                try:
                    for k in range(M):
                        status_text.text(f"Running simulation {k+1} of {M}...")
                        if _cat_labels:
                            # Select each categorical option for this step, then REBUILD
                            # the config from GUI so the signal-function choice takes effect.
                            for lbl in _cat_labels:
                                _cm = sweep_params[lbl]
                                st.session_state[_cm["session_key"]] = _cm["options"][series[lbl][k]]
                            c_k, ib_k, ip_k, is_k, mk_k = build_nominal_config_from_gui()
                        else:
                            c_k, ib_k, ip_k, is_k, mk_k = nominal_config, initial_B, initial_P, initial_S, model_kwargs
                        # Continuous params mutate the (possibly rebuilt) config.
                        for lbl in _cont_labels:
                            c_k, ib_k, ip_k, is_k, mk_k = apply_sweep_parameter(
                                series[lbl][k], sweep_params[lbl], c_k, ib_k, ip_k, is_k, mk_k)
                        t_prerun = mk_k.pop("_t_prerun_override", None)
                        if t_prerun is None:
                            t_prerun = st.session_state.get("int_t_prerun", 0.0)
                        if t_prerun > 0:
                            ic = stationary_phase_ic(c_k, t_prerun=t_prerun, B0=ib_k, initial_S=is_k)
                            ib_k = ic.B
                            is_k = max(float(ic.S), 0.0)
                            if ic.D is not None:
                                mk_k["initial_D"] = ic.D
                            if ic.Imm is not None:
                                mk_k["initial_Imm"] = ic.Imm
                            _carry_prerun_debris(ic, mk_k)
                        model = PBIModel(c_k, initial_B=ib_k, initial_P=ip_k, initial_S=is_k, **mk_k)
                        result = solve_ode(model, t_end=st.session_state.get("int_t_end", 48.0), dt=st.session_state.get("int_dt", 0.25), method=st.session_state.get("int_solver_method", "BDF"), extinction_threshold=st.session_state.get("int_extinction_threshold", 1.0) or None, extinction_check_interval=st.session_state.get("int_extinction_check_interval", 0.0) or None)
                        _cfu = result.sum_prefixes("B", "D")
                        _total_incl = result.sum_prefixes("B", "D", "I", "H")
                        _active = result.sum_prefixes("B")
                        total_bacteria = _cfu
                        _trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
                        t_clear = time_to_clearance(result, threshold=st.session_state.get("int_extinction_threshold", 1.0))
                        _row = {lbl: series[lbl][k] for lbl in labels}
                        _row.update({
                            "Nadir (cells/mL)": np.min(total_bacteria),
                            "AUC (cells·h/mL)": _trapz(total_bacteria, result.time),
                            "Clearance Time (h)": t_clear if t_clear is not None else np.nan,
                            "2-Log Red Time (h)": time_to_log_reduction(result, n_logs=2.0) or np.nan,
                        })
                        runs_outcomes.append(_row)
                        _lbl_txt = ", ".join(f"{_short_sweep_label(s)}={_fmt_sweep_value(series[s][k])}" for s in labels)
                        _lbl = f"Step {k+1}: {_lbl_txt}"
                        trajectories.append((result.time, _cfu, _lbl))
                        total_incl_trajectories.append((result.time, _total_incl, _lbl))
                        active_trajectories.append((result.time, _active, _lbl))
                        phage_trajectories.append((result.time, result.sum_prefixes("P"), _lbl))
                        if phage_pk_enabled(phages_gui):
                            central_phage_trajectories.append((result.time, central_phage_total(result, phages_gui), _lbl))
                        if _od_enabled:
                            _od = (_safe_od(result, total_bacteria))
                            od_trajectories.append((result.time, _od, _lbl))
                        progress_bar.progress((k + 1) / M)
                finally:
                    # Restore any categorical session keys the rebuild mutated.
                    for _sk, _sv in _saved_cat.items():
                        if _sv is _SWEEP_MISSING:
                            st.session_state.pop(_sk, None)
                        else:
                            st.session_state[_sk] = _sv

                status_text.text("Sweep completed!")
                st.session_state.param_sweep_result = {
                    "type": "coupled",
                    "labels": list(labels),
                    "summary": runs_outcomes,
                    "trajectories": [(np.asarray(t), np.asarray(b), lbl) for t, b, lbl in trajectories],
                    "total_incl_trajectories": [(np.asarray(t), np.asarray(b), lbl) for t, b, lbl in total_incl_trajectories],
                    "active_trajectories": [(np.asarray(t), np.asarray(b), lbl) for t, b, lbl in active_trajectories],
                    "phage_trajectories": [(np.asarray(t), np.asarray(p), lbl) for t, p, lbl in phage_trajectories],
                    "central_phage_trajectories": [(np.asarray(t), np.asarray(p), lbl) for t, p, lbl in central_phage_trajectories],
                    "od_trajectories": [(np.asarray(t), np.asarray(o), lbl) for t, o, lbl in od_trajectories],
                }

            else:
                # 2D Sweep
                if spacing == "Logarithmic":
                    if min_val <= 0 or max_val <= 0:
                        st.error("Logarithmic spacing requires positive bounds.")
                        st.stop()
                    sweep_values1 = np.logspace(np.log10(min_val), np.log10(max_val), int(steps))
                else:
                    sweep_values1 = np.linspace(min_val, max_val, int(steps))

                if spacing2 == "Logarithmic":
                    if min_val2 <= 0 or max_val2 <= 0:
                        st.error("Logarithmic spacing requires positive bounds.")
                        st.stop()
                    sweep_values2 = np.logspace(np.log10(min_val2), np.log10(max_val2), int(steps2))
                else:
                    sweep_values2 = np.linspace(min_val2, max_val2, int(steps2))

                total_sims = len(sweep_values1) * len(sweep_values2)
                sim_idx = 0

                grid_auc = np.zeros((len(sweep_values2), len(sweep_values1)))
                grid_nadir = np.zeros((len(sweep_values2), len(sweep_values1)))
                grid_clear = np.zeros((len(sweep_values2), len(sweep_values1)))

                for i2, val2 in enumerate(sweep_values2):
                    for i1, val1 in enumerate(sweep_values1):
                        status_text.text(f"Running simulation {sim_idx+1} of {total_sims}...")
                        
                        # Apply both parameters
                        c_k, ib_k, ip_k, is_k, mk_k = apply_sweep_parameter(
                            val1, meta1, nominal_config, initial_B, initial_P, initial_S, model_kwargs
                        )
                        c_k, ib_k, ip_k, is_k, mk_k = apply_sweep_parameter(
                            val2, meta2, c_k, ib_k, ip_k, is_k, mk_k
                        )

                        # equilibrate — carry the full stationary state (B, D, S, Imm).
                        t_prerun = mk_k.pop("_t_prerun_override", None)
                        if t_prerun is None:
                            t_prerun = st.session_state.get("int_t_prerun", 0.0)
                        if t_prerun > 0:
                            ic = stationary_phase_ic(c_k, t_prerun=t_prerun, B0=ib_k, initial_S=is_k)
                            ib_k = ic.B
                            is_k = max(float(ic.S), 0.0)
                            if ic.D is not None:
                                mk_k["initial_D"] = ic.D
                            if ic.Imm is not None:
                                mk_k["initial_Imm"] = ic.Imm
                            _carry_prerun_debris(ic, mk_k)

                        model = PBIModel(c_k, initial_B=ib_k, initial_P=ip_k, initial_S=is_k, **mk_k)
                        result = solve_ode(model, t_end=st.session_state.get("int_t_end", 48.0), dt=st.session_state.get("int_dt", 0.25), method=st.session_state.get("int_solver_method", "BDF"), extinction_threshold=st.session_state.get("int_extinction_threshold", 1.0) or None, extinction_check_interval=st.session_state.get("int_extinction_check_interval", 0.0) or None)

                        # 2D heatmap metrics use culturable CFU (B+D).
                        total_bacteria = result.sum_prefixes("B", "D")
                        nadir_val = np.min(total_bacteria)

                        _trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
                        auc_val = _trapz(total_bacteria, result.time)

                        t_clear = time_to_clearance(
                            result,
                            threshold=st.session_state.get("int_extinction_threshold", 1.0),
                        )
                        
                        grid_auc[i2, i1] = auc_val
                        grid_nadir[i2, i1] = nadir_val
                        grid_clear[i2, i1] = t_clear if t_clear is not None else st.session_state.get("int_t_end", 48.0)

                        sim_idx += 1
                        progress_bar.progress(sim_idx / total_sims)

                status_text.text("Sweep completed!")
                st.session_state.param_sweep_result = {
                    "type": "2D",
                    "param1_label": param1_label, "param2_label": param2_label,
                    "spacing": spacing, "spacing2": spacing2,
                    "sweep_values1": np.asarray(sweep_values1), "sweep_values2": np.asarray(sweep_values2),
                    "grid_auc": grid_auc, "grid_nadir": grid_nadir, "grid_clear": grid_clear,
                }

        # Render the (persisted) parameter-sweep result if one exists.
        _ps = st.session_state.get("param_sweep_result")
        if _ps:
            import plotly.graph_objects as go
            if _ps["type"] == "1D":
                _p1 = _ps["param1_label"]
                df_summary = pd.DataFrame(_ps["summary"])
                _sweep_summary_tiles(df_summary)
                st.markdown("#### Summary of Runs")
                st.dataframe(
                    df_summary.style.format({
                        _p1: "{:.2e}", "Nadir (cells/mL)": "{:.2e}", "AUC (cells·h/mL)": "{:.2e}",
                        "Clearance Time (h)": "{:.1f}", "2-Log Red Time (h)": "{:.1f}"
                    }),
                    width="stretch")
                _traces = {"CFU — culturable (B+D)": ("log", _ps["trajectories"])}
                if _ps.get("total_incl_trajectories"):
                    _traces["Total incl. infected (B+D+I+H)"] = ("log", _ps["total_incl_trajectories"])
                if _ps.get("active_trajectories"):
                    _traces["Active only (B)"] = ("log", _ps["active_trajectories"])
                if _ps.get("phage_trajectories"):
                    _traces["Free phage — infection site (PFU/mL)"] = ("log", _ps["phage_trajectories"])
                if _ps.get("central_phage_trajectories"):
                    _traces["Central phage — blood (PFU/mL)"] = ("log", _ps["central_phage_trajectories"])
                if _ps.get("od_trajectories"):
                    _traces["Optical density (AU)"] = ("linear", _ps["od_trajectories"])
                plot_sweep_traces(_traces, "ps1d_traj", title_suffix="across runs")
                st.markdown("#### Outcome Metrics vs Parameter Value")
                fig_metric = go.Figure()
                fig_metric.add_trace(go.Scatter(x=df_summary[_p1], y=df_summary["AUC (cells·h/mL)"], mode="lines+markers", name="Bacterial AUC", yaxis="y1"))
                fig_metric.add_trace(go.Scatter(x=df_summary[_p1], y=df_summary["Nadir (cells/mL)"], mode="lines+markers", name="Nadir", yaxis="y2"))
                fig_metric.update_layout(
                    xaxis=dict(title=_p1, type="log" if _ps["spacing"] == "Logarithmic" else "linear"),
                    yaxis=dict(title="AUC (cells·h/mL)", type="log"),
                    yaxis2=dict(title="Nadir (cells/mL)", type="log", overlaying="y", side="right"),
                    template="plotly_white" if theme_mode == "Light" else "plotly_dark")
                st.plotly_chart(fig_metric, width="stretch")
            elif _ps["type"] == "coupled":
                _labels = _ps["labels"]
                df_summary = pd.DataFrame(_ps["summary"])
                _sweep_summary_tiles(df_summary)
                st.markdown("#### Summary of Runs (linked parameters)")
                # Only numeric columns get the scientific format — a categorical
                # (signal-function) linked param holds option-name strings.
                _fmt = {c: "{:.2e}" for c in _labels if pd.api.types.is_numeric_dtype(df_summary[c])}
                _fmt.update({"Nadir (cells/mL)": "{:.2e}", "AUC (cells·h/mL)": "{:.2e}",
                             "Clearance Time (h)": "{:.1f}", "2-Log Red Time (h)": "{:.1f}"})
                st.dataframe(df_summary.style.format(_fmt), width="stretch")
                _traces = {"CFU — culturable (B+D)": ("log", _ps["trajectories"])}
                if _ps.get("total_incl_trajectories"):
                    _traces["Total incl. infected (B+D+I+H)"] = ("log", _ps["total_incl_trajectories"])
                if _ps.get("active_trajectories"):
                    _traces["Active only (B)"] = ("log", _ps["active_trajectories"])
                if _ps.get("phage_trajectories"):
                    _traces["Free phage — infection site (PFU/mL)"] = ("log", _ps["phage_trajectories"])
                if _ps.get("central_phage_trajectories"):
                    _traces["Central phage — blood (PFU/mL)"] = ("log", _ps["central_phage_trajectories"])
                if _ps.get("od_trajectories"):
                    _traces["Optical density (AU)"] = ("linear", _ps["od_trajectories"])
                plot_sweep_traces(_traces, "pscoupled_traj", title_suffix="across runs")
                st.markdown("#### Outcome Metrics vs Step Index")
                _step = list(range(1, len(df_summary) + 1))
                fig_metric = go.Figure()
                fig_metric.add_trace(go.Scatter(x=_step, y=df_summary["AUC (cells·h/mL)"], mode="lines+markers", name="Bacterial AUC", yaxis="y1"))
                fig_metric.add_trace(go.Scatter(x=_step, y=df_summary["Nadir (cells/mL)"], mode="lines+markers", name="Nadir", yaxis="y2"))
                fig_metric.update_layout(
                    xaxis=dict(title="Step index"),
                    yaxis=dict(title="AUC (cells·h/mL)", type="log"),
                    yaxis2=dict(title="Nadir (cells/mL)", type="log", overlaying="y", side="right"),
                    template="plotly_white" if theme_mode == "Light" else "plotly_dark")
                st.plotly_chart(fig_metric, width="stretch")
            elif _ps["type"] == "categorical":
                _p1 = _ps["param1_label"]
                df_summary = pd.DataFrame(_ps["summary"])
                _sweep_summary_tiles(df_summary)
                st.markdown(f"#### Summary of Runs — {_p1}")
                st.dataframe(
                    df_summary.style.format({
                        "Nadir (cells/mL)": "{:.2e}", "AUC (cells·h/mL)": "{:.2e}",
                        "Clearance Time (h)": "{:.1f}", "2-Log Red Time (h)": "{:.1f}"}),
                    width="stretch")
                _traces = {"CFU — culturable (B+D)": ("log", _ps["trajectories"])}
                if _ps.get("total_incl_trajectories"):
                    _traces["Total incl. infected (B+D+I+H)"] = ("log", _ps["total_incl_trajectories"])
                if _ps.get("active_trajectories"):
                    _traces["Active only (B)"] = ("log", _ps["active_trajectories"])
                if _ps.get("phage_trajectories"):
                    _traces["Free phage — infection site (PFU/mL)"] = ("log", _ps["phage_trajectories"])
                if _ps.get("central_phage_trajectories"):
                    _traces["Central phage — blood (PFU/mL)"] = ("log", _ps["central_phage_trajectories"])
                if _ps.get("od_trajectories"):
                    _traces["Optical density (AU)"] = ("linear", _ps["od_trajectories"])
                plot_sweep_traces(_traces, "pscat_traj", title_suffix="across options")
                st.markdown("#### Outcome Metrics by Option")
                _opts = df_summary["Option"].astype(str).tolist()
                fig_metric = go.Figure()
                fig_metric.add_trace(go.Bar(x=_opts, y=df_summary["AUC (cells·h/mL)"], name="Bacterial AUC", yaxis="y1"))
                fig_metric.add_trace(go.Scatter(x=_opts, y=df_summary["Nadir (cells/mL)"], mode="markers", marker=dict(size=11, symbol="diamond"), name="Nadir", yaxis="y2"))
                fig_metric.update_layout(
                    xaxis=dict(title=_p1, type="category"),
                    yaxis=dict(title="AUC (cells·h/mL)", type="log"),
                    yaxis2=dict(title="Nadir (cells/mL)", type="log", overlaying="y", side="right"),
                    barmode="group",
                    template="plotly_white" if theme_mode == "Light" else "plotly_dark")
                st.plotly_chart(fig_metric, width="stretch")
            else:
                _p1, _p2 = _ps["param1_label"], _ps["param2_label"]
                _xt = "log" if _ps["spacing"] == "Logarithmic" else "linear"
                _yt = "log" if _ps["spacing2"] == "Logarithmic" else "linear"
                st.markdown("#### Outcome Heatmaps (2D Sweep)")
                h1, h2 = st.columns(2)
                with h1:
                    fig_auc = go.Figure(data=go.Contour(z=_ps["grid_auc"], x=_ps["sweep_values1"], y=_ps["sweep_values2"], colorscale="Viridis", colorbar=dict(title="AUC")))
                    fig_auc.update_layout(title="Bacterial AUC Heatmap", xaxis=dict(title=_p1, type=_xt), yaxis=dict(title=_p2, type=_yt), template="plotly_white" if theme_mode == "Light" else "plotly_dark")
                    st.plotly_chart(fig_auc, width="stretch")
                with h2:
                    fig_nadir = go.Figure(data=go.Contour(z=_ps["grid_nadir"], x=_ps["sweep_values1"], y=_ps["sweep_values2"], colorscale="Magma", colorbar=dict(title="Nadir")))
                    fig_nadir.update_layout(title="Bacterial Nadir Heatmap", xaxis=dict(title=_p1, type=_xt), yaxis=dict(title=_p2, type=_yt), template="plotly_white" if theme_mode == "Light" else "plotly_dark")
                    st.plotly_chart(fig_nadir, width="stretch")
        else:
            st.info("Configure parameters and click **Run Sweep** to start the analysis.")

    with st.expander("View Python Reproduction Code"):
        st.caption("Standalone script that reproduces this sweep — the recorded base model "
                   "plus a loop calling the app's own apply_sweep_parameter per value. "
                   "(1D numeric sweeps only — categorical signal-function sweeps rebuild the "
                   "model per option through the full builder and aren't exported.)")
        try:
            st.code(generate_param_sweep_reproduction_code(), language="python")
        except Exception as _e:
            st.warning(f"Reproduction code unavailable: {_e}")

    # Persist the sweep controls so they survive navigation (see reseed above).
    save_widget_config("param_sweep_config", ("p1_", "p2_", "ps_", "pc_"))

# ── Interactive Simulator Page ────────────────────────────────────────────────
