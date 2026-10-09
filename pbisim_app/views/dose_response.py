"""Rendered by app.py when this page is selected."""
from pbisim_app.common import *  # noqa: F401,F403


def render():
    theme_mode = st.session_state.get("theme_mode", "Light")
    st.title("Dose-Response Simulator")
    st.caption("Perform multi-drug dose-response sweeps with MOI scaling, vector padding, and raw time-series visualization.")
    # Run the sweep against a CHOSEN Model (frozen snapshot or live draft).
    _sel = page_model_selector("dresp")
    with model_config_context(resolve_model_snapshot(_sel)):
        _render_body(theme_mode)


def _render_body(theme_mode):
    # Keep the sweep controls alive across navigation (re-seed before they render).
    reseed_widget_config("dr_sweep_config", ("dr_sweep_",))

    strains = st.session_state.get("int_strains", [])
    phages = st.session_state.get("int_phages", [])
    antibiotics = st.session_state.get("int_antibiotics", [])
    theme_mode = st.session_state.get("theme_mode", "Light")

    st.markdown(
        "<div class='info-banner'>Sweeps are based on the <b>Model selected above</b>. "
        "Change biological parameters (e.g. growth rates, adsorption) via the builder / a saved Model.</div>",
        unsafe_allow_html=True,
    )

    # 1. Sweep configurations
    swept_inputs = {}
    swept_units = {}
    swept_repeat_configs = {}

    col_setup, col_run = st.columns([1, 2])

    with col_setup:
        st.markdown("### Configure Sweeps")
        
        # Phages
        for j, p in enumerate(phages):
            st.markdown(f"#### Phage {j}: {p['name']}")
            do_sweep = st.checkbox(f"Sweep Phage {j}", key=f"dr_sweep_phg_en_{j}", value=False)
            if do_sweep:
                series_str = st.text_input(
                    "Dose series (comma-separated)",
                    value="0, 1e3, 1e5, 1e7, 1e9",
                    key=f"dr_sweep_phg_series_{j}"
                )
                unit = st.selectbox(
                    "Dose Unit",
                    ["PFU (absolute)", "MOI (relative to B(0))"],
                    key=f"dr_sweep_phg_unit_{j}"
                )
                swept_inputs[f"phage_{j}"] = series_str
                swept_units[f"phage_{j}"] = unit
                
                rep_en = st.checkbox("Configure as custom repeat dosing regimen", key=f"dr_sweep_phg_rep_en_{j}", value=False)
                if rep_en:
                    r_interval = st.number_input("Interval (hours)", min_value=1.0, value=12.0, key=f"dr_sweep_phg_rep_int_{j}")
                    r_count = st.number_input("Number of repeats", min_value=1, value=4, key=f"dr_sweep_phg_rep_count_{j}")
                    r_start = st.number_input("Start time (hours)", min_value=0.0, value=0.0, key=f"dr_sweep_phg_rep_start_{j}")
                    r_route = st.selectbox("Route", ["bolus", "infusion"], key=f"dr_sweep_phg_rep_route_{j}")
                    r_dur = 0.0
                    if r_route == "infusion":
                        r_dur = st.number_input("Duration (hours)", min_value=0.1, value=2.0, key=f"dr_sweep_phg_rep_dur_{j}")
                    swept_repeat_configs[f"phage_{j}"] = {
                        "interval": r_interval,
                        "count": int(r_count),
                        "start": r_start,
                        "route": r_route,
                        "duration": r_dur
                    }

        # Antibiotics
        for j, a in enumerate(antibiotics):
            st.markdown(f"#### Antibiotic {j}: {a['name']}")
            do_sweep = st.checkbox(f"Sweep Antibiotic {j}", key=f"dr_sweep_abx_en_{j}", value=False)
            if do_sweep:
                series_str = st.text_input(
                    "Dose series (comma-separated)",
                    value="0.5, 1.0, 2.0",
                    key=f"dr_sweep_abx_series_{j}"
                )
                swept_inputs[f"abx_{j}"] = series_str
                swept_units[f"abx_{j}"] = "absolute"
                
                rep_en = st.checkbox("Configure as custom repeat dosing regimen", key=f"dr_sweep_abx_rep_en_{j}", value=False)
                if rep_en:
                    r_interval = st.number_input("Interval (hours)", min_value=1.0, value=12.0, key=f"dr_sweep_abx_rep_int_{j}")
                    r_count = st.number_input("Number of repeats", min_value=1, value=4, key=f"dr_sweep_abx_rep_count_{j}")
                    r_start = st.number_input("Start time (hours)", min_value=0.0, value=0.0, key=f"dr_sweep_abx_rep_start_{j}")
                    r_route = st.selectbox("Route", ["bolus", "infusion"], key=f"dr_sweep_abx_rep_route_{j}")
                    r_dur = 0.0
                    if r_route == "infusion":
                        r_dur = st.number_input("Duration (hours)", min_value=0.1, value=2.0, key=f"dr_sweep_abx_rep_dur_{j}")
                    swept_repeat_configs[f"abx_{j}"] = {
                        "interval": r_interval,
                        "count": int(r_count),
                        "start": r_start,
                        "route": r_route,
                        "duration": r_dur
                    }

        st.markdown("<br>", unsafe_allow_html=True)
        run_sweep = st.button("Run Dose-Response Sweep", width="stretch", type="primary")

    with col_run:
        st.markdown("### Sweep Results")
        if run_sweep:
            # Parse vectors
            parsed_vectors = {}
            for k, val_str in swept_inputs.items():
                try:
                    vals = parse_comma_separated_series(val_str)
                    if not vals:
                        st.error(f"Sweep vector for '{k}' cannot be empty.")
                        st.stop()
                    parsed_vectors[k] = vals
                except ValueError:
                    st.error(f"Invalid comma-separated values for '{k}'. Please check formatting.")
                    st.stop()

            if not parsed_vectors:
                st.error("Please select at least one drug/phage to sweep.")
            else:
                # Perform padding
                padded, warnings = pad_vectors(parsed_vectors)
                for w in warnings:
                    st.warning(f"{w}")

                # Determine number of runs M
                first_key = list(padded.keys())[0]
                M = len(padded[first_key])
                
                # Baseline initial_B
                sum_initial_B = sum(s["initial_B"] for s in strains)

                # Warn once if the pre-run decimates the culture (death w/o dormancy).
                _pc_cfg, _pc_B0, _pc_iP, _pc_iS, *_ = build_nominal_config_from_gui()
                warn_if_prerun_collapses(_pc_cfg, _pc_B0, initial_S=_pc_iS)

                # Save original doses
                original_doses = list(st.session_state.get("int_doses", []))

                # When a phage's DOSE is swept, its baseline initial_P inoculum must not
                # leak in: otherwise a swept value of 0 still starts with the configured
                # P0 (default 1e6) and suppresses the bacteria. Zero the swept phages'
                # initial_P for the duration of the sweep — the dose delivers the phage —
                # and restore it afterwards. (Antibiotics have no such baseline: they are
                # delivered purely via dose events.)
                swept_phage_idx = [int(k.split("_")[1]) for k in padded if k.startswith("phage_")]
                original_initial_P = {j: phages[j]["initial_P"] for j in swept_phage_idx if j < len(phages)}
                for j in original_initial_P:
                    phages[j]["initial_P"] = 0.0

                runs_outcomes = []
                trajectories = [] # list of (time, cfu_b, label) — culturable CFU (B+D)
                total_incl_trajectories = [] # (time, B+D+I+H, label)
                active_trajectories = [] # (time, B, label)
                phage_trajectories = [] # list of (time, total_free_phage, label)
                central_phage_trajectories = []  # (time, central/blood conc. Σ Pc/Vc, label)
                od_trajectories = [] # list of (time, od, label) — only if OD/debris enabled
                _od_enabled = st.session_state.get("int_debris_enabled", False)
                _phage_pk_on = phage_pk_enabled(phages)

                # Progress bar
                progress_bar = st.progress(0)
                status_text = st.empty()
                
                try:
                    for k_idx in range(M):
                        status_text.text(f"Running simulation {k_idx+1} of {M}...")
                        
                        # Generate custom doses for this run
                        custom_doses = []
                        
                        # Add non-swept nominal doses
                        for nd in original_doses:
                            is_swept = False
                            if nd["target_type"] == "phage" and f"phage_{nd['target_idx']}" in padded:
                                is_swept = True
                            elif nd["target_type"] == "antibiotic" and f"abx_{nd['target_idx']}" in padded:
                                is_swept = True
                            
                            if not is_swept:
                                custom_doses.append(nd)

                        # Add swept doses
                        swept_label_parts = []
                        for key, vec in padded.items():
                            val = vec[k_idx]
                            
                            # Scale MOI
                            if swept_units.get(key) == "MOI (relative to B(0))":
                                val_absolute = val * sum_initial_B
                                swept_label_parts.append(f"{key}: {val} MOI ({val_absolute:.1e})")
                                val = val_absolute
                            else:
                                swept_label_parts.append(f"{key}: {val:.1e}")

                            # Check if repeat regimen or nominal dose overrides
                            if key in swept_repeat_configs:
                                rcfg = swept_repeat_configs[key]
                                target_type = "phage" if key.startswith("phage") else "antibiotic"
                                target_idx = int(key.split("_")[1])
                                for r in range(rcfg["count"]):
                                    custom_doses.append({
                                        "time": rcfg["start"] + r * rcfg["interval"],
                                        "amount": val,
                                        "target_type": target_type,
                                        "target_idx": target_idx,
                                        "route": rcfg["route"],
                                        "duration": rcfg["duration"]
                                    })
                            else:
                                # Overwrite nominal doses
                                target_type = "phage" if key.startswith("phage") else "antibiotic"
                                target_idx = int(key.split("_")[1])
                                target_nominal_doses = [d for d in original_doses if d["target_type"] == target_type and d["target_idx"] == target_idx]
                                if target_nominal_doses:
                                    for nd in target_nominal_doses:
                                        nd_copy = dict(nd)
                                        nd_copy["amount"] = val
                                        custom_doses.append(nd_copy)
                                else:
                                    # Create single dose event at t=0
                                    custom_doses.append({
                                        "time": 0.0,
                                        "amount": val,
                                        "target_type": target_type,
                                        "target_idx": target_idx,
                                        "route": "bolus",
                                        "duration": 0.0
                                    })

                        st.session_state.int_doses = custom_doses
                        
                        # Run
                        result, config = run_sim_from_gui_params()
                        
                        # Compute metrics. Summary metrics use culturable CFU (B+D); the
                        # other bacterial bases are captured too so the user can overlay them.
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
                        
                        run_label = ", ".join(swept_label_parts)
                        runs_outcomes.append({
                            "Run Index": k_idx + 1,
                            "Swept Doses": run_label,
                            "Nadir (cells/mL)": nadir_val,
                            "AUC (cells·h/mL)": auc_val,
                            "Clearance Time (h)": t_clear if t_clear is not None else np.nan,
                            "2-Log Red Time (h)": t_log_red if t_log_red is not None else np.nan
                        })
                        
                        trajectories.append((result.time, _cfu, f"Run {k_idx + 1}: {run_label}"))
                        total_incl_trajectories.append((result.time, _total_incl, f"Run {k_idx + 1}: {run_label}"))
                        active_trajectories.append((result.time, _active, f"Run {k_idx + 1}: {run_label}"))
                        phage_trajectories.append((result.time, result.sum_prefixes("P"), f"Run {k_idx + 1}: {run_label}"))
                        if _phage_pk_on:
                            central_phage_trajectories.append(
                                (result.time, central_phage_total(result, phages), f"Run {k_idx + 1}: {run_label}"))
                        if _od_enabled:
                            _od = (_safe_od(result, total_bacteria))
                            od_trajectories.append((result.time, _od, f"Run {k_idx + 1}: {run_label}"))
                        progress_bar.progress((k_idx + 1) / M)
                        
                    status_text.text("Sweep completed!")
                finally:
                    # Restore original doses + swept phages' baseline inoculum
                    st.session_state.int_doses = original_doses
                    for j, v in original_initial_P.items():
                        phages[j]["initial_P"] = v

                # Persist results so they survive navigation (rendered below, until
                # the sweep is re-run).
                st.session_state.dr_sweep_result = {
                    "summary": runs_outcomes,
                    "trajectories": [(np.asarray(t), np.asarray(b), lbl) for t, b, lbl in trajectories],
                    "total_incl_trajectories": [(np.asarray(t), np.asarray(b), lbl) for t, b, lbl in total_incl_trajectories],
                    "active_trajectories": [(np.asarray(t), np.asarray(b), lbl) for t, b, lbl in active_trajectories],
                    "phage_trajectories": [(np.asarray(t), np.asarray(p), lbl) for t, p, lbl in phage_trajectories],
                    "central_phage_trajectories": [(np.asarray(t), np.asarray(p), lbl) for t, p, lbl in central_phage_trajectories],
                    "od_trajectories": [(np.asarray(t), np.asarray(o), lbl) for t, o, lbl in od_trajectories],
                }

        # Render the (persisted) sweep result if one exists.
        _dr = st.session_state.get("dr_sweep_result")
        if _dr:
            import plotly.graph_objects as go
            df_summary = pd.DataFrame(_dr["summary"])
            _sweep_summary_tiles(df_summary)
            st.markdown("#### Summary of Runs")
            st.dataframe(
                df_summary.style.format({
                    "Nadir (cells/mL)": "{:.2e}",
                    "AUC (cells·h/mL)": "{:.2e}",
                    "Clearance Time (h)": "{:.1f}",
                    "2-Log Red Time (h)": "{:.1f}"
                }),
                width="stretch"
            )

            # Trajectories across doses — pick which subpopulation/metric to trace.
            _traces = {"CFU — culturable (B+D)": ("log", _dr["trajectories"])}
            if _dr.get("total_incl_trajectories"):
                _traces["Total incl. infected (B+D+I+H)"] = ("log", _dr["total_incl_trajectories"])
            if _dr.get("active_trajectories"):
                _traces["Active only (B)"] = ("log", _dr["active_trajectories"])
            if _dr.get("phage_trajectories") and st.session_state.get("int_phages"):
                _traces["Free phage — infection site (PFU/mL)"] = ("log", _dr["phage_trajectories"])
            if _dr.get("central_phage_trajectories"):
                _traces["Central phage — blood (PFU/mL)"] = ("log", _dr["central_phage_trajectories"])
            if _dr.get("od_trajectories"):
                _traces["Optical density (AU)"] = ("linear", _dr["od_trajectories"])
            plot_sweep_traces(_traces, "dr_traj", title_suffix="across doses")

            st.markdown("#### Outcome Metrics vs Run Index")
            fig_metrics = go.Figure()
            fig_metrics.add_trace(go.Scatter(x=df_summary["Run Index"], y=df_summary["AUC (cells·h/mL)"], mode="lines+markers", name="Bacterial AUC", yaxis="y1"))
            fig_metrics.add_trace(go.Scatter(x=df_summary["Run Index"], y=df_summary["Nadir (cells/mL)"], mode="lines+markers", name="Nadir", yaxis="y2"))
            fig_metrics.update_layout(
                title="Bacterial Efficacy Metrics Across Sweep Runs",
                xaxis=dict(title="Run Index"),
                yaxis=dict(title="AUC (cells·h/mL)", type="log"),
                yaxis2=dict(title="Nadir (cells/mL)", type="log", overlaying="y", side="right"),
                template="plotly_white" if theme_mode == "Light" else "plotly_dark")
            st.plotly_chart(fig_metrics, width="stretch")
        else:
            st.info("Configure the sweep on the left and click **Run Dose-Response Sweep** to view results.")

    with st.expander("View Python Reproduction Code"):
        st.caption("Standalone script that reproduces this sweep — the recorded base "
                   "model plus a loop rebuilding the per-run dose schedule exactly as the app does.")
        try:
            st.code(generate_dose_sweep_reproduction_code(), language="python")
        except Exception as _e:
            st.warning(f"Reproduction code unavailable: {_e}")

    # Persist the sweep controls so they survive navigation (see reseed above).
    save_widget_config("dr_sweep_config", ("dr_sweep_",))

# ── Parameter Sweeps Page ──────────────────────────────────────────────────────
