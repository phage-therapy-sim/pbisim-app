"""Rendered by app.py when this page is selected."""
from pbisim_app.common import *  # noqa: F401,F403


def render():
    theme_mode = st.session_state.get("theme_mode", "Light")
    st.title("Clinical Trials & Cohort Simulator")
    st.caption("Generate a virtual population (VPOP), apply statistical variability (IIV), and run matching parallel arms.")
    # The trial's nominal patient = a CHOSEN Model (frozen snapshot or live draft).
    _sel = page_model_selector("trial")
    with model_config_context(resolve_model_snapshot(_sel)):
        _render_body(theme_mode)


def _render_body(theme_mode):
    st.markdown(
        "<div class='info-banner'>Virtual Cohort simulations use the <b>Model selected above</b> as the "
        "baseline 'nominal patient'; IIV then varies it per patient. Arms/dosing are set below.</div>",
        unsafe_allow_html=True,
    )
    
    t_cols = st.columns([1, 2])
    
    with t_cols[0]:
        st.markdown("### Trial Settings")
        trial_patients = st.number_input("Cohort Size (N)", min_value=10, max_value=500, value=50, step=10)
        trial_seed = st.number_input("Cohort RNG Seed", value=42)
        trial_t_end = st.number_input("Trial Duration (hours)", min_value=12.0, max_value=336.0, value=72.0, step=12.0)
        trial_dt = st.number_input("Solver output step (dt)", min_value=0.05, max_value=1.0, value=0.25, step=0.05)
        trial_n_jobs = st.slider("Parallel workers (n_jobs)", min_value=1, max_value=16, value=1, help="Parallel patient simulation via joblib (loky). Keep at 1 on small/shared hosts (e.g. the free Render tier) — forked worker processes can segfault or OOM there; raise it on a beefier machine.")
        
        st.markdown("### Treatment Arms")
        st.caption("Define any number of arms (e.g. low-dose vs high-dose), each with its own phage / antibiotic regimen. The Control arm never receives doses.")
        trial_include_control = st.checkbox("Include Control arm (no doses)", value=True)

        _tphages = st.session_state.get("int_phages", [])
        _tabx = st.session_state.get("int_antibiotics", [])
        trial_arms = st.session_state.get("trial_arms", [])

        if not _tphages and not _tabx:
            st.info("Configure at least one phage or antibiotic in the Interactive Simulator to define dosed arms.")

        # Existing arms — editable in place
        for _ai, _arm in enumerate(list(trial_arms)):
            _arm.setdefault("_id", _next_uid("arm"))   # stable key across reorder/delete
            _aid = _arm["_id"]
            _lc, _dc = st.columns([6, 1])
            with _lc:
                st.markdown(f"**{_arm['name']}** — {arm_regimen_summary(_arm)}")
            with _dc:
                if st.button(":material/delete:", key=f"del_arm_{_aid}"):
                    trial_arms.pop(_ai)
                    st.session_state.trial_arms = trial_arms
                    st.rerun()
            with st.expander(f"Edit '{_arm['name']}'", expanded=False):
                _en = st.text_input("Arm name", value=_arm["name"], key=f"edit_arm_name_{_aid}")
                _ep, _ea = {"on": False}, {"on": False}
                if _tphages:
                    st.markdown("**Phage dosing**")
                    _ep = render_regimen_config(f"edit_arm_p_{_aid}", _tphages, "phage",
                                                1e9, "Amount (PFU)", initial=_arm.get("phage"))
                if _tabx:
                    st.markdown("**Antibiotic dosing**")
                    _ea = render_regimen_config(f"edit_arm_a_{_aid}", _tabx, "antibiotic",
                                                10.0, "Amount (mg)", initial=_arm.get("abx"))
                if st.button("Save changes", key=f"save_arm_{_aid}"):
                    _arm["name"], _arm["phage"], _arm["abx"] = _en, _ep, _ea
                    st.session_state.trial_arms = trial_arms
                    st.rerun()

        # Add-arm form
        with st.expander("+ Add treatment arm", expanded=not trial_arms):
            _new_name = st.text_input("Arm name", value=f"Arm {len(trial_arms) + 1}", key="new_arm_name")
            _pcfg, _acfg = {"on": False}, {"on": False}
            if _tphages:
                st.markdown("**Phage dosing**")
                _pcfg = render_regimen_config("new_arm_p", _tphages, "phage", 1e9, "Amount (PFU)", default_on=True)
            if _tabx:
                st.markdown("**Antibiotic dosing**")
                _acfg = render_regimen_config("new_arm_a", _tabx, "antibiotic", 10.0, "Amount (mg)", default_on=not _tphages)
            if st.button("+ Add arm", key="add_arm_btn"):
                trial_arms.append({"name": _new_name, "phage": _pcfg, "abx": _acfg})
                st.session_state.trial_arms = trial_arms
                st.rerun()

        st.markdown("### Parameter Variability (IIV)")
        
        # Active IIVs — editable in place
        trial_iivs = st.session_state.get("trial_iiv_inputs", [])

        for idx, iiv in enumerate(trial_iivs):
            iiv.setdefault("_id", _next_uid("iiv"))   # stable key across reorder/delete
            _iid = iiv["_id"]
            _pname = next((n for n, p in IIV_PARAMETERS.items() if p == iiv["path"]), iiv["path"])
            col_p, col_act = st.columns([6, 1])
            with col_p:
                st.markdown(f"**{_pname}** — {iiv['dist_type']} {iiv['params']} ({iiv['mode']})")
            with col_act:
                if st.button(":material/delete:", key=f"del_iiv_{_iid}"):
                    trial_iivs.pop(idx)
                    st.session_state.trial_iiv_inputs = trial_iivs
                    st.rerun()
            with st.expander("Edit", expanded=False):
                _edited = render_iiv_config(f"edit_iiv_{_iid}", initial=iiv)
                if st.button("Save changes", key=f"save_iiv_{_iid}"):
                    _edited["_id"] = _iid
                    trial_iivs[idx] = _edited
                    st.session_state.trial_iiv_inputs = trial_iivs
                    st.rerun()

        # Add IIV form
        with st.expander("+ Add Parameter Variability"):
            _new_iiv = render_iiv_config("new_iiv")
            if st.button("Add Parameter IIV"):
                trial_iivs.append(_new_iiv)
                st.session_state.trial_iiv_inputs = trial_iivs
                st.success("Added parameter variability.")
                st.rerun()

        # Run Button
        st.markdown("<br>", unsafe_allow_html=True)
        if st.button("Run Parallel Clinical Trial", width="stretch", type="primary"):
            with st.spinner("Generating cohort populations & simulating treatment arms..."):
                try:
                    # 1. Compile nominal base config
                    base_cfg, init_B, init_P, init_S, model_kwargs = build_nominal_config_from_gui()

                    # Phage is the *intervention*, delivered per-arm via each arm's
                    # dose schedule (Treatment Arms builder). In this crossover design
                    # all arms share initial_conditions and differ only by dose_schedule,
                    # so every arm starts with zero free phage and receives only its own
                    # configured doses — keeping the Control arm a true no-treatment arm.
                    init_P = np.asarray(init_P, dtype=float)
                    base_P = np.zeros_like(init_P)

                    from pbisim.trial.population import InitialConditions
                    base_cfg.initial_conditions = InitialConditions(
                        B=init_B,
                        P=base_P,
                        S=init_S,
                        D=model_kwargs.get("initial_D", None),
                        Imm=model_kwargs.get("initial_Imm", None),
                    )

                    # 2. Assemble arms from the Treatment Arms builder
                    arms = []
                    _used_names = set()

                    def _add_arm(name, doses):
                        # ClinicalTrial requires unique arm names
                        base = name.strip() or "Arm"
                        uniq, k = base, 2
                        while uniq in _used_names:
                            uniq = f"{base} ({k})"; k += 1
                        _used_names.add(uniq)
                        arms.append(TreatmentArm(name=uniq, dose_schedule=DoseSchedule(list(doses))))

                    if trial_include_control:
                        _add_arm("Control", [])

                    for _arm in trial_arms:
                        _doses = arm_dose_events(_arm)
                        if not _doses:
                            st.warning(f"Arm '{_arm['name']}' has no doses — it will behave like the Control arm.")
                        _add_arm(_arm["name"], _doses)

                    if not arms:
                        st.error("Add at least one treatment arm (or enable the Control arm) to run.")
                    else:
                        pretreatment_hours = st.session_state.get("int_t_prerun", 0.0)
                        
                        trial_result = run_trial_simulation(
                            base_cfg,
                            trial_iivs,
                            arms,
                            n_patients=int(trial_patients),
                            t_end=trial_t_end,
                            dt=trial_dt,
                            seed=int(trial_seed),
                            pretreatment_hours=pretreatment_hours,
                            n_jobs=int(trial_n_jobs),
                            base_initial_B=init_B,
                            base_initial_P=base_P,
                            base_initial_S=init_S,
                            inherit_debris=st.session_state.get("int_prerun_inherit_debris", True),
                            **model_kwargs
                        )
                        st.session_state.trial_result = trial_result
                        st.success("Clinical Trial cohort simulation completed successfully!")
                except Exception as e:
                    st.error(f"Trial Execution Error: {e}")
                    import traceback
                    st.code(traceback.format_exc())
                    
    with t_cols[1]:
        st.markdown("### Outcomes & Visualization")
        
        if st.session_state.trial_result is None:
            st.info("Run the clinical trial simulation on the left panel to display outcomes.")
        else:
            result = st.session_state.trial_result
            
            # Outcome selection
            _endpoint_labels = {
                "tte": "Time-to-Eradication (TTE)",
                "tt2lr": "Time-to-2-Log-Reduction (TT2LR)",
            }
            _metric_labels = {
                "max_log_reduction": "Maximum Log Reduction",
                "log_reduction_final": "Log Reduction (baseline → last obs)",
                "bacterial_auc": "Bacterial AUC",
                "nadir_count": "Nadir Count",
            }
            c_v1, c_v2 = st.columns(2)
            with c_v1:
                endpoint_choice = st.selectbox(
                    "Survival Endpoint (time-to-event)", ["tte", "tt2lr"], index=0,
                    format_func=lambda x: _endpoint_labels[x],
                )
            with c_v2:
                metric_choice = st.selectbox(
                    "Distribution Metric", list(_metric_labels), index=0,
                    format_func=lambda x: _metric_labels[x],
                )
                
            clearance_threshold = st.session_state.get("int_extinction_threshold", 100.0)

            # Cure-rate summary tiles (one per arm; eradication = reached clearance by t_end)
            try:
                _arm_names = list(result.arm_names)
            except Exception:
                _arm_names = []
            if _arm_names and len(_arm_names) <= 6:
                _tiles = st.columns(len(_arm_names))
                for _col, _arm in zip(_tiles, _arm_names):
                    try:
                        _pats = [r for r in result[_arm].results if r is not None]
                        _tt = [time_to_clearance(r, threshold=clearance_threshold) for r in _pats]
                        _cured = [t for t in _tt if t is not None]
                        _rate = (len(_cured) / len(_pats) * 100.0) if _pats else 0.0
                        _median = float(np.median(_cured)) if _cured else None
                        _sub = (f"{len(_cured)}/{len(_pats)} cured · median {_median:.0f} h"
                                if _median is not None else f"{len(_cured)}/{len(_pats)} cured")
                    except Exception:
                        _rate, _sub = 0.0, "n/a"
                    _col.markdown(
                        f"""
                        <div class="metric-container">
                            <div class="metric-label">Cure rate · {_arm}</div>
                            <div class="metric-value">{_rate:.0f}%</div>
                            <div class="metric-sub">{_sub}</div>
                        </div>
                        """,
                        unsafe_allow_html=True,
                    )
                st.markdown("<br>", unsafe_allow_html=True)

            # Raw PKPD time trajectories per arm — pick which outputs to show.
            st.markdown("#### PK/PD trajectories (median & IQR per arm)")
            _TRIAL_OUTPUTS = {
                "CFU — culturable (B+D)": (("B", "D"), "log₁₀ CFU/mL"),
                "Total incl. infected (B+D+I+H)": (("B", "D", "I", "H"), "log₁₀ cells/mL"),
                "Active only (B)": (("B",), "log₁₀ cells/mL"),
                "Free phage — infection site (PFU/mL)": (("P",), "log₁₀ PFU/mL"),
                "Immune effector": (("Imm",), "log₁₀ Imm"),
            }
            # Central-compartment (blood) phage is only present when phage PK is enabled.
            # Show it as a CONCENTRATION = Σ Pc / Vc (divide the summed central amount by
            # the central volume; uses the first PK phage's Vc — see caption if they differ).
            _pk_on = phage_pk_enabled(_tphages)
            _pk_vcs = [float(p.get("Vc", 5000.0)) for p in _tphages if p.get("pk_mode", "None") != "None"]
            _central_vc = _pk_vcs[0] if _pk_vcs else 1.0
            _central_name = "Central phage — blood (PFU/mL)"
            if _pk_on:
                _TRIAL_OUTPUTS[_central_name] = (("Pc",), "log₁₀ PFU/mL")
            if not _tphages:   # antibiotic-only model: no phage panel to offer
                _TRIAL_OUTPUTS.pop("Free phage — infection site (PFU/mL)")
            _outputs = st.multiselect(
                "Outputs", list(_TRIAL_OUTPUTS),
                default=["CFU — culturable (B+D)"]
                        + (["Free phage — infection site (PFU/mL)"] if _tphages else []),
                key="trial_pkpd_outputs",
            )
            if _pk_on and _central_name in _outputs and len(set(_pk_vcs)) > 1:
                st.caption("Central phage uses the first PK phage's Vc for the Σ Pc → concentration; "
                           "phages have differing Vc, so it's approximate for multi-phage cocktails.")
            for _name in _outputs:
                _prefixes, _ylab = _TRIAL_OUTPUTS[_name]
                _div = _central_vc if _prefixes == ("Pc",) else 1.0
                st.plotly_chart(
                    plot_pkpd_trajectories_plotly(result, prefixes=_prefixes, title=_name,
                                                  y_label=_ylab, divide_by=_div),
                    width="stretch",
                )

            # Step survival plot
            st.markdown("#### Step-Survival (Kaplan-Meier)")
            fig_km = plot_kaplan_meier_plotly(result, endpoint=endpoint_choice, t_end=trial_t_end, threshold=clearance_threshold, n_logs=2.0)
            st.plotly_chart(fig_km, width="stretch")

            # Metric distributions
            st.markdown("#### Distribution of outcomes")
            fig_dist = plot_metric_distributions_plotly(result, metric=metric_choice)
            st.plotly_chart(fig_dist, width="stretch")
            
            # Data Exports
            st.markdown("---")
            st.markdown("### Cohort Data Exports")
            
            cx1, cx2 = st.columns(2)
            with cx1:
                # Outcome Dataframe
                out_df = result.outcome_dataframe(endpoint=endpoint_choice, t_end=trial_t_end, threshold=clearance_threshold)
                csv_out = out_df.to_csv(index=False)
                st.download_button(
                    "Download Survival Outcomes DataFrame (CSV)",
                    data=csv_out,
                    file_name="pbisim_survival_outcomes.csv",
                    mime="text/csv",
                    width="stretch"
                )
            with cx2:
                # NLME Dataframe for pharmacometric models
                outputs_spec = {"DV_B": ("B", "D", "I", "H")}
                obs_times = np.linspace(0, trial_t_end, 10)
                nlme_df = result.nlme_dataframe(outputs_spec, times=obs_times)
                csv_nlme = nlme_df.to_csv(index=False)
                st.download_button(
                    "Download Pharmacometrics (NLME) DataFrame (CSV)",
                    data=csv_nlme,
                    file_name="pbisim_nlme_cohort.csv",
                    mime="text/csv",
                    width="stretch"
                )


# ── Dose-Response Sweeps Page ──────────────────────────────────────────────────
