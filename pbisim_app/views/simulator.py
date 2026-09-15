"""Rendered by app.py when this page is selected."""
from pbisim_app.common import *  # noqa: F401,F403


def render_model_builder(inoculum_mode="magnitude"):
    """Full model builder — builder-mode selector, growth/death signal
    functions, and the selected Direct / Binary-Genotypes / Custom-Strains body.
    Extracted from the Interactive Simulator's Tab 1 so the Calibration
    manual-tuning panel renders the exact same controls (no drift). Reads/writes
    the shared int_* session state and returns the active builder_mode string.

    ``inoculum_mode``:
      * ``"magnitude"`` (Simulator) — per-strain B₀/P₀ are absolute densities.
      * ``"ratio"`` (Calibration) — per-strain B₀ inputs express a *relative*
        ratio (the absolute total comes from the per-arm B₀ in the Calibration
        conditions), the absolute inoculum totals (BRG equilibrium total, phage
        P₀) are hidden, but the BRG equilibrium-IC checkbox is KEPT (it sets the
        resistant fraction from fitness cost + mutation).
    """
    _ratio_inoc = (inoculum_mode == "ratio")
    strains = st.session_state.get("int_strains", [])
    phages = st.session_state.get("int_phages", [])
    antibiotics = st.session_state.get("int_antibiotics", [])
    # Builder Mode Selector
    builder_mode = st.selectbox(
        "Bacterial Population Builder Mode",
        ["Direct (ModelBuilder)", "Binary Genotypes (BRG)", "Custom Strains & Graph (StrainSet)"],
        index=["Direct (ModelBuilder)", "Binary Genotypes (BRG)", "Custom Strains & Graph (StrainSet)"].index(
            st.session_state.get("int_builder_mode", "Direct (ModelBuilder)")
        ),
        key="widget_builder_mode"
    )
    if builder_mode != st.session_state.get("int_builder_mode", "Direct (ModelBuilder)"):
        st.session_state["int_builder_mode"] = builder_mode
        st.session_state.simulation_result = None
        st.session_state.simulation_config = None
        st.rerun()

    # ── Growth model (signal function + its half-saturation / carrying capacity) ──
    # Kept here in the model builder (not the Environment tab) since it defines the
    # growth kinetics. The nutrient environment (S0 / recycle / inflow / washout)
    # stays in Environment & Dosing.
    _gm1, _gm2 = st.columns([2, 1])
    with _gm1:
        _gs_cur_fn = st.session_state.get("int_growth_function", "monod_growth")
        _gs_labels = list(GROWTH_SIGNALS.keys())
        _gs_cur_label = next((L for L, (fn, _) in GROWTH_SIGNALS.items() if fn == _gs_cur_fn), _gs_labels[0])
        _gs_choice = st.selectbox(
            "Growth signal function", _gs_labels, index=_gs_labels.index(_gs_cur_label),
            help="How the per-strain growth rate is modulated:  constant = unlimited;  "
                 "nutrient = Monod S/(Ks+S);  density = logistic (1−ΣB/K);  "
                 "nutrient+density = Monod × logistic;  density-throttled = Monod × a soft "
                 "hyperbolic throttle 1/(1+ΣB/Kd) (never fully stops);  Gompertz = "
                 "exp(−exp(−k(S−S∞)));  sequential/diauxic = one nutrient pool consumed "
                 "through ordered Monod phases (set the phases below);  smooth two-efficiency "
                 "Monod = K blends from efficient (high S) to inefficient (low S) — a smooth, "
                 "better-fitting diauxie.",
        )
    _gs_fn, _gs_track = GROWTH_SIGNALS[_gs_choice]
    st.session_state["int_growth_function"] = _gs_fn
    st.session_state["int_track_nutrients"] = _gs_track
    with _gm2:
        # Monod Ks (nutrient half-saturation). For smooth_efficiency_monod this is the
        # EFFICIENT (high-S) K; the low-S K + transition are set in the block below.
        if _gs_fn in ("monod_growth", "monod_logistic_growth", "density_throttled_growth",
                      "smooth_efficiency_monod"):
            _is_smeff = _gs_fn == "smooth_efficiency_monod"
            st.session_state["int_monod_constant"] = st.number_input(
                "Efficient Monod K (high S)" if _is_smeff else "Monod constant (Ks)",
                value=float(st.session_state.get("int_monod_constant", 0.3)), step=0.05, format="%g",
                help=("High-substrate (efficient) Monod K — the small-K end." if _is_smeff
                      else "Nutrient half-saturation for Monod growth S/(Ks+S)."))
        # Carrying capacity K (logistic density ceiling).
        if _gs_fn in ("logistic_growth", "monod_logistic_growth"):
            st.session_state["int_carrying_capacity"] = st.number_input(
                "Carrying capacity K (CFU·mL⁻¹)",
                value=float(st.session_state.get("int_carrying_capacity", 1e9)),
                format="%.1e", help="Density ceiling for logistic growth (1 − ΣB/K).")
        # Density-throttle constant Kd (soft, hyperbolic — 0.5 at ΣB=Kd, never 0).
        if _gs_fn == "density_throttled_growth":
            st.session_state["int_density_growth_constant"] = st.number_input(
                "Density throttle Kd (CFU·mL⁻¹)",
                value=float(st.session_state.get("int_density_growth_constant", 1e9)),
                format="%.1e", help="Hyperbolic density throttle 1/(1+ΣB/Kd): 0.5 at ΣB=Kd, "
                                    "never fully stops (growth creeps past Kd toward the nutrient plateau).")
        # Gompertz nutrient-response sigmoid exp(−exp(−k(S−S∞))). k and S∞ are on the
        # NUTRIENT scale (S₀ default 1.0) — dedicated keys so S∞ doesn't inherit the
        # density-scale K (~1e9), which would overflow the exponent → flat 0-growth curve.
        if _gs_fn == "gompertz_growth":
            st.session_state["int_gompertz_k"] = st.number_input(
                "Gompertz shape k", min_value=0.0,
                value=float(st.session_state.get("int_gompertz_k", 10.0)), step=0.5, format="%g",
                help="Curvature of exp(−exp(−k(S−S∞))). Larger k = sharper nutrient response.")
            st.session_state["int_gompertz_sinf"] = st.number_input(
                "Gompertz inflection S∞ (nutrient units)",
                value=float(st.session_state.get("int_gompertz_sinf", 0.5)), step=0.05, format="%g",
                help="Nutrient level at the growth-rate inflection — on the S scale (S₀ default "
                     "1.0), NOT the density scale. A value like 1e9 makes growth ≈ 0 (flat curve).")
        # Smooth two-efficiency Monod: K blends from the efficient high-S K (above) to an
        # inefficient low-S K, via a Hill of S — a differentiable diauxie (better-conditioned
        # fit). K_low should be LARGER than the efficient K (poorer uptake when starved).
        if _gs_fn == "smooth_efficiency_monod":
            st.session_state["int_monod_K_low"] = st.number_input(
                "Inefficient Monod K (low S)", min_value=0.0,
                value=float(st.session_state.get("int_monod_K_low", 3.0)), step=0.5, format="%g",
                help="Low-substrate (inefficient) Monod K — typically larger than the efficient "
                     "K, so growth decelerates into a slow non-saturating creep as S depletes.")
            st.session_state["int_monod_efficiency_theta"] = st.number_input(
                "Efficiency transition θ (nutrient units)", min_value=0.0, max_value=1.0,
                value=float(st.session_state.get("int_monod_efficiency_theta", 0.5)),
                step=0.05, format="%g",
                help="Nutrient level where K is halfway between the efficient and inefficient "
                     "values (on the S scale, in (0,1)).")
            st.session_state["int_monod_efficiency_hill"] = st.number_input(
                "Efficiency transition sharpness (Hill)", min_value=1.0,
                value=float(st.session_state.get("int_monod_efficiency_hill", 4.0)),
                step=1.0, format="%g", help="Hill exponent for the K(S) transition; ≥ 1, larger = sharper.")

    # Diauxic (sequential_monod) per-phase editor — a single nutrient pool consumed
    # through X ordered phases, each with its own max-rate factor and Monod K. As S
    # drops below each threshold θ the next phase takes over (biphasic growth curves).
    if _gs_fn == "sequential_monod":
        with st.expander("Diauxic growth phases", expanded=True):
            st.caption(
                "One nutrient pool, X ordered Monod phases. Phase 1 is the reference "
                "(rate factor 1.0); each later phase scales the strain growth rate and "
                "has its own Monod K. Thresholds θ (strictly decreasing, in (0,1)) set "
                "the nutrient level at which each next phase begins.")
            _nph = int(st.number_input(
                "Number of phases (X)", min_value=2, max_value=4,
                value=int(st.session_state.get("int_growth_n_phases", 2) or 2), step=1,
                key="widget_growth_n_phases",
                help="Biphasic (2) covers most diauxie; up to 4 supported."))
            st.session_state["int_growth_n_phases"] = _nph
            st.markdown("**Per-phase max-rate factor & Monod K**")
            for _i in range(_nph):
                _pc1, _pc2 = st.columns(2)
                with _pc1:
                    if _i == 0:
                        st.number_input("Phase 1 rate factor", value=1.0, disabled=True,
                                        key="gp_rate_0_display",
                                        help="Phase 1 is the reference and is pinned to 1.0.")
                    else:
                        st.session_state[f"gp_rate_{_i}"] = st.number_input(
                            f"Phase {_i+1} rate factor", min_value=0.0,
                            value=float(st.session_state.get(f"gp_rate_{_i}", 0.5)),
                            step=0.05, format="%g", key=f"widget_gp_rate_{_i}")
                with _pc2:
                    _mdef = 0.3 if _i else float(st.session_state.get("int_monod_constant", 0.3))
                    st.session_state[f"gp_monod_{_i}"] = st.number_input(
                        f"Phase {_i+1} Monod K", min_value=0.0,
                        value=float(st.session_state.get(f"gp_monod_{_i}", _mdef)),
                        step=0.05, format="%g", key=f"widget_gp_monod_{_i}")
            if _nph > 1:
                st.markdown("**Phase thresholds θ (strictly decreasing, in (0,1))**")
                _tcols = st.columns(_nph - 1)
                for _i in range(_nph - 1):
                    with _tcols[_i]:
                        st.session_state[f"gp_thresh_{_i}"] = st.number_input(
                            f"θ{_i+1} (→ phase {_i+2})", min_value=0.0, max_value=1.0,
                            value=float(st.session_state.get(f"gp_thresh_{_i}", 0.3 / (_i + 2))),
                            step=0.01, format="%g", key=f"widget_gp_thresh_{_i}")
            _seq_err = validate_sequential_growth(*sequential_growth_phase_params())
            if _seq_err:
                st.error(f"Diauxic phases invalid: {_seq_err}")

    # Death signal (model-wide) — modulates the per-strain natural death rate dB.
    _dth_cur_fn = st.session_state.get("int_death_function", "constant_death")
    _dth_labels = list(DEATH_SIGNALS.keys())
    _dth_cur_label = next((L for L, fn in DEATH_SIGNALS.items() if fn == _dth_cur_fn), _dth_labels[0])
    _dth_choice = st.selectbox(
        "Death signal function", _dth_labels, index=_dth_labels.index(_dth_cur_label),
        help="How the per-strain natural death rate dB is modulated:  constant = flat rate d "
             "(the previous behaviour);  nutrient = starvation d·(1−S/(Ks+S)) (rises as nutrients "
             "deplete);  density = crowding d·min(1, ΣB/K). Note: starvation death only separates "
             "stationary from death phase when nutrients persist at a low plateau (recycling).")
    st.session_state["int_death_function"] = DEATH_SIGNALS[_dth_choice]

    # Lysis signal (model-wide) — how fast the infected→lysis (latent) chain progresses.
    _ly1, _ly2 = st.columns([2, 1])
    with _ly1:
        _ly_cur_fn = st.session_state.get("int_lysis_function", "constant_lysis")
        _ly_labels = list(LYSIS_SIGNALS.keys())
        _ly_cur_label = next((L for L, fn in LYSIS_SIGNALS.items() if fn == _ly_cur_fn), _ly_labels[0])
        _ly_choice = st.selectbox(
            "Lysis signal function", _ly_labels, index=_ly_labels.index(_ly_cur_label),
            help="How the infected→lysis (latent-chain) progression is modulated:  constant = "
                 "fixed rate (phi_lysis = 1, the previous behaviour);  nutrient (Monod) = "
                 "frac_lysis, phi_lysis = S/(Ks+S), which slows lysis when nutrients are scarce "
                 "(needs a nutrient-tracking growth signal).")
    _ly_fn = LYSIS_SIGNALS[_ly_choice]
    st.session_state["int_lysis_function"] = _ly_fn
    with _ly2:
        if _ly_fn == "frac_lysis":
            if growth_tracks_nutrients():   # any nutrient-tracking growth, not just Monod
                st.session_state["int_monod_constant_lysis"] = st.number_input(
                    "Lysis Ks (Ks_lysis)",
                    value=float(st.session_state.get("int_monod_constant_lysis", 0.3)),
                    step=0.05, format="%g",
                    help="Nutrient half-saturation for phi_lysis = S/(Ks_lysis + S). "
                         "Lower = lysis more sensitive to nutrient depletion.")
                st.session_state["int_lysis_floor"] = st.number_input(
                    "Lysis floor (φ_min)", min_value=0.0, max_value=1.0,
                    value=float(st.session_state.get("int_lysis_floor", 0.0) or 0.0),
                    step=0.05, format="%g",
                    help="Residual lysis efficacy retained as nutrients vanish: the signal is "
                         "rescaled φ = φ_min + (1−φ_min)·signal, so φ→1 at full growth and "
                         "φ→φ_min as S→0. 0 = efficacy vanishes at stationary phase (legacy); "
                         ">0 keeps adsorbed phage lysing slow/stationary cells at a reduced rate; "
                         "1 = constant lysis.")
    if _ly_fn == "frac_lysis" and not growth_tracks_nutrients():
        st.caption("⚠ Nutrient-coupled lysis needs a nutrient-tracking growth signal "
                   "(Monod, density-throttled, Gompertz, diauxic, …) — it will fall back to "
                   "constant lysis until you choose one.")

    # Dormancy signal functions (model-wide) — the engine's dormancy_function /
    # dormancy_monod_constant / dormancy_carrying_capacity are SINGLE fields for the whole
    # model, so these apply to every dormancy-enabled strain (per-strain RATES stay on each
    # strain card). Shown always, like the growth/death/lysis selectors above.
    _dm = st.columns(3)
    with _dm[0]:
        st.session_state["int_dormancy_signal"] = st.selectbox(
            "Dormancy signal function", SIGNAL_OPTIONS,
            index=SIGNAL_OPTIONS.index(canonical_signal(st.session_state.get("int_dormancy_signal", "nutrient"))),
            help="How the dormancy-ENTRY rate is modulated (applies to all strains with dormancy "
                 "enabled):  constant;  nutrient = Monod S/(Ks+S);  density = ΣB/K_dorm;  nutrient+density.")
    with _dm[1]:
        st.session_state["int_resuscitation_signal"] = st.selectbox(
            "Resuscitation signal function", SIGNAL_OPTIONS,
            index=SIGNAL_OPTIONS.index(canonical_signal(st.session_state.get("int_resuscitation_signal", "nutrient"))),
            help="How the wake-up (resuscitation) rate is modulated. Same four options.")
    with _dm[2]:
        st.session_state["int_diffusion_signal"] = st.selectbox(
            "Depth-diffusion signal function", SIGNAL_OPTIONS,
            index=SIGNAL_OPTIONS.index(canonical_signal(st.session_state.get("int_diffusion_signal", "constant"))),
            help="How cells ratchet between dormancy depth layers. constant = legacy symmetric; "
                 "nutrient/density drive cells deeper under prolonged stress.")
    _ds_now = canonical_signal(st.session_state["int_dormancy_signal"])
    _rs_now = canonical_signal(st.session_state["int_resuscitation_signal"])
    _dfs_now = canonical_signal(st.session_state["int_diffusion_signal"])
    _dm_needs_ks = any(s in ("nutrient", "nutrient+density") for s in (_ds_now, _rs_now, _dfs_now))
    _dm_needs_k = any(s in ("density", "nutrient+density") for s in (_ds_now, _rs_now, _dfs_now))
    if _dm_needs_ks or _dm_needs_k:
        _dmk = st.columns(2)
        if _dm_needs_ks:
            with _dmk[0]:
                st.session_state["int_dormancy_monod_constant"] = st.number_input(
                    "Dormancy nutrient half-saturation (Ks)",
                    value=float(st.session_state.get("int_dormancy_monod_constant", 0.0) or 0.0),
                    min_value=0.0, format="%g",
                    help="Half-saturation for the nutrient dormancy signal. 0 = inherit the growth Monod constant.")
        if _dm_needs_k:
            with _dmk[1]:
                st.session_state["int_dormancy_carrying_capacity"] = st.number_input(
                    "Dormancy density threshold (CFU·mL⁻¹)",
                    value=float(st.session_state.get("int_dormancy_carrying_capacity", 1e8) or 1e8),
                    min_value=0.0, format="%.2e",
                    help="Density threshold for the density dormancy signal (rate ∝ ΣB/K_dorm). "
                         "0 = inherit the growth carrying capacity.")

    # What "density" counts for ALL density signals (dormancy / resuscitation / death).
    st.session_state["int_density_total_cells"] = st.checkbox(
        "Density signals count all cell states (B + I + D + H)",
        value=bool(st.session_state.get("int_density_total_cells", False)),
        key="widget_density_total_cells",
        help="Off (default): density-dependent signals use ACTIVE bacteria only, "
             "min(1, ΣB/K). On: they use TOTAL cell density including infected (I), "
             "dormant (D) and hibernating (H) cells — so all compartments count toward "
             "crowding. Applies to every density / nutrient+density dormancy, "
             "resuscitation and death signal.")

    st.markdown("---")
    col1, col2 = st.columns(2)

    # ── DIRECT MODE ──
    if builder_mode == "Direct (ModelBuilder)":
        with col1:
            st.markdown("### Bacterial Strains")

            n_strains = counted_number_input("Number of strains", len(strains), "direct_n_strains", min_value=1, max_value=10)
            if n_strains != len(strains):
                st.session_state.simulation_result = None
                st.session_state.simulation_config = None
            # Adjust list size
            if n_strains > len(strains):
                for i in range(len(strains), n_strains):
                    strains.append(
                        {
                            "name": f"Strain {i}",
                            "initial_B": 1e7 if i == 0 else 0.0,
                            "initial_D": 0.0,
                            "growth_rate": 1.2,
                            "bacteria_to_resource_ratio": 1e9,
                            "death_rate_B": 0.0,
                            "death_rate_D": 0.0,
                            "dormancy_enabled": False,
                            "dormancy_depth": 1,
                            "dormancy_rate": 0.001,
                            "resuscitation_rate": 0.1,
                            "dormancy_diffusion_rate": 0.05,
                            "dormancy_signal": "nutrient",
                            "resuscitation_signal": "nutrient",
                        }
                    )
            elif n_strains < len(strains):
                strains = strains[:n_strains]
            st.session_state["int_strains"] = strains

            # Inputs per strain
            for i in range(n_strains):
                with st.expander(f"Strain {i}: {strains[i]['name']}", expanded=True):
                    strains[i]["name"] = st.text_input(
                        "Name", value=strains[i]["name"], key=f"str_name_{i}"
                    )
                    cc1, cc2 = st.columns(2)
                    with cc1:
                        strains[i]["initial_B"] = st.number_input(
                            "Initial ratio (relative)" if _ratio_inoc else "Initial Density (B0)",
                            value=float(strains[i]["initial_B"]),
                            format="%g" if _ratio_inoc else "%.1e",
                            key=f"str_init_{i}",
                            help=("Relative starting amount across strains; the absolute total "
                                  "comes from the per-arm B₀ in the Calibration conditions."
                                  if _ratio_inoc else None),
                        )
                    with cc2:
                        strains[i]["growth_rate"] = st.number_input(
                            "Growth rate (h⁻¹)",
                            value=float(strains[i]["growth_rate"]),
                            step=0.1,
                            key=f"str_growth_{i}",
                        )
                    strains[i]["bacteria_to_resource_ratio"] = st.number_input(
                        "Bacteria-to-resource ratio",
                        value=float(strains[i].get("bacteria_to_resource_ratio", 1e9)),
                        format="%.2e",
                        key=f"str_ratio_{i}",
                        help="Bacteria produced per unit resource consumed (yield). Governs how fast "
                             "growth depletes the substrate under nutrient-limited (Monod) growth.",
                    )
                    strains[i]["death_rate_B"] = st.number_input(
                        "Natural death rate (h⁻¹)",
                        value=float(strains[i].get("death_rate_B", 0.0)),
                        step=0.01,
                        key=f"str_death_{i}",
                    )
                    strains[i]["dormancy_enabled"] = st.checkbox(
                        "Enable Dormancy",
                        value=strains[i].get("dormancy_enabled", False),
                        key=f"str_dorm_en_{i}",
                    )
                    if strains[i]["dormancy_enabled"]:
                        cd1, cd2 = st.columns(2)
                        with cd1:
                            strains[i]["dormancy_depth"] = st.number_input(
                                "Depth layers (Q)",
                                min_value=1,
                                max_value=10,
                                value=int(strains[i].get("dormancy_depth", 1)),
                                key=f"str_depth_{i}",
                            )
                        with cd2:
                            strains[i]["dormancy_rate"] = st.number_input(
                                "Dormancy rate (sleep)",
                                value=float(strains[i].get("dormancy_rate", 0.001)),
                                step=0.05,
                                key=f"str_sleep_{i}",
                            )
                        cd3, cd4 = st.columns(2)
                        with cd3:
                            strains[i]["resuscitation_rate"] = st.number_input(
                                "Resuscitation rate (wake)",
                                value=float(strains[i].get("resuscitation_rate", 0.1)),
                                step=0.05,
                                key=f"str_wake_{i}",
                            )
                        with cd4:
                            strains[i]["death_rate_D"] = st.number_input(
                                "Dormant death rate (h⁻¹)",
                                value=float(strains[i].get("death_rate_D", 0.0)),
                                step=0.01,
                                key=f"str_death_d_{i}",
                            )
                        strains[i]["dormancy_diffusion_rate"] = st.number_input(
                            "Depth diffusion rate",
                            value=float(strains[i].get("dormancy_diffusion_rate", 0.05)),
                            step=0.01,
                            key=f"str_diff_{i}",
                        )
                        # Dormancy signal FUNCTIONS + their Ks/threshold are model-wide
                        # (set in the topmost panel) — only per-strain RATES live here.
                        strains[i]["initial_D"] = st.number_input(
                            "Initial dormant density (D0)",
                            value=float(strains[i].get("initial_D", 0.0)),
                            min_value=0.0,
                            format="%.1e",
                            key=f"str_init_d_{i}",
                            help="Dormant cells/mL at t=0. Defaults to 0. Distributed evenly across Q depth layers.",
                        )

        with col2:
            st.markdown("### Phage Strains")

            n_phages = counted_number_input("Number of phages", len(phages), "direct_n_phages", min_value=0, max_value=10)
            if n_phages != len(phages):
                st.session_state.simulation_result = None
                st.session_state.simulation_config = None
            # Adjust list size
            if n_phages > len(phages):
                for i in range(len(phages), n_phages):
                    phages.append(
                        {
                            "name": f"Phage {i}",
                            "initial_P": 1e6,
                            "adsorption_rates": 1e-8,
                            "adsorption_rates_dormant": 0.0,
                            "burst_sizes": 50.0,
                            "latent_periods": 0.5,
                            "phage_decay_rates": 0.1,
                            "pk_mode": "None",
                            "Vc": 5000.0,
                            "k_elim": 0.2,
                            "k_in": 0.1,
                            "k_out": 0.05,
                            "Vi": 10.0,
                            "Km_elim": 0.0,
                            "phage_decay_Km": 0.0,
                            "hibernation_rate_s": 0.0,
                            "hibernation_rate_r": 0.0,
                            "lytic_resumption_rate_s": 0.0,
                            "lytic_resumption_rate_r": 0.0,
                        }
                    )
            elif n_phages < len(phages):
                phages = phages[:n_phages]
            st.session_state["int_phages"] = phages

            # Inputs per phage
            for i in range(n_phages):
                with st.expander(f"Phage {i}: {phages[i]['name']}", expanded=True):
                    phages[i]["name"] = st.text_input(
                        "Name", value=phages[i]["name"], key=f"phg_name_{i}"
                    )
                    cc3, cc4 = st.columns(2)
                    with cc3:
                        if _ratio_inoc:
                            st.caption("Phage inoculum (P₀) comes from the per-arm dose "
                                       "(MOI / PFU) in Calibration.")
                        else:
                            phages[i]["initial_P"] = st.number_input(
                                "Initial Density (P0)",
                                value=float(phages[i]["initial_P"]),
                                format="%.1e",
                                key=f"phg_init_{i}",
                            )
                    with cc4:
                        phages[i]["burst_sizes"] = st.number_input(
                            "Burst size (PFU/cell)",
                            value=float(phages[i]["burst_sizes"]),
                            step=10.0,
                            key=f"phg_burst_{i}",
                        )
                    phages[i]["latent_periods"] = st.number_input(
                        "Latent period (hours)",
                        value=float(phages[i]["latent_periods"]),
                        step=0.1,
                        key=f"phg_latent_{i}",
                    )
                    phages[i]["phage_decay_rates"] = st.number_input(
                        "Phage decay rate (m)",
                        value=float(phages[i]["phage_decay_rates"]),
                        step=0.05,
                        key=f"phg_decay_{i}",
                    )

                    # Cross-resistance adsorption rates
                    st.markdown("**Adsorption Rates**")
                    for s_idx in range(n_strains):
                        ads_key = f"ads_{s_idx}_{i}"
                        ads_dorm_key = f"ads_dorm_{s_idx}_{i}"
                        # suscept ads
                        st.session_state[ads_key] = st.number_input(
                            f"Adsorption to {strains[s_idx]['name']} (mL·h⁻¹)",
                            value=float(st.session_state.get(ads_key, 1e-8 if s_idx == 0 else 0.0)),
                            format="%.1e",
                            key=f"ads_input_{s_idx}_{i}",
                        )
                        # dormant ads
                        st.session_state[ads_dorm_key] = st.number_input(
                            f"Adsorption to dormant {strains[s_idx]['name']} (mL·h⁻¹)",
                            value=float(st.session_state.get(ads_dorm_key, 0.0)),
                            format="%.1e",
                            key=f"ads_dorm_input_{s_idx}_{i}",
                        )

                    # Advanced biological and PK options
                    st.markdown("**Advanced Phage Kinetics**")
                    phages[i]["phage_decay_Km"] = st.number_input(
                        "Decay Km (Michaelis-Menten)",
                        value=float(phages[i].get("phage_decay_Km", 0.0)),
                        format="%.1e",
                        key=f"phg_decay_km_{i}",
                        help="Phage decay saturation. Set to 0 to disable."
                    )
                    phages[i]["attenuation_rate"] = st.number_input(
                        "Dormant adsorption attenuation (per depth layer)",
                        value=float(phages[i].get("attenuation_rate", 0.0)),
                        min_value=0.0,
                        step=0.1,
                        key=f"phg_atten_{i}",
                        help=(
                            "Exponential decay of adsorption to dormant cells with dormancy "
                            "depth: effective rate = adsorption_dormant × exp(−attenuation × "
                            "depth layer). 0 = phage penetrates all dormant layers equally."
                        ),
                    )
                    phages[i]["infected_nutrient_consumption"] = st.number_input(
                        "Infected-cell nutrient consumption (×)",
                        value=float(phages[i].get("infected_nutrient_consumption", 0.0)),
                        min_value=0.0, step=0.25, format="%g", key=f"phg_infnut_{i}",
                        help="Substrate drawn by cells infected by THIS phage, × the uninfected "
                             "per-capita uptake, while building phage. 0 = off (legacy); >1 = a "
                             "hijacked cell consumes more. Per-phage (a phage×host property) — "
                             "matters in cocktails; needs a nutrient-tracking growth signal.")

                    st.markdown("**Pseudolysogeny & Hibernation**")
                    phages[i]["hibernation_rate_s"] = st.number_input("Susceptible I->H rate", value=float(phages[i].get("hibernation_rate_s", 0.0)), step=0.05, key=f"phg_hib_s_{i}")
                    phages[i]["hibernation_rate_r"] = st.number_input("Resistant I->H rate", value=float(phages[i].get("hibernation_rate_r", 0.0)), step=0.05, key=f"phg_hib_r_{i}")
                    phages[i]["lytic_resumption_rate_s"] = st.number_input("Susceptible extra H->I rate", value=float(phages[i].get("lytic_resumption_rate_s", 0.0)), step=0.05, key=f"phg_res_s_{i}")
                    phages[i]["lytic_resumption_rate_r"] = st.number_input("Resistant extra H->I rate", value=float(phages[i].get("lytic_resumption_rate_r", 0.0)), step=0.05, key=f"phg_res_r_{i}")

                    st.markdown("**Pharmacokinetics (PK)**")
                    phages[i]["pk_mode"] = st.selectbox(
                        "Phage PK Mode",
                        ["None", "Effect Compartment", "Mass-Conserving"],
                        index=["None", "Effect Compartment", "Mass-Conserving"].index(phages[i].get("pk_mode", "None")),
                        key=f"phg_pk_{i}",
                    )
                    if phages[i]["pk_mode"] != "None":
                        pk1, pk2 = st.columns(2)
                        with pk1:
                            phages[i]["Vc"] = st.number_input("Central volume (Vc mL)", value=float(phages[i].get("Vc", 5000.0)), key=f"phg_vc_{i}")
                            phages[i]["k_elim"] = st.number_input("Elimination rate k_elim (h⁻¹)", value=float(phages[i].get("k_elim", 0.2)), key=f"phg_kelim_{i}")
                        with pk2:
                            phages[i]["k_in"] = st.number_input("Inflow rate k_in (h⁻¹)", value=float(phages[i].get("k_in", 0.1)), key=f"phg_kin_{i}")
                            phages[i]["k_out"] = st.number_input("Outflow rate k_out (h⁻¹)", value=float(phages[i].get("k_out", 0.05)), key=f"phg_kout_{i}")

                        phages[i]["Km_elim"] = st.number_input(
                            "Elimination Km",
                            value=float(phages[i].get("Km_elim", 0.0)),
                            format="%.1e",
                            key=f"phg_km_elim_{i}",
                            help="Central elimination saturation. Set to 0 to disable."
                        )
                        pc1, pc2 = st.columns(2)
                        with pc1:
                            phages[i]["k12"] = st.number_input(
                                "Central→peripheral k12 (h⁻¹)", value=float(phages[i].get("k12", 0.0)),
                                min_value=0.0, step=0.05, key=f"phg_k12_{i}",
                                help="Set > 0 to add a 2-compartment peripheral tissue compartment "
                                     "(Pp). 0 = one-compartment central model.")
                        with pc2:
                            phages[i]["k21"] = st.number_input(
                                "Peripheral→central k21 (h⁻¹)", value=float(phages[i].get("k21", 0.0)),
                                min_value=0.0, step=0.05, key=f"phg_k21_{i}",
                                help="Return transfer from the peripheral compartment. Used only when k12 > 0.")
                        if phages[i]["pk_mode"] == "Mass-Conserving":
                            phages[i]["Vi"] = st.number_input("Infection Site Volume (Vi mL)", value=float(phages[i].get("Vi", 10.0)), key=f"phg_vi_{i}")

            if n_phages > 0:
                st.markdown("---")
                st.markdown("### Bacterial Mutations (WT → R)")
                if n_strains == 2**n_phages:
                    st.caption("Per-phage-locus shortcut (binary-genotype layout, n_strains = 2^n_phages):")
                    phg_res_rates = []
                    _prev_mu = st.session_state.get("direct_phg_res_rates", [])
                    for j in range(n_phages):
                        # Seed from the persisted direct_phg_res_rates (a scenario data
                        # key) so a value of 0 survives navigation. The widget key
                        # direct_mu_{j} is dropped when the page isn't rendered, so
                        # value= must come from a key that actually persists — the old
                        # code read `direct_phg_mu_{j}`, which was never written, so it
                        # reverted to 1e-7 on return.
                        res_rate = st.number_input(
                            f"Mutation rate to {phages[j]['name']} resistance (mu)",
                            value=float(_prev_mu[j]) if j < len(_prev_mu) else 1e-7,
                            format="%.1e",
                            key=f"direct_mu_{j}"
                        )
                        phg_res_rates.append(res_rate)
                    st.session_state["direct_phg_res_rates"] = phg_res_rates

                with st.expander("Custom mutation network (any number of strains)",
                                 expanded=(n_strains != 2**n_phages)):
                    st.caption(
                        "Define arbitrary strain→strain mutation edges (per replication). Works "
                        "for any strain count (lifts the 2^n_phages requirement). **If any edge "
                        "is added here it overrides the per-locus shortcut above.** For "
                        "growth-independent phenotypic switching use the transitions section "
                        "below the builder columns."
                    )
                    render_mutation_graph_editor(strains, key_prefix="dir_trans")

    # ── BINARY RESISTANCE GENOTYPES (BRG) ──
    elif builder_mode == "Binary Genotypes (BRG)":
        with col1:
            st.markdown("### Base Bacteria (WT)")
            st.session_state["int_brg_base_growth"] = st.number_input(
                "Base growth rate (h⁻¹)", value=float(st.session_state.get("int_brg_base_growth", 1.2)), step=0.1
            )
            st.session_state["int_brg_base_ratio"] = st.number_input(
                "Resource consumption ratio", value=float(st.session_state.get("int_brg_base_ratio", 1e9)), format="%.1e"
            )
            st.session_state["int_brg_death_rate_B"] = st.number_input(
                "Natural death rate (h⁻¹)", value=float(st.session_state.get("int_brg_death_rate_B", 0.0)), step=0.01
            )
            st.session_state["int_brg_dormancy_enabled"] = st.checkbox(
                "Enable Dormancy", value=st.session_state.get("int_brg_dormancy_enabled", False)
            )
            if st.session_state["int_brg_dormancy_enabled"]:
                st.session_state["int_brg_dorm_rate"] = st.number_input(
                    "Dormancy rate (sleep)", value=float(st.session_state.get("int_brg_dorm_rate", 0.001)), step=0.05
                )
                st.session_state["int_brg_resus_rate"] = st.number_input(
                    "Resuscitation rate (wake)", value=float(st.session_state.get("int_brg_resus_rate", 0.1)), step=0.05
                )
                st.session_state["int_brg_diff_rate"] = st.number_input(
                    "Depth diffusion rate", value=float(st.session_state.get("int_brg_diff_rate", 0.05)), step=0.01
                )
                st.session_state["int_brg_n_depth"] = st.number_input(
                    "Dormancy depth layers (Q)",
                    min_value=1, max_value=10,
                    value=int(st.session_state.get("int_brg_n_depth", 1)),
                    help="Number of dormancy-depth compartments for all genotypes.",
                )
                st.session_state["int_brg_death_rate_D"] = st.number_input(
                    "Dormant death rate (h⁻¹)", value=float(st.session_state.get("int_brg_death_rate_D", 0.0)), step=0.01
                )
                st.caption("Dormancy/resuscitation/diffusion **signal functions** (+ their Ks / "
                           "density threshold) are model-wide — set them in the topmost builder panel.")

            # Renders the loci count
            st.markdown("---")
            st.markdown("### Phage Loci")
            n_phg_loci = counted_number_input("Number of phage species (loci)", len(phages), "brg_n_phg_loci", min_value=1, max_value=10)
            if n_phg_loci != len(phages):
                phages = phages[:n_phg_loci]
                while len(phages) < n_phg_loci:
                    phages.append({
                        "name": f"Phage {len(phages)}",
                        "initial_P": 1e6,
                        "burst_sizes": 50.0,
                        "latent_periods": 0.5,
                        "phage_decay_rates": 0.1,
                        "pk_mode": "None",
                        "Vc": 5000.0, "k_elim": 0.2, "k_in": 0.1, "k_out": 0.05, "Vi": 10.0,
                        "adsorption_s": 5e-8,
                        "adsorption_r": 0.0,
                        "adsorption_dormant_s": 0.0,
                        "adsorption_dormant_r": 0.0,
                        "fitness_cost": 0.05,
                        "mu": 1e-7,
                    })
                st.session_state["int_phages"] = phages

            for idx in range(n_phg_loci):
                with st.expander(f"Phage Locus {idx}: {phages[idx]['name']}", expanded=True):
                    phages[idx]["name"] = st.text_input("Locus name", value=phages[idx]["name"], key=f"brg_phg_name_{idx}")
                    if _ratio_inoc:
                        st.caption("Phage inoculum (P₀) comes from the per-arm dose in Calibration.")
                    else:
                        phages[idx]["initial_P"] = st.number_input("Initial count P₀ (PFU·mL⁻¹)", value=float(phages[idx]["initial_P"]), format="%.1e", key=f"brg_phg_init_{idx}")
                    phages[idx]["adsorption_s"] = st.number_input("Adsorption WT (mL·h⁻¹)", value=float(phages[idx].get("adsorption_s", 5e-8)), format="%.2e", key=f"brg_phg_ads_s_{idx}")
                    phages[idx]["adsorption_r"] = st.number_input("Adsorption Res (mL·h⁻¹)", value=float(phages[idx].get("adsorption_r", 0.0)), format="%.2e", key=f"brg_phg_ads_r_{idx}")
                    phages[idx]["adsorption_dormant_s"] = st.number_input(
                        "Adsorption to dormant WT (mL·h⁻¹)", value=float(phages[idx].get("adsorption_dormant_s", 0.0)),
                        format="%.2e", key=f"brg_phg_ads_dorm_s_{idx}",
                        help="Adsorption rate to dormant/persister WT cells (0 = phage cannot infect dormant cells).")
                    phages[idx]["adsorption_dormant_r"] = st.number_input(
                        "Adsorption to dormant Res (mL·h⁻¹)", value=float(phages[idx].get("adsorption_dormant_r", 0.0)),
                        format="%.2e", key=f"brg_phg_ads_dorm_r_{idx}")
                    phages[idx]["burst_sizes"] = st.number_input("Burst size (PFU/cell)", value=float(phages[idx]["burst_sizes"]), step=10.0, key=f"brg_phg_burst_{idx}")
                    phages[idx]["latent_periods"] = st.number_input("Latent period (h)", value=float(phages[idx]["latent_periods"]), step=0.1, key=f"brg_phg_latent_{idx}")
                    phages[idx]["phage_decay_rates"] = st.number_input("Phage decay rate (h⁻¹)", value=float(phages[idx]["phage_decay_rates"]), step=0.05, key=f"brg_phg_decay_{idx}")
                    phages[idx]["fitness_cost"] = st.number_input(
                        "Resistance fitness cost",
                        value=float(phages[idx].get("fitness_cost", 0.05)), step=0.01,
                        key=f"brg_phg_fit_{idx}",
                        help="Fractional growth-rate penalty of phage-resistant genotypes "
                             "(resistant growth = base × (1 − cost)). Drives the equilibrium "
                             "initial condition — with cost 0 the resistant mutants are neutral "
                             "and dominate the pre-treatment equilibrium.",
                    )
                    phages[idx]["mu"] = st.number_input("Mutation rate μ (per replication)", value=float(phages[idx].get("mu", 1e-7)), format="%.1e", key=f"brg_phg_mu_{idx}")
                    phages[idx]["attenuation_rate"] = st.number_input(
                        "Dormant adsorption attenuation (per depth layer)",
                        value=float(phages[idx].get("attenuation_rate", 0.0)),
                        min_value=0.0, step=0.1, key=f"brg_phg_atten_{idx}",
                        help="Exponential decay of dormant-cell adsorption with dormancy depth (0 = none).",
                    )
                    phages[idx]["infected_nutrient_consumption"] = st.number_input(
                        "Infected-cell nutrient consumption (×)",
                        value=float(phages[idx].get("infected_nutrient_consumption", 0.0)),
                        min_value=0.0, step=0.25, format="%g", key=f"brg_phg_infnut_{idx}",
                        help="Substrate drawn by cells infected by THIS phage, × the uninfected "
                             "per-capita uptake (0 = off; >1 = a hijacked cell consumes more). "
                             "Per-phage; needs a nutrient-tracking growth signal.")

                    st.markdown("**Advanced Phage Kinetics**")
                    phages[idx]["phage_decay_Km"] = st.number_input(
                        "Decay Km (Michaelis-Menten)",
                        value=float(phages[idx].get("phage_decay_Km", 0.0)),
                        format="%.1e",
                        key=f"brg_phg_decay_km_{idx}",
                        help="Phage decay saturation. Set to 0 to disable.",
                    )

                    st.markdown("**Pseudolysogeny & Hibernation**")
                    phages[idx]["hibernation_rate_s"] = st.number_input("Susceptible I->H rate", value=float(phages[idx].get("hibernation_rate_s", 0.0)), step=0.05, key=f"brg_phg_hib_s_{idx}")
                    phages[idx]["hibernation_rate_r"] = st.number_input("Resistant I->H rate", value=float(phages[idx].get("hibernation_rate_r", 0.0)), step=0.05, key=f"brg_phg_hib_r_{idx}")
                    phages[idx]["lytic_resumption_rate_s"] = st.number_input("Susceptible extra H->I rate", value=float(phages[idx].get("lytic_resumption_rate_s", 0.0)), step=0.05, key=f"brg_phg_res_s_{idx}")
                    phages[idx]["lytic_resumption_rate_r"] = st.number_input("Resistant extra H->I rate", value=float(phages[idx].get("lytic_resumption_rate_r", 0.0)), step=0.05, key=f"brg_phg_res_r_{idx}")

                    st.markdown("**Pharmacokinetics (PK)**")
                    phages[idx]["pk_mode"] = st.selectbox(
                        "Phage PK Mode",
                        ["None", "Effect Compartment", "Mass-Conserving"],
                        index=["None", "Effect Compartment", "Mass-Conserving"].index(phages[idx].get("pk_mode", "None")),
                        key=f"brg_phg_pkmode_{idx}",
                    )
                    if phages[idx]["pk_mode"] != "None":
                        bpk1, bpk2 = st.columns(2)
                        with bpk1:
                            phages[idx]["Vc"] = st.number_input("Central volume (Vc mL)", value=float(phages[idx].get("Vc", 5000.0)), key=f"brg_phg_vc_{idx}")
                            phages[idx]["k_elim"] = st.number_input("Elimination rate k_elim (h⁻¹)", value=float(phages[idx].get("k_elim", 0.2)), key=f"brg_phg_kelim_{idx}")
                        with bpk2:
                            phages[idx]["k_in"] = st.number_input("Inflow rate k_in (h⁻¹)", value=float(phages[idx].get("k_in", 0.1)), key=f"brg_phg_kin_{idx}")
                            phages[idx]["k_out"] = st.number_input("Outflow rate k_out (h⁻¹)", value=float(phages[idx].get("k_out", 0.05)), key=f"brg_phg_kout_{idx}")
                        phages[idx]["Km_elim"] = st.number_input(
                            "Elimination Km",
                            value=float(phages[idx].get("Km_elim", 0.0)),
                            format="%.1e",
                            key=f"brg_phg_km_elim_{idx}",
                            help="Central elimination saturation. Set to 0 to disable.",
                        )
                        bpc1, bpc2 = st.columns(2)
                        with bpc1:
                            phages[idx]["k12"] = st.number_input(
                                "Central→peripheral k12 (h⁻¹)", value=float(phages[idx].get("k12", 0.0)),
                                min_value=0.0, step=0.05, key=f"brg_phg_k12_{idx}",
                                help="Set > 0 to add a 2-compartment peripheral tissue compartment "
                                     "(Pp). 0 = one-compartment central model.")
                        with bpc2:
                            phages[idx]["k21"] = st.number_input(
                                "Peripheral→central k21 (h⁻¹)", value=float(phages[idx].get("k21", 0.0)),
                                min_value=0.0, step=0.05, key=f"brg_phg_k21_{idx}",
                                help="Return transfer from the peripheral compartment. Used only when k12 > 0.")
                        if phages[idx]["pk_mode"] == "Mass-Conserving":
                            phages[idx]["Vi"] = st.number_input("Infection Site Volume (Vi mL)", value=float(phages[idx].get("Vi", 10.0)), key=f"brg_phg_vi_{idx}")

        with col2:
            st.markdown("### Auto-generated Genotypes")
            # Show list of 2^(m+a) genotypes and initial conditions inputs
            n_abx = len(antibiotics)
            st.caption(f"Based on {n_phg_loci} phages and {n_abx} antibiotics, there are "
                       f"{2 ** (n_phg_loci + n_abx)} genotypes:")

            st.session_state["int_brg_use_eq_ic"] = st.checkbox(
                "Use equilibrium initial condition",
                value=st.session_state.get("int_brg_use_eq_ic", False),
                help="Compute B0 per genotype from replicator-dynamics equilibrium (BRG growth rates + mutation matrix). Overrides per-genotype inputs.",
            )
            if st.session_state["int_brg_use_eq_ic"]:
                if _ratio_inoc:
                    st.caption("Genotype **ratios** come from the equilibrium (fitness cost + "
                               "mutation); the absolute total B₀ is set per-arm in the "
                               "Calibration conditions.")
                else:
                    st.session_state["int_brg_eq_total_B"] = st.number_input(
                        "Total bacteria at t=0 (cells/mL)",
                        value=float(st.session_state.get("int_brg_eq_total_B", 1e7)),
                        format="%.1e",
                        key="brg_eq_total_B",
                    )
                    st.caption("Per-genotype B0 computed from `brg.equilibrium_initial_condition()` at run time.")
            else:
                brg_initial_B = st.session_state.get("int_brg_initial_B", {})
                for idx, lbl in enumerate(brg_genotype_labels(n_phg_loci, n_abx)):
                    brg_initial_B[lbl] = st.number_input(
                        (f"Initial ratio for genotype {lbl}" if _ratio_inoc
                         else f"Initial count for genotype {lbl}"),
                        value=float(brg_initial_B.get(lbl, 1e7 if idx == 0 else 0.0)),
                        format="%g" if _ratio_inoc else "%.1e",
                        key=f"brg_init_B_{lbl}"
                    )
                st.session_state.int_brg_initial_B = brg_initial_B

    # ── CUSTOM STRAINS & MUTATION GRAPH (StrainSet) ──
    elif builder_mode == "Custom Strains & Graph (StrainSet)":
        with col1:
            st.markdown("### Custom Bacterial Strains")

            n_strains = counted_number_input("Number of custom strains", len(strains), "ss_n_strains", min_value=1, max_value=10)
            if n_strains != len(strains):
                strains = strains[:n_strains]
                while len(strains) < n_strains:
                    strains.append({
                        "name": f"Strain {len(strains)}",
                        "initial_B": 1e7,
                        "growth_rate": 1.2,
                        "bacteria_to_resource_ratio": 1e9,
                        "death_rate_B": 0.0,
                        "death_rate_D": 0.0,
                        "dormancy_enabled": False,
                        "dormancy_rate": 0.001, "resuscitation_rate": 0.1, "dormancy_diffusion_rate": 0.05
                    })
                st.session_state["int_strains"] = strains

            # Renders expanders per strain
            for i in range(n_strains):
                with st.expander(f"Strain {i}: {strains[i]['name']}", expanded=True):
                    strains[i]["name"] = st.text_input("Strain name", value=strains[i]["name"], key=f"ss_str_name_{i}")
                    strains[i]["initial_B"] = st.number_input(
                        "Initial ratio (relative)" if _ratio_inoc else "Initial count B₀ (CFU·mL⁻¹)",
                        value=float(strains[i]["initial_B"]),
                        format="%g" if _ratio_inoc else "%.1e", key=f"ss_str_init_{i}",
                        help=("Relative amount across strains; absolute total from the per-arm B₀ "
                              "in Calibration." if _ratio_inoc else None))
                    strains[i]["growth_rate"] = st.number_input("Growth rate (h⁻¹)", value=float(strains[i]["growth_rate"]), step=0.1, key=f"ss_str_growth_{i}")
                    strains[i]["bacteria_to_resource_ratio"] = st.number_input(
                        "Bacteria-to-resource ratio", value=float(strains[i].get("bacteria_to_resource_ratio", 1e9)),
                        format="%.2e", key=f"ss_str_ratio_{i}",
                        help="Bacteria produced per unit resource consumed (yield). Governs how fast "
                             "growth depletes the substrate under nutrient-limited (Monod) growth.")
                    strains[i]["death_rate_B"] = st.number_input("Natural death rate (h⁻¹)", value=float(strains[i].get("death_rate_B", 0.0)), step=0.01, key=f"ss_str_death_{i}")

                    strains[i]["dormancy_enabled"] = st.checkbox("Enable Dormancy", value=strains[i].get("dormancy_enabled", False), key=f"ss_str_dorm_{i}")
                    if strains[i]["dormancy_enabled"]:
                        strains[i]["dormancy_depth"] = st.number_input("Depth layers (Q)", min_value=1, max_value=10, value=int(strains[i].get("dormancy_depth", 1)), key=f"ss_str_depth_{i}", help="Number of dormancy-depth compartments (max across strains sets the model n_depth).")
                        strains[i]["dormancy_rate"] = st.number_input("Dormancy rate (h⁻¹)", value=float(strains[i].get("dormancy_rate", 0.001)), key=f"ss_str_sleep_{i}")
                        strains[i]["resuscitation_rate"] = st.number_input("Resuscitation rate (h⁻¹)", value=float(strains[i].get("resuscitation_rate", 0.1)), key=f"ss_str_wake_{i}")
                        strains[i]["dormancy_diffusion_rate"] = st.number_input("Depth diffusion (h⁻¹)", value=float(strains[i].get("dormancy_diffusion_rate", 0.05)), key=f"ss_str_diff_{i}")
                        strains[i]["death_rate_D"] = st.number_input("Dormant death rate (h⁻¹)", value=float(strains[i].get("death_rate_D", 0.0)), step=0.01, key=f"ss_str_death_d_{i}")
                        # Signal FUNCTIONS + their Ks/threshold are model-wide (topmost panel);
                        # only per-strain rates live here.

                    if len(phages) > 0:
                        st.markdown("**Phage Adsorption Rates**")
                        for p_idx in range(len(phages)):
                            p_name = phages[p_idx]["name"]
                            st.session_state[f"ads_{i}_{p_idx}"] = st.number_input(
                                f"Adsorption of {p_name} (mL·h⁻¹)",
                                value=float(st.session_state.get(f"ads_{i}_{p_idx}", 1e-8 if i == 0 else 0.0)),
                                format="%.1e",
                                key=f"ss_ads_input_{i}_{p_idx}"
                            )
                            if strains[i]["dormancy_enabled"]:
                                st.session_state[f"ads_dorm_{i}_{p_idx}"] = st.number_input(
                                    f"Dormant adsorption of {p_name} (mL·h⁻¹)",
                                    value=float(st.session_state.get(f"ads_dorm_{i}_{p_idx}", 0.0)),
                                    format="%.1e",
                                    key=f"ss_ads_dorm_input_{i}_{p_idx}"
                                )

            # Mutation graph editor (division-coupled). Phenotypic transitions and the
            # phage generators live in the shared section below the builder columns.
            st.markdown("#### Mutation Graph")
            st.caption("Strain→strain mutation at division (rate per replication). Mass-conserving: "
                       "each origin's diagonal is rebalanced automatically.")
            render_mutation_graph_editor(strains, key_prefix="trans")

        with col2:
            st.markdown("### Phage Strains")
            n_phages = counted_number_input("Number of phages", len(phages), "ss_n_phages", min_value=0, max_value=10)
            if n_phages != len(phages):
                phages = phages[:n_phages]
                while len(phages) < n_phages:
                    phages.append({
                        "name": f"Phage {len(phages)}",
                        "initial_P": 1e6,
                        "burst_sizes": 50.0,
                        "latent_periods": 0.5,
                        "phage_decay_rates": 0.1,
                        "pk_mode": "None",
                        "Vc": 5000.0, "k_elim": 0.2, "k_in": 0.1, "k_out": 0.05, "Vi": 10.0,
                    })
                st.session_state["int_phages"] = phages

            for idx in range(n_phages):
                with st.expander(f"Phage {idx}: {phages[idx]['name']}", expanded=True):
                    phages[idx]["name"] = st.text_input("Phage name", value=phages[idx]["name"], key=f"ss_phg_name_{idx}")
                    if _ratio_inoc:
                        st.caption("Phage inoculum (P₀) comes from the per-arm dose in Calibration.")
                    else:
                        phages[idx]["initial_P"] = st.number_input("Initial count P₀ (PFU·mL⁻¹)", value=float(phages[idx]["initial_P"]), format="%.1e", key=f"ss_phg_init_{idx}")
                    phages[idx]["burst_sizes"] = st.number_input("Burst size (PFU/cell)", value=float(phages[idx].get("burst_sizes", 50.0)), step=10.0, key=f"ss_phg_burst_{idx}")
                    phages[idx]["latent_periods"] = st.number_input("Latent period (h)", value=float(phages[idx].get("latent_periods", 0.5)), step=0.1, key=f"ss_phg_latent_{idx}")
                    phages[idx]["phage_decay_rates"] = st.number_input("Phage decay rate (h⁻¹)", value=float(phages[idx]["phage_decay_rates"]), step=0.05, key=f"ss_phg_decay_{idx}")
                    phages[idx]["attenuation_rate"] = st.number_input(
                        "Dormant adsorption attenuation (per depth layer)",
                        value=float(phages[idx].get("attenuation_rate", 0.0)),
                        min_value=0.0, step=0.1, key=f"ss_phg_atten_{idx}",
                        help="Exponential decay of dormant-cell adsorption with dormancy depth (0 = none).",
                    )
                    phages[idx]["infected_nutrient_consumption"] = st.number_input(
                        "Infected-cell nutrient consumption (×)",
                        value=float(phages[idx].get("infected_nutrient_consumption", 0.0)),
                        min_value=0.0, step=0.25, format="%g", key=f"ss_phg_infnut_{idx}",
                        help="Substrate drawn by cells infected by THIS phage, × the uninfected "
                             "per-capita uptake (0 = off; >1 = a hijacked cell consumes more). "
                             "Per-phage; needs a nutrient-tracking growth signal.")

                    st.markdown("**Advanced Phage Kinetics**")
                    phages[idx]["phage_decay_Km"] = st.number_input(
                        "Decay Km (Michaelis-Menten)",
                        value=float(phages[idx].get("phage_decay_Km", 0.0)),
                        format="%.1e",
                        key=f"ss_phg_decay_km_{idx}",
                        help="Phage decay saturation. Set to 0 to disable.",
                    )

                    st.markdown("**Pseudolysogeny & Hibernation**")
                    phages[idx]["hibernation_rate_s"] = st.number_input("Susceptible I->H rate", value=float(phages[idx].get("hibernation_rate_s", 0.0)), step=0.05, key=f"ss_phg_hib_s_{idx}")
                    phages[idx]["hibernation_rate_r"] = st.number_input("Resistant I->H rate", value=float(phages[idx].get("hibernation_rate_r", 0.0)), step=0.05, key=f"ss_phg_hib_r_{idx}")
                    phages[idx]["lytic_resumption_rate_s"] = st.number_input("Susceptible extra H->I rate", value=float(phages[idx].get("lytic_resumption_rate_s", 0.0)), step=0.05, key=f"ss_phg_res_s_{idx}")
                    phages[idx]["lytic_resumption_rate_r"] = st.number_input("Resistant extra H->I rate", value=float(phages[idx].get("lytic_resumption_rate_r", 0.0)), step=0.05, key=f"ss_phg_res_r_{idx}")

                    st.markdown("**Pharmacokinetics (PK)**")
                    phages[idx]["pk_mode"] = st.selectbox(
                        "Phage PK Mode",
                        ["None", "Effect Compartment", "Mass-Conserving"],
                        index=["None", "Effect Compartment", "Mass-Conserving"].index(phages[idx].get("pk_mode", "None")),
                        key=f"ss_phg_pkmode_{idx}",
                    )
                    if phages[idx]["pk_mode"] != "None":
                        spk1, spk2 = st.columns(2)
                        with spk1:
                            phages[idx]["Vc"] = st.number_input("Central volume (Vc mL)", value=float(phages[idx].get("Vc", 5000.0)), key=f"ss_phg_vc_{idx}")
                            phages[idx]["k_elim"] = st.number_input("Elimination rate k_elim (h⁻¹)", value=float(phages[idx].get("k_elim", 0.2)), key=f"ss_phg_kelim_{idx}")
                        with spk2:
                            phages[idx]["k_in"] = st.number_input("Inflow rate k_in (h⁻¹)", value=float(phages[idx].get("k_in", 0.1)), key=f"ss_phg_kin_{idx}")
                            phages[idx]["k_out"] = st.number_input("Outflow rate k_out (h⁻¹)", value=float(phages[idx].get("k_out", 0.05)), key=f"ss_phg_kout_{idx}")
                        phages[idx]["Km_elim"] = st.number_input(
                            "Elimination Km",
                            value=float(phages[idx].get("Km_elim", 0.0)),
                            format="%.1e",
                            key=f"ss_phg_km_elim_{idx}",
                            help="Central elimination saturation. Set to 0 to disable.",
                        )
                        spc1, spc2 = st.columns(2)
                        with spc1:
                            phages[idx]["k12"] = st.number_input(
                                "Central→peripheral k12 (h⁻¹)", value=float(phages[idx].get("k12", 0.0)),
                                min_value=0.0, step=0.05, key=f"ss_phg_k12_{idx}",
                                help="Set > 0 to add a 2-compartment peripheral tissue compartment "
                                     "(Pp). 0 = one-compartment central model.")
                        with spc2:
                            phages[idx]["k21"] = st.number_input(
                                "Peripheral→central k21 (h⁻¹)", value=float(phages[idx].get("k21", 0.0)),
                                min_value=0.0, step=0.05, key=f"ss_phg_k21_{idx}",
                                help="Return transfer from the peripheral compartment. Used only when k12 > 0.")
                        if phages[idx]["pk_mode"] == "Mass-Conserving":
                            phages[idx]["Vi"] = st.number_input("Infection Site Volume (Vi mL)", value=float(phages[idx].get("Vi", 10.0)), key=f"ss_phg_vi_{idx}")

    render_generator_matrix_editors(builder_mode)
    return builder_mode


def render_generator_matrix_editors(builder_mode):
    """Phenotypic transitions (bacteria) + phage mutation / transition generators —
    the three engine generator matrices no builder column exposes (the bacterial
    MUTATION matrix is configured per mode above: per-locus μ / mutation graph).

    Mode-aware node labels: Direct & StrainSet use the strain names, BRG the
    auto-generated genotype labels; phages use their names in every mode. Wired into
    the build via ``generator_matrices_from_state`` (Direct → ``with_mutations``,
    StrainSet → ``set_transition_graph`` / ``set_phage_*_rates``, BRG → post-build
    ``cfg.<field>``)."""
    strains = st.session_state.get("int_strains", [])
    phages = st.session_state.get("int_phages", [])
    antibiotics = st.session_state.get("int_antibiotics", [])
    if builder_mode == "Binary Genotypes (BRG)":
        bact_names = brg_genotype_labels(len(phages), len(antibiotics))
    else:
        bact_names = [s["name"] for s in strains]
    phage_names = [p["name"] for p in phages]
    n_edges = sum(len(st.session_state.get(k, []) or [])
                  for k in ("int_bact_transitions", "int_phage_mutations", "int_phage_transitions"))
    if len(bact_names) < 2 and len(phage_names) < 2:
        return   # nothing to switch between

    st.markdown("---")
    with st.expander("Phenotypic transitions & phage evolution (generator matrices)",
                     expanded=n_edges > 0):
        st.caption(
            "Optional rate graphs the engine applies as mass-conserving generator matrices "
            "(`M[dest, origin]`, diagonal = −outflow). **Mutation** is coupled to replication "
            "(bacteria: per division; phage: per burst), so it vanishes when growth / lysis "
            "stops. **Phenotypic transitions** are first-order switches in h⁻¹ applied to the "
            "population directly (persister-like / phase-variation switching, phage "
            "state changes) — they act even in stationary phase."
        )
        if len(bact_names) >= 2:
            st.markdown("**Bacterial phenotypic transitions** (h⁻¹, growth-independent)")
            if builder_mode == "Binary Genotypes (BRG)":
                st.caption("Nodes are the BRG genotype labels (bacterial mutation itself is the "
                           "per-locus μ in the phage-locus cards).")
            render_rate_graph_editor(bact_names, "int_bact_transitions", "gen_btrans",
                                     rate_label="Rate (h⁻¹)", default_rate=0.01,
                                     add_label="+ Add bacterial transition")
        if len(phage_names) >= 2:
            st.markdown("**Phage mutation** (fraction of each burst switching phage type)")
            render_rate_graph_editor(phage_names, "int_phage_mutations", "gen_pmut",
                                     rate_label="Rate (per burst)", default_rate=1e-6,
                                     add_label="+ Add phage mutation")
            st.markdown("**Phage phenotypic transitions** (h⁻¹, applied to free phage)")
            render_rate_graph_editor(phage_names, "int_phage_transitions", "gen_ptrans",
                                     rate_label="Rate (h⁻¹)", default_rate=0.01,
                                     add_label="+ Add phage transition")


def render():
    theme_mode = st.session_state.get("theme_mode", "Light")
    st.title("Interactive Simulation Builder")
    st.caption("Configure custom variables, build mathematical parameters, and solve the ODE.")

    st.markdown(
        "<div class='info-banner'>Configure your biological layers using "
        "the tabs below, then scroll to the bottom and click <b>Run Simulation</b>!</div>",
        unsafe_allow_html=True,
    )

    # This page edits the live builder draft = the active Model. Loading a saved/demo
    # Model (sidebar) populates these widgets; edits here diverge from it until re-saved.
    _am = st.session_state.get("active_model", WORKING_DRAFT_LABEL)
    if _am != WORKING_DRAFT_LABEL:
        st.caption(f"Editing a working copy of model **{_am}** — save changes as a Model "
                   "in the sidebar to keep them; downstream tasks can then select it.")

    # Click-to-view snapshot of the whole configured model (the live draft being built
    # here) — every section resolved in one place, no tab-hunting. Toggle gates the build.
    if st.toggle("📋 Show current model config", key="sim_show_cfg",
                 help="A full snapshot of the model you're building — strains, phages, growth "
                      "& nutrient environment, immunity, OD/debris, initial conditions, solver."):
        render_model_snapshot()

    # 1. Retrieve current lists from state
    strains = st.session_state.get("int_strains", [])
    phages = st.session_state.get("int_phages", [])
    antibiotics = st.session_state.get("int_antibiotics", [])
    doses = st.session_state.get("int_doses", [])
    track_nutrients = st.session_state.get("int_track_nutrients", True)

    # 2. Main tabs for parameters configuration
    config_tabs = st.tabs(
        [
            "Strains & Phages",
            "Antibiotics & Immunity",
            "Environment & Dosing",
            "Solver Settings",
        ]
    )

    # ──── Tab 1: Strains & Phages ─────────────────────────────────────────────
    with config_tabs[0]:
        builder_mode = render_model_builder()

    # ──── Tab 2: Antibiotics & Immunity ───────────────────────────────────────
    with config_tabs[1]:
        col1, col2 = st.columns(2)

        # Antibiotics
        with col1:
            st.markdown("### Antibiotics")

            n_abx = st.number_input(
                "Number of antibiotics", min_value=0, max_value=6, value=len(antibiotics)
            )
            if n_abx != len(antibiotics):
                st.session_state.simulation_result = None
                st.session_state.simulation_config = None
            # Adjust list size
            if n_abx > len(antibiotics):
                for i in range(len(antibiotics), n_abx):
                    antibiotics.append(
                        {
                            "name": f"Drug {i}",
                            "Vc": 250.0,
                            "k_elim": 0.3,
                            "k12": 0.0,
                            "k21": 0.0,
                            "emax": 3.0,
                            "ec50": 0.2,
                            "hill": 1.5,
                            "f_lyse": 0.0,
                            "inoculum_effect_constant": 0.0,
                            "Km_elim": 0.0,
                            "emax_r": 0.3,
                            "ec50_r": 2.0,
                            "fitness_cost": 0.05,
                            "mu": 1e-7,
                        }
                    )
            elif n_abx < len(antibiotics):
                antibiotics = antibiotics[:n_abx]
            st.session_state["int_antibiotics"] = antibiotics

            # Inputs per antibiotic
            for i in range(n_abx):
                with st.expander(
                    f"Antibiotic {i}: {antibiotics[i]['name']}", expanded=True
                ):
                    antibiotics[i]["name"] = st.text_input(
                        "Name", value=antibiotics[i]["name"], key=f"abx_name_{i}"
                    )
                    c1, c2 = st.columns(2)
                    with c1:
                        antibiotics[i]["Vc"] = st.number_input(
                            "Volume (Vc L)",
                            value=float(antibiotics[i].get("Vc", 250.0)),
                            key=f"abx_vc_{i}",
                        )
                        antibiotics[i]["k_elim"] = st.number_input(
                            "Clearance k_elim (h⁻¹)",
                            value=float(antibiotics[i]["k_elim"]),
                            step=0.05,
                            key=f"abx_kelim_{i}",
                        )
                    with c2:
                        antibiotics[i]["k12"] = st.number_input(
                            "Transfer central->peripheral (k12)",
                            value=float(antibiotics[i].get("k12", 0.0)),
                            step=0.05,
                            key=f"abx_k12_{i}",
                        )
                        antibiotics[i]["k21"] = st.number_input(
                            "Transfer peripheral->central (k21)",
                            value=float(antibiotics[i].get("k21", 0.0)),
                            step=0.05,
                            key=f"abx_k21_{i}",
                        )

                    st.markdown("**Pharmacodynamics (PD)**")
                    
                    if builder_mode == "Binary Genotypes (BRG)":
                        # Susceptible vs Resistant inputs
                        antibiotics[i]["emax"] = st.number_input("Susceptible Max Efficacy (Emax_s)", value=float(antibiotics[i].get("emax", 3.0)), key=f"abx_emax_s_{i}")
                        antibiotics[i]["emax_r"] = st.number_input("Resistant Max Efficacy (Emax_r)", value=float(antibiotics[i].get("emax_r", 0.3)), key=f"abx_emax_r_{i}")
                        
                        antibiotics[i]["ec50"] = st.number_input("Susceptible EC50_s", value=float(antibiotics[i].get("ec50", 0.2)), key=f"abx_ec50_s_{i}")
                        antibiotics[i]["ec50_r"] = st.number_input("Resistant EC50_r", value=float(antibiotics[i].get("ec50_r", 2.0)), key=f"abx_ec50_r_{i}")
                        
                        antibiotics[i]["fitness_cost"] = st.number_input(
                            "Resistance fitness cost",
                            value=float(antibiotics[i].get("fitness_cost", 0.05)), step=0.01,
                            key=f"abx_fit_{i}",
                            help="Fractional growth-rate penalty of antibiotic-resistant genotypes. "
                                 "Drives the equilibrium initial condition; cost 0 → resistant "
                                 "mutants dominate the pre-treatment equilibrium.",
                        )
                        antibiotics[i]["mu"] = st.number_input("Mutation rate μ (per replication)", value=float(antibiotics[i].get("mu", 1e-7)), format="%.1e", key=f"abx_mu_{i}")
                    else:
                        # Direct parameters
                        c3, c4 = st.columns(2)
                        with c3:
                            antibiotics[i]["emax"] = st.number_input(
                                "Max efficacy (Emax)",
                                value=float(antibiotics[i]["emax"]),
                                step=0.5,
                                key=f"abx_emax_{i}",
                            )
                        with c4:
                            antibiotics[i]["ec50"] = st.number_input(
                                "Half-max conc (EC50)",
                                value=float(antibiotics[i]["ec50"]),
                                step=0.05,
                                key=f"abx_ec50_{i}",
                            )
                            
                    antibiotics[i]["hill"] = st.number_input(
                        "Hill coefficient (H)",
                        value=float(antibiotics[i].get("hill", 1.5)),
                        step=0.1,
                        key=f"abx_hill_{i}",
                    )
                    
                    # Advanced PD inoculum + lytic
                    st.markdown("**Advanced PD features**")
                    antibiotics[i]["f_lyse"] = st.slider(
                        "Bacteriolytic Fraction (f_lyse)",
                        min_value=0.0,
                        max_value=1.0,
                        value=float(antibiotics[i].get("f_lyse", 0.0)),
                        step=0.1,
                        key=f"abx_flyse_{i}",
                        help="0.0 = non-lytic (e.g. Cipro); 1.0 = fully lytic (e.g. beta-lactams). Affects OD debris."
                    )
                    antibiotics[i]["inoculum_effect_constant"] = st.number_input(
                        "Inoculum Effect Constant (K_inoc)",
                        value=float(antibiotics[i].get("inoculum_effect_constant", 0.0)),
                        format="%.1e",
                        key=f"abx_inoc_{i}",
                        help="Bacterial density (cells/mL) at which antibiotic EC50 doubles. Set to 0 to disable."
                    )
                    
                    antibiotics[i]["Km_elim"] = st.number_input(
                        "Nonlinear elimination Km",
                        value=float(antibiotics[i].get("Km_elim", 0.0)),
                        format="%.1e",
                        key=f"abx_km_elim_{i}",
                        help="Michaelis-Menten central elimination saturation. Set to 0 to disable."
                    )

        # Host Immunity
        with col2:
            st.markdown("### Host Immunity")

            st.session_state["int_immunity_enabled"] = st.checkbox(
                "Enable Immune System Module",
                value=st.session_state.get("int_immunity_enabled", False),
            )

            if st.session_state["int_immunity_enabled"]:
                _valid_modules = ["innate", "hill"]
                _cur_module = st.session_state.get("int_immune_module", "innate")
                if _cur_module not in _valid_modules:
                    _cur_module = "innate"
                st.session_state["int_immune_module"] = st.selectbox(
                    "Immune module",
                    _valid_modules,
                    index=_valid_modules.index(_cur_module),
                    help=(
                        "innate: Imm is integrated (dImm/dt = stim_rate·B/(stim50+B) "
                        "− decay·Imm) and killing ∝ imm_kill_rate·Imm.  |  "
                        "hill: Imm is frozen; killing = imm_max·B/(imm_kill50+B) directly "
                        "(imm_stim_rate / imm_kill_rate / imm_decay_rate / initial_Imm unused)."
                    ),
                )
                _module = st.session_state["int_immune_module"]

                if _module == "innate":
                    imm_col1, imm_col2 = st.columns(2)
                    with imm_col1:
                        st.session_state["int_imm_stim_rate"] = st.number_input(
                            "Stimulation rate (imm_stim_rate)",
                            value=float(st.session_state.get("int_imm_stim_rate", 0.1)),
                            format="%.2e",
                            help="Rate at which each bacterium recruits immune effectors.",
                        )
                        st.session_state["int_innate_kill_rate"] = st.number_input(
                            "Kill rate coefficient (imm_kill_rate)",
                            value=float(st.session_state.get("int_innate_kill_rate", 1e7)),
                            format="%.1e",
                            help="Per-bacterium immune killing coefficient.",
                        )
                        st.session_state["int_innate_decay_rate"] = st.number_input(
                            "Effector decay rate (imm_decay_rate)",
                            value=float(st.session_state.get("int_innate_decay_rate", 0.1)),
                            step=0.01,
                        )
                    with imm_col2:
                        st.session_state["int_imm_stim50"] = st.number_input(
                            "Stimulation half-sat. (imm_stim50)",
                            value=float(st.session_state.get("int_imm_stim50", 1e6)),
                            format="%.1e",
                            help="Bacterial density at half-max immune stimulation.",
                        )
                        st.session_state["int_innate_kill50"] = st.number_input(
                            "Killing half-sat. (imm_kill50)",
                            value=float(st.session_state.get("int_innate_kill50", 1e5)),
                            format="%.1e",
                            help="Bacterial density at half-max immune killing.",
                        )
                        st.session_state["int_imm_initial"] = st.number_input(
                            "Initial immune density (initial_Imm)",
                            value=float(st.session_state.get("int_imm_initial", 0.0)),
                            format="%.1e",
                            help="Starting immune effector level. Typically 0 — grows from bacterial stimulation.",
                        )
                else:  # hill
                    st.caption(
                        "**Hill module:** immune killing = `imm_max · B / (imm_kill50 + B_total)`. "
                        "The `Imm` state is *frozen* (not integrated), so `imm_stim_rate`, "
                        "`imm_kill_rate`, `imm_decay_rate` and `initial_Imm` have no effect here — "
                        "only the two parameters below (plus `imm_kill_rate_D`) apply."
                    )
                    imm_col1, imm_col2 = st.columns(2)
                    with imm_col1:
                        st.session_state["int_innate_max"] = st.number_input(
                            "Max clearance rate (imm_max)",
                            value=float(st.session_state.get("int_innate_max", 1e7)),
                            format="%.1e",
                            help=(
                                "Maximum immune clearance flux (cells/mL/h) at saturating bacterial "
                                "density: killing = imm_max·B/(imm_kill50+B). Set comparable to the "
                                "bacterial growth flux (~growth_rate × bacterial load)."
                            ),
                        )
                    with imm_col2:
                        st.session_state["int_innate_kill50"] = st.number_input(
                            "Killing half-sat. (imm_kill50)",
                            value=float(st.session_state.get("int_innate_kill50", 1e5)),
                            format="%.1e",
                            help="Bacterial density at half-max immune killing.",
                        )

                st.session_state["int_imm_kill_rate_D"] = st.number_input(
                    "Kill rate for dormant/hibernating cells (imm_kill_rate_D)",
                    value=float(st.session_state.get("int_imm_kill_rate_D", 0.0)),
                    format="%.1e",
                    help="Set > 0 to allow immune clearance of dormant compartments.",
                )

                # Immune-refuge warning: dormant cells are immune-privileged unless
                # imm_kill_rate_D > 0. With dormancy on, a resistant/persister
                # population can hide in the dormant compartment — the immune system
                # neither kills it (imm_kill_rate_D=0) nor is stimulated by it
                # (dormant cells are excluded from the immune signal). The infection
                # then never clears even though immunity is "on", which is a common
                # source of confusion.
                _dormancy_on = any(
                    s.get("dormancy_enabled", False)
                    for s in st.session_state.get("int_strains", [])
                )
                if _dormancy_on and st.session_state["int_imm_kill_rate_D"] <= 0:
                    st.warning(
                        "**Dormancy + immunity:** dormant/hibernating cells are "
                        "immune-privileged while `imm_kill_rate_D = 0` — the immune "
                        "system will not kill them and they do not stimulate it. A "
                        "phage-resistant (or persister) population can survive in the "
                        "dormant reservoir, so the infection may never clear and the "
                        "resistant fraction can stay near 100% even with immunity "
                        "enabled. Set `imm_kill_rate_D > 0` above to let immunity "
                        "clear dormant cells."
                    )

    # ──── Tab 3: Environment & Dosing ─────────────────────────────────────────
    with config_tabs[2]:
        col1, col2 = st.columns(2)

        # Environment & Debris
        with col1:
            st.markdown("### Nutrient environment")
            st.caption("The **growth signal function** (and its Monod constant / carrying capacity) "
                       "is set under **Strains & Phages → Growth model**. These are the medium/reactor "
                       "conditions for nutrient-tracking growth.")
            if st.session_state.get("int_track_nutrients", True):
                st.session_state["int_initial_S"] = st.number_input(
                    "Initial Resource Substrate (S0)",
                    value=float(st.session_state.get("int_initial_S", 1.0)), step=0.1,
                )
                st.session_state["int_recycle_fraction"] = st.number_input(
                    "Nutrient recycling fraction",
                    value=float(st.session_state.get("int_recycle_fraction", 0.0)),
                    min_value=0.0, max_value=1.0, step=0.1,
                )
                st.session_state["int_s_in"] = st.number_input(
                    "Continuous medium inflow (s_in)",
                    value=float(st.session_state.get("int_s_in", 0.0)), step=0.1,
                )
                st.session_state["int_s_out"] = st.number_input(
                    "Continuous Washout dilution (s_out)",
                    value=float(st.session_state.get("int_s_out", 0.0)), step=0.05,
                )
                st.caption("**Infected-cell nutrient consumption** is now a per-phage property "
                           "(a phage×host interaction) — set it on each phage, under *Strains & "
                           "Phages → Phage → Dormant adsorption attenuation*.")
            else:
                st.info("The selected growth signal is nutrient-independent (constant / density), "
                        "so nutrient substrate dynamics are inactive.")

            st.markdown("### Optical Density (OD) & Debris")
            st.session_state["int_debris_enabled"] = st.checkbox(
                "Track Bacteriolytic Cell Debris",
                value=st.session_state.get("int_debris_enabled", False),
            )

            if st.session_state["int_debris_enabled"]:
                st.session_state["int_debris_u"] = st.number_input(
                    "Scattering weight for intact dead cells (u)",
                    value=float(st.session_state.get("int_debris_u", 0.4)),
                    step=0.1,
                )
                st.session_state["int_debris_v"] = st.number_input(
                    "Scattering weight for lysed cell fragments (v)",
                    value=float(st.session_state.get("int_debris_v", 0.2)),
                    step=0.1,
                )
                st.session_state["int_debris_kdis"] = st.number_input(
                    "Debris dissolution rate (k_dis)",
                    value=float(st.session_state.get("int_debris_kdis", 0.01)),
                    step=0.05,
                )
                st.session_state["int_dormant_od_fraction"] = st.number_input(
                    "Dormant OD contribution fraction",
                    min_value=0.0, max_value=1.0,
                    value=float(st.session_state.get("int_dormant_od_fraction", 1.0)),
                    step=0.1,
                    help="Weight of dormant (D+H) cells in OD: 1.0 = counts like active cells "
                         "(default); lower if dormant cells scatter less light.",
                )
                st.session_state["int_od_to_cfu_conversion_factor"] = st.number_input(
                    "OD-to-CFU conversion factor",
                    value=float(st.session_state.get("int_od_to_cfu_conversion_factor", 2e8)),
                    format="%.1e",
                )
                # Nutrient-dependent OD absorptivity (Geng et al. 2024): live cells shrink as
                # nutrients deplete, so OD-per-cell falls into stationary phase (OD declines even
                # at constant CFU). OFF by default (equal asymptotes = constant OD/CFU).
                st.session_state["int_od_nutrient_enabled"] = st.checkbox(
                    "Nutrient-dependent OD (cell-size shrinkage)",
                    value=bool(st.session_state.get("int_od_nutrient_enabled", False)),
                    key="widget_od_nutrient_enabled",
                    help="Scale OD-per-cell by a Hill function of nutrient S: eps(S) goes from the "
                         "stationary-phase absorptivity (low S) to the exponential-phase absorptivity "
                         "(high S). Needs a nutrient-tracking growth signal. Equal asymptotes = off.")
                if st.session_state["int_od_nutrient_enabled"]:
                    _oc1, _oc2 = st.columns(2)
                    with _oc1:
                        st.session_state["int_od_abs_exp"] = st.number_input(
                            "Exponential-phase absorptivity", min_value=0.0,
                            value=float(st.session_state.get("int_od_abs_exp", 1.0)), step=0.1, format="%g",
                            help="OD-per-cell at high nutrient (exponential growth). Reference = 1.0.")
                        st.session_state["int_od_abs_S50"] = st.number_input(
                            "OD transition nutrient level (S50)", min_value=0.0,
                            value=float(st.session_state.get("int_od_abs_S50", 0.3)), step=0.05, format="%g",
                            help="Nutrient S at the half-transition of eps(S).")
                    with _oc2:
                        st.session_state["int_od_abs_stat"] = st.number_input(
                            "Stationary-phase absorptivity", min_value=0.0,
                            value=float(st.session_state.get("int_od_abs_stat", 0.4)), step=0.1, format="%g",
                            help="OD-per-cell at low nutrient (stationary). < exponential → OD falls as "
                                 "nutrients deplete even at constant CFU.")
                        st.session_state["int_od_abs_hill"] = st.number_input(
                            "OD transition sharpness (Hill)", min_value=0.1,
                            value=float(st.session_state.get("int_od_abs_hill", 2.0)), step=0.5, format="%g",
                            help="Hill exponent of the eps(S) sigmoid.")
                    if not growth_tracks_nutrients():
                        st.caption("⚠ Nutrient-dependent OD needs a nutrient-tracking growth signal "
                                   "(S must vary) — otherwise eps(S) is constant.")

        # Dosing Schedule
        with col2:
            st.markdown("### Dosing Schedule & Regimens")

            sub_col1, sub_col2 = st.columns(2)

            with sub_col1:
                st.markdown("#### Active Dosing Events")
                if doses:
                    st.caption("Edit time / amount / route inline; :material/delete: removes a row. "
                               "To change the target, delete and re-add.")
                _routes = ["bolus", "infusion"]
                for idx, dose in enumerate(doses):
                    dose.setdefault("_id", _next_uid("dose"))  # stable key across reorder/delete
                    _did = dose["_id"]
                    c_t1, c_t2, c_t3, c_t4 = st.columns([2, 3, 3, 1])
                    with c_t1:
                        dose["time"] = st.number_input(
                            "Time (h)", min_value=0.0, value=float(dose.get("time", 0.0)),
                            step=1.0, key=f"dose_time_{_did}", label_visibility="collapsed")
                    with c_t2:
                        dose["amount"] = st.number_input(
                            f"Amount → {dose['target_type']} #{dose['target_idx']}",
                            min_value=0.0, value=float(dose.get("amount", 0.0)), format="%.2e",
                            key=f"dose_amt_{_did}", label_visibility="collapsed")
                        st.caption(f"→ {dose['target_type']} #{dose['target_idx']}")
                    with c_t3:
                        dose["route"] = st.selectbox(
                            "Route", _routes, index=_routes.index(dose.get("route", "bolus")),
                            key=f"dose_route_{_did}", label_visibility="collapsed")
                        if dose["route"] == "infusion":
                            dose["duration"] = st.number_input(
                                "Infusion duration (h)", min_value=0.1,
                                value=float(dose.get("duration") or 2.0), step=0.5,
                                key=f"dose_dur_{_did}")
                    with c_t4:
                        if st.button(":material/delete:", key=f"del_dose_{_did}"):
                            doses.pop(idx)
                            st.session_state.int_doses = doses
                            st.rerun()
                st.session_state.int_doses = doses  # persist inline edits

                st.markdown("#### Add Single Dose Event")
                with st.expander("+ Define Single Dosing Event"):
                    # Build target options
                    target_ops = ["phage"]
                    if len(antibiotics) > 0:
                        target_ops.append("antibiotic")
                    target_ops.append("nutrient")

                    d_type = st.selectbox("Target Compartment", target_ops)

                    d_time = st.number_input("Time (hours)", min_value=0.0, value=0.0, step=1.0)
                    # Target-appropriate default amount (1e8 makes sense only for phage).
                    # Keyed by target so switching target suggests a sensible default.
                    d_amount = st.number_input(
                        DOSE_AMOUNT_LABELS[d_type],
                        min_value=0.0,
                        value=DOSE_AMOUNT_DEFAULTS[d_type],
                        format="%.1e",
                        key=f"single_dose_amount_{d_type}",
                        help="Phage in PFU (e.g. 1e8); antibiotic in mg (e.g. 10); nutrient in resource units.",
                    )

                    d_idx = 0
                    if d_type == "phage" and len(phages) > 1:
                        d_idx = st.selectbox("Phage Target index", list(range(len(phages))))
                    elif d_type == "antibiotic" and len(antibiotics) > 1:
                        d_idx = st.selectbox("Antibiotic Target index", list(range(len(antibiotics))))

                    d_route = st.selectbox("Administration Route", ["bolus", "infusion"])
                    d_dur = 0.0
                    if d_route == "infusion":
                        d_dur = st.number_input("Infusion Duration (hours)", min_value=0.1, value=2.0, step=0.5)

                    if st.button("+ Add Dose Event"):
                        doses.append(
                            {
                                "time": d_time,
                                "amount": d_amount,
                                "target_type": d_type,
                                "target_idx": d_idx,
                                "route": d_route,
                                "duration": d_dur,
                            }
                        )
                        st.session_state.int_doses = doses
                        st.success("Dose event appended successfully!")
                        st.rerun()

            with sub_col2:
                st.markdown("#### Add Repeat Dosing Regimen")
                with st.expander("+ Define Repeat Dosing Regimen"):
                    # Target options
                    target_ops_rep = ["phage"]
                    if len(antibiotics) > 0:
                        target_ops_rep.append("antibiotic")
                    target_ops_rep.append("nutrient")
                    r_type = st.selectbox("Target Compartment", target_ops_rep, key="rep_dose_type")

                    r_amount = st.number_input(
                        DOSE_AMOUNT_LABELS[r_type],
                        min_value=0.0,
                        value=DOSE_AMOUNT_DEFAULTS[r_type],
                        format="%.1e",
                        key=f"rep_dose_amount_{r_type}",
                        help="Phage in PFU (e.g. 1e8); antibiotic in mg (e.g. 10); nutrient in resource units.",
                    )
                    r_interval = st.number_input("Interdose Interval (hours)", min_value=1.0, value=12.0, step=1.0, key="rep_dose_interval")
                    r_count = st.number_input("Number of Repeats", min_value=1, value=4, step=1, key="rep_dose_count")
                    r_start = st.number_input("Start Time (hours)", min_value=0.0, value=0.0, step=1.0, key="rep_dose_start")

                    r_idx = 0
                    if r_type == "phage" and len(phages) > 1:
                        r_idx = st.selectbox("Phage Target index", list(range(len(phages))), key="rep_dose_phage_idx")
                    elif r_type == "antibiotic" and len(antibiotics) > 1:
                        r_idx = st.selectbox("Antibiotic Target index", list(range(len(antibiotics))), key="rep_dose_abx_idx")

                    r_route = st.selectbox("Administration Route", ["bolus", "infusion"], key="rep_dose_route")
                    r_dur = 0.0
                    if r_route == "infusion":
                        r_dur = st.number_input("Infusion Duration (hours)", min_value=0.1, value=2.0, step=0.5, key="rep_dose_duration")

                    if st.button("+ Add Repeat Regimen", key="rep_dose_add_btn"):
                        for k in range(int(r_count)):
                            doses.append({
                                "time": r_start + k * r_interval,
                                "amount": r_amount,
                                "target_type": r_type,
                                "target_idx": r_idx,
                                "route": r_route,
                                "duration": r_dur,
                            })
                        st.session_state.int_doses = doses
                        st.success(f"Added {r_count} repeat dose events successfully!")
                        st.rerun()


    # ──── Tab 4: Solver Settings ──────────────────────────────────────────────
    with config_tabs[3]:
        col1, col2 = st.columns(2)

        with col1:
            st.markdown("### Solver Time parameters")
            st.session_state["int_t_end"] = st.number_input(
                "Simulation end time (hours)",
                value=float(st.session_state.get("int_t_end", 48.0)),
                step=12.0,
            )
            st.session_state["int_dt"] = st.number_input(
                "ODE output interval (dt)",
                value=float(st.session_state.get("int_dt", 0.25)),
                step=0.05,
            )

            st.markdown("### Model structure")
            st.session_state["int_n_latent"] = st.number_input(
                "Latency compartments (n_latent)",
                min_value=1, max_value=20,
                value=int(st.session_state.get("int_n_latent", 5)),
                help="Number of latent-infection stages (Erlang-distributed latent period). "
                     "Applies to all builder modes (Direct / BRG / Custom Strains).",
            )

            st.markdown("### Advanced solver options")
            st.session_state["int_superinfection"] = st.checkbox(
                "Allow Phage Superinfection",
                value=st.session_state.get("int_superinfection", False),
                help="If active, phages will adsorb to latently-infected cells without producing extra burst."
            )
            
            st.session_state["int_t_prerun"] = st.number_input(
                "Pre-treatment prerun pre-growth (hours)",
                value=float(st.session_state.get("int_t_prerun", 0.0)),
                step=4.0,
                help="Let the bacteria/nutrient system equilibrate without treatments before t=0. Set to 0 to disable."
            )
            if st.session_state.get("int_t_prerun", 0.0) > 0 and st.session_state.get("int_debris_enabled", False):
                st.session_state["int_prerun_inherit_debris"] = st.checkbox(
                    "Inherit bacterial debris from the pre-run",
                    value=bool(st.session_state.get("int_prerun_inherit_debris", True)),
                    key="widget_prerun_inherit_debris",
                    help="On (default): the OD/debris that built up during the pre-run carries into "
                         "treatment (t=0 starts with a realistic dead-cell background). Off: the dead "
                         "cells are washed out — treatment starts with zero debris.")

        with col2:
            st.markdown("### ODE solver specifics")
            st.session_state["int_extinction_threshold"] = st.number_input(
                "Absorbing Extinction threshold",
                value=float(st.session_state.get("int_extinction_threshold", 1.0)),
                step=1.0,
                help="If density falls below this threshold, it is locked to 0 to prevent numerical recovery.",
            )
            st.session_state["int_extinction_check_interval"] = st.number_input(
                "Extinction check interval (hours)",
                value=float(st.session_state.get("int_extinction_check_interval", 0.0)),
                min_value=0.0,
                step=1.0,
                help=(
                    "Apply the extinction threshold every this many hours (not just at "
                    "dose events). Zeroes any sub-threshold strain before it can regrow "
                    "from a below-threshold pool. Set to 0 to check only at dose boundaries. "
                    "Ignored when the extinction threshold is 0."
                ),
            )
            st.session_state["int_solver_method"] = st.selectbox(
                "ODE Solver integration method",
                ["BDF", "Radau", "LSODA"],
                index=["BDF", "Radau", "LSODA"].index(
                    st.session_state.get("int_solver_method", "BDF")
                ),
            )
            st.session_state["int_phage_noise_floor"] = st.number_input(
                "Phage noise floor suppression",
                value=float(st.session_state.get("int_phage_noise_floor", 1e-9)),
                help="Suppress floating-point round-off of non-dosed phages. Set <= 0 to disable.",
            )

    # ──── Run Button ──────────────────────────────────────────────────────────
    st.markdown("<br>", unsafe_allow_html=True)
    if st.button("Run Simulation", width="stretch", type="primary"):
        with st.spinner("Assembling model equations & integrating..."):
            try:
                _pc_cfg, _pc_B0, _pc_iP, _pc_iS, *_ = build_nominal_config_from_gui()
                warn_if_prerun_collapses(_pc_cfg, _pc_B0, initial_S=_pc_iS)
                _t0 = time.perf_counter()
                result, config = run_sim_from_gui_params()
                st.session_state.sim_runtime = time.perf_counter() - _t0
                st.session_state.simulation_result = result
                st.session_state.simulation_config = config
                st.success("Simulation finished successfully!")
            except Exception as e:
                st.error(f"Solver Error: {e}")
                import traceback
                st.code(traceback.format_exc())

    # ──── Outputs / Results Section ───────────────────────────────────────────
    if st.session_state.simulation_result is not None:
        result = st.session_state.simulation_result
        config = st.session_state.simulation_config

        # 1. Calculate Metrics. "Total bacteria" here follows the user-chosen basis
        # (default culturable CFU = B+D); the selector lives in the Bacterial Dynamics tab.
        _metric_basis = st.session_state.get("sim_metric_basis", CFU_BASIS_LABEL)
        total_bacteria = bacterial_total(result, _metric_basis)
        nadir_val = np.min(total_bacteria)
        nadir_time = result.time[np.argmin(total_bacteria)]

        _trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
        auc_val = _trapz(total_bacteria, result.time)

        # Clearance & Log Reduction
        t_clear = time_to_clearance(
            result,
            threshold=st.session_state.get("int_extinction_threshold", 1.0),
        )
        t_log_red = time_to_log_reduction(result, n_logs=2.0)

        # Peak free-phage titre
        _phage_tot = np.asarray(result.sum_prefixes("P"), dtype=float)
        peak_phage = float(_phage_tot.max()) if _phage_tot.size else 0.0
        peak_phage_t = float(result.time[int(np.argmax(_phage_tot))]) if _phage_tot.size else 0.0

        # Outcome classification for the results header badge
        _b0 = float(total_bacteria[0]) if len(total_bacteria) else 0.0
        _bend = float(total_bacteria[-1]) if len(total_bacteria) else 0.0
        if t_clear is not None:
            _outcome, _obg, _ofg = "Cleared", "var(--teal)", "#fff"
        elif _b0 > 0 and _bend <= _b0 * 0.1:
            _outcome, _obg, _ofg = "Suppressed", "var(--teal-tint)", "var(--teal)"
        elif _b0 > 0 and nadir_val <= _b0 * 0.1 and _bend > nadir_val * 10:
            _outcome, _obg, _ofg = "Regrowth", "#f3e4cf", "#8a5a1a"
        else:
            _outcome, _obg, _ofg = "Uncontrolled", "#f4dedb", "#9b3b33"

        # Results header bar: title + solver/runtime meta + outcome badge
        _rt = st.session_state.get("sim_runtime")
        _meta = f"t = 0–{result.time[-1]:.0f} h · {st.session_state.get('int_solver_method', 'BDF')}"
        if _rt is not None:
            _meta += f" · {_rt:.2f} s"
        st.markdown(
            "<div style='display:flex;align-items:center;justify-content:space-between;margin:2px 0 14px'>"
            "<div><div style='font-size:1.25rem;font-weight:600;color:var(--ink)'>Simulation results</div>"
            f"<div class='section-label' style='margin-top:3px'>{_meta}</div></div>"
            f"<div style='background:{_obg};color:{_ofg};padding:6px 14px;border-radius:6px;"
            f"font-weight:600;font-size:13px'>{_outcome}</div>"
            "</div>",
            unsafe_allow_html=True,
        )

        # Render Metrics in Columns
        m_col1, m_col2, m_col3, m_col4, m_col5, m_col6 = st.columns(6)
        with m_col1:
            st.markdown(
                f"""
                <div class="metric-container">
                    <div class="metric-label">Bacterial Nadir</div>
                    <div class="metric-value">{nadir_val:.2e}</div>
                    <div class="metric-sub">cells/mL at t={nadir_time:.1f}h</div>
                </div>
                """,
                unsafe_allow_html=True,
            )
        with m_col2:
            st.markdown(
                f"""
                <div class="metric-container">
                    <div class="metric-label">Bacterial AUC</div>
                    <div class="metric-value">{auc_val:.1e}</div>
                    <div class="metric-sub">cells·h/mL</div>
                </div>
                """,
                unsafe_allow_html=True,
            )
        with m_col3:
            clear_lbl = f"{t_clear:.1f}h" if t_clear is not None else "Never"
            st.markdown(
                f"""
                <div class="metric-container">
                    <div class="metric-label">Time to Clearance</div>
                    <div class="metric-value">{clear_lbl}</div>
                    <div class="metric-sub">threshold={st.session_state.get('int_extinction_threshold', 1.0)}</div>
                </div>
                """,
                unsafe_allow_html=True,
            )
        with m_col4:
            log_lbl = f"{t_log_red:.1f}h" if t_log_red is not None else "Never"
            st.markdown(
                f"""
                <div class="metric-container">
                    <div class="metric-label">2-Log Reduction Time</div>
                    <div class="metric-value">{log_lbl}</div>
                    <div class="metric-sub">from baseline B(0)</div>
                </div>
                """,
                unsafe_allow_html=True,
            )
        with m_col5:
            # show final resistant fraction
            total_active = result.sum_prefixes("B")
            final_res_frac = (total_active[-1] - result.get("B0")[-1]) / total_active[-1] if (total_active[-1] > 0 and len(strains) > 1) else 0.0
            st.markdown(
                f"""
                <div class="metric-container">
                    <div class="metric-label">Resistant fraction</div>
                    <div class="metric-value">{final_res_frac*100:.1f}%</div>
                    <div class="metric-sub">at t_end={result.time[-1]}h</div>
                </div>
                """,
                unsafe_allow_html=True,
            )
        with m_col6:
            st.markdown(
                f"""
                <div class="metric-container">
                    <div class="metric-label">Peak Phage Titre</div>
                    <div class="metric-value">{peak_phage:.2e}</div>
                    <div class="metric-sub">PFU/mL at t={peak_phage_t:.1f}h</div>
                </div>
                """,
                unsafe_allow_html=True,
            )

        st.markdown("<br>", unsafe_allow_html=True)

        # 2. Plot Tabs
        plot_tabs = st.tabs(
            [
                "Bacterial Dynamics",
                "Phage Dynamics",
                "Nutrients & OD",
                "Antibiotics & Immunity",
            ]
        )

        t = result.time

        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        def _sim_layout(fig, title, ylab):
            fig.update_layout(title=title, xaxis_title="Time (hours)", yaxis_title=ylab,
                              template="plotly_white", height=440, margin=dict(t=48, b=40),
                              legend=dict(orientation="h", yanchor="bottom", y=-0.28, x=0))
            return fig

        # Plottable-series registry (single source of truth for entity choices).
        _series = build_series(result, config=config, strains=strains, phages=phages,
                               antibiotics=antibiotics, builder_mode=builder_mode)

        # Bacterial Dynamics Plot — pick which compartments to show.
        with plot_tabs[0]:
            _bact = [s for s in _series if s.category == "Bacteria"]
            st.selectbox(
                "Outcome metric basis", list(BACTERIAL_BASES.keys()), index=0,
                key="sim_metric_basis",
                help="Which compartments count as 'bacteria' for the metrics and outcome badge "
                     "above. CFU = culturable (active B + dormant D); infected (I) and "
                     "hibernating (H) cells don't form colonies on a plate.")
            st.caption("Series to plot:")
            _chosen = series_selector("sim_bact_sel", _bact)
            if _chosen:
                fig = plot_series(result, _chosen, plot_axis_controls("sim_bact", default_y="Log"),
                                  title="Bacterial Population Trajectories", y_label="Density (cells/mL)")
                st.plotly_chart(fig, width="stretch")
            else:
                st.info("Select at least one series to plot.")

        # Phage Dynamics Plot — pick which phage series to show.
        with plot_tabs[1]:
            _phg = [s for s in _series if s.category == "Phage"]
            if _phg:
                st.caption("Series to plot:")
                _chosen = series_selector("sim_phage_sel", _phg)
                if _chosen:
                    fig = plot_series(result, _chosen, plot_axis_controls("sim_phage", default_y="Log"),
                                      title="Phage Population Trajectories", y_label="Density (phages/mL)")
                    st.plotly_chart(fig, width="stretch")
                else:
                    st.info("Select at least one series to plot.")
            else:
                st.info("No phages were configured in this simulation.")

        # Nutrients & OD Plot
        with plot_tabs[2]:
            col_nut1, col_nut2 = st.columns(2)
            with col_nut1:
                if track_nutrients:
                    fig = go.Figure()
                    fig.add_trace(go.Scatter(x=t, y=result.get("S"), mode="lines",
                                             name="Substrate (S)", line=dict(color="#c1873a", width=2)))
                    _sim_layout(fig, "Nutrient Resource Depletion", "Substrate concentration")
                    st.plotly_chart(fig, width="stretch")
                else:
                    st.info("Nutrient tracking is disabled (constant/logistic growth).")
            with col_nut2:
                if st.session_state.get("int_debris_enabled", False):
                    cfu_od = total_bacteria / st.session_state.get("int_od_to_cfu_conversion_factor", 2e8)
                    fig = go.Figure()
                    fig.add_trace(go.Scatter(x=t, y=_safe_od(result, total_bacteria), mode="lines",
                                             name="OD (AU)", line=dict(color="#0d7a68", width=2.5)))
                    fig.add_trace(go.Scatter(x=t, y=cfu_od, mode="lines", name="Live-only OD",
                                             line=dict(color="#0d7a68", dash="dash"), opacity=0.6))
                    _sim_layout(fig, "Simulated Optical Density (Live + Debris)", "Optical Density (AU)")
                    st.plotly_chart(fig, width="stretch")
                else:
                    st.info("Bacterial debris & Optical Density (OD) tracking was not enabled.")

        # Antibiotics & Host Immunity Plot
        with plot_tabs[3]:
            abx_present = len(antibiotics) > 0
            imm_present = st.session_state.get("int_immunity_enabled", False)
            if abx_present or imm_present:
                fig = make_subplots(specs=[[{"secondary_y": True}]])
                if abx_present:
                    for j, abx in enumerate(antibiotics):
                        Vc = abx.get("Vc", 1.0)
                        fig.add_trace(go.Scatter(x=t, y=result.get(f"Ac{j}") / Vc, mode="lines",
                                                 name=f"{abx['name']} (Blood)", line=dict(color="#3b6fb5")),
                                      secondary_y=False)
                        if abx.get("k12", 0.0) > 0:
                            fig.add_trace(go.Scatter(x=t, y=result.get(f"Ap{j}") / Vc, mode="lines",
                                                     name=f"{abx['name']} (Peripheral)",
                                                     line=dict(color="#3b6fb5", dash="dash")),
                                          secondary_y=False)
                # When an antibiotic is present, put Imm on the secondary y-axis; otherwise on the
                # primary (so an immunity-only run is never drawn on an empty/throwaway axis).
                if imm_present:
                    fig.add_trace(go.Scatter(x=t, y=result.get("Imm"), mode="lines",
                                             name="Immune Effector", line=dict(color="#b5487f")),
                                  secondary_y=bool(abx_present))
                fig.update_xaxes(title_text="Time (hours)")
                fig.update_yaxes(title_text="Antibiotic concentration", secondary_y=False)
                if imm_present:
                    fig.update_yaxes(title_text="Immune Effector Cells (Imm)", secondary_y=bool(abx_present))
                fig.update_layout(title="Pharmacokinetics & Host Defense Dynamics",
                                  template="plotly_white", height=440, margin=dict(t=48, b=40),
                                  legend=dict(orientation="h", yanchor="bottom", y=-0.28, x=0))
                st.plotly_chart(fig, width="stretch")
            else:
                st.info("No antibiotics or immune modules were configured.")

        # 3. Export Code & Data
        st.markdown("### Export & Reproducibility")

        c_down1, c_down2 = st.columns(2)
        with c_down1:
            # Generate CSV data
            df_export = pd.DataFrame({"time": t})
            for col in result.state_names:
                df_export[col] = result.get(col)
            df_export["CFU_total"] = total_bacteria

            csv_buffer = io.StringIO()
            df_export.to_csv(csv_buffer, index=False)
            csv_str = csv_buffer.getvalue()

            st.download_button(
                "Download Simulation Trajectories (CSV)",
                data=csv_str,
                file_name="pbisim_simulation_results.csv",
                mime="text/csv",
                width="stretch",
            )

        with c_down2:
            rep_code = generate_reproduction_code()
            st.download_button(
                "Download Python Script",
                data=rep_code,
                file_name="pbisim_run.py",
                mime="text/x-python",
                width="stretch",
            )

        with st.expander("View Python Reproduction Code"):
            st.code(rep_code, language="python")
