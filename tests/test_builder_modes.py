"""Tests for pbisim-app builder modes (BRG, StrainSet) and cohort trials."""

from __future__ import annotations

import numpy as np
import pytest
from pbisim import PBIModel, solve_ode
from pbisim.strains import StrainDefinition, StrainSet
from pbisim.strains.genotypes import BinaryResistanceGenotypes, BacterialStrain, PhageStrain, Antibiotic
from pathlib import Path as _Path

# Absolute path — Streamlit >=1.59 resolves a relative AppTest.from_file() path against
# the CALLER's directory (tests/), not the CWD, so a relative "pbisim_app/app.py" breaks.
APP = str(_Path(__file__).resolve().parents[1] / "pbisim_app" / "app.py")


def test_brg_build_and_run():
    """Verify BinaryResistanceGenotypes constructs config and runs successfully."""
    # 1. Setup BRG
    b = BacterialStrain(base_growth_rate=1.2)
    p = [PhageStrain(name="Phage0", burst_size_s=50.0)]
    abx = [Antibiotic(name="Drug0", emax_s=3.0, ec50_s=0.2, emax_r=0.1, ec50_r=2.0)]
    
    brg = BinaryResistanceGenotypes.from_strains(p, bacteria=b, antibiotics=abx)
    config = brg.to_config(n_latent=5, n_depth=1)
    
    # 2 phages/abx loci -> 4 genotypes: '00', '10', '01', '11'
    initial_B = np.array([1e7, 10.0, 10.0, 0.0])
    initial_P = np.array([1e6])
    initial_S = 1.0
    
    model = PBIModel(config, initial_B=initial_B, initial_P=initial_P, initial_S=initial_S)
    result = solve_ode(model, t_end=24.0, dt=1.0)
    assert result is not None
    assert len(result.time) > 0


def test_strainset_build_and_run():
    """Verify StrainSet with graph mutation dictionary constructs config and runs successfully."""
    # 1. Setup StrainSet
    ss = StrainSet(n_phages=1)
    from pbisim.pk.antibiotic import AntibioticDefinition, AntibioticSensitivity
    ss.add_antibiotic(AntibioticDefinition("cipro", k_elim=0.5))
    
    ss.add_strain(StrainDefinition(
        name="WT",
        growth_rate=1.2,
        adsorption_rates=np.array([1e-7]),
        adsorption_rates_dormant=np.array([0.0]),
        burst_sizes=np.array([50.0]),
        latent_periods=np.array([0.5]),
        latent_periods_dormant=np.array([0.5]),
        bacteria_to_resource_ratio=1e9,
        dormancy_rate=0.0,
        resuscitation_rate=0.0,
        dormancy_diffusion_rate=0.0,
        imm_stim_rate=0.0,
        imm_kill_rate=0.0,
        attenuation_rate=np.zeros(1),
        antibiotic_sensitivity={"cipro": AntibioticSensitivity(emax=5.0, ec50=0.5)}
    ))
    
    ss.add_strain(StrainDefinition(
        name="R_cipro",
        growth_rate=1.1,
        adsorption_rates=np.array([1e-7]),
        adsorption_rates_dormant=np.array([0.0]),
        burst_sizes=np.array([50.0]),
        latent_periods=np.array([0.5]),
        latent_periods_dormant=np.array([0.5]),
        bacteria_to_resource_ratio=1e9,
        dormancy_rate=0.0,
        resuscitation_rate=0.0,
        dormancy_diffusion_rate=0.0,
        imm_stim_rate=0.0,
        imm_kill_rate=0.0,
        attenuation_rate=np.zeros(1),
        antibiotic_sensitivity={"cipro": AntibioticSensitivity(emax=0.5, ec50=5.0)}
    ))
    
    ss.set_mutation_graph({"WT": {"R_cipro": 1e-7}})
    config = ss.to_config(
        n_latent=3,
        n_depth=1,
        phage_decay_rates=np.array([0.03]),
        imm_decay_rate=0.1,
        imm_stim50=1e6,
        imm_kill50=1e6,
        monod_constant=1.0,
        recycle_fraction=0.5
    )
    
    initial_B = np.array([1e7, 10.0])
    initial_P = np.array([1e6])
    initial_S = 1.0
    
    model = PBIModel(config, initial_B=initial_B, initial_P=initial_P, initial_S=initial_S)
    result = solve_ode(model, t_end=24.0, dt=1.0)
    assert result is not None
    assert len(result.time) > 0


def test_clinical_trial_with_pretreatment():
    """Verify that a ClinicalTrial running with a PretreatmentPhase succeeds with InitialConditions attached."""
    from pbisim.trial.clinical import TreatmentArm
    from pbisim.trial.population import InitialConditions
    from pbisim.builder import ModelBuilder
    from pbisim_app.trial_helper import run_trial_simulation
    
    cfg = ModelBuilder(n_bacteria=1, n_phages=0).with_growth_rates(1.0).build()
    
    init_B = np.array([1e6])
    init_P = np.zeros(0)
    init_S = 1.0
    
    cfg.initial_conditions = InitialConditions(
        B=init_B,
        P=init_P,
        S=init_S
    )
    
    iiv_inputs = [
        {"path": "growth_rates", "dist_type": "LogNormal", "params": {"cv": 0.2}, "mode": "multiplicative"}
    ]
    
    arms = [TreatmentArm(name="Control", dose_schedule=None)]
    
    res = run_trial_simulation(
        cfg,
        iiv_inputs,
        arms,
        n_patients=2,
        t_end=5.0,
        dt=1.0,
        seed=42,
        pretreatment_hours=2.0,
        n_jobs=1,
        base_initial_B=init_B,
        base_initial_P=init_P,
        base_initial_S=init_S
    )
    
    assert res is not None
    assert "Control" in res.arm_names


def test_trial_pretreatment_carries_dormant_reservoir():
    """A trial PretreatmentPhase must carry its dormant reservoir into treatment.

    Same class of bug as the interactive-simulator pre-run: PretreatmentPhase replaces
    the patient config's initial_conditions with the full stationary-phase state, but
    the model factory previously took initial_D/initial_Imm from the GUI base kwargs,
    discarding the (dominant) dormant population. A long pretreatment on a dormancy
    model must still start the treatment near the stationary carrying capacity.
    """
    from pbisim.trial.clinical import TreatmentArm
    from pbisim.trial.population import InitialConditions
    from pbisim.builder import ModelBuilder
    from pbisim_app.trial_helper import run_trial_simulation

    b = ModelBuilder(n_bacteria=1, n_phages=0, n_latent=5, n_depth=3)
    b = b.with_growth_rates([1.2], bacteria_to_resource_ratio=[1e9])
    b = b.with_dormancy(dormancy_rate=np.array([0.2]), resuscitation_rate=np.array([0.1]),
                        dormancy_diffusion_rate=np.array([0.05]))
    b = b.with_nutrient(track_nutrients=True, monod_constant=0.3)
    cfg = b.build()

    init_B = np.array([1e7])
    init_P = np.zeros(0)
    init_S = 1.0
    cfg.initial_conditions = InitialConditions(B=init_B, P=init_P, S=init_S)

    res = run_trial_simulation(
        cfg, [], [TreatmentArm(name="Control", dose_schedule=None)],
        n_patients=2, t_end=12.0, dt=0.5, seed=1,
        pretreatment_hours=48.0,  # long prerun -> population mostly dormant
        n_jobs=1, base_initial_B=init_B, base_initial_P=init_P, base_initial_S=init_S,
    )

    for r in res["Control"].results:
        assert r is not None
        total = r.sum_prefixes("B", "D", "I", "H")
        assert total[0] > 1e8, "treatment must start from the full stationary population"


def test_trial_control_arm_has_no_phage():
    """Regression: the phage inoculum must not leak into non-phage arms.

    Crossover trials share initial conditions across arms and differ only by
    dose schedule, so seeding the phage via ``initial_P`` would contaminate the
    Control (and Antibiotic-Only) arms and eradicate the bacteria there, making
    the control indistinguishable from phage therapy.  The app instead starts
    every arm at zero free phage and delivers the inoculum as a t=0 bolus only
    in phage-containing arms.  This test locks in that invariant: Control keeps
    phage at zero and bacteria grow, while the phage arm eradicates them.
    """
    from pbisim.trial.clinical import TreatmentArm
    from pbisim.trial.population import InitialConditions
    from pbisim.builder import ModelBuilder
    from pbisim.pk.dosing import DoseEvent, DoseSchedule
    from pbisim_app.trial_helper import run_trial_simulation

    builder = ModelBuilder(n_bacteria=1, n_phages=1, n_latent=5, n_depth=1)
    builder = builder.with_growth_rates([0.6], bacteria_to_resource_ratio=[1e9])
    builder = builder.with_phage_params(
        adsorption_rates=np.array([[1e-8]]),
        adsorption_rates_dormant=np.array([[0.0]]),
        burst_sizes=np.array([[50.0]]),
        latent_periods=np.array([[0.5]]),
        phage_decay_rates=np.array([0.1]),
    )
    builder = builder.with_nutrient(track_nutrients=True, monod_constant=0.3)
    base_cfg = builder.build()

    init_B = np.array([1e7])
    init_P = np.array([1e6])
    init_S = 1.0

    # Replicate the fixed app arm-assembly: zero the shared phage baseline and
    # deliver the inoculum only as a t=0 phage bolus in the phage arm.
    base_P = np.zeros_like(init_P)
    phage_inoculum_doses = [
        DoseEvent(time=0.0, amount=float(init_P[i]), target="phage",
                  index=i, route="bolus", duration=0.0)
        for i in range(len(init_P)) if float(init_P[i]) > 0.0
    ]
    base_cfg.initial_conditions = InitialConditions(B=init_B, P=base_P, S=init_S)

    arms = [
        TreatmentArm(name="Control", dose_schedule=DoseSchedule([])),
        TreatmentArm(name="Phage-Only", dose_schedule=DoseSchedule(phage_inoculum_doses)),
    ]

    res = run_trial_simulation(
        base_cfg,
        [],  # no IIV — deterministic
        arms,
        n_patients=3,
        t_end=48.0,
        dt=0.5,
        seed=1,
        pretreatment_hours=0.0,
        n_jobs=1,
        base_initial_B=init_B,
        base_initial_P=base_P,
        base_initial_S=init_S,
    )

    # Control: no phage present at any time, bacteria grow (no eradication).
    control = res["Control"]
    for r in control.results:
        assert r is not None
        assert np.max(r.sum_prefixes("P")) < 1.0, "phage leaked into the Control arm"
        total_b = r.sum_prefixes("B", "D", "I", "H")
        assert total_b[-1] > total_b[0], "Control bacteria should grow untreated"

    # Phage-Only: phage delivered, bacteria eradicated.
    phage_arm = res["Phage-Only"]
    for r in phage_arm.results:
        assert r is not None
        assert np.max(r.sum_prefixes("P")) > init_P[0], "phage inoculum not delivered"
        total_b = r.sum_prefixes("B", "D", "I", "H")
        assert total_b[-1] < 1.0, "phage arm should eradicate bacteria"


def test_dormancy_creates_immune_refuge():
    """Dormant cells are immune-privileged unless imm_kill_rate_D > 0.

    This locks in the mechanism behind the app's dormancy+immunity warning: with
    dormancy enabled, a resistant strain survives in the dormant reservoir that
    immunity neither kills (imm_kill_rate_D=0) nor is stimulated by, so the
    infection never clears — whereas imm_kill_rate_D > 0 lets immunity clear it.
    """
    from pbisim.builder import ModelBuilder

    def build(kill_rate_D):
        b = ModelBuilder(n_bacteria=2, n_phages=1, n_latent=5, n_depth=3)
        b = b.with_growth_rates([1.2, 1.2], bacteria_to_resource_ratio=[1e9, 1e9])
        b = b.with_dormancy(
            dormancy_rate=np.array([0.2, 0.2]),
            resuscitation_rate=np.array([0.1, 0.1]),
            dormancy_diffusion_rate=np.array([0.05, 0.05]),
        )
        b = b.with_phage_params(
            adsorption_rates=np.array([[1e-8], [0.0]]),
            adsorption_rates_dormant=np.array([[0.0], [0.0]]),
            burst_sizes=np.array([[50.0], [50.0]]),
            latent_periods=np.array([[0.5], [0.5]]),
            phage_decay_rates=np.array([0.1]),
        )
        b = b.with_mutations(phage_resistance_rates=[1e-7])
        b = b.with_nutrient(track_nutrients=True, monod_constant=0.3)
        b = b.with_immunity(
            imm_stim_rate=np.full(2, 1.0), imm_stim50=1e6,
            imm_kill_rate=np.full(2, 1e7), imm_kill50=1e8,
            imm_decay_rate=0.05, immune_module="innate",
            imm_kill_rate_D=(np.array([kill_rate_D, kill_rate_D])
                             if kill_rate_D > 0 else None),
        )
        m = PBIModel(b.build(), initial_B=np.array([1e7, 0.0]),
                     initial_P=np.array([1e6]), initial_S=1.0, initial_Imm=0.0)
        return solve_ode(m, t_end=48.0, dt=0.5, method="BDF", extinction_threshold=1.0)

    refuge = build(0.0).sum_prefixes("B", "D", "I", "H")[-1]
    cleared = build(1e7).sum_prefixes("B", "D", "I", "H")[-1]

    # With no dormant killing, a large reservoir persists (>1e7); enabling it
    # drops the burden by orders of magnitude.
    assert refuge > 1e7, "dormant reservoir should persist when imm_kill_rate_D=0"
    assert cleared < refuge / 100.0, "imm_kill_rate_D>0 should clear the reservoir"


def test_hill_immune_module_active_at_app_defaults():
    """The hill module must actually kill bacteria at the app's default parameters.

    Hill killing is imm_max·B/(imm_kill50+B); with the app defaults (imm_max=1e7,
    imm_kill50=1e5) it should control an otherwise-unchecked bloom. Also verifies the
    engine ignores the innate-only fields in hill mode (imm_stim_rate/imm_kill_rate/
    imm_decay_rate), which is why the UI hides them — changing them must not move the
    hill result.
    """
    from pbisim.builder import ModelBuilder

    def build(immune, **imm):
        b = ModelBuilder(n_bacteria=1, n_phages=0, n_latent=5, n_depth=1)
        b = b.with_growth_rates([1.2], bacteria_to_resource_ratio=[1e9])
        b = b.with_nutrient(track_nutrients=True, monod_constant=0.3)
        if immune:
            b = b.with_immunity(immune_module="hill", **imm)
        m = PBIModel(b.build(), initial_B=np.array([1e7]),
                     initial_P=np.zeros(0), initial_S=1.0)
        return solve_ode(m, t_end=48.0, dt=0.5, method="BDF",
                         extinction_threshold=1.0).sum_prefixes("B", "D", "I", "H")

    defaults = dict(imm_max=1e7, imm_kill50=1e5, imm_stim_rate=np.full(1, 0.1),
                    imm_kill_rate=np.full(1, 1e7), imm_decay_rate=0.1, imm_stim50=1e6)

    off = build(False)[-1]
    on = build(True, **defaults)[-1]
    assert off > 1e8, "control should bloom to carrying capacity"
    assert on < off / 1000.0, "hill immunity must control the bloom at app defaults"

    # Innate-only fields are inert in hill mode: perturbing them by orders of
    # magnitude must not change the outcome.
    inert = {**defaults, "imm_stim_rate": np.full(1, 1e3),
             "imm_kill_rate": np.full(1, 1e13), "imm_decay_rate": 5.0}
    assert np.isclose(build(True, **inert)[-1], on, atol=1.0), \
        "hill output must not depend on imm_stim_rate/imm_kill_rate/imm_decay_rate"


def test_prerun_carries_dormant_reservoir():
    """A long stationary-phase prerun must not collapse the treatment population.

    Guards the run_sim_from_gui_params prerun fix: stationary_phase_ic returns the
    active B *and* the dormant reservoir D (which dominates at stationary phase).
    Keeping only B (the old behaviour) discards most of the culture, so a longer
    prerun starts the treatment with fewer and fewer cells until it reads as ~0.
    Carrying ic.D forward must preserve the full population regardless of prerun
    length.
    """
    from pbisim.builder import ModelBuilder
    from pbisim.analysis import stationary_phase_ic

    b = ModelBuilder(n_bacteria=1, n_phages=0, n_latent=5, n_depth=3)
    b = b.with_growth_rates([1.2], bacteria_to_resource_ratio=[1e9])
    b = b.with_dormancy(dormancy_rate=np.array([0.2]), resuscitation_rate=np.array([0.1]),
                        dormancy_diffusion_rate=np.array([0.05]))
    b = b.with_nutrient(track_nutrients=True, monod_constant=0.3)
    cfg = b.build()

    def treat(t_prerun, keep_dormant):
        ic = stationary_phase_ic(cfg, t_prerun=t_prerun, B0=np.array([1e7]))
        kw = {}
        if keep_dormant and ic.D is not None:
            kw["initial_D"] = ic.D
        m = PBIModel(cfg, initial_B=ic.B, initial_P=np.zeros(0),
                     initial_S=max(float(ic.S), 0.0), **kw)
        return solve_ode(m, t_end=24.0, dt=0.5, method="BDF",
                         extinction_threshold=1.0).sum_prefixes("B", "D", "I", "H")

    dropped = treat(48.0, keep_dormant=False)   # old behaviour
    kept = treat(48.0, keep_dormant=True)        # fixed behaviour
    assert dropped.max() < 1e7, "dropping ic.D should lose most of the stationary culture"
    assert kept.max() > 1e8, "carrying ic.D must preserve the stationary population"
    assert kept.max() > 100 * dropped.max(), "fix must retain orders of magnitude more biomass"


def test_bacteria_to_resource_ratio_editable_in_all_modes():
    """bacteria_to_resource_ratio has an editable widget in Direct + StrainSet
    (it was previously read-only there); pseudolysogeny is now editable in
    BRG + StrainSet."""
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(APP, default_timeout=200)
    at.run()
    keys = lambda a: {n.key for n in a.number_input if n.key}
    assert "str_ratio_0" in keys(at)  # Direct

    at.session_state["widget_builder_mode"] = "Custom Strains & Graph (StrainSet)"
    at.run(); at.run()
    k = keys(at)
    assert "ss_str_ratio_0" in k
    assert all(f"ss_phg_{s}_0" in k for s in ("hib_s", "hib_r", "res_s", "res_r"))

    at.session_state["widget_builder_mode"] = "Binary Genotypes (BRG)"
    at.run(); at.run()
    k = keys(at)
    assert all(f"brg_phg_{s}_0" in k for s in ("hib_s", "hib_r", "res_s", "res_r"))
    assert len(at.exception) == 0


def test_direct_mutation_rate_survives_navigation():
    """A mutation rate of 0 in Direct mode must survive navigating to the parameter
    sweep and back (the widget value= previously read an unwired key and reverted
    to 1e-7)."""
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(APP, default_timeout=180)
    at.run()
    # 2 strains, 1 phage -> the 2^m mutation shortcut renders direct_mu_0
    wt = dict(at.session_state["int_strains"][0])
    at.session_state["int_strains"] = [wt, {"name": "R", "initial_B": 0.0, "growth_rate": 1.2,
                                            "bacteria_to_resource_ratio": 1e9, "death_rate_B": 0.0,
                                            "dormancy_enabled": False}]
    at.run()
    at.session_state["direct_mu_0"] = 0.0
    at.run()
    assert at.session_state["direct_phg_res_rates"] == [0.0]
    at.session_state["current_page_radio"] = "Parameter Sweeps"
    at.run()
    at.session_state["current_page_radio"] = "Interactive Simulator"
    at.run()
    muw = [n for n in at.number_input if n.key == "direct_mu_0"]
    assert muw and muw[0].value == 0.0
    assert at.session_state["direct_phg_res_rates"] == [0.0]


def test_dormancy_signals_and_growth_signals():
    """nutrient+density dormancy no longer errors (was a string mismatch), the dormancy
    Ks flows through, and all four growth signals build."""
    from streamlit.testing.v1 import AppTest

    # nutrient+density dormancy with a custom Ks runs (the reported bug). Drive the
    # widgets (not the dict) so keyed widget state doesn't override the edits.
    at = AppTest.from_file(APP, default_timeout=200)
    at.run()
    at.session_state["int_carrying_capacity"] = 1e9
    at.session_state["str_dorm_en_0"] = True   # enable dormancy (per-strain rate control)
    at.run()
    # Signal FUNCTIONS + Ks/threshold are now model-wide (topmost panel).
    at.session_state["int_dormancy_signal"] = "nutrient+density"
    at.session_state["int_resuscitation_signal"] = "nutrient+density"
    at.run()
    at.session_state["int_dormancy_monod_constant"] = 0.05   # model-wide dormancy Ks
    at.run()
    [b for b in at.button if "Run Simulation" in (b.label or "")][0].click().run()
    assert len(at.exception) == 0 and len(at.error) == 0
    assert at.session_state["simulation_config"].dormancy_monod_constant == 0.05

    # all four growth signals build without error
    for name, track in [("constant_growth", False), ("monod_growth", True),
                        ("logistic_growth", False), ("monod_logistic_growth", True)]:
        a = AppTest.from_file(APP, default_timeout=200)
        a.run()
        a.session_state["int_growth_function"] = name
        a.session_state["int_track_nutrients"] = track
        a.run()
        [b for b in a.button if "Run Simulation" in (b.label or "")][0].click().run()
        assert len(a.error) == 0, f"{name}: {[e.value for e in a.error]}"
        assert "simulation_result" in a.session_state and a.session_state["simulation_result"] is not None


def test_dormancy_signal_config_in_brg_and_strainset():
    """Dormancy signal FUNCTIONS are model-wide (single selector in the topmost panel) and
    map to the right engine dormancy function (+ density threshold) in BRG and StrainSet."""
    from streamlit.testing.v1 import AppTest

    at = AppTest.from_file(APP, default_timeout=200)
    at.run()
    at.session_state["widget_builder_mode"] = "Binary Genotypes (BRG)"
    at.run(); at.run()
    at.session_state["int_brg_dormancy_enabled"] = True
    # model-wide dormancy signal
    at.session_state["int_dormancy_signal"] = "nutrient+density"
    at.session_state["int_resuscitation_signal"] = "nutrient+density"
    at.run()
    [b for b in at.button if "Run Simulation" in (b.label or "")][0].click().run()
    cfg = at.session_state["simulation_config"]
    assert cfg.dormancy_function.__name__ == "nutrient_and_density_dormancy"
    assert cfg.dormancy_carrying_capacity == 1e8  # default density threshold
    assert len(at.error) == 0

    a = AppTest.from_file(APP, default_timeout=200)
    a.run()
    a.session_state["widget_builder_mode"] = "Custom Strains & Graph (StrainSet)"
    a.run(); a.run()
    a.session_state["ss_str_dorm_0"] = True
    a.session_state["int_dormancy_signal"] = "density"
    a.session_state["int_resuscitation_signal"] = "density"
    a.run()
    [b for b in a.button if "Run Simulation" in (b.label or "")][0].click().run()
    assert a.session_state["simulation_config"].dormancy_function.__name__ == "density_dependent_dormancy"
    assert len(a.error) == 0


def test_phage_two_compartment_pk_wired_into_config():
    """A phage with PK enabled and k12 > 0 builds a 2-compartment PhagePKConfig (peripheral
    tissue): the app's build sets k12/k21 on phage_pk_config and the solved result carries
    both the central (Pc) and peripheral (Pp) compartment states."""
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(APP, default_timeout=200)
    at.run()  # Direct mode by default
    # enable phage PK on phage 0, then set the peripheral-compartment rates via their widgets
    at.session_state["phg_pk_0"] = "Effect Compartment"
    at.run()
    at.session_state["phg_k12_0"] = 0.3
    at.session_state["phg_k21_0"] = 0.15
    at.run()
    [b for b in at.button if "Run Simulation" in (b.label or "")][0].click().run()
    assert len(at.error) == 0, at.error
    pk = at.session_state["simulation_config"].phage_pk_config
    assert pk is not None
    assert pk.k12 is not None and float(pk.k12[0]) == 0.3
    assert pk.k21 is not None and float(pk.k21[0]) == 0.15
    r = at.session_state["simulation_result"]
    assert r.get("Pc0") is not None and r.get("Pp0") is not None   # both compartments present


def test_reset_environment_clears_brg_and_strainset():
    """Reset Environment must clear BRG / StrainSet config, not just the Direct
    builder — the mode returns to Direct and the mode-specific keys are gone."""
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(APP, default_timeout=160)
    at.run()
    # configure BRG + a StrainSet mutation graph
    at.session_state["widget_builder_mode"] = "Binary Genotypes (BRG)"
    at.run(); at.run()
    at.session_state["int_brg_base_growth"] = 0.55
    at.session_state["int_brg_dormancy_enabled"] = True
    at.run()
    at.session_state["int_transitions"] = [{"from": "WT", "to": "R", "rate": 1e-7}]
    at.run()
    # Reset Environment (button lives on the AI Assistant page)
    at.session_state["current_page_radio"] = "AI Assistant"
    at.run()
    [b for b in at.button if "Reset Environment" in (b.label or "")][0].click().run()
    assert at.session_state["int_builder_mode"] == "Direct (ModelBuilder)"
    assert "int_brg_base_growth" not in at.session_state
    assert "int_brg_dormancy_enabled" not in at.session_state
    # the mutation graph is cleared (re-initialised to empty by session init)
    assert at.session_state["int_transitions"] == []
    assert "widget_builder_mode" not in at.session_state
    assert len(at.exception) == 0


def test_death_signal_function_selector():
    """The death-signal selector sets the model's death_function (default constant),
    and density death gets a carrying capacity so it doesn't crash under Monod growth."""
    from streamlit.testing.v1 import AppTest
    DEATH = {"constant": "constant_death",
             "nutrient (starvation)": "nutrient_dependent_death",
             "density (crowding)": "density_dependent_death",
             "nutrient + density": "nutrient_and_density_death"}
    for label, fn in DEATH.items():
        a = AppTest.from_file(APP, default_timeout=160)
        a.run()
        a.session_state["str_death_0"] = 0.05   # nonzero death rate via the widget
        a.run()
        [s for s in a.selectbox if "Death signal function" in (s.label or "")][0].set_value(label).run()
        [b for b in a.button if "Run Simulation" in (b.label or "")][0].click().run()
        cfg = a.session_state["simulation_config"]
        assert cfg.death_function.__name__ == fn, f"{label}: {cfg.death_function}"
        assert len(a.error) == 0


def test_lysis_signal_function_selector():
    """The lysis-signal selector wires frac_lysis (+ Ks_lysis) into the config across all
    three builder modes; constant leaves the engine default (None); and nutrient-coupled
    lysis falls back to constant when nutrients aren't tracked."""
    from streamlit.testing.v1 import AppTest

    # Direct, default = constant lysis → lysis_progression_function stays None.
    a = AppTest.from_file(APP, default_timeout=160); a.run()
    [b for b in a.button if "Run Simulation" in (b.label or "")][0].click().run()
    assert a.session_state["simulation_config"].lysis_progression_function is None

    # Direct, nutrient (Monod) lysis → frac_lysis + Ks_lysis flow through (Monod growth default).
    a = AppTest.from_file(APP, default_timeout=160); a.run()
    [s for s in a.selectbox if "Lysis signal function" in (s.label or "")][0].set_value("nutrient (Monod)").run()
    a.session_state["int_monod_constant_lysis"] = 0.15
    a.run()
    [b for b in a.button if "Run Simulation" in (b.label or "")][0].click().run()
    cfg = a.session_state["simulation_config"]
    assert cfg.lysis_progression_function.__name__ == "frac_lysis"
    assert cfg.monod_constant_lysis == 0.15
    assert len(a.error) == 0

    # Coercion: nutrient lysis + constant growth (S frozen) → falls back to constant lysis.
    a = AppTest.from_file(APP, default_timeout=160); a.run()
    a.session_state["int_growth_function"] = "constant_growth"
    a.session_state["int_track_nutrients"] = False
    a.session_state["int_lysis_function"] = "frac_lysis"
    a.run()
    [b for b in a.button if "Run Simulation" in (b.label or "")][0].click().run()
    assert a.session_state["simulation_config"].lysis_progression_function is None
    assert len(a.error) == 0

    # BRG + StrainSet: frac_lysis wired via to_config.
    for mode in ("Binary Genotypes (BRG)", "Custom Strains & Graph (StrainSet)"):
        b = AppTest.from_file(APP, default_timeout=200); b.run()
        b.session_state["widget_builder_mode"] = mode
        b.run(); b.run()
        b.session_state["int_lysis_function"] = "frac_lysis"
        b.session_state["int_monod_constant_lysis"] = 0.2
        b.run()
        [btn for btn in b.button if "Run Simulation" in (btn.label or "")][0].click().run()
        cfg = b.session_state["simulation_config"]
        assert cfg.lysis_progression_function.__name__ == "frac_lysis", mode
        assert cfg.monod_constant_lysis == 0.2, mode
        assert len(b.error) == 0, mode


def test_lysis_floor_reaches_config_all_modes():
    """The lysis floor φ_min (residual phage efficacy as nutrients deplete) reaches the
    config as ``lysis_floor`` under frac_lysis in all three builder modes; default is 0."""
    from streamlit.testing.v1 import AppTest
    for mode in ("Direct (ModelBuilder)", "Binary Genotypes (BRG)",
                 "Custom Strains & Graph (StrainSet)"):
        a = AppTest.from_file(APP, default_timeout=200); a.run()
        if not mode.startswith("Direct"):
            a.session_state["widget_builder_mode"] = mode; a.run(); a.run()
        a.session_state["int_lysis_function"] = "frac_lysis"
        a.session_state["int_lysis_floor"] = 0.3
        a.run()
        [b for b in a.button if "Run Simulation" in (b.label or "")][0].click().run()
        cfg = a.session_state["simulation_config"]
        assert cfg.lysis_progression_function.__name__ == "frac_lysis", mode
        assert cfg.lysis_floor == 0.3, (mode, cfg.lysis_floor)
        assert len(a.error) == 0, (mode, [e.value for e in a.error])

    # default (constant lysis) leaves the floor at 0
    a = AppTest.from_file(APP, default_timeout=160); a.run()
    [b for b in a.button if "Run Simulation" in (b.label or "")][0].click().run()
    assert a.session_state["simulation_config"].lysis_floor == 0.0


def test_new_growth_signals_build():
    """density_throttled_growth and gompertz_growth (new pbisim growth functions) build and
    reach the config; density-throttled forwards its Kd; gompertz reuses monod_constant/K."""
    from streamlit.testing.v1 import AppTest

    # Direct — each new signal builds and is set on the config.
    for fn, label in [("density_throttled_growth", "density-throttled (Monod × 1/(1+ΣB/Kd))"),
                      ("gompertz_growth", "Gompertz (nutrient)")]:
        a = AppTest.from_file(APP, default_timeout=160); a.run()
        [s for s in a.selectbox if "Growth signal function" in (s.label or "")][0].set_value(label).run()
        [b for b in a.button if "Run Simulation" in (b.label or "")][0].click().run()
        cfg = a.session_state["simulation_config"]
        assert cfg.growth_function.__name__ == fn, (label, cfg.growth_function)
        assert len(a.error) == 0, (label, [e.value for e in a.error])


def test_gompertz_uses_nutrient_scale_defaults_and_grows():
    """Gompertz k/S∞ come from dedicated NUTRIENT-scale keys (not the density-scale
    carrying_capacity), so with defaults the culture actually grows rather than sitting
    at a flat 0-growth curve (the exp overflow when S∞≈1e9)."""
    import numpy as np
    from streamlit.testing.v1 import AppTest
    a = AppTest.from_file(APP, default_timeout=160); a.run()
    # a density-scale carrying capacity must NOT leak into Gompertz S∞
    a.session_state["int_carrying_capacity"] = 1e9
    [s for s in a.selectbox if "Growth signal function" in (s.label or "")][0].set_value(
        "Gompertz (nutrient)").run()
    [b for b in a.button if "Run Simulation" in (b.label or "")][0].click().run()
    cfg = a.session_state["simulation_config"]
    assert cfg.monod_constant == 10.0 and cfg.carrying_capacity == 0.5  # k, S∞ from gompertz keys
    assert len(a.error) == 0
    # growth actually happens (not a flat line pinned at the inoculum)
    cfu = a.session_state["simulation_result"].sum_prefixes("B", "D")
    assert np.max(cfu) > 2.0 * cfu[0], (cfu[0], np.max(cfu))

    # density-throttled forwards density_growth_constant (Kd).
    a = AppTest.from_file(APP, default_timeout=160); a.run()
    [s for s in a.selectbox if "Growth signal function" in (s.label or "")][0].set_value(
        "density-throttled (Monod × 1/(1+ΣB/Kd))").run()
    a.session_state["int_density_growth_constant"] = 5e8
    a.run()
    [b for b in a.button if "Run Simulation" in (b.label or "")][0].click().run()
    assert a.session_state["simulation_config"].density_growth_constant == 5e8

    # BRG mode also accepts a new growth signal via to_config.
    b = AppTest.from_file(APP, default_timeout=200); b.run()
    b.session_state["widget_builder_mode"] = "Binary Genotypes (BRG)"; b.run(); b.run()
    b.session_state["int_growth_function"] = "density_throttled_growth"; b.run()
    [btn for btn in b.button if "Run Simulation" in (btn.label or "")][0].click().run()
    assert b.session_state["simulation_config"].growth_function.__name__ == "density_throttled_growth"
    assert len(b.error) == 0


def test_infected_nutrient_consumption_is_per_phage():
    """infected_nutrient_consumption is now a PER-PHAGE (phage×host) property: a per-phage
    input reaches the config as a (n_phages,) array in all three builder modes."""
    import numpy as np
    from streamlit.testing.v1 import AppTest
    for mode, key in [("Direct (ModelBuilder)", "phg_infnut_0"),
                      ("Binary Genotypes (BRG)", "brg_phg_infnut_0"),
                      ("Custom Strains & Graph (StrainSet)", "ss_phg_infnut_0")]:
        a = AppTest.from_file(APP, default_timeout=200); a.run()
        if not mode.startswith("Direct"):
            a.session_state["widget_builder_mode"] = mode; a.run(); a.run()
        _w = [n for n in a.number_input if n.key == key]
        assert _w, (mode, f"per-phage input {key} missing")
        _w[0].set_value(2.0).run()
        [b for b in a.button if "Run Simulation" in (b.label or "")][0].click().run()
        cfg = a.session_state["simulation_config"]
        inc = np.atleast_1d(np.asarray(cfg.infected_nutrient_consumption, dtype=float))
        assert inc.shape == (1,) and inc[0] == 2.0, (mode, inc)   # per-phage array on the config
        assert len(a.error) == 0, (mode, [e.value for e in a.error])

    # default is off — a zero per-phage array (or scalar 0), backward-compatible
    a = AppTest.from_file(APP, default_timeout=160); a.run()
    [b for b in a.button if "Run Simulation" in (b.label or "")][0].click().run()
    assert not np.any(a.session_state["simulation_config"].infected_nutrient_consumption)


def test_smooth_efficiency_monod_all_modes():
    """smooth_efficiency_monod builds in all three modes, forwarding the efficient K
    (monod_constant), the inefficient low-S K, and the θ/hill transition to the config."""
    from streamlit.testing.v1 import AppTest
    for mode in ("Direct (ModelBuilder)", "Binary Genotypes (BRG)",
                 "Custom Strains & Graph (StrainSet)"):
        a = AppTest.from_file(APP, default_timeout=200); a.run()
        if not mode.startswith("Direct"):
            a.session_state["widget_builder_mode"] = mode; a.run(); a.run()
        a.session_state["int_growth_function"] = "smooth_efficiency_monod"
        a.session_state["int_monod_constant"] = 0.2       # efficient (high-S) K
        a.session_state["int_monod_K_low"] = 5.0          # inefficient (low-S) K
        a.session_state["int_monod_efficiency_theta"] = 0.4
        a.session_state["int_monod_efficiency_hill"] = 3.0
        a.run()
        [b for b in a.button if "Run Simulation" in (b.label or "")][0].click().run()
        cfg = a.session_state["simulation_config"]
        assert cfg.growth_function.__name__ == "smooth_efficiency_monod", mode
        assert cfg.monod_constant == 0.2 and cfg.monod_K_low == 5.0, mode
        assert cfg.monod_efficiency_theta == 0.4 and cfg.monod_efficiency_hill == 3.0, mode
        assert len(a.error) == 0, (mode, [e.value for e in a.error])


def test_sequential_diauxic_growth_all_modes():
    """sequential_monod (diauxic growth) builds in all three modes, forwarding the three
    per-phase arrays to the config; the app validation helper mirrors the engine's rules."""
    from streamlit.testing.v1 import AppTest
    import numpy as np
    from pbisim_app.common import validate_sequential_growth

    for mode in ("Direct (ModelBuilder)", "Binary Genotypes (BRG)",
                 "Custom Strains & Graph (StrainSet)"):
        a = AppTest.from_file(APP, default_timeout=200); a.run()
        if not mode.startswith("Direct"):
            a.session_state["widget_builder_mode"] = mode; a.run(); a.run()
        a.session_state["int_growth_function"] = "sequential_monod"
        a.session_state["int_growth_n_phases"] = 2
        a.session_state["gp_monod_0"] = 0.5
        a.session_state["gp_rate_1"] = 0.4
        a.session_state["gp_monod_1"] = 0.1
        a.session_state["gp_thresh_0"] = 0.2
        a.run()
        [b for b in a.button if "Run Simulation" in (b.label or "")][0].click().run()
        cfg = a.session_state["simulation_config"]
        assert cfg.growth_function.__name__ == "sequential_monod", mode
        assert np.allclose(cfg.growth_phase_rate_factors, [1.0, 0.4]), mode
        assert np.allclose(cfg.growth_phase_monod, [0.5, 0.1]), mode
        assert np.allclose(cfg.growth_phase_thresholds, [0.2]), mode
        assert len(a.error) == 0, (mode, [e.value for e in a.error])

    # validation helper: strictly-decreasing thresholds in (0,1); matching lengths
    assert validate_sequential_growth([1.0, 0.5], [0.3, 0.3], [0.2]) is None
    assert validate_sequential_growth([1.0, 0.5], [0.3, 0.3], [1.5]) is not None   # θ ≥ 1
    assert validate_sequential_growth([1.0, 0.5, 0.2], [0.3, 0.3, 0.3],
                                      [0.2, 0.4]) is not None                       # not decreasing
    assert validate_sequential_growth([1.0, 0.5], [0.3], [0.2]) is not None         # length mismatch


def test_nutrient_dependent_od_config_all_modes():
    """Enabling nutrient-dependent OD sets od_absorptivity_function (+ asymptotes) on the
    config across all three builder modes, and get_od() stays finite."""
    from streamlit.testing.v1 import AppTest
    import numpy as np
    for mode in ("Direct (ModelBuilder)", "Binary Genotypes (BRG)",
                 "Custom Strains & Graph (StrainSet)"):
        a = AppTest.from_file(APP, default_timeout=200); a.run()
        if not mode.startswith("Direct"):
            a.session_state["widget_builder_mode"] = mode; a.run(); a.run()
        a.session_state["int_debris_enabled"] = True
        a.session_state["int_od_nutrient_enabled"] = True
        a.session_state["int_od_abs_exp"] = 1.0
        a.session_state["int_od_abs_stat"] = 0.4
        a.run()
        [b for b in a.button if "Run Simulation" in (b.label or "")][0].click().run()
        cfg = a.session_state["simulation_config"]
        assert cfg.od_absorptivity_function is not None, mode
        assert cfg.od_abs_stat == 0.4 and cfg.od_abs_exp == 1.0, mode
        assert len(a.error) == 0, (mode, [e.value for e in a.error])
        assert np.all(np.isfinite(a.session_state["simulation_result"].get_od())), mode


def test_prerun_inherit_debris_checkbox():
    """The 'Inherit bacterial debris' checkbox appears only when a pre-run AND the debris
    ODE are both on, and a run with it works."""
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(APP, default_timeout=160)
    at.run()

    # debris off -> no checkbox even with a pre-run
    at.session_state["int_t_prerun"] = 12.0
    at.session_state["int_debris_enabled"] = False
    at.run()
    assert not any("Inherit bacterial debris" in (c.label or "") for c in at.checkbox)

    # debris on + pre-run -> checkbox appears (default ticked) and a run succeeds
    at.session_state["int_debris_enabled"] = True
    at.session_state["int_prerun_inherit_debris"] = True
    at.run()
    box = [c for c in at.checkbox if "Inherit bacterial debris" in (c.label or "")]
    assert box and box[0].value is True
    [b for b in at.button if "Run Simulation" in (b.label or "")][0].click().run()
    assert len(at.exception) == 0, at.exception
    assert at.session_state["simulation_result"] is not None


def test_count_increment_sticks_no_revert():
    """Increasing the number of phages / bacteria persists — it must not blink back
    to the previous value on the follow-up rerun (the keyless value=len() self-fight)."""
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(APP, default_timeout=140)
    at.run()
    for label, state_key in (("Number of phages", "int_phages"),
                             ("Number of strains", "int_strains")):
        ni = [n for n in at.number_input if n.label == label]
        assert ni, f"{label} input missing"
        start = len(at.session_state[state_key])
        ni[0].set_value(start + 1).run()
        assert len(at.session_state[state_key]) == start + 1, f"{label} did not increase"
        at.run()  # the follow-up rerun (the 'blink') must not revert it
        assert len(at.session_state[state_key]) == start + 1, f"{label} reverted"
        assert [n for n in at.number_input if n.label == label][0].value == start + 1


def test_brg_dormant_adsorption_wired():
    """BRG exposes an adsorption-to-dormant widget and carries it into the built
    config (previously missing — dormant adsorption was silently 0 in BRG)."""
    import numpy as np
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(APP, default_timeout=160)
    at.run()
    [s for s in at.selectbox if s.label == "Bacterial Population Builder Mode"][0] \
        .set_value("Binary Genotypes (BRG)").run()
    ni = [n for n in at.number_input if n.label == "Adsorption to dormant WT (mL·h⁻¹)"]
    assert ni, "BRG dormant-adsorption widget missing"
    ni[0].set_value(5e-9).run()
    [b for b in at.button if "Run Simulation" in (b.label or "")][0].click().run()
    assert len(at.exception) == 0, at.exception
    cfg = at.session_state["simulation_config"]
    assert np.any(np.asarray(cfg.adsorption_rates_dormant) > 0), "dormant adsorption not applied"


def test_dormant_od_fraction_reaches_config_all_modes():
    """dormant_od_fraction (optical weight of dormant cells in OD) must reach the built
    config in EVERY builder mode. Previously only Direct passed it via with_od_debris;
    BRG/StrainSet omitted it from to_config's extra kwargs → engine default 1.0 silently
    ignored the user's setting."""
    from streamlit.testing.v1 import AppTest
    for mode in ("Direct (ModelBuilder)", "Binary Genotypes (BRG)",
                 "Custom Strains & Graph (StrainSet)"):
        at = AppTest.from_file(APP, default_timeout=200)
        at.run()
        at.session_state["widget_builder_mode"] = mode
        at.run(); at.run()
        at.session_state["int_debris_enabled"] = True
        at.session_state["int_dormant_od_fraction"] = 0.3
        at.run()
        [b for b in at.button if "Run Simulation" in (b.label or "")][0].click().run()
        cfg = at.session_state["simulation_config"]
        assert abs(float(cfg.dormant_od_fraction) - 0.3) < 1e-9, (mode, cfg.dormant_od_fraction)
        assert len(at.error) == 0, (mode, at.error)


def test_frac_lysis_survives_under_any_nutrient_growth():
    """Nutrient-coupled frac_lysis (and nutrient dormancy signals) must be kept under ANY
    nutrient-tracking growth signal — not silently coerced to constant unless growth is
    literally Monod. Regression: the check hardcoded (monod_growth, monod_logistic_growth)."""
    from streamlit.testing.v1 import AppTest
    for glabel, gfn in [("Gompertz (nutrient)", "gompertz_growth"),
                        ("smooth two-efficiency Monod", "smooth_efficiency_monod"),
                        ("density-throttled (Monod × 1/(1+ΣB/Kd))", "density_throttled_growth")]:
        a = AppTest.from_file(APP, default_timeout=160); a.run()
        [s for s in a.selectbox if "Growth signal function" in (s.label or "")][0].set_value(glabel).run()
        [s for s in a.selectbox if "Lysis signal function" in (s.label or "")][0].set_value("nutrient (Monod)").run()
        [b for b in a.button if "Run Simulation" in (b.label or "")][0].click().run()
        cfg = a.session_state["simulation_config"]
        assert cfg.lysis_progression_function.__name__ == "frac_lysis", (glabel, cfg.lysis_progression_function)
        assert len(a.error) == 0, (glabel, [e.value for e in a.error])

    # a NON-nutrient growth (constant) still coerces frac_lysis away (→ engine default
    # constant lysis, which the config represents as lysis_progression_function=None)
    a = AppTest.from_file(APP, default_timeout=160); a.run()
    [s for s in a.selectbox if "Growth signal function" in (s.label or "")][0].set_value("constant (unlimited)").run()
    [s for s in a.selectbox if "Lysis signal function" in (s.label or "")][0].set_value("nutrient (Monod)").run()
    [b for b in a.button if "Run Simulation" in (b.label or "")][0].click().run()
    _lfn = a.session_state["simulation_config"].lysis_progression_function
    assert _lfn is None or _lfn.__name__ != "frac_lysis"   # coerced away from frac_lysis


# ── Generator matrices: phenotypic transitions (bacteria) + phage mutation / transitions ──

def _sel(at, label):
    return [s for s in at.selectbox if label in (s.label or "")][0]


@pytest.mark.parametrize("mode", [
    "Direct (ModelBuilder)",
    "Binary Genotypes (BRG)",
    "Custom Strains & Graph (StrainSet)",
])
def test_generator_matrices_reach_config_and_repro_in_all_modes(mode):
    """The engine's transition_rates (bacterial phenotypic switching), mutation_rates_phage
    and transition_rates_phage generators used to be unreachable from every builder mode.
    Each mode now wires the shared edge lists (int_bact_transitions / int_phage_mutations /
    int_phage_transitions) into a mass-conserving M[dest, origin] matrix — via
    with_mutations (Direct), set_*_graph / set_phage_*_rates (StrainSet) or post-build
    cfg.<field> (BRG, whose to_config hard-codes zeros) — and the repro script mirrors it."""
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(APP, default_timeout=200)
    at.run()
    if mode != "Direct (ModelBuilder)":
        _sel(at, "Bacterial Population Builder Mode").set_value(mode)
        at.run()
    if mode == "Direct (ModelBuilder)":
        at.session_state["direct_n_strains"] = 2; at.session_state["direct_n_phages"] = 2
    elif mode == "Custom Strains & Graph (StrainSet)":
        at.session_state["ss_n_strains"] = 2; at.session_state["ss_n_phages"] = 2
    else:
        at.session_state["brg_n_phg_loci"] = 2
    at.run()
    strains, phages = at.session_state["int_strains"], at.session_state["int_phages"]
    bn = ["00", "01", "10", "11"] if mode == "Binary Genotypes (BRG)" else [s["name"] for s in strains]
    pn = [p["name"] for p in phages]
    at.session_state["int_bact_transitions"] = [{"from": bn[0], "to": bn[1], "rate": 0.02}]
    at.session_state["int_phage_mutations"] = [{"from": pn[0], "to": pn[1], "rate": 1e-5}]
    at.session_state["int_phage_transitions"] = [{"from": pn[1], "to": pn[0], "rate": 0.03}]
    at.run()
    assert len(at.exception) == 0, at.exception
    # the editors render (the expander auto-opens when edges exist)
    assert any("Add bacterial transition" in (b.label or "") for b in at.button)
    assert any("Add phage mutation" in (b.label or "") for b in at.button)
    [b for b in at.button if "Run Simulation" in (b.label or "")][0].click().run()
    assert len(at.exception) == 0, at.exception

    cfg = at.session_state["simulation_config"]
    n = len(bn)
    assert cfg.transition_rates.shape == (n, n)
    assert cfg.transition_rates[1, 0] == 0.02 and cfg.transition_rates[0, 0] == -0.02
    assert np.allclose(cfg.transition_rates.sum(axis=0), 0.0)      # mass-conserving
    assert cfg.mutation_rates_phage[1, 0] == 1e-5 and cfg.mutation_rates_phage[0, 0] == -1e-5
    assert cfg.transition_rates_phage[0, 1] == 0.03 and cfg.transition_rates_phage[1, 1] == -0.03
    # the reproduction script builds the same matrices
    ns = {}
    exec(compile(at.session_state["_last_repro_code"], "<repro>", "exec"), ns)
    for fld in ("transition_rates", "mutation_rates_phage", "transition_rates_phage"):
        assert np.allclose(getattr(ns["cfg"], fld), getattr(cfg, fld)), fld
    # the snapshot lists the switching section
    at.session_state["sim_show_cfg"] = True
    at.run()
    assert any("phenotypic switching" in (m.value or "") for m in at.markdown)


def test_phenotypic_transition_acts_without_growth():
    """transition_rates is a first-order switch on B itself (growth-independent), unlike
    mutation_rates which is coupled to division: with growth = 0 a WT→R transition still
    moves cells, a WT→R mutation does not."""
    from pbisim import ModelBuilder, PBIModel, solve_ode

    def _run(**mut):
        cfg = (ModelBuilder(n_bacteria=2, n_phages=1)
               .with_growth_rates([0.0, 0.0])
               .with_phage_params(adsorption_rates=[[0.0], [0.0]], burst_sizes=50,
                                  latent_periods=0.5, phage_decay_rates=0.1)
               .with_mutations(**mut).build())
        m = PBIModel(cfg, initial_B=np.array([1e7, 0.0]), initial_P=np.array([0.0]))
        r = solve_ode(m, t_end=10.0, dt=0.5)
        return float(r.get("B1")[-1]), float(r.sum_prefixes("B")[-1])

    b1_tr, tot_tr = _run(transition_rates=[[-0.1, 0.0], [0.1, 0.0]])
    b1_mu, _ = _run(mutation_rates=[[-0.1, 0.0], [0.1, 0.0]])
    assert b1_tr > 1e6                     # ~63% switched after 10 h at 0.1 h⁻¹
    assert np.isclose(tot_tr, 1e7, rtol=1e-3)   # conserved, just redistributed
    assert b1_mu < 1.0                     # no division → no mutation flux


def test_rate_graph_helpers_roundtrip():
    from pbisim_app.common import (rate_matrix_from_graph, edges_from_rate_matrix,
                                   graph_dict_from_edges, brg_genotype_labels)
    names = ["A", "B", "C"]
    edges = [{"from": "A", "to": "B", "rate": 0.1}, {"from": "A", "to": "C", "rate": 0.2},
             {"from": "C", "to": "A", "rate": 0.05},
             {"from": "A", "to": "A", "rate": 9.0},          # self-loop ignored
             {"from": "A", "to": "ZZZ", "rate": 1.0},        # unknown node ignored
             {"from": "B", "to": "A", "rate": 0.0}]          # zero rate ignored
    M = rate_matrix_from_graph(edges, names)
    assert M[1, 0] == 0.1 and M[2, 0] == 0.2 and M[0, 2] == 0.05
    assert np.isclose(M[0, 0], -0.3) and np.isclose(M[2, 2], -0.05) and M[1, 1] == 0
    assert np.allclose(M.sum(axis=0), 0.0)
    assert edges_from_rate_matrix(M, names) == [
        {"from": "A", "to": "B", "rate": 0.1}, {"from": "A", "to": "C", "rate": 0.2},
        {"from": "C", "to": "A", "rate": 0.05}]
    assert graph_dict_from_edges(edges, names) == {"A": {"B": 0.1, "C": 0.2}, "C": {"A": 0.05}}
    assert rate_matrix_from_graph([], names) is None
    assert brg_genotype_labels(2, 0) == ["00", "01", "10", "11"]
    assert brg_genotype_labels(1, 1) == ["phi0_abx0", "phi0_abx1", "phi1_abx0", "phi1_abx1"]
    assert brg_genotype_labels(0, 1) == ["abx0", "abx1"]
    # matches the engine's own labelling
    from pbisim.strains.genotypes import BinaryResistanceGenotypes, BacterialStrain, PhageStrain, Antibiotic
    brg = BinaryResistanceGenotypes.from_strains(
        [PhageStrain(name="p")], bacteria=BacterialStrain(base_growth_rate=1.2),
        antibiotics=[Antibiotic(name="a", emax_s=3.0, ec50_s=0.2, emax_r=0.1, ec50_r=2.0)])
    assert brg.strain_labels == brg_genotype_labels(1, 1)
