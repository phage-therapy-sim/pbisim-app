"""pbisim-fit NLS integration (Calibration §5c).

The nls_fit helpers must (a) offer only free params that fit the model's
dimensions, (b) build a pbisim-fit dataset from the app's aggregated calibration
data, and (c) recover ground-truth parameters from the committed tutorial CSV via
`refine_nls`. The fit is torch-free (lazy `import pbisim_fit`), so these run without
the SBI extras. Skipped if pbisim-fit is not installed in the env.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("pbisim_fit")

from pbisim import ModelBuilder
from pbisim.growth.signals import monod_growth

from pbisim_app import nls_fit as nls
from pbisim_app.fit_helper import normalize_fit_dataframe, aggregate_observations

CSV = "pbisim_app/examples/tutorial_synthetic_brg.csv"


def _agg_from_csv():
    df = pd.read_csv(CSV)
    long, _conds = normalize_fit_dataframe(df, "time", "value", "observable", ["arm"], None)
    return aggregate_observations(long, "raw", None)


def _monod_base(ratio=1e9):
    return (ModelBuilder(n_bacteria=1, n_phages=1)
            .with_growth_rates([1.2], bacteria_to_resource_ratio=[ratio])
            .with_growth_function(monod_growth)
            .with_nutrient(monod_constant=0.3)).build()


def test_available_free_params_respects_dimensions():
    cfg = _monod_base()
    labels = [l for l, *_ in nls.available_free_params(cfg)]
    assert "Growth rate — strain 0" in labels
    assert "Burst size — strain 0 × phage 0" in labels     # 1 strain × 1 phage
    assert "Growth rate — strain 1" not in labels          # only 1 strain
    assert "Monod Ks (global)" in labels


def test_build_dataset_dose_unit_pfu_vs_moi():
    """A dose value can be an absolute PFU/mL titre (unit='pfu') or an MOI multiplier
    (unit='moi'); build_dataset emits the right DoseRecord unit."""
    agg = _agg_from_csv()
    cond = {"1e+09 PFU": {"moi": 1e9}}
    ds_pfu = nls.build_dataset(agg, ["1e+09 PFU"], ["cfu"], cond, dose_unit="pfu")
    ds_moi = nls.build_dataset(agg, ["1e+09 PFU"], ["cfu"], cond, dose_unit="moi")
    d_pfu = ds_pfu.arms[0].dose_events[0]
    d_moi = ds_moi.arms[0].dose_events[0]
    assert d_pfu.unit == "pfu" and d_pfu.amount == 1e9   # absolute titre
    assert d_moi.unit == "moi" and d_moi.amount == 1e9   # multiplier of B0


def test_build_dataset_shapes_and_moi_dose():
    agg = _agg_from_csv()
    ds = nls.build_dataset(agg, ["control"], ["cfu", "od"],
                           {"control": {"moi": 2.0}}, od_to_cfu=2e8)
    assert len(ds.arms) == 1
    arm = ds.arms[0]
    assert arm.cfu is not None and arm.od is not None
    assert len(arm.cfu) == len(ds.time)
    # MOI condition became a phage dose
    assert any(d.target == "phage" and d.unit == "moi" and d.amount == 2.0
               for d in arm.dose_events)


def test_build_dataset_emits_bacteria_dose_additive_b0():
    """Additive-B0 model: a KNOWN inoculum (b0_is_dose) becomes a t=0 bacteria dose;
    a pre-run arm carries it as pretreatment_inoculum instead; first-observation mode
    (b0_is_dose=False) records neither, so pbisim-fit falls back to cfu[0]."""
    agg = _agg_from_csv()
    ds = nls.build_dataset(agg, ["control"], ["cfu"],
                           {"control": {"b0": 5e6, "b0_is_dose": True, "prerun": 0.0}})
    bd = [d for d in ds.arms[0].dose_events if d.target == "bacteria"]
    assert len(bd) == 1 and bd[0].amount == 5e6 and bd[0].unit == "cfu"
    assert ds.arms[0].pretreatment_inoculum is None

    ds_pr = nls.build_dataset(agg, ["control"], ["cfu"],
                              {"control": {"b0": 5e6, "b0_is_dose": True, "prerun": 24.0}})
    assert not [d for d in ds_pr.arms[0].dose_events if d.target == "bacteria"]
    assert ds_pr.arms[0].pretreatment_inoculum == 5e6

    ds_fo = nls.build_dataset(agg, ["control"], ["cfu"],
                              {"control": {"b0": 1e6, "b0_is_dose": False, "prerun": 0.0}})
    assert not [d for d in ds_fo.arms[0].dose_events if d.target == "bacteria"]
    assert ds_fo.arms[0].pretreatment_inoculum is None


def test_build_dataset_uses_imported_dose_records():
    """Imported NONMEM/Monolix dose rows (arm_doses) are emitted verbatim and override
    the manual per-arm fields for the targets they specify (inoculum + phage)."""
    agg = _agg_from_csv()
    arm_doses = {"control": [
        {"time": 0.0, "target": "bacteria", "amount": 3e6, "unit": "cfu"},
        {"time": 0.0, "target": "phage", "amount": 1e8, "unit": "pfu"}]}
    cond = {"control": {"moi": 5.0, "b0": 9e9, "b0_is_dose": True, "prerun": 0.0}}
    ds = nls.build_dataset(agg, ["control"], ["cfu"], cond, arm_doses=arm_doses)
    doses = ds.arms[0].dose_events
    bac = [d for d in doses if d.target == "bacteria"]
    ph = [d for d in doses if d.target == "phage"]
    assert len(bac) == 1 and bac[0].amount == 3e6           # data dose, not the manual 9e9
    assert len(ph) == 1 and ph[0].amount == 1e8 and ph[0].unit == "pfu"  # not the manual moi=5


def test_estimate_od_to_cfu_recovers_data_ratio():
    agg = _agg_from_csv()
    r = nls.estimate_od_to_cfu(agg, ["control"])
    assert 1.5e8 < r < 2.6e8   # tutorial truth od_to_cfu ~ 2e8


def test_build_param_spec_sharing_reduces_param_count():
    """Reparameterization: sharing ties several paths to ONE estimated value, so the
    spec has fewer parameters than paths."""
    from pbisim_fit.synthetic import reference_config
    cfg = reference_config()
    shared = [{"paths": ["growth_rates[0]", "growth_rates[1]"], "lo": 0.1, "hi": 3.0, "log": False}]
    _c, spec = nls.build_param_spec(cfg, [], shared_groups=shared)
    assert spec.n_params == 1                       # 2 paths → 1 shared theta


def test_nls_fit_sharing_recovers_shared_growth():
    """A shared growth theta across the two BRG strains recovers ~truth from control-arm
    CFU, with a single estimated parameter."""
    from pbisim_fit.synthetic import reference_config
    agg = _agg_from_csv()
    cfg = reference_config()
    # hi=2.0, not 3.0: the multi-start midpoint of the search box is the first start,
    # and the control CFU curve saturates by ~12 h, so growth is near-unidentifiable
    # above the saturating rate. With hi=3.0 the midpoint is 1.55 — in that flat region
    # AND above this test's own 1.5 ceiling — so the fit returns a barely-moved start.
    # Measured: hi=3.0 -> 1.550 (2 restarts) / 1.369 (3); hi=2.0 -> 1.193 at 2 restarts,
    # same wall-clock. Truth is 1.2. See "known product fragility" note in CLAUDE.md.
    shared = [{"paths": ["growth_rates[0]", "growth_rates[1]"], "lo": 0.1, "hi": 2.0, "log": False},
              {"paths": ["bacteria_to_resource_ratio[0]", "bacteria_to_resource_ratio[1]"],
               "lo": 1e6, "hi": 1e10, "log": True}]
    ds = nls.build_dataset(agg, ["control"], ["cfu"], {}, od_to_cfu=None)
    fp = nls.run_nls_fit(cfg, [], ds, ["cfu"], shared_groups=shared, n_restarts=3, max_nfev=300)
    assert fp.n_params == 2                         # 4 paths → 2 shared thetas
    m = fp.map()
    assert 0.9 < m["shr0"] < 1.5                    # shared growth ~ 1.2
    assert 5e7 < m["shr1"] < 2e8                    # shared ratio ~ 1e8


def test_build_param_spec_fitness_cost_derives_target():
    """A fitness-cost link binds target = source × (1 − cost); the fitted config's
    target path equals source × (1 − cost_MAP)."""
    from pbisim_fit.synthetic import reference_config
    cfg = reference_config()
    indiv = [("WT growth", "growth_rates[0]", 0.1, 3.0, False)]
    links = [{"source": "growth_rates[0]", "target": "growth_rates[1]",
              "kind": "fitness_cost", "lo": 0.0, "hi": 0.9}]
    _c, spec = nls.build_param_spec(cfg, indiv, fitness_links=links)
    assert spec.n_params == 2                        # WT growth + one cost
    ds = nls.build_dataset(_agg_from_csv(), ["control"], ["cfu"], {}, od_to_cfu=None)
    fp = nls.run_nls_fit(cfg, indiv, ds, ["cfu"], fitness_links=links, n_restarts=2, max_nfev=150)
    m = fp.map()
    fc = fp.to_config()
    g0, g1 = float(np.atleast_1d(fc.growth_rates)[0]), float(np.atleast_1d(fc.growth_rates)[1])
    assert abs(g1 - g0 * (1.0 - m["link0"])) < 1e-6  # derivation actually bound


def test_build_param_spec_scale_link():
    """A scale link binds target = source × factor."""
    from pbisim_fit.synthetic import reference_config
    cfg = reference_config()
    indiv = [("WT ads", "adsorption_rates[0,0]", 1e-10, 1e-7, True)]
    links = [{"source": "adsorption_rates[0,0]", "target": "adsorption_rates[1,0]",
              "kind": "scale", "lo": 0.0, "hi": 1.0}]
    _c, spec = nls.build_param_spec(cfg, indiv, fitness_links=links)
    assert spec.n_params == 2                        # WT adsorption + one scale factor


def test_estimate_b0_wires_free_initial_conditions():
    """The B0-source 'Estimate' modes add an estimated additive-B0 parameter via
    free_initial_conditions; a role-table fit_initial_cfu target is skipped so the two
    don't double-wire."""
    from pbisim_fit.synthetic import reference_config
    cfg = reference_config()
    ds = nls.build_dataset(_agg_from_csv(), ["control"], ["cfu"], {}, od_to_cfu=None)
    tg = [{"path": "growth_rates[0]", "free": True, "value": 1.0,
           "lo": 0.1, "hi": 3.0, "log": False, "prior_mu": None, "prior_sd": None}]
    _c0, spec0 = nls.build_param_spec_v2(cfg, tg)
    _c1, spec1 = nls.build_param_spec_v2(cfg, tg, dataset=ds, estimate_b0="shared")
    assert spec1.n_params == spec0.n_params + 1          # one shared B0 θ added
    # a fit_initial_cfu free target is dropped when estimating (no double-wire)
    tg2 = tg + [{"path": "fit_initial_cfu", "free": True, "value": 1e6,
                 "lo": 1e3, "hi": 1e11, "log": True, "prior_mu": None, "prior_sd": None}]
    _c2, spec2 = nls.build_param_spec_v2(cfg, tg2, dataset=ds, estimate_b0="shared")
    assert spec2.n_params == spec1.n_params


def test_estimate_b0_fit_runs_and_sets_initial_cfu():
    """A fit with estimate_b0='shared' completes and the fitted config carries an
    estimated fit_initial_cfu (the additive B0 offset)."""
    from pbisim_fit.synthetic import reference_config
    cfg = reference_config()
    ds = nls.build_dataset(_agg_from_csv(), ["control"], ["cfu"], {}, od_to_cfu=None)
    tg = [{"path": "growth_rates[0]", "free": True, "value": 1.0,
           "lo": 0.1, "hi": 3.0, "log": False, "prior_mu": None, "prior_sd": None}]
    fp = nls.run_nls_fit_v2(cfg, tg, [], [], ds, ["cfu"],
                            n_restarts=1, max_nfev=100, estimate_b0="shared")
    fc = fp.to_config()
    _ic = getattr(fc, "fit_initial_cfu", None)
    assert _ic is not None and float(np.atleast_1d(_ic).ravel()[0]) > 0


def test_arm_covariate_values_matches_arm_labels():
    """Per-arm covariate values are keyed by the same ' | '-joined arm label as the
    normaliser, so they line up with the dataset arms."""
    from pbisim_app.fit_helper import arm_covariate_values
    df = pd.DataFrame([{"PHAGE": "A", "MOI": 0.1, "temp": 30},
                       {"PHAGE": "A", "MOI": 0.1, "temp": 30},
                       {"PHAGE": "B", "MOI": 1.0, "temp": 37}])
    out = arm_covariate_values(df, ["PHAGE", "MOI"], ["MOI", "temp"])
    assert out == {"A | 0.1": {"MOI": 0.1, "temp": 30.0},
                   "B | 1.0": {"MOI": 1.0, "temp": 37.0}}


def test_build_dataset_attaches_per_arm_covariates():
    agg = _agg_from_csv()
    ds = nls.build_dataset(agg, ["control"], ["cfu"], {},
                           arm_covariates={"control": {"temp": 37.0}})
    assert ds.arms[0].covariates == {"temp": 37.0}


def test_build_param_spec_v2_wires_covariate_effect():
    """A covariate effect adds one estimated β to the spec and attaches the
    CovariateEffect structure to the config (consumed at solve time by pbisim-fit)."""
    from pbisim_fit.synthetic import reference_config
    cfg = reference_config()
    tg = [{"path": "growth_rates[0]", "free": True, "value": 1.0,
           "lo": 0.1, "hi": 3.0, "log": False, "prior_mu": None, "prior_sd": None}]
    _c0, spec0 = nls.build_param_spec_v2(cfg, tg)
    cov = [{"path": "growth_rates[0]", "covariate": "temp", "form": "power",
            "ref": 37.0, "beta_lo": -2.0, "beta_hi": 2.0, "beta_init": 0.0}]
    c1, spec1 = nls.build_param_spec_v2(cfg, tg, covariate_effects=cov)
    assert spec1.n_params == spec0.n_params + 1            # one β θ added
    effs = getattr(c1, "covariate_effects", None)
    assert effs and effs[0].path == "growth_rates[0]" and effs[0].covariate == "temp"
    assert effs[0].form == "power" and effs[0].ref == 37.0


def test_covariate_fit_runs_and_estimates_beta():
    """A short fit with a covariate link completes; the fitted config carries the
    covariate structure so the app can apply it per arm in the overlay."""
    from pbisim_fit.synthetic import reference_config
    cfg = reference_config()
    ds = nls.build_dataset(_agg_from_csv(), ["control"], ["cfu"],
                           {}, arm_covariates={"control": {"temp": 37.0}})
    tg = [{"path": "growth_rates[0]", "free": True, "value": 1.0,
           "lo": 0.1, "hi": 3.0, "log": False, "prior_mu": None, "prior_sd": None}]
    cov = [{"path": "growth_rates[0]", "covariate": "temp", "form": "power",
            "ref": 37.0, "beta_lo": -2.0, "beta_hi": 2.0, "beta_init": 0.0}]
    fp = nls.run_nls_fit_v2(cfg, tg, [], [], ds, ["cfu"], n_restarts=1, max_nfev=60,
                            covariate_effects=cov)
    fc = fp.to_config()
    assert getattr(fc, "covariate_effects", None)          # structure survives to_config


def test_available_targets_reflects_builder_mode():
    """The catalog is model-aware so there is one control per quantity: mutation/debris/
    monod are always estimable, B₀ is NEVER a table row (it's the per-arm B₀-source
    control), and the resistant-genotype parameterization depends on the builder mode —
    BRG exposes fitness cost / resistant fraction and HIDES the derived resistant growth
    rate; Direct/StrainSet exposes independent growth rates and HIDES those virtuals."""
    from pbisim_fit.synthetic import reference_config
    cfg = reference_config()   # 2-strain, debris on, mutation matrix present

    def paths(mode):
        return {p for (_l, p, *_r) in nls.available_targets(cfg, builder_mode=mode)}

    direct = paths("Direct (ModelBuilder)")
    assert "growth_rates[0]" in direct
    assert "mutation_rates[1,0]" in direct                    # mutation estimable
    assert {"debris_u", "debris_v", "debris_kdis"} <= direct  # debris estimable
    assert "monod_constant" in direct
    assert "fit_initial_pfu" in direct                        # phage titre still estimable
    assert "fit_initial_cfu" not in direct                    # B₀ is the per-arm B₀-source control
    # Direct/StrainSet: independent resistant growth, no BRG selection virtuals
    assert "growth_rates[1]" in direct
    assert "fitness_cost" not in direct and "init_resistant_fraction" not in direct
    ss = paths("Custom Strains & Graph (StrainSet)")
    assert "growth_rates[1]" in ss and "fitness_cost" not in ss

    # BRG: resistant growth = fitness cost; the raw resistant growth rate is hidden
    brg = paths("Binary Genotypes (BRG)")
    assert "fitness_cost" in brg and "init_resistant_fraction" in brg
    assert "growth_rates[0]" in brg and "growth_rates[1]" not in brg
    assert "fit_initial_cfu" not in brg


def test_available_targets_are_growth_signal_aware():
    """The fit-parameter catalog exposes the estimables specific to the ACTIVE growth
    signal (Gompertz S∞, smooth-efficiency low-S K / θ / Hill, diauxic phase params) plus
    the nutrient-environment params — and each builds a valid v2 fit spec."""
    from pbisim import (ModelBuilder, gompertz_growth, logistic_growth)

    def _paths(cfg):
        return {p for (_l, p, *_r) in nls.available_targets(cfg, initial_cfu=1e7)}

    base = lambda: ModelBuilder(n_bacteria=1, n_phages=1).with_growth_rates([1.0])

    smeff = base().with_smooth_efficiency_growth(3.0).with_nutrient(
        monod_constant=0.2, infected_nutrient_consumption=np.array([0.0])).build()
    p = _paths(smeff)
    assert {"monod_K_low", "monod_efficiency_theta", "monod_efficiency_hill"} <= p
    assert {"s_in", "s_out"} <= p                          # nutrient env now estimable
    assert "infected_nutrient_consumption[0]" in p         # per-phage (phage×host) estimable

    gom = base().with_growth_function(gompertz_growth).with_nutrient(monod_constant=10.0, carrying_capacity=0.5).build()
    assert "carrying_capacity" in _paths(gom)   # = Gompertz S∞

    seq = base().with_sequential_growth([1.0, 0.5], [0.1, 0.5], [0.3]).with_nutrient(monod_constant=0.1).build()
    ps = _paths(seq)
    assert "growth_phase_monod[0]" in ps and "growth_phase_thresholds[0]" in ps
    assert "growth_phase_rate_factors[1]" in ps and "growth_phase_rate_factors[0]" not in ps  # phase 1 pinned

    log = base().with_growth_function(logistic_growth).with_nutrient(carrying_capacity=1e9, monod_constant=0.3).build()
    assert "carrying_capacity" in _paths(log)

    # a nutrient-frozen (track_nutrients=False) model must NOT offer the nutrient-env estimables
    import dataclasses
    con = dataclasses.replace(base().build(), track_nutrients=False)
    assert "s_in" not in _paths(con)

    # each new estimable builds a valid v2 fit spec (incl. the per-phage indexed path)
    def _t(path, lo, hi, log_, val):
        return {"path": path, "lo": lo, "hi": hi, "log": log_, "value": val, "free": True}
    _, spec = nls.build_param_spec_v2(
        smeff, [_t("monod_K_low", 0.01, 50.0, True, 3.0),
                _t("infected_nutrient_consumption[0]", 0.0, 5.0, False, 0.0)], [], [])
    assert spec.n_params == 2
    _, spec2 = nls.build_param_spec_v2(seq, [_t("growth_phase_monod[1]", 0.01, 5.0, True, 0.5)], [], [])
    assert spec2.n_params == 1


def test_estimable_fitness_cost_and_initial_cfu():
    """fitness_cost and fit_initial_cfu are freeable via the v2 spec and recover
    sensible values (fit_initial_cfu ≈ the control-arm B₀)."""
    from pbisim_fit.synthetic import reference_config
    cfg = reference_config()
    tgts = [
        # hi=2.0 not 3.0 — see the saturation/midpoint note on the sharing test above.
        {"path": "growth_rates[0]", "free": True, "value": 1.0, "lo": 0.1, "hi": 2.0, "log": False},
        {"path": "fitness_cost", "free": True, "value": 0.0, "lo": 0.0, "hi": 0.9, "log": False},
        {"path": "fit_initial_cfu", "free": True, "value": 5e6, "lo": 1e3, "hi": 1e11, "log": True},
    ]
    ds = nls.build_dataset(_agg_from_csv(), ["control"], ["cfu"], {}, od_to_cfu=None)
    fp = nls.run_nls_fit_v2(cfg, tgts, [], [], ds, ["cfu"], n_restarts=2, max_nfev=200)
    m = fp.map()
    assert 0.9 < m["free0"] < 1.5                    # growth
    assert 1e6 < m["free2"] < 5e7                     # fit_initial_cfu ≈ control B₀ (5e6)


def test_virtual_params_not_set_when_fixed():
    """A fixed fitness_cost must NOT be pinned (setting fitness_cost=0 would wipe a
    BRG's baked-in resistant growth) — it acts only when freed."""
    from pbisim_fit.synthetic import reference_config
    import numpy as np
    cfg = reference_config()
    g_before = np.array(cfg.growth_rates, dtype=float).copy()
    tgts = [{"path": "growth_rates[0]", "free": True, "value": 1.0, "lo": 0.1, "hi": 3.0, "log": False},
            {"path": "fitness_cost", "free": False, "value": 0.0, "lo": 0.0, "hi": 1.0, "log": False}]
    out_cfg, _spec = nls.build_param_spec_v2(cfg, tgts, [], [])
    # fitness_cost must be unset on the built config (so _apply_fitness_cost is a no-op)
    assert getattr(out_cfg, "fitness_cost", None) is None


def test_validate_expr_guards_bad_input():
    ok, _ = nls.validate_expr("theta1*(1-theta2)", ["theta1", "theta2"])
    assert ok
    bad, _ = nls.validate_expr("theta1 + nope", ["theta1"])         # unknown name
    assert not bad
    bad2, _ = nls.validate_expr("__import__('os')", ["theta1"])     # no calls/attrs
    assert not bad2


def test_build_param_spec_v2_free_and_mapping():
    """v2 spec: a freed target plus a theta/mapping produce the right parameter count
    and bind the mapped path."""
    from pbisim_fit.synthetic import reference_config
    cfg = reference_config()
    targets = [{"path": p, "free": (p == "bacteria_to_resource_ratio[0]"), "value": v,
                "lo": lo, "hi": hi, "log": log}
               for (lab, p, v, lo, hi, log) in nls.available_targets(cfg)]
    # theta1 hi=2.0 not 3.0 — see the saturation/midpoint note on the sharing test above.
    thetas = [{"name": "theta1", "lo": 0.1, "hi": 2.0, "log": False, "initial": 1.0},
              {"name": "theta2", "lo": 0.0, "hi": 0.9, "log": False, "initial": 0.1}]
    mappings = [{"path": "growth_rates[0]", "expr": "theta1"},
                {"path": "growth_rates[1]", "expr": "theta1*(1-theta2)"}]
    _c, spec = nls.build_param_spec_v2(cfg, targets, thetas, mappings)
    assert spec.n_params == 3         # 1 free target + theta1 + theta2
    ds = nls.build_dataset(_agg_from_csv(), ["control"], ["cfu"], {}, od_to_cfu=None)
    fp = nls.run_nls_fit_v2(cfg, targets, thetas, mappings, ds, ["cfu"], n_restarts=2, max_nfev=200)
    m = fp.map()
    g = np.atleast_1d(np.asarray(fp.to_config().growth_rates))
    assert abs(g[1] - g[0] * (1.0 - m["theta2"])) < 1e-6   # mapping bound
    assert 0.9 < m["theta1"] < 1.5


def test_unbounded_and_one_sided_bounds():
    """Blank bounds → ±inf; run_nls_fit_v2 caps to a single start and still recovers."""
    import numpy as np
    from pbisim import ModelBuilder
    from pbisim.growth.signals import monod_growth
    cfg = (ModelBuilder(n_bacteria=1, n_phages=1)
           .with_growth_rates([1.2], bacteria_to_resource_ratio=[1e8])
           .with_growth_function(monod_growth).with_nutrient(monod_constant=0.3)).build()
    # growth fully unbounded (both inf); ratio one-sided (lower finite, upper inf)
    targets = [
        {"path": "growth_rates[0]", "free": True, "value": 1.0,
         "lo": float("-inf"), "hi": float("inf"), "log": False},
        {"path": "bacteria_to_resource_ratio[0]", "free": True, "value": 1e8,
         "lo": 1e6, "hi": float("inf"), "log": True},
    ]
    assert nls.has_unbounded(targets, [])
    ds = nls.build_dataset(_agg_from_csv(), ["control"], ["cfu"], {}, od_to_cfu=None)
    fp = nls.run_nls_fit_v2(cfg, targets, [], [], ds, ["cfu"], n_restarts=5, max_nfev=200)
    m = fp.map()
    assert 0.9 < m["free0"] < 1.5
    assert 5e7 < m["free1"] < 2e8


def test_prior_regularizes_map_estimate():
    """A Gaussian prior (prior_mu/prior_sd) turns the NLS into a MAP estimate: a tight
    prior pulls the estimate toward its mean, away from the data-only value (~1.2)."""
    from pbisim_fit.synthetic import reference_config
    cfg = reference_config()
    ds = nls.build_dataset(_agg_from_csv(), ["control"], ["cfu"], {}, od_to_cfu=None)

    def _fit(prior):
        # hi=2.5, not 3.0: the multi-start's first point is the box midpoint, and the
        # control CFU curve saturates by ~12 h, so growth is near-unidentifiable above
        # the saturating rate. hi=3.0 put that first point at 1.55 — inside the flat
        # region — so the fit returned a barely-moved start. hi=2.5 puts it at 1.30,
        # in the informative region, while keeping the 2.0 prior below comfortably
        # interior (hi=2.0 would sit that optimum exactly on the bound).
        # Measured: 3 restarts -> 1.1932 / 0.6011 / 1.9999, all three assertions, 54s.
        t = [{"path": "growth_rates[0]", "free": True, "value": 1.0,
              "lo": 0.1, "hi": 2.5, "log": False}]
        if prior:
            t[0]["prior_mu"], t[0]["prior_sd"] = prior
        return nls.run_nls_fit_v2(cfg, t, [], [], ds, ["cfu"], n_restarts=3, max_nfev=200).map()["free0"]

    assert 1.0 < _fit(None) < 1.4              # data-only ≈ 1.2
    assert abs(_fit((0.6, 0.03)) - 0.6) < 0.05  # tight prior at 0.6 dominates
    assert abs(_fit((2.0, 0.03)) - 2.0) < 0.05  # tight prior at 2.0 dominates


def test_fit_spec_dsl_parse_and_roundtrip():
    """The statement-based DSL parses into the table schemas, validates paths, and
    round-trips (serialize → parse gives the same thetas/maps/freed)."""
    from pbisim_fit.synthetic import reference_config
    cat = nls.available_targets(reference_config(), initial_cfu=5e6)
    spec = (
        "theta g bounds=0.1..3.0 prior=1.2,0.3\n"
        "theta cost bounds=0..0.9 init=0.05\n"
        "map growth_rates[0] = g\n"
        "map growth_rates[1] = g * (1 - cost)\n"
        "free bacteria_to_resource_ratio[0] bounds=1e6..1e10 log\n"
        "fix death_rate_B[0] = 0.0\n"
    )
    tdf, thdf, errs = nls.parse_fit_spec(spec, cat)
    assert not errs
    assert list(thdf["name"]) == ["g", "cost"]
    # mappings now live on the target rows as role='Derived' + an expression
    _derived = {r["path"]: r["expression"] for _, r in tdf.iterrows() if r["role"] == "Derived"}
    assert set(_derived) == {"growth_rates[0]", "growth_rates[1]"}
    assert "bacteria_to_resource_ratio[0]" in [r["path"] for _, r in tdf.iterrows() if r["role"] == "Free"]

    # unknown path is reported, not silently dropped
    _t, _th, e2 = nls.parse_fit_spec("free not_a_param bounds=0..1", cat)
    assert any("unknown parameter" in x for x in e2)

    # serialize the structures back and re-parse — thetas/maps preserved
    targets = [{"path": r["path"], "free": (r["role"] == "Free"),
                "value": float(r["value"]) if str(r["value"]).strip() else 0.0,
                "lo": (float(r["lower"]) if str(r["lower"]).strip() else 0.0),
                "hi": (float(r["upper"]) if str(r["upper"]).strip() else float("inf")),
                "log": bool(r["log"]), "prior_mu": None, "prior_sd": None}
               for _, r in tdf.iterrows()]
    thetas = [{"name": r["name"], "lo": float(r["lower"]), "hi": float(r["upper"]),
               "log": bool(r["log"]), "initial": (float(r["initial"]) if str(r["initial"]).strip() else None),
               "prior_mu": None, "prior_sd": None} for _, r in thdf.iterrows()]
    maps = [{"path": p, "expr": e} for p, e in _derived.items()]
    txt = nls.serialize_fit_spec(targets, thetas, maps, cat)
    _t2, _th2, e3 = nls.parse_fit_spec(txt, cat)
    assert not e3
    _derived2 = {r["path"] for _, r in _t2.iterrows() if r["role"] == "Derived"}
    assert list(_th2["name"]) == ["g", "cost"] and len(_derived2) == 2


def test_nls_fit_recovers_growth_from_tutorial_csv():
    """The headline check: CSV → refine_nls recovers the control-arm growth rate
    and yield (Monod base config, CFU + OD, few restarts)."""
    agg = _agg_from_csv()
    cfg = _monod_base()
    free = nls.available_free_params(cfg)
    by = {l: (l, p, lo, hi, g) for l, p, lo, hi, g in free}
    picks = [by["Growth rate — strain 0"], by["Bacteria/resource ratio — strain 0"]]
    od_link = nls.estimate_od_to_cfu(agg, ["control"])
    ds = nls.build_dataset(agg, ["control"], ["cfu", "od"], {}, od_to_cfu=od_link)
    # Uses the APP DEFAULT bounds from available_free_params (growth 0.1-2.0). 2 restarts
    # suffice now that the default upper bound puts the multi-start midpoint at 1.05,
    # inside the informative region — see the note on that default in nls_fit.py.
    fp = nls.run_nls_fit(cfg, picks, ds, ["cfu", "od"],
                         od_to_cfu=od_link, n_restarts=2, max_nfev=200)
    m = fp.map()
    assert 0.9 < m["growth_rates[0]"] < 1.5              # truth 1.2
    assert 5e7 < m["bacteria_to_resource_ratio[0]"] < 2e8   # truth 1e8


# ── Curve stripping (analytic warm-start) ────────────────────────────────────

def _od_strip_agg():
    """An `_agg`-shaped OD dataset: a no-phage control (MOI 0, logistic growth) plus
    two phage arms (growth → lysis dip → regrowth). od_to_cfu=1.0 in the test, so the
    OD values are read directly as CFU."""
    t = np.arange(0, 25.0)
    ctrl = np.minimum(1e7 * np.exp(0.38 * t), 1e9)
    up = np.minimum(1e7 * np.exp(0.38 * t[:8]), 1e9)
    lysis = np.concatenate([up, np.geomspace(up[-1], 1e3, 9), np.geomspace(1.5e3, 5e6, 8)])
    rows = []
    for arm, y in (("ctrl", ctrl), ("moi_0p1", lysis), ("moi_1", lysis)):
        for ti, v in zip(t, y):
            rows.append({"arm": arm, "observable": "od", "time": float(ti),
                         "value": float(v), "lo": float(v), "hi": float(v)})
    conds = {"ctrl": {"moi": 0.0}, "moi_0p1": {"moi": 0.1}, "moi_1": {"moi": 1.0}}
    return pd.DataFrame(rows), conds


def test_curve_strip_and_sanity_check():
    """strip_curves reads growth/capacity off the OD control and returns config-path
    initials; strip_sanity_check flags fitted values against the analytic bands."""
    agg, conds = _od_strip_agg()
    R = nls.strip_curves(agg, conds, b_fixed=50.0, od_to_cfu=1.0, monod_constant=0.3)

    assert np.isfinite(R.value("g_max")) and R.value("g_max") > 0
    assert np.isfinite(R.value("cap")) and R.value("cap") > 0
    inits = R.initials
    assert inits and "growth_rates[0]" in inits and "bacteria_to_resource_ratio[0]" in inits
    assert isinstance(R.report(), str) and R.report()

    # feed the stripped g_max back as the "fitted" value → ratio 1.0 → in band, ok
    rep = nls.strip_sanity_check(R, {"growth_rates[0]": R.value("g_max")})
    g_rows = [row for row in rep.rows if row[0] == "g_max"]
    assert g_rows and g_rows[0][5] is True          # (name, stripped, fitted, ratio, band, ok)
    assert rep.ok and not rep.flagged

    # a wildly-off fitted value is flagged
    rep2 = nls.strip_sanity_check(R, {"growth_rates[0]": R.value("g_max") * 100.0})
    assert "g_max" in rep2.flagged and not rep2.ok


def test_curve_strip_requires_od_control_and_moi_arm():
    """Clear ValueErrors when the data can't support stripping."""
    agg, conds = _od_strip_agg()
    # no OD observable
    with pytest.raises(ValueError):
        nls.strip_curves(agg.assign(observable="cfu"), conds,
                         b_fixed=50.0, od_to_cfu=1.0, monod_constant=0.3)
    # no MOI-0 control
    no_ctrl = {"ctrl": {"moi": 1.0}, "moi_0p1": {"moi": 0.1}, "moi_1": {"moi": 1.0}}
    with pytest.raises(ValueError):
        nls.strip_curves(agg, no_ctrl, b_fixed=50.0, od_to_cfu=1.0, monod_constant=0.3)
    # no phage arm
    all_ctrl = {"ctrl": {"moi": 0.0}, "moi_0p1": {"moi": 0.0}, "moi_1": {"moi": 0.0}}
    with pytest.raises(ValueError):
        nls.strip_curves(agg, all_ctrl, b_fixed=50.0, od_to_cfu=1.0, monod_constant=0.3)


def test_curve_strip_manual_control_arm():
    """A control with blank/NaN MOI isn't auto-detected, but an explicit control_arm works."""
    agg, _ = _od_strip_agg()
    nan_conds = {"ctrl": {"moi": float("nan")}, "moi_0p1": {"moi": 0.1}, "moi_1": {"moi": 1.0}}
    with pytest.raises(ValueError):
        nls.strip_curves(agg, nan_conds, b_fixed=50.0, od_to_cfu=1.0, monod_constant=0.3)
    R = nls.strip_curves(agg, nan_conds, b_fixed=50.0, od_to_cfu=1.0, monod_constant=0.3,
                         control_arm="ctrl")
    assert np.isfinite(R.value("g_max")) and "growth_rates[0]" in R.initials
    # a control_arm with no OD data is a clear error
    with pytest.raises(ValueError):
        nls.strip_curves(agg, nan_conds, b_fixed=50.0, od_to_cfu=1.0, monod_constant=0.3,
                         control_arm="nope")


def test_curve_strip_fit_f0_refinement():
    """fit_f0=True runs pbisim-fit's 1-D refine_f0 → the f0 estimate is tagged 'fit'."""
    agg, conds = _od_strip_agg()
    R = nls.strip_curves(agg, conds, b_fixed=50.0, od_to_cfu=1.0, monod_constant=0.3, fit_f0=True)
    assert R.estimates["f0"][1] == "fit"
    assert np.isfinite(R.value("f0"))


def test_set_config_path_and_skips_unknown():
    """set_config_path writes by path grammar and returns False (skips) for paths with no field."""
    cfg = (ModelBuilder(n_bacteria=2, n_phages=1).with_growth_rates([1.0, 0.9])
           .with_phage_params(burst_sizes=[[50], [50]], latent_periods=[[0.5], [0.5]],
                              adsorption_rates=[[1e-8], [1e-8]]).build())
    assert nls.set_config_path(cfg, "growth_rates[0]", 1.7) and nls._get_path(cfg, "growth_rates[0]") == 1.7
    assert nls.set_config_path(cfg, "adsorption_rates[1,0]", 2e-9)
    assert nls._get_path(cfg, "adsorption_rates[1,0]") == 2e-9
    assert nls.set_config_path(cfg, "monod_constant", 0.25) and nls._get_path(cfg, "monod_constant") == 0.25
    # unsettable field / out-of-range index → False (safely skipped, not an error)
    assert not nls.set_config_path(cfg, "init_resistant_fraction", 1e-3)
    assert not nls.set_config_path(cfg, "growth_rates[9]", 1.0)
    assert not nls.set_config_path(cfg, "adsorption_rates[0,9]", 1.0)


def test_strip_growth_signal_mapping():
    """App growth-function names map to the three signals curve stripping supports."""
    assert nls.strip_growth_signal("monod_growth") == "monod"
    assert nls.strip_growth_signal("logistic_growth") == "logistic"
    assert nls.strip_growth_signal("constant_growth") == "constant"
    for fn in ("gompertz_growth", "density_throttled_growth", "smooth_efficiency_monod",
               "monod_logistic_growth", "sequential_monod", "something_new"):
        assert nls.strip_growth_signal(fn) == "monod"


def test_curve_strip_growth_signal_and_debris_passthrough():
    """growth_signal picks the forward model; fit_f0 + a debris dict runs the direct refine_f0."""
    agg, conds = _od_strip_agg()
    R = nls.strip_curves(agg, conds, b_fixed=50.0, od_to_cfu=1.0, monod_constant=0.3,
                         growth_signal="logistic")
    assert np.isfinite(R.value("g_max"))
    R2 = nls.strip_curves(agg, conds, b_fixed=50.0, od_to_cfu=1.0, monod_constant=0.3,
                          fit_f0=True, debris={"v": 0.25, "kdis": 0.02, "u": 0.1}, f0_latent=0.8)
    assert R2.estimates["f0"][1] == "fit" and np.isfinite(R2.value("f0"))


def test_broadcast_growth_maps_g_to_all_species():
    """_broadcast_growth spreads g_max to every strain's growth_rates[i]; other params untouched."""
    from pbisim_app.views.calibration import _broadcast_growth
    out = _broadcast_growth({"growth_rates[0]": 0.5, "bacteria_to_resource_ratio[1]": 3e8}, 3)
    assert out["growth_rates[0]"] == out["growth_rates[1]"] == out["growth_rates[2]"] == 0.5
    assert out["bacteria_to_resource_ratio[1]"] == 3e8
    # no growth estimate → unchanged
    assert _broadcast_growth({"adsorption_rates[0,0]": 1e-8}, 2) == {"adsorption_rates[0,0]": 1e-8}


def test_split_initial_B_by_f0():
    """f0 maps onto the WT/resistant initial_B split (strain 1 = resistant), total preserved."""
    from pbisim_app.views.calibration import _split_initial_B_by_f0
    strains = [{"initial_B": 6e6}, {"initial_B": 2e6}]        # total 8e6
    assert _split_initial_B_by_f0(strains, 0.25)
    assert abs(strains[0]["initial_B"] - 6e6) < 1e-3          # (1-0.25)*8e6
    assert abs(strains[1]["initial_B"] - 2e6) < 1e-3          # 0.25*8e6
    # guards: <2 strains, zero total, None
    assert not _split_initial_B_by_f0([{"initial_B": 1e7}], 0.1)
    assert not _split_initial_B_by_f0([{"initial_B": 0}, {"initial_B": 0}], 0.1)
    assert not _split_initial_B_by_f0(strains, None)
    # clamp to [0,1]
    _s2 = [{"initial_B": 1e6}, {"initial_B": 1e6}]
    assert _split_initial_B_by_f0(_s2, 1.5)
    assert _s2[0]["initial_B"] == 0.0 and _s2[1]["initial_B"] == 2e6


def test_strip_curves_picker_guards():
    """require_growth_peak + min_virulence thread through to propose_initials (opt-in, no crash)."""
    agg, conds = _od_strip_agg()
    R = nls.strip_curves(agg, conds, b_fixed=50.0, od_to_cfu=1.0, monod_constant=0.3,
                         require_growth_peak=True, min_virulence=0.05)
    assert np.isfinite(R.value("g_max"))          # control-derived, unaffected by the guards


def _amort_dataset():
    import pandas as pd
    t = np.arange(0, 25.0, 1.0)
    ctrl = np.minimum(0.02 * np.exp(0.32 * t), 1.2)
    up = np.minimum(0.02 * np.exp(0.32 * t[:8]), 1.2)
    lysis = np.concatenate([up, np.geomspace(up[-1], 0.03, 9), np.geomspace(0.04, 0.6, 8)])
    rows = []
    for moi, y in ((0.0, ctrl), (0.1, lysis), (1.0, lysis)):
        for ti, v in zip(t, y):
            rows.append({"arm": f"m{moi}", "observable": "od", "time": float(ti),
                         "value": float(v), "lo": v, "hi": v})
    agg = pd.DataFrame(rows)
    ac = {f"m{m}": {"b0": 4e6, "moi": m} for m in (0.0, 0.1, 1.0)}
    return nls.build_dataset(agg, list(ac), ["od"], ac, od_to_cfu=2.4e9, dose_unit="moi")


def test_amortized_available_and_nets():
    """list_amortized_nets works WITHOUT importing onnxruntime (lazy) — the page can list nets."""
    nets = nls.list_amortized_nets()
    assert isinstance(nets, list)
    if nets:
        assert any("MXP1001" in n for n in nets)


def test_amortized_fit_torch_free():
    """The full amortized ONNX path is torch-free and returns map/ci/initials/config."""
    import sys
    pytest.importorskip("onnxruntime")            # skips on the dev box unless conda lib is on LD path
    if not nls.list_amortized_nets():
        pytest.skip("no amortized nets available")
    ds = _amort_dataset()
    res = nls.amortized_fit(ds, nls.list_amortized_nets()[0], burst=50.0, latent=0.5)
    assert set(res["map"]) >= {"g", "cap", "cap_r", "k", "f0"}
    assert res["initials"]["growth_rates[0]"] == pytest.approx(res["map"]["g"])
    assert res["initials"]["latent_periods[0,0]"] == 0.5      # frozen context in initials
    assert res["config"] is not None and float(res["config"].latent_periods[0, 0]) == 0.5
    assert "torch" not in sys.modules


# ── Generator matrices: fit targets + mass-conserving diagonal ────────────────────────

def _generator_cfg():
    return (ModelBuilder(n_bacteria=2, n_phages=2)
            .with_growth_rates([1.0, 0.9])
            .with_phage_params(adsorption_rates=[[1e-8, 1e-8], [0.0, 0.0]], burst_sizes=50,
                               latent_periods=0.5, phage_decay_rates=[0.1, 0.1])
            .with_mutations(mutation_rates=[[-1e-7, 0], [1e-7, 0]],
                            transition_rates=[[-0.02, 0], [0.02, 0]],
                            mutation_rates_phage=[[-1e-6, 0], [1e-6, 0]],
                            transition_rates_phage=[[0, 0.03], [0, -0.03]])
            .build())


def test_available_targets_include_configured_generator_edges():
    """Phenotypic-transition + phage generator entries are offered ONLY for edges the
    model configures (nonzero), while mutation_rates keeps offering every off-diagonal."""
    paths = {p for (_l, p, *_r) in nls.available_targets(_generator_cfg())}
    assert {"mutation_rates[1,0]", "mutation_rates[0,1]"} <= paths
    assert "transition_rates[1,0]" in paths and "transition_rates[0,1]" not in paths
    assert "mutation_rates_phage[1,0]" in paths and "mutation_rates_phage[0,1]" not in paths
    assert "transition_rates_phage[0,1]" in paths and "transition_rates_phage[1,0]" not in paths
    assert not any(p.endswith("[0,0]") or p.endswith("[1,1]") for p in paths
                   if p.split("[")[0] in nls._GENERATOR_FIELDS)     # diagonals never freeable


def test_freed_generator_entry_keeps_column_mass_conserving():
    """Freeing / re-fixing an off-diagonal generator rate must rebalance that column's
    diagonal at solve time — otherwise the fit would create or destroy cells instead of
    moving them (this leaked mass for the pre-existing mutation_rates targets too)."""
    cfg = _generator_cfg()
    targets = [
        {"path": "mutation_rates[1,0]", "free": True, "value": 1e-7, "lo": 1e-9, "hi": 1e-4, "log": True},
        {"path": "transition_rates[1,0]", "free": False, "value": 0.05, "lo": 1e-4, "hi": 10, "log": True},
        {"path": "transition_rates_phage[0,1]", "free": True, "value": 0.03, "lo": 1e-4, "hi": 10, "log": True},
        {"path": "mutation_rates[0,1]", "free": False, "value": 0.0, "lo": 1e-9, "hi": 1e-4, "log": True},
    ]
    base, spec = nls.build_param_spec_v2(cfg, targets)
    assert spec.paths == ["free0", "free2"]
    x = spec.initial_vector()
    c2 = spec.apply(x, base)
    vals = spec.as_dict(x)
    assert np.isclose(c2.mutation_rates[1, 0], vals["free0"])
    assert np.isclose(c2.mutation_rates[0, 0], -vals["free0"])
    assert np.isclose(c2.transition_rates_phage[0, 1], vals["free2"])
    assert np.isclose(c2.transition_rates_phage[1, 1], -vals["free2"])
    assert c2.transition_rates[1, 0] == 0.05 and c2.transition_rates[0, 0] == -0.05   # re-fixed
    for fld in nls._GENERATOR_FIELDS:
        assert np.allclose(getattr(c2, fld).sum(axis=0), 0.0), fld
    # a mapped (Derived) generator entry is rebalanced too
    thetas = [{"name": "k", "lo": 1e-4, "hi": 1.0, "log": True, "initial": 0.1}]
    base, spec = nls.build_param_spec_v2(
        cfg, [{"path": "transition_rates[1,0]", "free": False, "value": 0.02,
               "lo": 1e-4, "hi": 10, "log": True}],
        thetas=thetas, mappings=[{"path": "transition_rates[1,0]", "expr": "2*k"}])
    c3 = spec.apply(spec.initial_vector(), base)
    assert np.isclose(c3.transition_rates[1, 0], 0.2) and np.isclose(c3.transition_rates[0, 0], -0.2)
