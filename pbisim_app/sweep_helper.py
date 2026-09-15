"""
sweep_helper.py — Utility functions for mapping, applying, and executing parameter sweeps
and dose-response simulations in the pbisim Streamlit app.
"""

from __future__ import annotations

import dataclasses
import numpy as np
from pbisim import PBIModel, solve_ode, stationary_phase_ic

def get_sweep_parameters(config, strains=None, phages=None, antibiotics=None) -> dict:
    """
    Dynamically constructs a dictionary of sweepable parameters based on the current model configuration.
    Maps user-friendly labels to target fields, types, and indices.
    """
    params = {}

    # 1. Global scalars (always available)
    params["Monod Constant (Ks)"] = {
        "type": "scalar",
        "field": "monod_constant",
    }
    params["Carrying Capacity (K)"] = {
        "type": "scalar",
        "field": "carrying_capacity",
    }
    params["Recycle Fraction"] = {
        "type": "scalar",
        "field": "recycle_fraction",
    }
    params["Medium Inflow (s_in)"] = {
        "type": "scalar",
        "field": "s_in",
    }
    params["Medium Washout (s_out)"] = {
        "type": "scalar",
        "field": "s_out",
    }
    params["Immune Decay Rate"] = {
        "type": "scalar",
        "field": "imm_decay_rate",
    }
    params["Immune Kill 50 (K_kill)"] = {
        "type": "scalar",
        "field": "imm_kill50",
    }
    params["Immune Capacity (Imm_max)"] = {
        "type": "scalar",
        "field": "imm_max",
    }
    params["OD-to-CFU Conversion Factor"] = {
        "type": "scalar",
        "field": "od_to_cfu_conversion_factor",
    }
    params["Debris Dissolution Rate (k_dis)"] = {
        "type": "scalar",
        "field": "debris_kdis",
    }
    params["Immune Stim 50 (K_stim)"] = {
        "type": "scalar",
        "field": "imm_stim50",
    }
    # Pre-run duration is an app-level step (stationary_phase_ic before treatment),
    # NOT a ModelConfig field — apply_sweep_parameter stashes the swept value and the
    # sweep loop runs the pre-run per point (still honouring the debris-inheritance
    # choice, since each duration accumulates a different amount of debris).
    params["Pre-run duration (hours)"] = {
        "type": "prerun",
        "default": 24.0,
    }
    # Dormancy sensing thresholds (None → inherit growth Ks / K; give a numeric
    # `default` so the 1D UI has a sensible nominal even when the field is None).
    params["Dormancy Nutrient Half-saturation (Ks_dorm)"] = {
        "type": "scalar", "field": "dormancy_monod_constant", "default": 0.3,
    }
    params["Dormancy Density Threshold (K_dorm)"] = {
        "type": "scalar", "field": "dormancy_carrying_capacity", "default": 1e8,
    }
    params["Lysis Monod Constant (Ks_lysis)"] = {
        "type": "scalar", "field": "monod_constant_lysis", "default": 0.3,
    }
    params["Lysis Floor (φ_min)"] = {
        "type": "scalar", "field": "lysis_floor", "default": 0.0,
    }
    if getattr(config, "debris_u", None) is not None:
        params["Debris Yield · deaths (u)"] = {"type": "scalar", "field": "debris_u"}
    if getattr(config, "debris_v", None) is not None:
        params["Debris Yield · lysis (v)"] = {"type": "scalar", "field": "debris_v"}

    # Dimensions (Q & L)
    params["Dormancy Depth Layers (Q)"] = {
        "type": "dimension",
        "field": "n_depth",
    }
    params["Phage Latent Stages (L)"] = {
        "type": "dimension",
        "field": "n_latent",
    }

    # 2. Strain-specific parameters
    for i in range(config.n_bacteria):
        strain_name = f"Strain {i}"
        if strains and i < len(strains):
            strain_name = f"Strain {i} ({strains[i]['name']})"
        elif hasattr(config, "strain_labels") and config.strain_labels and i < len(config.strain_labels):
            strain_name = f"Strain {i} ({config.strain_labels[i]})"

        params[f"Growth Rate - {strain_name}"] = {
            "type": "array1d",
            "field": "growth_rates",
            "index": i,
        }
        params[f"Bacteria-Resource Ratio - {strain_name}"] = {
            "type": "array1d",
            "field": "bacteria_to_resource_ratio",
            "index": i,
        }
        params[f"Dormancy Rate (sleep) - {strain_name}"] = {
            "type": "array1d",
            "field": "dormancy_rate",
            "index": i,
        }
        params[f"Resuscitation Rate (wake) - {strain_name}"] = {
            "type": "array1d",
            "field": "resuscitation_rate",
            "index": i,
        }
        params[f"Dormancy Diffusion Rate - {strain_name}"] = {
            "type": "array1d",
            "field": "dormancy_diffusion_rate",
            "index": i,
        }
        params[f"Natural Death Rate (dB) - {strain_name}"] = {
            "type": "array1d_or_none",
            "field": "death_rate_B",
            "index": i,
        }
        params[f"Dormant Death Rate (dD) - {strain_name}"] = {
            "type": "array1d_or_none",
            "field": "death_rate_D",
            "index": i,
        }
        params[f"Immune Kill Rate (Dormant) - {strain_name}"] = {
            "type": "array1d_or_none",
            "field": "imm_kill_rate_D",
            "index": i,
        }
        params[f"Immune Stim Rate - {strain_name}"] = {
            "type": "array1d",
            "field": "imm_stim_rate",
            "index": i,
        }
        params[f"Immune Kill Rate (Active) - {strain_name}"] = {
            "type": "array1d",
            "field": "imm_kill_rate",
            "index": i,
        }

    # 2b. Strain→strain generator rates: mutation (per division) and phenotypic
    # transitions (h⁻¹, growth-independent). Mass-conserving: sweeping origin→dest
    # rebalances the origin column's diagonal (see the "mutation" handler in
    # apply_sweep_parameter, which serves every generator matrix via meta["field"]).
    if config.n_bacteria > 1:
        def _sname(k):
            if strains and k < len(strains):
                return strains[k]["name"]
            if getattr(config, "strain_labels", None) and k < len(config.strain_labels):
                return config.strain_labels[k]
            return f"Strain {k}"
        for _fld, _lab, _dflt in (("mutation_rates", "Mutation Rate", 1e-7),
                                  ("transition_rates", "Phenotypic Transition Rate", 0.01)):
            if getattr(config, _fld, None) is None:
                continue
            for _o in range(config.n_bacteria):
                for _d in range(config.n_bacteria):
                    if _o == _d:
                        continue
                    params[f"{_lab} - {_sname(_o)} → {_sname(_d)}"] = {
                        "type": "mutation", "field": _fld, "origin": _o, "dest": _d,
                        "default": _dflt,
                    }

    # 3. Phage-specific parameters
    def _pname(j):
        return f"Phage {j} ({phages[j]['name']})" if phages and j < len(phages) else f"Phage {j}"
    # 3a. Phage→phage generator rates: mutation (per burst) and phenotypic transitions
    # (h⁻¹ on free phage) — same mass-conserving "mutation" handler, keyed by field.
    if config.n_phages > 1:
        for _fld, _lab, _dflt in (("mutation_rates_phage", "Phage Mutation Rate", 1e-6),
                                  ("transition_rates_phage", "Phage Transition Rate", 0.01)):
            if getattr(config, _fld, None) is None:
                continue
            for _o in range(config.n_phages):
                for _d in range(config.n_phages):
                    if _o == _d:
                        continue
                    params[f"{_lab} - {_pname(_o)} → {_pname(_d)}"] = {
                        "type": "mutation", "field": _fld, "origin": _o, "dest": _d,
                        "default": _dflt,
                    }
    for j in range(config.n_phages):
        phage_name = _pname(j)

        params[f"Phage Decay Rate - {phage_name}"] = {
            "type": "array1d",
            "field": "phage_decay_rates",
            "index": j,
        }
        if getattr(config, "phage_decay_Km", None) is not None:
            params[f"Phage Decay Km (MM saturation) - {phage_name}"] = {
                "type": "array1d", "field": "phage_decay_Km", "index": j,
            }
        # infected_nutrient_consumption is per-phage (an array on the config); only
        # sweepable per-phage when the model actually tracks nutrients (else it's a scalar).
        _inc = getattr(config, "infected_nutrient_consumption", 0.0)
        if np.ndim(_inc) > 0 and j < np.size(_inc):
            params[f"Infected-cell nutrient consumption - {phage_name}"] = {
                "type": "array1d", "field": "infected_nutrient_consumption", "index": j,
            }

        # 2D arrays: adsorption, burst, latent (bacteria x phage)
        for i in range(config.n_bacteria):
            strain_name = f"Strain {i}"
            if strains and i < len(strains):
                strain_name = strains[i]['name']
            elif hasattr(config, "strain_labels") and config.strain_labels and i < len(config.strain_labels):
                strain_name = config.strain_labels[i]

            params[f"Adsorption - {phage_name} on {strain_name}"] = {
                "type": "array2d",
                "field": "adsorption_rates",
                "index_row": i,
                "index_col": j,
            }
            params[f"Dormant Adsorption - {phage_name} on {strain_name}"] = {
                "type": "array2d",
                "field": "adsorption_rates_dormant",
                "index_row": i,
                "index_col": j,
            }
            params[f"Burst Size - {phage_name} on {strain_name}"] = {
                "type": "array2d",
                "field": "burst_sizes",
                "index_row": i,
                "index_col": j,
            }
            params[f"Latent Period - {phage_name} on {strain_name}"] = {
                "type": "array2d",
                "field": "latent_periods",
                "index_row": i,
                "index_col": j,
            }
            params[f"Dormant Latent Period - {phage_name} on {strain_name}"] = {
                "type": "array2d",
                "field": "latent_periods_dormant",
                "index_row": i,
                "index_col": j,
            }
            params[f"Dormant Adsorption Attenuation - {phage_name} on {strain_name}"] = {
                "type": "array2d",
                "field": "attenuation_rate",
                "index_row": i,
                "index_col": j,
            }
            # Pseudolysogeny (only present when enabled): I→H hibernation and H→I
            # lytic resumption, per strain × phage.
            if getattr(config, "hibernation_rate", None) is not None:
                params[f"Hibernation Rate (I→H) - {phage_name} on {strain_name}"] = {
                    "type": "array2d", "field": "hibernation_rate", "index_row": i, "index_col": j,
                }
            if getattr(config, "lytic_resumption_rate", None) is not None:
                params[f"Lytic Resumption Rate (H→I) - {phage_name} on {strain_name}"] = {
                    "type": "array2d", "field": "lytic_resumption_rate", "index_row": i, "index_col": j,
                }

    # 4. Antibiotics-specific parameters
    if config.n_antibiotics > 0 and config.pk_config is not None:
        for j in range(config.n_antibiotics):
            abx_name = f"Antibiotic {j}"
            if antibiotics and j < len(antibiotics):
                abx_name = f"Antibiotic {j} ({antibiotics[j]['name']})"

            params[f"Clearance (k_elim) - {abx_name}"] = {
                "type": "pk_array1d",
                "field": "k_elim",
                "index": j,
            }
            params[f"Volume (Vc) - {abx_name}"] = {
                "type": "pk_array1d",
                "field": "Vc",
                "index": j,
            }

            if config.pd_config is not None:
                for i in range(config.n_bacteria):
                    strain_name = f"Strain {i}"
                    if strains and i < len(strains):
                        strain_name = strains[i]['name']
                    elif hasattr(config, "strain_labels") and config.strain_labels and i < len(config.strain_labels):
                        strain_name = config.strain_labels[i]

                    params[f"Emax - {abx_name} on {strain_name}"] = {
                        "type": "pd_array2d",
                        "field": "abx_emax",
                        "index_row": i,
                        "index_col": j,
                    }
                    params[f"EC50 - {abx_name} on {strain_name}"] = {
                        "type": "pd_array2d",
                        "field": "abx_ec50",
                        "index_row": i,
                        "index_col": j,
                    }
                    params[f"Hill (H) - {abx_name} on {strain_name}"] = {
                        "type": "pd_array2d",
                        "field": "abx_hill",
                        "index_row": i,
                        "index_col": j,
                    }

    # 5. Initial conditions
    for i in range(config.n_bacteria):
        strain_name = f"Strain {i}"
        if strains and i < len(strains):
            strain_name = f"Strain {i} ({strains[i]['name']})"
        elif hasattr(config, "strain_labels") and config.strain_labels and i < len(config.strain_labels):
            strain_name = f"Strain {i} ({config.strain_labels[i]})"

        params[f"Initial Density (B0) - {strain_name}"] = {
            "type": "initial_B",
            "index": i,
        }

    for j in range(config.n_phages):
        phage_name = f"Phage {j}"
        if phages and j < len(phages):
            phage_name = f"Phage {j} ({phages[j]['name']})"

        params[f"Initial Density (P0) - {phage_name}"] = {
            "type": "initial_P",
            "index": j,
        }

    params["Initial Resource Substrate (S0)"] = {
        "type": "initial_S",
    }

    # 6. Broadcast parameters — apply one value to ALL strains / phages at once.
    # Useful when strains share a parameter (e.g. WT and mutant have the same growth
    # rate) so the user sweeps a single control instead of one per strain.
    if config.n_bacteria > 1:
        for label, field, typ in [
            ("Growth Rate (ALL strains)", "growth_rates", "array1d_broadcast"),
            ("Bacteria-Resource Ratio (ALL strains)", "bacteria_to_resource_ratio", "array1d_broadcast"),
            ("Dormancy Rate (ALL strains)", "dormancy_rate", "array1d_broadcast"),
            ("Resuscitation Rate (ALL strains)", "resuscitation_rate", "array1d_broadcast"),
            ("Dormancy Diffusion Rate (ALL strains)", "dormancy_diffusion_rate", "array1d_broadcast"),
            ("Natural Death Rate dB (ALL strains)", "death_rate_B", "array1d_broadcast_or_none"),
            ("Dormant Death Rate dD (ALL strains)", "death_rate_D", "array1d_broadcast_or_none"),
            ("Immune Stim Rate (ALL strains)", "imm_stim_rate", "array1d_broadcast"),
            ("Immune Kill Rate (ALL strains)", "imm_kill_rate", "array1d_broadcast"),
            ("Immune Kill Rate Dormant (ALL strains)", "imm_kill_rate_D", "array1d_broadcast_or_none"),
            ("Initial Density B0 (ALL strains)", "growth_rates", "initial_B_broadcast"),
        ]:
            params[label] = {"type": typ, "field": field}
    if config.n_phages > 1:
        params["Phage Decay Rate (ALL phages)"] = {"type": "array1d_broadcast", "field": "phage_decay_rates"}
        if getattr(config, "phage_decay_Km", None) is not None:
            params["Phage Decay Km (ALL phages)"] = {"type": "array1d_broadcast", "field": "phage_decay_Km"}
        if np.ndim(getattr(config, "infected_nutrient_consumption", 0.0)) > 0:
            params["Infected-cell nutrient consumption (ALL phages)"] = {
                "type": "array1d_broadcast", "field": "infected_nutrient_consumption"}

    # Strain × phage matrix params — one value across every strain AND phage.
    if config.n_phages >= 1 and config.n_bacteria * config.n_phages > 1:
        for label, field in [
            ("Adsorption (ALL strains × phages)", "adsorption_rates"),
            ("Dormant Adsorption (ALL strains × phages)", "adsorption_rates_dormant"),
            ("Burst Size (ALL strains × phages)", "burst_sizes"),
            ("Latent Period (ALL strains × phages)", "latent_periods"),
            ("Dormant Latent Period (ALL strains × phages)", "latent_periods_dormant"),
            ("Dormant Adsorption Attenuation (ALL strains × phages)", "attenuation_rate"),
        ]:
            params[label] = {"type": "array2d_broadcast", "field": field}
        for field, lbl in [("hibernation_rate", "Hibernation Rate (I→H)"),
                           ("lytic_resumption_rate", "Lytic Resumption Rate (H→I)")]:
            if getattr(config, field, None) is not None:
                params[f"{lbl} (ALL strains × phages)"] = {"type": "array2d_broadcast", "field": field}

    # Antibiotic PD matrix params — one value across every strain × antibiotic.
    if config.n_antibiotics > 0 and config.pd_config is not None \
            and config.n_bacteria * config.n_antibiotics > 1:
        for label, field in [
            ("Emax (ALL strains × abx)", "abx_emax"),
            ("EC50 (ALL strains × abx)", "abx_ec50"),
            ("Hill (ALL strains × abx)", "abx_hill"),
        ]:
            params[label] = {"type": "pd_array2d_broadcast", "field": field}

    # ── Categorical (signal-function) sweeps ──────────────────────────────────
    # Model-wide signal selectors are discrete choices, not numerics. A categorical entry
    # carries its option set + the session key that selects it; param_sweeps sweeps by
    # setting that key and REBUILDING per option (the signal is resolved at build time).
    from pbisim_app.common import GROWTH_SIGNALS, DEATH_SIGNALS, LYSIS_SIGNALS, SIGNAL_OPTIONS
    _sig_opts = {o: o for o in SIGNAL_OPTIONS}
    params["Growth signal function"] = {
        "type": "categorical", "session_key": "int_growth_function",
        "options": {L: fn for L, (fn, _t) in GROWTH_SIGNALS.items()}}
    params["Death signal function"] = {
        "type": "categorical", "session_key": "int_death_function", "options": dict(DEATH_SIGNALS)}
    params["Lysis signal function"] = {
        "type": "categorical", "session_key": "int_lysis_function", "options": dict(LYSIS_SIGNALS)}
    params["Dormancy signal function"] = {
        "type": "categorical", "session_key": "int_dormancy_signal", "options": _sig_opts}
    params["Resuscitation signal function"] = {
        "type": "categorical", "session_key": "int_resuscitation_signal", "options": _sig_opts}
    params["Depth-diffusion signal function"] = {
        "type": "categorical", "session_key": "int_diffusion_signal", "options": _sig_opts}

    return params


_SWEEP_CATEGORY_ORDER = [
    "Bacterial", "Phage", "Immune", "Nutrient & environment", "Antibiotic",
    "OD & debris", "Initial conditions", "Structure & pre-run", "Signal functions", "Other",
]


def sweep_category(label: str, meta: dict) -> str:
    """Assign a sweepable parameter to a builder-style category so the UI can group
    them (mirrors the Interactive Simulator's parameter sections)."""
    t = meta.get("type", "")
    field = meta.get("field", "")
    if t == "categorical":
        return "Signal functions"
    if t in ("initial_B", "initial_B_broadcast", "initial_P", "initial_S"):
        return "Initial conditions"
    if t in ("prerun", "dimension"):
        return "Structure & pre-run"
    if t == "mutation":
        return "Phage" if field.endswith("_phage") else "Bacterial"
    if field.startswith("imm_") or "Immune" in label:
        return "Immune"
    if field in ("monod_constant", "carrying_capacity", "recycle_fraction", "s_in", "s_out",
                 "monod_constant_lysis", "lysis_floor",
                 "dormancy_monod_constant", "dormancy_carrying_capacity"):
        return "Nutrient & environment"
    if field in ("od_to_cfu_conversion_factor", "debris_kdis", "debris_u", "debris_v"):
        return "OD & debris"
    if field in ("growth_rates", "bacteria_to_resource_ratio", "dormancy_rate", "resuscitation_rate",
                 "dormancy_diffusion_rate", "death_rate_B", "death_rate_D"):
        return "Bacterial"
    if field in ("phage_decay_rates", "phage_decay_Km", "adsorption_rates", "adsorption_rates_dormant",
                 "burst_sizes", "latent_periods", "latent_periods_dormant", "attenuation_rate",
                 "hibernation_rate", "lytic_resumption_rate", "infected_nutrient_consumption"):
        return "Phage"
    if t in ("pk_array1d", "pd_array2d", "pd_array2d_broadcast") \
            or field in ("k_elim", "Vc", "abx_emax", "abx_ec50", "abx_hill"):
        return "Antibiotic"
    return "Other"


def categorize_sweep_params(params: dict) -> dict:
    """Group sweep-param labels by category, in a stable category order; each
    category's labels are sorted. Empty categories are dropped."""
    from collections import OrderedDict
    groups = OrderedDict((c, []) for c in _SWEEP_CATEGORY_ORDER)
    for label, meta in params.items():
        groups.setdefault(sweep_category(label, meta), []).append(label)
    return OrderedDict((c, sorted(ls)) for c, ls in groups.items() if ls)


def apply_sweep_parameter(val: float, meta: dict, config, initial_B, initial_P, initial_S, model_kwargs) -> tuple:
    """
    Applies a sweep value to the model configuration or initial conditions.
    Returns a new tuple of (config, initial_B, initial_P, initial_S, model_kwargs).
    Resizes initial_D if n_depth changes.
    """
    # Clone everything to prevent side-effects
    config = dataclasses.replace(config)
    initial_B = np.copy(initial_B)
    initial_P = np.copy(initial_P)
    model_kwargs = dict(model_kwargs)
    if "initial_D" in model_kwargs and model_kwargs["initial_D"] is not None:
        model_kwargs["initial_D"] = np.copy(model_kwargs["initial_D"])

    param_type = meta["type"]

    if param_type == "scalar":
        setattr(config, meta["field"], val)

    elif param_type == "dimension":
        # Compartment counts (n_depth, n_latent) must be integers >= 1. A sweep can
        # hand us a fractional or sub-1 value (e.g. a default log range starting at
        # 0.1); round to the nearest integer and clamp to 1 so the model never gets
        # 0 compartments (which crashes with an IndexError).
        field = meta["field"]
        new_val = max(1, int(round(float(val))))
        if field == "n_depth":
            new_depth = new_val
            config = dataclasses.replace(config, n_depth=new_depth)
            # Resize initial_D if present
            if "initial_D" in model_kwargs and model_kwargs["initial_D"] is not None:
                init_D = model_kwargs["initial_D"]
                total_D = np.sum(init_D, axis=0) # shape (n_bacteria,)
                new_init_D = np.zeros((new_depth, config.n_bacteria))
                for q in range(new_depth):
                    new_init_D[q, :] = total_D / new_depth
                model_kwargs["initial_D"] = new_init_D
        elif field == "n_latent":
            config = dataclasses.replace(config, n_latent=new_val)

    elif param_type == "array1d":
        arr = np.copy(getattr(config, meta["field"]))
        arr[meta["index"]] = val
        setattr(config, meta["field"], arr)

    elif param_type == "array1d_or_none":
        arr = getattr(config, meta["field"])
        if arr is None:
            arr = np.zeros(config.n_bacteria)
        else:
            arr = np.copy(arr)
        arr[meta["index"]] = val
        setattr(config, meta["field"], arr)

    elif param_type == "array1d_broadcast":
        arr = np.copy(getattr(config, meta["field"]))
        arr[:] = val  # same value for every strain / phage
        setattr(config, meta["field"], arr)

    elif param_type == "array1d_broadcast_or_none":
        n = config.n_phages if meta["field"].startswith("phage") else config.n_bacteria
        arr = getattr(config, meta["field"])
        arr = np.full(n, val, dtype=float) if arr is None else np.copy(arr)
        arr[:] = val
        setattr(config, meta["field"], arr)

    elif param_type == "initial_B_broadcast":
        initial_B = np.full_like(initial_B, val)

    elif param_type == "array2d":
        arr = np.copy(getattr(config, meta["field"]))
        arr[meta["index_row"], meta["index_col"]] = val
        setattr(config, meta["field"], arr)

    elif param_type == "array2d_broadcast":
        arr = np.copy(getattr(config, meta["field"]))
        arr[:, :] = val   # same value for every strain × phage
        setattr(config, meta["field"], arr)

    elif param_type == "pk_array1d":
        # pk_array1d always targets the antibiotic PKConfig, not PhagePKConfig.
        # Using phage_pk_config here when both are set would corrupt phage PK params.
        pk_config = dataclasses.replace(config.pk_config)
        arr = np.copy(getattr(pk_config, meta["field"]))
        arr[meta["index"]] = val
        setattr(pk_config, meta["field"], arr)
        config = dataclasses.replace(config, pk_config=pk_config)

    elif param_type == "pd_array2d":
        pd_config = dataclasses.replace(config.pd_config)
        arr = np.copy(getattr(pd_config, meta["field"]))
        arr[meta["index_row"], meta["index_col"]] = val
        setattr(pd_config, meta["field"], arr)
        config = dataclasses.replace(config, pd_config=pd_config)

    elif param_type == "pd_array2d_broadcast":
        pd_config = dataclasses.replace(config.pd_config)
        arr = np.copy(getattr(pd_config, meta["field"]))
        arr[:, :] = val
        setattr(pd_config, meta["field"], arr)
        config = dataclasses.replace(config, pd_config=pd_config)

    elif param_type == "initial_B":
        initial_B[meta["index"]] = val

    elif param_type == "initial_P":
        initial_P[meta["index"]] = val

    elif param_type == "initial_S":
        initial_S = val

    elif param_type == "mutation":
        # Origin→dest rate of a generator matrix (mutation_rates / transition_rates /
        # mutation_rates_phage / transition_rates_phage — meta["field"], default the
        # bacterial mutation matrix). Each column sums to 0 (mass conservation), so set
        # the off-diagonal entry and rebalance that column's diagonal.
        fld = meta.get("field", "mutation_rates")
        M = np.copy(getattr(config, fld))
        o, d = meta["origin"], meta["dest"]
        M[d, o] = val
        col = M[:, o].copy()
        col[o] = 0.0
        M[o, o] = -float(np.sum(col))
        config = dataclasses.replace(config, **{fld: M})

    elif param_type == "prerun":
        # Not a ModelConfig field — the app applies the pre-run (stationary_phase_ic)
        # before the treatment solve. Stash the swept duration for the sweep loop to
        # pick up (it pops "_t_prerun_override" before building the PBIModel).
        model_kwargs["_t_prerun_override"] = float(val)

    return config, initial_B, initial_P, initial_S, model_kwargs


def parse_comma_separated_series(text: str) -> list[float]:
    """
    Parses a comma-separated string of values into a list of floats.
    Supports linear, scientific notation, and spaces.
    """
    if not text.strip():
        return []
    parts = text.split(",")
    vals = []
    for p in parts:
        p_clean = p.strip()
        if p_clean:
            vals.append(float(p_clean))
    return vals


def pad_vectors(vectors: dict[str, list[float]]) -> tuple[dict[str, list[float]], list[str]]:
    """
    Pads all non-empty lists in vectors to the maximum length of any list.
    Padding is done by repeating the last value in the list.
    Returns the padded vectors and a list of warning messages detailing which lists were padded.
    """
    non_empty = {k: v for k, v in vectors.items() if len(v) > 0}
    if not non_empty:
        return {}, []

    max_len = max(len(v) for v in non_empty.values())
    padded_vectors = {}
    warnings = []

    for k, v in non_empty.items():
        if len(v) < max_len:
            pad_val = v[-1]
            padded_v = list(v) + [pad_val] * (max_len - len(v))
            padded_vectors[k] = padded_v
            warnings.append(f"Dose series for '{k}' was padded from length {len(v)} to {max_len} using the last value {pad_val:.1e}.")
        else:
            padded_vectors[k] = list(v)

    return padded_vectors, warnings
