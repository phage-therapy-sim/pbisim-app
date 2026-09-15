"""
test_system_prompt_sync.py — guard against drift between the system prompt and
the real pbisim API.

`prompts/system_prompt.md` hand-documents pbisim signatures so the agent
generates valid code. If pbisim renames or removes any of that surface, the
prompt silently goes stale and the agent emits broken code. These tests fail
loudly when that happens, so the prompt (and executor namespace) can be updated
in lockstep.

If a test here fails after a pbisim upgrade, update BOTH this file and
`prompts/system_prompt.md` (and `pbisim_app/executor.py` if the namespace
changed).
"""

import inspect

import pytest


# ── 1. Every import path the prompt + executor rely on must resolve ───────────

def test_prompt_import_paths_resolve():
    """The deep import paths quoted in system_prompt.md still exist."""
    from pbisim.builder import ModelBuilder            # §1, §7
    from pbisim.core.model import PBIModel             # §3, §7
    from pbisim.core.solver import solve_ode           # §4, §7
    from pbisim.pk.dosing import DoseSchedule, DoseEvent          # §6
    from pbisim.strains.builder import StrainSet, StrainDefinition  # §8
    from pbisim.pk.antibiotic import (                 # §8
        AntibioticDefinition,
        AntibioticSensitivity,
    )
    from pbisim.strains.genotypes import BinaryResistanceGenotypes  # §8


# ── 2. ModelBuilder constructor + fluent methods named in §1–§2 ───────────────

def _params(func):
    return set(inspect.signature(func).parameters)


def test_modelbuilder_constructor_signature():
    from pbisim.builder import ModelBuilder
    params = _params(ModelBuilder.__init__)
    for name in ("n_bacteria", "n_phages", "n_latent", "n_depth"):
        assert name in params, f"ModelBuilder.__init__ lost parameter '{name}'"


@pytest.mark.parametrize(
    "method",
    [
        "with_growth_rates",
        "with_phage_params",
        "with_mutations",
        "with_nutrient",
        "with_antibiotic",
        "with_growth_function",
        "with_sequential_growth",
        "with_smooth_efficiency_growth",
        "with_transition_function",
        "build",
    ],
)
def test_modelbuilder_methods_exist(method):
    from pbisim.builder import ModelBuilder
    assert callable(getattr(ModelBuilder, method, None)), (
        f"system_prompt.md uses ModelBuilder.{method}() but it no longer exists"
    )


def test_with_mutations_signature():
    """The prompt documents all four generator matrices + the per-locus shortcut."""
    from pbisim.builder import ModelBuilder
    params = _params(ModelBuilder.with_mutations)
    for name in ("mutation_rates", "transition_rates", "mutation_rates_phage",
                 "transition_rates_phage", "phage_resistance_rates"):
        assert name in params, f"with_mutations lost parameter '{name}'"


def test_strainset_generator_setters_exist():
    from pbisim import StrainSet
    for m in ("set_mutation_graph", "set_transition_graph",
              "set_phage_mutation_rates", "set_phage_transition_rates"):
        assert callable(getattr(StrainSet, m, None)), f"system_prompt.md uses StrainSet.{m}()"


def test_with_antibiotic_signature():
    from pbisim.builder import ModelBuilder
    params = _params(ModelBuilder.with_antibiotic)
    for name in ("name", "k_elim", "Vc", "emax", "ec50", "hill"):
        assert name in params, f"with_antibiotic lost parameter '{name}'"


# ── 3. PBIModel / solve_ode signatures named in §3–§4 ─────────────────────────

def test_pbimodel_signature():
    from pbisim.core.model import PBIModel
    params = _params(PBIModel.__init__)
    for name in ("config", "initial_B", "initial_P", "initial_S"):
        assert name in params, f"PBIModel.__init__ lost parameter '{name}'"


def test_solve_ode_signature():
    from pbisim.core.solver import solve_ode
    params = _params(solve_ode)
    for name in ("t_end", "dt"):
        assert name in params, f"solve_ode lost parameter '{name}'"


# ── 4. SimulationResult accessors named in §5 ─────────────────────────────────

def test_simulation_result_accessors():
    """`result.get(...)` and `result.sum_prefixes(...)` are the §5 contract."""
    from pbisim.core.solver import SimulationResult
    for method in ("get", "sum_prefixes"):
        assert callable(getattr(SimulationResult, method, None)), (
            f"system_prompt.md uses result.{method}() but it no longer exists"
        )


# ── 5. §13 stationary_phase_ic — B0 is required (the prompt example passes it) ─

def test_stationary_phase_ic_requires_B0():
    from pbisim import stationary_phase_ic
    params = _params(stationary_phase_ic)
    assert "B0" in params, "system_prompt.md §13 passes B0= to stationary_phase_ic"


# ── 6. §18 StrainSet surface — required StrainDefinition fields + to_config args ─

def test_strain_definition_required_fields():
    """The §18 StrainDefinition example sets exactly the required fields; if pbisim
    adds/removes a required field the example goes stale."""
    import dataclasses as dc
    from pbisim.strains.builder import StrainDefinition
    required = {f.name for f in dc.fields(StrainDefinition)
               if f.default is dc.MISSING and f.default_factory is dc.MISSING}
    documented = {
        "name", "growth_rate", "adsorption_rates", "adsorption_rates_dormant",
        "burst_sizes", "latent_periods", "latent_periods_dormant",
        "bacteria_to_resource_ratio", "dormancy_rate", "resuscitation_rate",
        "dormancy_diffusion_rate", "imm_stim_rate", "imm_kill_rate", "attenuation_rate",
    }
    assert required == documented, (
        f"StrainDefinition required fields changed — update system_prompt.md §18.\n"
        f"missing from prompt: {required - documented}\n"
        f"stale in prompt: {documented - required}"
    )


def test_strainset_to_config_required_args():
    """§18 documents these to_config args; a rename would break generated StrainSet code."""
    from pbisim.strains.builder import StrainSet
    params = _params(StrainSet.to_config)
    for name in ("n_latent", "n_depth", "phage_decay_rates", "imm_decay_rate",
                 "imm_stim50", "imm_kill50", "monod_constant", "recycle_fraction"):
        assert name in params, f"StrainSet.to_config lost parameter '{name}' (see §18)"


# ── 7. §2/§3 shape contracts the model kept getting wrong (eval-driven) ────────

def test_with_phage_params_decay_is_per_phage():
    """§2: phage_decay_rates has shape (n_phages,), NOT (n_bacteria, n_phages).
    The example teaches this; guard that the documented shape actually builds."""
    import numpy as np
    from pbisim.builder import ModelBuilder
    (ModelBuilder(n_bacteria=2, n_phages=1)
        .with_growth_rates([1.2, 1.1])
        .with_phage_params(
            adsorption_rates=np.array([[1e-8], [0.0]]),
            burst_sizes=np.array([[100], [100]]),
            latent_periods=np.array([[0.5], [0.5]]),
            phage_decay_rates=np.array([0.02]),   # (n_phages,)
        ))


def test_pbimodel_no_phage_takes_empty_initial_P():
    """§3: initial_P is required even with no phages — pass np.array([])."""
    import numpy as np
    from pbisim.builder import ModelBuilder
    from pbisim.core.model import PBIModel
    cfg = ModelBuilder(n_bacteria=1, n_phages=0).with_growth_rates(1.2).build()
    PBIModel(cfg, initial_B=np.array([1e7]), initial_P=np.array([]), initial_S=1.0)
