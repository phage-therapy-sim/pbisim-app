# CLAUDE.md — pbisim-app session context

This file is read at the start of every Claude Code session in this repository.

> **Workspace context:** `pbisim-app` is one of three packages in the
> phage-therapy-sim workspace. **Read the root
> [`../CLAUDE.md`](../CLAUDE.md) first** for workspace-wide conventions, and
> [`../ECOSYSTEM.md`](../ECOSYSTEM.md) for architecture and the cross-package API
> contracts. This file covers only pbisim-app specifics.

---

## Project identity

**pbisim-app v0.1.0** — AI-powered (Anthropic Claude) Streamlit interface for the
`pbisim` phage-bacteria simulation engine. Users can explore the engine interactively,
run parameter sweeps, design clinical trials, and ask an AI assistant to build and
explain simulations in natural language.

**Status:** active development, **deployed on Render** (Standard instance, auto-deploy
from `main`). **300 tests passing.** Depends on `pbisim>=1.0` (engine is now **2.0.0**)
**and `pbisim-fit>=0.1`** (now **1.0.0**)
— the fit integration is **LIVE** (the Calibration page runs pbisim-fit's NLS; lazy
import keeps the app torch-free). See §5.3 in ECOSYSTEM.md.

---

## Repository layout

```
pbisim-app/
├── pbisim_app/
│   ├── __init__.py          package docstring / launch instructions
│   ├── app.py     (~410)    THIN ENTRY: favicon/page_config, load styles.css,
│   │                        session init, sidebar, and a dispatch that imports the
│   │                        matching views/*.py and calls render(). (Was ~6050 lines;
│   │                        split 2026-07-20.)
│   ├── common.py  (~2300)   ALL shared helpers, the _ReproRecorder class, constants,
│   │                        and the st.number_input precision monkeypatch. Re-exports
│   │                        everything via __all__ (119 names) → a view needs only
│   │                        `from pbisim_app.common import *`.
│   ├── views/               One render() module per page (NOT Streamlit's reserved
│   │                        pages/): library, calibration, assistant, trials,
│   │                        dose_response, param_sweeps, simulator.
│   ├── static/              Bundled IBM Plex woff2 (served at /app/static/) + styles.css.
│   ├── viz_helper.py        plot_axis_controls + apply_axis_mpl/apply_axis_plotly
│   │                        (log/linear + axis-limit controls, both backends).
│   ├── agent.py   (~126)    SimulationAgent (wraps anthropic.Anthropic),
│   │                        AgentResponse, _parse_response() (regex extraction).
│   ├── executor.py(~155)    execute_code(): sandboxed namespace, exec()s
│   │                        generated code, captures stdout + matplotlib figures.
│   ├── sweep_helper.py      ModelConfig mutation helpers for 1D/2D parameter sweeps
│   │                        and dose-response sweeps; vector padding for MOI sweeps.
│   ├── trial_helper.py      IIV distribution factory, ClinicalTrial orchestration,
│   │                        Plotly KM + box plots, metric dataframe helpers.
│   └── fit_helper.py        Calibration/Phase-A: observable registry (CFU/PFU/OD/lum),
│                            data ingestion (normalize→pbisim-fit long format), overlay/RMSE.
├── prompts/
│   └── system_prompt.md     ~400-line hand-maintained pbisim API reference for the
│                            AI assistant. Sections 1–11 original; §12–17 added
│                            2026-06-23 covering new pbisim features.
├── tests/
│   ├── test_agent.py        5 tests — response parsing
│   ├── test_executor.py     10 tests — sandbox, capture, security, errors
│   ├── test_builder_modes.py 14 tests — BRG, StrainSet, cohort, phage-leak + immune + prerun + death-signal guards
│   ├── test_sweeps.py       14 tests — sweep_helper params, broadcast/coupled sweeps, OD trajectories, persistence
│   ├── test_self_healing.py 2 tests — self-healing loop and history rollback
│   ├── test_trial_features.py 5 tests — dose regimens, multi-arm, metrics, PK/PD trajectories
│   ├── test_scenarios.py     2 tests — scenario save/load round-trip (AppTest)
│   ├── test_parts.py         2 tests — parts save/load/export, host-tag (AppTest)
│   ├── test_calibration.py   4 tests — Calibration page config persistence, overlay, tuning (AppTest)
│   ├── test_fit_helper.py    9 tests — observable registry, ingestion, overlay/RMSE, tuning keys
│   └── test_repro_code.py    7 tests — reproduction script execs + full-config parity vs build (AppTest)
│   (test_system_prompt_sync.py  12 tests — sync guard, run after pbisim API changes)
└── pyproject.toml           entry point: pbisim-app = "pbisim_app.app:main"
```

---

## App pages

| Page | Description |
|---|---|
| Interactive Simulator | Three builder modes (Direct/ModelBuilder, Binary Genotypes/BRG, Custom Strains/StrainSet). Repeat-dosing regimen builder. Run button → plots + metrics. |
| Dose-Response Sweeps | Log/Lin dose range per agent, MOI scaling, vector padding warnings, color-coded trajectories. |
| Parameter Sweeps | 1D/2D sweeps over any ModelConfig field. Contour maps for 2D. n_depth resizing guard. |
| Clinical Trials & Cohorts | Full ClinicalTrial API integration: IIV, PretreatmentPhase, parallel arms, KM plots, metric distributions, CSV/NLME export. |
| Calibration | Upload experimental CSV → normalize to pbisim-fit long format (auto-detect + Monolix column-map) → **filter rows**, **regroup** by chosen variables, aggregate replicates (**raw / mean / median + percentile band**) → overlay the current model vs observations (group multiselect) with live RMSE; run the pbisim-fit NLS fit; **compare candidate models by AIC/BIC** (ΔAIC/ΔBIC parsimony panel). Extensible observable registry (CFU/PFU/OD/luminescence). |
| AI Assistant | Natural-language → pbisim code. Self-healing loop (up to 3 retries with history rollback). Dynamic model listing from `/v1/models`. |
| Library | Two sections: **💾 Scenarios** (save/load full-config snapshots) and **🧬 Parts** (composable bacteria/phages/antibiotics — save a current entity, load into config, host-tagged phages); each export/import as versioned JSON. (Tutorial presets + `presets.py`/`test_presets.py` removed 2026-07-10 — they tracked the pbisim tutorials, which change independently.) |
| Help | Curated in-app orientation (quick start, per-page guide, key concepts, troubleshooting) + the bundled `USER_GUIDE.md` rendered from disk at runtime (single source of truth). `views/help.py`; imported as `help_view` in the dispatch to avoid shadowing the builtin. |

---

## How the AI assistant works (data flow)

```
Streamlit UI (app.py)
   → SimulationAgent.ask(user_message)                     agent.py
   → Claude (model dynamically selected), system = prompts/system_prompt.md
   → _parse_response() → AgentResponse(code, narrative, assumptions, raw_text)
   → executor.execute_code(code) → ExecutionResult(success, figures, stdout, error)
     ↓ on failure (up to 3 times):
     → agent.ask(traceback) → corrected code → execute_code()
     → on final failure: history rolled back to pre-request state
   → UI renders narrative + figures (+ optional code/assumptions/stdout)
   → conversation history kept in agent.history (multi-turn)
```

- **No tool-use, no streaming.** Claude returns markdown; code is pulled from a
  ```` ```python ```` block by regex (`agent.py:_parse_response`).
- **API key:** `ANTHROPIC_API_KEY` env var, else Streamlit sidebar password input.
- **Model:** dynamically selected in sidebar; default pin is `claude-sonnet-4-6`.

---

## Package-specific conventions

1. **The system prompt is the contract with the engine.** `prompts/system_prompt.md`
   hard-codes `pbisim` signatures — it is NOT generated. If the `pbisim` public
   API changes, update this prompt in lockstep (coordinate via the integration
   role). Never let the prompt invent method names.
2. **Generated code may use either pre-loaded sandbox names OR `import` statements.**
   The executor sandbox (`executor.py`) pre-loads: `np/numpy`, `plt/matplotlib`,
   `ModelBuilder`, `PBIModel`, `solve_ode`, `DoseSchedule`, `DoseEvent`, `StrainSet`,
   `StrainDefinition`, `BinaryResistanceGenotypes`, `BacterialStrain`, `PhageStrain`,
   `Antibiotic`, `PhagePKConfig`, `AntibioticDefinition`, `AntibioticSensitivity`,
   `stationary_phase_ic`, `time_to_clearance`, `time_to_log_reduction`, and (if
   pbisim.trial is installed) `TrialRunner`, `IIVSpec`, `VirtualPopulation`,
   `default_metrics`, `LogNormal`, `Normal`, `Uniform`, `Fixed`, `ClinicalTrial`,
   `TreatmentArm`, `PretreatmentPhase`. Add to the namespace AND the prompt
   together when new surface is needed.
3. **Honor the shared engine contracts** (see ECOSYSTEM.md §3.2): the prompt must
   keep instructing the model to use `result.sum_prefixes("B","D")` for **CFU**
   (culturable = active B + dormant D; infected `I` / hibernating `H` cells lyse rather
   than plating — use `("B","D","I","H")` only for *total live load*), `result.get('B0')`
   for individual series, `np.maximum(x, 1.0)` before `log10`, and to set `initial_S` on
   `PBIModel` (not `solve_ode`). The app plots/metrics default to CFU=B+D (with total and
   active selectable); Calibration's observation model defaults CFU to B+D and threads the
   chosen compartments to pbisim-fit via `NLSConfig.cfu_compartments`. **Engine-side
   `default_metrics`/clinical clearance still default to `B+D+I+H`** pending the ECOSYSTEM
   §3.2 integration ruling — a known residual gap on the Clinical Trials page.
4. **Resistance seeding rule** (in the prompt): never start the resistant strain
   at exactly 0 — use `initial_B[resistant] = max(mutation_rate * initial_B[0], 10.0)`,
   because mutation flux → 0 when growth → 0 under nutrient limitation.
5. **Sandbox is research-grade only.** It strips `open/exec/eval/compile` and
   captures stdout/figures, but is not safe for untrusted/public input. Do not
   weaken it; do not deploy publicly without real isolation (Docker /
   RestrictedPython).
6. **Never commit API keys.** Keys come from env or UI at runtime only.
7. **`sweep_helper.py` pk_array1d** always targets the antibiotic `PKConfig` — never
   `phage_pk_config`. Both are on `ModelConfig` and must not be confused.

---

## Test commands

Activate the project env first (`conda activate pbisim`), then:

```bash
# Full suite — expected: 86 passed
python -m pytest tests/ -q

# By file:
python -m pytest tests/test_executor.py -q        # 10 tests
python -m pytest tests/test_agent.py -q           # 5 tests
python -m pytest tests/test_builder_modes.py -q   # 14 tests (BRG, StrainSet, cohort, phage-leak + immune + prerun + death-signal)
python -m pytest tests/test_sweeps.py -q          # 14 tests (params, broadcast/coupled sweeps, OD trajectories, persistence)
python -m pytest tests/test_self_healing.py -q    # 2 tests
python -m pytest tests/test_trial_features.py -q  # 5 tests (dose regimens, multi-arm, metrics, PK/PD trajectories)
python -m pytest tests/test_calibration.py -q     # 4 tests (config persistence, overlay, manual tuning)
python -m pytest tests/test_fit_helper.py -q      # 9 tests (observable registry, ingestion, overlay/RMSE)
python -m pytest tests/test_repro_code.py -q      # 7 tests (repro script execs + full-config parity vs build)
```

After any pbisim API change, also run:
```bash
python -m pytest tests/test_system_prompt_sync.py -v
```

## Running the app

```bash
python -m streamlit run pbisim_app/app.py
# or, once installed:  pbisim-app
```

**Desktop window (optional):** `pbisim_app/desktop.py` (console script `pbisim-app-desktop`,
extra `.[desktop]` → pywebview) starts the *same* local Streamlit server and wraps it in a
native OS window — additive, no app-logic or Render-deploy changes (browser mode unchanged;
updates flow normally, nothing is frozen). Groundwork for a future local "power-user mode"
(in-app scripting + local-LLM assistant), only safe when everything runs locally.
- **Verified 2026-07-28** on this Linux box: launcher works end-to-end (server up → native
  Qt window → clean shutdown). **BUT the webview's browser engine must be modern enough for
  current Streamlit** (needs `Object.hasOwn`, ES2022 / Chromium ≥ ~93). This dev box's old
  PyQt5 QtWebEngine (Ubuntu 20.04) is too old → the window opens but the UI stays blank (a JS
  error in the page, NOT a Python exception, so it can't auto-fall-back). Windows (WebView2 =
  evergreen), macOS (WKWebView), and up-to-date Linux (recent WebKit2GTK/Qt) are fine.
- Escape hatches: `PBISIM_APP_BROWSER=1` forces the browser; a missing pywebview / native
  backend also falls back to the browser automatically. Linux native window needs a system
  webview lib (`apt install gir1.2-webkit2-4.1 python3-gi`, or a recent PyQt6-WebEngine).

## Access gate (deployment)

`pbisim_app/auth.py` adds an optional shared-credential sign-in (`require_login()`,
called early in `app.py`). It is **active only when `APP_PASSWORD` or
`APP_PASSWORD_HASH` is set in the environment** — set it in the Render dashboard
(Environment vars; never commit), leave it unset locally. Optional `APP_USERNAME`.
Basic gate to keep the public out, NOT enterprise SSO; the AI exec sandbox is still
research-grade, so only share the password with trusted people.

**Persistent "stay logged in" (2026-08-03).** The gate used to store auth only in
`st.session_state`, which is in-memory + tied to the WebSocket session — so every
reconnect (laptop sleep, network blip, backgrounded tab), OOM restart, or redeploy
bounced the user back to the password page (a tier upgrade doesn't fix this; it only
removed free-tier spin-down). Fixed with a **signed, expiring cookie**: on sign-in a
`pbisim_auth` cookie = `<expiry>.<hmac>` is set (via `streamlit-cookies-controller`,
now a core dep); each fresh session reads it back **server-side** with
`st.context.cookies` (reliable, no JS-component load flash) and restores the login.
The HMAC secret is derived from the configured credential (rotating the password
invalidates all cookies) or `APP_AUTH_SECRET`; lifetime is `APP_AUTH_TTL_HOURS`
(default 168 = 7 days). Sign-out removes the cookie **and** sets `_cookie_suppressed`
for the session (the connection's request cookies are fixed for its lifetime, so the
stale cookie would otherwise re-login until the next reconnect). Cookie set/remove is
queued (`_pending_cookie`/`_forget_cookie`) and applied on the next clean render so
`st.rerun()` can't abort the component. Degrades to session-only if the component is
absent. The cookie only proves "this browser passed the password" — not an identity
token. Tests in `tests/test_auth.py` (token sign/verify/expiry/tamper/rotate, cookie
restore, invalid-cookie, sign-out suppression).

---

## Known gaps / next steps

- **pbisim-fit integration is LIVE** (Calibration page): CSV ingest → overlay →
  NLS fit via `refine_nls` (lazy import, torch-free); role-based fit-parameter
  table + MAP priors + statement DSL; **additive-B0** (bacteria-dose /
  `free_initial_conditions` estimate / `cfu[0]` fallback); **NONMEM/Monolix
  dose-row import**; real-builder reuse for tuning. *Remaining (non-urgent):* route
  ingestion through pbisim-fit's `EventTable.from_csv` to retire the app's
  hand-rolled dose parser + get covariate pass-through (see
  `../pbisim-fit/APP_INTEGRATION_NOTES.md`); optional `to_model_config`.
- No streaming, no tool-use, brittle regex parsing, no rate limiting, no UI tests.
- Model pin lives in `agent.py` `_MODEL` (currently `claude-sonnet-4-6`); keep it
  current and re-run `tests/test_system_prompt_sync.py` after any pbisim upgrade.
- Preset `script_code` strings for type="single" presets (01–10, 13) are reference
  only — they are not executed. Any API mismatch there is cosmetic but should be fixed.

## Done this session (2026-09-15) — generator matrices: phenotypic transitions + phage mutation/transitions

Owner: "the app doesn't incorporate transition matrices applied to bacteria or phages." Confirmed:
the engine has FOUR mass-conserving generators (`M[dest, origin]`, diagonal = −outflow) —
`mutation_rates` (bacteria, applied to the DIVISION flux → vanishes when growth stops),
`transition_rates` (bacteria, first-order switch on `B` itself, h⁻¹, growth-independent),
`mutation_rates_phage` (applied to the lysis yield) and `transition_rates_phage` (on free `P`) —
and the app only ever wired the first (per-locus μ / `int_transitions` mutation graph). The other
three were unreachable from every builder mode. Now:
- **Shared edge lists** (all `int_*` → captured by Scenarios): `int_bact_transitions`,
  `int_phage_mutations`, `int_phage_transitions` (`GENERATOR_GRAPH_KEYS` in common.py; the
  historically-named `int_transitions` stays the bacterial MUTATION graph). Generic
  `render_rate_graph_editor(names, state_key, prefix)` replaces the two duplicated inline
  editors (`render_mutation_graph_editor` is a thin wrapper); `rate_matrix_from_graph` /
  `edges_from_rate_matrix` / `graph_dict_from_edges` / `brg_genotype_labels` (mirrors the
  engine's `strain_labels`, asserted equal in a test) are the pure helpers.
- **UI:** one mode-aware expander below the builder columns (`render_generator_matrix_editors`
  in simulator.py, so Calibration's embedded builder gets it too): bacterial nodes = strain
  names (Direct/StrainSet) or genotype labels (BRG), phage nodes = phage names; shown when
  ≥2 strains or ≥2 phages, auto-expanded when edges exist.
- **Build wiring:** Direct → ONE `with_mutations(...)` call carrying all four kwargs;
  StrainSet → `set_transition_graph` / `set_phage_mutation_rates` / `set_phage_transition_rates`;
  BRG → post-build `cfg.<field> = M` via `apply_generator_matrices` (BRG `to_config` HARD-CODES
  zeros for these three AND forwards `**extra_config_kwargs` to `ModelConfig`, so passing them
  would raise "multiple values" — post-build assignment is the app-side path, recorded in the
  repro via the new `rec.assign`). Repro parity asserted for all three modes.
- **Sweeps:** the `"mutation"` sweep type gained a `field` (default `mutation_rates`); entries
  for all four matrices (phage ones under the Phage category), mass-conserving apply.
- **Fit table:** targets for configured (nonzero) transition/phage-generator edges (bounds
  families added); **and a real fix** — `build_param_spec_v2` now derives every touched
  generator column's diagonal (`_rebalance_generator_diagonals`: `M[j,j] = −Σ_{i≠j} M[i,j]`,
  dynamic when any entry is freed/mapped). The pre-existing `mutation_rates[i,j]` targets freed
  only the off-diagonal, so a fit created/destroyed cells instead of moving them.
  `_apply_config_to_state` writes fitted generators back into the edge lists
  (`write_back_generator_matrices`; only graphs already configured; BRG skips mutation).
- **Snapshot:** "Evolution & phenotypic switching" section (nonzero matrices only).
- **AI:** `configure_simulator` accepts `transition_graph` / `phage_mutation_graph` /
  `phage_transition_graph` (any mode). `system_prompt.md` `with_mutations` block rewritten — it
  was WRONG (said `[i,j] = i→j` and "diagonal ignored / auto-normalised"; the engine is
  `M[dest, origin]` with an explicit negative diagonal, `_to_nn` does no normalisation) — and now
  documents all four matrices + `with_transition_function` + the StrainSet setters; prompt-sync
  guard extended (`with_mutations` signature, StrainSet setters).
- Docs: USER_GUIDE §4.4b + parameter table row; Help "Resistance" concept. Tests +12
  (builder_modes ×5 incl. a growth=0 test proving transition ≠ mutation, sweeps, nls_fit ×2,
  ai_configure, prompt_sync ×2). **300 passing.** `with_transition_function` (a callable)
  stays scripting-only.

## Done this session (2026-08-07) — phage 2-compartment PK + central-compartment display

Owner: (1) show central-compartment phage concentration wherever phage is plotted, (2) wire
pbisim's 2-compartment phage PK (peripheral tissue) which the app lacked (antibiotic PK already
had k12/k21). Engine recap: phage PK states are `P` (infection site, plotted), `Pc` (central/
blood, amount → conc = Pc/Vc), `Pp` (peripheral tissue amount, no volume in engine); 2-comp is
`PhagePKConfig(k12=, k21=)`.

- **2-compartment (A):** added **k12** (central→peripheral) + **k21** (peripheral→central) inputs
  to the phage-PK section in all 3 builder modes (Direct/BRG/StrainSet; shown when PK mode ≠ None,
  default 0 = 1-compartment). Build wires them via new `_phage_pk_k12_k21(phages)` → all three
  `PhagePKConfig(...)` calls in common.py (arrays when any PK phage has k12>0, else None; engine
  requires both-or-neither). `build_series` adds a **`Pp{j}`** peripheral series (opt-in) when k12>0,
  and relabelled the central series "(blood)" → "(central/blood conc.)".
- **Central display throughout (B):** the Simulator already showed `Pc{j}` (central conc.). Added it
  to the other surfaces: **Dose-Response** + **Parameter Sweeps** (all 3 loops) collect a
  `central_phage_trajectories` (Σ Pc/Vc) and offer a **"Central phage — blood (PFU/mL)"** trace
  (infection-site relabelled "Free phage — infection site"); **Clinical Trials** PK/PD outputs gained
  a central option (`prefixes=("Pc",)` with new `divide_by=Vc` on `plot_pkpd_trajectories_plotly`;
  uses the first PK phage's Vc, caption if phages' Vc differ).
- Helpers (common.py, exported): `phage_pk_enabled`, `central_phage_total` (Σ Pc/Vc conc.),
  `peripheral_phage_total` (Σ Pp amount), `phage_uses_two_compartment`, `_phage_pk_k12_k21`.
- k12/k21 default via `.get(...,0.0)` (no default-dict churn; captured in scenarios once the PK
  section renders). Tests: `test_redesign::test_phage_pk_central_and_peripheral_series_and_helpers`
  (series+helpers), `test_builder_modes::test_phage_two_compartment_pk_wired_into_config` (app build
  → config.phage_pk_config.k12/k21 + Pc/Pp states), `test_sweeps::
  test_dose_response_shows_central_phage_when_pk_enabled`. Suite **270 green**. AI system_prompt has
  no phage-PK section at all (pre-existing) — left out of scope. **NOT committed** (UI change).

## Done this session (2026-08-06) — model-aware fit-parameter table (one knob per quantity)

Owner flagged real internal inconsistencies in the Calibration fit table: it exposed BOTH
coupled/redundant controls with invisible precedence — e.g. `growth_rates[1]` AND `fitness_cost`
(a BRG model derives the resistant growth from the cost), and `fit_initial_cfu` in the table AND
the per-arm B₀-source radio. Confirmed precedence in pbisim-fit: `_apply_fitness_cost` sets
`growth_rates[1:] = growth_rates[0]·(1−fitness_cost)` at solve time **only when fitness_cost is
freed**, OVERWRITING a freed `growth_rates[1]` (fitness_cost wins; a fixed fitness_cost value is
ignored). A freed `fit_initial_cfu` overrides the "First observation" B₀ source. Nothing was
broken, but which control wins was invisible → confusing.

Fix (owner chose "model-aware, one knob each" via AskUserQuestion): `available_targets` gained a
`builder_mode` arg (calibration passes the CHOSEN model's `int_builder_mode`):
- **B₀ is never a table row** — dropped `fit_initial_cfu`. B₀ is the per-arm **B₀-source** radio
  whose "Estimate (shared/per-arm)" options ARE how you free it (`free_initial_conditions`).
  (`fit_initial_pfu` stays — phage-titre estimate has no per-arm analog.)
- **BRG:** exposes WT `growth_rates[0]` + `fitness_cost`/`init_resistant_fraction`; HIDES the
  derived per-genotype `growth_rates[i≥1]` (its knob is the fitness cost).
- **Direct/StrainSet:** exposes each independent `growth_rates[i]`; HIDES the BRG-only
  `fitness_cost`/`init_resistant_fraction` (they'd be a second, overwriting control).
UI caption under the fit table states which model's params are shown + that B₀ is the §4 control.
Note: mappings/DSL still bind any path (a mapped path is applied regardless of table membership),
so advanced reparameterization is unaffected in Direct/StrainSet; in BRG the resistant growth is
intentionally the fitness cost (independent resistant growth ⇒ use Direct/StrainSet). Tests:
`test_available_targets_reflects_builder_mode` (per-mode catalog); the two B₀ tests now use the
"Estimate (shared)" radio instead of freeing a table row. Full suite 267 green. **NOT committed**
— held for owner local testing (standing calibration-UI rule).

## Done this session (2026-08-06) — CI test-suite fix (absolute AppTest paths for Streamlit ≥1.59)

The GitHub Actions **test-and-deploy** workflow (`.github/workflows/deploy.yml`) has been
failing at "Run test suite" on every push for a while — while the Render deploy succeeded
(autoDeploy is independent of Actions). Root cause: CI does `pip install -e .` (no pins), so
it gets the **latest** streamlit; a recent version (≈1.59+) changed `AppTest.from_file(<relative
path>)` to resolve against the **calling test file's directory** (`tests/`), not the CWD — so
`AppTest.from_file("pbisim_app/app.py")` looked for `tests/pbisim_app/app.py` → `FileNotFoundError`
in **all ~100 AppTest tests**. Reproduced locally in a throwaway venv with latest deps
(streamlit 1.61.1 / numpy 2.5 / pandas 3.0): 107 failed → after the fix, **267 passed** (so the
newer numpy/pandas are otherwise fine).

Fix: every test now passes an **absolute** path (from_file uses an absolute path verbatim on all
versions) — `APP = str(Path(__file__).resolve().parents[1] / "pbisim_app" / "app.py")` — replacing
the module-level `APP = "pbisim_app/app.py"` and every inline `AppTest.from_file("pbisim_app/app.py"…)`
across all test files. (No dependency pinning — the app itself runs fine on latest streamlit on
Render; only the test path was version-fragile.) Verified on both streamlit 1.58 (dev) and 1.61 (CI-like).

## Done this session (2026-08-06) — Calibration fit runs in a SEPARATE PROCESS (real nav-freeze fix)

The freeze-on-navigate persisted even with no auto-poll, because the fit ran in a **thread**:
measured, a threaded pbisim fit leaves the main thread at only ~22% CPU, and the *real*
pbisim-fit NLS is heavier — the GIL starves Streamlit, so navigating freezes until the fit
ends. No Streamlit-side polling trick fixes a GIL-bound thread. Owner chose (AskUserQuestion)
to keep navigation, so the fit now runs in a **separate process**.

- `pbisim_app/fit_process.py` (NEW, Streamlit-free): `run_fit_in_process(queue, payload)` runs
  `run_nls_fit_v2` and puts `("done", {map, ci, fitted_config})` or `("error", msg)` on the
  queue. It EXTRACTS results to picklable primitives in the child (the posterior may not
  pickle; a `ModelConfig` does — verified).
- `_start_fit_process(payload)` (calibration.py) forks via `mp.get_context("fork")` +
  `Process(..., daemon=True)`. **Fork, not spawn/forkserver**: spawn re-imports the Streamlit
  app script in the child (breaks it); fork inherits memory, so the payload (config, dataset,
  covariate resolvers) is NOT pickled to pass — only the *result* is pickled back through the
  queue. Parent pre-imports `nls_fit`/`pbisim_fit.refinement.nls` before forking so the child
  runs no imports (import lock = the usual fork-in-a-thread deadlock source). n_arm_jobs=1 → no
  nested workers (a daemon process can't have children). Fork is the Linux/Render default
  (deploy = python:3.11-slim, no fork DeprecationWarning; 3.13 dev box warns but works).
- `_monitor_fit_job` polls the queue non-blocking each render (harvest done/error; detect a
  process that exited with no result = crash/OOM). Still **no auto-poll loop** — a running fit
  shows a STATIC banner + manual **🔄 Refresh status** (+ Stop, which `terminate()`s the
  process = real cancellation). Navigation is a plain rerun; the process has its own GIL so the
  app stays fully responsive. Result harvested on Refresh/return. `_harvest_fit_done` reads
  `_job["result"]` (the queue dict) instead of a live posterior.
- Removed the thread worker `_fit_worker` + `import threading`.
- Tests: monitor tests now inject a plain `queue.Queue` with a `("done"/"error", …)` item
  (deterministic); the 9 real fit tests fork a real subprocess and are driven by `_settle_fit`.
  Full suite green. **Risk noted:** fork-in-a-multithreaded-process can (rarely) deadlock;
  pre-importing in the parent mitigates the common case. If it ever bites on deploy, the
  fallback is a synchronous spinner.

## Done this session (2026-08-06) — Calibration fit progress: no auto-poll (superseded same day)

Saga of the running-fit progress display (all on the Calibration page; the fit itself always
ran fine in a daemon thread — this was only about how the UI *observes* it):
1. `sleep(1)+st.rerun()` loop at page top → monopolised Streamlit's script runner → app froze
   ~1 min, nav radio **desynced to the wrong page**, finished fit never harvested. Also blew
   CI's per-`.run()` wall-clock timeout (one `.run()` spanned the whole fit).
2. `st.fragment(run_every=1.0)` → fixed on-page ticking, but the fragment timer kept firing /
   re-triggering `st.rerun(scope="app")` after the user navigated away → **froze again**.

Measured the actual cost: a threaded pbisim fit leaves the main thread at ~172 ticks/s (vs ~930
idle, worst stall 40 ms) — i.e. the **thread does NOT freeze the app**; every freeze was the
*polling mechanism* fighting Streamlit's execution model. So: **removed all auto-poll.** A
running fit now shows a **STATIC banner + a manual `🔄 Refresh status` button** (+ Stop).
Navigation is a plain rerun (nothing to starve → no freeze, no desync); the result is harvested
by `_monitor_fit_job` on the next render — when the user clicks Refresh or simply returns to the
page. Trade-off: the banner doesn't tick live and the result isn't *pushed* — the user pulls it
(Refresh/return). `_harvest_fit_done` extracted; records `calib_fit_error` instead of throwing
if a converged result can't be read.

**Test impact:** the fit is async, so a single `.run()` after "Run NLS fit" returns while the
thread runs. `_settle_fit(at)` (test_models.py) re-runs the app until `fit_job` clears — added
after each of the 9 fit-run clicks (each `at.run()` acts like a Refresh). New `test_calibration
.py::test_running_fit_does_not_freeze_page_or_navigation` (running job shows banner, page still
renders, nav away works). Suite green.

Possible future upgrade if live auto-update is wanted without the freeze: run the fit in a
SUBPROCESS (frees the main thread fully) — but the polling-vs-Streamlit fragility is the same
regardless of thread/process, so the manual-Refresh model is the robust baseline.

## Done this session (2026-08-06) — coupled sweep can mix categorical + continuous

Bug: a **Coupled (linked)** parameter sweep that paired a **categorical** signal-function
param (e.g. Lysis signal function) with a **continuous** one (e.g. Initial Phage Density)
didn't work — the coupled loop applied EVERY selected param via `apply_sweep_parameter`,
which mutates an already-built config, but a signal-function choice resolves at BUILD time
(session key + `*_kwargs()` helpers), so the categorical one was silently ignored. (The 1D
categorical sweep already handled this by set-key-then-rebuild; coupled didn't.)

Fix (`views/param_sweeps.py`): split coupled labels into categorical vs continuous. Per
step, if any categorical is present, set each categorical **session key** to that step's
option and **rebuild** via `build_nominal_config_from_gui()`; then apply the continuous
params via `apply_sweep_parameter` on the rebuilt config; solve. Categorical session keys
saved before / restored after in a `finally` (no residue on the live model). UI: a categorical
coupled param uses **one dropdown PER STEP** — a "number of steps" `st.number_input`
(`pc_catn_{i}`) sets how many `st.selectbox`es (`pc_catopt_{i}_{j}`) render, so an option can
**REPEAT** across steps (e.g. constant×3 then nutrient×3 paired with 0,1e3,1e5,0,1e3,1e5 =
the full 2×3 factorial as a linked series). A single `st.multiselect` was tried first but it
can't repeat an option (each disappears once chosen) — that's why per-step boxes. Each keyed
selectbox is sanitized (pop a stale/invalid persisted value before render so it can't raise).
Compute reads `[pc_catopt_{i}_{j} for j in range(pc_catn_{i})]`, drops invalids. Continuous
params stay comma-text (`pc_series_{i}`); the "same length" check pairs them. Trajectory-label + summary formatting made
value-type-aware (`_short_sweep_label` / `_fmt_sweep_value`; the coupled summary applies the
`{:.2e}` numeric format only to numeric columns — a categorical column holds option-name
strings). Pure-numeric coupled sweeps unchanged (no rebuild, no restore).

**2D + categorical:** NOT supported — 2D axes are both numeric (heatmap grid). Added a UI
guard: picking a signal-function param for either 2D axis shows a warning and **disables**
Run 2D, pointing to 1D (compare options) or Coupled (pair with a continuous param). (Repro
path untouched — still the "1D only" message for coupled/2D.) Test: `test_sweeps.py::
test_coupled_sweep_mixes_categorical_and_continuous` (mix lysis-signal + P0 via the dropdown;
asserts option strings + densities in the summary, session key restored, empty selection
rejected). Suite green.

## Done this session (2026-08-05) — Calibration fit-monitor survives navigation

Bug: run an NLS fit, navigate away from Calibration and back → the progress banner
and the **Run NLS fit** button vanish; the fit *looks* dead. Diagnosis: the fit was
**never killed** — it runs in a `daemon` thread (`_fit_worker`) writing only into the
`fit_job` holder (a plain `session_state` dict), so it keeps running across navigation
(and its completion lands in the holder). But the running banner **and** the
done→result harvest lived deep inside section 5c's parameter-table `try` block, which
rebuilds from arm-selection / `fit_targets*` state that is deliberately *not* persisted
across navigation → on return that block could raise before reaching the poll code,
hiding the banner and leaving a *completed* fit un-harvested until 5c rendered cleanly.
Fix: **decoupled the monitor from 5c.** New module-level `_monitor_fit_job()` (calibration.py)
runs at the TOP of `render()` (right after the fit-config re-seed, before the dataset
gate): harvests done/error into `calib_fit_result`/`calib_fitted_config`/`calib_overlay_result`,
and while a fit is **running** renders ONLY the banner + **Stop** and self-polls
(`sleep(1)`+`rerun`) FROM THE TOP — it does NOT fall through to rebuild the rest of the
(heavy) page each cycle. Section 5c now only renders the Run button (disabled while
running) + the result table; its old poll/harvest block (~70 lines) deleted.

Why top-poll (not a bottom-of-page poll): the page below the banner is static stored
data during a fit, so re-rendering it every second buys nothing and is the main
per-cycle cost — rebuilding the fit table (`available_targets`, `build_config_from_model`).
`AppTest.run()` measures a single wall-clock `timeout` across the WHOLE rerun chain of
one `.run()`, so a fit's many 1-s poll cycles × a full-page rebuild on a slow 2-core CI
runner blew that budget → **the CI test suite failed (deploy gate) even though it passed
locally.** Rendering only the banner per cycle makes each poll cheap and keeps CI under
budget. It's also strictly more race-safe than a bottom poll: the running branch reruns
unconditionally each cycle, so the instant the thread flips status the next top pass
harvests it (a completed fit can't be stranded). Tests: `test_calibration.py` +2
(`test_background_fit_result_harvested_by_top_monitor` injects a stub `status='done'`
job with **no dataset loaded** → harvested + `fit_job` cleared;
`test_background_fit_error_surfaced_by_top_monitor`). Full suite 265 green.

Caveat (separate issue): a full **WebSocket reconnect / OOM restart / redeploy** starts a
*new* session (or a new process), so the in-memory `fit_job`/thread genuinely can't be
recovered there — that's the cookie-auth/disconnect class of problem (2026-08-03), not
this navigation fix.

## Done this session (2026-08-04) — lysis floor φ_min (residual efficacy)

pbisim `45e389c` added `config.lysis_floor` (φ_min ∈ [0,1], default 0) + `with_lysis_function(
…, phi_min=)`: the lysis signal is rescaled `φ = φ_min + (1−φ_min)·signal`, so frac_lysis keeps
some phage efficacy as nutrients deplete (φ→φ_min at S→0) instead of vanishing at stationary
phase. φ_min=0 = legacy; 1 = constant lysis; no-op for constant lysis. App:
- `lysis_kwargs()` returns `lysis_floor` (read from `int_lysis_floor`, only when frac_lysis).
- Direct build → `with_lysis_function(fn, Ks_lysis=, phi_min=)`; BRG/StrainSet → `extra_kwargs
  ["lysis_floor"]` (to_config). UI: "Lysis floor (φ_min)" input next to Ks_lysis (shown for
  frac_lysis + nutrient-tracking growth). Defaults + preset/scenario load + snapshot(n/a — lysis
  not shown). Sweep: scalar "Lysis Floor (φ_min)" (Nutrient & environment). Fit table: estimable
  "Lysis floor (φ_min)" when frac_lysis (`_FAMILY` (0,1)). Test: floor reaches config all modes.

## Done this session (2026-08-03) — infected_nutrient_consumption is now per-phage

pbisim `196d1d8` made `infected_nutrient_consumption` a **per-phage** `(n_phages,)` array
(a phage×host property) instead of a system-wide scalar (scalar still broadcasts). App
updated to match:
- **UI:** per-phage input in each phage card (Direct/BRG/StrainSet, keys `phg_infnut_{i}` /
  `brg_phg_infnut_{idx}` / `ss_phg_infnut_{idx}`), next to the dormant-adsorption attenuation.
  Removed the old single nutrient-env scalar input (simulator + calibration compact block →
  captions pointing to the phage cards).
- **Build:** `growth_nutrient_kwargs` assembles the `(n_phages,)` array from the phage dicts
  and passes it to `with_nutrient` (Direct) / `to_config` (BRG/StrainSet); ALWAYS an array
  when phages exist (so `config.infected_nutrient_consumption` is per-phage-sweepable). Legacy
  `int_infected_nutrient_consumption` scalar still seeds any phage without its own value.
- **Sweep:** per-phage `array1d` entries + an ALL-phages broadcast, recategorized to **Phage**.
- **Fit table:** per-phage `infected_nutrient_consumption[j]` targets (was one scalar), gated on
  the config field being an array. Snapshot row relabeled "× per phage". Scenario/preset phage
  seeding carries the field.
- Tests updated: per-phage build (all modes), per-phage sweep, per-phage fit target. 262 passing.

## Done this session (2026-08-03) — real code editor for the Scripting page

- **`streamlit-code-editor` (Ace) replaces `st.text_area` in the Scripting cells** — Tab-indent,
  auto-indent, syntax highlighting, line numbers, and a ▶ Run button bound to **Ctrl/Cmd+Enter**
  (native `st.text_area` can do none of these: Tab moves focus, no auto-indent/highlight).
  `views/scripting.py`: `_HAS_CODE_EDITOR` import guard; `code_editor(..., response_mode="blur",
  buttons=_EDITOR_BUTTONS)` per cell; edits persist on blur. `_apply_editor_response(cid, resp)`
  folds the buffer into the plain `script_src_{cid}` key and runs **once per submit `id`** (the
  component re-emits its last value every rerun, so a bare `type=='submit'` check would re-execute
  on unrelated reruns; empty-`type` responses never clobber the source).
- **Graceful fallback**: absent the component, cells fall back to `st.text_area` (functional, no
  niceties). Source is `script_src_{cid}` in both paths, so run/add/delete logic is identical.
- **Packaging**: new `[scripting]` optional extra (`streamlit-code-editor>=0.1.20`, frontend-only,
  lazy-imported); **Dockerfile** now `pip install -e '.[scripting]'` so the editor is present on
  the deploy wherever scripting is enabled. Auto-deploys from `main` (Docker, `autoDeploy: true`);
  the changed install line + app source invalidate that layer → the component installs
  automatically (no manual cache clear; the pbisim/pbisim-fit git-ref layers stay cached).
- Tests: `test_scripting.py` rewritten editor-agnostic (source via `script_src_{cid}` /
  `script_cell_ids`, not `at.text_area`) + `_apply_editor_response` submit-dedup unit test +
  text_area-fallback test. **257 passing.** Caveat: AppTest can't drive real keystrokes, so
  Tab/auto-indent/Ctrl-Enter were verified by owner in a live browser.

## Done this session (2026-08-03) — upstream sync: new growth signal, infected-nutrient draw, AIC/BIC, prompt

Reflected recent pbisim / pbisim-fit updates in the app (A→F sequence):
- **A — `smooth_efficiency_monod`** growth signal (Monod K blends efficient→inefficient via a Hill
  of S; a differentiable, better-conditioned diauxie). `GROWTH_SIGNALS` + `_GROWTH_NUTRIENT`;
  `growth_nutrient_kwargs` emits `monod_K_low`/`θ`/`hill`; Direct build → `with_smooth_efficiency_growth`
  (recorded), BRG/StrainSet via `to_config`; dedicated UI inputs; snapshot rows; auto in the categorical
  sweep. `_GROWTH_SMOOTH_EFF` const.
- **B — `infected_nutrient_consumption`** (latent-I cells draw down substrate; default 0 = legacy).
  UI input in the nutrient-environment section; forwarded via `growth_nutrient_kwargs` → `with_nutrient`
  (Direct) / `to_config` (BRG/StrainSet); snapshot row.
- **C — AIC/BIC model comparison on Calibration.** `_compute_overlay` now returns the pooled log10
  `residuals` (+ `n_resid`); `fit_helper.model_information_criteria` / `compare_fit_models` delegate to
  pbisim-fit's `information_criteria`/`compare_models` (local fallback for older installs). New
  **"Compare models (AIC / BIC)"** expander: snapshot a candidate (residuals + free-param count k) →
  ranked table with ΔAIC/ΔBIC; warns when candidates were overlaid on different n. Buttons in
  `_FIT_NOPERSIST`; store `fit_model_comparison`.
- **D — NLS solver default.** Confirmed the app builds `NLSConfig(...)` WITHOUT `fit_method`/`x_scale`,
  so it inherits pbisim-fit's new BDF + `x_scale='jac'` defaults automatically; added a caption on the
  §5c fit control.
- **E — AI prompt sync.** Added `infected_nutrient_consumption` + a growth-signal-functions block
  (`with_growth_function`/`with_sequential_growth`/`with_smooth_efficiency_growth`, Gompertz nutrient
  note) to `prompts/system_prompt.md`; extended `test_system_prompt_sync.py` method guard.
- **F — docs/memory:** Help key-concepts + Calibration row; this note; memory `upstream-sync-2026-08`.
- Tests: +A (build+repro), +B (all-modes+default-off), +C (helpers + overlay residuals + panel e2e),
  +E (guard). Full suite green (run before commit).

## Done this session (2026-07-27) — EventTable dose-parser migration + model-config snapshot

- **EventTable migration (anti-drift).** `fit_helper.parse_dose_rows` now delegates the
  obs/dose split + canonical event model to pbisim-fit's **`EventTable`** (builds a
  canonical events frame from the app's mapped columns → reads `et.doses`), instead of the
  app's hand-rolled EVID loop. Dose-target vocabulary is sourced from pbisim-fit via
  `fit_helper.dose_targets()` (`DOSE_TARGETS` kept as a synced module constant). Used the
  df constructor `EventTable(events=df)` — the in-memory equivalent of `from_csv` (the app
  already holds a dataframe + its own column-map UI, so `from_csv`'s file/dialect handling
  isn't needed). Observation ingestion (normalize/filter/group/aggregate) is app-specific
  and stays. Verified byte-identical behavior (all dose/NONMEM/nls tests pass).
- **Model-config snapshot.** New `common.render_model_snapshot(container, *, snapshot=None)`
  renders a sectioned, mode-agnostic summary of the **fully-resolved** config (builds it
  via `build_nominal_config_from_gui` inside `model_config_context`, so Direct/BRG/StrainSet
  all show the same fields, no drift): growth & nutrient environment, death & dormancy,
  phage, immunity, OD/debris, ICs, solver. Click-to-view **`📋 Show model config`** toggle
  (gates the build) added to the Simulator (live draft) and `page_model_selector` (chosen
  frozen Model). No more tab-hunting to see a model's config.
- Tests: `test_models.py::test_model_snapshot_renders_on_simulator_and_pages`;
  `parse_dose_rows` migration covered by existing tests. **204 tests passing.**

## Done this session (2026-07-26) — B0/dose calibration overhaul (pbisim-fit additive-B0)

pbisim-fit finished the additive-B0 / EventTable (NONMEM/Monolix) work (its
`APP_INTEGRATION_NOTES.md` is the hand-off); app side implemented in 3 pieces:

- **B — ratio-mode builder.** `render_model_builder(inoculum_mode="magnitude"|"ratio")`.
  Calibration renders `"ratio"`: per-strain B0 → "Initial ratio (relative)", absolute
  inoculum totals hidden (phage P0, BRG `brg_eq_total_B`), BRG equilibrium checkbox KEPT
  (sets the resistant fraction). Simulator unchanged (`"magnitude"`).
- **C — B0 source + additive dose.** Per-arm conditions gained a B0 radio (First
  observation [+ noise warning] / Shared / Per-arm). Fixed modes → `build_dataset` emits
  `DoseRecord(target="bacteria", unit="cfu")` (or `pretreatment_inoculum` for pre-run
  arms); first-obs emits neither (pbisim-fit's `cfu[0]` fallback + warning).
- **A — NONMEM dose-row import.** Column-map gained EVID/AMT/unit selectors;
  `fit_helper.parse_dose_rows()` parses dose rows (EVID=1; observable = target
  compartment) into per-arm records that gate the manual per-arm fields and are emitted
  verbatim by `build_dataset(arm_doses=)` (overriding the manual fields for covered
  targets). Overlay respects per-arm `moi_unit`.

Tests: parse_dose_rows, additive-dose + imported-dose override (test_nls_fit), ratio
labels + B0-mode + NONMEM gating (test_calibration). **201 tests passing.** See memory
[[calibration-additive-b0-dose]]. Committed + pushed this session (engine initial_S fix
`abf3b5f` too).

## Done this session (2026-07-25) — Calibration manual-tuning = the real model builder

Owner audit: the Calibration **Manual parameter tuning** panel (a hand-curated
`STRAIN_TUNABLES`/`PHAGE_TUNABLES`/per-mode list in `fit_helper.py`) had drifted —
**BRG mode showed no dormancy kinetics at all** (strain block skipped in BRG; dormant
attenuation still showed because it's a *phage* tunable), and `dormant_od_fraction`
(+ signals, pseudolysogeny, …) were missing. Root cause = the curated list is
structurally drift-prone. **Fix (owner-approved): reuse the real builder.**

- **Extracted** `simulator.py`'s Tab-1 body (mode selector + growth/death signals +
  Direct/BRG/StrainSet) verbatim (pure move, 0 content diff) into module-level
  **`render_model_builder()`** (seeds strains/phages/antibiotics from `int_*`, returns
  `builder_mode`). The Simulator's Tab 1 now calls it.
- **Calibration** manual-tuning renders the SAME `render_model_builder()` under a
  `fit_show_builder` toggle (outside any expander — the builder has its own), so every
  parameter and all 3 builder modes are present and can never drift. A compact "Global &
  structural" block keeps the non-Tab-1 params (n_latent, nutrient env, OD/debris) and
  now includes the **`dormant_od_fraction`** input (`fit_edit_dorm_od`).
- **Deleted** the curated `*_TUNABLES` + `entity_param_key` (fit_helper + common
  `__all__`/imports + app.py + their 2 tests). fit-apply now also pops the builder
  widget keys (`_BUILDER_WIDGET_PREFIXES`) so applied fits re-seed the builder inputs.
- Tests: `test_calibration.py::test_calibration_embeds_full_builder_all_modes` (all 3
  modes render; Direct dormancy exposes str_sleep/str_wake); updated
  `test_manual_tuning_edits_model_directly` + `test_globals_and_debris_*`. **193 passing.**
- **NOT committed** — held for owner local testing (standing "test locally first" on
  calibration UI). Note: `render.yaml plan: standard` from the prior turn is also staged.

## Done this session (2026-07-24) — unified role-based fit-parameter table

Owner-approved redesign of the Calibration §5c fit UI (the layered free?-checkbox +
separate **Shared parameters** + **Reparameterization** panels confused users: hidden
precedence — "if a param is set free in the mapping but not free in the first table,
what wins?"). **Now ONE table, one Role per row.**

- **Unified parameter table** (`fit_targets_df`): each model parameter has a **role**
  SelectboxColumn — **Fixed** (held at value) · **Free** (estimated 1:1, uses
  bounds/prior) · **Derived** (`= expression` of a θ). Plus an `expression` column for
  Derived rows. Deleted the separate "Shared parameters (quick)" expander and the
  Reparameterization **mappings** table (`fit_map_df`) entirely — mappings now live on
  the target rows. Kept the **Custom parameters (θ)** table (`fit_thetas_df`).
- **Share helper** (one-click): multiselect rows → **Share →** sets them to Derived
  with a common auto-named θ (`shared1`…) and appends that θ. = the way to tie
  parameters together. Self-clears (`_share_clear`).
- **Assembly** (`views/calibration.py`): role → `free=(role=='Free')` target; Derived +
  expression → a `{path, expr}` mapping (validated vs θ names). Backend
  `build_param_spec_v2`/`run_nls_fit_v2` UNCHANGED (targets still carry `free` internally).
- **DSL** (`nls_fit.parse_fit_spec`): now returns **(targets_df, thetas_df, errors)** —
  the separate map_df is gone; `map <path> = <expr>` sets that row's role=Derived +
  expression. `serialize_fit_spec` prefixes a `# available parameters` comment header so
  paths are always to hand (addresses "user doesn't know the <path> names").
- Tests updated to the role model (`_free_targets`/`_derive_targets` helpers) + new
  `test_share_helper_ties_rows_to_one_theta`. **194 tests passing.**
- **NOT committed** — held for owner local testing (per the standing "test locally
  first" on the fit-spec/redesign work). See memory [[nls-od-link-on-config]] (v3 section).

## Done this session (2026-07-20) — visual redesign (branch `feature/redesign`)

Adopting the owner-approved "scientific instrument" mockup (Claude Design). Decided
to **stay in Streamlit** for the public launch — the real launch blockers (the AI
`exec()` sandbox is research-grade/unsafe for public input, plus auth/storage and
concurrency/cost) are framework-independent; keep the compute/AI core UI-agnostic so a
future frontend swap stays cheap. Work is on `feature/redesign` (pushed), **143 tests**.

- **Pass A — identity/theme:** self-host IBM Plex Sans+Mono woff2 in `pbisim_app/static/
  fonts/`, served via `enableStaticServing` at `/app/static/` (verified 200 font/woff2 on
  a headless boot); `config.toml` + a full rewrite of the injected CSS block (light+dark)
  to muted-teal `#0d7a68` / warm-paper `#faf9f5` / 6px bordered cards / mono data values;
  dropped the rainbow gradient-text headings and gradient buttons; de-emoji'd titles/
  headers/buttons (kept one 🐍); φ-mark sidebar brand.
- **Pass B — structure:** primary/secondary button hierarchy (Run = `type="primary"`);
  Interactive Simulator results header bar (solver + runtime meta + outcome badge) and a
  Peak Phage Titre tile; moved the `show_code`/`show_assumptions` toggles out of the
  sidebar onto the AI page; units on all core physical-parameter labels; expanders styled
  as cards.
- **Pass C — visualization:** new `pbisim_app/viz_helper.py` (`plot_axis_controls` +
  `apply_axis_mpl`/`apply_axis_plotly`, log/linear/log-log + axis limits, plotly log10-range
  gotcha handled). Wired into the sim bacterial+phage (mpl) and Dose-Response + Param-Sweep
  trajectories (plotly). Hardcoded `semilogy`/`yaxis_type="log"` removed from those.
- Tests: `tests/test_redesign.py`. **Merged to `main`** (feature/redesign deleted). Metric
  tiles now on every result page (Sim results header+badge+peak-phage, Trials per-arm cure
  rate, Calibration pooled RMSE+R², Dose-Response + Param-Sweep summaries).
- **Post-merge follow-ups:** AI settings collapsed into a sidebar `st.expander`; unit labels
  standardized (`h⁻¹`/`mL·h⁻¹`, dropped confusing "(r)"/"dB"/"dD"/"Y" symbols); **all
  app-owned plots migrated to Plotly** for consistency (repro script + AI-generated figures
  stay matplotlib by design). **Dark mode deferred** (light-only for launch; toggle hidden,
  dark CSS dormant).
- Deploy note: auto-deploys from `main`, but the earlier pbisim engine change still needs a
  **Clear build cache & deploy** on Render (cached pip layer).

## Done this session (2026-06-23 continued)

- **Adsorption default 2e-9 → 1e-8** for WT phage in all preset parameter dicts and
  script_code strings (`presets.py`); list first-elements updated per-strain; collateral-
  sensitivity preset ([1e-9, 5e-9]) left intact; IIV distribution lines left intact.
- **`system_prompt.md`** example updated: `adsorption_s=2e-9` → `adsorption_s=1e-8`.
- **`tests/test_presets.py`** fallback default updated: `2e-9` → `1e-8`. 48 tests pass.
- **`USER_GUIDE.md`** written: task-oriented GUI guide (12 sections, 5 worked workflows,
  parameter reference tables, troubleshooting). See `USER_GUIDE.md` at repo root.
- **`README.md`** written: repo front door with features, install, run/stop, pages table,
  doc links, test command. See `README.md` at repo root.
- **Immunity module full audit and fix** (commit `a397edd`):
  - Preset Tutorial 04: removed `adaptive_stimulation_rate/decay_rate/max/delay`
    (invented by scaffold, not in pbisim); added `immune_module`, `imm_stim_rate`,
    `imm_stim50`, `innate_decay_rate`, `imm_initial` with correct pbisim names.
  - `load_preset_to_state()`: backward-compat translation of `adaptive_decay_rate` →
    `innate_decay_rate`, `"adaptive"` module → `"innate"`; new keys `int_imm_stim_rate`,
    `int_imm_stim50`, `int_imm_initial`.
  - Direct / BRG / StrainSet simulation paths: all now pass `imm_stim_rate`,
    `imm_stim50`, correct `immune_module`, and `initial_Imm` (was inadvertently set
    to `imm_max` value).
  - Immunity UI (Tab 2): module selector fixed to `["innate", "hill"]` (removed invalid
    `"adaptive"`); `imm_stim_rate`, `imm_stim50`, `initial_Imm` widgets added; `imm_max`
    made conditional on hill module.
  - Repro code (Direct): updated with all new immunity parameters.
  - Repro code (BRG): was entirely missing immunity; now emits `brg.to_config(...)` with
    all immunity kwargs when enabled.
  - Repro code (StrainSet): was hardcoded `imm_stim_rate=0.0, imm_kill_rate=0.0`; now
    reads session state; `ss.to_config()` now includes immunity params, actual phage
    decay rates, and correct `n_depth`.
- **Clinical trial Combo arm guard** (same commit): `st.warning` emitted when Combo arm
  contains only phage or only antibiotic doses (Combo = monotherapy → identical KM
  curves). Explains the "strange outputs" the user observed.

## Done this session (2026-07-17) — signal functions (dormancy/resus/growth) + dormancy Ks

- **BUG: nutrient+density dormancy raised a ValueError.** The app sent the dormancy/
  resuscitation signal as `"nutrient_and_density"`, but pbisim's key is `"nutrient+density"`
  (`_DORMANCY_SIGNALS`). Fixed via `canonical_signal()` (also translates legacy stored
  values) + `SIGNAL_OPTIONS = ["constant","nutrient","density","nutrient+density"]` — the
  **"constant"** option was also missing. Direct-mode dormancy/resuscitation selectors now
  offer all four; build path canonicalises before calling `with_dormancy`.
- **Configurable dormancy nutrient half-saturation (engine + app).** Added
  `dormancy_monod_constant` to pbisim (`ModelConfig` field, `_monod_signal` uses it, else
  falls back to `monod_constant`; `with_dormancy(dormancy_monod_constant=…)`). App exposes
  a "Dormancy nutrient half-saturation (Ks)" input (Direct mode, shown for nutrient signals;
  0 = inherit growth Ks). Engine: **1144 tests pass** (+1).
- **Comprehensive growth signals.** New "Growth signal" selector (constant / nutrient (Monod)
  / density (logistic) / nutrient+density) mapped to pbisim growth functions via
  `growth_nutrient_kwargs()`, wired into all three builder modes (replaces the bare
  track-nutrients checkbox). Compatibility handled: when a growth signal freezes S,
  `compat_dormancy_signal()` / `dormancy_compat_kwargs()` pin nutrient-independent dormancy
  functions (the engine rejects nutrient dormancy with a frozen S); density-based dormancy
  now gets a `dormancy_carrying_capacity` even under Monod growth. Verified all 4 signals ×
  3 modes build. Repro code emits the growth function + nutrient config.
- Tests: `test_builder_modes.py::test_dormancy_signals_and_growth_signals`,
  `pbisim/tests/test_dormancy.py::test_dormancy_monod_constant_overrides_growth_ks`.
  **App: 74 tests pass.**
- **Follow-up (docs discipline):** the `dormancy_monod_constant` engine change was made
  without asking + went undocumented — corrected: documented in pbisim `API_REFERENCE.md`
  + `tutorial_03_dormancy.ipynb` (markdown-only; owner re-executes/redeploys), noted in the
  workspace CLAUDE.md. Also exposed the **pre-existing** `dormancy_carrying_capacity`
  (density dormancy threshold) as a Direct-mode input (`str_dcc_{i}`, shown for density
  signals; 0 = inherit growth K) — it was previously only auto-inherited. See memory
  `ask-before-changing-pbisim-engine`.

## Done this session (2026-07-16) — OD trajectories in parameter sweep

- **Parameter sweep now plots OD trajectories** when the OD/debris module is enabled
  (previously only the dose-response sweep did). Added `od_trajectories` collection to
  the 1D and Coupled compute loops (`result.get_od()`), stored in `param_sweep_result`,
  and an "Optical Density" chart rendered after the CFU trajectories in both. Omitted
  when debris is off; 2D (heatmaps only) is unaffected. Test in `test_sweeps.py`.
  **73 tests pass.**

## Done this session (2026-07-16) — mutation-rate persistence + coupled/broadcast sweeps

- **Direct-mode mutation rate reverted to 1e-7 on navigation.** The 2^m-shortcut
  mutation widget (`direct_mu_{j}`) read its `value=` from `direct_phg_res_mu_{j}`
  (sic — `direct_phg_mu_{j}`), a key that was NEVER written; the actual value lives in
  `direct_phg_res_rates`. So when the widget key was dropped on navigation the input
  re-defaulted to 1e-7 and overwrote a user's 0. Fixed to seed `value=` from the
  persisted `direct_phg_res_rates` list. Regression test in `test_builder_modes.py`.
- **Sweeping shared / linked parameters** (parameter-sweep feature):
  - **Broadcast params** (`sweep_helper`): `get_sweep_parameters` now emits
    `… (ALL strains)` / `(ALL phages)` entries (types `array1d_broadcast` /
    `array1d_broadcast_or_none` / `initial_B_broadcast`) that apply one value to every
    strain/phage — e.g. sweep a growth rate shared by WT + mutant with one control.
  - **Coupled (linked) sweep**: a third `ps_sweep_type` mode. Pick several parameters
    (multiselect `pc_labels`) and give each a value series of equal length
    (`pc_series_{i}`); at step k, value[k] of every parameter is applied together —
    e.g. (dormancy, resuscitation) = (0,1),(0.5,0.5),(1,0). Summary has a column per
    parameter; trajectories labelled by the tuple; metrics vs step index. Controls +
    result persist across navigation (added `pc_` to the persist prefixes).
  - Tests: `test_sweeps.py` (broadcast unit + coupled AppTest + mismatch guard).
- **72 tests pass.**

## Done this session (2026-07-15) — sweep CONFIG persistence + zero-dose default

- **Sweep controls now survive navigation** (previously only the *results* did). The
  dose-response `dr_sweep_*` widgets and the parameter-sweep `p1_*`/`p2_*`/`ps_*`
  widgets got dropped on navigation (Streamlit discards un-rendered widget keys).
  Added reusable `reseed_widget_config()` / `save_widget_config()` (shadow the
  selections into `dr_sweep_config` / `param_sweep_config` plain dicts, re-seed before
  render). For the parameter sweep the previously-unkeyed `sweep_type` (radio) and the
  1D `ps_1d_min/max/steps/spacing` inputs were keyed so they can persist; the 1D range
  widgets re-autoscale when the swept parameter changes (`_ps_1d_last_param` guard pops
  them so a persisted key doesn't pin a stale range).
- **Dose-response default phage series now `0, 1e3, 1e5, 1e7, 1e9`** (adds the zero-dose
  control). Test `test_dose_response_shows_od_trajectories_when_enabled` updated.
- Test: `test_sweeps.py::test_sweep_configs_survive_navigation`. **69 pass.**

## Done this session (2026-07-15) — pre-run collapse, sweep persistence, BRG calibration clash

- **Pre-run made CFU/OD scale very low (reported as a dose-response OD bug).** Root
  cause: `stationary_phase_ic` applies `death_rate_B` throughout the pre-run, but once
  nutrients exhaust growth stops — so with a death rate and no dormancy the culture
  declines for the whole pre-run and the treatment starts from a decimated inoculum
  (reproduced: death 0.5 → OD 1.16 no-prerun vs 0.03 after 24 h pre-run). It is the
  engine's documented behavior, not an app bug, but was silent. Added
  `prerun_collapse_fraction()` / `warn_if_prerun_collapses()` and surfaced the warning
  in the main sim run, the dose-response sweep and the parameter sweep when the pre-run
  leaves <10% of the inoculum (B+D). (Dormancy avoids it — persisters survive.)
- **Dose-response + parameter sweep results now survive navigation.** Both drew results
  inside the Run button's `if`, so navigating away wiped them. Split compute-from-render:
  the button stores results in `dr_sweep_result` / `param_sweep_result`; a separate block
  renders from the stored data every run (alive until re-run).
- **Calibration BRG clash + missed params** (see the 2026-07-15 tuning entry below for
  the panel details): BRG strain kinetics live on `int_brg_base_*` (not the per-strain
  dicts), so BRG strain edits were silent no-ops — the panel is now mode-aware. Mutation
  rate (μ) + fitness cost added per phage in BRG; dormant adsorption + adsorption_r added.
- Tests: `test_sweeps.py` (+2 persistence, +1 collapse warning), plus the earlier
  calibration/fit_helper additions. **68 pass.**

## Done this session (2026-07-15) — OD-debris propagation, structural/global tuning, save-calibrated

- **OD/debris now propagates into the Calibration overlay.** The overlay always
  simulated *with* debris (build_nominal includes it), but the OD *observable* was
  computed as biomass/link and ignored the debris state. `predicted_observable()`
  gained `use_model_od`: when the OD/debris module is on, OD comes from the model's
  debris-inclusive `result.get_od()` (uses `od_to_cfu_conversion_factor`) instead of
  the simple biomass/link scaling. The overlay's link input is then replaced by a note
  pointing to the OD/debris params in the tuning panel.
- **Manual calibration now exposes global & structural parameters** (a "Global &
  structural" block in the tuning panel): n_latent (latent compartments), K, Ks;
  nutrient S₀ / recycle / s_in / s_out (when nutrient tracking is on); the OD/debris
  params (od_to_cfu, u, v, k_dis) when the module is on; and per-strain dormancy depth
  (Q) added to the dormancy row. All edit the live `int_*` session keys directly.
- **Save the calibrated model** (section 6): edits are already live in the Interactive
  Simulator (same dicts), so "apply to builder" is automatic; added a **💾 Save
  calibrated config as Scenario** button (reuses `dump_state_to_scenario`) to persist
  the whole config to the Library for reload. `fit_save_scenario` (a button) is in
  `_FIT_NOPERSIST` — persisting a button key pre-sets it and makes `st.button` raise;
  the fit_config re-seed loop now also scrubs stale non-persistable keys.
- Tests: `test_fit_helper.py` (+1 debris-OD), `test_calibration.py` (+1 globals/debris/
  save-scenario). **66 pass.**

## Done this session (2026-07-15) — overlay persistence, builder-mode param parity, fuller tuning

- **Calibration overlay now survives navigation** (item 1): the overlay was drawn
  inside `if st.button(...)`, so leaving the page and returning (without re-clicking)
  wiped the plot. Split compute-from-render: the button computes and stores the plot
  data in `st.session_state.calib_overlay_result`; a separate block renders it on
  every run while present. Stays alive until re-run or the dataset is cleared.
- **Builder-mode parameter parity** (item 2): closed the "read-but-not-editable" gaps
  found by auditing all three builder UIs vs `build_nominal_config_from_gui`:
  - `bacteria_to_resource_ratio` — was read (default 1e9) but had no widget in
    **Direct** (`str_ratio_{i}`) or **StrainSet** (`ss_str_ratio_{i}`); added to both
    (BRG already had `int_brg_base_ratio`).
  - Pseudolysogeny (`hibernation_rate_s/r`, `lytic_resumption_rate_s/r`) — read into
    `PhageStrain` but no widgets in **BRG** (`brg_phg_{hib,res}_{s,r}_{idx}`) or
    **StrainSet** (`ss_phg_...`); added to both (Direct already had them).
- **Manual calibration exposes the fuller parameter set** (item 3): `STRAIN_TUNABLES`
  now growth / bacteria_to_resource_ratio / death_rate_B / initial_B; a conditional
  `STRAIN_DORMANCY_TUNABLES` row (dormancy/resuscitation/diffusion/dormant-death) shown
  per strain when dormancy is on; `PHAGE_TUNABLES` adds phage_decay_Km + attenuation_rate
  (on top of burst/latent/decay + the mode-aware adsorption editor).
- Tests: `test_calibration.py` (+1 overlay persistence), `test_builder_modes.py`
  (+1 all-mode ratio/pseudolysogeny), `test_fit_helper.py` registry update. **64 pass.**

## Done this session (2026-07-15) — dose-response zero-dose phage leak

- **Dose-Response Sweeps: a swept phage dose of 0 still suppressed the bacteria.**
  Same class of bug as the clinical-trial phage-leak fix: the model always seeds
  phage from the phage config's `initial_P` (default 1e6) at t=0, independent of the
  dose. The sweep overrides only the dose *events*, so `dose=0` still started with
  1e6 PFU and crushed the culture. Fix (`app.py`, Dose-Response Sweeps page): for the
  duration of the sweep, zero the swept phages' `initial_P` (saved/restored in the
  `finally` alongside `int_doses`) so the swept dose fully controls phage exposure;
  `dose=0` now means no phage. Antibiotics need no analog (no baseline inoculum — all
  via dose events). Verified: `dose=0` nadir stays ~1e7, `dose=1e8` eradicates.
  Regression test `tests/test_sweeps.py::test_dose_response_zero_phage_dose_does_not_suppress`.
  **61 tests pass.**

## Done this session (2026-07-14) — Calibration Phase B + config persistence

- **Calibration config now survives navigation** (commit `07315d1`): Streamlit
  drops a widget's key from `session_state` when the widget isn't rendered on a
  rerun, so leaving Calibration to change the model reset all the filter/grouping/
  statistic selections. Fix: shadow the `fit_*` widget selections into a plain
  `fit_config` session key (survives) and re-seed the widgets from it at the top of
  the page, before they render. Buttons + the file-uploader are excluded from the
  shadow (they can't be re-seeded). Regression test `tests/test_calibration.py`.
- **Phase B — manual parameter tuning** on the Calibration page (section 5, a
  "🎛 Manual parameter tuning" expander). Edits the model's **actual absolute
  parameter values** (like the ModelBuilder), *not* multipliers — the widgets read
  from / write to the live `int_strains`/`int_phages` dicts and session keys, so an
  edit **is** the model (no separate apply step) and is savable as a Part.
  - Exposed: global K / Ks; per strain growth_rate + initial_B; per phage burst /
    latent / decay; and **adsorption**, which is builder-mode specific — Direct &
    Custom-Strains store it in the pairwise `ads_{i}_{j}` session keys (edited per
    strain×phage pair), Binary-Genotypes on the phage dict `adsorption_s`. The link
    factor (od_to_cfu / rlu_per_cell) is keyed per-observable (`fit_link_{obs}`).
  - `fit_helper.STRAIN_TUNABLES`/`PHAGE_TUNABLES`/`ADSORPTION_PHAGE_KEYS` +
    `entity_param_key()` (mode-aware key resolution). Pure/no-Streamlit → unit-tested.
  - `fit_edit_*` widgets are excluded from the `fit_config` shadow (the dicts are
    authoritative + already persistent; re-seeding a stale copy would clobber edits).
  - Tests: `test_fit_helper.py`, `test_calibration.py`. **60 tests pass.** (Earlier
    in this session an initial ×multiplier version shipped, then was replaced per
    owner feedback: multipliers were inconvenient and silently missed Direct-mode
    adsorption, which lives in `adsorption_rates`/`ads_{i}_{j}`, not `adsorption_s`.)

**Known gaps / still open after this session:**
- CS/CR (26 `cr_*` fields on `PhageStrain`/`BacterialStrain`) not exposed in BRG UI —
  complex feature, deferred. Document in USER_GUIDE as "advanced scripting only".
- `immune_module="custom"` has no UI support — deferred.
- StrainSet repro code: `n_depth` calculation uses `dormancy_depth` key which may
  differ from actual session state key; edge case, not a blocker.
- Test count: 48.

## Done this session (2026-07-07) — clinical trial phage-leak fix

- **CRITICAL BUG FIXED — Control arm was secretly receiving the phage inoculum**
  (`app.py`, Clinical Trials & Cohorts page). Root cause: the crossover
  `ClinicalTrial` shares `initial_conditions` across all arms and differentiates
  them only by `dose_schedule` (`ClinicalTrial._apply_arm` never touches ICs). The
  app seeded phage via `initial_P` on the shared `base_cfg.initial_conditions.P`, so
  **every** arm — Control and Antibiotic-Only included — started with the full phage
  inoculum (default `1e6`/mL), which eradicated the bacteria in ~6 h. Every arm
  looked identical and the untreated control "cured" all patients (the strange
  outputs the user reported).
  - **Fix:** phage is the intervention, so it is now delivered per-arm. The trial
    starts every arm at **zero free phage** (`base_P = np.zeros_like(init_P)`) and
    injects the configured inoculum as a **t=0 phage bolus** only into the Phage-Only
    and Combo arms (verified numerically identical to seeding `initial_P`). Control
    and Antibiotic-Only are now genuine no-phage arms.
  - Arm-existence guards + Combo "identical to monotherapy" warning updated to treat
    the inoculum as phage presence (`_has_phage = inoculum or nominal phage dose`).
  - `base_initial_P` passed to `run_trial_simulation` also zeroed (model-factory
    fallback consistency).
  - Verified post-fix (30 patients, IIV on growth): Control 0/30 eradicated
    (bacteria → ~1e9), Antibiotic-Only 1/30 (weak monotx), Phage/Combo 30/30.
- **Regression test added**: `test_builder_modes.py::test_trial_control_arm_has_no_phage`
  — asserts Control keeps `P≈0` throughout and bacteria grow, while the phage arm
  delivers phage and eradicates. **49 tests passing** (was 48).

**Known limitation introduced:** the IIV option "Initial Phage Density" (`ic.P`) now
acts on a zeroed baseline in trials, so it no longer varies the inoculum — phage comes
from the dose. If per-patient phage-dose variability is wanted, wire IIV to the dose
amount (small follow-up).

## Done this session (2026-07-07) — immune module investigation

- **Investigated "immune module doesn't work" reports. The immune module is NOT
  broken** — innate immunity kills active bacteria (both strains) correctly in all
  three builder modes and the exact app solve path (BDF + extinction_threshold). Two
  real issues found:
  1. **Dormancy is an immune refuge (root cause of the user's 2-strain report).**
     Reproduced exactly (2 strains, dormancy on, defaults, immunity on → nadir 1.57e6,
     never clears, resistant fraction ~99.7% — matches user's numbers). Mechanism:
     immunity crushes the *active* resistant strain (B1: 4.7e8→2.9e6) but the resistant
     population mass-converts to dormancy when nutrients crash, and dormant/hibernating
     cells are **immune-privileged** — the engine neither kills them (`imm_kill_rate_D`
     defaults to 0/None) nor lets them stimulate immunity (excluded from
     `bac_total_active`). So a ~3.6e8 dormant reservoir persists and the infection never
     clears. Confirmed fix path: `imm_kill_rate_D > 0` drops the burden to 1.8e5. This
     is correct model biology (persister immune evasion), but the app gave no hint.
     **→ Added a UI warning** in the immunity tab (Tab 2) when dormancy + immunity are
     both on and `imm_kill_rate_D <= 0`, explaining the refuge and pointing to the
     control. No biology/default changes. Regression test
     `test_builder_modes.py::test_dormancy_creates_immune_refuge`. **51 tests passing.**
  2. **The `hill` immune module — FIXED this session.** Hill killing is
     `imm_max·B/(imm_kill50+B)` with `Imm` frozen, so it **ignores** `imm_kill_rate`,
     `imm_stim_rate`, `imm_decay_rate`, `initial_Imm` (verified: 1e7→1e13 on
     `imm_kill_rate` changes nothing). Previously `imm_max=1e7` did nothing at
     `imm_kill50=1e8`; the owner's `imm_kill50` 1e8→1e5 change (see below) now makes
     hill effective at the default `imm_max=1e7` (OFF→1e9, ON→0). Remaining defects
     were UX/docs: (a) the module-selector help described a stimulation formula the
     engine doesn't implement — corrected; (b) the immunity tab showed the inert
     innate-only fields in hill mode — now hill renders only `imm_max` + `imm_kill50`
     (+ `imm_kill_rate_D`) with a caption explaining `Imm` is frozen. Regression test
     `test_hill_immune_module_active_at_app_defaults` (hill controls the bloom at
     defaults; output is invariant to the inert fields).

- **Fixed blank "Antibiotics & Host Immunity" results graph.** When immunity was on
  but no antibiotic was configured, the plot did `ax2 = ... else plt.subplots(...)[1]`,
  drawing `Imm` onto a throwaway figure while `st.pyplot(fig)` showed the empty
  original — so the graph was blank. Now plots `Imm` on `ax1` (the displayed figure)
  in the no-antibiotic case. (User's reported symptom.)
- **Updated default innate immune parameters** (per owner): `imm_stim_rate` 1.0→**0.1**,
  `imm_decay_rate` 0.05→**0.1**, `imm_kill50` 1e8→**1e5**; `imm_kill_rate` (1e7) and
  `imm_stim50` (1e6) unchanged. Applied across session defaults, all three builder
  paths, UI widgets, preset-load, and repro-code.
- **Exposed `extinction_check_interval`** (pbisim solver arg) as a new "Extinction
  check interval (hours)" UI control (0 = check only at dose boundaries). Wired into
  the main sim, 1D/2D sweep solves, repro-code, and preset-load. Zeroes a sub-threshold
  strain at the chosen cadence so it can't regrow from a below-threshold pool (verified:
  a 3.4-cell residual regrows to 9.35e8 with interval off, stays extinct at interval=6h).

- **Fixed the `hill` immune-module UX/docs** (see #2 above): corrected the
  module-selector help text and hid the inert innate-only fields in hill mode.

**Still open (immune):** none outstanding.

## Done this session (2026-07-07) — stationary-phase pre-run fix

- **Fixed: long pre-run collapsed the treatment to a flat 0 CFU curve.**
  `run_sim_from_gui_params` (and the 1D/2D sweeps + repro-code) equilibrated with
  `stationary_phase_ic` but kept only `ic.B` and `ic.S`, **discarding the dormant
  reservoir `ic.D` and immune priming `ic.Imm`**. At stationary phase most of the
  culture is dormant, so the longer the pre-run the more cells were silently thrown
  away — the surviving active `B` shrinks with pre-run length until it drops below the
  extinction floor and the treatment plots as flat 0 (reproduced: with dormancy, active
  `ic.B` falls 2.6e8→2.0e5 as t_prerun 12→48 while `ic.D` grows to ~1e9; treatment peak
  collapses 2.6e8→2.0e5). Fix: carry the full stationary state — `initial_D = ic.D`,
  `initial_Imm = ic.Imm` — into the treatment `PBIModel`; treatment now starts at the
  full ~1e9 population for any pre-run length.
- Also: removed the bogus `S0=` kwarg passed to `stationary_phase_ic` (no such
  parameter — it was silently forwarded to scipy and ignored with a warning), and
  clamped the carried nutrient `initial_S = max(ic.S, 0)` (the pre-run can leave S
  slightly negative numerically). Applied in the main sim, 1D + 2D sweeps, and the
  auto-generated reproduction code. Regression test
  `test_prerun_carries_dormant_reservoir`. **52 tests passing.**

- **Fixed the same drop in the clinical-trial `PretreatmentPhase` path**
  (`trial_helper.create_model_factory`). The factory read `ic.B/P/S` from the
  per-patient config but took `initial_D`/`initial_Imm` from the GUI base kwargs, so a
  trial pretreatment discarded its dormant reservoir + immune priming. Now the factory
  prefers `config.initial_conditions.D`/`.Imm` (set by `PretreatmentPhase`) and clamps
  `initial_S = max(ic.S, 0)`. Verified: a 48 h pretreatment on a dormancy model now
  starts every patient's treatment at ~1e9 (was collapsing). Regression test
  `test_trial_pretreatment_carries_dormant_reservoir`. **53 tests passing.**

**Still open (pre-run):** none outstanding.

## Done this session (2026-07-08) — trial dosing, PK/PD outputs, builder consistency

Six owner-requested feature updates (all UI-verified via streamlit AppTest: every page
+ builder mode renders and a full trial runs with 0 exceptions):

1. **Antibiotic-aware default dose** (Environment & Dosing): dose-amount default is now
   target-specific — phage 1e8 PFU, antibiotic 10 mg, nutrient 1.0 — via
   `DOSE_AMOUNT_DEFAULTS`/`DOSE_AMOUNT_LABELS`; target selector moved before the amount
   (keyed per target). Applies to the single and repeat-regimen forms.
2. **`attenuation_rate` exposed** (phage config, all 3 builders): engine param
   "phage penetration decay with dormancy depth" — effective dormant adsorption =
   `adsorption_dormant × exp(−attenuation × depth)`. Per-phage input on the shared
   `phages` dict; wired into Direct (`with_phage_params`), BRG (`to_config`, broadcast),
   StrainSet (`StrainDefinition`), and all repro-code.
3. **Trial output dropdowns** (item 1c/1d): survival endpoint keeps tte + tt2lr;
   distribution-metric options are now Maximum Log Reduction, Log Reduction (baseline→last
   obs), Bacterial AUC, Nadir Count (time_to_clearance removed — it's a survival endpoint).
   New per-patient metrics `max_log_reduction` / `log_reduction_final` added to
   `trial_helper.trial_metric_fns()` and passed to `ClinicalTrial(metric_fns=...)`.
4. **Configurable n_latent + n_depth** (item 4): global **n_latent** control (Solver
   Settings → Model structure, `int_n_latent`) used by all three builders (replaced the
   hardcoded 5). **n_depth**: BRG gained a "Dormancy depth layers (Q)" control
   (`int_brg_n_depth`); StrainSet gained a per-strain "Depth layers (Q)" input (Direct
   already had one). All wired into build + repro.
5. **Trial dedicated dose editor + regimen** (item 1a): the Clinical Trials page now has
   its own "Trial Dosing Regimen" section — per-agent (phage/antibiotic) amount + start
   time + regimen (single, or repeat qX h × N) via `trial_helper.build_regimen_doses`.
   Arms (Control/Phage/Antibiotic/Combo) derive from it instead of the simulator's
   `int_doses`; the old init_P-inoculum injection is gone (the editor's t=0 dose is the
   inoculum). Base config still starts every arm at zero free phage.
6. **Trial raw PK/PD trajectories** (item 1b): per-arm median + IQR-band Plotly plots of
   total bacteria (CFU) and free phage (PFU) via `trial_helper.plot_pkpd_trajectories_plotly`,
   shown at the top of the trial outputs.

Regression tests in new `tests/test_trial_features.py` (regimen builder, distribution
metrics, PK/PD trajectory plots). **56 tests passing.**

**Note:** the Clinical Trials page no longer reads the simulator's Environment & Dosing
`int_doses` — trial doses come solely from the new Trial Dosing editor.

## Done this session (2026-07-08 cont.) — multi-arm trials + BRG fitness cost

- **Clinical trial now supports arbitrary named dose arms** (item follow-up: the fixed
  Control/Phage/Antibiotic/Combo checkboxes only gave one dose level). The Clinical
  Trials page has a **Treatment Arms builder**: add any number of named arms (e.g.
  "Low dose" / "High dose"), each with its own phage + antibiotic regimen (single or
  repeat qX h × N) via `render_regimen_config` / `arm_dose_events` (app.py) →
  `trial_helper.build_regimen_doses`. Arms stored in `st.session_state.trial_arms`; a
  separate "Include Control arm" toggle. Arm assembly builds one `TreatmentArm` per
  config with de-duplicated names; every arm still starts at zero free phage. Verified
  via AppTest (arms = Control / Low dose / High dose, 0 exceptions) and regression test
  `test_multiple_dose_arms_produce_distinct_outcomes` (higher dose → lower bacterial AUC,
  larger phage peak).
- **BRG fitness-cost default 0.0 → 0.05** (phage loci + antibiotics, append + preset-load
  defaults) with help text noting it drives the equilibrium IC. The engine already wired
  `fitness_cost` correctly (fc=0 → resistant neutral → resistant-dominated equilibrium;
  fc=0.05 → WT-dominated `[1e7, 20]`); the input existed but defaulted to 0, which
  produced the fully-resistant equilibrium the owner observed. **56 tests passing.**

**Note:** the old fixed-arm checkboxes (`run_control`/`run_phage`/`run_abx`/`run_combo`)
and the single Trial Dosing editor are gone — replaced by the Treatment Arms builder.

## Done this session (2026-07-10) — bugfixes + Tier-1 scenario library

Bugfixes (each committed separately, all verified via streamlit AppTest):
- **OD/debris crash (Direct mode):** the app passed debris params to
  `builder.build(**extra_kwargs)`, but `ModelBuilder.build()` takes no kwargs →
  `TypeError` whenever "Track Bacteriolytic Cell Debris" was enabled. Fixed to use
  `builder.with_od_debris(u, v, kdis, od_to_cfu_conversion_factor)`. BRG/StrainSet
  were fine (their `to_config` accepts `**extra_config_kwargs`).
- **Number-input precision:** Streamlit infers `"%.2f"` from the step, so 0.001 was
  rounded to 0.00 on entry. Global wrapper injects `format="%g"` for any float
  `st.number_input` without an explicit format (auto-skips ints / scientific).
- **BRG + StrainSet phage PK/advanced:** those two modes lacked the Advanced Phage
  Kinetics (`phage_decay_Km`) + Phage PK (`pk_mode`/Vc/k_elim/k_in/k_out/Vi/Km_elim)
  widgets Direct had (build paths already read them; UI was missing). Added widgets +
  `phage_decay_Km` forwarding for both.

**Tier-1 Scenario library (presets rework):** a "scenario" = the *entire* input
configuration (builder mode, strains/phages/antibiotics, pairwise adsorption, dosing,
nutrient, immune, debris, solver, prerun, trial arms/IIV). Implemented as a **session
snapshot** (`dump_state_to_scenario` captures all `int_*` + `ads_<s>_<p>` +
`direct_phg_res_rates` + `trial_*` keys; `load_scenario_to_state` clears widget keys
then restores) rather than an inverse of `load_preset_to_state` — so new params are
captured automatically and every builder mode is covered. UI in the **Presets &
Tutorials** page: Save current config, list/Load/Delete, and **Export/Import the whole
library as versioned JSON** (`schema_version=1`) — the portable "personal DB" that works
on the stateless deploy. Also fixed a **pre-existing navigation bug**: the keyed nav
radio overrode `current_page`, so Load buttons (scenario *and* tutorial) didn't switch
pages — now routed through a pending-`_nav_to` hop applied before the radio instantiates.
Round-trip verified (config + adsorption + builder mode restored, navigates, re-runs).
Tests: `tests/test_scenarios.py` (AppTest round-trip + JSON export).
- **Tutorial presets REMOVED** (owner's call): deleted `pbisim_app/presets.py` +
  `tests/test_presets.py`, dropped the tutorial-card browser, and renamed the page
  **"Presets & Tutorials" → "Scenarios"** (nav string + routing). A new inline
  `DEFAULT_SCENARIO` constant in `app.py` replaces `TUTORIALS[0]` as the fresh-session /
  Reset startup config, so the app no longer depends on the pbisim tutorials.
  **48 tests passing** (was 60; −12 removed preset tests).

**Design direction (agreed with owner):** presets are a two-tier model — **Tier 1 =
Scenarios (full-config snapshots, DONE)**; **Tier 2 = a composable Parts library**
(bacteria / phages / antibiotics). Phage kinetics (burst/latent/adsorption) are NOT
phage-intrinsic — they're phage×host pair properties — so Tier-2 phage parts will carry
a **reference-host tag** + soft "verify for this strain" flag, and map directly onto a
future `pbisim-fit` "fit → save part" pipeline. Persistence stays JSON export/import
(a real per-user DB needs auth+storage, deferred).

## Done this session (2026-07-13) — Tier-2 Parts library

- **Parts library** (composable building blocks) added to the **Library** page
  (renamed from "Scenarios"; now has two sections: 💾 Scenarios + 🧬 Parts).
  A part = one reusable entity (bacterium / phage / antibiotic) = its param dict +
  provenance (`source`: educated guess | literature | pbisim-fit | experimental) +
  annotation. Three tabs; each: **save a current entity as a part**, list with
  **Load** (appends to the shared `int_strains`/`int_phages`/`int_antibiotics`, so it
  composes across all pages, capped at the builder max) and **Delete**, plus
  **Export/Import the whole library as versioned JSON** (`PARTS_SCHEMA_VERSION=1`).
- **Host-tagged phages** (agreed design): phage kinetics (burst/latent/adsorption) are
  phage×host properties, not phage-intrinsic — so phage parts carry a `reference_host`
  (selected from current strains) and loading one whose host isn't among the current
  strains raises a soft "verify kinetics for this host" flash. This maps directly onto
  a future `pbisim-fit` "fit → save part" pipeline (fit output = host-tagged phage part).
- Helpers in `app.py`: `PART_CATEGORIES` / `PART_SOURCES`, `export_parts_json` /
  `import_parts_json`, `clear_entity_widgets` (pops entity widget keys so appended
  entities re-read from data), and a cross-rerun `_flash` mechanism displayed after the
  sidebar (used for part-load feedback + host warning + navigation).
- Tests: `tests/test_parts.py` (save/load/append + host-tag + JSON export). **50 passing.**
- Note: AppTest raises a spurious `KeyError('part_pick_bacteria')` if you do an *extra*
  `.run()` after a part-load navigation — a harness artifact (leaving a page with a keyed
  widget). The real app is fine (the load render itself is 0-exception); tests avoid the
  post-nav extra run.

**Design status:** both preset tiers now DONE — **Tier 1 Scenarios** + **Tier 2 Parts**.
Future: `pbisim-fit` → save-part pipeline; a real per-user DB (needs auth+storage) if the
JSON export/import backbone is outgrown.

## Done this session (2026-07-14) — mutation-graph fix, segfault hardening, Calibration (Phase A)

- **Direct-mode custom mutation network**: the 2^m-strains rule was a pbisim-app
  limitation (only used `with_mutations(phage_resistance_rates=...)`); added a
  strain→strain→rate graph editor (shared `int_transitions`) → builds the (n,n)
  mass-conserving matrix (`mutation_matrix_from_transitions`) → `with_mutations(
  mutation_rates=...)`. Any strain count now supports mutation.
- **SIGSEGV (exit 139) hardening** for Render: forced `matplotlib.use("Agg")` before
  pyplot (GUI backend from Streamlit's thread segfaults headless); Dockerfile
  `MPLBACKEND=Agg` + `OMP/OPENBLAS/MKL/NUMEXPR_NUM_THREADS=1` (OpenBLAS reads the host
  core count and over-spawns threads in the container → classic segfault); clinical-
  trial `n_jobs` default 4→1 (joblib loky forking); enabled `faulthandler`. Root cause
  was environmental (host placement / thread over-subscription), exposed by frequent
  auto-deploys — the thread cap fixes it deterministically. Crashes stopped.
- **Calibration page (Phase A of pbisim-fit integration)** — `pbisim_app/fit_helper.py`:
  - **Observable registry** (`OBSERVABLES`): CFU / PFU / OD / **luminescence**, each
    declaring the model compartments it reflects + a **link** (None | ÷param | ×param).
    OD = biomass ÷ `od_to_cfu`; luminescence = *active* biomass (`B` only) × `rlu_per_cell`.
    Adding a signal = one entry (extensible, per owner's bioluminescence use case).
  - **Ingestion** (`normalize_fit_dataframe`): any CSV → canonical **pbisim-fit long
    format** (`time, arm, observable, value` + per-arm MOI conditions). Auto-detects the
    pbisim-fit format; a column-mapping UI handles Monolix (`ID,TIME,DV,MOI,PHAGE,EXPERI`,
    DV=OD, PHAGE×MOI=arm). Validated against the real `monophage_data`/`ck_data` (33 arms).
  - **Overlay**: arm **multiselect** → simulate the current model per arm (`initial_P =
    MOI × B0`) → overlay predicted vs observed + per-arm **RMSE** (log for CFU/PFU).
  - The app does NOT import pbisim_fit yet (kept out of the deploy); it mirrors the schema
    so the ingested dataset feeds pbisim-fit directly when wired.
  - Tests: `tests/test_fit_helper.py` (registry, links, Monolix ingestion, RMSE). **54 passing.**

**Design plan (agreed):** Phase A DONE. **Phase B** = manual parameter tuning (focused
sliders + the link factors `od_to_cfu`/`rlu_per_cell`, re-overlay, save tuned params as a
Part). **Phase C** = pbisim-fit hand-off (manual tune → `NLSRefiner` warm-start → NPE
posterior → host-tagged Part → posterior→IIV for trials). Manual tuning is kept as the
human-in-the-loop / warm-start front-end, not deleted. Luminescence fitting later needs
pbisim-fit to add the observable (cross-package coordination item).
