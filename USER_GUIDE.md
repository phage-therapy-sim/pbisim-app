# pbisim-app User Guide

**Version:** 0.1.0  
**Engine:** pbisim 2.0.0 · **Estimation:** pbisim-fit 1.0.0  
**Last updated:** 2026-08-28

---

## Contents

1. [Introduction](#1-introduction)
2. [Getting Started](#2-getting-started)
3. [App Layout](#3-app-layout)
4. [Interactive Simulator](#4-interactive-simulator)
5. [Dose-Response Sweeps](#5-dose-response-sweeps)
6. [Parameter Sweeps](#6-parameter-sweeps)
7. [Clinical Trials & Cohorts](#7-clinical-trials--cohorts)
8. [Calibration (fit to data)](#8-calibration-fit-to-data)
9. [AI Assistant](#9-ai-assistant)
10. [Library (Scenarios & Parts)](#10-library-scenarios--parts)
11. [Key Parameter Reference](#11-key-parameter-reference)
12. [Common Workflows](#12-common-workflows)
13. [Troubleshooting](#13-troubleshooting)

---

## 1. Introduction

**pbisim-app** is an interactive web interface for the `pbisim` phage–bacteria ODE simulation
engine, with parameter estimation via `pbisim-fit`. It lets you:

- Simulate phage therapy, antibiotic treatment, and combination regimens without writing code.
- Explore dose-response relationships and sweep any biological parameter.
- Design virtual clinical trials with patient-level variability.
- **Calibrate the model against your own experimental data** (CFU / PFU / OD / luminescence)
  and fit parameters with pbisim-fit's non-linear least squares.
- Save configurations as reusable **Models**, **Scenarios**, and **Parts**.
- Ask an AI assistant to build simulations from plain-language descriptions.

The app is aimed at researchers who understand the biology of phage–bacteria interactions
but want a rapid-exploration tool. Familiarity with basic PK/PD concepts helps for the
clinical-trial and calibration features. No programming experience is required for most
workflows; the AI assistant handles code generation on demand.

For deep documentation of the underlying ODE model — compartment structure, resistance
genetics, PK/PD sub-models — see the pbisim
[API Reference](https://phage-therapy-sim.github.io/pbisim-docs/API_REFERENCE.html) and the
13 pbisim tutorials published at
[phage-therapy-sim.github.io/pbisim-docs](https://phage-therapy-sim.github.io/pbisim-docs).
This guide focuses on navigating the app itself.

---

## 2. Getting Started

### Prerequisites

- Python ≥ 3.10 with the `pbisim` conda environment active (has `pbisim`, `pbisim-fit`, and
  this package installed). See the [README](README.md) for setup.
- The Anthropic API key is **only** needed for the AI Assistant page; all simulation, sweep,
  and calibration features run entirely offline.

### Launching the server (browser)

```bash
conda activate pbisim
cd /path/to/pbisim-app
python -m streamlit run pbisim_app/app.py     # or, once installed:  pbisim-app
```

The terminal prints the local URL (default **http://localhost:8501**). Open it in a browser.

### Desktop window (optional)

For a native-window experience instead of a browser tab, install the desktop extra and run
the desktop launcher:

```bash
pip install -e '.[desktop]'
pbisim-app-desktop
```

It runs the same local server inside a native OS window. If no native webview engine is
available it opens your browser instead; `PBISIM_APP_BROWSER=1` forces the browser. (Windows
and macOS work out of the box; Linux needs a current system webview library.)

### Hosted version

A deployed instance runs on Render and may be gated by a shared password — no local install
needed, but the AI code sandbox is research-grade, so treat the hosted instance as trusted-use
only.

### Stopping the server

Press **Ctrl+C** in the terminal where Streamlit is running (or close the desktop window). If
the process persists: `pkill -f "streamlit run"`.

### Session state and cache

Each browser tab has its own independent session. Closing a tab and reopening it starts
fresh. If the app feels stale after a code update, the fastest reset is to open a new
tab — or use the **Reset Environment** button in the sidebar (see §3).

---

## 3. App Layout

The sidebar on the left holds the navigation, the active-Model selector, and the AI settings.

### Navigation radio buttons

Select one of seven pages:

| Page | Purpose |
|---|---|
| **Interactive Simulator** | Build and run a single simulation manually |
| **Calibration** | Fit the model to experimental data (pbisim-fit NLS) |
| **Dose-Response Sweeps** | Sweep dose range for one or two agents |
| **Parameter Sweeps** | Sweep any model parameter (1D or 2D) |
| **Clinical Trials & Cohorts** | Virtual patient cohorts, KM curves, NLME export |
| **AI Assistant** | Natural-language simulation builder |
| **Library** | Save/load full-config Scenarios and composable Parts |

### Models — reusable frozen configurations

The **Interactive Simulator** edits a *live builder draft*. You can **freeze** the current
organism/kinetics/environment as a named **Model** (sidebar), and every task page (sweeps,
trials, calibration) has a **Model selector** to run against either the live draft or a chosen
frozen Model — so builder edits don't silently change a task you set up earlier. A Model
captures organism/kinetics/environment/solver only, not dosing or trial design. On the
Simulator and each task page, a **📋 Show model config** toggle reveals a full sectioned
snapshot of the selected model (growth & nutrient environment, dormancy, phage, immunity,
OD/debris, initial conditions, solver) — no tab-hunting.

### AI Settings (sidebar expander)

- **Anthropic API Key** — masked; needed only for the AI Assistant, never stored between
  sessions.
- **Model** — which Claude model powers the AI. With a valid key the list is fetched from your
  account (`/v1/models`); otherwise a curated default is shown.
- **Test API Key & List Models** — verifies the key and lists authorized models.
- **Show generated code / assumptions** toggles live on the AI Assistant page itself.

### Reset Environment button

At the bottom of the sidebar. Clears session state and reloads the default configuration —
a fresh start without restarting the server.

---

## 4. Interactive Simulator

> **Screenshot placeholder** — full simulator page with tabs visible.

The Interactive Simulator is organized into four tabs:

| Tab | Content |
|---|---|
| 🧫 Strains & Phages | Bacterial and phage strain parameters |
| 🧪 Antibiotics & Immunity | Antibiotic PK/PD and immune response |
| 📅 Environment & Dosing | Nutrients, dose schedule, repeat-dose builder |
| ⚙️ Solver Settings | ODE solver, time span, extinction threshold |

Configure each tab, then scroll to the bottom of the page and click **Run Simulation**.

---

### 4.1 Builder Mode

At the top of the **Strains & Phages** tab, choose a *Builder Mode*:

| Mode | When to use |
|---|---|
| **Direct (ModelBuilder)** | Most simulations — any number of bacteria/phage strains with explicit per-strain parameters |
| **Binary Genotypes (BRG)** | Systematic resistance evolution — enumerate all 2ⁿ combinations of n resistance loci automatically |
| **Custom Strains & Graph (StrainSet)** | Named strains with custom mutation graphs — full control over topology |

Switching modes resets all strain/phage configuration and clears any cached simulation result.

---

### 4.2 Direct Mode

**Bacterial Strains**

Click **+ Add Strain** (or change the "Number of strains" counter) to add strains.
Parameters per strain:

| Parameter | Meaning | Typical range |
|---|---|---|
| Name | Display label | — |
| Initial density (B0) | Starting CFU/mL | 10⁶–10⁸ |
| Growth rate (μ) | Maximum growth rate (h⁻¹) | 0.5–2.0 |
| Bacteria-to-resource ratio | CFU per unit resource | 10⁹ |
| Fitness cost | Fractional reduction in μ for resistant strains | 0–0.3 |
| Dormancy enabled | Adds dormant compartment (D) | toggle |
| Q-depth | Number of dormancy layers | 1–5 |
| Initial dormant density (D0) | Starting dormant CFU/mL | 0 (default) |

> **Multi-strain tip:** The first strain (index 0) is treated as wild-type. Additional
> strains default to B0 = 0 (seeded by mutation). Only change B0 for strains that should
> be pre-existing in the inoculum.

**Phage Strains**

| Parameter | Meaning | Typical range |
|---|---|---|
| Name | Display label | — |
| Initial density (P0) | Starting PFU/mL | 10⁶–10⁸ |
| Burst size (Y) | Progeny phage per lysis event | 10–200 |
| Latent period (h) | Time from adsorption to lysis | 0.3–1.0 |
| Phage decay rate (m) | First-order decay h⁻¹ | 0.05–0.3 |
| Adsorption Rates | Per-strain matrix (mL h⁻¹) — one value per bacterial strain. WT default 1×10⁻⁸, resistant strains default 0. | 10⁻¹⁰–10⁻⁷ |
| Dormant adsorption | Adsorption rate to dormant compartment | 0 (usually) |
| PK mode | None / 1-compartment / 2-compartment | toggle |

> **Adsorption rates are per-strain:** each row is a bacterial strain, each column a phage.
> A resistant mutant should have a much lower adsorption rate (e.g. 10⁻¹¹) or zero. The
> phage cannot kill what it cannot adsorb.

---

### 4.3 Binary Genotypes (BRG) Mode

BRG mode enumerates all 2ⁿ genotype combinations across *n* resistance loci (phage and/or
antibiotic loci). The wild-type strain is the all-zero genotype; each locus adds a
resistance allele.

**Configuration panels:**

- **Bacterial base strain** — growth rate, fitness cost per locus, mutation rates, dormancy.
- **Phage definitions** — one panel per phage; set `adsorption_s` (to sensitive bacteria)
  and `adsorption_r` (to resistant bacteria at this phage's locus).
- **Antibiotic definitions** — EC50, Emax for sensitive and resistant genotypes.
- **Initial conditions** — either set per-genotype B0 manually, or enable
  **Equilibrium IC** to compute the mutation-selection balance distribution automatically
  (recommended for multi-locus models; total bacteria is configurable).

> The BRG generates `2ⁿ` strains internally. With 4 loci you get 16 strains; with 5 you
> get 32. Simulation time grows accordingly.

---

### 4.4 Custom Strains (StrainSet) Mode

Define named strains and specify their mutation graph explicitly. Each strain has
independent growth rates, adsorption rates, and sensitivity parameters. Mutation paths
are specified as directed edges (source strain → target strain, rate).

Use this mode when the resistance topology is not a simple binary lattice — for example,
sequential mutational pathways, or strains with heterogeneous fitness landscapes.

---

### 4.4b Phenotypic transitions & phage evolution (all modes)

Below the builder columns, the **Phenotypic transitions & phage evolution** expander
exposes the engine's three remaining *generator matrices* — the mutation matrix above is
the fourth. Each is a mass-conserving rate graph (`M[dest, origin]`, diagonal = −outflow)
edited as directed edges (from → to, rate):

| Graph | Engine field | Units | Coupled to |
|---|---|---|---|
| Bacterial mutation (per mode, above) | `mutation_rates` | per division | replication — vanishes when growth stops |
| **Bacterial phenotypic transitions** | `transition_rates` | h⁻¹ | nothing — a first-order switch on the population itself (persister / phase-variation switching; acts in stationary phase) |
| **Phage mutation** | `mutation_rates_phage` | fraction of each burst | lysis |
| **Phage phenotypic transitions** | `transition_rates_phage` | h⁻¹ | nothing — applied to free phage |

Nodes are the strain names (Direct / StrainSet), the auto-generated genotype labels
(BRG: `00`, `01`, … or `phi01_abx0`) and the phage names. The section appears once the
model has ≥ 2 strains or ≥ 2 phages. Configured edges are sweepable (Parameter Sweeps,
grouped under Bacterial / Phage), estimable on the Calibration page (only the edges you
configured are offered; the diagonal is always derived so the fit can't create or destroy
cells), captured by Scenarios, and mirrored in the reproduction script. The
state-dependent `with_transition_function` (a callable) is scripting-only.

---

### 4.5 Antibiotics & Immunity Tab

**Antibiotics**

Add one or more antibiotics. Each has:

| Parameter | Meaning |
|---|---|
| PK: Vc, k_elim, k12, k21 | Central volume, elimination rate, inter-compartment transfer (2-cpt) |
| PD: EC50, Emax, Hill | Sigmoidal kill curve parameters |
| Lytic fraction (f_lyse) | Fraction of antibiotic-killed cells that lyse (vs. die intact). Affects OD measurements and debris. |
| Inoculum effect constant | Shifts EC50 with log bacterial density (inoculum effect) |

**Immune Response**

Toggle **Enable adaptive immunity** to activate a T-cell-like effector compartment (Imm).
Parameters: stimulation rate, kill rate, kill50 (half-saturation), decay rate.

---

### 4.6 Environment & Dosing Tab

**Nutrients (Monod kinetics)**

Toggle **Track nutrients** to enable nutrient-limited growth (recommended for most
simulations). Parameters: Monod constant (K_s), nutrient inflow (s_in), washout rate
(s_out), initial nutrient (S₀), and **infected-cell nutrient consumption** (a multiple of
the uninfected per-capita uptake by which latent-infected cells draw down the substrate
while building phage; 0 = off).

> Leaving **Track nutrients** off gives unbounded exponential growth regardless of CFU —
> only appropriate for very short simulations or model-fitting contexts.

**Growth signal function** (set under *Strains & Phages → Growth model*)

Choose how the per-strain growth rate is modulated: constant, Monod (nutrient), logistic
(density), Monod×logistic, density-throttled, Gompertz (a nutrient-scale sigmoid), or two
approaches to *biphasic / diauxic* growth — **sequential (diauxic) Monod** (ordered phases
over one nutrient pool) and **smooth two-efficiency Monod** (a differentiable diauxie whose
Monod K blends from efficient at high nutrient to inefficient at low nutrient — the
better-conditioned choice when fitting growth curves). Any of these can be swept
categorically on the Parameter Sweeps page.

**Dose Schedule**

Add individual dose events:
- *Time* — when the dose is applied (hours post-inoculation)
- *Amount* — dose size (CFU/mL for bacteria, PFU/mL for phage, mg/L for antibiotics)
- *Target* — phage index, antibiotic index, or bacteria (rarely used directly)

**Repeat Dosing Regimen Builder**

For multi-dose schedules (e.g. "1×10⁸ PFU q12h × 3 doses starting at hour 2"):
1. Set start time, dose amount, interval, number of doses, and target.
2. Click **Append Repeat Doses** — the computed events are added to the dose table.

---

### 4.7 Solver Settings Tab

| Setting | Effect |
|---|---|
| Simulation end time (h) | How long to integrate |
| Time step dt (h) | Output resolution (not solver step size) |
| Solver method | BDF (default, stiff), LSODA (auto-stiff detection), Radau (stiff, slower) |
| Extinction threshold | Below this CFU/mL, a strain is zeroed out and treated as extinct. `0` = disabled. |

> Use **BDF** for most simulations involving phage + dormancy. **LSODA** is faster for
> simple models. **Radau** is sometimes more accurate for very stiff problems but slower.

---

### 4.8 Running the Simulation

Click **▶ Run Simulation** at the bottom of the page. Results appear below:

**Plot: CFU and PFU over time**

- **Solid lines** — bacterial populations (B + dormant D + infected I + hibernating H),
  i.e. total CFU/mL on a log₁₀ scale.
- **Dashed lines** — phage populations (PFU/mL).
- A horizontal grey line at the detection limit may be shown.

> CFU includes all viable bacteria: `B + D + I + H`. Dormant cells (D) are
> phage-invisible but can resuscitate; infected cells (I) are committed to lysis.

**Metrics table**

Below the plot, a summary table shows:
- Peak CFU and time of peak
- Minimum CFU and time of nadir
- Final CFU at t_end
- Time to clearance (first crossing below 1 CFU/mL, if any)

**Reproduction code**

Click **Show reproduction code** to expand a Python code block that exactly reproduces the
current simulation using the pbisim API. This is useful for:
- Archiving a configuration
- Running the simulation in a notebook
- Feeding the code to the AI assistant for further modification

---

### 4.9 Interpreting Results

| Observation | Likely cause |
|---|---|
| Bacteria immediately rebound after phage | Resistant mutants with B0 > 0 seeded at start, or phage extinct |
| No bacterial killing despite phage | Adsorption rate too low, or phage dose too small relative to bacteria |
| Simulation looks identical to no-phage run | Check solver method — BDF recommended; also check extinction threshold |
| Bacteria grow without bound | Nutrient tracking disabled; enable it or reduce simulation time |
| Dormant fraction persists after clearance | Dormant bacteria are phage-invisible; use a higher burst size or second dosing regimen |

---

## 5. Dose-Response Sweeps

> **Screenshot placeholder** — dose-response page showing log-spaced sweep controls and
> resulting trajectories.

Sweeps one or two agents over a range of doses and plots outcomes.

### Setting up a sweep

1. Configure a baseline simulation in the **Interactive Simulator** first (strains, phages,
   antibiotics, dosing schedule). The sweep will hold everything fixed except the swept
   doses.
2. Navigate to **Dose-Response Sweeps**.
3. For each agent to sweep, select: agent type (phage / antibiotic), index, dose spacing
   (Log or Linear), number of points, min/max dose.
4. Optionally enable **MOI scaling** — doses are expressed as MOI relative to initial
   bacterial density and converted automatically.
5. Click **Run Dose-Response Sweep**.

### Reading the output

- **Trajectories panel** — one time-series line per dose level, color-coded from low
  (light) to high (dark). Select which population to display (total CFU, individual strains,
  phage).
- **Metrics table** — tabulates peak CFU, nadir CFU, and clearance time for each dose level.

> If your sweep has mismatched vector lengths between strains and phages (e.g. one phage
> with two bacterial strains), a **padding warning** appears. The app pads missing values
> with zeros; verify this matches your intent.

---

## 6. Parameter Sweeps

> **Screenshot placeholder** — 1D sweep showing metrics vs. parameter value and 2D contour.

Sweeps any numeric field of the model configuration.

### 1D Sweep

1. Select **Parameter** from the dropdown (all `ModelConfig` fields are listed).
2. Set min, max, number of points, spacing (Log/Linear).
3. Select **Output metric** (e.g. minimum CFU, clearance time).
4. Click **Run 1D Sweep**. Output: a metric-vs-parameter line chart and a trajectory panel.

### 2D Sweep

1. Select two independent parameters.
2. Set ranges for each.
3. Click **Run 2D Sweep**. Output: a filled contour map of the metric in the 2D parameter space.

### n_depth sweeps (dormancy depth)

When sweeping **n_depth** (number of dormancy layers), the app automatically resizes the
model's state vector to match the new depth. This is handled transparently — no manual
adjustment needed.

---

## 7. Clinical Trials & Cohorts

> **Screenshot placeholder** — Kaplan-Meier curve panel with multiple treatment arms.

Simulates a cohort of virtual patients, each with individually drawn parameters from
log-normal or other distributions.

### Workflow overview

1. **Define patient variability (IIV)** — set coefficient-of-variation (CV%) for each
   parameter that varies across patients (e.g. bacterial inoculum, growth rate, phage
   adsorption). Parameters with CV = 0 are fixed across all patients.
2. **Configure treatment arms** — add one or more arms (e.g. "Phage alone", "Antibiotic
   alone", "Combination"). Each arm specifies a dose schedule applied on top of the shared
   IIV patient pool.
3. **Pre-treatment phase** — optionally run each patient's bacteria for a fixed period
   before treatment to reach a physiological steady state (e.g. 24 h stationary phase
   pre-inoculation).
4. **Set cohort size and endpoint** — number of patients per arm; endpoint: time-to-event
   (TTE), time to 2-log reduction (tt2lr), or custom.
5. Click **Run Clinical Trial**.

### Reading the output

- **Kaplan-Meier curves** — survival (non-clearance) probability over time per arm,
  with 95% confidence bands.
- **Metric distributions** — box plots of clearance times, nadir CFU, or other metrics
  per arm.
- **Export CSV** — full per-patient outcome table for offline analysis.
- **Export NLME** — longitudinal data in NONMEM/Monolix format for population PK/PD
  modelling (columns: ID, TIME, DV, DVID, EVID, MDV, AMT).

---

## 8. Calibration (fit to data)

The **Calibration** page fits the model to your own experimental data using `pbisim-fit`'s
non-linear least squares (`refine_nls`). It runs against the Model chosen in the page's Model
selector (live draft or a frozen Model).

### Workflow

1. **Upload & column mapping (§1).** Upload a CSV. Map the columns: **time**, **value** (the
   measurement / NONMEM `DV`), **observable** (fixed, or a column of `cfu` / `pfu` / `od` /
   `lum`), and one or more **arm-defining** columns. If the file interleaves NONMEM/Monolix
   **dose rows**, also map the **dose-event (EVID)** and **amount (AMT)** columns (see *Dose-row
   import* below).
2. **Filter · group · statistics (§2–3).** Filter rows, choose the grouping variables that
   define your curves/arms, and aggregate replicates (raw points, mean, or median + a
   percentile band).
3. **Overlay (§4).** Pick the arms and observables to overlay. The model is simulated once per
   arm and each observable is projected from that trajectory; you get per-arm **RMSE / R²** and
   a pooled objective. For **OD** data, set the single **od_to_cfu** factor here (see below).
4. **Per-arm conditions.** Each arm can carry its own initial density **B₀**, growth phase
   (**pre-run** = 0 for log phase, > 0 to equilibrate toward stationary phase), and phage dose.
5. **Manual tuning.** A **Show model builder** toggle reveals the *same* builder as the
   Interactive Simulator (all modes, every parameter, in **ratio** form for the inoculum) plus a
   compact global/structural block — so you can hand-tune and re-overlay before fitting.
6. **Fit parameters (§5).** In the parameter table set each parameter's **role** (see below),
   optionally add custom **θ** and priors, then **Run NLS fit**. The fit runs in a background
   thread; the fitted curves overlay automatically. **Apply** writes the fitted values back into
   the model, and you can **save the calibrated configuration** as a Scenario.
7. **Compare models (AIC / BIC).** Under the overlay, the **Compare models** panel snapshots the
   current model's pooled residuals plus its free-parameter count *k*, and ranks the candidates
   by AIC/BIC with ΔAIC/ΔBIC. Because a larger model always fits the training data at least as
   well (lower RSS), the information criteria add a per-parameter penalty so a richer model — say
   diauxic vs plain Monod growth — must *earn* its complexity. Only compare candidates overlaid on
   the **same** data (the panel warns otherwise). ΔAIC < 2 ≈ comparable support; > 10 ≈ decisively
   worse.

### Observables & the `od_to_cfu` factor

Supported observables: **CFU/mL**, **PFU/mL**, **OD**, and **luminescence (RLU)**. OD and
luminescence need a *link* to model state. For OD there is **one** factor, `od_to_cfu` (CFU per
OD unit), edited once in the §4 Overlay OD field. That single value drives everything: the
initial density (B₀ = first OD × od_to_cfu), the overlay (model OD = biomass ÷ od_to_cfu), the
debris `get_od()` if the OD/debris module is on, and the fit — so there are no competing OD
factors.

### Initial density B₀ (source modes)

B₀ follows pbisim-fit's **additive model** — a known inoculum is a *dose*, exactly like a phage
dose. The per-arm **B₀** selector offers:

- **First observation** — each arm starts at its earliest CFU (or OD × od_to_cfu). Simple, but
  carries measurement noise (a warning is shown).
- **Shared value** / **Per-arm values** — a *known* inoculum, recorded as a `bacteria` dose.
- **Estimate (shared / per-arm)** — B₀ is fit as an additive offset (`free_initial_conditions`).

If the dataset carries a `bacteria` dose row for an arm, that overrides the manual field.

### Fit-parameter roles

Each model parameter has a **role**: **Fixed** (held at its value), **Free** (estimated 1:1 with
its bounds/prior), or **Derived** (`= expression` of a custom θ). Blank bounds mean
unconstrained; **log** = log-space search; optional Gaussian **priors** (μ, σ) turn the fit into
MAP. Tie parameters together with the **Share** helper (sets rows to Derived = one θ). A
two-way **fit-spec text** panel round-trips the whole spec as a compact DSL
(`free … / fix … / theta … / map …`).

### Dose-row import (NONMEM/Monolix)

When you map an EVID column, dose rows (`EVID = 1`) are parsed — the *observable* column names
the target compartment (`bacteria` / `phage` / `antibiotic` / `nutrient`) and `AMT` its amount —
into per-arm doses that **gate the manual fields** (an arm with a data-provided dose hides its
manual B₀ / MOI) and feed the fit. Parsing is delegated to pbisim-fit's `EventTable`, so it
stays consistent with the fitter.

### Covariate effects (per-arm parameter links)

By default the fit builds **one shared model** for all arms (only the doses differ). The
**Covariate effects — per-arm parameter links (advanced)** panel (in the Fit-parameters section)
lets a parameter instead vary **per arm** by a covariate, NONMEM/Monolix-style:

- **θᵢ = θ_ref · (cov / ref)^β** (`power`), or `1 + β·(cov/ref − 1)` (`linear`), or
  `exp(β·(cov/ref − 1))` (`exponential`). The parameter's own row in the fit table sets **θ_ref**;
  the slope **β is estimated** as an extra fit parameter (one per link).
- A **covariate** is a **numeric grouping column** (e.g. temperature, inoculum) or **MOI** (the
  phage dose). Add the column as a grouping variable (§3) to make it available in the panel.
- Fill one row per link: *parameter · covariate · form · **ref** (the covariate value where the
  multiplier = 1) · β lower / upper / init*. The multiplier is exactly 1 at `cov = ref`, so θ_ref
  keeps its plain meaning. Control arms (no covariate) are left unmodulated.
- After the fit, the **fitted overlay applies each arm's link**, so every arm's model curve
  reflects its own covariate value.

*Example:* link **adsorption** to **MOI** (power form, ref = 1) to test whether apparent adsorption
scales with the multiplicity of infection — the fit returns the exponent β. The panel is
opt-in: an empty table leaves the fit as a single shared model.

---

## 9. AI Assistant

> **Screenshot placeholder** — chat interface with code block and plot output.

The AI Assistant translates natural-language requests into pbisim code, executes it, and
returns plots and narrative explanations.

### Getting started

1. Enter your Anthropic API key in the sidebar (see §3).
2. Select a Claude model (default: `claude-sonnet-4-6`).
3. Type your request in the chat box.

### Example prompts

```
Simulate 1 wild-type strain starting at 1e7 CFU/mL treated with a phage 
dose of 1e8 PFU/mL at t=0 and t=12. Show total CFU over 48 h.

Compare phage monotherapy vs phage + ciprofloxacin combination therapy. 
Ciprofloxacin: EC50=0.2, Emax=3, dosed at 1 mg/L at t=0.

Model the emergence of phage resistance with mutation rate 1e-7. 
Start with 1e7 WT bacteria and no resistant cells.
```

### How it works

1. Your message is sent to Claude along with the full pbisim API reference as system context.
2. Claude writes a Python code block using the pbisim API.
3. The app executes the code in a sandboxed namespace and displays the figures.
4. If the code raises an error, the app sends the traceback back to Claude automatically
   (up to 3 retries). On persistent failure, it rolls back to the pre-request conversation
   state.

### Tips for better results

- **Be specific about initial conditions** — mention starting CFU/mL and PFU/mL.
- **Mention doses explicitly** — time, amount, and whether it's phage or antibiotic.
- **Ask follow-up questions** — the conversation is multi-turn. After a plot appears, you
  can say "now add a second phage dose at t=24" and the AI will update the code.
- **Inspect the generated code** — enable "Show generated code" in the sidebar to see
  exactly what pbisim calls the AI made. You can copy this to a notebook for further work.
- **Use a saved configuration as a starting point** — load a Scenario from the Library (or set
  up a model in the Simulator), then switch to the AI Assistant and ask it to modify the loaded
  configuration.

### Limitations

- The AI generates code; the code runs in a lightweight sandbox that restricts file access
  and subprocess calls. It is not suitable for untrusted input.
- Complex multi-panel figures sometimes require prompting the AI to break them into steps.
- Very long simulation requests (e.g. large cohort trials) may time out in the chat
  interface; use the Clinical Trials page directly for those.

---

## 10. Library (Scenarios & Parts)

The **Library** page is your personal store of reusable configurations. It has two sections,
each with **Export/Import as versioned JSON** so a library is portable across machines (useful
since the hosted deploy is stateless).

### 💾 Scenarios — full-config snapshots

A **Scenario** captures the *entire* input configuration in one snapshot: builder mode,
strains/phages/antibiotics, pairwise adsorption, dosing, nutrient environment, immunity,
OD/debris, solver, pre-run, and clinical-trial arms/IIV. **Save current config** stores it;
**Load** restores the whole configuration (and navigates to the Simulator); **Delete** removes
it. Use scenarios to snapshot and reload complete experiments.

### 🧬 Parts — composable building blocks

A **Part** is a single reusable entity — a **bacterium**, **phage**, or **antibiotic** — with
its parameter dict plus provenance (educated guess / literature / pbisim-fit / experimental) and
an annotation. Save the current entity as a part; **Load** appends it to the shared strain/phage/
antibiotic lists so it composes across all pages. Because phage kinetics (burst / latent /
adsorption) are properties of a phage × *host* pair, not the phage alone, phage parts carry a
**reference-host tag**; loading one whose host isn't among the current strains raises a soft
"verify kinetics for this host" note. This maps directly onto a future pbisim-fit "fit → save
part" pipeline.

> Note: the older tutorial-preset browser was removed — the app no longer tracks the pbisim
> tutorials (which evolve independently). Fresh sessions start from a built-in default
> configuration; use Scenarios/Parts to build your own library.

---

## 11. Key Parameter Reference

### Bacterial strain

| Symbol | Parameter name | Unit | Notes |
|---|---|---|---|
| μ | `base_growth_rate` | h⁻¹ | Maximum (nutrient-saturated) growth rate |
| K_s | `monod_constant` | (resource units) | Half-saturation for Monod growth |
| φ | `dorm_frac` | — | Fraction of net growth entering dormancy |
| Q | `n_depth` | integer | Number of dormancy layers |
| γ | `resus_rates` | h⁻¹ | Dormant-to-active resuscitation rate (per layer) |
| μ_mut | `mutation_rate` | — | Per-division probability of resistance mutation |
| — | `transition_rates` | h⁻¹ | Growth-independent phenotypic switching between strains (§4.4b) |
| fc | `fitness_cost` | — | Fractional reduction in μ for each resistance allele |

### Phage strain

| Symbol | Parameter name | Unit | Notes |
|---|---|---|---|
| φ | `adsorption_s` / `adsorption_r` | mL h⁻¹ | Adsorption rate constant; s = sensitive bacteria, r = resistant |
| Y | `burst_size_s` | PFU/lysis | Progeny released per lysis |
| τ | `latent_period_s` | h | Intracellular development time |
| m | `decay_rate` | h⁻¹ | Free phage first-order decay |

### Antibiotic

| Symbol | Parameter name | Unit | Notes |
|---|---|---|---|
| E_max | `emax_s` / `emax_r` | × μ | Maximum kill rate as multiple of growth rate |
| EC50 | `ec50_s` / `ec50_r` | mg/L | Concentration at half-maximum kill |
| H | `hill` | — | Hill coefficient; > 1 = cooperative dose-response |
| f_l | `f_lyse` | — | Lytic fraction of antibiotic-killed cells (affects OD/debris) |

### PK parameters (1-compartment)

| Symbol | Parameter name | Unit | Notes |
|---|---|---|---|
| V_c | `Vc` | mL | Central volume of distribution |
| k_e | `k_elim` | h⁻¹ | First-order elimination rate |
| k_e0 | `ke0` | h⁻¹ | Effect-site equilibration rate (optional) |

---

## 12. Common Workflows

### Workflow A: First simulation — phage monotherapy

1. Press **Reset Environment** in the sidebar (reloads the default configuration).
2. Go to **Interactive Simulator → Strains & Phages**.
3. Confirm: 1 bacterial strain (B0 = 1×10⁷), 1 phage (P0 = 1×10⁶).
4. Go to **Environment & Dosing**: confirm a dose event at t=0 for 1×10⁶ phage.
5. Click **Run Simulation**.
6. Observe bacterial clearance around 24–36 h, phage peak then decline.

### Workflow B: Adding antibiotic to a failing phage therapy

1. Run Workflow A first to establish baseline.
2. Go to **Antibiotics & Immunity** tab, click **+ Add Antibiotic**.
3. Set EC50 = 0.2, Emax = 2.5, Hill = 1.5, k_elim = 0.3, Vc = 250.
4. Go to **Environment & Dosing**, add a dose event: time = 0, amount = 1.0 (mg/L),
   target = Antibiotic 0.
5. Click **Run Simulation**. Compare clearance time with and without the antibiotic.

### Workflow C: Modelling resistance emergence

1. Go to **Strains & Phages**, switch **Builder Mode** to **Binary Genotypes (BRG)**.
2. Set 1 phage locus, mutation rate = 1×10⁻⁷.
3. Leave phage adsorption_s = 1×10⁻⁸, adsorption_r = 1×10⁻¹¹.
4. Enable **Equilibrium IC** (recommended for BRG — sets realistic initial strain frequencies).
5. Dose phage at t=0 and t=12.
6. Click **Run Simulation**. The resistant strain (adsorption_r) should emerge and rebound
   after approximately 24–48 h.

### Workflow D: Systematic dose optimization (dose-response sweep)

1. Configure a baseline in the Interactive Simulator (Workflow A or B).
2. Navigate to **Dose-Response Sweeps**.
3. Set Agent = Phage 0, spacing = Log, 20 points, min = 10⁴, max = 10¹⁰ PFU/mL.
4. Click **Run Dose-Response Sweep**.
5. Identify the minimum dose achieving bacterial clearance from the metrics table.

### Workflow E: Virtual clinical trial

1. Configure a single treatment arm in the Interactive Simulator.
2. Navigate to **Clinical Trials & Cohorts**.
3. Set IIV: bacterial inoculum CV = 0.5 (50%), growth rate CV = 0.2.
4. Add Treatment Arm 1 (Phage) and Treatment Arm 2 (Combination).
5. Set cohort size = 50 patients per arm, endpoint = time-to-clearance.
6. Optionally enable **Pre-treatment phase** (24 h stationary pre-run).
7. Click **Run Clinical Trial**.
8. Inspect KM curves. Export CSV for survival analysis in R or Stata.

### Workflow F: Calibrate the model to experimental data

1. Build (or select) a Model whose structure matches your assay (e.g. Monod growth for a
   plateauing OD curve).
2. Navigate to **Calibration** and upload your CSV; map time / value / observable / arm
   columns (and EVID/AMT if the file has dose rows).
3. Group by the variables that define your curves; choose an aggregation (mean/median + band).
4. Select the arms and observables to overlay. For OD, set **od_to_cfu** in §4.
5. Pick a **B₀ source** (first-observation, a known shared/per-arm dose, or *estimate*).
6. Optionally **Show model builder** to hand-tune, then click **Overlay** and read the RMSE/R².
7. In the parameter table, set roles (**Free** the parameters to estimate; **Derived**/θ for
   shared or reparameterized ones), add priors if desired, then **Run NLS fit**.
8. The fitted curves overlay automatically. **Apply** the fit to the model and **save** it as a
   Scenario (or a host-tagged phage Part) for reuse in sweeps/trials.

---

## 13. Troubleshooting

### App won't start

- Confirm you are in the correct conda environment (`conda activate pbisim`).
- Check port 8501 is free: `ss -tlnp | grep 8501`. If occupied, kill the old process
  before restarting.

### App started but port is wrong

Streamlit auto-increments the port if 8501 is occupied (moves to 8502, etc.). Always read
the URL printed in the terminal on startup.

### Adsorption rate or other parameter shows an unexpected value after update

Streamlit preserves session state between hot-reloads. Use the **Reset Environment** button
or open a new browser tab to get a clean state.

### "Address already in use" when restarting

```bash
pkill -f "streamlit run"
sleep 2
python -m streamlit run pbisim_app/app.py
```

### Simulation result looks wrong (shallow killing, no clearance)

- Ensure **Solver Method = BDF** in the Solver Settings tab — LSODA can miss stiff dynamics
  in phage+dormancy models.
- Check that the phage's adsorption rate to the target bacteria is ≥ 10⁻⁹ mL h⁻¹. A value
  of 0 means the phage cannot infect that strain.
- Check the **Extinction Threshold** — if set too high, strains are zeroed prematurely.
- Confirm the dose schedule: make sure the phage dose has **target = phage** and the right
  phage index.

### AI Assistant code fails repeatedly

- Check your API key is valid (sidebar → "Test API Key & List Models").
- For complex simulations, break the request into smaller steps (build the model first,
  then add dosing, then add resistance).
- If the session history seems confused, click **Reset Environment** to clear history and
  start fresh.

### Clinical trial is very slow

- Reduce cohort size (start with 10–20 patients per arm for exploration).
- Shorten the simulation end time.
- The trial uses all available CPU cores by default (`n_jobs=-1`); run on a machine with
  more cores for large cohorts.

---

*For the underlying pbisim API documentation, see the
[pbisim API Reference](https://phage-therapy-sim.github.io/pbisim-docs/API_REFERENCE.html).*
