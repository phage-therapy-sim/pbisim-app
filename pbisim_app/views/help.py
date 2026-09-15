"""Rendered by app.py when the Help page is selected.

A curated in-app orientation (quick start, what each page does, key concepts,
troubleshooting) plus the full bundled USER_GUIDE.md as the deep reference — read
from disk at runtime so it stays a single source of truth (no hand-duplication).
"""
from pathlib import Path

from pbisim_app.common import *  # noqa: F401,F403


# Repo root holds README.md / USER_GUIDE.md (one level above the package dir).
_REPO_ROOT = Path(__file__).resolve().parents[2]
_DOCS_URL = "https://phage-therapy-sim.github.io/pbisim-docs"
_API_REF = f"{_DOCS_URL}/API_REFERENCE.html"


def _read_doc(name):
    """Return the text of a repo-root markdown doc, or None if it can't be found."""
    p = _REPO_ROOT / name
    try:
        return p.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return None


def render():
    st.title("Help & User Guide")
    st.markdown(
        "<div class='info-banner'>New here? Start with <b>Quick start</b> below, then skim "
        "<b>What each page does</b>. The complete reference is in <b>Full user guide</b> at the "
        "bottom. Every simulation, sweep, and calibration runs locally — only the AI Assistant "
        "needs an API key.</div>",
        unsafe_allow_html=True,
    )

    # ── Quick start ───────────────────────────────────────────────────────────
    st.markdown("### Quick start")
    st.markdown(
        "1. **Build a model** on the **Interactive Simulator** page — pick a builder mode, set "
        "the bacteria / phage / antibiotic parameters, and configure the growth, death, dormancy, "
        "and dosing settings.\n"
        "2. **Run it** with the **Run Simulation** button → CFU / phage / OD trajectories and "
        "outcome metrics appear on the right.\n"
        "3. **Explore** — sweep a dose (**Dose-Response Sweeps**) or any parameter / signal "
        "function (**Parameter Sweeps**), or design a virtual trial (**Clinical Trials & "
        "Cohorts**).\n"
        "4. **Fit to your data** on the **Calibration** page — upload CFU / PFU / OD / luminescence, "
        "overlay the model, and run the non-linear least-squares fit.\n"
        "5. **Save & reuse** — freeze the current build as a **Model** (sidebar), or store a full "
        "config (**Scenario**) or a single organism (**Part**) in the **Library**.\n\n"
        "Stuck with a stale value after editing? Use **Reset Environment** (sidebar) for a clean "
        "session."
    )

    # ── What each page does ───────────────────────────────────────────────────
    st.markdown("### What each page does")
    _pages = [
        ("Interactive Simulator",
         "The core page. Three builder modes — **Direct** (ModelBuilder), **Binary Genotypes** "
         "(BRG, resistance loci), and **Custom Strains** (StrainSet, arbitrary strain/mutation "
         "graph). Set organism kinetics, growth/death/lysis/dormancy signal functions, immunity, "
         "OD/debris, nutrient environment, dosing, and solver settings; run to see trajectories + "
         "metrics. A **Show model config** toggle prints the fully-resolved config, and every run "
         "emits a standalone **reproduction script**."),
        ("Dose-Response Sweeps",
         "Sweep a phage or antibiotic **dose** across a series (log/linear, or MOI-scaled) and "
         "compare trajectories + summary metrics. A dose of 0 is a genuine no-treatment control."),
        ("Parameter Sweeps",
         "Sweep **any model parameter** in 1D or 2D (contour maps for 2D), or a **signal function** "
         "categorically (compare Monod vs logistic vs diauxic growth, dormancy signals, etc.). "
         "Coupled mode varies several linked parameters together."),
        ("Clinical Trials & Cohorts",
         "Virtual in-silico trials: define named treatment arms, add patient-to-patient variability "
         "(IIV), an optional stationary-phase pre-run, and run parallel cohorts → Kaplan–Meier "
         "curves, per-arm metric distributions, and CSV / NLME export."),
        ("Calibration",
         "Fit the model to your data. Upload a CSV of CFU / PFU / OD / luminescence, map the "
         "columns (auto-detect + Monolix/NONMEM support), filter / regroup / aggregate replicates, "
         "overlay the current model vs observations with live RMSE, then run the pbisim-fit NLS "
         "fit (role-based parameter table: Fixed / Free / Derived, MAP priors, per-arm covariates). "
         "A **Compare models (AIC / BIC)** panel ranks candidate models with a parsimony penalty, so "
         "a richer model must earn its extra parameters."),
        ("AI Assistant",
         "Describe a simulation in plain language; Claude writes and runs the pbisim code and "
         "explains the result (with a self-healing retry loop). Needs an Anthropic API key "
         "(sidebar or `ANTHROPIC_API_KEY`)."),
        ("Library",
         "Two stores: **Scenarios** (whole-config snapshots) and **Parts** (reusable bacteria / "
         "phages / antibiotics, phages host-tagged). Export/import either as versioned JSON — your "
         "portable personal database."),
    ]
    if scripting_enabled():
        _pages.append(
            ("Scripting",
             "A notebook-style Python scratchpad with a shared kernel (opt-in power-user page). "
             "Cells use a real code editor (Tab-indent, syntax highlighting, line numbers) — run a "
             "cell with its ▶ Run button or Ctrl/Cmd+Enter. The execution sandbox is "
             "**research-grade, not a security boundary** — only for trusted local/authenticated use."))
    for _name, _desc in _pages:
        with st.expander(_name):
            st.markdown(_desc)

    st.markdown(
        "**Models vs Scenarios vs Parts.** A **Model** (sidebar) freezes the organism/kinetics so "
        "downstream tasks (sweeps, trials, fitting) run against a fixed config instead of your live "
        "edits. A **Scenario** is a full snapshot of *everything* (dosing, solver, trial setup). A "
        "**Part** is one reusable entity (a bacterium, phage, or antibiotic)."
    )

    # ── Key concepts ──────────────────────────────────────────────────────────
    st.markdown("### Key concepts")
    with st.expander("CFU vs total load (what the plots count)"):
        st.markdown(
            "**CFU (colony-forming units) = active `B` + dormant `D`.** Infected (`I`) and "
            "hibernating (`H`) cells lyse and don't form colonies, so they're excluded from the "
            "plate count. The plots let you also show **total live load** (`B+D+I+H`) or "
            "**active only** (`B`). **OD** includes scattering from infected cells and (optionally) "
            "bacterial debris, and — if enabled — a nutrient-dependent cell-size factor.")
    with st.expander("MOI and the phage inoculum"):
        st.markdown(
            "Multiplicity of infection (MOI) = phage : bacteria ratio at *t = 0*. In sweeps and "
            "calibration you can specify the phage dose as an absolute titre (PFU) or as an MOI "
            "that scales with the bacterial inoculum. In clinical trials the phage is delivered "
            "as a dose event, so control / antibiotic-only arms start at genuinely zero phage.")
    with st.expander("Signal functions (growth / death / lysis / dormancy)"):
        st.markdown(
            "Several processes are modulated by a selectable **signal function** (set model-wide in "
            "the topmost builder panel):\n\n"
            "- **Growth** — constant, Monod (nutrient), logistic (density), Monod×logistic, "
            "density-throttled, Gompertz, **sequential / diauxic** (one nutrient pool consumed "
            "through ordered Monod phases), or **smooth two-efficiency Monod** (a differentiable "
            "diauxie — the Monod K blends from efficient at high nutrient to inefficient at low "
            "nutrient; the better-conditioned choice for fitting growth curves).\n"
            "- **Death** — constant, nutrient (starvation), density (crowding), or nutrient+density.\n"
            "- **Lysis progression** — constant, or nutrient-coupled (`frac_lysis`).\n"
            "- **Dormancy / resuscitation / depth-diffusion** — constant, nutrient, density, or "
            "nutrient+density.\n\n"
            "You can also **sweep across signal functions categorically** on the Parameter Sweeps "
            "page to compare them side by side.")
    with st.expander("Resistance & cross-resistance"):
        st.markdown(
            "Use **Binary Genotypes (BRG)** for a small number of resistance loci (phage and/or "
            "antibiotic) with cross-resistance / collateral-sensitivity, or **Custom Strains** for "
            "an arbitrary strain + mutation graph. Never start a resistant strain at exactly 0 — "
            "mutation flux vanishes when growth stops, so seed it at "
            "`max(mutation_rate × B₀, ~10)`. Growth-**independent** switching (persisters, "
            "phase variation) is a *phenotypic transition* (h⁻¹, applied to the population "
            "directly) — set it, and phage mutation / transition graphs, in the "
            "**Phenotypic transitions & phage evolution** expander below the builder columns.")

    # ── Troubleshooting ───────────────────────────────────────────────────────
    st.markdown("### Troubleshooting")
    with st.expander("A result looks wrong — shallow killing, no clearance, or a flat 0 curve"):
        st.markdown(
            "- Set **Solver Method = BDF** (Solver Settings) — LSODA can miss stiff phage+dormancy "
            "dynamics.\n"
            "- Check the phage **adsorption rate** to the target strain is ≥ 10⁻⁹ mL·h⁻¹ (0 = the "
            "phage can't infect that strain).\n"
            "- Check the **Extinction Threshold** isn't so high it zeroes strains prematurely.\n"
            "- Confirm each dose has the right **target** (phage vs antibiotic) and index.\n"
            "- A long stationary-phase **pre-run with a death rate but no dormancy** decimates the "
            "inoculum before treatment — the app warns when the pre-run leaves < 10 % of the "
            "culture.")
    with st.expander("A parameter reverted to an old value after I edited it"):
        st.markdown(
            "Streamlit keeps session state across hot-reloads. Use **Reset Environment** (sidebar) "
            "or open a fresh tab for a clean state. Very small decimals: the inputs use `%g` "
            "formatting, so values like 1e-3 are preserved.")
    with st.expander("The AI Assistant keeps failing"):
        st.markdown(
            "- Verify the API key (sidebar → **Test API Key & List Models**).\n"
            "- Break complex requests into steps (build the model, then dosing, then resistance).\n"
            "- If the conversation seems confused, **Reset Environment** clears the history.")
    with st.expander("A clinical trial is slow"):
        st.markdown(
            "Start with 10–20 patients per arm, shorten the end time, and run on more CPU cores "
            "(the trial parallelises across arms/patients).")

    # ── Full user guide (single source of truth) ──────────────────────────────
    st.markdown("### Full user guide")
    _guide = _read_doc("USER_GUIDE.md")
    if _guide:
        with st.expander("Open the complete USER_GUIDE.md", expanded=False):
            st.markdown(_guide)
    else:
        st.caption("USER_GUIDE.md not found in this deployment — see the online docs below.")

    # ── Links & about ─────────────────────────────────────────────────────────
    st.markdown("### More documentation")
    st.markdown(
        f"- **pbisim engine** — [API Reference]({_API_REF}) and "
        f"[13 tutorials]({_DOCS_URL}) (compartment structure, resistance genetics, PK/PD).\n"
        "- **Calibration / fitting** — the Calibration page ingests NONMEM/Monolix-style long "
        "data; see the in-page help for column mapping.\n"
    )
    st.caption("pbisim-app 0.1.0 · engine pbisim 1.0 · estimation pbisim-fit 0.1")
