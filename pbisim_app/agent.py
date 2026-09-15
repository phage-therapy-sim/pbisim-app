"""
agent.py — Claude API integration for pbisim-app.

Translates natural-language simulation requests into pbisim Python code,
executes the code in a sandboxed namespace, and returns the results.
"""

from __future__ import annotations

import os
import re
import textwrap
from pathlib import Path
from typing import NamedTuple

import anthropic

from pbisim_app.knowledge import select_cards


# ── Load system prompt ────────────────────────────────────────────────────────
_PROMPT_PATH = Path(__file__).parent.parent / "prompts" / "system_prompt.md"
_SYSTEM_PROMPT = _PROMPT_PATH.read_text(encoding="utf-8")


def _system_blocks(user_message: str, *, extra: list | None = None) -> list:
    """Build the ``system`` blocks for a request.

    The base prompt (API contract + general reasoning) is a single **cached** block so it
    is re-used across turns and self-healing retries. Keyword-gated **knowledge cards** for
    this specific query are appended as a small, uncached tail (:func:`knowledge.select_cards`),
    so the always-cached prefix stays roughly constant-size as the knowledge library grows.
    """
    blocks = [{"type": "text", "text": _SYSTEM_PROMPT, "cache_control": {"type": "ephemeral"}}]
    blocks += extra or []
    cards = select_cards(user_message or "")
    if cards:
        blocks.append({"type": "text", "text": cards})
    return blocks

# ── Claude model ──────────────────────────────────────────────────────────────
# Default to the strongest code model — one-shot accuracy matters more here than
# per-call cost, since a wrong first guess triggers a full extra round-trip. The
# sidebar dropdown lets the user switch to a faster/cheaper model when desired.
_MODEL = "claude-opus-4-8"


class AgentResponse(NamedTuple):
    """Structured response from the simulation agent."""
    code: str            # generated Python code
    narrative: str       # plain-English interpretation
    assumptions: str     # bullet list of assumptions
    raw_text: str        # full model response (for debugging)


class AgentRun(NamedTuple):
    """Result of an agentic (tool-using) generation.

    Unlike AgentResponse (single blind generation), this is produced by the model
    running its own code via the ``run_pbisim_code`` tool and correcting itself within
    one turn. ``result`` is the ExecutionResult of the last code the model ran (holds
    the figures/stdout to display); ``tool_calls`` counts how many times it ran code.
    """
    narrative: str       # the model's final explanation (after the code worked)
    code: str            # the last code the model executed
    result: object       # executor.ExecutionResult of that code (figures, stdout, error)
    success: bool        # did the final executed code run cleanly
    tool_calls: int      # number of run_pbisim_code (code) executions (1 = ran code once)
    transcript: tuple = ()  # per-run ({code, success, error, figures}); [0] = the FIRST execution
    configured: bool = False  # the model called configure_simulator (populated the widgets)


# Tool the model uses to run and iterate on its code inside a single turn.
_RUN_TOOL = {
    "name": "run_pbisim_code",
    "description": (
        "Execute Python code in the pbisim sandbox and return its stdout, the number of "
        "matplotlib figures it created, and the full traceback if it raised. Use this to "
        "TEST and DEBUG your simulation code before answering: call it, read the result, "
        "and if it errored or produced no figure, fix the code and call it again. The "
        "sandbox already has the full pbisim API, numpy (np) and matplotlib (plt) loaded. "
        "When the code runs cleanly and produces the requested plot, stop calling the tool "
        "and give your final plain-English explanation."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "code": {"type": "string", "description": "Complete, self-contained pbisim Python code."}
        },
        "required": ["code"],
    },
}

# Built-in tool: look up the REAL signature + docstring of a pbisim symbol from the
# installed package (ground truth, drift-proof) so the model doesn't guess.
_LOOKUP_TOOL = {
    "name": "pbisim_api_lookup",
    "description": (
        "Look up the REAL signature and docstring of a pbisim API symbol from the installed "
        "package (always current, never guessed). Call this BEFORE writing code whenever you "
        "are unsure of an exact signature, argument name, or array shape — e.g. "
        "'ModelBuilder.with_phage_params', 'PBIModel', 'stationary_phase_ic', "
        "'BinaryResistanceGenotypes.from_strains', 'StrainSet.to_config'. Pass a dotted name."
    ),
    "input_schema": {
        "type": "object",
        "properties": {"name": {"type": "string",
                                "description": "dotted pbisim name, e.g. ModelBuilder.with_phage_params"}},
        "required": ["name"],
    },
}


def _pbisim_api_lookup(name: str) -> str:
    """Resolve a dotted pbisim name and return its real signature + docstring."""
    import inspect
    import pbisim

    name = (name or "").strip()
    if name.startswith("pbisim."):
        name = name[len("pbisim."):]
    if not name:
        return "Give a pbisim name, e.g. 'ModelBuilder.with_phage_params' or 'stationary_phase_ic'."

    obj = pbisim
    try:
        for part in name.split("."):
            obj = getattr(obj, part)
    except AttributeError:
        return (f"'{name}' is not in the pbisim public API. Check the name, or run "
                "`dir(pbisim)` / `inspect.signature(...)` via run_pbisim_code to explore.")

    out = [f"{name}"]
    try:
        out[0] += str(inspect.signature(obj))
    except (TypeError, ValueError):
        out[0] += f"  (type: {type(obj).__name__})"
    doc = inspect.getdoc(obj)
    if doc:
        out.append(doc[:2000])
    return "\n".join(out)


# Tool the model uses to populate the app's Interactive Simulator (Direct mode) instead
# of writing throwaway code. Its schema is a bounded subset of the app config; anything it
# can't express is a signal to fall back to run_pbisim_code.
_CONFIGURE_TOOL = {
    "name": "configure_simulator",
    "description": (
        "Populate the app's Interactive Simulator from a structured configuration so the "
        "user sees the widgets filled in and can tweak & run it. Covers ALL THREE builder "
        "modes — set builder_mode:\n"
        "- 'direct' (default): an explicit list of strains + phages + antibiotics + doses. "
        "adsorption_rates is per-strain (length = #strains; 0 = resistant).\n"
        "- 'strainset': explicitly named strains with a custom mutation_graph "
        "(from→to→rate edges between strain names).\n"
        "- 'brg': binary resistance genotypes — ONE base strain (strains[0]) plus phages "
        "and antibiotics acting as resistance loci; the 2^n genotypes are generated "
        "automatically. Give each phage adsorption_s (susceptible) / adsorption_r "
        "(resistant) and each antibiotic emax_s/emax_r/ec50_s/ec50_r; optionally set "
        "equilibrium_ic + total_bacteria.\n"
        "Use this whenever the request fits the app (any of the three modes). Only fall "
        "back to run_pbisim_code for parameter/dose SWEEPS, clinical-trial cohorts, or "
        "bespoke analysis/plots the Interactive Simulator can't express."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "builder_mode": {
                "type": "string", "enum": ["direct", "brg", "strainset"],
                "description": "which Interactive Simulator mode to set up (default 'direct')",
            },
            "strains": {
                "type": "array", "minItems": 1,
                "items": {"type": "object", "properties": {
                    "name": {"type": "string"},
                    "growth_rate": {"type": "number", "description": "h^-1, e.g. 1.2"},
                    "initial_B": {"type": "number", "description": "CFU/mL at t=0, e.g. 1e7"},
                    "death_rate_B": {"type": "number"},
                    "dormancy_enabled": {"type": "boolean"},
                    "dormancy_rate": {"type": "number"},
                    "resuscitation_rate": {"type": "number"},
                }, "required": ["name"]},
                "description": "strain list (direct/strainset); for 'brg', strains[0] is the base strain",
            },
            "phages": {
                "type": "array",
                "items": {"type": "object", "properties": {
                    "name": {"type": "string"},
                    "burst_sizes": {"type": "number", "description": "progeny per lysis, e.g. 50"},
                    "latent_periods": {"type": "number", "description": "hours, e.g. 0.5"},
                    "phage_decay_rates": {"type": "number", "description": "h^-1, e.g. 0.1"},
                    "initial_P": {"type": "number", "description": "PFU/mL at t=0"},
                    "adsorption_rates": {
                        "type": "array", "items": {"type": "number"},
                        "description": "DIRECT/STRAINSET: per-strain adsorption (length = #strains); 0 = resistant",
                    },
                    "adsorption_s": {"type": "number", "description": "BRG: adsorption of susceptible genotype"},
                    "adsorption_r": {"type": "number", "description": "BRG: adsorption of resistant genotype (often 0)"},
                    "fitness_cost": {"type": "number", "description": "BRG: fitness cost of resistance to this phage"},
                    "mu": {"type": "number", "description": "BRG: mutation rate to resistance for this phage"},
                }, "required": ["name"]},
            },
            "antibiotics": {
                "type": "array",
                "items": {"type": "object", "properties": {
                    "name": {"type": "string"},
                    "emax": {"type": "number"}, "ec50": {"type": "number"},
                    "k_elim": {"type": "number", "description": "elimination rate h^-1"},
                    "emax_r": {"type": "number", "description": "BRG: Emax on the resistant genotype"},
                    "ec50_r": {"type": "number", "description": "BRG: EC50 on the resistant genotype"},
                    "fitness_cost": {"type": "number", "description": "BRG: fitness cost of resistance"},
                    "mu": {"type": "number", "description": "BRG: mutation rate to resistance"},
                }, "required": ["name"]},
            },
            "mutation_graph": {
                "type": "array",
                "description": "STRAINSET: strain→strain mutation edges (names must match strains)",
                "items": {"type": "object", "properties": {
                    "from": {"type": "string"}, "to": {"type": "string"}, "rate": {"type": "number"},
                }, "required": ["from", "to", "rate"]},
            },
            "transition_graph": {
                "type": "array",
                "description": ("Bacterial PHENOTYPIC transitions (growth-independent switching, "
                                "h^-1, e.g. persister / phase variation) as from→to→rate edges. "
                                "Node names = strain names (direct / strainset) or BRG genotype "
                                "labels ('00','01',... / 'phi01_abx0')."),
                "items": {"type": "object", "properties": {
                    "from": {"type": "string"}, "to": {"type": "string"}, "rate": {"type": "number"},
                }, "required": ["from", "to", "rate"]},
            },
            "phage_mutation_graph": {
                "type": "array",
                "description": "Phage→phage mutation at lysis (fraction of each burst), edges between phage names",
                "items": {"type": "object", "properties": {
                    "from": {"type": "string"}, "to": {"type": "string"}, "rate": {"type": "number"},
                }, "required": ["from", "to", "rate"]},
            },
            "phage_transition_graph": {
                "type": "array",
                "description": "Phage→phage phenotypic transitions of free phage (h^-1), edges between phage names",
                "items": {"type": "object", "properties": {
                    "from": {"type": "string"}, "to": {"type": "string"}, "rate": {"type": "number"},
                }, "required": ["from", "to", "rate"]},
            },
            "equilibrium_ic": {"type": "boolean", "description": "BRG: seed genotype densities from the replicator equilibrium"},
            "total_bacteria": {"type": "number", "description": "BRG: total CFU/mL for the equilibrium IC"},
            "doses": {
                "type": "array",
                "items": {"type": "object", "properties": {
                    "time": {"type": "number"}, "amount": {"type": "number"},
                    "target_type": {"type": "string", "enum": ["phage", "antibiotic", "nutrient"]},
                    "target_idx": {"type": "integer"},
                    "route": {"type": "string", "enum": ["bolus", "infusion"]},
                }, "required": ["time", "amount", "target_type"]},
            },
            "t_end": {"type": "number", "description": "simulation end time (h), default 48"},
            "dt": {"type": "number"},
            "track_nutrients": {"type": "boolean"},
            "immunity_enabled": {"type": "boolean"},
            "debris_enabled": {"type": "boolean", "description": "track optical density via the debris ODE"},
        },
        "required": ["strains"],
    },
}


# Tool the model uses to read the user's current simulation results before interpreting.
_SUMMARY_TOOL = {
    "name": "get_simulation_summary",
    "description": (
        "Read a compact summary of the user's CURRENT simulation results (the run showing "
        "in the Interactive Simulator): outcome metrics such as initial/final/nadir "
        "bacteria, time to clearance, 2-log-reduction time, per-strain final densities and "
        "resistant fraction, free phage, immune level, and OD. Call this WHENEVER the user "
        "asks you to interpret or discuss their results ('why did resistance win?', 'did it "
        "clear?', 'what's the outcome?', 'explain this'). Returns text; no input needed."
    ),
    "input_schema": {"type": "object", "properties": {}},
}


_TOOL_INSTRUCTION = (
    "Choose the right response for the user's request:\n"
    "- INTERPRET / discuss the user's CURRENT results ('did it clear?', 'why did resistance "
    "win?', 'explain the outcome') → call get_simulation_summary to read the metrics, then "
    "explain in plain English. Do NOT re-run a simulation just to interpret existing results.\n"
    "- QUESTION / CHAT / EXPLAIN (general, not about the current run) → answer directly in "
    "plain text; call NO tool.\n"
    "- SET UP a model the app's Interactive Simulator supports → call configure_simulator. "
    "This covers ALL THREE builder modes: 'direct' (explicit strain list), 'strainset' "
    "(named strains + a custom mutation graph), and 'brg' (binary resistance genotypes from "
    "a base strain + phage/antibiotic loci). It populates the app's own widgets so the user "
    "can see, tweak, and run it (editable + reproducible). Prefer this over writing code "
    "whenever the request fits the simulator. After configuring, tell the user it's ready in "
    "the Interactive Simulator.\n"
    "- Something the simulator CANNOT express (parameter/dose SWEEPS, clinical-trial cohorts, "
    "comparing many configurations at once, or bespoke analysis/plots), OR the user "
    "explicitly asks to run it now → use run_pbisim_code: write "
    "and verify the code, fix any traceback and re-run (a few attempts). If you are unsure "
    "of any pbisim signature, argument name, or array shape, call pbisim_api_lookup FIRST to "
    "check the real API — this prevents shape/name mistakes. Your LAST "
    "run_pbisim_code call must be the COMPLETE final script that produces the requested plot "
    "(and any printed values) — do not follow a working script with a separate diagnostic-"
    "only run, or the plot will be lost. Then reply with a concise explanation; do not paste "
    "the code in your final message.\n"
    "When unsure whether configuring covers the request, prefer run_pbisim_code (it is fully "
    "general). If configure_simulator reports it cannot express the request, switch to "
    "run_pbisim_code."
)


class SimulationAgent:
    """
    Wraps the Claude API to translate natural-language queries into
    pbisim simulation code.

    Parameters
    ----------
    api_key : str or None
        Anthropic API key.  If None, reads from ``ANTHROPIC_API_KEY`` env var.
    model : str
        Claude model identifier.
    max_tokens : int
        Maximum tokens in the response.
    history : list
        Conversation history (for multi-turn refinement).
    """

    def __init__(
        self,
        api_key: str | None = None,
        model: str = _MODEL,
        max_tokens: int = 4096,
    ) -> None:
        self.client = anthropic.Anthropic(
            api_key=api_key or os.environ.get("ANTHROPIC_API_KEY")
        )
        self.model = model
        self.max_tokens = max_tokens
        self.history: list[dict] = []

    def ask(self, user_message: str) -> AgentResponse:
        """
        Send a simulation request and return structured response.

        Parameters
        ----------
        user_message : str
            Natural-language simulation request from the user.

        Returns
        -------
        AgentResponse
            Contains generated code, narrative, and assumptions.
        """
        self.history.append({"role": "user", "content": user_message})
        try:
            response = self.client.messages.create(
                model=self.model,
                max_tokens=self.max_tokens,
                # Cached API-reference block + any query-relevant knowledge cards.
                system=_system_blocks(user_message),
                messages=self.history,
            )
        except Exception:
            self.history.pop()  # Keep history clean
            raise

        raw_text = response.content[0].text
        self.history.append({"role": "assistant", "content": raw_text})
        # Expose token usage of the most recent call (used by the eval harness for
        # cost/latency reporting; harmless for the app).
        self.last_usage = getattr(response, "usage", None)

        return _parse_response(raw_text)

    def generate(self, user_message: str, execute, configure=None, summarize=None,
                 max_tool_calls: int = 6) -> AgentRun:
        """Agentic generation. The model decides what the request needs and, within one
        turn: answers directly (no tool), reads the current results via
        ``get_simulation_summary`` (``summarize`` handler) to interpret them, populates the
        app's Interactive Simulator via ``configure_simulator`` (``configure`` handler), or
        writes and runs code via ``run_pbisim_code`` (``execute``), self-correcting.

        ``execute`` maps ``code:str -> ExecutionResult``. ``configure`` maps a config dict
        -> summary (applies as a side effect). ``summarize`` maps ``{} -> results text``.
        Optional handlers, when None, simply aren't offered. Returns an :class:`AgentRun`.
        """
        self._trim_history()             # keep the conversation (and memory/cost) bounded
        _entry_len = len(self.history)   # roll back to here if this turn fails
        self.history.append({"role": "user", "content": user_message})
        system = _system_blocks(
            user_message, extra=[{"type": "text", "text": _TOOL_INSTRUCTION}])
        tools = ([_RUN_TOOL, _LOOKUP_TOOL]
                 + ([_CONFIGURE_TOOL] if configure is not None else [])
                 + ([_SUMMARY_TOOL] if summarize is not None else []))

        last_code, last_result, tool_calls, transcript = "", None, 0, []
        results_hist = []
        configured = False
        try:
            for _ in range(max_tool_calls):
                response = self.client.messages.create(
                    model=self.model, max_tokens=self.max_tokens,
                    system=system, tools=tools, messages=self.history,
                )
                self.last_usage = getattr(response, "usage", None)
                self.history.append({"role": "assistant", "content": response.content})

                tool_uses = [b for b in response.content if getattr(b, "type", None) == "tool_use"]
                if not tool_uses:
                    narrative = "".join(getattr(b, "text", "") for b in response.content
                                        if getattr(b, "type", None) == "text").strip()
                    display = _pick_display(last_result, results_hist)
                    ok = bool(display and display.success)
                    return AgentRun(narrative, last_code, display, ok, tool_calls,
                                    tuple(transcript), configured)

                tool_results = []
                for tu in tool_uses:
                    name = getattr(tu, "name", "run_pbisim_code")
                    if name == "configure_simulator" and configure is not None:
                        summary = configure(tu.input or {})
                        is_err = isinstance(summary, str) and summary.strip().upper().startswith("ERROR")
                        configured = configured or not is_err
                        tool_results.append({
                            "type": "tool_result", "tool_use_id": tu.id,
                            "content": summary, "is_error": is_err,
                        })
                    elif name == "get_simulation_summary" and summarize is not None:
                        tool_results.append({
                            "type": "tool_result", "tool_use_id": tu.id,
                            "content": summarize(tu.input or {}),
                        })
                    elif name == "pbisim_api_lookup":
                        tool_results.append({
                            "type": "tool_result", "tool_use_id": tu.id,
                            "content": _pbisim_api_lookup((tu.input or {}).get("name", "")),
                        })
                    else:  # run_pbisim_code
                        code = (tu.input or {}).get("code", "")
                        result = execute(code) if code else _no_code_execution_result()
                        last_code, last_result, tool_calls = code, result, tool_calls + 1
                        results_hist.append(result)
                        transcript.append({"code": code, "success": bool(result.success),
                                           "error": result.error, "figures": len(result.figures)})
                        tool_results.append({
                            "type": "tool_result", "tool_use_id": tu.id,
                            "content": _format_tool_result(result),
                            "is_error": not result.success,
                        })
                self.history.append({"role": "user", "content": tool_results})

            # Ran out of tool budget — return the best (last) attempt.
            display = _pick_display(last_result, results_hist)
            ok = bool(display and display.success)
            return AgentRun(
                "Reached the maximum number of code-execution attempts."
                + ("" if ok else " The last attempt still errored."),
                last_code, display, ok, tool_calls, tuple(transcript), configured,
            )
        except Exception:
            # Roll the whole turn back — leaving a partial tool exchange (an assistant
            # tool_use with no matching tool_result) would break the next API call.
            del self.history[_entry_len:]
            raise

    def _trim_history(self, max_turns: int = 16) -> None:
        """Keep the API conversation bounded so a very long chat doesn't grow cost/latency
        without limit (the tool loop appends code, tracebacks, and tool results every turn).
        The cap is generous — normal conversations keep full context; only very long ones
        drop their oldest turns. Trims whole turns from the front at user-prompt boundaries,
        preserving each kept turn's tool_use/tool_result pairing and the leading user msg."""
        starts = [i for i, m in enumerate(self.history)
                  if m.get("role") == "user" and isinstance(m.get("content"), str)]
        if len(starts) > max_turns:
            del self.history[:starts[len(starts) - max_turns]]

    def reset(self) -> None:
        """Clear conversation history (start a new simulation session)."""
        self.history.clear()


def _pick_display(last_result, history):
    """Which execution's figures/stdout to surface. Normally the last run, but if it
    produced no figure (e.g. the model ran a follow-up diagnostic that didn't re-plot),
    fall back to the most recent successful run that DID make a figure, so the user (and
    the eval) still see the plot."""
    if last_result is not None and getattr(last_result, "figures", None):
        return last_result
    for r in reversed(history):
        if r.success and r.figures:
            return r
    return last_result


def _no_code_execution_result():
    from pbisim_app.executor import ExecutionResult
    return ExecutionResult(success=False, figures=[], stdout="",
                           error="the tool call contained no code")


def _format_tool_result(result) -> str:
    """Compact textual feedback the model reads after running code."""
    parts = [f"figures_created: {len(result.figures)}"]
    out = (result.stdout or "").strip()
    if out:
        parts.append("stdout:\n" + out[:2000])
    if result.success:
        parts.append("status: SUCCESS")
    else:
        err = (result.error or "").strip()
        parts.append("status: ERROR (fix the code and run it again)\n" + err[-2500:])
    return "\n".join(parts)


# Any fenced block: captures an optional language tag and the body.
_FENCE_RE = re.compile(r"```[ \t]*([A-Za-z0-9_+.-]*)[ \t]*\r?\n(.*?)```", re.DOTALL)
_PBISIM_TOKENS = ("solve_ode", "ModelBuilder", "PBIModel", "BinaryResistanceGenotypes",
                  "StrainSet", "DoseSchedule", "import ", "plt.", "np.")


def _extract_code(text: str) -> str:
    """Pull the Python code out of a model response, robustly.

    Tolerates ```py / ```Python / bare ``` fences and multiple blocks: prefers
    Python-tagged (or untagged) blocks, and among those picks the one that most
    looks like a pbisim script (falling back to the longest). Returns "" if none.
    """
    blocks = [(lang.lower(), body) for lang, body in _FENCE_RE.findall(text)]
    if not blocks:
        return ""
    py = [b for lang, b in blocks if lang in ("python", "py", "python3", "")]
    candidates = py or [b for _, b in blocks]

    def score(body: str) -> tuple:
        hits = sum(tok in body for tok in _PBISIM_TOKENS)
        return (hits, len(body))

    return textwrap.dedent(max(candidates, key=score)).strip()


def _parse_response(text: str) -> AgentResponse:
    """Extract code block, narrative, and assumptions from model response."""
    code = _extract_code(text)

    # Remove ALL fenced blocks for narrative extraction (any language).
    text_no_code = _FENCE_RE.sub("", text).strip()

    # Extract assumptions bullet list (after "assumption" keyword)
    assumptions_match = re.search(
        r"(?i)(assumptions?.*?)(?:\n\n|\Z)", text_no_code, re.DOTALL
    )
    assumptions = assumptions_match.group(1).strip() if assumptions_match else ""

    # Narrative: everything that isn't the assumptions block
    narrative = text_no_code.replace(assumptions, "").strip()

    return AgentResponse(
        code=code,
        narrative=narrative,
        assumptions=assumptions,
        raw_text=text,
    )
