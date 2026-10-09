"""An antibiotic-only (phage-free) model must work end to end in the app.

Setting "Number of phages" to 0 used to crash the Simulator after a successful
solve: the results panel summed free phage with ``result.sum_prefixes("P")``,
which raised ``KeyError`` when the model had no phage states (fixed in pbisim:
a known-but-absent compartment now sums to zeros).  The dosing UI also offered a
"phage" target that the engine would silently ignore.
"""
from pathlib import Path

APP = str(Path(__file__).resolve().parents[1] / "pbisim_app" / "app.py")


def _phage_free_app():
    import matplotlib
    matplotlib.use("Agg")
    from streamlit.testing.v1 import AppTest

    at = AppTest.from_file(APP, default_timeout=300)
    at.run()
    [w for w in at.number_input if w.label == "Number of phages"][0].set_value(0).run()
    [w for w in at.number_input if w.label == "Number of antibiotics"][0].set_value(1).run()
    at.session_state["int_doses"] = [
        {"time": 0.0, "amount": 10.0, "target_type": "antibiotic", "target_idx": 0,
         "route": "bolus", "duration": 0.0}]
    at.run()
    return at


def test_antibiotic_only_simulation_runs():
    at = _phage_free_app()
    assert len(at.exception) == 0, at.exception
    [b for b in at.button if b.label == "Run Simulation"][0].click().run()
    assert len(at.exception) == 0, at.exception
    assert len(at.error) == 0, [e.value for e in at.error]
    # the dose target dropdown no longer offers a compartment the model lacks
    targets = [s for s in at.selectbox if s.label == "Target Compartment"]
    assert targets and all("phage" not in s.options for s in targets)
