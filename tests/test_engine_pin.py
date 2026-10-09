"""The deployed image installs the engine packages at commits pinned in
deploy/*.ref (see Dockerfile).  Installing ``@main`` with a constant build arg let
Docker reuse a cached layer, so deploys silently kept a stale engine while the
app code updated (2026-10: the phage-free fix shipped without its engine half)."""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_refs_are_full_commit_shas():
    for name in ("pbisim.ref", "pbisim-fit.ref"):
        ref = (ROOT / "deploy" / name).read_text().strip()
        assert re.fullmatch(r"[0-9a-f]{40}", ref), f"deploy/{name} must be a 40-char SHA, got {ref!r}"


def test_dockerfile_installs_the_pinned_refs():
    df = (ROOT / "Dockerfile").read_text()
    assert "COPY deploy/pbisim.ref deploy/pbisim-fit.ref" in df
    assert "pbisim.git@${REF}" in df and "pbisim-fit.git@${FIT_REF}" in df
    code = "\n".join(l for l in df.splitlines() if not l.lstrip().startswith("#"))
    assert "@main" not in code and "PBISIM_REF=main" not in code


def test_help_reports_build_info():
    from pbisim_app.views.help import engine_build_info
    info = engine_build_info()
    assert set(info) == {"pbisim", "pbisim-fit"}
    assert info["pbisim"][0] is not None          # version resolvable
