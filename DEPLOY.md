# Deploying pbisim-app on Render

`pbisim-app` is **public**, but it depends on the **private** `pbisim` engine. The
only real constraint when hosting is: *the build must authenticate to the private
`pbisim` repo, and that credential must live only in the host's secret store —
never in this public repo.*

This directory ships everything needed:

| File | Purpose |
|---|---|
| `Dockerfile` | Multi-stage build; installs private `pbisim` in a throwaway builder stage (token never reaches the final image), editable-installs the app so `prompts/` resolves. |
| `.dockerignore` | Keeps the build context small; blocks `.streamlit/secrets.toml`. |
| `render.yaml` | Render Blueprint (service + env var declarations; secrets are `sync: false`). |
| `.github/workflows/deploy.yml` | Runs the test suite, then triggers a Render deploy **only if tests pass**. |
| `deploy/pbisim-repo-notify.yml` | **Template** to copy into the private `pbisim` repo so engine pushes redeploy the app. |
| `DEPLOY_TROUBLESHOOTING.md` | **When a working deploy breaks:** build-error decision table, why a valid PAT can go missing, token rotation checklist, post-deploy verification. |

---

## 1. Create a read-only token for the private engine

GitHub → **Settings → Developer settings → Fine-grained tokens → Generate**:

- **Resource owner:** `phage-therapy-sim`
- **Repository access:** *Only select repositories* → `pbisim`
- **Permissions:** `Contents: Read-only`
- Copy the token (starts with `github_pat_…`). This is `PBISIM_TOKEN`.

> Use a fine-grained, single-repo, read-only token — not your account PAT. Rotate it
> if it is ever exposed. Never paste it into the repo, a Dockerfile, or a log.

---

## 2. Deploy to Render

### Option A — Blueprint (recommended)

1. Render Dashboard → **New → Blueprint** → connect this repo. Render reads `render.yaml`.
2. When prompted, set the `sync: false` values:
   - `PBISIM_TOKEN` = the token from step 1 (**secret**)
   - `PBISIM_FIT_TOKEN` = a READ-ONLY PAT for the **pbisim-fit** repo (powers the
     Calibration NLS fit). **Omit it** if `PBISIM_TOKEN` already has read access to
     pbisim-fit — the build falls back to it. Set it only for a separate per-repo token.
   - `ANTHROPIC_API_KEY` = leave **empty** for public use (users enter their own key);
     set it only for a private/gated deployment (see §5).
   - Engine versions are **not** env vars: they are pinned to commit SHAs in
     `deploy/pbisim.ref` / `deploy/pbisim-fit.ref` (see §4). A leftover `PBISIM_REF` /
     `PBISIM_FIT_REF` in the dashboard is ignored — delete it.
     (pbisim-fit installs its torch-free core only — the fitting import is lazy.)
3. Deploy. Render exposes these env vars to the Docker build as build ARGs, so the
   Dockerfile's `ARG PBISIM_TOKEN` is populated automatically.
   Confirm in the build logs that the `pip install git+https://***@github.com/...`
   step succeeds (the token is masked).

### Option B — Manual service

New → **Web Service** → this repo → **Runtime: Docker** → add the same env vars →
Health check path `/_stcore/health`.

---

## 3. Gate deploys on passing tests (recommended)

`render.yaml` sets `autoDeploy: false` so a push doesn't ship a broken build. Instead
CI runs the suite first:

1. Render service → **Settings → Deploy Hook** → copy the URL.
2. This repo → **Settings → Secrets and variables → Actions** → add:
   - `RENDER_DEPLOY_HOOK_URL` = that URL
   - `PBISIM_TOKEN` = the read-only token (CI needs it to install the engine for tests)
3. Now `.github/workflows/deploy.yml` runs on every push to `main`: install engine +
   app → `pytest` → **on green**, POST the deploy hook → Render rebuilds & deploys.

Prefer zero-CI simplicity? Set `autoDeploy: true` in `render.yaml` and delete the
workflow — every push deploys immediately (no test gate).

---

## 4. Keep the app current when EITHER repo changes

- **`pbisim-app` changes:** you push → CI tests → deploy hook → Render rebuilds. ✔ automatic.
- **`pbisim` (engine) changes:** Render watches this repo, *not* the engine, so an
  engine push must fan out a trigger. Install the companion workflow:
  1. Copy `deploy/pbisim-repo-notify.yml` into the **pbisim** repo at
     `.github/workflows/notify-pbisim-app.yml`.
  2. In the **pbisim** repo, add secret `PBISIM_APP_DISPATCH_TOKEN` — a fine-grained PAT
     with `Contents: Read` + `Actions: Read/Write` on `pbisim-app`.
  3. Now a push to `pbisim` fires a `repository_dispatch` here → the same test-gated
     deploy runs against the new engine commit (passed through as `client_payload.ref`).

**The deployed engine is pinned.** The Dockerfile installs `pbisim` and `pbisim-fit` at the
40-character commit SHAs in `deploy/pbisim.ref` and `deploy/pbisim-fit.ref`. The ref files
are `COPY`ed in right before the install steps, so changing a SHA rebuilds exactly those
layers, and an unchanged SHA reuses them. The image therefore provably contains the named
commits; the Help page footer shows the installed short hashes.

**To ship an engine change:** push the engine, write its SHA into `deploy/pbisim.ref`
(`git -C ../pbisim rev-parse HEAD > deploy/pbisim.ref`), and push the app. CI tests against
the same pinned refs.

> **Why not `main`?** Until 2026-10 the build installed `@main` with a constant build arg.
> Docker's cache key for that step never changed, so Render reused a cached engine layer
> indefinitely: the app code updated on every deploy while the engine stayed at whatever
> commit was current when the layer was first built. A fix that spanned both repos (the
> 0-phage crash) deployed "successfully" and still failed in production.

---

## 5. Security — decide BEFORE exposing it publicly

Two app-specific risks (both flagged in `CLAUDE.md`):

1. **The AI Assistant `exec()`s Claude-generated code.** The executor sandbox is
   *research-grade only* — not safe for untrusted public input. A public AI Assistant
   page is effectively remote code execution.
2. **The Anthropic API key.** A server-side key means every anonymous visitor spends
   your money. Ship with **no** `ANTHROPIC_API_KEY` so users supply their own in the
   sidebar, and/or gate access.

Recommended for a first deployment — **put the whole app behind auth**, which neutralizes
both at once:

- **Cloudflare Access** (free tier): add your Render URL as a Cloudflare app, allow only
  specific emails/domains. No app changes. Best option.
- **Quick password gate:** add a small `st.text_input(type="password")` check reading a
  value from Render env/`st.secrets` at the top of `app.py` — weak but fine for a private demo.
- **Harden later** for true public use: disable the AI Assistant page, or isolate the
  executor (per-session container / RestrictedPython / gVisor).

---

## 6. Build & run locally (sanity check)

```bash
# BuildKit keeps the token out of the final image; --secret avoids it in `docker history`
export PBISIM_TOKEN=github_pat_xxx
docker build \
  --build-arg PBISIM_TOKEN="$PBISIM_TOKEN" \
  # engine refs come from deploy/*.ref; for a one-off: --build-arg PBISIM_REF_OVERRIDE=<ref>
  -t pbisim-app .

docker run --rm -p 8501:8501 \
  -e PORT=8501 \
  -e ANTHROPIC_API_KEY="$ANTHROPIC_API_KEY" \
  pbisim-app
# open http://localhost:8501
```

Verify the final image has no token: `docker history pbisim-app` should not reveal it
(it lives only in the discarded builder stage).

---

## 7. Hygiene

- The workspace git remotes currently embed a PAT — never let that value reach a
  committed file, CI log, or the public repo. Deploy uses a *separate* fine-grained
  read-only token.
- Rotate `PBISIM_TOKEN` periodically and immediately if exposed; update it in Render and
  in this repo's Actions secrets.
