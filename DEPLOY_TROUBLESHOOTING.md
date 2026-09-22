# Deploy troubleshooting — Render build failures & token rotation

Companion to [`DEPLOY.md`](DEPLOY.md), which covers **first-time setup**. This file is for
when a deploy that used to work stops working. Work top-down: §0 tells you which section
you need.

The build depends on two credentials-gated steps that run **before** any app code is
touched, so most "build failed" emails are about those two lines in the `Dockerfile`, not
about a code change you just pushed.

---

## 0. Read the error string before changing anything

Render's email says only `exit status 1`. Open the service → **Logs** (or the failed
deploy → **Build logs**) and find the first red line. Then:

| Build-log line | What it means | What it is **not** |
|---|---|---|
| `ERROR: PBISIM_TOKEN build arg is required` | Render handed the build an **empty** value for `PBISIM_TOKEN`. | Not an expired/invalid token — see §1. |
| `Failed to connect to github.com port 443` | The build host could not open a TCP connection to GitHub. | Not an auth problem — see §2. |
| `Authentication failed` / `403` / `could not read Username for 'https://github.com'` | The token was **present but rejected**: expired, revoked, or lacking access to that repo. | Not a missing env var — see §3. |
| `no open ports detected` | Streamlit didn't bind `$PORT` inside Render's scan window. | Not a token problem — see §4. |

> **Do not mint a new PAT reflexively.** Three of those four rows are unaffected by
> minting a new token, and rotating one costs you four places to update (§6).

### Run this test first: retry the deploy with no changes

Render Dashboard → **Manual Deploy → Deploy latest commit**. It costs one build and it
partitions the problem cleanly:

| Retry result | What it proves |
|---|---|
| **Green** | The configuration is fine — `PBISIM_TOKEN` is set on the service and the PAT is valid, because the same guard would otherwise fire again. The red build hit a **transient build-host fault**; go to §1.0 and stop investigating tokens. |
| **Red, same error** | A real, persistent misconfiguration. Continue to the matching section below. |

Observed 2026-09-22: a build failed with the empty-token guard **and** a port-443 failure,
and the next deploy went green with no change to the repo, the service, or the token. Two
different failures of the build environment in one build, neither caused by anything in
this repo.

---

## 1. `ERROR: PBISIM_TOKEN build arg is required`

### What produces it

`Dockerfile`, the private-engine step:

```dockerfile
ARG PBISIM_TOKEN
RUN test -n "$PBISIM_TOKEN" || (echo "ERROR: PBISIM_TOKEN build arg is required" >&2; exit 1) \
    && pip install "git+https://${PBISIM_TOKEN}@github.com/phage-therapy-sim/pbisim.git@${PBISIM_REF}"
```

`test -n` tests for a **non-empty string**. It says nothing about whether the token is
valid. An expired, revoked, or wrong-scope PAT is still a non-empty string: it sails past
this guard and fails a few seconds later inside `git`/`pip` with one of the §3 messages
instead.

**So if you see this message, the token's expiry date is irrelevant.** Render supplied
nothing. Your PAT on GitHub may be perfectly healthy.

### Why it can appear suddenly when nothing changed

**0. A transient build-host fault — check this first.** The guard reports an empty
variable, but "empty" can also mean *the build ran without the service's environment
injected*, which is an infrastructure failure rather than a configuration one. The tell is
decisive: **a no-change retry succeeds.** If the same red build is accompanied by a §2
port-443 error, treat the two as one degraded build host and retry before changing
anything. This is what happened on 2026-09-22.

If a retry does **not** fix it, the variable really is missing, and these are the reasons —
ranked by how often each explains a *persistent* failure:

1. **The Docker layer cache was hiding it.** `ARG PBISIM_TOKEN` participates in the cache
   key of the `RUN` that follows. While that layer stayed cached, Render never executed
   the guard — so the variable could have been missing for weeks of green builds. Any
   cache-busting event (Render evicting the build cache, a build scheduled on a fresh
   builder, a plan/region change, a changed ARG value) makes that step run for the first
   time in a long while and it fails immediately.
   **Corollary worth acting on:** if that layer was cached, the running image's `pbisim`
   may be far older than `main`. Verify with §7.
2. **A Blueprint sync or a service re-creation.** `render.yaml` declares both tokens
   `sync: false`, which means *the Blueprint deliberately carries no value*. A service
   created — or re-created — from the Blueprint starts with those variables **empty**
   until someone types them into the dashboard. Re-syncing a Blueprint after editing
   `render.yaml` is the usual trigger.
3. **The variable is set somewhere the build can't see it.** Check whether it lives
   directly on the service or was moved into a shared **Environment Group**. Group and
   secret-file values are dependable at *runtime*; confirm in the build log that they are
   also injected as Docker **build args** before assuming they are.
4. **It is set on a different service.** Easy to do after a plan upgrade, a region move,
   or when a preview/staging service exists alongside the production one. Confirm you are
   looking at the service whose logs you just read.
5. **Key typo or a stray character.** `PBISIM_TOKEN ` (trailing space in the *key*) is a
   different variable. Note a trailing space in the *value* would still pass `test -n` and
   fail in §3 instead.

GitHub-side revocation is **not** on this list, because it cannot produce this message.
It is still worth knowing that a fine-grained PAT can die before its expiry date:

- GitHub secret-scanning auto-revokes a token it finds published anywhere.
- For an **org-owned** resource (`phage-therapy-sim` is an org), an org policy change —
  e.g. enforcing a maximum token lifetime, or requiring owner approval for fine-grained
  tokens — revokes or suspends existing tokens regardless of their expiry.
- An org owner can revoke any token granted against org repos.

Those all land you in §3, not here.

### Diagnose

1. Render Dashboard → the service → **Environment**. Is `PBISIM_TOKEN` listed, on **this**
   service, with a non-empty value? (Render masks it; "•••" means set.)
2. Check whether it is inherited from a linked **Environment Group** rather than set
   directly.
3. Scroll the failed **build log** to the step header — it names the Dockerfile line, which
   confirms you are in step 1 (engine) and not step 1b (pbisim-fit).

### Fix

- **Variable missing or empty** → re-enter it (§5 if you also need a new token), then
  redeploy with **Clear build cache** (§7). Do not skip clearing the cache: it is what
  makes the failing step actually re-run and what guarantees the image picks up the
  current engine.
- **Variable present and you believe the token is good** → you are probably in §3; verify
  the token with the command in §3 before rotating it.

---

## 2. `Failed to connect to github.com port 443`

This is a TCP-level failure — the build host could not reach GitHub at all. Authentication
never happened, so it is unrelated to tokens, scopes, and expiry.

1. **Retry the deploy.** Transient egress failures on the build host are the common cause.
   If this appears *alongside* the §1 empty-token error in the same build, that combination
   is strong evidence of a degraded build host rather than two separate bugs — retry before
   touching any configuration.
2. If it recurs across several attempts, it is Render-side infrastructure — check Render's
   status page, then open a support ticket. Nothing in this repo will fix it.
3. If the URL in the log looks malformed (e.g. a repo name that isn't
   `pbisim.git` / `pbisim-fit.git`), that is a different bug — capture the exact line, as
   both the `Dockerfile` and `.github/workflows/deploy.yml` spell those URLs literally.

---

## 3. `Authentication failed` / `403` / `could not read Username`

The token was supplied but GitHub rejected it. Check the token itself — run this
**locally**, in a terminal, with the token pasted into the variable:

```bash
TOKEN=github_pat_...        # paste; do not commit, do not echo into a file
curl -sS -o /dev/null -w '%{http_code}\n' \
  -H "Authorization: Bearer $TOKEN" \
  https://api.github.com/repos/phage-therapy-sim/pbisim
```

| Response | Meaning |
|---|---|
| `200` | Token is alive and can read `pbisim`. The problem is elsewhere — re-read §1. |
| `401` | Token is expired, revoked, or mistyped → mint a new one (§5). |
| `404` | Token is alive but has **no access to that repository** — the usual cause is a fine-grained token whose *Repository access* list omits the repo, or one still **pending org approval**. Fix the grant (§5) rather than minting another. |

Repeat with `.../repos/phage-therapy-sim/pbisim-fit` — one token must cover **both** repos
unless you set `PBISIM_FIT_TOKEN` separately. The `Dockerfile` falls back to
`PBISIM_TOKEN` when `PBISIM_FIT_TOKEN` is unset, so a token scoped to `pbisim` **only**
builds the engine fine and then fails on the `pbisim-fit` step.

---

## 4. `no open ports detected`

Streamlit didn't bind `$PORT` before Render's scan timed out. The `Dockerfile` already
mitigates this by pre-compiling bytecode (`compileall` with `--invalidation-mode
unchecked-hash`) so cold start isn't spent compiling. If it reappears: confirm the service
is still on the **standard** plan (`render.yaml` pins it; the 512 MB free tier was
OOM-unusable), and that `healthCheckPath` is `/_stcore/health`.

---

## 5. Mint a new fine-grained PAT

Only if §3 showed `401`, or you are rotating deliberately.

1. GitHub → click your avatar → **Settings**.
2. Left sidebar, bottom → **Developer settings**.
3. **Personal access tokens → Fine-grained tokens** → **Generate new token**.
4. Fill in:
   - **Token name:** something traceable, e.g. `render-pbisim-app-deploy-2026-09`.
   - **Resource owner:** `phage-therapy-sim` (the **org**, not your personal account).
     If the org requires approval, the token is created in **pending** state and returns
     `404` on everything until an org owner approves it under
     *Org → Settings → Personal access tokens → Pending requests*.
   - **Expiration:** custom; 90 days is a reasonable cadence. Put the date in your
     calendar — a build that dies on expiry day is otherwise a surprise.
   - **Repository access:** *Only select repositories* → select **`pbisim`** *and*
     **`pbisim-fit`** (both, so one token covers both build steps).
   - **Permissions → Repository permissions → Contents: `Read-only`**. Nothing else. Do
     not grant write, and do not use a classic PAT.
5. **Generate token** and copy it (`github_pat_…`). It is shown **once**.
6. Verify it immediately with the `curl` check in §3 — expect `200` for both repos.

---

## 6. Where the credential lives — rotation checklist

Rotating means updating **every** copy. Miss one and you get a half-broken pipeline whose
failure mode is silent (see the last row).

| # | Location | How to update | Symptom if stale |
|---|---|---|---|
| 1 | Render service env var `PBISIM_TOKEN` | Dashboard → service → **Environment** → edit → Save | Build fails at the engine step (§1/§3) |
| 2 | Render service env var `PBISIM_FIT_TOKEN` | same — **only if** you use a separate per-repo token; otherwise leave unset | Build fails at the pbisim-fit step |
| 3 | GitHub Actions secret `PBISIM_TOKEN` in the **`pbisim-app`** repo | Repo → Settings → Secrets and variables → **Actions** | **Silent.** `deploy.yml`'s first step sets `configured=false` when the secret is empty and skips the whole job **green** — you lose the test gate with no red ✗. If the secret is present but invalid, the job fails loudly at the pip step instead. |
| 4 | Local git remotes on any machine whose remote URL embeds a PAT | switch that machine to SSH (preferred), or re-set the remote URL | `git push`/`pull` fails with auth errors on that machine only |

Also revoke the **old** token afterwards (GitHub → Fine-grained tokens → the old entry →
Delete), so a leaked copy stops being useful.

Never paste a token into a tracked file, a log, a commit message, or a Dockerfile. The
`Dockerfile` deliberately confines it to the **builder** stage and then deletes
`direct_url.json` files, so it never reaches the shipped image or `docker history`.

---

## 7. Redeploy and verify

1. Render Dashboard → service → **Manual Deploy** → **Clear build cache & deploy**.
   Always clear the cache after a token fix. Because `PBISIM_REF` is `main` (an unchanged
   ARG value), a cached engine layer can otherwise survive indefinitely and pin the image
   to a **stale `pbisim`** even though the build "succeeded".
2. Watch the build log for the engine step succeeding — the token appears masked:
   `pip install git+https://***@github.com/phage-therapy-sim/pbisim.git@main`.
3. When the service is live, open **Shell** on the service and confirm what actually
   shipped:

   ```bash
   python -c "import pbisim; print('pbisim', pbisim.__version__)"
   python -c "import pbisim_fit; print('pbisim-fit', pbisim_fit.__version__)"
   ```

4. Confirm the image has the engine commit you expect, not a cached older one. Check for a
   known-recent change rather than the version string (the version does not move on every
   fix). For example, the 2026-09 removal of the `min(1, ·)` cap on rising density signals:

   ```bash
   python -c "
   import inspect, pbisim.dormancy.transitions as t
   src = inspect.getsource(t.density_dependent_dormancy)
   print('capped (STALE engine)' if 'minimum(1.0' in src else 'uncapped (current)')"
   ```

5. `curl -fsS https://<your-service>.onrender.com/_stcore/health` → `ok`.

---

## 8. Recommended hardening

Not yet applied — each is a small change to this repo:

- **Pin `PBISIM_REF` to a commit SHA** in `render.yaml` instead of `main`. It busts the
  pip layer cache by construction, so every deploy provably contains a known engine
  instead of depending on whether the cache happened to be cleared. Cost: you bump the SHA
  when you want a new engine — which is the point.
- **Make the guard message name its source**, e.g. "set PBISIM_TOKEN in the Render
  dashboard (Environment tab) — see DEPLOY_TROUBLESHOOTING.md §1". Saves a debugging cycle
  next time.
- **Calendar the token expiry** at mint time (§5), so rotation is planned rather than
  discovered from a failed deploy.
