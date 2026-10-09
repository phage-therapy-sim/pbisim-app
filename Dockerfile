# syntax=docker/dockerfile:1
#
# pbisim-app container image.
#
# The private engine `pbisim` is installed from its private GitHub repo using a
# READ-ONLY token supplied at build time. The token lives only in the *builder*
# stage (a build ARG), so it never appears in the final image's layers or
# `docker history`. NEVER hard-code the token here or commit it — set it as a
# secret environment variable on the host (see DEPLOY.md).
#
# Build args:
#   PBISIM_TOKEN  fine-grained GitHub PAT with read-only access to the pbisim repo
#
# Engine versions are PINNED to exact commits in deploy/pbisim.ref and
# deploy/pbisim-fit.ref (one 40-char SHA each). Those files are COPYed in right
# before the install steps, so bumping a SHA invalidates exactly those layers and
# the image provably contains the named commits. (Previously the install used
# `@main` with a constant PBISIM_REF build arg, so Docker reused a cached layer
# and deploys silently kept a stale engine while the app code updated.)
# To deploy a new engine: push it, then put its SHA in deploy/pbisim.ref and push
# the app. PBISIM_REF_OVERRIDE / PBISIM_FIT_REF_OVERRIDE (unset by default) force a
# different ref for a one-off local build.

########################  builder  ########################
FROM python:3.11-slim AS builder

ENV PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_NO_CACHE_DIR=1

RUN apt-get update \
    && apt-get install -y --no-install-recommends git \
    && rm -rf /var/lib/apt/lists/*

RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# --- 1. Private engine, pinned (layer rebuilds whenever deploy/pbisim.ref changes) ---
COPY deploy/pbisim.ref deploy/pbisim-fit.ref /tmp/refs/
ARG PBISIM_TOKEN
ARG PBISIM_REF_OVERRIDE=
RUN test -n "$PBISIM_TOKEN" || (echo "ERROR: PBISIM_TOKEN build arg is required — set it in the Render dashboard (Environment tab); see DEPLOY_TROUBLESHOOTING.md §1" >&2; exit 1) \
    && REF="${PBISIM_REF_OVERRIDE:-$(tr -d '[:space:]' < /tmp/refs/pbisim.ref)}" \
    && echo "Installing pbisim @ ${REF}" \
    && pip install "git+https://${PBISIM_TOKEN}@github.com/phage-therapy-sim/pbisim.git@${REF}" \
    && echo "${REF}" > /opt/venv/pbisim.installed-ref

# --- 1b. Private estimation package `pbisim-fit` (powers the Calibration NLS fit).
#         Its CORE is torch-free (torch lives only in the [sbi]/[neural_ode] extras,
#         which we do NOT install), so this stays lightweight and fits the free tier.
#         Reuses PBISIM_TOKEN when PBISIM_FIT_TOKEN is unset — set a separate token
#         only if your fine-grained PAT is scoped to a single repo. ---
ARG PBISIM_FIT_TOKEN=
ARG PBISIM_FIT_REF_OVERRIDE=
RUN FIT_TOKEN="${PBISIM_FIT_TOKEN:-$PBISIM_TOKEN}" \
    && test -n "$FIT_TOKEN" || (echo "ERROR: PBISIM_FIT_TOKEN or PBISIM_TOKEN is required" >&2; exit 1) \
    && FIT_REF="${PBISIM_FIT_REF_OVERRIDE:-$(tr -d '[:space:]' < /tmp/refs/pbisim-fit.ref)}" \
    && echo "Installing pbisim-fit @ ${FIT_REF}" \
    && pip install "git+https://${FIT_TOKEN}@github.com/phage-therapy-sim/pbisim-fit.git@${FIT_REF}" \
    && echo "${FIT_REF}" > /opt/venv/pbisim-fit.installed-ref

# --- 2. App source + remaining deps (editable so prompts/ resolves at runtime) ---
# [scripting] adds the code-editor component for the opt-in Scripting page; [onnx] adds
# onnxruntime for the Calibration page's instant amortized fit (torch-free). Both are
# lazy-imported, so they cost nothing when unused and degrade gracefully if absent.
WORKDIR /app
COPY . /app
RUN pip install -e '.[scripting,onnx]'

# --- 3. Scrub any token left in VCS install metadata ---
RUN find /opt/venv -name direct_url.json -delete || true

# --- 4. Pre-compile bytecode so the runtime container starts fast ---
# Without this the image ships no .pyc, so every boot compiles all of
# Streamlit's + pbisim's bytecode on first import — far too slow on the free
# plan's fractional CPU, so Streamlit misses Render's port-scan window and the
# deploy fails with "no open ports detected". Compiling here (on the build
# instance's full CPU) bakes the .pyc into the image → fast cold start.
# `unchecked-hash` makes the .pyc always valid regardless of source mtime, so the
# cross-stage COPY can't invalidate them and force a runtime recompile.
RUN python -m compileall -q -j 0 --invalidation-mode unchecked-hash /opt/venv /app || true

########################  runtime  ########################
FROM python:3.11-slim AS runtime

ENV PYTHONUNBUFFERED=1 \
    PATH="/opt/venv/bin:$PATH" \
    # headless matplotlib backend (belt-and-suspenders with matplotlib.use in app.py)
    MPLBACKEND=Agg \
    # cap native (BLAS/OMP) thread pools — multi-threaded BLAS forked by joblib on a
    # small shared instance is a common SIGSEGV (exit 139) source
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    NUMEXPR_NUM_THREADS=1

RUN useradd --create-home --uid 10001 appuser

# venv + source (the editable install links pbisim_app -> /app)
COPY --from=builder /opt/venv /opt/venv
COPY --from=builder /app /app

WORKDIR /app
USER appuser

# Render/most hosts inject $PORT; default to 8501 for local `docker run`.
ENV PORT=8501
EXPOSE 8501

CMD streamlit run pbisim_app/app.py \
    --server.port "${PORT}" \
    --server.address 0.0.0.0 \
    --server.headless true \
    --browser.gatherUsageStats false
