# MERIDIAN containerization — design spec

**Date:** 2026-09-17
**Status:** Implemented (files below), build not yet verified — 2026-09-23
**Branch context:** authored on `validation/wetlands-flyway`

## Implementation notes (2026-09-23)

All four files landed. Two deviations from the spec text, both deliberate:

1. **`.dockerignore` keeps `projects/`.** The spec listed `projects/*/` under
   exclusions, but in this repo those dirs hold only small `config.yaml` / prep
   scripts (no study data — data lives under the mounted `project_root`), and the
   "reproduce the manuscript" invocation bakes `projects/chicken_batch1/config.yaml`
   into the image. Excluding them would break that invocation.
2. **Dockerfile splits the renv `COPY`.** The spec's single
   `COPY renv.lock renv/activate.R renv/settings.json ./renv/` would land the
   lockfile at `./renv/renv.lock` (Docker multi-source COPY puts every source in
   the dest dir), but `renv::restore(lockfile='renv.lock')` and `renv/activate.R`
   expect it at the project root. Split into two `COPY` lines.

Still pending: Tier-2 build reality (`docker build`) + Tier-3 HPC full run — both
require Docker/Apptainer, done by Julio. Quarto pin left at ARG default 1.5.57.

## Goal

Package the MERIDIAN analysis pipeline (the R + Quarto analysis/reporting
layer in this repo) as a single, version-tagged container image so that:

1. **Manuscript reproducibility** — the image tag is a citable artifact that
   re-runs the exact analysis behind the Applications Note.
2. **HPC / collaborator portability** — the image converts cleanly to an
   Apptainer/Singularity `.sif` and runs against arbitrary study directories.
3. **Escape Windows toolchain friction** — a stable Linux environment removes
   the recurring `RSTUDIO_PANDOC` / Python-collision pain of the native setup.

## Scope

**In scope:** the core pipeline — `run_pipeline.R` → the `R/NN_*.R` stage
modules → the Quarto HTML report render. This is exactly what
`run_pipeline.R` does today.

**Out of scope:**
- Upstream classification tools (basecalling, wf-metagenomics, Kraken2/Bracken,
  Abricate, Re-centrifuge). These are *not* in this repo; the pipeline consumes
  their pre-computed CSV/text outputs. The README must state this boundary
  explicitly so the image is not oversold as "the whole metagenomics pipeline."
- Manuscript Figure-1 Python/matplotlib toolchain (`generate_figure1.py`).
- Word / COM tooling (Windows-only; tracked-diff, manuscript review).

## Non-goals

- Replacing the native development workflow. The `/tune-panels` ~45s loop,
  Word COM, and native `Rscript run_pipeline.R` all keep working untouched. The
  container is an **added run/publish target**, not a replacement dev
  environment.
- Byte-exact cross-OS output parity. Native runs happen on Windows; exact
  floating-point byte-parity across OSes is neither guaranteed nor required for
  the reproducibility claim.

## Architecture

A single `Dockerfile` at the repo root produces `meridian:<version>`, where the
tag mirrors the manuscript/release version (e.g. `meridian:v1.1.0`). That image
bundles, pinned together:

- **Base:** `bioconductor/bioconductor_docker:RELEASE_3_22` — ships R 4.5.0 and
  the system libraries the 57 Bioconductor packages (ALDEx2 + core Bioc chain)
  and spatial/graphics R packages compile against. Chosen over `rocker/r-ver`
  to avoid iterative "missing libX" build failures across 283 packages.
- **Quarto** (pinned via `ARG`) — not present in the Bioc base;
  `run_pipeline.R` shells out to `quarto render`.
- **Pipeline code** — baked in at build via `COPY`. Rebuilding the image is how
  code is updated; the tag pins code + all 283 renv packages + toolchain as one
  artifact.

Data and per-study configs are **not** baked in — they are supplied at runtime
via bind mounts + an environment variable (see Runtime contract).

HPC: Docker is the single source of truth. Apptainer users convert with
`apptainer build meridian_<version>.sif docker-daemon://meridian:<version>`.
No separate Singularity definition file is maintained.

## Build (Dockerfile)

```dockerfile
FROM bioconductor/bioconductor_docker:RELEASE_3_22

# 1. Quarto (pinned) — not in the Bioc base; run_pipeline.R calls `quarto render`
ARG QUARTO_VERSION=1.5.57
RUN curl -sL https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb \
      -o /tmp/quarto.deb && dpkg -i /tmp/quarto.deb && rm /tmp/quarto.deb

# 2. renv restore — copy ONLY the lockfile first so this layer caches across
#    code edits and only re-runs when renv.lock actually changes.
WORKDIR /opt/meridian
COPY renv.lock renv/activate.R renv/settings.json ./renv/
RUN R -e "install.packages('renv'); renv::restore(lockfile='renv.lock', prompt=FALSE)"

# 3. Copy the rest of the repo (code edits invalidate only from here down).
COPY . /opt/meridian

# 4. Build-time smoke check: toolchain + keystone packages load.
RUN quarto --version && \
    R -e "stopifnot(nzchar(Sys.which('quarto'))); library(ALDEx2); library(quarto)"

ENTRYPOINT ["Rscript", "run_pipeline.R"]
```

Key build decisions:

- **Layer ordering for cache reuse.** `renv.lock` is copied and restored before
  the full `COPY .`. Editing an R stage rebuilds only from step 3; the expensive
  ~283-package restore stays cached until the lockfile changes.
- **Quarto version pin.** Default to a recent 1.5.x; at implementation time,
  pin to whatever `quarto --version` reports in the environment that produced
  the current manuscript figures, if recoverable.
- **ENTRYPOINT `Rscript run_pipeline.R`** so `docker run meridian:<v> <config>`
  reads naturally — the config path is the container argument.
- **System-lib gaps** discovered during the real `renv::restore()` are fixed by
  adding `apt-get install` lines before step 2 (see Testing — build reality).

## Runtime contract

Exactly one code change is required, in `load_config()` (`R/00_setup.R`): an
environment-variable override applied **before** the `dir.exists()` assertion,
so the config's baked-in Windows host path is never validated in-container.

```r
load_config <- function(path) {
  cfg <- yaml::read_yaml(path)
  # Runtime override: MERIDIAN_PROJECT_ROOT wins over the config's host path,
  # so the SAME config.yaml runs unchanged natively and in-container.
  env_root <- Sys.getenv("MERIDIAN_PROJECT_ROOT", unset = "")
  if (nzchar(env_root)) cfg$project_root <- env_root
  stopifnot(
    !is.null(cfg$project_root),
    dir.exists(cfg$project_root),
    !is.null(cfg$metadata$file),
    !is.null(cfg$metadata$sample_id_col)
  )
  cfg
}
```

This is backward-compatible: with the env var unset, native behaviour is
identical. Because `inputs$*`, `metadata$file`, and `outputs$*` already resolve
*relative* to `project_root` via `file.path()`, overriding the root is
sufficient — no other path surgery is needed. Outputs land under the mounted
root, so results appear on the host.

**Contract in one line:** mount the study dir somewhere, point
`MERIDIAN_PROJECT_ROOT` at that container path, pass a config path as the
argument. The config's own `project_root:` value is irrelevant in-container.

### Invocation — reproduce the manuscript (config baked into the image)

```bash
docker run --rm \
  -e MERIDIAN_PROJECT_ROOT=/data \
  -v /host/path/to/study_dir:/data \
  meridian:v1.1.0 projects/chicken_batch1/config.yaml
```

### Invocation — new study (collaborator mounts their own config too)

```bash
docker run --rm \
  -e MERIDIAN_PROJECT_ROOT=/data \
  -v /host/study_dir:/data \
  -v /host/my_config.yaml:/project/config.yaml \
  meridian:v1.1.0 /project/config.yaml
```

## Testing & validation

**Tier 1 — build-time smoke (in Dockerfile).** `quarto --version` +
`library(ALDEx2)` + `nzchar(Sys.which('quarto'))`. Fails the build early if the
toolchain or a keystone package is broken.

**Tier 2 — build reality.** The real `renv::restore()` against a fresh Bioc base
is where a clean-looking Dockerfile can diverge from reality (a package needing
a system lib the base lacks). Implementation is not "done" until the image
**actually builds**; any missing system libs are added as `apt-get install`
lines and re-verified.

**Tier 3 — full-pipeline validation on HPC (performed by Julio).** No small
input fixture exists yet, so there is no automated functional test in v1.
Instead, validation is a documented full-run on HPC against a real study
directory, following the steps below.

### HPC validation steps (Apptainer/Singularity)

Prerequisites on the HPC login node: `apptainer` (or `singularity`) available;
the study directory (containing `Datasets/`, `Metadata/`, `Reports/`) staged on
a filesystem the compute nodes can read; enough scratch for the image.

1. **Build the Docker image** (on a machine with Docker — your workstation):
   ```bash
   docker build -t meridian:v1.1.0 .
   ```

2. **Export the image to a tar** so it can move to the cluster:
   ```bash
   docker save meridian:v1.1.0 -o meridian_v1.1.0.tar
   ```
   Copy `meridian_v1.1.0.tar` to the HPC (e.g. `scp`/`rsync`).

3. **Convert to a `.sif` on the HPC** (Docker is usually unavailable there;
   Apptainer reads the tar directly):
   ```bash
   apptainer build meridian_v1.1.0.sif docker-archive://meridian_v1.1.0.tar
   ```
   (If your workstation can push to a registry instead, you may skip steps 2–3
   and run `apptainer build meridian_v1.1.0.sif docker://<registry>/meridian:v1.1.0`.)

4. **Run the pipeline against a real study** — bind the study dir to `/data`
   and point the env var at it. Wrap in your scheduler's job script
   (SLURM example shown; adjust partition/resources):
   ```bash
   #!/bin/bash
   #SBATCH --job-name=meridian
   #SBATCH --cpus-per-task=8
   #SBATCH --mem=32G
   #SBATCH --time=02:00:00

   apptainer run \
     --env MERIDIAN_PROJECT_ROOT=/data \
     --bind /scratch/$USER/study_dir:/data \
     meridian_v1.1.0.sif projects/chicken_batch1/config.yaml
   ```
   For a **new study**, also bind your own config and pass its container path:
   ```bash
   apptainer run \
     --env MERIDIAN_PROJECT_ROOT=/data \
     --bind /scratch/$USER/study_dir:/data \
     --bind /scratch/$USER/my_config.yaml:/project/config.yaml \
     meridian_v1.1.0.sif /project/config.yaml
   ```

5. **Confirm outputs** appear under the bound study dir on the host
   (`/scratch/$USER/study_dir`): the manifest JSON, panel PNGs under the
   configured figures dir, and the Quarto HTML report. A zero exit code plus
   these artifacts constitutes a passing full-run.

Notes for HPC:
- Apptainer mounts `$HOME` and the CWD by default; `--bind` adds the study dir.
  If the site runs `--contain`-style isolation, add `--bind` for any other
  paths the config references.
- No root is needed for `apptainer run`. Building the `.sif` from a tar
  (step 3) also works rootless on most clusters.

## Housekeeping / docs

- **`.dockerignore`** (new): exclude `renv/library/` (host Windows binaries —
  must not leak into the Linux image), `projects/*/` data, `docs/manuscript/*.docx`,
  `docs/manuscript/_backup/`, `scratch/`, `.git/`, `*.log`, `.RData`, `.Rhistory`,
  `Rplots.pdf`.
- **README** — add a "Running in a container" section: the invocation patterns
  above, the Apptainer conversion, and the HPC validation steps.
- **CI** — out of scope for v1 (needs an input fixture, which does not yet
  exist). Revisit once a trimmed reference dataset is available.

## Files touched

| File | Change |
| --- | --- |
| `Dockerfile` | new — build definition (above) |
| `.dockerignore` | new — exclude data/host-arch/scratch |
| `R/00_setup.R` | edit — `MERIDIAN_PROJECT_ROOT` override in `load_config()` |
| `README.md` | edit — "Running in a container" + HPC steps |

## Open items (resolve at implementation time, not blocking)

1. **Quarto version pin** — recover the exact version behind current manuscript
   figures; else default recent 1.5.x.
2. **System-lib gaps** — surface via the real build; add `apt-get` lines.
3. **Input fixture** — none yet; deferred. Blocks future CI, not v1.
