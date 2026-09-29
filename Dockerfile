# MERIDIAN — containerised analysis/reporting layer.
#
# Packages run_pipeline.R -> R/NN_*.R stage modules -> Quarto HTML report as a
# single version-tagged image. Upstream classification tools (basecalling,
# wf-metagenomics, Kraken2/Bracken, ABRicate, Re-centrifuge) are NOT included:
# the pipeline consumes their pre-computed CSV/text outputs. See README
# "Running in a container".
#
# Build:  docker build -t meridian:v1.1.0 .
# Run:    docker run --rm -e MERIDIAN_PROJECT_ROOT=/data \
#           -v /host/study_dir:/data meridian:v1.1.0 projects/chicken_batch1/config.yaml

FROM bioconductor/bioconductor_docker:RELEASE_3_22

# 1. Quarto CLI (pinned) — not in the Bioc base; run_pipeline.R shells out to
#    `quarto render` via system2() (the CLI, not the quarto R package).
ARG QUARTO_VERSION=1.5.57
RUN curl -sL https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb \
      -o /tmp/quarto.deb && dpkg -i /tmp/quarto.deb && rm /tmp/quarto.deb

# 2. renv restore. Copy .Rprofile + the renv bootstrap files FIRST so that
#    starting R auto-activates the project (via renv/activate.R) and restore
#    installs into the project library (renv/library) that the pipeline uses at
#    runtime. Without .Rprofile here, restore runs un-activated and the runtime
#    library ends up empty. Copying only these before the full source keeps this
#    expensive ~283-package layer cached until renv.lock changes.
#
#    CACHE DISABLED ON PURPOSE: by default renv symlinks the project library into
#    a global cache under $HOME. In a read-only Apptainer/Singularity container
#    those symlinks break at run time (rootfs read-only + $HOME differs from build),
#    so renv thinks packages "need reinstalling" and fails writing to the read-only
#    library. Disabling the cache makes restore copy real package files into the
#    library, so the image is self-contained and load() never writes.
ENV RENV_CONFIG_CACHE_ENABLED=FALSE \
    RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE
WORKDIR /opt/meridian
COPY .Rprofile renv.lock ./
COPY renv/activate.R renv/settings.json ./renv/
RUN R -e "renv::restore(prompt = FALSE)"

# 3. Copy the rest of the repo (code edits invalidate only from here down).
COPY . /opt/meridian

# 4. Build-time smoke check: Quarto CLI on PATH + keystone renv package loads
#    through the activated project library. Run via Rscript to mirror the
#    ENTRYPOINT exactly. The pipeline uses the quarto CLI, not the quarto R
#    package, so we deliberately do not library(quarto).
RUN quarto --version && \
    Rscript -e "cat('libPaths:', .libPaths(), sep='\n'); stopifnot(nzchar(Sys.which('quarto'))); library(ALDEx2); cat('smoke check OK\n')"

ENTRYPOINT ["Rscript", "run_pipeline.R"]
