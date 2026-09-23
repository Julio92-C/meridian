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

# 1. Quarto (pinned) — not in the Bioc base; run_pipeline.R calls `quarto render`.
ARG QUARTO_VERSION=1.5.57
RUN curl -sL https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb \
      -o /tmp/quarto.deb && dpkg -i /tmp/quarto.deb && rm /tmp/quarto.deb

# 2. renv restore — copy ONLY the lockfile first so this layer caches across
#    code edits and only re-runs when renv.lock actually changes.
WORKDIR /opt/meridian
COPY renv.lock ./
COPY renv/activate.R renv/settings.json ./renv/
RUN R -e "install.packages('renv'); renv::restore(lockfile='renv.lock', prompt=FALSE)"

# 3. Copy the rest of the repo (code edits invalidate only from here down).
COPY . /opt/meridian

# 4. Build-time smoke check: toolchain + keystone packages load.
RUN quarto --version && \
    R -e "stopifnot(nzchar(Sys.which('quarto'))); library(ALDEx2); library(quarto)"

ENTRYPOINT ["Rscript", "run_pipeline.R"]
