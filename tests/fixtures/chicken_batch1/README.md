# chicken_batch1 — self-contained test fixture

Committed pipeline inputs for the **published, open** chicken_batch1 cohort
(NCBI BioProject **PRJNA1406192**), so MERIDIAN can be run end-to-end from a
fresh clone or a container with **no external data staging**. These are the
pre-computed upstream outputs the pipeline consumes (Kraken2/Bracken, ABRicate,
Re-centrifuge, wf-metagenomics diversity) plus sample metadata — not raw reads.

Total size ~5.6 MB. This is the manuscript's primary validation cohort, so a
run here reproduces the paper's chicken_batch1 numbers.

## Run it

From the repository root (outputs land in `test_run/`, which is gitignored):

```bash
Rscript run_pipeline.R tests/fixtures/chicken_batch1/config.yaml
```

Container (image baked with the repo, so the fixture is already inside):

```bash
docker run --rm -v "$PWD/out:/opt/meridian/tests/fixtures/chicken_batch1/test_run" \
  ghcr.io/julio92-c/meridian:v1.0.0 tests/fixtures/chicken_batch1/config.yaml
```

Apptainer on HPC (UCL Myriad — build the .sif on a login node):

```bash
module load apptainer
export APPTAINER_CACHEDIR=~/Scratch/.apptainer
apptainer build meridian_v1.0.0.sif docker://ghcr.io/julio92-c/meridian:v1.0.0
apptainer run meridian_v1.0.0.sif tests/fixtures/chicken_batch1/config.yaml
```

## Notes

- `config.yaml` uses a **repo-relative** `project_root`, resolved from the
  current working directory. Override it to point at your own bind-mounted
  study directory with the `MERIDIAN_PROJECT_ROOT` environment variable.
- No `taxid_fixes.csv` ships with this cohort (none exists); `taxid_fixes_file`
  is left null.
