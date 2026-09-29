# chicken_batch1 — self-contained test fixture

Committed pipeline inputs for the **published, open** chicken_batch1 cohort
(NCBI BioProject **PRJNA1406192**), so MERIDIAN can be run end-to-end from a
fresh clone or a container with **no external data staging**. These are the
pre-computed upstream outputs the pipeline consumes (Kraken2/Bracken, ABRicate,
Re-centrifuge, wf-metagenomics diversity) plus sample metadata — not raw reads.

Total size ~5.6 MB. This is the manuscript's primary validation cohort, so a
run here reproduces the paper's chicken_batch1 numbers.

## Layout

Canonical project_root structure — two input folders (`Metadata/` + `Reports/`);
`Datasets/` is an OUTPUT dir generated under `test_run/`, never an input.

```
Metadata/chicken_metadata1.csv
Reports/Abricate/summary_reportClean.csv
Reports/Bracken/bracken_arranged.csv
Reports/Kraken2/kraken2_db1_combined_reports.txt
Reports/Re-centrifuge/*.rcf.data.csv
Reports/wf-metagenomics/wf-metagenomics-diversity.csv
```

## Run it

Native, from the repository root (outputs land in `test_run/`, gitignored):

```bash
Rscript run_pipeline.R tests/fixtures/chicken_batch1/config.yaml
```

Apptainer on HPC (UCL Myriad). The container rootfs is **read-only**, so the
inputs are staged to a writable Scratch dir and `MERIDIAN_PROJECT_ROOT` points
the pipeline there (that env var overrides the config's `project_root`):

```bash
module load apptainer
export APPTAINER_CACHEDIR=~/Scratch/.apptainer                 # keep layers off home quota
apptainer build ~/Scratch/meridian_v1.0.0.sif docker://ghcr.io/julio92-c/meridian:v1.0.0

# stage a writable copy of the baked-in fixture, then run against it
mkdir -p ~/Scratch/ck_test
apptainer exec ~/Scratch/meridian_v1.0.0.sif \
  cp -a /opt/meridian/tests/fixtures/chicken_batch1/. ~/Scratch/ck_test/
apptainer run --pwd /opt/meridian \
  --env MERIDIAN_PROJECT_ROOT=/data -B ~/Scratch/ck_test:/data \
  ~/Scratch/meridian_v1.0.0.sif tests/fixtures/chicken_batch1/config.yaml
# outputs -> ~/Scratch/ck_test/test_run/
```

`--pwd /opt/meridian` is required so `run_pipeline.R` and `.Rprofile` (renv
activation) resolve; without it Apptainer runs in the host `$PWD`.

## Notes

- `config.yaml` uses a **repo-relative** `project_root`; override it with the
  `MERIDIAN_PROJECT_ROOT` env var to point at any bind-mounted study directory.
- No `taxid_fixes.csv` ships with this cohort (none exists); `taxid_fixes_file`
  is left null.
