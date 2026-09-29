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

**Container — `meridian.sh` (recommended).** The `meridian.sh` wrapper at the
repo root makes the container behave exactly like the native call — you pass
only a config. It reads `project_root` from the config, binds that host
directory into the container unchanged, and auto-detects Apptainer (HPC) or
Docker (laptop).

```bash
# 1. get the image once
apptainer build meridian_v1.0.0.sif docker://ghcr.io/julio92-c/meridian:v1.0.0

# 2. stage this fixture as a study directory (or use your own)
apptainer exec meridian_v1.0.0.sif cp -a /opt/meridian/tests/fixtures/chicken_batch1/. ./study/

# 3. point the config at that directory (absolute host path), then run
sed -i "s|^project_root:.*|project_root: \"$PWD/study\"|" ./study/config.yaml
./meridian.sh ./study/config.yaml
# outputs -> ./study/test_run/
```

`meridian.sh` env overrides: `MERIDIAN_SIF` (path to the `.sif`) and
`MERIDIAN_IMAGE` (Docker image ref, default `ghcr.io/julio92-c/meridian:v1.0.0`).

<details><summary>Manual invocation (what the wrapper runs)</summary>

```bash
apptainer run --pwd /opt/meridian -B "$PWD/study:$PWD/study" \
  meridian_v1.0.0.sif "$PWD/study/config.yaml"
# docker equivalent:
docker run --rm -w /opt/meridian -v "$PWD/study:$PWD/study" \
  ghcr.io/julio92-c/meridian:v1.0.0 "$PWD/study/config.yaml"
```

`--pwd /opt/meridian` (Apptainer) / `-w /opt/meridian` (Docker) is required so
`run_pipeline.R` and `.Rprofile` (renv activation) resolve inside the image.
</details>

## Notes

- The committed `config.yaml` uses a **repo-relative** `project_root` for native
  runs; for the container, set it to the **absolute** host study path (as above)
  or override with the `MERIDIAN_PROJECT_ROOT` env var.
- No `taxid_fixes.csv` ships with this cohort (none exists); `taxid_fixes_file`
  is left null.
