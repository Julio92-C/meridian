#!/usr/bin/env bash
# Run the pipeline once per grouping variable (one config each), sequentially.
# Run from the repo root after combine_inputs.R:  bash projects/lung_microbiome_final/run_all.sh [slug ...]
# Logs: scripts/_run_lungfinal_<slug>_<timestamp>.log
set -u
RSCRIPT="${RSCRIPT:-/c/Program Files/R/R-4.5.0/bin/Rscript.exe}"
export RSTUDIO_PANDOC="${RSTUDIO_PANDOC:-C:/Users/julio/AppData/Local/Programs/Quarto/bin/tools}"
DIR=projects/lung_microbiome_final
SLUGS=("$@")
[[ ${#SLUGS[@]} -gt 0 ]] || SLUGS=(sample_type cohort immunodeficiency immunocompromised_type sex bronchiectasis macrolide age_band)
TS=$(date +%Y%m%d_%H%M%S)
for s in "${SLUGS[@]}"; do
    cfg=$DIR/config_$s.yaml; [[ $s == sample_type ]] && cfg=$DIR/config.yaml
    log=scripts/_run_lungfinal_${s}_$TS.log
    echo "[$(date '+%F %T')] start $s ($cfg)"
    "$RSCRIPT" run_pipeline.R "$cfg" > "$log" 2>&1; rc=$?
    echo "[$(date '+%F %T')] end   $s exit=$rc -> $log ($(grep -m1 'Total' "$log" | tr -s ' '))"
done
