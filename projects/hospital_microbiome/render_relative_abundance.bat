@echo off
REM Standalone hospital relative-abundance figures (species top30 / genus top10 per group).
REM Lives in projects\hospital_microbiome; runs from the repo root two levels up.
pushd "%~dp0..\.."
> projects\hospital_microbiome\_relab_polish.log 2>&1 (
  echo === START %DATE% %TIME%
  "C:\Program Files\R\R-4.5.0\bin\Rscript.exe" projects/hospital_microbiome/polish_relative_abundance.R
  echo === DONE
)
popd
