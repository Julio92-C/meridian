@echo off
REM Redraw UpSet PNGs, then rebuild composite panels for hospital_microbiome.
REM Output is logged to scripts\_render_hospital.log
pushd "%~dp0"
set RSTUDIO_PANDOC=C:/Users/julio/AppData/Local/Programs/Quarto/bin/tools
set RS="C:\Program Files\R\R-4.5.0\bin\Rscript.exe"
> scripts\_render_hospital.log 2>&1 (
  echo === START %DATE% %TIME%
  %RS% scripts/upset_only.R projects/hospital_microbiome/config.yaml
  %RS% scripts/panels_only.R projects/hospital_microbiome/config.yaml
  echo === DONE %DATE% %TIME%
)
popd
