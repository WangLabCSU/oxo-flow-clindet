#!/usr/bin/env bash
set -euo pipefail

echo "[INFO] Post-deploy started at $(date)"
echo "[INFO] CONDA_PREFIX=${CONDA_PREFIX}"

export R_REMOTES_NO_ERRORS_FROM_WARNINGS=true

Rscript --vanilla - <<'RSCRIPT'
options(warn = 0)

cat("[INFO] R version:", R.version.string, "\n")
cat("[INFO] warn option:", getOption("warn"), "\n")
cat("[INFO] R_REMOTES_NO_ERRORS_FROM_WARNINGS:",
    Sys.getenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"), "\n")

cat("[INFO] Installing gpgr from GitHub...\n")
remotes::install_github(
  "umccr/gpgr",
  ref = "main",
  upgrade = "never"
)

cat("[INFO] gpgr installation finished.\n")
RSCRIPT

echo "[INFO] Post-deploy finished at $(date)"