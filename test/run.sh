#!/usr/bin/env bash
# Acceptance test for the oxo-flow-clindet port (zyllifeworld/clindet @ 582a9131).
# Single-entry workflow (upstream Snakefile form): run_type selects
# wes / wgs / rna; paired vs unpaired WES derives from the sample sheet
# (wildcard.control), so each run type validates + dry-runs here.
# Usage: ./test/run.sh            (uses ./main.oxoflow)
#        OXO=/path/to/oxo-flow ./test/run.sh
set -euo pipefail
cd "$(dirname "$0")/.."
OXO=${OXO:-oxo-flow}

# Only run/dry-run take --arg overrides; validate/lint/debug are static
# structure checks (run_type-gated behavior is exercised by the dry-run loop).
"$OXO" validate main.oxoflow
"$OXO" lint main.oxoflow > /tmp/oxo-lint-$$.txt 2>&1
"$OXO" debug main.oxoflow 2>&1 | grep -E '\{config\.|\{pair_id\}' > /dev/null && {
  echo "unexpanded wildcards in debug output"
  exit 1
}

for ARGS in "" "--arg run_type=wgs" "--arg run_type=rna"; do
  LABEL=${ARGS:-default-wes}
  echo "==> [$LABEL] dry-run"
  "$OXO" dry-run main.oxoflow $ARGS > /tmp/oxo-dryrun-$$.txt 2>&1
  grep -q "would execute" /tmp/oxo-dryrun-$$.txt
done

echo "PASS"
