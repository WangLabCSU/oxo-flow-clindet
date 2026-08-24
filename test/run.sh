#!/usr/bin/env bash
# Acceptance test for the oxo-flow-clindet port (zyllifeworld/clindet @ 582a9131).
# Usage: ./test/run.sh            (uses ./main.oxoflow)
#        OXO=/path/to/oxo-flow ./test/run.sh
set -euo pipefail
cd "$(dirname "$0")/.."
OXO=${OXO:-oxo-flow}

for WF in main.oxoflow main_unpaired.oxoflow; do
  echo "==> [$WF] validate"
  "$OXO" validate "$WF"

  echo "==> [$WF] lint (warnings are acceptable, errors are not)"
  "$OXO" lint "$WF" > /tmp/oxo-lint-$$.txt 2>&1

  echo "==> [$WF] dry-run with default config"
  "$OXO" dry-run "$WF" > /tmp/oxo-dryrun-$$.txt 2>&1
  grep -q "would execute" /tmp/oxo-dryrun-$$.txt

  echo "==> [$WF] debug: expanded commands contain no literal {config./{pair_id} placeholders"
  if "$OXO" debug "$WF" 2>&1 | grep -E '\{config\.|\{pair_id\}' > /dev/null; then
    echo "unexpanded wildcards in debug output"
    exit 1
  fi
done

echo "PASS"
