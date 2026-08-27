#!/usr/bin/env bash
# Acceptance test for the oxo-flow-clindet RNA sub-workflow port
# (zyllifeworld/clindet workflow/RNA). Usage: ./test/run_rna.sh
#        OXO=/path/to/oxo-flow ./test/run_rna.sh
set -euo pipefail
cd "$(dirname "$0")/.."
OXO=${OXO:-oxo-flow}

# The upstream mini-test default stages (arriba + call_mut) as rule targets:
# the RSEM/kallisto/salmon/TRUST4 quant rules and isofox are NOT in the
# default stages upstream.
TARGETS=(
  -t fastp_trim -t STAR_1_pass -t STAR_arriba_map -t STAR_mut_map
  -t arriba_fusion -t link_bam -t SplitNCigarReads
  -t mutect2 -t M2_filter -t unpaired -t call_variants -t lofreq
  -t varscan2 -t norm_filter
)

# validate/lint/debug are static structure checks and do not accept --arg
# or -t (run_type-gated behavior is exercised by the targeted dry-run below;
# run.sh covers the same wildcard-placeholder check for the default WES).
echo "==> validate"
"$OXO" validate main.oxoflow

echo "==> lint (warnings are acceptable, errors are not)"
"$OXO" lint main.oxoflow

echo "==> dry-run with the upstream default stages"
"$OXO" dry-run main.oxoflow --arg run_type=rna "${TARGETS[@]}" > /tmp/oxo-dryrun-rna-$$.txt 2>&1
grep -q "would execute" /tmp/oxo-dryrun-rna-$$.txt

echo "PASS"
