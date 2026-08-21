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

echo "==> validate"
"$OXO" validate main_rna.oxoflow

echo "==> lint (warnings are acceptable, errors are not)"
"$OXO" lint main_rna.oxoflow

echo "==> dry-run with the upstream default stages"
"$OXO" dry-run main_rna.oxoflow "${TARGETS[@]}" > /tmp/oxo-dryrun-rna-$$.txt 2>&1
grep -q "would execute" /tmp/oxo-dryrun-rna-$$.txt

echo "==> debug: expanded commands contain no literal {config./{pair_id} placeholders"
if "$OXO" debug main_rna.oxoflow "${TARGETS[@]}" 2>&1 | grep -E '\{config\.|\{pair_id\}' > /dev/null; then
  echo "unexpanded wildcards in debug output"
  exit 1
fi

echo "PASS"
