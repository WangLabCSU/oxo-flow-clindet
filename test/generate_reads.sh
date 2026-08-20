#!/usr/bin/env bash
# Regenerate the synthetic mini reads under fixtures/reads/.
#
# 10k pairs per sample, 150bp, ~1% error, 500bp mean insert — sized so that
# Manta's GetAlignmentStats has >=100 high-confidence read pairs per read
# group (the original 2-pair fixtures died there live:
# "Too few high-confidence read pairs (0) to determine pair orientation").
#
# Requires wgsim: conda install -c bioconda wgsim
set -euo pipefail
cd "$(dirname "$0")/fixtures/reads"
REF=../refs/sequence/Homo_sapiens_assembly38_chr21.fasta
command -v wgsim >/dev/null || { echo "wgsim not found (conda install -c bioconda wgsim)"; exit 1; }
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
wgsim -N 10000 -1 150 -2 150 -e 0.01 -d 500 -s 50 -S 42  "$REF" "$TMP/mini-T_R1.fq"  "$TMP/mini-T_R2.fq"  >/dev/null
wgsim -N 10000 -1 150 -2 150 -e 0.01 -d 500 -s 50 -S 1337 "$REF" "$TMP/mini-NC_R1.fq" "$TMP/mini-NC_R2.fq" >/dev/null
for f in mini-T_R1 mini-T_R2 mini-NC_R1 mini-NC_R2; do
  gzip -9 -c "$TMP/$f.fq" > "$f.fq.gz"
done
echo "regenerated: $(pwd)"
ls -la *.fq.gz
