#!/usr/bin/env python3
"""Append synthetic chimeric (fusion) read pairs to the RNA fixture.

The RNA sub-workflow's arriba_fusion rule requires chimeric alignments:
on the tiny fixture reference, wgsim reads can never produce split
reads, so arriba dies with "no split reads or discordant mates found".
This script builds INTER-CHROMOSOMAL fusion reads (chr21 donor ->
chrX acceptor) directly from the reference: two perfect-match segments
(100 bp donor + 50 bp acceptor) spliced at a GT/AG genomic junction
(so STAR's --chimScoreJunctionNonGTAG 0 does not zero the score).
Intra-chromosomal junctions do NOT work: arriba filters close
same-chromosome events as deletions (observed live: exit 0 but zero
rows in the fusion TSV).

Output: test/fixtures/reads/mini-T_RNA_R{1,2}.fq.gz — the existing
mini-T wgsim reads (kept for realism) plus 25 fusion pairs spread over
5 distinct junctions.  The DNA-side fixture files are untouched; the
#21 live-PASS evidence on mini-T_*.fq.gz stays valid.

Usage: python3 test/generate_rna_fusions.py
"""
import gzip
import random
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REFS = ROOT / "test/fixtures/refs/sequence/Homo_sapiens_assembly38_chr21.fasta"
READS = ROOT / "test/fixtures/reads"

SEG1, SEG2 = 100, 50          # donor / acceptor segment lengths
FUSION_COUNT = 5              # distinct junctions
READS_PER_JUNCTION = 5        # fusion pairs per junction

COMPL = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq: str) -> str:
    return seq.translate(COMPL)[::-1]


def load_contigs() -> dict[str, str]:
    contigs, name, chunks = {}, None, []
    with REFS.open() as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    contigs[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
        if name:
            contigs[name] = "".join(chunks)
    return contigs


def find_junctions(donor: str, acceptor: str) -> list[tuple[int, int]]:
    """Return (i, j) with donor GT at donor[i+99:i+101] and acceptor AG
    at acceptor[j-2:j]. Inter-chromosomal: no distance constraint."""
    found = []
    for i in range(len(donor) - SEG1):
        if donor[i + SEG1 - 1 : i + SEG1 + 1] != "GT":
            continue
        for j in range(2, len(acceptor) - 150):
            if acceptor[j - 2 : j] != "AG":
                continue
            found.append((i, j))
    return found


def main() -> None:
    contigs = load_contigs()
    assert "chr21" in contigs and "chrX" in contigs, contigs.keys()
    donor, acceptor = contigs["chr21"], contigs["chrX"]
    junctions = find_junctions(donor, acceptor)
    if len(junctions) < FUSION_COUNT:
        print(f"only {len(junctions)} GT/AG junction pairs found — "
              f"need {FUSION_COUNT}", file=sys.stderr)
        sys.exit(1)
    rng = random.Random(2026)
    picks = rng.sample(junctions, FUSION_COUNT)

    r1_out = READS / "mini-T_RNA_R1.fq.gz"
    r2_out = READS / "mini-T_RNA_R2.fq.gz"
    # Start from copies of the wgsim reads; the chimeric pairs are appended.
    shutil.copyfile(READS / "mini-T_R1.fq.gz", r1_out)
    shutil.copyfile(READS / "mini-T_R2.fq.gz", r2_out)

    n = 0
    with gzip.open(r1_out, "ab") as r1, gzip.open(r2_out, "ab") as r2:
        for jx, (i, j) in enumerate(picks):
            seg1 = donor[i : i + SEG1]
            seg2 = acceptor[j : j + SEG2]
            # The mate maps at the ACCEPTOR locus: a mate sitting next to
            # the donor makes STAR resolve the pair as a proper pair and
            # never attempt chimera detection (observed live: 100M50S with
            # zero chimeric reads).  A discordant mate is also the real
            # fusion biology (mates flank the junction from both sides).
            mate = revcomp(acceptor[j : j + 150])
            for k in range(READS_PER_JUNCTION):
                n += 1
                name = f"fusion_{jx}_{k}"
                q1 = "I" * (SEG1 + SEG2)
                q2 = "I" * 150
                r1.write(f"@{name}/1\n{seg1 + seg2}\n+\n{q1}\n".encode())
                r2.write(f"@{name}/2\n{mate}\n+\n{q2}\n".encode())
    print(f"appended {n} fusion pairs over {FUSION_COUNT} junctions:")
    for jx, (i, j) in enumerate(picks):
        print(f"  fusion_{jx}: chr21@{i}+{SEG1} -> chrX@{j}+{SEG2}")
    print(f"wrote {r1_out.name} / {r2_out.name}")


if __name__ == "__main__":
    main()
