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
READS_PER_JUNCTION = 20       # fusion pairs per junction (arriba filters
                              # low-support events; 5/junction produced
                              # zero reported rows live)

COMPL = str.maketrans("ACGTN", "TGCAN")
TRANS = str.maketrans("ACGT", "GTCA")  # transversion for per-read uniqueness


def revcomp(seq: str) -> str:
    return seq.translate(COMPL)[::-1]


def unique_snp(seq: str, pos: int) -> str:
    """One transversion at `pos` — makes otherwise-identical synthetic
    reads distinct (arriba dedups exact duplicates: 20 -> 1, live)."""
    return seq[:pos] + seq[pos].translate(TRANS) + seq[pos + 1 :]


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
    # The furthest mate slice = j + 12000 + (N-1)*25 + 150 — keep it
    # inside the contig (an out-of-range slice silently produced short
    # reads and a malformed fastq, live). The donor must lie inside
    # mini_gene1/2 (chr21:50-250 / 400-700) and the acceptor inside
    # mini_gene3 (chrX:400-700): intergenic or boundary breakpoints get
    # classified as in-vitro/end-to-end artifacts and filtered (live).
    max_j = len(acceptor) - 12150 - (READS_PER_JUNCTION - 1) * 25
    donor_spans = [(50, 250), (400, 700)]
    for lo, hi in donor_spans:
        for i in range(max(lo, 2), min(hi - SEG1, len(donor) - 150)):
            if donor[i + SEG1 - 1 : i + SEG1 + 1] != "GT":
                continue
            for j in range(400, min(700 - SEG2, max_j)):
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
            for k in range(READS_PER_JUNCTION):
                n += 1
                name = f"fusion_{jx}_{k}"
                # The mate maps on the ACCEPTOR contig, >=12 kb away from
                # the junction: arriba discards chimeric reads whose mate
                # sits <=10 kb from the breakpoint as read-through
                # fragments (live: 1605 -> 4 alignments, zero rows), and
                # chrX is 20 kb exactly so this structure can exist at
                # all. A donor-side mate is worse: STAR resolves the pair
                # as proper and never attempts chimera detection (live:
                # 0 chimeric reads).
                #
                # Per-read variation in TWO dimensions, because arriba
                # dedups by alignment signature, not sequence: (a) the
                # mate offset drifts by k*25 bp (insert-size spread,
                # real-library-like; identical mate positions collapsed
                # 20 reads to 1 live), (b) one transversion SNP in the
                # DONOR segment (a SNP in the 50 bp acceptor tail costs
                # the chimera: 15/20 survived with SNPs spread over the
                # full read, live).
                mate = revcomp(
                    acceptor[j + 12000 + k * 25 : j + 12150 + k * 25])
                r1_seq = unique_snp(seg1 + seg2, (k * 7) % SEG1)
                r2_seq = unique_snp(mate, (k * 7) % 150)
                q1 = "I" * (SEG1 + SEG2)
                q2 = "I" * 150
                r1.write(f"@{name}/1\n{r1_seq}\n+\n{q1}\n".encode())
                r2.write(f"@{name}/2\n{r2_seq}\n+\n{q2}\n".encode())
    print(f"appended {n} fusion pairs over {FUSION_COUNT} junctions:")
    for jx, (i, j) in enumerate(picks):
        print(f"  fusion_{jx}: chr21@{i}+{SEG1} -> chrX@{j}+{SEG2}")
    print(f"wrote {r1_out.name} / {r2_out.name}")


if __name__ == "__main__":
    main()
