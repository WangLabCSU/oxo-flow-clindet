#!/usr/bin/env python3
"""Generate the CNV mini-fixture reference files for the 922 bp chr21
fixture reference (test/fixtures/refs/sequence/Homo_sapiens_assembly38_chr21.fasta).

Outputs (committed, tiny):
  test/fixtures/cnv/ascat_loci/21        ASCAT loci file   (chr pos rsID)
  test/fixtures/cnv/ascat_alleles/21     ASCAT alleles file (chr pos rsID A B)
  test/fixtures/cnv/ascat_gc.txt         ASCAT GC-content file (chr pos GC)
  test/fixtures/cnv/freec_chrlen.txt     Control-FREEC chrLenFile
  test/fixtures/cnv/sequenza_gc.wig      sequenza-utils-style GC wiggle (50 bp)

These are SYNTHETIC smoke-test data: loci are spaced across the 922 bp
reference with invented rsIDs and reference-derived alleles. Real runs use
the upstream reference trees (ASCAT G1000 loci/alleles, freec chrLenFile,
sequenza gc50 wiggle) — see README.
"""
import gzip
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
REF = ROOT / "test/fixtures/refs/sequence/Homo_sapiens_assembly38_chr21.fasta"
OUT = ROOT / "test/fixtures/cnv"
CHROM = "chr21"

COMP = {"A": "G", "C": "T", "G": "A", "T": "C", "N": "A"}


def read_refs():
    """Parse the multi-contig fixture fasta into {contig: sequence}."""
    seqs, name, parts = {}, None, []
    with open(REF) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(parts).upper()
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
        if name is not None:
            seqs[name] = "".join(parts).upper()
    return seqs


def gc_window(seq, start, size=50):
    seg = seq[start - 1:start - 1 + size]
    if not seg:
        return 0.0
    gc = sum(1 for b in seg if b in "GC")
    return round(gc / len(seg), 6)


def main():
    seqs = read_refs()
    seq = seqs[CHROM]
    length = len(seq)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "ascat_loci").mkdir(exist_ok=True)
    (OUT / "ascat_alleles").mkdir(exist_ok=True)

    # loci every ~90 bp inside the chr21 contig (the other contigs of the
    # fixture reference — the RNA fusion test chrX — carry no loci)
    positions = list(range(90, length, 90))
    loci_lines, alleles_lines, gc_lines = [], [], []
    for i, pos in enumerate(positions):
        rsid = f"rsMINI{i + 1:03d}"
        ref_base = seq[pos - 1]
        alt = COMP[ref_base]
        loci_lines.append(f"{CHROM}\t{pos}\t{rsid}")
        alleles_lines.append(f"{CHROM}\t{pos}\t{rsid}\t{ref_base}\t{alt}")
        gc_lines.append(f"{CHROM}\t{pos}\t0.50")

    (OUT / "ascat_loci" / "21.txt").write_text("\n".join(loci_lines) + "\n")
    (OUT / "ascat_alleles" / "21.txt").write_text("\n".join(alleles_lines) + "\n")
    (OUT / "ascat_gc.txt").write_text("\n".join(gc_lines) + "\n")

    # exomedepth regions bed: the 2-region exome bed is degenerate for
    # ExomeDepth's reference-set logic (somatic.CNV.call errors with <3
    # targets); 10 synthetic regions across the chr21 contig make the
    # mini test mechanically valid
    bed_lines = [
        f"{CHROM}\t{start}\t{start + 50}"
        for start in range(50, length - 60, 80)
    ]
    (OUT / "exomedepth_regions.bed").write_text("\n".join(bed_lines) + "\n")

    # freec chrLenFile (tab-separated chr<TAB>len, like a .fai without offsets)
    (OUT / "freec_chrlen.txt").write_text(
        "".join(f"{name}\t{len(s)}\n" for name, s in seqs.items())
    )

    # sequenza-style fixedStep wiggle, 50 bp windows on chr21 (same shape as
    # sequenza-utils gc_wiggle -w 50 output)
    wig = [f"fixedStep chrom={CHROM} start=1 step=50 span=50"]
    for start in range(1, length + 1, 50):
        wig.append(str(gc_window(seq, start, 50)))
    (OUT / "sequenza_gc.wig").write_text("\n".join(wig) + "\n")

    print(f"wrote {len(positions)} loci for a {length} bp {CHROM} contig "
          f"(reference contigs: {list(seqs)}) -> {OUT}")


if __name__ == "__main__":
    main()
