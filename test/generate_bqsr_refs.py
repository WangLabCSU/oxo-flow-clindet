#!/usr/bin/env python3
"""Generate SYNTHETIC mini known-sites VCFs for the BQSR rules.

The upstream-known-sites files (1000G phase1 snps / Mills gold-standard
indels, hg38_chr21 subsets) carry the full-chromosome contig header
(##contig=<ID=chr21,length=46709983>), which GATK BaseRecalibrator rejects
against the 900 bp chr21 fixture reference ("contig features" mismatch).
This generator synthesizes two mini known-sites VCFs whose contig headers
match the fixture reference and whose sites sit inside the reference
contig — same role as the arriba mini DBs (documented synthetic fixtures;
real runs point known_sites1/2 at the real varanno files).

Outputs (committed, tiny):
  test/fixtures/refs/annotations/known_sites1.mini.vcf.gz   (SNPs)
  test/fixtures/refs/annotations/known_sites2.mini.vcf.gz   (indels)
"""
import gzip
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
REF = ROOT / "test/fixtures/refs/sequence/Homo_sapiens_assembly38_chr21.fasta"
OUT = ROOT / "test/fixtures/refs/annotations"
COMP = {"A": "G", "C": "T", "G": "A", "T": "C", "N": "A"}


def read_ref():
    seqs, name, parts = {}, None, []
    for line in REF.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(parts).upper()
            name = line[1:].split()[0]
            parts = []
        else:
            parts.append(line.strip())
    if name is not None:
        seqs[name] = "".join(parts).upper()
    return seqs


def main():
    seqs = read_ref()
    seq = seqs["chr21"]
    length = len(seq)

    def header(source):
        return (
            "##fileformat=VCFv4.2\n"
            f"##source={source}\n"
            f"##contig=<ID=chr21,length={length}>\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        )

    snp_lines, indel_lines = [], []
    for i, pos in enumerate(range(100, length, 100)):
        ref_base = seq[pos - 1]
        alt = COMP[ref_base]
        snp_lines.append(
            f"chr21\t{pos}\trsMINI_SNP{i + 1:03d}\t{ref_base}\t{alt}\t50\tPASS\t."
        )
        indel_lines.append(
            f"chr21\t{pos}\trsMINI_INDEL{i + 1:03d}\t{ref_base}\t{ref_base}{alt}\t50\tPASS\t."
        )

    with gzip.open(OUT / "known_sites1.mini.vcf.gz", "wt") as fh:
        fh.write(header("synthetic-mini-snps") + "\n".join(snp_lines) + "\n")
    with gzip.open(OUT / "known_sites2.mini.vcf.gz", "wt") as fh:
        fh.write(header("synthetic-mini-indels") + "\n".join(indel_lines) + "\n")

    print(f"wrote {len(snp_lines)} SNP + {len(indel_lines)} indel mini sites "
          f"for a {length} bp chr21 -> {OUT}")


if __name__ == "__main__":
    main()
