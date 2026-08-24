#!/usr/bin/env python3
"""Add reference contigs missing from a BAM header's @SQ list.

The fixture reference gained a chrX contig with the RNA fusion fixture
(test/extend_fixture_ref.py), while the DNA mini reads are chr21-only.
Strict tools (Manta: "Reference genome mismatch: Tumor BAM/CRAM file is
missing a chromosome found in the reference fasta file") reject such BAMs.
This helper re-writes the BAM with the full reference contig list in @SQ
(zero-coverage contigs are legal in BAM).

Usage: python3 scripts/patch_bam_contigs.py IN.bam OUT.bam REF.fai
"""
import sys

import pysam


def main():
    if len(sys.argv) != 4:
        sys.exit("usage: patch_bam_contigs.py IN.bam OUT.bam REF.fai")
    in_bam, out_bam, fai = sys.argv[1:4]

    ref_contigs = []
    with open(fai) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                ref_contigs.append((fields[0], int(fields[1])))

    bam = pysam.AlignmentFile(in_bam, "rb")
    header = bam.header.to_dict()
    have = {c["SN"] for c in header.get("SQ", [])}
    for name, length in ref_contigs:
        if name not in have:
            header["SQ"].append({"SN": name, "LN": length})

    out = pysam.AlignmentFile(out_bam, "wb", header=header)
    for rec in bam:
        out.write(rec)
    out.close()
    bam.close()


if __name__ == "__main__":
    main()
