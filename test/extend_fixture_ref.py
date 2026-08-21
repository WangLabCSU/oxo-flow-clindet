#!/usr/bin/env python3
"""Extend the mini fixture reference with a second contig (chrX).

arriba filters intra-chromosomal events at ~300-700 bp distance as
deletions, so the single-contig 922 bp fixture can never produce a
reported fusion (observed live: exit 0 but zero rows in the fusion
TSV).  A second contig lets generate_rna_fusions.py build genuine
inter-chromosomal chimeric reads (chr21 donor -> chrX acceptor).

Regenerates, in place:
  - test/fixtures/refs/sequence/Homo_sapiens_assembly38_chr21.fasta
  - .../Homo_sapiens_assembly38_chr21.fasta.fai   (offsets recomputed)
  - .../Homo_sapiens_assembly38_chr21.dict        (SAM header rewritten)
  - test/fixtures/refs/annotations/mini_chr21.gtf (mini_gene3 on chrX)
  - test/fixtures/refs/annotations/arriba_protein_domains_mini.gff3
    (domain on mini_gene3)

The chrX sequence is the reverse of chr21, which keeps the two contigs
distinguishable but equally "genomic-looking".

Usage: python3 test/extend_fixture_ref.py
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SEQ_DIR = ROOT / "test/fixtures/refs/sequence"
ANN_DIR = ROOT / "test/fixtures/refs/annotations"
FASTA = SEQ_DIR / "Homo_sapiens_assembly38_chr21.fasta"

COMPL = str.maketrans("ACGTN", "TGCAN")


def read_fasta(path: Path) -> list[tuple[str, str]]:
    seqs, name = [], None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name:
                seqs.append((name, "".join(chunks)))
            name, chunks = line[1:].split()[0], []
        else:
            chunks.append(line.strip())
    if name:
        seqs.append((name, "".join(chunks)))
    return seqs


def main() -> None:
    seqs = read_fasta(FASTA)
    assert [n for n, _ in seqs] == ["chr21"], seqs
    chr21 = seqs[0][1]
    chrX = chr21.translate(COMPL)[::-1]  # reverse complement: distinct, valid

    # --- fasta + fai + dict -------------------------------------------------
    lines = []
    for n, s in (("chr21", chr21), ("chrX", chrX)):
        lines.append(f">{n}")
        lines.extend(s[i : i + 60] for i in range(0, len(s), 60))
    FASTA.write_text("\n".join(lines) + "\n")
    fai_lines, offset = [], 0
    for n, s in (("chr21", chr21), ("chrX", chrX)):
        fai_lines.append(f"{n}\t{len(s)}\t{offset}\t60\t61")
        offset += len(f">{n}\n") + len(s) + 1
    (SEQ_DIR / "Homo_sapiens_assembly38_chr21.fasta.fai").write_text(
        "\n".join(fai_lines) + "\n")
    (SEQ_DIR / "Homo_sapiens_assembly38_chr21.dict").write_text(
        "@HD\tVN:1.6\tSO:unsorted\n"
        f"@SQ\tSN:chr21\tLN:{len(chr21)}\tM5:00000000000000000000000000000000\n"
        f"@SQ\tSN:chrX\tLN:{len(chrX)}\tM5:00000000000000000000000000000000\n")

    # --- GTF: mini_gene3 on chrX --------------------------------------------
    gtf = ANN_DIR / "mini_chr21.gtf"
    text = gtf.read_text()
    assert "mini_gene3" not in text
    gtf.write_text(text + (
        'chrX\tmini\tgene\t50\t250\t.\t+\t.\tgene_id "mini_gene3"; '
        'gene_name "MINI3";\n'
        'chrX\tmini\ttranscript\t50\t250\t.\t+\t.\tgene_id "mini_gene3"; '
        'transcript_id "mini_gene3_t1"; gene_name "MINI3";\n'
        'chrX\tmini\texon\t50\t250\t.\t+\t.\tgene_id "mini_gene3"; '
        'transcript_id "mini_gene3_t1"; exon_number "1"; gene_name "MINI3";\n'))

    # --- protein domains: a domain on mini_gene3 -----------------------------
    pd = ANN_DIR / "arriba_protein_domains_mini.gff3"
    text = pd.read_text()
    assert "mini_gene3" not in text
    pd.write_text(text + (
        'chrX\tpfam\tprotein_domain\t50\t250\t0\t+\t.\t'
        'Name=mini_domain3;color=#00FF00;gene_id=mini_gene3;'
        'gene_name=MINI3;protein_domain_id=PF00003\n'))

    print("chr21:", len(chr21), "chrX:", len(chrX))
    print("fasta/fai/dict/gtf/protein_domains updated")


if __name__ == "__main__":
    main()
