# oxo-flow-clindet — Clinical WES tumor/normal variant calling

[![CI](https://github.com/WangLabCSU/oxo-flow-clindet/actions/workflows/ci.yml/badge.svg)](https://github.com/WangLabCSU/oxo-flow-clindet/actions/workflows/ci.yml)

> ☆ Community · ⇄ Official port of [`zyllifeworld/clindet`](https://github.com/zyllifeworld/clindet) @ `582a9131`
> — same tools, same versions, same commands. Part of the
> [oxo-flow-community catalog](https://oxo-flow-community.github.io/).

Clinical WES tumor/normal pipeline: fastp QC → bwa+fixmate+samtools sort →
GATK MarkDuplicates → five SNV callers (Mutect2, VarDict, VarScan2, MuSE,
HaplotypeCaller) plus germline Strelka2+Manta and CaVEMan → bcftools
normalization → vcf2maf with VEP annotation → merged MAF → region-based
mutation flagging → cancer case report (Rmd/knitr) and MultiQC.

## Installation

### 1. Install oxo-flow

Requires **oxo-flow ≥ 0.12.0**. Release binary (recommended):

```bash
curl -fL -o oxo-flow.tar.gz \
  https://github.com/Traitome/oxo-flow/releases/latest/download/oxo-flow-latest-x86_64-unknown-linux-gnu.tar.gz
tar xzf oxo-flow.tar.gz
sudo mv oxo-flow /usr/local/bin/
```

Alternative: `conda install -c bioconda oxo-flow-cli` (the bioconda package
may lag behind releases; other platform binaries are on the releases page).

### 2. Get this workflow

```bash
git clone https://github.com/WangLabCSU/oxo-flow-clindet.git
```

### 3. Requirements

- **Reference data** (paths in `[config]`, edit `main.oxoflow`): genome FASTA
  with `.fai` and `.dict`, target BED, dbSNP / Mills gold-standard indels /
  gnomAD VCFs, and a BWA index of the reference FASTA; annotation VCFs need
  tabix indexes (runtime, when running the pipeline for real).
- **VEP cache**: GRCh38, cache version 110, at the configured `vep_data`
  path (upstream default: `resources/ref_genome/hg38/vep`).
- **Compute**: up to 30 threads / 10 GB per rule (mapping rules request
  30 threads; the cancer report 10 GB).
- **Tools**: conda envs in `envs/` (pinned: varscan 2.4.6, vcf2maf 1.6.22,
  ensembl-vep 114.2, libboost 1.85.0, strelka 2.9.7 — the 2.9.10 conda
  build's strelka2 binary segfaults in the loader on glibc 2.39, see
  `envs/strelka.yaml`) plus pinned Singularity/Docker images
  for GATK 4.6.2.0, VarDict 1.8.3, MuSE 2.1.2 and CaVEMan 1.15.3 (see the
  `[rules.environment]` entries).

## Usage

```bash
# 1. prepare data (see test/fixtures for the expected layout; fill in your
#    own reads/reference/annotation paths in main.oxoflow or override on the
#    command line)
# 2. preview the plan
oxo-flow dry-run main.oxoflow
# 3. run
oxo-flow run main.oxoflow -j 8
```

The samplesheet (`samplesheet.csv`) holds one `pair_id` per tumor/normal
pair: columns `pair_id,experiment,control,experiment_type`.

## Fidelity

| Upstream process/rule | oxo-flow rule | Tool (version) | Notes |
|---|---|---|---|
| fastp (T/N) | `fastp_tumor_sample` / `fastp_normal_sample` | fastp | identical flags (`-w 8 -Q -c -L`) |
| bam flagstat (T/N) | `bam_flagstat_tumor` / `bam_flagstat_normal` | samtools | identical command |
| map_reads (T/N) | `map_reads_tumor` / `map_reads_normal` | bwa >=0.7.18, samtools | `bwa mem -MR` + `fixmate` + `sort` |
| mark_duplicates (T/N) | `mark_duplicates_tumor` / `mark_duplicates_normal` | gatk4 4.6.2.0 (container) | `MarkDuplicates --CREATE_INDEX true`, `VALIDATION_STRINGENCY SILENT` |
| recal_link (T/N) | `recal_link_tumor` / `recal_link_normal` | ln -s | BQSR off (upstream mini-test default) |
| bed_to_interval_list | `bed_to_interval_list` | gatk4 (container) | `BedToIntervalList --SORT true` |
| picard_collect_wes (T/N) | `picard_collect_wes_tumor` / `picard_collect_wes_normal` | gatk4 (container) | `CollectHsMetrics` |
| call_variants_HaplotypeCaller | `call_variants_HaplotypeCaller` | gatk4 (container) | identical `-A` annotation flags |
| vardict_paired_mode | `vardict_paired_mode` | vardict-java 1.8.3 (container) | `vardict-java` + `testsomatic.R` + `var2vcf_paired.pl` |
| vardict_filter_somatic | `vardict_filter_somatic` | bcftools >=1.22 | StrongSomatic/LikelySomatic + `SSF <= 0.05` |
| varscan2 mpileup/call/processSomatic/filter | `varscan2_mpileup` … `merge_somatic` | varscan 2.4.6, samtools, bcftools | `--strand-filter 1`, `--output-vcf 1`, concat chain |
| muse_call / muse_sump | `muse_call` / `muse_sump` | muse 2.1.2 (container) | `sump -E -n {threads} -D {dbsnp_gz}` |
| mutect2 chain | `M2_ST`/`M2_SNC`/`M2_contam`/`mutect2`/`M2_filter` | gatk4 (container) | `-pon` only when `wes_pon` set (upstream hg38_chr21 has none) |
| call_config_strelka | `call_config_strelka` | manta | `configManta.py --exome --callRegions` |
| call_strelka_manta_germline | `call_strelka_manta_germline` | strelka2, manta | `configureStrelkaGermlineWorkflow` + `runWorkflow.py` |
| merge_strelka_manta | `merge_strelka_manta` | bcftools | upstream `{params.indel}` bug fixed (single concat+sort) |
| strelka somatic (via manta) | `call_strelka_somatic_manta` + `merge_strelka_somatic_manta` | strelka2, manta, bcftools | same pipeline on somatic config |
| CM_cnv | `CM_cnv` | touch | empty tumour/normal CNV beds (upstream default) |
| CM_call / CM_flag | `CM_call` / `CM_flag` | caveman 1.15.3 (container) | `-td 2 -nd 2 -seqType WGS -no-flagging`, flag with `-umv .`; `-ignore-file` fed a one-region bed on a contig absent from the reference (upstream passes `""`, which caveman 1.15.3 rejects — same "no ignore regions" semantics); flagger gets real `-c`/`-v` configs from `test/fixtures/flag` (GRCh38 params verbatim, bed-based flags dropped — no chr21 flag data) plus empty `-b`/`-ab` dirs and `-t genomic` (upstream's `""`/`""`/`"genome"` are rejected by cgpFlagCaVEMan 1.15.3) |
| CM_germ_flag | `CM_germ_flag` | bcftools | `-e 'DP<=30' -s LowDP --mode x` |
| vcf_norm (per caller) | `vcf_norm_{Mutect2,vardict,varscan2,muse,HaplotypeCaller,germline_strelkamanta,germline_caveman}` | bcftools >=1.22 | verbatim per-caller FILTER rules incl. vardict contig-header branch |
| loop_vcf2maf_paired | `vcf2maf_{Mutect2,vardict,varscan2,muse,HaplotypeCaller}` | vcf2maf 1.6.22, ensembl-vep 114.2 | verbatim tumor/normal IDs per `get_vcf_name` |
| loop_vcf2maf_germ_paired | `vcf2maf_germ_strelkamanta` / `vcf2maf_germ_caveman` | vcf2maf 1.6.22 | TUMOUR/NORMAL for CaVEMan |
| merge_loop (somatic) | `merge_paired_maf` | merge_maf.R (verbatim) | driven via scripts/smk.R shim |
| merge_loop_germline | `merge_paired_germ_maf` | merge_maf.R (verbatim) | |
| make_region_bed_list + flag_mutation_pairead_maf | `make_region_bed_list` + `flag_mutation_pairead_maf` | flag_mutation_maf.R (verbatim) | empty bed_list = header-only TSV |
| run_cancer_report | `run_cancer_report` | R >=4.4 (knitr, gpgr via post-deploy) | only MAF/panel/Rmd params; CNV/QC params unset (NULL) as in upstream default path |
| combined_multiqc | `prep_multiqc_data` + `combined_multiqc_prep_multiqc_data` + `combined_multiqc` | multiqc | conpair/purple inputs out of scope |

**Not ported** (upstream branches outside the default paired WES path, with
reasons): CNV branch (purple/amber/cobalt/ASCAT/FACETS/sequenza/freec/
exomedepth — Broad-only, unbuildable without commercial licenses and heavy
reference data), extended SV branch (delly/gridss/svaba/BRASS/linx/
igv-caller/jasmine — reference-data heavy), RNA branch (arriba/TRUST4/
isofox/RNA SNV — separate sample type), unpaired (single-sample) mode.

## Source

Ported from **[zyllifeworld/clindet](https://github.com/zyllifeworld/clindet)**,
version `582a9131` (MIT). Created 2026-08-15; this workflow **may lag
upstream releases**. Attribution in `NOTICE.md`; upstream license in
`LICENSE.upstream`.

## Test

```bash
bash test/run.sh   # validate + lint + dry-run, exits 0
```

## License

Apache-2.0. Copyright (c) 2026 oxo-flow-community.

## Community

https://oxo-flow-community.github.io/
