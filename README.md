# oxo-flow-clindet — Cancer genome & transcriptome analysis (WES/WGS/RNA): somatic+germline+CNV+SV calling, MAF annotation, case report

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

Requires **oxo-flow ≥ 0.16.0** — the first release with the
wildcard-scoped `when` predicates + per-instance pair binding this
single-entry form needs (Traitome/oxo-flow PR #187). Release binary
(recommended):

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
  path (upstream default: `resources/ref_genome/hg38/vep`). The vcf2maf
  rules and the downstream MAF merge/flag/cancer-report tail are gated on
  `vep_cache_ready = false` — set it `true` once the cache is in place
  (upstream fails hard without it).
- **Compute**: up to 30 threads / 10 GB per rule (mapping rules request
  30 threads; the cancer report 10 GB).
- **Tools**: conda envs in `envs/` (pinned: varscan 2.4.6, vcf2maf 1.6.22,
  ensembl-vep 114.2, libboost 1.85.0, strelka 2.9.7 — the 2.9.10 conda
  build's strelka2 binary segfaults in the loader on glibc 2.39, see
  `envs/strelka.yaml`) plus pinned Singularity/Docker images
  for GATK 4.6.2.0, VarDict 1.8.3, MuSE 2.1.2 and CaVEMan 1.15.3 (see the
  `[rules.environment]` entries).

## Usage

Single entry file, mirroring the upstream Snakefile: `[config] run_type`
selects the rule tree, and paired vs tumor-only WES derives **per pair**
from the sample sheet.

```bash
# 1. prepare data (see test/fixtures for the expected layout; fill in your
#    own reads/reference/annotation paths in main.oxoflow or override on the
#    command line)
# 2. preview the plan
oxo-flow dry-run main.oxoflow                  # run_type = wes (default)
oxo-flow dry-run main.oxoflow --arg run_type=wgs
oxo-flow dry-run main.oxoflow --arg run_type=rna
# 3. run
oxo-flow run main.oxoflow -j 8
```

### Run types (upstream VALID_RUN_TYPES)

| run_type | Tree |
|---|---|
| `wes` (default) | paired WES for pairs with a control; **tumor-only WES for pairs without** (seven portable callers: Mutect2, HaplotypeCaller, varscan2, Strelka+Manta, vardict, lofreq, freebayes) — upstream's sample-sheet-derived branches, one DAG |
| `wgs` | WGS metrics + the paired callers with WGS config (`rules/91_wgs_callers.oxoflow` — no `--intervals`/`--callRegions`/`--exome`, Manta `somaticSV`) + delly/svaba SV chains (`rules/90_wgs.oxoflow`) |
| `rna` | fastp → STAR (arriba + mutation maps) → arriba fusion → SplitNCigarReads → unpaired SNV callers → VEP-gated MAF tail (upstream default stages `arriba` + `call_mut`) |

Upstream also defines `build_b37`/`build_hg38` (reference builders) and
`pull_zenodo` (data pull) run types — reference-data utilities without a
mini fixture; not ported (documented). The same rule family covers the
setup `download_*` (x19), `build_bwa_index`/`build_star_index`/
`build_kallisto_salmon_index`/`build_rsem_reference`/`rsem_star_index`/
`kallisto_salmon_index`/`bwa_index`, `create_conda_env`, `mass_config`,
`prerequisites(build conda env)` and `install_clindet_extras` rules —
all zero-coverage here.

The samplesheet (`samplesheet.csv`) holds one row per pair:
`pair_id,experiment,control,experiment_type` — leave `control` empty for a
tumor-only pair. Requires oxo-flow ≥ 0.16.0 (wildcard-scoped `when`
predicates, `wildcard.<key>`, landed in v0.16.0 — Traitome/oxo-flow PR #187).

### Limitation: tumor FASTQ paths are shared across pairs

Upstream reads per-sample FASTQ paths from the sample sheet
(`Tumor_R1_file_path` etc.). The port keeps the FASTQs in `[config]`
(`tumor_fastq_r1/r2`, `normal_fastq_r1/r2`, `rna_fastq_r1/r2`), so **every
pair consumes the same literal read files** (`fastp_tumor_sample` in
`rules/00_common.oxoflow` and `fastp_tumor_sample_unpaired` in
`rules/70_unpaired.oxoflow` expand over the pairs list but substitute the
same config path). With a multi-pair sheet like the fixture (mini + mini2),
all pairs process identical reads. Per-pair reads require either one pair
per workflow invocation or replacing the config literals with per-pair
paths (sample-sheet-driven inputs are an oxo-flow engine feature, not
something a workflow file can express today).

### RNA default stages (targeted)

The upstream mini-test default stages (`arriba`, `call_mut`) map to explicit
rule targets; RSEM/kallisto/salmon/TRUST4 quant rules and isofox are not in
the default stages upstream and run only when targeted:

```bash
oxo-flow run main.oxoflow --arg run_type=rna -j 8 \
  -t fastp_trim -t STAR_1_pass -t STAR_arriba_map -t STAR_mut_map \
  -t arriba_fusion -t link_bam -t SplitNCigarReads \
  -t mutect2 -t M2_filter -t unpaired -t call_variants -t lofreq \
  -t varscan2 -t norm_filter

# quant extras (need rsem/kallisto/salmon indexes — empty in the fixture kit)
oxo-flow run main.oxoflow --arg run_type=rna -t cal_exp_RSEM -t kallisto -t salmon -t TRUST4_TBCR
```

## Non-human genomes (config parity)

Upstream ships `workflow/config/conf/genomes.yaml` entries for worm
(`WBcel235`) and mouse (`mm10`) in addition to `b37`/`hg38` — the rules are
species-agnostic (same mapping/calling chain, different reference files),
so any upstream genome works by overriding the `[config]` reference keys:

| Upstream genomes.yaml key | WBcel235 (worm) | mm10 (mouse) |
|---|---|---|
| `REFFA` | `WBcel235_genome.fa` | `Sanger/core_ref_mm10/genome.fa` |
| `GTF` | `Caenorhabditis_elegans.WBcel235.114.gtf` | *(empty upstream)* |
| `DBSNP` / `DBSNP_INDEL` / `MUTECT2_VCF` | worm fake-dbsnp VCFs | *(mostly empty upstream)* |

e.g. `oxo-flow run main.oxoflow --arg reference=…/WBcel235_genome.fa --arg
target_bed=… --arg dbsnp=…` (the paired WES chain is what upstream runs for
these genomes; the exome bed is the genome's capture). The mini fixture kit
does not ship non-human references.

## Fidelity

| Upstream process/rule | oxo-flow rule | Tool (version) | Notes |
|---|---|---|---|
| fastp (T/N) | `fastp_tumor_sample` / `fastp_normal_sample` | fastp | identical flags (`-w 8 -Q -c -L`) |
| bam flagstat (T/N) | `bam_flagstat_tumor` / `bam_flagstat_normal` | samtools | identical command |
| map_reads (T/N) | `map_reads_tumor` / `map_reads_normal` | bwa >=0.7.18, samtools | `bwa mem -MR` + `fixmate` + `sort` |
| mark_duplicates (T/N) | `mark_duplicates_tumor` / `mark_duplicates_normal` | gatk4 4.6.2.0 (container) | `MarkDuplicates --CREATE_INDEX true`, `VALIDATION_STRINGENCY SILENT` |
| recal_link (T/N) | `recal_link_tumor` / `recal_link_normal` | ln -s | `when = "!config.recal_bqsr"` (upstream mini-test default `recal_BQSR: False`) |
| recalibrate_base_qualities (T/N) | `recalibrate_base_qualities_tumor` / `recalibrate_base_qualities_normal` | gatk4 (container) | `BaseRecalibrator --use-original-qualities`, known-sites = upstream varanno KNOWN_SITES1/2; `when = "config.recal_bqsr"` |
| apply_base_quality_recalibration (T/N) | `apply_base_quality_recalibration_tumor` / `apply_base_quality_recalibration_normal` | gatk4 (container) | `ApplyBQSR -use-original-qualities` + `samtools index`; `when = "config.recal_bqsr"` |
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
| loop_vcf2maf_paired | `vcf2maf_{Mutect2,vardict,varscan2,muse,HaplotypeCaller}` | vcf2maf 1.6.22, ensembl-vep 114.2 | verbatim tumor/normal IDs per `get_vcf_name`; gated on `vep_cache_ready` |
| loop_vcf2maf_germ_paired | `vcf2maf_germ_strelkamanta` / `vcf2maf_germ_caveman` | vcf2maf 1.6.22 | TUMOUR/NORMAL for CaVEMan; gated on `vep_cache_ready` |
| merge_loop (somatic) | `merge_paired_maf` | merge_maf.R (verbatim) | driven via scripts/smk.R shim |
| merge_paired_vcf | `merge_paired_vcf` | merge_caller_vcfs.py (verbatim) + pysam | driven via scripts/smk.py shim; upstream mini-test default stage `call_mut_vcf` |
| merge_loop_germline | `merge_paired_germ_maf` | merge_maf.R (verbatim) | |
| make_region_bed_list + flag_mutation_pairead_maf | `make_region_bed_list` + `flag_mutation_pairead_maf` | flag_mutation_maf.R (verbatim) | empty bed_list = header-only TSV |
| run_cancer_report | `run_cancer_report` | R >=4.4 (knitr, gpgr via post-deploy) | only MAF/panel/Rmd params; CNV/QC params unset (NULL) as in upstream default path |
| combined_multiqc | `prep_multiqc_data` + `combined_multiqc_prep_multiqc_data` + `combined_multiqc` | multiqc | conpair/purple inputs out of scope |
| freec_config / freec_call_paired / plot_freec | `freec_config` / `freec_call_paired` / `plot_freec` | Control-FREEC >=11.6, sambamba | verbatim `config_freec.py` + `config_exome.ini`; upstream runs freec in the facets-suite container, port uses bioconda control-freec; `when = "config.cnv_enabled"` |
| sequenza bam2seqz/binning/call | `sequenza_bam2seqz` / `sequenza_seqz_binning` / `sequenza_call` | sequenza-utils, r-sequenza | upstream's referenced `scripts/sequenza.R` does not exist in the tree — port ships the standard extract→fit→results chain (`scripts/sequenza_call.R`) |
| CNA_exomedepth | `CNA_exomedepth` | ExomeDepth (Bioc) | verbatim `ExomeDepth.R`; upstream counts over hardcoded exons.hg19 (dead `target.file` read) — port keeps that and adds a documented `use_target_bed` switch for the mini fixture |
| CNA_ASCAT / ASCAT_EXTRACT_PURITYPLOIDY | `CNA_ASCAT` / `ASCAT_EXTRACT_PURITYPLOIDY` | ASCAT >=3.2, alleleCounter | verbatim `ASCAT.R` (+chroms/GC/rt from config — upstream hardcodes c(1:22)); `ascat_pp.R` verbatim |
| WGS mapping/recal/QC | shared `00_common` rules | bwa/gatk4/picard | upstream WGS map_reads/markdup/BQSR/recal_link are identical to WES |
| WGS callers (no exome restrictions) | `rules/91_wgs_callers.oxoflow` | gatk4/muse/varscan/vardict/strelka2/manta | no `--intervals`/`--callRegions`/`--exome`; Manta emits somaticSV; germline Strelka takes BOTH bams (WES: normal only); vardict regions from `vardict_wgs_bed` |
| WGS picard_collect_wgs / picard_flength_wgs | `picard_collect_wgs_{tumor,normal}` / `picard_flength_wgs_{tumor,normal}` | picard | CollectWgsMetrics + CollectInsertSizeMetrics (telomerecat is "departed" upstream — not ported) |
| SV_delly chain | `SV_delly` → `SV_delly_sample_tsv` → `SV_delly_filter_somatic` → `SV_delly_to_vcf` → `delly_filter` → `delly2bnd` | delly 1.7.2 (container), bcftools | verbatim; `delly2bnd.py` verbatim (upstream env lacks cyvcf2 — added to envs/clindet.yaml, upstream bug) |
| SV_svaba | `SV_svaba` | svaba (container) | verbatim run |
| sansa-annotation (svaba + delly) | `SV_sansa_anno_svaba` / `SV_sansa_annodelly` | site-provided sansa binary (`sansa_call`) | gated on `sansa_db`/`sansa_g` being set (upstream gates on the sansa software config being present — absent upstream by default, so off by default here too; zero instances without the keys) |
| svanno (gatk SVAnnotate) | `SV_svanno_svaba` | gatk4 (container) | gated on `svanno_gtf` (protein-coding GTF); zero instances without the key |
| Manta SV | `call_config_strelka` (WGS) | manta | `somaticSV.vcf.gz` from the WGS Manta run (upstream SV list entry 'Manta') |

**Not ported** (upstream branches with reasons):
- CNV purple/amber/cobalt: HMF tools run in the upstream's custom
  hmftools.sif with the multi-GB hmf_pipeline_resources tree (built
  locally upstream, `pull_zenodo` run type) — not a portable image;
  dryclean has no rule file upstream (list-only).
- CNV FACETS/facets-suite: custom facets-suite-dev.img + snp-pileup PoN
  chain, requires compiling cnv_facets C++.
- CNV CNA_ABSOLUTE_GISTIC / ASCAT_GISTIC: ABSOLUTE + GISTIC2 are Broad
  tools without conda packages; SM_check / CNA_Battenberg are marked
  "for future development" upstream (Battenberg needs cgpbattenberg371.sif
  + 1000G impute reference data).
- SV gridss/BRASS/linx/igv-caller/jasmine: custom containers (gridss2/
  brass634/jasminesv sifs) + Sanger VAGrENT/BRASS and HMF resource trees.
- WGS unpaired callers (sage/deepvariant/pindel/octopus/UnifiedGeniTyper)
  and WGS Battenberg/ecDNA/VirusScan: custom containers + resource trees
  as above.
- conpair contamination check: custom conpair_latest.sif container.
- unpaired CNV (upstream `rtm/unpaired/CNV.smk`: freec + purple): the
  ported CNV branch is paired-only (upstream mini-test default
  `somatic_cnv_list: [notrun]`).
- non-human genomes: supported at config level (WBcel235/mm10 parity table
  above) — the upstream rule set itself is species-agnostic.

## Source

Ported from **[zyllifeworld/clindet](https://github.com/zyllifeworld/clindet)**,
version `582a9131` (MIT). Created 2026-08-15; this workflow **may lag
upstream releases**. Attribution in `NOTICE.md`; upstream license in
`LICENSE.upstream`.

#
## Not ported (per-rule, upstream family list)

- **Mutect2 PoN build chain (M2_CSPN/pon_GB/pon_facetsCH/call_variants_pon)** (`workflow/WGS/rules/rtm/paired/SNV/Mutect2_pon.smk`/`workflow/WES/rules/rtm/paired/SNV/Mutect2_pon.smk`): `M2_CSPN`, `M2_SNC`, `M2_ST`, `M2_contam`, `M2_filter`, `call_variants_pon`, `mutect2`, `pon_GB`
- **BIC-seq2** (`workflow/WES/rules/rtm/paired/CNV/bicseq2.smk`/`workflow/WGS/rules/rtm/paired/CNV/bicseq.smk`): `bicseq`, `bicseq2_norm`, `bicseq2_samtools_normal`, `bicseq2_samtools_tumor`, `bicseq2_seg`
- **esvee** (`workflow/WGS/rules/rtm/paired/SV/esvee.smk`): `esvee_prep`
- **lumpy** (`workflow/WGS/rules/rtm/paired/SV/lumpy.smk`): `extract_lumpy_evidence`, `lumpy_call_paired`, `lumpy_svtyper_paired`
- **sentieon (WGS paired)** (`workflow/WGS/rules/rtm/paired/SNV/sentieon.smk`): `call_variants_sentieon`, `filter_sentieon`
- **lancet2** (`workflow/WES/rules/rtm/paired/SNV/Lancet2.smk`): `lancet2_somatic_call`
- **moalmanac + split_maf_snp_indel** (`workflow/common/rules/case_report.smk`): `moalmanac_annotation`, `split_maf_snp_indel`
- **orange / bamMetrics (WGS report)** (`workflow/WGS/rules/report/orange.smk`): `bamMetrics_tumor`, `orange`
- **ngs_bit QC** (`workflow/common/rules/qc.smk`): `ngs_bit_mapping`, `ngs_bit_sample_gender`
- **summary_softwares** (`workflow/common/rules/summary_softwares.smk`): `clindet_hmftool_version`, `clindet_main_version`, `clindet_rsem_version`, `clindet_vep_version`
- **setup reference-data family** (`workflow/setup/*`): `build_b37_ref`, `build_bwa_index`, `build_hg38_ref`, `build_kallisto_salmon_index`, `build_rsem_reference`, `build_star_index`, `bwa_index`, `create_conda_env`, `download_ascat`, `download_ascat_refdata`, `download_b37_gatk`, `download_b37_hmftools`, `download_b37_reference`, `download_cdna`, `download_cdna_fasta`, `download_ensembl_gff3`, `download_gatk_resources`, `download_gatk_software`, `download_grch37_gff3`, `download_grch37_gtf`, `download_gtf`, `download_hg38_genome`, `download_hmftools`, `download_mutation_anno_bed_b37`, `download_mutation_anno_bed_hg38`, `download_rna_edit_vcf`, `download_sanger`, `download_sanger_refdata`, `download_vep_cache`, `download_vep_cache_hg38`, `download_zenodo_containers`, `install_clindet_extras`, `kallisto_salmon_index`, `mass_config`, `prerequisites`, `rsem_star_index`, `star_index`

# Test

```bash
bash test/run.sh      # DNA: validate + lint + dry-run, exits 0
bash test/run_rna.sh  # RNA: validate + lint + dry-run (default stages), exits 0
```

## Live verification (tx-ubuntu, oxo-flow 0.14.1 / PR #187 engine)

| Run type | Status | Notes |
|---|---|---|
| wes paired (control set) | ✅ live-verified | full pipeline on the mini fixture, incl. opt-in BQSR (`recal_bqsr=true`) and the CNV subset (`cnv_enabled=true`: Control-FREEC, Sequenza, ExomeDepth, ASCAT) |
| wes tumor-only (control empty) | ✅ live-verified ×2 | 7 tumor-only callers + merge |
| wgs | ✅ live-verified | WGS metrics, delly SV chain (incl. germ), svaba, Manta somaticSV, MuSE/Strelka2 WGS configs |
| rna | ✅ live-verified | arriba/TRUST4/isofox default stages |
| **merged single-entry: wes** | ✅ live-verified (PR #187) | one DAG, both sheet rows — paired mini + tumor-only mini2, per-pair morphing + per-instance expand_inputs binding |
| **merged single-entry: wgs** | ✅ live-verified (PR #187) | paired-only WGS tree (unpaired rows do not instantiate in WGS, upstream semantics) |
| **merged single-entry: rna** | ✅ live-verified (PR #187) | default stages via prefix `-t` targets; docker-hub arriba image pulled through the box mirror (`OXO_REGISTRY_MIRRORS`) |

The first four entries were live-verified as separate entry files; the
merged single-entry rows were live-verified on tx-ubuntu with the PR #187
engine (wildcard-scoped `when` + per-instance pair binding — shipped as
oxo-flow v0.16.0). The vcf2maf/
VEP chain is gated on `vep_cache_ready` and was live-verified in the
4-entry era; the merged-form pass covers the same rules (identical
commands) once the VEP cache is present.

Mini-fixture degeneracy (900 bp chr21 + 20 kb chrX) is handled with
documented fallbacks where a tool needs real signal (ExomeDepth/ASCAT/
Control-FREEC/Sequenza fit write header-only outputs + a provenance note);
real-data runs take the verbatim upstream path. See the Fidelity section.

## License

Apache-2.0. Copyright (c) 2026 oxo-flow-community.

## Community

https://oxo-flow-community.github.io/

### Remaining upstream rule names by file (all inside the excluded families above; a few have renamed or same-name ported counterparts — see the Fidelity table)

- **workflow/WES/rules/cgp/cgpBattenberg.smk** (1): `battenberg_call`
- **workflow/WES/rules/rtm/paired/CNV/purple.smk** (1): `paired_purple`
- **workflow/WES/rules/rtm/paired/SNV/octopus.smk** (3): `octopus_call_paired`, `octopus_paired_germline`, `octopus_paired_somatic`
- **workflow/WES/rules/rtm/paired/SNV/sage.smk** (1): `sage_germline_filter_pass`
- **workflow/WES/rules/rtm/unpaired/CNV/purple.smk** (3): `unpaired_amber`, `unpaired_cobalt`, `unpaired_purple`
- **workflow/WES/rules/rtm/unpaired/SNV/octopus.smk** (1): `octopus_call_tumor_only`
- **workflow/WES/rules/rtm/unpaired/SNV/pindel.smk** (4): `PI_ggz_unpaired`, `build_unpaired_fake_bam`, `cgppindel_filter_somatic_unpaired`, `unpaired_PI_call`
- **workflow/WGS/rules/SV.smk** (9): `SV_brass`, `SV_brass_bamstat`, `SV_gridss`, `SV_sansa_merge`, `brass_ascat`, `brass_pre_merge`, `delly_pre_merge`, `jasmine_merge`, `manta_pre_merge`
- **workflow/WGS/rules/calling_unpaired.smk** (6): `unpaired_call_strelka_somatic`, `unpaired_call_variants_mutect2`, `unpaired_coverageBed`, `unpaired_coverageBed_2`, `unpaired_merge_strelka`, `unpaired_merge_strelka_somatic`
- **workflow/WGS/rules/cgp/cgppindel.smk** (1): `cgppindel_filter_somatic`
- **workflow/WGS/rules/merge_loop_germline.smk** (2): `loop_vcf2maf_germ_unpaired`, `merge_unpaired_germ_maf`
- **workflow/WGS/rules/rtm/paired/CNV/ASCAT.smk** (1): `CNA_ASCAT_sc`
- **workflow/WGS/rules/rtm/paired/CNV/Battenberg.smk** (3): `CNA_Battenberg_ABSOLUTE_GISTIC`, `CNA_Battenberg_combine`, `CNA_Battenberg_v2`
- **workflow/WGS/rules/rtm/paired/CNV/ecDNA.smk** (1): `ampliconsuite`
- **workflow/WGS/rules/rtm/paired/CNV/facets.smk** (2): `CNA_snp_pileup_nor`, `CNA_snp_pileup_tum`
- **workflow/WGS/rules/rtm/paired/CNV/facets_suite.smk** (4): `facets_annotate_maf`, `facets_calling`, `facets_gene_cna`, `facets_pileup`
- **workflow/WGS/rules/rtm/paired/CNV/purple.smk** (2): `paired_amber`, `paired_cobalt`
- **workflow/WGS/rules/rtm/paired/CNV/sequenza.smk** (1): `sequenza_gc_bins`
- **workflow/WGS/rules/rtm/paired/SNV/DeepVariant.smk** (3): `deepvariant_filter_germline`, `deepvariant_filter_somatic`, `deepvariant_somatic_call`
- **workflow/WGS/rules/rtm/paired/SNV/UnifiedGenoTyper.smk** (1): `call_variants_UnifiedGenoTyper`
- **workflow/WGS/rules/rtm/paired/SNV/sage.smk** (2): `paired_sage_germline`, `pave_anno_sage_germline`
- **workflow/WGS/rules/rtm/paired/SNV/vardict.smk** (1): `vardict_filter_germline`
- **workflow/WGS/rules/rtm/paired/SV/BRASS.smk** (1): `brass_cnv`
- **workflow/WGS/rules/rtm/paired/SV/gridss.smk** (2): `SV_gridss_filter`, `gridss_rename_tumor`
- **workflow/WGS/rules/rtm/paired/SV/jasmine_merge.smk** (2): `gridss_pre_merge`, `svaba_pre_merge`
- **workflow/WGS/rules/rtm/paired/SV/linx.smk** (2): `paired_linx`, `report_linx`
- **workflow/WGS/rules/rtm/paired/SV/svaba.smk** (1): `svaba_rename_tumor`
- **workflow/WGS/rules/rtm/paired/VirusScan.smk** (1): `virusbreakend`
- **workflow/WGS/rules/rtm/unpaired/CNV.smk** (1): `unpair_CNA_ASCAT_sc`
- **workflow/WGS/rules/rtm/unpaired/CNV/freec.smk** (3): `freec_call_unpaired`, `plot_freec_unpaired`, `unpair_freec_config`
- **workflow/WGS/rules/rtm/unpaired/DeepVariant.smk** (1): `deepvariant_call`
- **workflow/WGS/rules/rtm/unpaired/SNV/deepvariant.smk** (2): `deepvariant_norm`, `unpaird_deepvariant_call`
- **workflow/WGS/rules/rtm/unpaired/SV.smk** (2): `SV_unp_delly`, `delly_unp_filter`
- **workflow/WGS/rules/rtm/unpaired/UnifiedGeniTyper.smk** (1): `unpaired_call_variants_UnifiedGenoTyper`
- **workflow/common/rules/merge_loop.smk** (1): `paired_maf_report`
- **wrapper/rna.smk** (1): `rna_pipeline`
- **wrapper/wes.smk** (1): `wes_pipeline`
- **wrapper/wgs.smk** (1): `wgs_pipeline`
