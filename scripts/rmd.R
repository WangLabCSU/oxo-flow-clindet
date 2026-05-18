#' Generate UMCCR Cancer Report
#'
#' Generates a UMCCR Cancer Report. It does so with the following steps:
#' 1. move the img_dir into 'tmp/img_dir'
#' 2. copy the rmd into 'tmp/cancer_report.Rmd'
#' 3. render the rmd inside 'tmp/'
#' 4. return the path to the output HTML
#'
# af_global Path to `af_tumor.txt` file.
# af_keygenes Path to `af_tumor_keygenes.txt` file.
# batch_name Name of batch sample.
# conda_list Path to `conda_pkg_list.txt` file.
# img_dir Path to directory containing PURPLE plots.
# key_genes Path to UMCCR cancer gene file.
# oncokb_genes Path to OncoKB database file.
# virusbreakend_tsv Path to VIRUSBreakend summary file.
# virusbreakend_vcf Path to VIRUSBreakend VCF file.
# purple_purity Path to `purple.purity.tsv`.
# purple_qc Path to `purple.qc`.
# purple_som_cnv Path to `purple.cnv.somatic.tsv`.
# purple_som_cnv_ann Path to annotated and prioritised `purple.cnv.somatic.tsv`.
# purple_som_gene_cnv Path to `purple.cnv.gene.tsv`.
# purple_som_snv_vcf Path to `purple.somatic.vcf.gz`.
# result_outdir Path to directory to write tidy JSON/TSV results.
# somatic_snv_vcf Path to `somatic-PASS.vcf.gz` SNV VCF.
# somatic_snv_summary Path to `somatic_snv_summary.json` JSON.
# somatic_sv_tsv Path to SV TSV file.
# somatic_sv_vcf Path to SV VCF file.
# tumor_name Name of tumor sample.
# out_file Path to output HTML file (needs '.html' suffix) (def: `{tumor_name}_cancer_report.html`).
# quiet Suppress log printing during rendering.

report_rmd <- snakemake@params[['report_rmd']]
output_file <- snakemake@params[['output_file']]
out_dir <- dirname(output_file)
result_outdir <- dirname(output_file)
sv_available = FALSE # will not report SVs in WES and target panel
pars <- list(
  batch_name = snakemake@params[['batch_name']],
  af_keygenes = snakemake@params[['af_keygenes']],
  genome_version = snakemake@params[['genome_version']],
  sv_available = sv_available,
  # bcftools_stats = snakemake@params[['bcftools_stats']],
  # key_genes = snakemake@params[['key_genes']],
  # somatic_snv_vcf = snakemake@params[['somatic_snv_vcf']],
  # somatic_snv_summary = snakemake@params[['somatic_snv_summary']],
  # somatic_sv_tsv = snakemake@params[['somatic_sv_tsv']],
  # somatic_sv_vcf = snakemake@params[['somatic_sv_vcf']],
  ### CNV results and plot dir
  ascat_res_dir = snakemake@params[['ascat_res_dir']],
  purple_res_dir = snakemake@params[['purple_res_dir']],
  ##
  purple_som_gene_cnv = snakemake@params[['purple_som_gene_cnv']],
  # purple_som_cnv_ann = snakemake@params[['purple_som_cnv_ann']],
  purple_som_cnv = snakemake@params[['purple_som_cnv']],
  purple_purity = snakemake@params[['purple_purity']],
  purple_qc = snakemake@params[['purple_qc']],
  purple_som_snv_vcf = snakemake@params[['purple_som_snv_vcf']],
  virusbreakend_tsv = snakemake@params[['virusbreakend_tsv']],
  virusbreakend_vcf = snakemake@params[['virusbreakend_vcf']],
  result_outdir = result_outdir,
  tumor_name = snakemake@params[['tumor_name']]
)

# suppress DT large size warning
quiet = FALSE
options(DT.warn.size = FALSE)
rmarkdown::render(
  input = report_rmd,
  params = pars,
  output_dir = out_dir,
  output_file = I(output_file),
  run_pandoc = TRUE,
  quiet = quiet
)
