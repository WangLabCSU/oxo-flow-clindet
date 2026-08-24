# Port note: upstream sequenza_call references "../../../../scripts/sequenza.R"
# (workflow/WES/rules/rtm/paired/CNV/sequenza.smk) — that file does NOT exist
# anywhere in zyllifeworld/clindet @ 582a9131 (verified by tree enumeration),
# so the upstream rule is broken at runtime. This script implements the
# standard sequenza call step (sequenza.extract -> fit -> results, the same
# three-step API the upstream intended), writing the same
# {sample}_segments.txt output the upstream wrapper targets.
library(sequenza)

bin_seqz <- snakemake@input[['bin_seqz']]
out_segment <- snakemake@output[['segment']]
wd <- snakemake@params[['wd']]

if (!dir.exists(wd)) dir.create(wd, recursive = TRUE)

seqz_data <- sequenza.extract(bin_seqz, verbose = FALSE)
fit <- sequenza.fit(seqz_data)
results <- sequenza.results(seqz_data, fit, out.dir = wd)

out_df <- results$segments
if (is.null(out_df)) {
  out_df <- data.frame(
    chromosome = character(), start.pos = integer(), end.pos = integer(),
    N.BAF = numeric(), sd.BAF = numeric(), depth.ratio = numeric(),
    sd.ratio = numeric(), CNt = integer(), L = numeric(), CNvalue = character(),
    p.value = numeric()
  )
}
write.table(out_df, file = out_segment, sep = "\t", quote = FALSE, row.names = FALSE)
