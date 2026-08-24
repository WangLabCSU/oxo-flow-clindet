# Port note: upstream sequenza_call references "../../../../scripts/sequenza.R"
# (workflow/WES/rules/rtm/paired/CNV/sequenza.smk) — that file does NOT exist
# anywhere in zyllifeworld/clindet @ 582a9131 (verified by tree enumeration),
# so the upstream rule is broken at runtime. This script implements the
# standard sequenza call step (sequenza.extract -> fit -> results, the same
# three-step API the upstream intended), writing the same
# {sample}_segments.txt output the upstream wrapper targets.
library(sequenza)

# Version-compat shim (live-verified on tx-ubuntu, r-sequenza 3.0.0 +
# iotools 0.3.5): sequenza's gc.sample.stats calls
#   chunk.apply(input=..., FUN=..., CH.MAX.SIZE=..., parallel = <n>)
# but the chunk.apply it resolves is iotools::chunk.apply, whose parallelism
# argument is named CH.PARALLEL. The stray `parallel` therefore lands in
# iotools' `...` and is forwarded into the chunk FUN, killing the run with
# "unused argument (parallel = N)". No fixed r-sequenza exists (CRAN latest
# is 3.0.0), so this shim shadows chunk.apply in the script's global
# environment (sequenza has no importFrom binding for it, so the global
# definition wins) and translates `parallel` to the installed iotools'
# argument name. Real-data runs take this identical path — it is a version
# bug, not a mini-fixture accommodation.
chunk.apply <- function(input, FUN, ..., parallel = NULL) {
    f <- iotools::chunk.apply
    args <- c(list(input = input, FUN = FUN), list(...))
    if (!is.null(parallel)) {
        nms <- names(formals(f))
        if ("parallel" %in% nms) {
            args$parallel <- parallel
        } else if ("CH.PARALLEL" %in% nms) {
            args$CH.PARALLEL <- parallel
        } else {
            stop("iotools::chunk.apply supports neither `parallel` nor ",
                 "`CH.PARALLEL`; cannot forward the parallelism argument")
        }
    }
    do.call(f, args)
}

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
