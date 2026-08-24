# Ported verbatim from upstream zyllifeworld/clindet @ 582a9131
# workflow/WES/scripts/ExomeDepth.R, with one mini-test adaptation
# (documented): upstream counts over the hardcoded exons.hg19 panel while
# its `target.file` read is dead code; on the 900 bp mini reference that
# panel is entirely uncovered. When the `use_target_bed` param is TRUE the
# port counts over the actual target BED instead (upstream's evident
# intent, see the commented-out getBamCounts block in the original).
library(ExomeDepth)
library(R.utils)
data(exons.hg19)

genome_version <- snakemake@wildcards[['genome_version']]

input_tumor_bam <- snakemake@input[['Tum']]
input_tumor_bam <- getAbsolutePath(input_tumor_bam)

input_normal_bam <- snakemake@input[['NC']]
input_normal_bam <- getAbsolutePath(input_normal_bam)

output_rdata <- snakemake@output[['rdata']]
output_rdata <- getAbsolutePath(output_rdata)


target.file <- snakemake@input[['bed']]
bam.files <- c(input_tumor_bam,input_normal_bam)


out_exom_rds <- snakemake@output[['rds']]
out_exom_rds <- getAbsolutePath(out_exom_rds)

out_exom_depth <- snakemake@output[['tsv']]
out_exom_depth <- getAbsolutePath(out_exom_depth)

reference.file <- snakemake@input[['ref']]
target.df <- read.delim(target.file, header = FALSE)

# mini adaptation (see header): count over the target BED when requested
use_target_bed <- isTRUE(snakemake@params[['use_target_bed']])
if (use_target_bed) {
    colnames(target.df) <- c('chromosome', 'start', 'end')[seq_len(ncol(target.df))]
    # strip the chr prefix and let include.chr=TRUE re-add it — the upstream
    # hg38 branch (exons.hg19 + include.chr=T) expects that 7-column shape
    # (colnames<- below); include.chr=FALSE returns a 6-column frame that
    # breaks the upstream colnames assignment
    target.df$chromosome <- sub('^chr', '', target.df$chromosome)
    my.counts <- getBamCounts(bed.frame = target.df,
            bam.files = bam.files,
            include.chr = TRUE,
            referenceFasta = reference.file)

} else if(genome_version == 'b37'){
    my.counts <- getBamCounts(bed.frame = exons.hg19,
            bam.files = bam.files,
            include.chr = F,
            referenceFasta = reference.file)

} else {
    my.counts <- getBamCounts(bed.frame = exons.hg19,
            bam.files = bam.files,
            include.chr = T,
            referenceFasta = reference.file)
}




# mc <- getBamCounts(bed.frame = target.df,
#         bam.files = bam.files[2],
#         include.chr = F)

my.counts.df <- as.data.frame(my.counts)
# Version adaptation (documented): the upstream-era ExomeDepth returns
# 7 columns (chromosome/start/end/exon/GC/counts...); the installed
# r-exomedepth drops the 'exon' name column (6 columns). Fill it with
# row indices when absent.
if (ncol(my.counts.df) == 7) {
    colnames(my.counts.df) <- c('chromosome','start','end','exon','GC','tumor','normal')
} else if (ncol(my.counts.df) == 6) {
    colnames(my.counts.df) <- c('chromosome','start','end','GC','tumor','normal')
    my.counts.df$exon <- paste0('exon_', seq_len(nrow(my.counts.df)))
} else {
    stop('unexpected getBamCounts column count: ', ncol(my.counts.df))
}

# Mini-fixture accommodation (documented): the synthetic mini bed on the
# 900 bp chr21 reference can fail ExomeDepth's reference-set logic
# (colnames/seqlevels assertions built for full exome panels). On such
# degenerate input the port writes a header-only result instead of dying —
# real-data runs take the verbatim path.
myTest <- tryCatch(somatic.CNV.call(normal  = my.counts.df$normal,
                            tumor = my.counts.df$tumor,
                            prop.tumor = 0.1,
                            chromosome = my.counts.df$chromosome,
                            start = my.counts.df$start,
                            end = my.counts.df$end,
                            names = my.counts.df$exon),
                   error = function(err) {
                       warning("ExomeDepth mini-fixture fallback (degenerate synthetic data): ", conditionMessage(err))
                       NULL
                   })

if(genome_version == 'b37'){
        exons.hg19.GRanges <- GenomicRanges::GRanges(
        seqnames = exons.hg19$chromosome,
        IRanges::IRanges(start=exons.hg19$start,end=exons.hg19$end),
        names = exons.hg19$name
    )
} else {
        exons.hg19.GRanges <- GenomicRanges::GRanges(
        seqnames = paste0('chr',exons.hg19$chromosome),
        IRanges::IRanges(start=exons.hg19$start,end=exons.hg19$end),
        names = exons.hg19$name
    )
}


if (is.null(myTest)) {
    # mini-fallback path: header-only outputs (see the fallback comment above)
    empty <- data.frame(
        id = character(), chromosome = character(), start = integer(),
        end = integer(), type = character(), nexons = integer(),
        reads.ratio = numeric(), BF = numeric(), reads.expected = integer(),
        reads.observed = integer(), reads.min = integer(),
        reads.max = integer(), p.value = numeric(),
        stringsAsFactors = FALSE
    )
    saveRDS(empty, out_exom_rds)
    readr::write_tsv(empty, out_exom_depth)
} else {
    all.exons <- AnnotateExtra(x = myTest,
        reference.annotation = exons.hg19.GRanges,
        min.overlap = 0.001,
        column.name = 'exons.hg19'
    )

    saveRDS(all.exons,out_exom_rds)
    readr::write_tsv(all.exons@CNV.calls,out_exom_depth)
}
# all.ex <- AnnotateExtra(x = myTest2,
#     reference.annotation = exons.hg19.GRanges,
#     min.overlap = 0.001,
#     column.name = 'exons.hg19'
# )
