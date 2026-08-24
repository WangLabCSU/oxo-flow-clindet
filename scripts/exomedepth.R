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
    my.counts <- getBamCounts(bed.frame = target.df,
            bam.files = bam.files,
            include.chr = FALSE,  # the target BED already carries chr-prefixed names
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
colnames(my.counts.df) <- c('chromosome','start','end','exon','GC','tumor','normal')

myTest <- somatic.CNV.call(normal  = my.counts.df$normal,
                            tumor = my.counts.df$tumor,
                            prop.tumor = 0.1,
                            chromosome = my.counts.df$chromosome,
                            start = my.counts.df$start,
                            end = my.counts.df$end,
                            names = my.counts.df$exon)

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


all.exons <- AnnotateExtra(x = myTest,
    reference.annotation = exons.hg19.GRanges,
    min.overlap = 0.001,
    column.name = 'exons.hg19'
)

saveRDS(all.exons,out_exom_rds)
readr::write_tsv(all.exons@CNV.calls,out_exom_depth)
# all.ex <- AnnotateExtra(x = myTest2,
#     reference.annotation = exons.hg19.GRanges,
#     min.overlap = 0.001,
#     column.name = 'exons.hg19'
# )
