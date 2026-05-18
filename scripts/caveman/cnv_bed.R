library(tidyverse)

cnv_rdata <- snakemake@input[['rdata']]
tcnv_bed <- snakemake@output[['Tcnv']]
nccnv_bed <- snakemake@output[['NCcnv']]
fai_file <- snakemake@params[['fasta_fai']]

load(cnv_rdata)
fai_contigs <- read_tsv(fai_file,col_names = c('chr','endpos','d1','d2','d3'))
chr_index <- c(1:22,'X','Y')
fai_contigs <- fai_contigs %>% 
    mutate(startpos = 1) %>% 
    select(chr,startpos,endpos)
if (is.null(ascat.output$segments)){
    tumor_seg <- fai_contigs
    tumor_seg$cn <- 2
    tumor_seg <- tumor_seg %>% dplyr::select(chr,startpos,endpos,cn)
    write_tsv(tumor_seg,tcnv_bed,col_names = F)
    write_tsv(tumor_seg,nccnv_bed,col_names = F)
} else {
    tumor_seg <- ascat.output$segments
    tumor_seg$cn <- tumor_seg$nMajor + tumor_seg$nMinor
    tumor_seg$chr <- paste0('chr',tumor_seg$chr)
    tumor_seg <- tumor_seg %>% dplyr::select(chr,startpos,endpos,cn)
    write_tsv(tumor_seg,tcnv_bed,col_names = F)

    normal_seg <- ascat.output$segments
    normal_seg$cn <- 2
    normal_seg$chr <- paste0('chr',normal_seg$chr)
    normal_seg <- normal_seg %>% dplyr::select(chr,startpos,endpos,cn)
    write_tsv(normal_seg,nccnv_bed,col_names = F)

}

