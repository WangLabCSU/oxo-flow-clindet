library(tidyverse)
library(GenomicRanges)

input_cnv <- snakemake@input[['cnv_rdata']]
input_cnv <- R.utils::getAbsolutePath(input_cnv)


threads = snakemake@threads

output_seg <- snakemake@output[['seg']]

out_absolute_dir <- snakemake@output[['absolute_dir']]

sample_name <- basename(input_cnv) %>% str_remove('_ASCAT.rdata')

load(input_cnv)

snp_df <- ascat.bc$SNPpos
snp_df$end <- snp_df$Position + 1

snp_gr <- makeGRangesFromDataFrame(
  snp_df,seqnames.field = 'Chromosome',start.field ='Position',
  end.field = 'end'
)

segment_gr <- makeGRangesFromDataFrame(
  ascat.output$segments,
  seqnames.field = 'chr',start.field = 'startpos',end.field = 'endpos'
)

seg <- ascat.output$segments
ov_reg <- findOverlaps(snp_gr,segment_gr)

ov_df <- as.data.frame(ov_reg) %>% distinct(queryHits,.keep_all=T)

snp_logr <- cbind(ov_df ,ascat.bc$Tumor_LogR)
#### as safe way maybe
# snp_logr <- cbind(ov_df ,ascat.bc$Tumor_LogR[ov_df$queryHits,])
snp_sum <- snp_logr %>% group_by(subjectHits) %>% summarise(
    n = n(),
    # seg.mean = log2(2*mean(2^Tumor)) - 1
    seg.mean = mean(Tumor,na.rm = T) # ref :https://www.biostars.org/p/306180/ https://www.biostars.org/p/315164/ and ASCAT raw data
)

# snp_sum <- as.data.frame(ov_reg) %>% count(subjectHits)
# snp_sum

# seg$Segment_Mean <- log2(seg$nMajor + seg$nMinor) -1

seg$Segment_Mean <- snp_sum$seg.mean
seg$Num_Probes <- 0
seg$Num_Probes[snp_sum$subjectHits] <- snp_sum$n
seg$Sample = sample_name

Seg <- seg %>% dplyr::select(Sample,chr,startpos,endpos,Num_Probes,Segment_Mean)
colnames(Seg) <- c("Sample","Chromosome","Start","End", "Num_Probes","Segment_Mean")
write_tsv(Seg,output_seg)


