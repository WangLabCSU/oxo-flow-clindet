library(VariantAnnotation)
library(tidyverse)
library(stringi)

input_vcf <- snakemake@input[['vcf']]
output_vcf <- snakemake@output[['vcf']]
# output_vcf <- snakemake@params[['vcf']]
sample_name <- snakemake@wildcards[['sample']]
x = readVcf(input_vcf)
colnames(x) <- str_c(sample_name,c('_NC',''))
writeVcf(x, output_vcf,index = F)