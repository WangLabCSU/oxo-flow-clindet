library(R.utils)
library(glue)

input_rdata <- snakemake@input[['rdata']]
input_rdata <- getAbsolutePath(input_rdata)


# output a tsv which column names as same as purple output, so user can choose which result to use
output_tsv <- snakemake@output[['tsv']]
output_tsv <- getAbsolutePath(output_tsv)

load(input_rdata)

pp_df <- data.frame(
    ploidy = ascat.output[['ploidy']],
    purity = ascat.output[['purity']],
    goodnessOfFit = ascat.output[['goodnessOfFit']]
)

write.table(pp_df, file = output_tsv, row.names=FALSE, sep="\t",quote=FALSE)