library(tidyverse)
library(glue)

input_cnv <- snakemake@input[['cnv']]
output_ascat <- snakemake@output[['ascat']]
brass_cnv <- snakemake@params[['brass_cnv']]
if (brass_cnv == 'ascat'){
    load(input_cnv)
    sink(output_ascat,append=FALSE,split = FALSE)
    print(glue('rho ',unname(ascat.output$purity)))
    print(glue('Ploidy ',unname(ascat.output$ploidy)))
    print(glue('GenderChr X'))
    print(glue('GenderChrFound N'))
    sink()
} else if (brass_cnv == 'purple') {
    df <- read_tsv(input_cnv)
    gender <- if_else(df$gender == 'MALE','Y','X')
    sink(output_ascat,append=FALSE,split = FALSE)
    print(glue('rho ',df$purity[1]))
    print(glue('Ploidy ',df$ploidy[1]))
    print(glue('GenderChr {gender}'))
    print(glue('GenderChrFound N'))
    sink()
}
