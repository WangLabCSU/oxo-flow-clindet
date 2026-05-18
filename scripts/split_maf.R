library(readr)
library(dplyr)

input_maf <- snakemake@input[['maf']]
output_snp <- snakemake@output[['snp']]
output_indel <- snakemake@output[['indel']]

lines <- readLines(if (grepl("\\.gz$", input_maf)) gzfile(input_maf, "rt") else input_maf)

comment_idx <- grepl("^#", lines)
first_data_idx <- which(!comment_idx)[1]

comment_lines <- if (first_data_idx > 1) lines[1:(first_data_idx - 1)] else character()
data_text <- paste(lines[first_data_idx:length(lines)], collapse = "\n")

maf <- read_tsv(I(data_text), show_col_types = FALSE)

if (!"Variant_Type" %in% colnames(maf)) {
    stop("Variant_Type column not found in MAF header.")
}

maf_snp <- maf %>% filter(Variant_Type == "SNP")
maf_indel <- maf %>% filter(Variant_Type %in% c("INS", "DEL"))

write_with_comments <- function(df, outfile, comments = character()) {
  if (length(comments) > 0) {
    writeLines(comments, outfile, sep = "\n")
    write("\n", file = outfile, append = TRUE)
  } else {
    file.create(outfile)
  }

  readr::write_tsv(df, outfile, append = TRUE, col_names = TRUE)
}

write_with_comments(maf_snp, output_snp, comment_lines)
write_with_comments(maf_indel, output_indel, comment_lines)