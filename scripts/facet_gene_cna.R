library(tidyverse)
library(R.utils)
genome_version <- snakemake@wildcards[['genome_version']]

input_gene_level <- snakemake@input[['gene_level']]
input_gene_level <- getAbsolutePath(input_gene_level)

output_cna <- snakemake@output[['cna']]

gene_cna <- read_tsv(input_gene_level)

results <- gene_cna %>% filter(filter == 'PASS',cn_state %in% c('GAIN','HETLOSS'),cf.em >= 0.5,!spans_segs) %>%
  mutate(
    call = case_when(
      # 强缺失
      cn_state == "HOMDEL" ~ "Deletion",
      
      # LOH型缺失
      cn_state == "HETLOSS" & lcn.em == 0 ~ "Deletion",
      
      # 扩增（高拷贝优先）
      tcn.em >= 6 ~ "Amplification",
      
      # 中等扩增（可选）
      cn_state == "GAIN" & tcn.em >= 3 ~ "Amplification",
      
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(gene) %>%
  summarise(
    call = case_when(
      any(call == "Deletion") ~ "Deletion",
      any(call == "Amplification") ~ "Amplification",
      TRUE ~ "Wild type"
    ),
    .groups = "drop"
  )

write_tsv(results,output_cna)