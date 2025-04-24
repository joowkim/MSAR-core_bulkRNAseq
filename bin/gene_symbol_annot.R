args <- commandArgs(TRUE)

# star count matrix
star_mat <- args[1]

# gene_id, symbol annot table
annot_mat <- args[2]

library(readr)
library(dplyr)

cnt_df <- read_tsv(star_mat)
annot_df <- read_csv(annot_mat)

merged_cnt_df <- left_join(cnt_df, annot_df, by = c('gene_id')) |> relocate(gene_name, .after=gene_id) |> rename(gene_symbol=gene_name)

write_tsv(merged_cnt_df, "gene_cnt_table.tsv")