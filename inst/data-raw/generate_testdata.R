library(tidyverse)
library(monocle3)
library(blaseRtools)
set.seed(42)

counts_a <- matrix(
  rpois(100 * 40, lambda = 2),
  nrow = 100,
  dimnames = list(
    paste0("gene_", seq_len(100)),
    paste0("cell_", seq_len(40))
  )
)

cell_data_a <- DataFrame(row.names = paste0("cell_", seq_len(40)))
row_data_a <- DataFrame(row.names = paste0("gene_", seq_len(100)))
cds_bpcells_a <- monocle3::new_cell_data_set(expression_data = counts_a, cell_metadata = cell_data_a, gene_metadata = row_data_a)
save_monocle_disk(cds_bpcells_a, data_directory = "data", extdata_directory = "inst/extdata")

set.seed(24)

counts_b <- matrix(
  rpois(80 * 30, lambda = 2),
  nrow = 80,
  dimnames = list(
    paste0("gene_", seq_len(80)),
    paste0("cell_", seq_len(30))
  )
)

cell_data_b <- DataFrame(row.names = paste0("cell_", seq_len(30)))
row_data_b <- DataFrame(row.names = paste0("gene_", seq_len(80)))
cds_bpcells_b <- monocle3::new_cell_data_set(expression_data = counts_b, cell_metadata = cell_data_b, gene_metadata = row_data_b)
save_monocle_disk(cds_bpcells_b, data_directory = "data", extdata_directory = "inst/extdata")

ordinary_data <- data.frame(letters = letters[1:10])
save(ordinary_data, file = "data/ordinary_data.rda")
