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
projectData.monocle3::save_monocle_cds(
  cds = cds_bpcells_a,
  name = "cds_bpcells_a",
  package_dir = ".",
  overwrite = TRUE,
  convert_to_bpcells = TRUE
)


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
projectData.monocle3::save_monocle_cds(
  cds = cds_bpcells_b,
  name = "cds_bpcells_b",
  package_dir = ".",
  overwrite = TRUE,
  convert_to_bpcells = TRUE
)

ordinary_data <- data.frame(letters = letters[1:10])
save(ordinary_data, file = "data/ordinary_data.rda")
