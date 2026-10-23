library(tidyverse)
library(monocle3)
library(blaseRtools)
library(SummarizedExperiment)
library(Matrix)
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


set.seed(42)

# Dimensions
n_genes <- 100
n_cells <- 40

# Generate count matrix
counts <- matrix(
  rpois(n_genes * n_cells, lambda = 3),
  nrow = n_genes,
  dimnames = list(
    paste0("gene_", seq_len(n_genes)),
    paste0("cell_", seq_len(n_cells))
  )
)

# Cell metadata
cell_metadata <- data.frame(
  condition = rep(c("control", "treated"), each = 20),
  row.names = colnames(counts)
)

# Gene metadata
gene_metadata <- data.frame(
  gene_short_name = rownames(counts),
  row.names = rownames(counts)
)

# Create Monocle CDS
cds_multiassay <- monocle3::new_cell_data_set(
  expression_data = Matrix(counts, sparse = TRUE),
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

# Additional assay 1: log-normalized counts
normalized <- log1p(
  sweep(counts, 2, pmax(colSums(counts), 1), "/") * 10000
)

# Additional assay 2: simulated spliced counts
spliced <- matrix(
  rbinom(
    length(counts),
    size = as.vector(counts),
    prob = 0.7
  ),
  nrow = n_genes,
  dimnames = dimnames(counts)
)

# Add assays
SummarizedExperiment::assay(
  cds_multiassay, "logcounts"
) <- normalized

SummarizedExperiment::assay(
  cds_multiassay, "spliced"
) <- Matrix(spliced, sparse = TRUE)

# Add a dimensionality reduction for additional testing
set.seed(123)

SingleCellExperiment::reducedDim(
  cds_multiassay, "PCA"
) <- matrix(
  rnorm(n_cells * 5),
  nrow = n_cells,
  dimnames = list(colnames(counts), paste0("PC", 1:5))
)

# Inspect
cds_multiassay
SummarizedExperiment::assayNames(cds_multiassay)

library(SingleCellExperiment)
library(SummarizedExperiment)
library(Matrix)

set.seed(123)

n_cells <- ncol(cds_multiassay)
cell_names <- colnames(cds_multiassay)

# ---- Alternative experiment 1: ADT ----

adt_counts <- matrix(
  rpois(12 * n_cells, lambda = 5),
  nrow = 12,
  dimnames = list(
    paste0("ADT_", seq_len(12)),
    cell_names
  )
)

adt <- SingleCellExperiment(
  assays = list(
    counts = Matrix(adt_counts, sparse = TRUE),
    logcounts = log1p(adt_counts)
  ),
  rowData = DataFrame(
    protein_name = paste0("CD", seq_len(12))
  )
)

metadata(adt)$assay_type <- "Antibody-derived tags"

altExp(cds_multiassay, "ADT") <- adt


# ---- Alternative experiment 2: CRISPR ----

crispr_counts <- matrix(
  rpois(8 * n_cells, lambda = 2),
  nrow = 8,
  dimnames = list(
    paste0("guide_", seq_len(8)),
    cell_names
  )
)

crispr <- SingleCellExperiment(
  assays = list(
    counts = Matrix(crispr_counts, sparse = TRUE)
  ),
  rowData = DataFrame(
    target_gene = paste0("GeneTarget_", seq_len(8))
  )
)

metadata(crispr)$assay_type <- "CRISPR guide capture"

altExp(cds_multiassay, "CRISPR") <- crispr

cds_multiassay

altExpNames(cds_multiassay)
assayNames(altExp(cds_multiassay, "ADT"))
assayNames(altExp(cds_multiassay, "CRISPR"))

projectData.monocle3::save_monocle_cds(
  cds = cds_multiassay,
  name = "cds_multiassay",
  package_dir = ".",
  convert_to_bpcells = TRUE,
  overwrite = TRUE
)

ordinary_data <- data.frame(letters = letters[1:10])
save(ordinary_data, file = "data/ordinary_data.rda")
