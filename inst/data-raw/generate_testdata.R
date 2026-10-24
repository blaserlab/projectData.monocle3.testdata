library(monocle3)
library(SummarizedExperiment)
library(Matrix)

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


# Start with the original in-memory CDS
cds_bpcells <- monocle3::convert_counts_matrix(
  cds_multiassay,
  matrix_control = list(
    matrix_class = "BPCells",
    matrix_path = file.path(
      tempdir(),
      "test_bpcells_conversion"
    )
  )
)

# Confirm that counts are now BPCells-backed
stopifnot(
  methods::is(
    SummarizedExperiment::assay(cds_bpcells, "counts"),
    "IterableMatrix"
  )
)

# Inspect the assays and alternative experiments
assayNames(cds_bpcells)
altExpNames(cds_bpcells)

# Save without converting again
projectData.monocle3::save_monocle_cds(
  cds = cds_bpcells,
  name = "cds_already_bpcells",
  package_dir = ".",
  convert_to_bpcells = FALSE
)

projectData.monocle3::save_monocle_cds(
  cds = cds_multiassay,
  name = "cds_multiassay",
  package_dir = ".",
  convert_to_bpcells = TRUE
)
