# Generates the fixtures of projectData.monocle3.testdata:
#
#   cds_bpcells_a   100 genes x 40 cells, BPCells-backed counts
#   cds_bpcells_b    80 genes x 30 cells, independent BPCells-backed counts
#   cds_multiassay  100 genes x 40 cells, BPCells-backed counts plus
#                   `logcounts` and `spliced` assays and `ADT`/`CRISPR` altExps
#   ordinary_data   data.frame(letters = letters[1:10]), a conventional dataset
#
# All data are synthetic. This script is never run by installation, loading or
# R CMD check. Run it deliberately from the package root:
#
#   Sys.setenv(PDM3_TESTDATA_REGENERATE = "true")   # only to replace fixtures
#   source("inst/data-raw/generate_testdata.R")
#
# Without PDM3_TESTDATA_REGENERATE=true the script refuses to touch existing
# fixtures. Needs projectData, projectData.monocle3, monocle3 and BPCells.
# Counts are drawn with fixed seeds, so the matrices are reproducible; the
# BPCells/Monocle files on disk are written by those packages and are not
# promised to be byte-identical across versions (see inst/provenance).

library(Matrix)
library(SummarizedExperiment)
library(SingleCellExperiment)

stopifnot(
  "run from the package root" = identical(
    tryCatch(read.dcf("DESCRIPTION", "Package")[[1]], error = \(e) NA_character_),
    "projectData.monocle3.testdata"
  )
)
for (pkg in c("projectData", "projectData.monocle3", "monocle3", "BPCells")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("generating the fixtures requires the '", pkg, "' package", call. = FALSE)
  }
}

cds_names <- c("cds_bpcells_a", "cds_bpcells_b", "cds_multiassay")
# Fixtures of earlier package versions that the current contract no longer has
retired <- "cds_already_bpcells"

regenerate <- identical(Sys.getenv("PDM3_TESTDATA_REGENERATE"), "true")
existing <- c(
  file.path("inst", "extdata", c(cds_names, retired)),
  file.path("data", paste0(c(cds_names, retired, "ordinary_data"), ".rda"))
)
existing <- existing[file.exists(existing)]
if (length(existing) && !regenerate) {
  stop(
    "fixtures already exist:\n  ", paste(existing, collapse = "\n  "),
    "\nSet PDM3_TESTDATA_REGENERATE=true to replace them deliberately.",
    call. = FALSE
  )
}

seed_rng <- function(seed) {
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
}

as_csc <- function(m) methods::as(Matrix(m, sparse = TRUE), "CsparseMatrix")

# Sparse nonnegative integer counts
make_counts <- function(n_genes, n_cells, seed, gene_prefix, cell_prefix) {
  seed_rng(seed)
  as_csc(matrix(
    rpois(n_genes * n_cells, lambda = 2),
    nrow = n_genes,
    dimnames = list(
      sprintf("%s%03d", gene_prefix, seq_len(n_genes)),
      sprintf("%s%03d", cell_prefix, seq_len(n_cells))
    )
  ))
}

make_cds <- function(counts, seed) {
  seed_rng(seed)
  n_cells <- ncol(counts)
  cells <- data.frame(
    condition = rep(c("control", "treated"), length.out = n_cells),
    batch = rep(c("b1", "b2"), each = ceiling(n_cells / 2))[seq_len(n_cells)],
    row.names = colnames(counts)
  )
  genes <- data.frame(
    gene_short_name = paste0("sym_", rownames(counts)),
    row.names = rownames(counts)
  )
  cds <- monocle3::new_cell_data_set(counts, cell_metadata = cells, gene_metadata = genes)
  reducedDim(cds, "PCA") <- matrix(
    rnorm(n_cells * 5), nrow = n_cells,
    dimnames = list(colnames(counts), paste0("PC", 1:5))
  )
  cds
}

# Small in-memory experiment with its own identifiers and metadata
make_alt_exp <- function(cds, n_rows, prefix, seed) {
  seed_rng(seed)
  m <- matrix(
    rpois(n_rows * ncol(cds), lambda = 5), nrow = n_rows,
    dimnames = list(sprintf("%s_%02d", prefix, seq_len(n_rows)), colnames(cds))
  )
  se <- SummarizedExperiment(
    assays = list(counts = m),
    rowData = DataFrame(feature_type = rep(prefix, n_rows), row.names = rownames(m)),
    colData = DataFrame(row.names = colnames(cds))
  )
  metadata(se) <- list(origin = "synthetic", feature_set = prefix)
  se
}

counts_a <- make_counts(100, 40, seed = 4817, gene_prefix = "gA", cell_prefix = "cA")
counts_b <- make_counts(80, 30, seed = 7251, gene_prefix = "gB", cell_prefix = "cB")
counts_m <- make_counts(100, 40, seed = 2963, gene_prefix = "gM", cell_prefix = "cM")

cds_bpcells_a <- make_cds(counts_a, seed = 5120)
cds_bpcells_b <- make_cds(counts_b, seed = 6034)

cds_multiassay <- make_cds(counts_m, seed = 8396)
dense_m <- as.matrix(counts_m)
assay(cds_multiassay, "logcounts") <- as_csc(
  log1p(t(t(dense_m) / pmax(colSums(dense_m), 1)) * 1e4)
)
seed_rng(3571)
assay(cds_multiassay, "spliced") <- as_csc(matrix(
  rbinom(length(dense_m), size = as.vector(dense_m), prob = 0.7),
  nrow = nrow(dense_m), dimnames = dimnames(dense_m)
))
altExp(cds_multiassay, "ADT") <- make_alt_exp(cds_multiassay, 10, "ADT", seed = 1409)
altExp(cds_multiassay, "CRISPR") <- make_alt_exp(cds_multiassay, 6, "CRISPR", seed = 9182)

fixtures <- list(
  cds_bpcells_a = cds_bpcells_a,
  cds_bpcells_b = cds_bpcells_b,
  cds_multiassay = cds_multiassay
)

if (regenerate) {
  # Retired resources are removed here; the manifest is reconciled below
  unlink(file.path("inst", "extdata", retired), recursive = TRUE)
  unlink(file.path("data", paste0(retired, ".rda")))
}

# Counts are converted to BPCells by save_monocle_cds()
for (nm in names(fixtures)) {
  projectData.monocle3::save_monocle_cds(
    fixtures[[nm]], name = nm, package_dir = ".", overwrite = regenerate
  )
}

projectData.monocle3::write_monocle_manifest(package_dir = ".")

# Ordinary dataset: deliberately no manifest entry
ordinary_data <- data.frame(letters = letters[1:10])
dir.create("data", showWarnings = FALSE)
save(ordinary_data, file = file.path("data", "ordinary_data.rda"),
     version = 2, compress = "xz")

dir.create(file.path("inst", "provenance"), recursive = TRUE, showWarnings = FALSE)
versions <- vapply(
  c("projectData", "projectData.monocle3", "monocle3", "BPCells",
    "SingleCellExperiment", "SummarizedExperiment", "Matrix"),
  \(p) as.character(utils::packageVersion(p)), character(1)
)
write.dcf(
  data.frame(
    Script = "inst/data-raw/generate_testdata.R",
    Seeds = "counts 4817/7251/2963; cds 5120/6034/8396; spliced 3571; altExp 1409/9182",
    R = paste(R.version$major, R.version$minor, sep = "."),
    Platform = R.version$platform,
    t(versions),
    check.names = FALSE
  ),
  file = file.path("inst", "provenance", "fixtures.dcf")
)

message("Fixtures written. Run devtools::document() and the tests next.")
