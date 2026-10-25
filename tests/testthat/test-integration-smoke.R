# One end-to-end smoke test through projectData in a clean child process.
# Detailed adapter behaviour is covered by the projectData.monocle3 tests.

test_that("the installed fixtures open through projectData and the Monocle adapter", {
  skip_on_cran()
  for (p in c("projectData", "projectData.monocle3", "monocle3", "BPCells",
              "SummarizedExperiment", "SingleCellExperiment")) {
    skip_if_not_installed(p)
  }
  built <- built_fixture_lib()
  out <- rscript_lines(c(
    sprintf("projectData::data_attach('%s')", pkg_name),
    sprintf("for (nm in c(%s)) {", paste0("'", cds_resources, "'", collapse = ", ")),
    "  cds <- get(nm)",
    "  cat('CDS', nm, class(cds)[1], dim(cds), '\\n')",
    "}",
    "cds <- cds_multiassay",
    "cat('ASSAYS', sort(SummarizedExperiment::assayNames(cds)), '\\n')",
    "cat('ALTEXPS', sort(SingleCellExperiment::altExpNames(cds)), '\\n')",
    "cat('BACKED', methods::is(SummarizedExperiment::assay(cds, 'counts'), 'MatrixDir'), '\\n')"
  ), built$lib)
  expect_true(any(out == "CDS cds_bpcells_a cell_data_set 100 40"), info = paste(out, collapse = "\n"))
  expect_true(any(out == "CDS cds_bpcells_b cell_data_set 80 30"))
  expect_true(any(out == "CDS cds_multiassay cell_data_set 100 40"))
  expect_true(any(out == "ASSAYS counts counts_row_order logcounts spliced") ||
                any(out == "ASSAYS counts logcounts spliced"))
  expect_true(any(out == "ALTEXPS ADT CRISPR"))
  expect_true(any(out == "BACKED TRUE"))
})
