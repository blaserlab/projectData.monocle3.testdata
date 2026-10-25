# Built-tarball installation. Slow, so only from a source checkout and never on CRAN.

test_that("the package builds and installs from a source tarball without the Monocle stack", {
  skip_on_cran()
  built <- built_fixture_lib()
  expect_true(file.exists(built$tarball))
  expect_lt(file.size(built$tarball), 1e6)
  root <- file.path(built$lib, pkg_name)
  expect_true(file.exists(file.path(root, "projectData", "manifest.dcf")))
  expect_setequal(list.files(file.path(root, "extdata")), cds_resources)
  expect_false(any(grepl("data-raw/git_commands|library_catalogs",
                         list.files(root, recursive = TRUE))))

  # Only base packages and testthat-free code are needed to use the data
  out <- rscript_lines(c(
    sprintf("loadedNamespaces_before <- loadedNamespaces()"),
    sprintf("utils::data('ordinary_data', package = '%s')", pkg_name),
    sprintf("utils::data('cds_multiassay', package = '%s')", pkg_name),
    "cat('ORDINARY', identical(ordinary_data, data.frame(letters = letters[1:10])), '\\n')",
    "cat('PLACEHOLDER', class(cds_multiassay), '\\n')",
    "cat('HEAVY', any(c('monocle3', 'BPCells', 'projectData') %in% loadedNamespaces()), '\\n')"
  ), built$lib)
  expect_true(any(out == "ORDINARY TRUE"), info = paste(out, collapse = "\n"))
  expect_true(any(out == "PLACEHOLDER projectData_monocle3_placeholder"))
  expect_true(any(out == "HEAVY FALSE"))
})
