# Helpers for the fixture-integrity tests. Everything here is base R (plus
# testthat/withr) so that structure checks run without the Monocle stack.

pkg_name <- "projectData.monocle3.testdata"
cds_resources <- c("cds_bpcells_a", "cds_bpcells_b", "cds_multiassay")
cds_dims <- list(cds_bpcells_a = c(100L, 40L), cds_bpcells_b = c(80L, 30L),
                 cds_multiassay = c(100L, 40L))
bpcells_files <- c("col_names", "idxptr", "index", "row_names", "shape",
                   "storage_order", "val", "version")

# Root of the package being tested: the installed package, or `inst/` of the
# source tree under devtools::test().
fixture_root <- function() {
  root <- system.file(package = pkg_name)
  if (!nzchar(root)) stop("package not found: ", pkg_name)
  root
}

# Source checkout root, or NULL when testing an installed package.
source_root <- function() {
  candidate <- normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = FALSE)
  desc <- file.path(candidate, "DESCRIPTION")
  if (file.exists(desc) && dir.exists(file.path(candidate, "inst")) &&
      identical(unname(read.dcf(desc, "Package")[1, 1]), pkg_name)) {
    candidate
  }
}

# Dimensions of a saved BPCells matrix directory, read from its own files.
bpcells_dims <- function(dir) {
  con <- file(file.path(dir, "shape"), "rb")
  on.exit(close(con))
  readBin(con, "raw", 8L)
  readBin(con, "integer", 2L, size = 4L)
}

# Fresh library holding the package built from the source tarball. Built once
# per test run.
.build_cache <- new.env()
built_fixture_lib <- function() {
  if (!is.null(.build_cache$lib)) return(.build_cache)
  src <- source_root()
  if (is.null(src)) testthat::skip("not running from a source checkout")
  work <- tempfile("fixture-build-")
  dir.create(work)
  lib <- file.path(work, "lib")
  dir.create(lib)
  old <- setwd(work)
  on.exit(setwd(old), add = TRUE)
  r <- file.path(R.home("bin"), "R")
  build <- suppressWarnings(system2(r, c("CMD", "build", shQuote(src)),
                                    stdout = TRUE, stderr = TRUE))
  tarball <- list.files(work, pattern = "\\.tar\\.gz$", full.names = TRUE)
  if (length(tarball) != 1L) stop("R CMD build failed:\n", paste(build, collapse = "\n"))
  install <- suppressWarnings(system2(
    r, c("CMD", "INSTALL", "-l", shQuote(lib), shQuote(tarball)),
    stdout = TRUE, stderr = TRUE
  ))
  if (!file.exists(file.path(lib, pkg_name, "DESCRIPTION"))) {
    stop("R CMD INSTALL failed:\n", paste(install, collapse = "\n"))
  }
  .build_cache$lib <- lib
  .build_cache$tarball <- tarball
  .build_cache$install_log <- install
  withr::defer(unlink(work, recursive = TRUE), testthat::teardown_env())
  .build_cache
}

# Runs R code in a clean child process that sees `lib` first.
rscript_lines <- function(lines, lib) {
  script <- tempfile(fileext = ".R")
  on.exit(unlink(script))
  writeLines(lines, script)
  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)),
    stdout = TRUE, stderr = TRUE,
    env = paste0("R_LIBS=", paste(c(lib, .libPaths()), collapse = .Platform$path.sep))
  ))
  trimws(out)
}

# DESCRIPTION of the package under test (source tree, or installed copy).
description_file <- function() {
  src <- source_root()
  if (is.null(src)) file.path(fixture_root(), "DESCRIPTION") else file.path(src, "DESCRIPTION")
}
