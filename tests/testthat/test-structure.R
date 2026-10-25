# Structure checks on the package as installed (or its source `inst/`).
# These need only testthat: no projectData, monocle3 or BPCells.

test_that("the saved CDS directories are complete", {
  extdata <- file.path(fixture_root(), "extdata")
  expect_setequal(list.files(extdata), cds_resources)
  for (nm in cds_resources) {
    dir <- file.path(extdata, nm)
    expect_true(file.exists(file.path(dir, "file_index.rds")), info = nm)
    expect_true(file.exists(file.path(dir, "cds_object.rds")), info = nm)
    expect_setequal(list.files(dir), c("bpcells_matrix_dir", "cds_object.rds", "file_index.rds"))
    expect_setequal(list.files(file.path(dir, "bpcells_matrix_dir")), bpcells_files)
  }
})

test_that("saved BPCells matrices have the expected shape and valid identifiers", {
  for (nm in cds_resources) {
    mdir <- file.path(fixture_root(), "extdata", nm, "bpcells_matrix_dir")
    expect_identical(bpcells_dims(mdir), cds_dims[[nm]], info = nm)
    for (f in c("row_names", "col_names")) {
      ids <- readLines(file.path(mdir, f))
      expect_length(ids, if (f == "row_names") cds_dims[[nm]][1] else cds_dims[[nm]][2])
      expect_false(anyNA(ids))
      expect_true(all(nzchar(ids)))
      expect_false(anyDuplicated(ids) > 0L, info = paste(nm, f))
    }
  }
})

test_that("stored counts are finite, nonnegative and whole numbers", {
  for (nm in cds_resources) {
    val <- file.path(fixture_root(), "extdata", nm, "bpcells_matrix_dir", "val")
    con <- file(val, "rb")
    readBin(con, "raw", 8L)
    x <- readBin(con, "double", file.size(val) / 8)
    close(con)
    expect_gt(length(x), 0L)
    expect_true(all(is.finite(x)) && all(x >= 0) && all(x == round(x)), info = nm)
  }
})

test_that("the manifest is valid and declares the adapter loader", {
  path <- file.path(fixture_root(), "projectData", "manifest.dcf")
  expect_true(file.exists(path))
  records <- read.dcf(path, all = TRUE)
  expect_identical(unname(as.character(records$ManifestVersion[[1]])), "1")
  # The version record is its own paragraph, followed by one per resource
  res <- records[!is.na(records$Resource), , drop = FALSE]
  expect_setequal(res$Resource, cds_resources)
  expect_false(anyDuplicated(res$Resource) > 0L)
  expect_true(all(res$Loader == "projectData.monocle3::cds"))
  expect_identical(res$Source, paste0("extdata/", res$Resource))
  expect_identical(res$Objects, res$Resource)
  expect_true(all(res$RecordedClass == "cell_data_set"))
  expect_true(all(res$AdapterPackage == "projectData.monocle3"))
  expect_true(all(dir.exists(file.path(fixture_root(), res$Source))))
  expect_false("ordinary_data" %in% res$Resource)
})

test_that("placeholders exist, are consistent, and do not restore the CDS", {
  index <- utils::data(package = pkg_name)$results[, "Item"]
  expect_setequal(index, c(cds_resources, "ordinary_data"))
  for (nm in cds_resources) {
    env <- new.env()
    utils::data(list = nm, package = pkg_name, envir = env)
    x <- get(nm, envir = env)
    expect_type(x, "character")
    expect_identical(class(x), "projectData_monocle3_placeholder")
    expect_identical(attr(x, "object_name"), nm)
    expect_identical(attr(x, "package"), pkg_name)
    expect_true(any(grepl("placeholder", x, fixed = TRUE)))
    expect_true(any(grepl(paste0("rm(", nm, ")"), x, fixed = TRUE)))
    expect_true(any(grepl("data_attach", x, fixed = TRUE)))
    expect_false(methods::is(x, "SummarizedExperiment"))
  }
})

test_that("ordinary_data is an ordinary data frame", {
  env <- new.env()
  utils::data("ordinary_data", package = pkg_name, envir = env)
  expect_identical(env$ordinary_data, data.frame(letters = letters[1:10]))
})

test_that("DESCRIPTION keeps installation free of the Monocle stack", {
  desc <- read.dcf(description_file())
  hard <- unlist(desc[1, intersect(c("Depends", "Imports", "LinkingTo"), colnames(desc))])
  hard <- trimws(gsub("\\(.*?\\)", "", unlist(strsplit(hard, ","))))
  expect_identical(setdiff(hard, "R"), character())
  expect_identical(unname(desc[1, "LazyData"]), "true")
})

test_that("generation provenance is recorded", {
  prov <- read.dcf(file.path(fixture_root(), "provenance", "fixtures.dcf"))
  expect_true(all(c("R", "monocle3", "BPCells", "projectData", "projectData.monocle3",
                    "Seeds") %in% colnames(prov)))
})

test_that("no development-machine paths are recorded in readable files", {
  files <- list.files(fixture_root(), recursive = TRUE, full.names = TRUE)
  abs_path <- "(^|[\"' =(])(/home/|/tmp/|/workspace/|/Users/|/var/folders/)|(^|[\"' (=])[A-Za-z]:[/\\\\]"
  # Personal notes in source checkouts; excluded from the built package
  files <- files[basename(files) != "git_commands.R"]
  opaque <- grepl("\\.rds$|\\.rda$|/(val|index|idxptr)$", files)
  for (f in files[!opaque]) {
    lines <- readLines(f, warn = FALSE)
    expect_false(any(grepl(abs_path, lines, perl = TRUE, useBytes = TRUE)),
                 info = paste("absolute path in", sub(fixture_root(), "", f, fixed = TRUE)))
  }
})

test_that("absolute paths inside Monocle's native files are informational only", {
  # monocle3 records the staging and temporary BPCells directories of the
  # machine that saved the CDS. They cannot be edited out of the native files
  # and the loader never reads them; the fixture must not depend on them, so
  # none of them may point at anything that still exists.
  #
  # cds_object.rds holds an S4 cell_data_set, so it is scanned as a decompressed
  # byte stream: dput()/str() would need the Bioconductor class definitions
  # (e.g. GenomeInfoDb), which are not available under R CMD check.
  serialized_paths <- function(f) {
    con <- gzfile(f, "rb")
    on.exit(close(con))
    bytes <- readBin(con, "raw", file.size(f) * 100)
    printable <- bytes >= as.raw(32) & bytes < as.raw(127)
    txt <- rawToChar(replace(bytes, !printable, as.raw(10)))
    m <- gregexpr("/[A-Za-z0-9_.@%+-]+(/[A-Za-z0-9_.@%+-]+)+", txt)
    unique(regmatches(txt, m)[[1]])
  }
  for (nm in cds_resources) {
    dir <- file.path(fixture_root(), "extdata", nm)

    # file_index.rds is a plain list: inspect it directly
    index <- readRDS(file.path(dir, "file_index.rds"))
    expect_false(file.exists(index$directory), info = nm)
    expect_false(any(grepl("^(/|[A-Za-z]:)", index$files$file_path)), info = nm)
    expect_true(all(file.exists(file.path(dir, index$files$file_path))), info = nm)

    recorded <- serialized_paths(file.path(dir, "cds_object.rds"))
    expect_false(any(file.exists(recorded)), info = nm)
    # Anything recorded must be a temporary location, never a project path
    expect_false(any(grepl(fixture_root(), recorded, fixed = TRUE)), info = nm)
    expect_false(any(grepl("/inst/extdata/", recorded, fixed = TRUE)), info = nm)
  }
})
