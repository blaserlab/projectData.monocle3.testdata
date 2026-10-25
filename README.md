# projectData.monocle3.testdata

Small synthetic fixtures for testing the integration of
[projectData](https://github.com/blaserlab/projectData), its Monocle 3 adapter
[projectData.monocle3](https://github.com/blaserlab/projectData.monocle3) and
BPCells-backed `cell_data_set` objects. No patient or proprietary data.

| Resource | Contents |
|---|---|
| `cds_bpcells_a` | 100 genes x 40 cells, BPCells counts |
| `cds_bpcells_b` | 80 genes x 30 cells, independent BPCells counts |
| `cds_multiassay` | 100 x 40, BPCells counts, `logcounts`, `spliced`, altExps `ADT` and `CRISPR` |
| `ordinary_data` | `data.frame(letters = letters[1:10])`, no manifest entry |

The CDS resources are saved Monocle directories in `inst/extdata/<name>/`, registered
in `inst/projectData/manifest.dcf` (loader `projectData.monocle3::cds`), with
placeholder objects in `data/`. `utils::data()` loads only a placeholder; use
`projectData::data_attach("projectData.monocle3.testdata")` to get the real CDS.

## Install

Installing needs none of projectData, monocle3 or BPCells.

``` r
remotes::install_github("blaserlab/projectData.monocle3.testdata")
# or from a local checkout / built tarball
install.packages("projectData.monocle3.testdata_0.0.0.9004.tar.gz", repos = NULL)
```

Using the CDS resources needs `projectData`, `projectData.monocle3`, `monocle3` and `BPCells`
(all in `Suggests`, with `Remotes`).

## Validate

``` r
devtools::test()   # structure checks always; Monocle checks skip if the stack is missing
```

Structure checks need only testthat. A test that builds a tarball and installs it into a
temporary library runs from a source checkout when `NOT_CRAN=true`. The CDS behaviour itself is
covered by the integration tests in `projectData.monocle3`.

## Regenerate

Fixtures are committed and the generator refuses to overwrite them. From the package root,
with the full stack installed:

``` r
Sys.setenv(PDM3_TESTDATA_REGENERATE = "true")
source("inst/data-raw/generate_testdata.R")
devtools::document()
```

Run it in a fresh R session (`Rscript --vanilla`) so that no other copy of the adapter is
picked up. Seeds are fixed; the versions used are recorded in `inst/provenance/fixtures.dcf`.
The generator is never run by installation, loading or `R CMD check`.

Monocle's native files record some absolute paths from the generating machine (a staging
directory in `file_index.rds` and the temporary BPCells directory in `cds_object.rds`). They
are informational only: the loader copies the matrix from the installed files.

## Version coupling

The `projectData.monocle3` integration tests (`helper-integration.R`) mention fixture
version 0.0.0.9002 in comments and in a skip message; they check resources, not the version.
