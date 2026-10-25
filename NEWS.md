# projectData.monocle3.testdata 0.0.0.9004

* Fixtures now match the contract of the `projectData.monocle3` 0.0.0.9007 integration tests:
  `cds_bpcells_a` (100 x 40), `cds_bpcells_b` (80 x 30), `cds_multiassay` (100 x 40,
  with `logcounts`, `spliced`, and `ADT`/`CRISPR` alternative experiments) and
  `ordinary_data`. The 0.0.0.9003 resource `cds_already_bpcells` is retired.
* `inst/data-raw/generate_testdata.R` generates everything with fixed seeds through
  `projectData.monocle3::save_monocle_cds()` and refuses to overwrite fixtures unless
  `PDM3_TESTDATA_REGENERATE=true`. Generation versions are recorded in
  `inst/provenance/fixtures.dcf`.
* Fixtures (`data/*.rda`, `inst/extdata/`) are now committed, so the package installs
  from GitHub or a tarball without the Monocle stack.
* Added a testthat suite; structural checks need only testthat.
* Package metadata and `R/data.R` documentation completed; `License` is `MIT + file LICENSE`.

# projectData.monocle3.testdata 0.0.0.9003

* Complex CDS objects (superseded).
