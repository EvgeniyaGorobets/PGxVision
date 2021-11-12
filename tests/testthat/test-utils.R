# Test selectExperiment function
test_that("selectExperiment doesn't accept faulty experiment vector", {
  dt <- setDT(Biomarkers)
  notVector <- list("Lung", "Panobinostat", "rna")
  expect_error(selectExperiment(dt, experiment=notVector))

  badVector <- setNames(c("Lung", "Panobinostat", "rna"),
                        c("tissue", "drug", "mDataType"))
  expect_error(selectExperiment(dt, experiment=badVector))
})

test_that("selectExperiment warns when it returns empty data.table", {
  dt <- setDT(Biomarkers)
  experiment <- setNames(c("Lung", "unknownDrug", "rna"),
                         c("tissue", "compound", "mDataType"))
  expect_warning(selectExperiment(dt, experiment))
})

# TODO: write a test that checks proper subsetting occurs
