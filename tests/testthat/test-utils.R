# Test selectExperiment function
test_that("selectExperiment doesn't accept faulty user input", {
  # Accepts data.table only, not data.frame
  expect_error(selectExperiment(data.frame(a=c(1, 2), b=c(3, 4))))

  df <- copy(Biomarkers)
  dt <- setDT(df)

  # only accepts strings for other params
  expect_error(selectExperiment(dt, tissue = 3))
  expect_error(selectExperiment(dt, compound = list("notAString")))
  expect_error(selectExperiment(dt, mDataType = F))
})

test_that("selectExperiment warns when it returns empty data.table", {
  df <- copy(Biomarkers)
  dt <- setDT(df)
  expect_warning(selectExperiment(dt, "Lung", "unknownDrug", "rna"))
})

test_that("selectExperiment returns correct rows", {
  df <- copy(Biomarkers)
  dt <- setDT(df)
  result <- selectExperiment(dt, "Lung", "Trametinib", "rna")

  expect_true(all(result$tissue == "Lung"))
  expect_true(all(result$compound == "Trametinib"))
  expect_true(all(result$mDataType == "rna"))
  expect_equal(nrow(result), 1106)
})
