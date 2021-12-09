# Test that buildVolcanoPlot returns an error if the user input is wrong
test_that("buildVolcanoPlot doesn't accept faulty user input", {
  notDf <- list(3, "c", data.frame(3, 4, 5))
  expect_error(buildVolcanoPlot(biomarkerDf=notDf))

  df <- Biomarkers
  badDf <- data.frame(Biomarkers)
  names(badDf)[2] <- "drug"
  expect_error(buildVolcanoPlot(biomarkerDF=badDf))

  experiment <- setNames(c("Lung", "Panobinostat", "rna"),
                         c("tissue", "compound", "mDataType"))
  bigPval <- 2
  negPval <- -1
  expect_error(buildVolcanoPlot(df, experiment, bigPval))
  expect_error(buildVolcanoPlot(df, experiment, negPval))

})

