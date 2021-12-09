# Test buildNetworkPlot function
test_that("buildNetworkPlot doesn't accept faulty user input", {
  df <- data.frame("gs1" = c("a", "b", "c"), "gs2" = c("b", "c", "a"),
                   "similarity" = c(0, 0.5, 0.8))

  notADf <- list(df)
  expect_error(buildNetworkPlot(notADf))

  badDf <- df
  names(badDf) <- c("gs1", "gs2", "Similarity")
  expect_error(buildNetworkPlot(badDf))

  badCutoff <- "1"
  expect_error(buildNetworkPlot(gsSimilarityDf, similarityCutoff=badCutoff))

  badCutoff <- 3
  expect_error(buildNetworkPlot(gsSimilarityDf, similarityCutoff=badCutoff))

  badTitle <- 123
  expect_error(buildNetworkPlot(gsSimilarityDf, title=badTitle))
})


test_that("buildNetworkPlot sets edge weights to gene set similarity", {
  similarityScores <- c(0.1, 0.2, 0.3)
  df <- data.frame("gs1" = c("a", "b", "c"), "gs2" = c("b", "c", "a"),
                   "similarity" = similarityScores)
  graph <- buildNetworkPlot(df, similarityCutoff=0)
  expect_equal(graph$x$edges$value, similarityScores)
})


# Test geneSetAnalysis function
# Since this is a wrapper around other functions which are already thoroughly
# tested, I only test that errors are handled correctly
test_that("geneSetAnalysis aborts if getGeneSets returns 0 gene sets", {
  # getGeneSets("ENSG00000043143", "GO:MF") returns empty data.table
  expect_error(geneSetAnalysis("ENSG00000043143", "GO:MF"))
})
