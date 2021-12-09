# Test getGeneSets function
test_that("getGeneSets doesn't accept faulty user input", {
  badGene <- 100
  expect_error(getGeneSets(badGene, "GO:CC"))

  badSubtype <- "GO:BF"
  expect_error(getGeneSets("geneId", badSubtype))
})

test_that("getGeneSets returns correct gene set IDs and columns", {
  geneSets <- getGeneSets("ENSG00000012124", "GO:CC")

  expectedCols <- c("gs_id", "gs_name", "gs_exact_source", "gs_url",
                    "gs_description")
  expect_equal(colnames(geneSets), expectedCols)

  expectedGeneSets <- c("M17635", "M25848", "M11182", "M9432", "M5343", "M2251",
                        "M9232", "M17771", "M17514", "M17678")
  expect_equal(geneSets$gs_id, expectedGeneSets)
})


# Test expandGeneSet function
test_that("expandGeneSets doesn't accept faulty user input", {
  badIds <- c(1, 2, 3)
  expect_error(expandGeneSets(badIds))

  badSubtype <- "GO:BF"
  expect_error(getGeneSets("geneId", badSubtype))
})

test_that("expandGeneSets correctly expands gene sets", {
  testGene <- "ENSG00000012124"
  geneSetIds <- getGeneSets(testGene, "GO:CC")$gs_id
  geneSets <- expandGeneSets(geneSetIds)
  # This test may need to be updated if the MSigDb changes
  expect_equal(geneSets, TestGeneSets)

  # Test that the test gene shows up in each gene set!
  expect_equal(geneSets[geneSets$ensembl_gene == testGene, gs_id], geneSetIds)
})


# Test overlapDistance helper function
test_that("overlapDistance computes overlap correctly", {
  expect_equal(overlapDistance(TestGeneSets, "M17635", "M25848"), 0.055267703)
  expect_equal(overlapDistance(TestGeneSets, "fakeGS1", "M25848"), 0)
  expect_equal(overlapDistance(TestGeneSets, "M17635", "fakeGS2"), 0)
})


# Test computeGeneSetSimilarity function
test_that("computeGeneSetSimilarity doesn't accept faulty user input", {
  notADf <- list(TestGeneSets)
  expect_error(computeGeneSetSimilarity(notADf))

  badDf <- TestGeneSets
  names(badDf) <- c("gene_set_id", "gene_id")
  expect_error(computeGeneSetSimilarity(badDf))
})

test_that("computeGeneSetSimilarity returns data.frame with correct columns", {
  df <- computeGeneSetSimilarity(TestGeneSets)
  expect_equal(names(df), c("gs1", "gs2", "similarity"))
  expect_true(is.data.frame(df))
})
