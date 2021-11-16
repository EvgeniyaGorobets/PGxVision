# Test queryGene function
test_that("queryGene doesn't accept faulty user input", {
  badGene <- 100
  expect_error(queryGene(badGene, "GO:CC"))

  badSubtype <- "GO:BF"
  expect_error(queryGene("geneId", badSubtype))
})

test_that("queryGene returns corrects gene set IDs", {
  geneSetIds <- queryGene("ENSG00000012124", "GO:CC")
  expected <- c("M17635", "M25848", "M11182", "M9432", "M5343", "M2251",
                "M9232", "M17771", "M17514", "M17678")
  expect_equal(geneSetIds, expected)
})


# Test expandGeneSet function
test_that("expandGeneSets doesn't accept faulty user input", {
  badIds <- c(1, 2, 3)
  expect_error(expandGeneSets(badIds))

  badSubtype <- "GO:BF"
  expect_error(queryGene("geneId", badSubtype))
})

test_that("expandGeneSets correctly expands gene sets", {
  testGene <- "ENSG00000012124"
  geneSetIds <- queryGene(testGene, "GO:CC")
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
