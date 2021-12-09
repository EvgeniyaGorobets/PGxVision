# Test that buildManhattanPlot returns an error if the user input is wrong
test_that("buildManhattanPlot doesn't accept faulty user input", {
  notDf <- list(3, "c", data.frame(3, 4, 5))
  expect_error(buildManhattanPlot(biomarkerDf=notDf))

  df <- Biomarkers
  badDf <- data.frame(Biomarkers)
  names(badDf)[2] <- "drug"
  expect_error(buildManhattanPlot(biomarkerDF=badDf))

  notDf <- setNames(c("a", 10), c("chrName", "chrLength"))
  expect_error(buildManhattanPlot(biomarkerDf=df, chromosomeDf=notDf))

  chrDf <- GRCh38.p13.Assembly
  badDf <- data.frame(GRCh38.p13.Assembly)
  names(badDf)[1] <- "chromosomeName"
  expect_error(buildManhattanPlot(biomarkerDf=df, chromosomeDf=badDf))

  experiment <- setNames(c("Lung", "Panobinostat", "rna"),
                         c("tissue", "compound", "mDataType"))
  bigPval <- 2
  negPval <- -1
  expect_error(buildManhattanPlot(biomarkerDf=df, chromosomeDf=chrDf,
                                  experiment=experiment, pValCutoff=bigPval))
  expect_error(buildManhattanPlot(biomarkerDf=df, chromosomeDf=chrDf,
                                  experiment=experiment, pValCutoff=negPval))

})

# Test absolutizeGenomicCoords (helper of buildManhattanPlot)
test_that("absolutizeGenomicCoords computes genomic coords correctly", {
  df <- Biomarkers
  chrDf <- GRCh38.p13.Assembly[2:25,]
  data.table::setDT(df)
  data.table::setDT(chrDf)
  result <- absolutizeGenomicCoords(df, chrDf)
  shiftedBiomrks <- result[[1]]
  shiftedChr <- result[[2]]

  # Check that all the shifts are correct
  expect_equal(shiftedChr$chrLength,
               shiftedChr$seq_end - shiftedChr$seq_start + 1)
  expect_equal(shiftedBiomrks[chr == "chrX", abs_gene_seq_start],
               shiftedChr[chrName == "X", seq_start] +
               shiftedBiomrks[chr == "chrX", gene_seq_start] - 1)

  # Check for NAs
  expect_false(is.na(shiftedChr[, sum(seq_start) + sum(seq_end)]))
  expect_false(is.na(shiftedBiomrks[, sum(abs_gene_seq_start)]))

  #  Test that there are no warnings about truncations, etc.
  expect_warning(absolutizeGenomicCoords(df, chrDf), regexp = NA)

})
