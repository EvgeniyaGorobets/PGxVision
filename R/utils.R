# A helper function for buildManhattanPlot and buildVolcanoPlot
# Select out all gene drug sensitivity results matching a chosen experiment
# (An experiment is defined by a tissue, a compound, and an mDataType)
selectExperiment <- function(biomarkerDt, tissueName="", compoundName="",
                             molecularType="") {
  # Local bindings to satisfy check() and DT syntax
  tissue <- compound <- mDataType <- NULL

  # Validate user input
  if (!data.table::is.data.table(biomarkerDt)) {
    stop("biomarkerDt must be of type data.table.")
  }
  checkmate::assertString(tissueName)
  checkmate::assertString(compoundName)
  checkmate::assertString(molecularType)

  # Select only biomarkers from the chosen experiment
  selectedBiomrks <- biomarkerDt[tissue == tissueName, ]
  selectedBiomrks <- selectedBiomrks[compound == compoundName, ]
  selectedBiomrks <- selectedBiomrks[mDataType == molecularType, ]

  # Warn if resulting data.table is empty
  if (nrow(selectedBiomrks) == 0) {
    warning(paste0("There were no experiments with the combination ",
      sprintf("tissue=%s, compound=%s, mDataType=%s.\n", tissueName,
      compoundName, molecularType), "Plots will be empty."))
  }

  return(selectedBiomrks)
}

absolutizeGenomicCoords <- function(selectedBiomrks, chromosomeDf) {
  # Local bindings to satisfy check() and DT syntax
  seq_end <- seq_start <- abs_gene_seq_start <- gene_seq_start <- chrLength <-
    chrName <- chr <- NULL

  # Make a copy of the data.tables so you're not modifying user data
  selectedBiomrks <- data.table::copy(selectedBiomrks)
  chromosomeDf <- data.table::copy(chromosomeDf)

  # Add columns to track absolute position in genome (TODO: can maybe remove this)
  chromosomeDf$seq_start <- 1
  chromosomeDf$seq_end <- 1
  selectedBiomrks$abs_gene_seq_start <- 1

  # Shift genomic coords to be relative to genome instead of each chromosome
  chromEnd <- 0  # the end of the previous chromosome
  for (row in 1:nrow(chromosomeDf)) {
    # Set the absolute seq_start & seq_end of each chromosome
    chromosomeDf[row, seq_start := chromEnd + 1]
    chromosomeDf[row, seq_end := chromEnd + chrLength]

    # Update all biomarkers in that chromosome with absolute gene_seq_start
    chromName <- paste("chr", chromosomeDf[row, chrName], sep="") #TODO: not generic enough
    selectedBiomrks[chr == chromName,
                    abs_gene_seq_start := gene_seq_start + chromEnd]

    # Update chromEnd to end of current chromosome
    chromEnd <- chromosomeDf[row, seq_end]
  }

  # Return the modified dfs
  return(list(selectedBiomrks, chromosomeDf))
}

# [END]
