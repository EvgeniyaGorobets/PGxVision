# Possible extra packages
install.packages('dplyr')   # tidyverse data manipulation
install.packages('ggrepel') # for overlapping labels in ggplot2
install.packages("jsonlite") # for processing json files

# Constants
chromosomes <- c(1:22, "X", "Y")

# TODO: example data should be saved as an R obj? why?
biomarkerFile <- "data/pharmacodb_biomarkers_by_tissue.csv"
biomarkerDf <- read.csv(biomarkerFile)

# TODO: determine what part of this should be within the fxn
chromosomeData <- jsonlite::fromJSON("data/chromosome_lengths.json",
                                     simplifyDataFrame = TRUE)
chromosomeDf <- chromosomeData$`Chromosome Info`
colnames(chromosomeDf)[colnames(chromosomeDf) == "value"] <- "total-length"

experiment <- c("Lung", "Trametinib", "rna")
names(experiment) <- c("tissue", "compound", "mDataType")


buildManhattanPlot <- function(biomarkerDf=NULL,
                               chromosomeDf=NULL,
                               experiment=NULL,
                               relativeGenomeCoords=TRUE) {
  # Check user input
  # TODO: parameterize the column names
  # use a renaming map to reduce the # of parameters for a fxn (i.e., a named character vector or list).
  if (is.null(biomarkerDf) | !is.data.frame(biomarkerDf)) {
    stop("Please provide a biomarkerDf of type data.frame")
  }

  # Check that appropriate columns are provided
  # TODO: check columns for other df too (tissue, compoound, mDataType, gene_seq_start, pvalue/fdr)
  if (!is.data.frame(chromosomeDf) |
      !('molecule-name' %in% colnames(chromosomeDf) &
      'total-length' %in% colnames(chromosomeDf))) {
    stop("chromosomeDf must be a data.frame with columns 'molecule-name'
         (corresponding to different chromosomes) and 'total-length'.")
  }

  # Pre-process data
  selectedBiomrks <- selectExperiment(biomarkerDf, experiment)
  # If the provided genomic coordinates are relative, absolute-ize them
  if (relativeGenomeCoords) {
    result <- absolutizeGenomicCoords(selectedBiomrks, chromosomeDf)
    selectedBiomrks <- result[[1]]
    chromosomeDf <- result[[2]]
  }

  # Build the Manhattan plot
  plot <- ggplot(selectedBiomrks, aes(x=abs_gene_seq_start, y=-log10(pvalue)))
  plot <- plot + geom_point()

  # Customize x-axis ticks and labels
  # TODO: add theme!
  midChromosome <- (chromosomeDf$seq_start + chromosomeDf$seq_end)/2
  totalGenomeLen <- chromosomeDf$`total-length`[chromosomeDf$`molecule-name` == "all"]
  plot <- plot + # theme(axis.text.x=element_text(family, face, colour, size)) +
    scale_x_continuous("Chromosome", breaks=midChromosome[-1],
        minor_breaks=c(chromosomeDf$seq_start, chromosomeDf$seq_end),
        labels=chromosomes, limits=c(1, totalGenomeLen)) +
    guides(x = guide_prism_minor())

  return(plot)
}

selectExperiment <- function(biomarkerDf, experiment) {
  # Validate the experiment vector
  if (is.null(experiment) | !("mDataType" %in% names(experiment) &
      "tissue" %in% names(experiment) & "compound" %in% names(experiment))) {
    stop("You must choose an experiment to plot. Please provide a vector with
         a tissue, a compound, and a molecular data type (mDataType).")
  }

  # Select only biomarkers from the chosen experiment
  tissue <- experiment["tissue"]
  selectedBiomrks <- biomarkerDf[biomarkerDf$tissue == tissue, ]
  compound <- experiment["compound"]
  selectedBiomrks <- selectedBiomrks[selectedBiomrks$compound == compound, ]
  mDataType <- experiment["mDataType"]
  selectedBiomrks <- selectedBiomrks[selectedBiomrks$mDataType == mDataType, ]

  # Check that the resulting df is not empty
  if (nrow(selectedBiomrks) == 0) {
    warning(paste0("There were no experiments with the combination ",
                   sprintf("tissue=%s, compound=%s, mDataType=%s.\n",
                           tissue, compound, mDataType),
                   "Manhattan plot will be empty."))
  }

  return(selectedBiomrks)
}


absolutizeGenomicCoords <- function(selectedBiomrks, chromosomeDf) {
  # TODO: update this to use a reference genome?

  # Add columns to track absolute position in genome
  chromosomeDf$seq_start <- 1
  chromosomeDf$seq_end <- chromosomeDf$`total-length`[chromosomeDf$`molecule-name` == "all"]
  selectedBiomrks$abs_gene_seq_start <- selectedBiomrks$gene_seq_start

  # For each chromosome (except chr1), compute its absolute position in the
  # genome and shift biomarker positions accordingly
  for (chr in (seq_along(chromosomes[-1])+1)) {
    # Get the absolute seq_start and seq_end of the chromosome
    prevChrom <- chromosomeDf$`molecule-name` == chromosomes[chr-1]
    prevChromEnd <- chromosomeDf$seq_start[prevChrom] + chromosomeDf$`total-length`[prevChrom]
    chromosomeDf$seq_end[prevChrom] <- prevChromEnd - 1

    currentChrom <- chromosomeDf$`molecule-name` == chromosomes[chr]
    chromosomeDf$seq_start[currentChrom] <- prevChromEnd

    # Update all biomarkers in that chromosome with absolute gene_seq_start
    # TODO: I am modifying in place? should make a copy instead?
    chromName <- paste("chr", chromosomes[chr], sep="")
    relativeSeqStart <- selectedBiomrks$gene_seq_start[selectedBiomrks$chr == chromName]
    absoluteSeqStart <- relativeSeqStart + prevChromEnd - 1
    selectedBiomrks$abs_gene_seq_start[selectedBiomrks$chr == chromName] <- absoluteSeqStart
  }

  # Return the modified dfs
  return(list(selectedBiomrks, chromosomeDf))
}
