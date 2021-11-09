# For development only, load sample data
if (sys.nframe() == 0) {
  # development code here
  # TODO: example data should be saved as an R obj? why?
  biomarkerFile <- "data/pharmacodb_biomarkers_by_tissue.csv"
  biomarkerDf <- read.csv(biomarkerFile)

  # TODO: determine what part of this should be within the fxn
  chromosomeData <- jsonlite::fromJSON("data/chromosome_lengths.json",
                                       simplifyDataFrame = TRUE)
  chromosomeDf <- chromosomeData$`Chromosome Info`
  colnames(chromosomeDf)[colnames(chromosomeDf) == "value"] <- "chrLength"
  colnames(chromosomeDf)[colnames(chromosomeDf) == "molecule-name"] <- "chrName"

  experiment <- c("Lung", "Trametinib", "rna")
  names(experiment) <- c("tissue", "compound", "mDataType")
}

# TODO: move to zzz.R
chromosomes <- c(1:22, "X", "Y")

#' @importFrom data.table setDT
#' @importFrom ggplot2 ggplot, geom_point, scale_x_continuous, guides, theme
#' @importFrom ggprism guide_prism_minor
# TODO: convert internally to data.table
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
      !('chrName' %in% colnames(chromosomeDf) &
      'chrLength' %in% colnames(chromosomeDf))) {
    stop("chromosomeDf must be a data.frame with columns 'chrName'
         (corresponding to different chromosomes) and 'chrLength'.")
  }

  # Convert dfs to data.table by reference
  setDT(biomarkerDf, keep.rownames=TRUE)
  setDT(chromosomeDf, keep.rownames=TRUE)

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
  totalGenomeLen <- chromosomeDf[chrName == "all", chrLength]
  plot <- plot + # theme(axis.text.x=element_text(family, face, colour, size)) +
    scale_x_continuous("Chromosome", breaks=midChromosome[-1],
        minor_breaks=c(chromosomeDf$seq_start, chromosomeDf$seq_end),
        labels=chromosomes, limits=c(1, totalGenomeLen)) +
    guides(x = guide_prism_minor())

  return(plot)
}


absolutizeGenomicCoords <- function(selectedBiomrks, chromosomeDf) {
  # TODO: update this to use a reference genome?

  # Add columns to track absolute position in genome
  chromosomeDf$seq_start <- 1
  chromosomeDf$seq_end <- chromosomeDf[chrName == "all", chrLength]
  selectedBiomrks$abs_gene_seq_start <- selectedBiomrks$gene_seq_start

  # For each chromosome (except chr1), compute its absolute position in the
  # genome and shift biomarker positions accordingly
  for (chr in (seq_along(chromosomes[-1])+1)) {
    # Get the absolute seq_start and seq_end of the chromosome
    prevChrom <- chromosomeDf$chrName == chromosomes[chr-1]
    prevChromEnd <- chromosomeDf$seq_start[prevChrom] + chromosomeDf$chrLength[prevChrom]
    chromosomeDf$seq_end[prevChrom] <- prevChromEnd - 1

    currentChrom <- chromosomeDf$chrName == chromosomes[chr]
    chromosomeDf$seq_start[currentChrom] <- prevChromEnd

    # Update all biomarkers in that chromosome with absolute gene_seq_start
    # TODO: I am modifying in place? should make a copy instead?
    chromName <- paste("chr", chromosomes[chr], sep="")
    relativeSeqStart <- selectedBiomrks[chr == chromName, gene_seq_start]
    absoluteSeqStart <- relativeSeqStart + prevChromEnd - 1
    selectedBiomrks[chr == chromName, abs_gene_seq_start] <- absoluteSeqStart
  }

  # Return the modified dfs
  return(list(selectedBiomrks, chromosomeDf))
}
