#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

# Chromosomes
chromosomes <- c(1:22, "X", "Y")


# Volcano Plot Tutorial:
# https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/rna-seq-viz-with-volcanoplot-r/tutorial.html
install.packages('dplyr')   # tidyverse data maniuplation
install.packages('ggplot2') # tidyverse plotting
library('ggplot2')
install.packages('ggprism')
library('ggprism')
install.packages('ggrepel') # for overlapping labels in ggplot2


buildVolcacoPlot <- function() {
  print("Hello, world!")
  # waiting for rnaseq data
}

buildForestPlot <- function() {

}


biomarkerFile <- "data/pharmacodb_biomarkers_by_tissue.csv"

# TODO: determine what part of this should be within the fxn
install.packages("jsonlite")
chromosomeData <- jsonlite::fromJSON("data/chromosome_lengths.json",
                                     simplifyDataFrame = TRUE)
chromosomeDf <- chromosomeData$`Chromosome Info`
colnames(chromosomeDf)[colnames(chromosomeDf) == "value"] <- "total-length"

experiment <- c("Lung", "Trametinib", "rna")
names(experiment) <- c("Tissue", "Compound", "mDataType")

buildManhattanPlot <- function(biomarkerDf=NULL,
                               biomarkerFile=NULL,
                               chromosomeDf=NULL,
                               experiment=NULL,
                               relativeGenomeCoords=TRUE) {
  # Check user input
  # TODO: parameterize the column names
  # use a renaming map to reduce the # of parameters for a fxn (i.e., a named character vector or list).
  if (is.null(biomarkerDf) & is.null(biomarkerFile)) {
    stop("You must provide a data frame or a data file with columns tissue,
         compound, mDataType, seq_start, pvalue/fdr, ...")
  }

  if (is.null(biomarkerDf) & !is.character(biomarkerFile)) {
    stop("biomarkerFile must be a string (the path to your data file)")
  }

  if (is.null(biomarkerFile) & !is.data.frame(biomarkerDf)) {
    stop("biomarkerDf must be of type data.frame")
  }

  # If a biomarker file is provided, read it into a data.frame
  if (is.character(biomarkerFile)) {
    biomarkerDf = read.csv(biomarkerFile)
  }

  # Check that appropriate columns are provided
  # TODO: check columns for other df too; check her pattern for this
  if (!is.data.frame(chromosomeDf) |
      !('molecule-name' %in% colnames(chromosomeDf)) |
      !('total-length' %in% colnames(chromosomeDf))) {
    stop("chromosomeDf must be a data.frame with columns 'molecule-name'
         (corresponding to different chromosomes) and 'total-length'.")
  }

  # Check that an experiment was chosen
  if (is.null(experiment)) {
    stop("You must choose an experiment to plot. Please provide a vector with
         a tissue, a compount, and a molecular data type.")
  }

  # Select only biomarkers from the chosen experiment
  tissue <- experiment["Tissue"]
  selectedBiomrks <- biomarkerDf[biomarkerDf$tissue == tissue, ]
  compound <- experiment["Compound"]
  selectedBiomrks <- selectedBiomrks[selectedBiomrks$compound == compound, ]
  mDataType <- experiment["mDataType"]
  selectedBiomrks <- selectedBiomrks[selectedBiomrks$mDataType == mDataType, ]
  if (nrow(selectedBiomrks) == 0) {
    warning(paste0("There were no experiments with the combination ",
                   sprintf("tissue=%s, compound=%s, mDataType=%s.\n",
                           tissue, compound, mDataType),
                   "Manhattan plot will be empty."))
  }

  # TODO: update this to use a reference genome?
  # If the provided genomic coordinates are relative, absolute-ize them
  # TODO: separate into helper fxn and write a test for it
  if (relativeGenomeCoords) {
    # Add genomic sequence start/end to chromosomeDf
    chromosomeDf$seq_start <- 1
    chromosomeDf$seq_end <- chromosomeDf$`total-length`[1] # FOR NOW
    # Add absoluteSeqStart to biomarkerDf
    selectedBiomrks$abs_gene_seq_start <- selectedBiomrks$gene_seq_start
    for (chr in seq_along(chromosomes)) {
      if (chr > 1) {
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
    }
  }


  # Build the Manhattan plot
  plot <- ggplot(selectedBiomrks, aes(x=abs_gene_seq_start, y=-log10(pvalue))) + geom_point()

  # Customize x-axis ticks and labels
  # TODO: add theme!
  midChromosome <- (chromosomeDf$seq_start + chromosomeDf$seq_end)/2
  totalGenomeLen <- chromosomeDf$`total-length`[chromosomeDf$`molecule-name` == "all"]
  plot <- plot + # theme(axis.text.x=element_text(family, face, colour, size)) +
    scale_x_continuous("Chromosome", breaks=midChromosome[-1],
        minor_breaks=c(chromosomeDf$seq_start, chromosomeDf$seq_end),
        labels=chromosomes, limits=c(1, totalGenomeLen)) +
    guides(x = guide_prism_minor())


}

# So for the x-axis, you will sum all the lengths, then use seq_len(genome_length) to make your x-axis. We want to add minor ticks at the chromosome ends, and major ticks in half way through the chromosome with the name as a label. (edited)
