#' Build a Manhattan plot of biomarkers
#'
#' @param biomarkerDf A data.frame of drug sensitivity biomarkers that you want
#' to plot, including the compound, tissue, and molecular data type used in the
#' experiment, the genomic coordinates of the genes which were tested on, and
#' the resulting p-value/fdr
#' @param chromosomeDf A data.frame of the reference genome, including the
#' lengths of all the chromosomes. The row order should match the chromosome
#' order in the reference genome. If there is an "all" row, it should be first.
#' @param experiment A named character vector representing the experiment for
#' which you want to plot biomarkers; an experiment is defined by a "tissue",
#' "compound", and "mDataType" (molecular data type)
#' @param relativeGenomeCoords A boolean that indicates whether the genomic
#' coordinates given in biomarkerDf are relative to the chromosome (TRUE) or
#' relative to the entire genome (absolute, FALSE). Default is TRUE.
#' @param genomeName The name of the reference genome (used in plot title).
#' Default value is "GRCh28.p13"
#' @return A ggplot2 plot object mapping the biomarkers of the experiment
#' (x-axis = genome position; y-axis = -log10(p-value or fdr))
#'
#' @examples
#' data(Biomarkers)
#' data(GRCh38.p13.Assembly)
#' experiment <- c("Lung", "Trametinib", "rna")
#' names(experiment) <- c("tissue", "compound", "mDataType")
#' buildManhattanPlot(Biomarkers, GRCh38.p13.Assembly, experiment, TRUE)
#'
#' @importFrom data.table setDT copy
#' @importFrom ggplot2 ggplot geom_point scale_x_continuous guides theme aes
#' scale_color_manual ggtitle element_text
#' @importFrom ggprism guide_prism_minor
#' @export
buildManhattanPlot <- function(biomarkerDf=NULL,
                               chromosomeDf=NULL,
                               experiment=NULL,
                               relativeGenomeCoords=TRUE,
                               genomeName="GRCh38.p13") {
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

  # Select drug sensitivity results from the experiment
  selectedBiomrks <- selectExperiment(biomarkerDf, experiment)

  # If there is an "all" row in the chromosomeDf, remove it
  if (chromosomeDf[1, chrName] == "all") {
    chromosomeDf <- chromosomeDf[2:nrow(chromosomeDf)]
  }

  # If the provided genomic coordinates are relative to the chromosome,
  # shift them to be relative to the genome and compute chromosome coords
  if (relativeGenomeCoords) {
    result <- absolutizeGenomicCoords(selectedBiomrks, chromosomeDf)
    selectedBiomrks <- result[[1]]
    chromosomeDf <- result[[2]]
  }

  # Build the Manhattan plot
  totalGenomeLen <- chromosomeDf[nrow(chromosomeDf), seq_start + chrLength]
  plot <- ggplot(selectedBiomrks, aes(x=abs_gene_seq_start, y=-log10(pvalue),
                                      xmin=1, xmax=totalGenomeLen, color=chr))

  # Add title and colors
  title <- paste0(genomeName, " gene response to ", experiment["compound"],
                  " in ", experiment["tissue"],
                  " tissue\n(drug sensitivity measured using -log10(p-value))")
  plot <- plot + geom_point() + ggtitle(title) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values=rainbow(nrow(chromosomeDf)))

  # Customize x-axis ticks and labels
  # TODO: add theme!
  midChromosome <- (chromosomeDf$seq_start + chromosomeDf$seq_end)/2
  chromosomeNames <- chromosomeDf$chrName
  plot <- plot + # theme(axis.text.x=element_text(family, face, colour, size)) +
    scale_x_continuous("Chromosome", breaks=midChromosome,
        minor_breaks=c(1, chromosomeDf$seq_end), labels=chromosomeNames) +
    guides(x = guide_prism_minor())

  return(plot)
}


absolutizeGenomicCoords <- function(selectedBiomrks, chromosomeDf) {
  # Make a copy of the data.tables so you're not modifying user data
  selectedBiomrks <- copy(selectedBiomrks)
  chromosomeDf <- copy(chromosomeDf)

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
